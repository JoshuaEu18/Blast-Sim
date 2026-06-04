"""
Parametric building geometry generator.

Produces Panels (walls, floors, windows), Columns, and Rooms.
Each Panel gets SDOF structural properties computed from plate theory
(Timoshenko, simply-supported boundary conditions with a fixed-edge factor).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple
import numpy as np

from .materials import Material, GLASS


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Panel:
    id: int
    panel_type: str          # 'ext_wall' | 'int_wall' | 'floor' | 'roof' | 'window'
    corners: np.ndarray      # (4, 3) CCW from outside view
    center: np.ndarray       # (3,)
    normal: np.ndarray       # (3,) outward unit normal
    width: float             # panel shorter horizontal/span dimension [m]
    height: float            # panel taller / span dimension [m]
    area: float              # [m²]
    thickness: float         # [m]
    material: Material
    floor_idx: int
    room_inside: int         # room ID on interior side, -1 if exterior
    room_outside: int        # room ID on exterior side, -1 if exterior

    # SDOF structural properties (computed by _sdof_properties)
    omega_n: float = 0.0     # natural angular frequency [rad/s]
    Ke: float = 0.0          # equivalent stiffness [N/m]
    Me: float = 0.0          # equivalent mass [kg]
    Ce: float = 0.0          # equivalent damping [N·s/m]
    yield_disp: float = 0.0  # mid-panel displacement at first yield [m]
    fail_disp: float = 0.0   # displacement at failure [m]

    # Simulation state (written by simulation.py)
    max_disp: float = 0.0
    damage_index: float = 0.0   # max_disp / yield_disp  (>1 = yielded, >μ = failed)
    failed: bool = False


@dataclass
class Column:
    id: int
    base: np.ndarray         # (3,) bottom [m]
    top: np.ndarray          # (3,) top [m]
    height: float            # [m]
    material: Material
    b: float                 # cross-section width [m]
    d: float                 # cross-section depth [m]

    # Capacities
    M_cap: float = 0.0       # plastic moment capacity [N·m]
    V_cap: float = 0.0       # shear capacity [N]
    P_cap: float = 0.0       # axial capacity [N]
    EI: float = 0.0          # flexural rigidity [N·m²]

    # State
    damage_index: float = 0.0
    failed: bool = False


@dataclass
class Room:
    id: int
    floor_idx: int
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    z_bot: float
    z_top: float
    n_occupants: int
    panel_ids: List[int] = field(default_factory=list)

    # Simulation results
    interior_peak_pressure: float = 0.0  # kPa
    casualties_fatal: float = 0.0
    casualties_severe: float = 0.0
    casualties_minor: float = 0.0


# ---------------------------------------------------------------------------
# Plate theory SDOF properties
# ---------------------------------------------------------------------------

# Timoshenko Table: simply-supported plate, b >= a (short side a)
# b/a :  1.0     1.2     1.4     1.6     1.8     2.0     3.0     ∞
_AR_TABLE   = [1.0,   1.2,    1.4,    1.6,    1.8,    2.0,    3.0,   10.0]
_ALPHA_DEFL = [0.00406, 0.00564, 0.00698, 0.00830, 0.00931, 0.01013, 0.01223, 0.01302]
_BETA_MOM   = [0.0479, 0.0627,  0.0755,  0.0862,  0.0948,  0.1017,  0.1189,  0.1250]

# UFC 3-340-02 Table 3-12 equivalent KLM for simply-supported plate, uniform load
_KLM_SS = 0.58   # 2-way spanning
_KLM_FF = 0.52   # fixed on all 4 edges

# Boundary condition stiffness correction factor (fixed vs SS):
# fixed plate has ~4× higher peak stiffness than SS plate
_BC_STIFFNESS_FACTOR = 2.5   # intermediate (partly fixed)
_BC_MOMENT_FACTOR    = 0.6   # reduces yield moment demand at centre


def _sdof_properties(panel: Panel) -> Panel:
    """Compute and attach SDOF structural parameters to a Panel in-place."""
    mat = panel.material
    t   = panel.thickness
    a   = min(panel.width, panel.height)   # short side
    b   = max(panel.width, panel.height)   # long side
    ar  = min(b / max(a, 1e-4), 10.0)

    D = mat.E * t ** 3 / (12.0 * (1.0 - mat.nu ** 2))   # flexural rigidity [N·m]

    alpha_s = float(np.interp(ar, _AR_TABLE, _ALPHA_DEFL))
    beta_m  = float(np.interp(ar, _AR_TABLE, _BETA_MOM))

    # Equivalent stiffness: k = F/w where F = p*A, w = alpha*p*a^4/D
    # k_equiv = D * A / (alpha * a^4)  [N/m] (force = p*A, displacement = p*a^4*alpha/D)
    A = panel.area
    k_ss = D * A / (alpha_s * a ** 4)
    Ke   = k_ss * _BC_STIFFNESS_FACTOR    # corrected for boundary conditions

    # Equivalent mass
    m_total = mat.rho * A * t
    KLM = (_KLM_SS + _KLM_FF) / 2.0      # average BC factor
    Me  = KLM * m_total

    # Damping (viscous)
    Ce = 2.0 * mat.damping_ratio * np.sqrt(Ke * Me)

    omega_n = np.sqrt(Ke / max(Me, 1e-12))

    # Yield displacement: mid-panel reaches yield strain
    # Max moment at centre: M = beta * p * a^2
    # Yield pressure: p_y = fy * t^2 / (6 * beta * a^2) * BC_factor
    beta_eff = beta_m * _BC_MOMENT_FACTOR
    p_yield = mat.fy * t ** 2 / max(6.0 * beta_eff * a ** 2, 1e-20)
    w_yield = alpha_s * p_yield * a ** 4 / D / _BC_STIFFNESS_FACTOR
    w_yield = max(w_yield, 1e-6)

    panel.Ke         = Ke
    panel.Me         = Me
    panel.Ce         = Ce
    panel.omega_n    = omega_n
    panel.yield_disp = w_yield
    panel.fail_disp  = w_yield * mat.ductility
    return panel


# ---------------------------------------------------------------------------
# Column capacity
# ---------------------------------------------------------------------------

def _column_capacities(col: Column) -> Column:
    mat = col.material
    A   = col.b * col.d
    Ix  = col.b * col.d ** 3 / 12.0
    col.EI    = mat.E * Ix
    col.P_cap = mat.fy * A
    col.V_cap = 0.6 * mat.fy * A       # shear yielding
    # Plastic moment: Z = b*d^2/4
    Z_pl = col.b * col.d ** 2 / 4.0
    col.M_cap = mat.fy * Z_pl
    return col


# ---------------------------------------------------------------------------
# Building generator
# ---------------------------------------------------------------------------

def create_building(
    width: float = 12.0,          # x direction [m]
    depth: float = 10.0,          # y direction [m]
    n_floors: int = 3,
    floor_height: float = 3.0,    # [m]
    wall_thickness: float = 0.20, # exterior walls [m]
    floor_thickness: float = 0.25,# [m]
    col_size: float = 0.40,       # square column side [m]
    wall_mat: Material = None,
    floor_mat: Material = None,
    col_mat: Material = None,
    glass_mat: Material = None,
    n_rooms_x: int = 2,
    n_rooms_y: int = 2,
    window_frac: float = 0.30,    # fraction of exterior wall area glazed
    total_occupants: int = 100,
) -> Tuple[List[Panel], List[Column], List[Room]]:
    from .materials import CONCRETE, STEEL_STRUCTURAL, GLASS as _GLASS

    wall_mat  = wall_mat  or CONCRETE
    floor_mat = floor_mat or CONCRETE
    col_mat   = col_mat   or STEEL_STRUCTURAL
    glass_mat = glass_mat or _GLASS

    panels: List[Panel]  = []
    columns: List[Column] = []
    rooms: List[Room]    = []

    pid = 0  # panel counter
    cid = 0  # column counter
    rid = 0  # room counter

    total_rooms = n_floors * n_rooms_x * n_rooms_y
    occ_per_room = max(1, total_occupants // total_rooms)

    # Grid lines in x and y for room subdivision
    xs = np.linspace(0.0, width, n_rooms_x + 1)
    ys = np.linspace(0.0, depth, n_rooms_y + 1)

    for f in range(n_floors):
        z_bot = f * floor_height
        z_top = (f + 1) * floor_height

        # ---- Rooms --------------------------------------------------------
        floor_rooms: dict[tuple, int] = {}  # (ix, iy) -> room_id
        for ix in range(n_rooms_x):
            for iy in range(n_rooms_y):
                rm = Room(
                    id=rid, floor_idx=f,
                    x_min=xs[ix], x_max=xs[ix+1],
                    y_min=ys[iy], y_max=ys[iy+1],
                    z_bot=z_bot, z_top=z_top,
                    n_occupants=occ_per_room,
                )
                rooms.append(rm)
                floor_rooms[(ix, iy)] = rid
                rid += 1

        # ---- Exterior walls -----------------------------------------------
        def _ext_wall(corners4, nrm, w, h, thick, mat,
                      room_in, win_frac, wall_type='ext_wall'):
            nonlocal pid
            center = corners4.mean(axis=0)
            area_total = w * h

            # Solid wall panel (wall_frac of area, same geometry — area used for force)
            wall_area = area_total * (1.0 - win_frac)
            p = Panel(
                id=pid, panel_type=wall_type,
                corners=corners4.copy(), center=center.copy(),
                normal=nrm.copy(), width=w, height=h,
                area=wall_area, thickness=thick,
                material=mat, floor_idx=f,
                room_inside=room_in, room_outside=-1,
            )
            _sdof_properties(p)
            panels.append(p)
            pid += 1

            # Window panel
            if win_frac > 0.0:
                win_area = area_total * win_frac
                win_w = w * np.sqrt(win_frac)
                win_h = h * np.sqrt(win_frac)
                pw = Panel(
                    id=pid, panel_type='window',
                    corners=corners4.copy(), center=center.copy(),
                    normal=nrm.copy(), width=win_w, height=win_h,
                    area=win_area, thickness=0.006,   # 6 mm glazing
                    material=glass_mat, floor_idx=f,
                    room_inside=room_in, room_outside=-1,
                )
                _sdof_properties(pw)
                panels.append(pw)
                pid += 1

        # South wall (y=0, normal=[0,-1,0]): rooms iy=0
        for ix in range(n_rooms_x):
            w = xs[ix+1] - xs[ix]
            corners = np.array([
                [xs[ix],   0.0, z_bot],
                [xs[ix+1], 0.0, z_bot],
                [xs[ix+1], 0.0, z_top],
                [xs[ix],   0.0, z_top],
            ])
            _ext_wall(corners, np.array([0.,-1.,0.]), w, floor_height,
                      wall_thickness, wall_mat,
                      floor_rooms[(ix, 0)], window_frac)

        # North wall (y=depth, normal=[0,1,0]): rooms iy=n_rooms_y-1
        for ix in range(n_rooms_x):
            w = xs[ix+1] - xs[ix]
            corners = np.array([
                [xs[ix+1], depth, z_bot],
                [xs[ix],   depth, z_bot],
                [xs[ix],   depth, z_top],
                [xs[ix+1], depth, z_top],
            ])
            _ext_wall(corners, np.array([0.,1.,0.]), w, floor_height,
                      wall_thickness, wall_mat,
                      floor_rooms[(ix, n_rooms_y-1)], window_frac)

        # West wall (x=0, normal=[-1,0,0]): rooms ix=0
        for iy in range(n_rooms_y):
            h_span = ys[iy+1] - ys[iy]
            corners = np.array([
                [0., ys[iy],   z_bot],
                [0., ys[iy+1], z_bot],
                [0., ys[iy+1], z_top],
                [0., ys[iy],   z_top],
            ])
            _ext_wall(corners, np.array([-1.,0.,0.]), h_span, floor_height,
                      wall_thickness, wall_mat,
                      floor_rooms[(0, iy)], window_frac)

        # East wall (x=width, normal=[1,0,0]): rooms ix=n_rooms_x-1
        for iy in range(n_rooms_y):
            h_span = ys[iy+1] - ys[iy]
            corners = np.array([
                [width, ys[iy+1], z_bot],
                [width, ys[iy],   z_bot],
                [width, ys[iy],   z_top],
                [width, ys[iy+1], z_top],
            ])
            _ext_wall(corners, np.array([1.,0.,0.]), h_span, floor_height,
                      wall_thickness, wall_mat,
                      floor_rooms[(n_rooms_x-1, iy)], window_frac)

        # ---- Interior walls -----------------------------------------------
        for ix in range(1, n_rooms_x):   # vertical partitions (x = xs[ix])
            for iy in range(n_rooms_y):
                h_span = ys[iy+1] - ys[iy]
                corners = np.array([
                    [xs[ix], ys[iy],   z_bot],
                    [xs[ix], ys[iy+1], z_bot],
                    [xs[ix], ys[iy+1], z_top],
                    [xs[ix], ys[iy],   z_top],
                ])
                center = corners.mean(axis=0)
                p = Panel(
                    id=pid, panel_type='int_wall',
                    corners=corners, center=center,
                    normal=np.array([1.,0.,0.]),
                    width=h_span, height=floor_height,
                    area=h_span * floor_height,
                    thickness=wall_thickness * 0.8,
                    material=wall_mat, floor_idx=f,
                    room_inside=floor_rooms[(ix, iy)],
                    room_outside=floor_rooms[(ix-1, iy)],
                )
                _sdof_properties(p)
                panels.append(p)
                pid += 1

        for iy in range(1, n_rooms_y):   # horizontal partitions (y = ys[iy])
            for ix in range(n_rooms_x):
                w = xs[ix+1] - xs[ix]
                corners = np.array([
                    [xs[ix],   ys[iy], z_bot],
                    [xs[ix+1], ys[iy], z_bot],
                    [xs[ix+1], ys[iy], z_top],
                    [xs[ix],   ys[iy], z_top],
                ])
                center = corners.mean(axis=0)
                p = Panel(
                    id=pid, panel_type='int_wall',
                    corners=corners, center=center,
                    normal=np.array([0.,1.,0.]),
                    width=w, height=floor_height,
                    area=w * floor_height,
                    thickness=wall_thickness * 0.8,
                    material=wall_mat, floor_idx=f,
                    room_inside=floor_rooms[(ix, iy)],
                    room_outside=floor_rooms[(ix, iy-1)],
                )
                _sdof_properties(p)
                panels.append(p)
                pid += 1

        # ---- Floor slab (skip ground floor bottom) ------------------------
        if f > 0:
            p = Panel(
                id=pid, panel_type='floor',
                corners=np.array([[0,0,z_bot],[width,0,z_bot],
                                  [width,depth,z_bot],[0,depth,z_bot]],dtype=float),
                center=np.array([width/2, depth/2, z_bot]),
                normal=np.array([0.,0.,-1.]),
                width=width, height=depth,
                area=width*depth, thickness=floor_thickness,
                material=floor_mat, floor_idx=f,
                room_inside=-1, room_outside=-1,
            )
            _sdof_properties(p)
            panels.append(p)
            pid += 1

        # ---- Roof ----------------------------------------------------------
        if f == n_floors - 1:
            p = Panel(
                id=pid, panel_type='roof',
                corners=np.array([[0,0,z_top],[width,0,z_top],
                                  [width,depth,z_top],[0,depth,z_top]],dtype=float),
                center=np.array([width/2, depth/2, z_top]),
                normal=np.array([0.,0.,1.]),
                width=width, height=depth,
                area=width*depth, thickness=floor_thickness,
                material=floor_mat, floor_idx=f,
                room_inside=-1, room_outside=-1,
            )
            _sdof_properties(p)
            panels.append(p)
            pid += 1

        # ---- Columns -------------------------------------------------------
        for xi in xs:
            for yi in ys:
                col = Column(
                    id=cid,
                    base=np.array([xi, yi, z_bot]),
                    top=np.array([xi, yi, z_top]),
                    height=floor_height,
                    material=col_mat,
                    b=col_size, d=col_size,
                )
                _column_capacities(col)
                columns.append(col)
                cid += 1

    # Assign panels to rooms
    for p in panels:
        if p.room_inside >= 0:
            rooms[p.room_inside].panel_ids.append(p.id)
        if p.room_outside >= 0:
            rooms[p.room_outside].panel_ids.append(p.id)

    return panels, columns, rooms
