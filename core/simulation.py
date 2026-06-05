"""
Main simulation engine.

Panel dynamics: each panel solved as an independent SDOF system
  KLM·M·ÿ + C·ẏ + K·y = KL·P(t)·A
  with explicit 4th-order Runge-Kutta via scipy.integrate.solve_ivp.

Building frame: shear-building FEM (n_floors lateral DOFs).
  [M]{ÿ} + [C]{ẏ} + [K]{y} = {F(t)}
  Solved by the same RK45 integrator.

Blast interior propagation:
  When exterior panels of a room fail, the interior blast pressure equals
  the exterior Pso attenuated by a simple transmission factor.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict, Any
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from .blast import BlastSource, friedlander
from .geometry import Panel, Column, Room
from .materials import Material


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

@dataclass
class PanelResult:
    panel_id: int
    time: np.ndarray
    displacement: np.ndarray      # mid-panel out-of-plane [m]
    velocity: np.ndarray
    blast_pressure: np.ndarray    # applied load [kPa]
    peak_pressure: float          # kPa
    peak_displacement: float      # m
    damage_index: float           # peak_disp / yield_disp
    failed: bool
    failure_time: float           # s  (inf if no failure)
    meta: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ColumnResult:
    column_id: int
    shear_demand: float   # N  (from tributary area blast)
    shear_capacity: float # N
    damage_index: float   # demand / capacity
    failed: bool


@dataclass
class BuildingResult:
    time: np.ndarray
    story_disp: np.ndarray   # (n_floors, n_t)  lateral displacement [m]
    story_drift: np.ndarray  # (n_floors,)  peak inter-story drift ratio
    failed_floors: List[int]


@dataclass
class SimulationResult:
    time: np.ndarray
    panel_results: Dict[int, PanelResult]
    column_results: Dict[int, ColumnResult]
    building: BuildingResult
    # Populated after propagation
    room_peak_pressure: Dict[int, float]     # room_id -> kPa
    casualties: Dict[int, Dict[str, float]]  # room_id -> {fatal, severe, minor}
    total_casualties: Dict[str, float]


# ---------------------------------------------------------------------------
# SDOF panel integration
# ---------------------------------------------------------------------------

def _integrate_panel(panel: Panel, blast: BlastSource,
                     t_arr: np.ndarray) -> PanelResult:
    """Solve SDOF equation for one panel."""
    pressure_arr, meta = blast.panel_loading(panel.center, panel.normal, t_arr)
    # Only panels facing toward the blast receive load
    # (panels on sheltered side get Pso/10 approximation)
    r = panel.center - blast.pos
    n = panel.normal
    facing = float(np.dot(-r / max(np.linalg.norm(r), 0.1), n))

    KL = 0.53    # load transformation factor (fixed-fixed)
    Ke = panel.Ke
    Me = panel.Me    # equivalent mass already includes KLM (computed in geometry.py)
    Ce = panel.Ce
    A  = panel.area

    # Build pressure interpolant. For near-field blasts the positive-phase
    # duration td can be shorter than dt_out (0.5 ms), meaning the entire
    # blast pulse falls between two coarse time steps and is never sampled.
    # When td < 10·dt, insert a fine sub-grid over the positive-phase window
    # (td/20 resolution) so the peak is captured correctly.
    facing_factor = 0.1 if facing < -0.1 else 1.0
    pressure_arr  = pressure_arr * facing_factor   # apply once; used by both paths
    dt_out = float(t_arr[1] - t_arr[0]) if len(t_arr) > 1 else 5e-4
    td_s   = meta['td_ms'] * 1e-3
    if td_s < 10.0 * dt_out:
        toa_s   = meta['toa_s']
        dt_fine = max(td_s / 20.0, dt_out / 10.0)
        t_w0    = max(t_arr[0], toa_s)
        t_w1    = min(t_arr[-1], toa_s + td_s * 1.5)
        t_fine  = np.arange(t_w0, t_w1 + dt_fine, dt_fine)
        t_comb  = np.unique(np.concatenate([t_arr, t_fine]))
        t_sh    = (t_comb - toa_s) * 1e3          # shifted to ms for Friedlander
        p_comb  = friedlander(t_sh, meta['Pr'], meta['td_ms'], meta['b']) * facing_factor
        P_interp = interp1d(t_comb, p_comb * 1e3,
                            kind='linear', bounds_error=False, fill_value=0.0)
        peak_p = float(np.max(p_comb))
    else:
        P_interp = interp1d(t_arr, pressure_arr * 1e3,
                            kind='linear', bounds_error=False, fill_value=0.0)
        peak_p = float(np.max(pressure_arr))

    def ode(t, y):
        x, v = y
        F = KL * P_interp(t) * A
        a = (F - Ce * v - Ke * x) / max(Me, 1e-15)
        return [v, a]

    sol = solve_ivp(ode, [t_arr[0], t_arr[-1]], [0.0, 0.0],
                    method='RK45', t_eval=t_arr,
                    rtol=1e-4, atol=1e-7, dense_output=False)

    disp = sol.y[0]
    vel  = sol.y[1]
    peak_d = float(np.max(np.abs(disp)))
    di     = peak_d / max(panel.yield_disp, 1e-12)
    failed = peak_d >= panel.fail_disp

    t_fail = np.inf
    if failed:
        mask = np.abs(disp) >= panel.fail_disp
        if mask.any():
            t_fail = float(t_arr[np.argmax(mask)])

    return PanelResult(
        panel_id=panel.id,
        time=t_arr.copy(),
        displacement=disp,
        velocity=vel,
        blast_pressure=pressure_arr,
        peak_pressure=peak_p,
        peak_displacement=peak_d,
        damage_index=di,
        failed=failed,
        failure_time=t_fail,
        meta=meta,
    )


# ---------------------------------------------------------------------------
# Shear-building FEM
# ---------------------------------------------------------------------------

def _shear_building_fem(columns: List[Column], panels: List[Panel],
                         blast: BlastSource, t_arr: np.ndarray,
                         n_floors: int, floor_height: float,
                         floor_mass_per_m2: float = 500.0,
                         building_width: float = 12.0,
                         building_depth: float = 10.0) -> BuildingResult:
    """
    Lateral shear-building model.
    DOF: horizontal (x) displacement at each floor level.
    """
    if n_floors < 1:
        return BuildingResult(t_arr, np.zeros((0, len(t_arr))),
                              np.array([]), [])

    # Inter-story stiffness: sum of 12·EI/h^3 for columns in that story
    def story_stiffness(story: int) -> float:
        # Columns belonging to this story: base.z within floor band
        z_lo = story * floor_height - 0.01
        z_hi = (story + 1) * floor_height + 0.01
        k = 0.0
        for col in columns:
            if z_lo <= col.base[2] < z_hi:
                k += 12.0 * col.EI / col.height ** 3
        return max(k, 1e6)  # floor if no columns found

    # Floor masses
    floor_area = building_width * building_depth
    M_floor = floor_mass_per_m2 * floor_area   # kg per floor (approx)

    m_vec = np.full(n_floors, M_floor)
    k_vec = np.array([story_stiffness(f) for f in range(n_floors)])

    # Assemble tridiagonal stiffness matrix
    K = np.zeros((n_floors, n_floors))
    for i in range(n_floors):
        K[i, i] += k_vec[i]
        if i < n_floors - 1:
            K[i, i]     += k_vec[i+1]
            K[i, i+1]   -= k_vec[i+1]
            K[i+1, i]   -= k_vec[i+1]
    M = np.diag(m_vec)

    # Rayleigh damping: 5% at 1st and 3rd modes (approx)
    # Generalised eigenvalue via M^(-1/2) K M^(-1/2) (M is diagonal)
    m_inv_sqrt = 1.0 / np.sqrt(np.maximum(m_vec, 1e-12))
    K_tilde = (K * m_inv_sqrt).T * m_inv_sqrt
    omegas = np.sqrt(np.maximum(np.linalg.eigvalsh(K_tilde), 0.0))
    omegas = np.clip(omegas, 0.1, None)
    w1, w3 = omegas[0], omegas[min(2, len(omegas)-1)]
    alpha = 2*w1*w3*0.05 / (w1 + w3)
    beta  = 2*0.05 / (w1 + w3)
    C = alpha * M + beta * K

    # Lateral blast force at each floor: sum of horizontal blast forces on
    # exterior panels at that floor, projected onto x-direction
    F_floors = np.zeros((n_floors, len(t_arr)))
    for p in panels:
        if p.panel_type in ('ext_wall', 'window') and abs(p.normal[0]) > 0.5:
            prs, _ = blast.panel_loading(p.center, p.normal, t_arr)
            force = prs * 1e3 * p.area * abs(p.normal[0])   # N
            f_idx = p.floor_idx
            if 0 <= f_idx < n_floors:
                F_floors[f_idx] += force

    F_interp = interp1d(t_arr, F_floors, axis=1,
                        kind='linear', bounds_error=False, fill_value=0.0)

    M_inv = np.diag(1.0 / np.diag(M))   # diagonal mass

    def ode(t, y):
        u  = y[:n_floors]
        du = y[n_floors:]
        F  = F_interp(t)
        a  = M_inv @ (F - C @ du - K @ u)
        return np.concatenate([du, a])

    y0  = np.zeros(2 * n_floors)
    sol = solve_ivp(ode, [t_arr[0], t_arr[-1]], y0,
                    method='RK45', t_eval=t_arr,
                    rtol=1e-4, atol=1e-6)

    story_disp = sol.y[:n_floors, :]  # (n_floors, n_t)

    # Inter-story drift ratio
    drift = np.zeros(n_floors)
    for i in range(n_floors):
        if i == 0:
            delta = np.max(np.abs(story_disp[0]))
        else:
            delta = np.max(np.abs(story_disp[i] - story_disp[i-1]))
        drift[i] = delta / floor_height

    failed_floors = [i for i, dr in enumerate(drift) if dr > 0.05]

    return BuildingResult(
        time=t_arr.copy(),
        story_disp=story_disp,
        story_drift=drift,
        failed_floors=failed_floors,
    )


# ---------------------------------------------------------------------------
# Interior blast propagation
# ---------------------------------------------------------------------------

def _interior_pressure(rooms: List[Room], panels: List[Panel],
                       panel_results: Dict[int, PanelResult],
                       blast: BlastSource) -> Dict[int, float]:
    """
    Estimate peak interior blast pressure for each room.
    If an exterior panel of a room has failed, assume blast penetration.
    Transmission factor for failed wall: 0.5 (reflected off interior).
    Transmission factor for intact exterior wall: 0.05 (leakage / cracks).
    """
    room_pressure: Dict[int, float] = {}
    for room in rooms:
        p_in = 0.0
        for pid in room.panel_ids:
            pr = panel_results.get(pid)
            if pr is None:
                continue
            panel = panels[pid]
            is_ext = panel.panel_type in ('ext_wall', 'window')
            if is_ext:
                # Use incident side-on Pso (not reflected Pr) for interior pressure
                pso_val = pr.meta.get('Pso', pr.peak_pressure)
                factor = 0.5 if pr.failed else 0.05
                p_in = max(p_in, pso_val * factor)
            else:
                # Interior panel – propagate from adjacent room if it failed
                # (simplified: just carry 80% of neighbour's interior pressure)
                # This is handled in a second pass below
                pass
        room_pressure[room.id] = p_in

    # Second pass: propagate through failed interior walls
    for panel in panels:
        if panel.panel_type != 'int_wall':
            continue
        pr = panel_results.get(panel.id)
        if pr is None or not pr.failed:
            continue
        ri = panel.room_inside
        ro = panel.room_outside
        if ri >= 0 and ro >= 0:
            p_max = max(room_pressure.get(ri, 0.0), room_pressure.get(ro, 0.0))
            room_pressure[ri] = max(room_pressure.get(ri, 0.0), p_max * 0.8)
            room_pressure[ro] = max(room_pressure.get(ro, 0.0), p_max * 0.8)

    return room_pressure


# ---------------------------------------------------------------------------
# Column checks
# ---------------------------------------------------------------------------

def _check_columns(columns: List[Column], panels: List[Panel],
                   blast: BlastSource, t_arr: np.ndarray,
                   floor_height: float) -> Dict[int, ColumnResult]:
    results: Dict[int, ColumnResult] = {}
    for col in columns:
        # Tributary area: area of adjacent exterior wall panels
        # Approximate by collecting panels on same floor
        trib_area = 0.0
        for p in panels:
            if (p.panel_type == 'ext_wall' and p.floor_idx == 0):
                trib_area += p.area / max(len([c for c in columns]), 1)
                break
        # Conservative: use peak pressure × tributary area as shear demand
        r = col.base - blast.pos
        R = max(np.linalg.norm(r), 0.1)
        from .blast import scaled_distance, peak_overpressure
        from .materials import CONCRETE
        Z = np.clip(scaled_distance(R, blast.W_eff), 0.05, 40.0)
        Pso = peak_overpressure(Z) * 1e3   # Pa
        # Tributary area approximation
        trib = floor_height * min(col.b * 5, 3.0)   # ~5 col widths
        V_demand = Pso * trib
        DI = V_demand / max(col.V_cap, 1.0)
        results[col.id] = ColumnResult(
            column_id=col.id,
            shear_demand=V_demand,
            shear_capacity=col.V_cap,
            damage_index=DI,
            failed=DI > 1.0,
        )
    return results


# ---------------------------------------------------------------------------
# Top-level runner
# ---------------------------------------------------------------------------

def run_simulation(
    blast: BlastSource,
    panels: List[Panel],
    columns: List[Column],
    rooms: List[Room],
    n_floors: int,
    floor_height: float,
    building_width: float,
    building_depth: float,
    t_end: float = 0.25,      # simulation duration [s]
    dt: float = 5e-4,         # output time step [s]
) -> SimulationResult:
    from .casualties import compute_room_casualties

    t_arr = np.arange(0.0, t_end + dt, dt)

    # 1. Panel SDOF
    panel_results: Dict[int, PanelResult] = {}
    for p in panels:
        panel_results[p.id] = _integrate_panel(p, blast, t_arr)
        # Write back to panel state
        pr = panel_results[p.id]
        p.max_disp    = pr.peak_displacement
        p.damage_index = pr.damage_index
        p.failed      = pr.failed

    # 2. Building frame FEM
    bldg = _shear_building_fem(columns, panels, blast, t_arr,
                                n_floors, floor_height,
                                building_width=building_width,
                                building_depth=building_depth)

    # 3. Column checks
    col_results = _check_columns(columns, panels, blast, t_arr, floor_height)
    for col in columns:
        if col.id in col_results:
            col.damage_index = col_results[col.id].damage_index
            col.failed       = col_results[col.id].failed

    # 4. Interior blast
    room_pressure = _interior_pressure(rooms, panels, panel_results, blast)

    # 5. Casualties
    from .casualties import compute_room_casualties
    casualties: Dict[int, Dict[str, float]] = {}
    for room in rooms:
        prom = room_pressure.get(room.id, 0.0)
        room.interior_peak_pressure = prom
        c = compute_room_casualties(prom, room.n_occupants)
        room.casualties_fatal  = c['fatal']
        room.casualties_severe = c['severe']
        room.casualties_minor  = c['minor']
        casualties[room.id] = c

    totals = {
        'fatal':  sum(c['fatal']  for c in casualties.values()),
        'severe': sum(c['severe'] for c in casualties.values()),
        'minor':  sum(c['minor']  for c in casualties.values()),
    }

    return SimulationResult(
        time=t_arr,
        panel_results=panel_results,
        column_results=col_results,
        building=bldg,
        room_peak_pressure=room_pressure,
        casualties=casualties,
        total_casualties=totals,
    )
