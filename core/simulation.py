"""
Main simulation engine.

Panel dynamics: each panel solved as an independent SDOF system
  KLM·M·ÿ + C·ẏ + K·y = P(t)·A
  with explicit 4th-order Runge-Kutta via scipy.integrate.solve_ivp.
  (Biggs equivalent system: KM·M·ÿ + KL·K·y = KL·F; dividing by KL gives
  KLM·M·ÿ + K·y = F with KLM = KM/KL, so the load is NOT factored again.)

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

from .blast import BlastSource
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
    """Solve SDOF equation for one panel.

    The response is computed in three exact phases so the short blast pulse
    can never be skipped by the adaptive ODE integrator:

      1. Before arrival (t < toa): force and response are identically zero.
      2. During the pulse (toa .. toa+td): solve_ivp with max_step = td/20
         and the Friedlander force evaluated directly (no interpolant).
      3. After the pulse: closed-form damped free vibration from the
         end-of-pulse state — this is where the peak usually occurs and the
         analytic form captures it exactly.

    An earlier implementation integrated the full window with default
    (unbounded) max_step: RK45 grew its step size over the quiet pre-arrival
    phase and frequently leapt over the entire pulse, so whether a panel got
    loaded at all depended on where the adaptive steps happened to land.
    """
    pressure_arr, meta = blast.panel_loading(panel.center, panel.normal, t_arr)
    # Only panels facing toward the blast receive load
    # (panels on sheltered side get Pso/10 approximation)
    r = panel.center - blast.pos
    n = panel.normal
    facing = float(np.dot(-r / max(np.linalg.norm(r), 0.1), n))
    facing_factor = 0.1 if facing < -0.1 else 1.0
    pressure_arr  = pressure_arr * facing_factor

    # Biggs load-mass form: Me = KLM·M with KLM = KM/KL (set in geometry.py),
    # so the applied force is the *unfactored* P(t)·A. Multiplying the load
    # by KL here as well would double-count the transformation.
    Ke = panel.Ke
    Me = max(panel.Me, 1e-15)
    Ce = panel.Ce
    A  = panel.area

    toa   = float(meta['toa_s'])
    td_ms = float(meta['td_ms'])
    td_s  = td_ms * 1e-3
    t_end = float(t_arr[-1])
    Pr_eff = float(meta['Pr']) * facing_factor

    disp   = np.zeros_like(t_arr)
    vel    = np.zeros_like(t_arr)
    peak_d = 0.0
    peak_p = Pr_eff if toa < t_end else 0.0   # Friedlander peak is at arrival
    pulse_end = min(toa + td_s, t_end)

    if toa < t_end and td_s > 0.0:
        # ── Phase 2: forced response over the pulse window ──────────────
        b_wave = meta['b']

        def ode(t, y):
            x, v = y
            tau_ms = (t - toa) * 1e3
            p_kPa = Pr_eff * (1.0 - tau_ms / td_ms) * np.exp(-b_wave * tau_ms / td_ms) \
                    if 0.0 <= tau_ms <= td_ms else 0.0
            F = p_kPa * 1e3 * A
            return [v, (F - Ce * v - Ke * x) / Me]

        sol = solve_ivp(ode, [toa, pulse_end], [0.0, 0.0],
                        method='RK45', max_step=max(td_s / 20.0, 1e-7),
                        rtol=1e-5, atol=1e-9, dense_output=True)
        x1 = float(sol.y[0, -1])
        v1 = float(sol.y[1, -1])
        peak_d = float(np.max(np.abs(sol.y[0])))

        in_pulse = (t_arr >= toa) & (t_arr <= pulse_end)
        if in_pulse.any():
            y_p = sol.sol(t_arr[in_pulse])
            disp[in_pulse] = y_p[0]
            vel[in_pulse]  = y_p[1]

        # ── Phase 3: analytic damped free vibration after the pulse ─────
        wn   = np.sqrt(Ke / Me)
        zeta = min(Ce / max(2.0 * np.sqrt(Ke * Me), 1e-15), 0.999)
        wd   = wn * np.sqrt(1.0 - zeta ** 2)
        C1   = x1
        C2   = (v1 + zeta * wn * x1) / max(wd, 1e-12)

        def free_resp(tau):
            e = np.exp(-zeta * wn * tau)
            x = e * (C1 * np.cos(wd * tau) + C2 * np.sin(wd * tau))
            v = e * ((C2 * wd - zeta * wn * C1) * np.cos(wd * tau)
                     - (C1 * wd + zeta * wn * C2) * np.sin(wd * tau))
            return x, v

        after = t_arr > pulse_end
        if after.any():
            xf, vf = free_resp(t_arr[after] - pulse_end)
            disp[after] = xf
            vel[after]  = vf

        # Peak of the free phase occurs within the first natural period;
        # sample it finely so peak displacement is not limited by dt_out.
        if t_end > pulse_end and wd > 0.0:
            tau_f = np.linspace(0.0, min(2.0 * np.pi / wd, t_end - pulse_end), 256)
            xf, _ = free_resp(tau_f)
            peak_d = max(peak_d, float(np.max(np.abs(xf))))

    di     = peak_d / max(panel.yield_disp, 1e-12)
    failed = peak_d >= panel.fail_disp

    t_fail = np.inf
    if failed:
        mask = np.abs(disp) >= panel.fail_disp
        # Fine-grid peak may cross the threshold between coarse samples;
        # fall back to end-of-pulse as the failure time in that case.
        t_fail = float(t_arr[np.argmax(mask)]) if mask.any() else float(pulse_end)

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
    dt_out = float(t_arr[1] - t_arr[0]) if len(t_arr) > 1 else 5e-4
    # max_step bounds the adaptive stepper so it cannot leap over the blast
    # pulse during the quiet pre-arrival phase (same failure mode as the
    # panel SDOF solver had).
    sol = solve_ivp(ode, [t_arr[0], t_arr[-1]], y0,
                    method='RK45', t_eval=t_arr, max_step=dt_out,
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
    from .blast import scaled_distance, peak_overpressure

    results: Dict[int, ColumnResult] = {}
    for col in columns:
        # Conservative: use peak pressure × tributary area as shear demand
        r = col.base - blast.pos
        R = max(np.linalg.norm(r), 0.1)
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
