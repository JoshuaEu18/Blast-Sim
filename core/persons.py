"""
Individual person injury assessment for blast loading.

Computes:
  1. Effective blast loading at each person's position (with wall shielding
     via a ray-panel intersection test)
  2. Primary blast injuries  (overpressure-driven, Baker 1983 Probit)
  3. Tertiary injuries       (whole-body displacement from impulse)
  4. Secondary injuries      (debris / fragments from nearby failed panels)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple
import numpy as np
from scipy.special import ndtr

from .blast import (BlastSource, scaled_distance,
                    peak_overpressure, positive_duration, specific_impulse)

_FRONTAL_AREA = 0.7   # m²  representative frontal area of standing person


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Person:
    id: int
    name: str
    x: float
    y: float
    z: float           # absolute Z in building coords (floor_idx * fh + ~1 m standing)
    mass_kg: float = 70.0


@dataclass
class InjuryDetail:
    name: str
    probability: float  # 0–1
    severity: str       # None / Low / Mild / Moderate / Severe / Fatal
    mechanism: str


@dataclass
class PersonInjury:
    person_id: int
    name: str
    position: Tuple[float, float, float]

    standoff_m: float = 0.0
    raw_pressure_kPa: float = 0.0
    effective_pressure_kPa: float = 0.0
    impulse_kPa_ms: float = 0.0
    shielding_factor: float = 1.0
    walls_between: int = 0

    injuries: List[InjuryDetail] = field(default_factory=list)
    throw_velocity_ms: float = 0.0
    nearest_failed_m: float = float('inf')

    overall_severity: str = 'Uninjured'
    survival_probability: float = 1.0


# ---------------------------------------------------------------------------
# Shielding  (ray–panel intersection)
# ---------------------------------------------------------------------------

# Fraction of overpressure transmitted through one wall of each material type
_TRANSMISSION = {
    'Reinforced Concrete': 0.12,
    'Brick Masonry':       0.22,
    'Steel Plate':         0.10,
    'Annealed Glass':      0.70,
    'Structural Steel':    0.10,
}
_FAILED_TRANSMISSION = 0.88   # failed panel ≈ open gap


def _ray_triangle(origin: np.ndarray, direction: np.ndarray,
                  v0: np.ndarray, v1: np.ndarray, v2: np.ndarray):
    """Möller–Trumbore ray-triangle test. Returns distance t or None."""
    EPS = 1e-9
    e1 = v1 - v0
    e2 = v2 - v0
    h  = np.cross(direction, e2)
    a  = float(np.dot(e1, h))
    if abs(a) < EPS:
        return None
    f = 1.0 / a
    s = origin - v0
    u = f * float(np.dot(s, h))
    if not (0.0 <= u <= 1.0):
        return None
    q = np.cross(s, e1)
    v = f * float(np.dot(direction, q))
    if v < 0.0 or u + v > 1.0:
        return None
    t = f * float(np.dot(e2, q))
    return t if t > EPS else None


def _ray_hits_panel(origin, direction, ray_len, panel) -> bool:
    """True if the ray origin→(origin + ray_len*direction) passes through panel."""
    c = panel.corners
    for tri in [(c[0], c[1], c[2]), (c[0], c[2], c[3])]:
        t = _ray_triangle(origin, direction, *tri)
        if t is not None and 0.0 < t < ray_len:
            return True
    return False


def compute_shielding(person_pos: np.ndarray, blast_pos: np.ndarray,
                      panels, panel_results: dict):
    """
    Returns (attenuation_factor, n_walls_between).
    attenuation_factor ∈ [0, 1]  —  1 = no shielding.
    """
    ray = person_pos - blast_pos
    ray_len = float(np.linalg.norm(ray))
    if ray_len < 0.01:
        return 1.0, 0

    ray_dir = ray / ray_len
    attenuation = 1.0
    n_walls = 0

    for panel in panels:
        if panel.panel_type in ('floor', 'roof'):
            continue   # ignore horizontal surfaces for horizontal blast path

        if not _ray_hits_panel(blast_pos, ray_dir, ray_len, panel):
            continue

        n_walls += 1
        pr = panel_results.get(panel.id)
        failed = (pr.failed if pr is not None else getattr(panel, 'failed', False))

        if failed:
            attenuation *= _FAILED_TRANSMISSION
        else:
            mat_name = getattr(panel.material, 'name', '')
            attenuation *= _TRANSMISSION.get(mat_name, 0.20)

    return float(attenuation), n_walls


# ---------------------------------------------------------------------------
# Injury assessment
# ---------------------------------------------------------------------------

def _probit_prob(Y: float) -> float:
    return float(np.clip(ndtr(Y - 5.0), 0.0, 1.0))


def assess_person_injuries(person: Person, blast: BlastSource,
                            panels, panel_results: dict) -> PersonInjury:
    """
    Full injury assessment for one person.

    Parameters
    ----------
    panels        : list of Panel (or duck-typed objects with .corners/.center/.normal/.failed/.material/.panel_type)
    panel_results : dict[panel_id -> result_object_with_.failed]
    """
    pos  = np.array([person.x, person.y, person.z], dtype=float)
    R    = max(float(np.linalg.norm(pos - blast.pos)), 0.1)
    W    = blast.W_eff
    Z    = float(np.clip(scaled_distance(R, W), 0.05, 40.0))

    raw_Pso = float(peak_overpressure(Z))
    Is_raw  = float(specific_impulse(Z, W))          # kPa·ms

    shielding, n_walls = compute_shielding(pos, blast.pos, panels, panel_results)
    eff_Pso = raw_Pso * shielding
    eff_Is  = Is_raw  * shielding

    injuries: List[InjuryDetail] = []
    Pso_Pa = eff_Pso * 1e3   # Pa for Probit formulas

    # ── Eardrum rupture (Baker 1983) ──────────────────────────────────────
    p_ear = _probit_prob(-15.6 + 1.93 * np.log(max(Pso_Pa, 1.0)))
    injuries.append(InjuryDetail(
        name='Eardrum Rupture',
        probability=p_ear,
        severity=('None' if p_ear < 0.05 else
                  'Mild' if p_ear < 0.50 else 'Severe'),
        mechanism='Overpressure differential ruptures tympanic membrane',
    ))

    # ── Pulmonary barotrauma / lung haemorrhage (Baker 1983) ──────────────
    Pr_lung     = -77.1 + 6.91 * np.log(max(Pso_Pa, 1.0))
    p_lung_f    = _probit_prob(Pr_lung)
    p_lung_s    = max(0.0, _probit_prob(Pr_lung + 2.0) - p_lung_f)
    p_lung_m    = max(0.0, _probit_prob(Pr_lung + 4.0) - p_lung_s - p_lung_f)
    lung_total  = p_lung_f + p_lung_s + p_lung_m
    lung_sev    = ('Fatal'    if p_lung_f > 0.10 else
                   'Severe'   if p_lung_s > 0.20 else
                   'Mild'     if p_lung_m > 0.10 else
                   'Low Risk' if eff_Pso > 40  else 'None')
    injuries.append(InjuryDetail(
        name='Pulmonary Barotrauma',
        probability=min(lung_total, 1.0),
        severity=lung_sev,
        mechanism='Rapid pressure wave causes alveolar rupture and haemorrhage',
    ))

    # ── Abdominal / GI barotrauma ─────────────────────────────────────────
    p_gi = float(np.clip((eff_Pso - 150.0) / 200.0, 0.0, 1.0)) if eff_Pso > 150 else 0.0
    injuries.append(InjuryDetail(
        name='Abdominal / GI Injury',
        probability=p_gi,
        severity=('None' if p_gi < 0.05 else
                  'Mild' if p_gi < 0.35 else 'Severe'),
        mechanism='Gas-filled organs (bowel, stomach) rupture from pressure wave',
    ))

    # ── Blast TBI ─────────────────────────────────────────────────────────
    p_tbi = float(np.clip((eff_Pso - 80.0) / 350.0, 0.0, 1.0)) if eff_Pso > 80 else 0.0
    injuries.append(InjuryDetail(
        name='Blast Traumatic Brain Injury',
        probability=p_tbi,
        severity=('None'     if p_tbi < 0.05 else
                  'Mild'     if p_tbi < 0.30 else
                  'Moderate' if p_tbi < 0.65 else 'Severe'),
        mechanism='Pressure wave coupling through skull causes concussive brain injury',
    ))

    # ── Tertiary: whole-body displacement ────────────────────────────────
    # Throw velocity: Is (Pa·s ≡ kPa·ms) / specific mass (kg/m²)
    v_throw = float(eff_Is * _FRONTAL_AREA / person.mass_kg)  # m/s
    p_tert  = float(np.clip((v_throw - 2.0) / 8.0, 0.0, 1.0)) if v_throw > 2 else 0.0
    tert_sev = ('None'     if v_throw < 2.0  else
                'Minor'    if v_throw < 5.0  else
                'Moderate' if v_throw < 10.0 else
                'Severe'   if v_throw < 16.0 else 'Fatal')
    injuries.append(InjuryDetail(
        name='Whole-Body Displacement',
        probability=p_tert,
        severity=tert_sev,
        mechanism=f'Impulse accelerates body to ~{v_throw:.1f} m/s; fractures / head impact on landing',
    ))

    # ── Secondary: debris / fragments ────────────────────────────────────
    nearest_failed = float('inf')
    n_failed_near  = 0
    pos_arr = np.asarray(pos)
    for panel in panels:
        pr = panel_results.get(panel.id)
        failed_p = (pr.failed if pr is not None else getattr(panel, 'failed', False))
        if failed_p:
            d = float(np.linalg.norm(np.asarray(panel.center) - pos_arr))
            nearest_failed = min(nearest_failed, d)
            if d < 5.0:
                n_failed_near += 1

    if nearest_failed < 1.0:
        deb_sev, p_deb = 'Severe',   0.90
        deb_mech = f'Structural failure {nearest_failed:.1f} m away — crush/impact likely'
    elif nearest_failed < 3.0:
        deb_sev, p_deb = 'Moderate', 0.55
        deb_mech = f'Failed panel at {nearest_failed:.1f} m — glass/fragment lacerations'
    elif nearest_failed < 6.0:
        deb_sev, p_deb = 'Mild',     0.25
        deb_mech = f'Failed panel at {nearest_failed:.1f} m — secondary debris possible'
    else:
        deb_sev, p_deb = 'None',     0.00
        deb_mech = 'No nearby panel failures'
    injuries.append(InjuryDetail(
        name='Fragment / Debris',
        probability=p_deb,
        severity=deb_sev,
        mechanism=deb_mech,
    ))

    # ── Overall severity ─────────────────────────────────────────────────
    _rank = {'None': 0, 'Low Risk': 0, 'Low': 0,
             'Mild': 1, 'Minor': 1,
             'Moderate': 2,
             'Severe': 3, 'Fatal': 4}
    max_rank = max(_rank.get(inj.severity, 0) for inj in injuries)
    overall  = ['Uninjured', 'Minor', 'Moderate', 'Severe', 'Fatal'][max_rank]

    p_fatal_total = max(p_lung_f,
                        1.0 if v_throw > 16 else p_tert * 0.5,
                        1.0 if nearest_failed < 0.5 else 0.0)
    survival = float(np.clip(1.0 - p_fatal_total, 0.0, 1.0))

    return PersonInjury(
        person_id=person.id,
        name=person.name,
        position=(person.x, person.y, person.z),
        standoff_m=R,
        raw_pressure_kPa=raw_Pso,
        effective_pressure_kPa=eff_Pso,
        impulse_kPa_ms=eff_Is,
        shielding_factor=shielding,
        walls_between=n_walls,
        injuries=injuries,
        throw_velocity_ms=v_throw,
        nearest_failed_m=nearest_failed,
        overall_severity=overall,
        survival_probability=survival,
    )
