"""
Blast loading model based on Kinney & Graham (1985) empirical relations
and the Friedlander modified exponential waveform.

All pressures in kPa, times in seconds, distances in metres,
TNT equivalent mass in kg.
"""

import numpy as np
from scipy.optimize import brentq

P_ATM = 101.325   # kPa  atmospheric pressure
C_SOUND = 340.0   # m/s  ambient speed of sound


# ---------------------------------------------------------------------------
# Scaled-distance relations (Kinney & Graham 1985)
# ---------------------------------------------------------------------------

def scaled_distance(standoff_m: float, tnt_kg: float) -> float:
    return standoff_m / max(tnt_kg, 1e-6) ** (1.0 / 3.0)


def peak_overpressure(Z: float) -> float:
    """Free-air side-on peak overpressure [kPa].  Z in m/kg^(1/3)."""
    Z = np.maximum(Z, 0.05)
    num = 808.0 * (1.0 + (Z / 4.5) ** 2)
    den = (np.sqrt(1.0 + (Z / 0.048) ** 2) *
           np.sqrt(1.0 + (Z / 0.32) ** 2) *
           np.sqrt(1.0 + (Z / 1.35) ** 2))
    return P_ATM * num / den


def positive_duration(Z: float, W: float) -> float:
    """Positive phase duration [ms].  Z in m/kg^(1/3), W in kg."""
    Z = np.maximum(Z, 0.05)
    num = 980.0 * (1.0 + (Z / 0.54) ** 10)
    den = ((1.0 + (Z / 0.02) ** 3) *
           (1.0 + (Z / 0.74) ** 6) *
           np.sqrt(1.0 + (Z / 6.9) ** 2))
    return (num / den) * W ** (1.0 / 3.0)


def specific_impulse(Z: float, W: float) -> float:
    """Specific positive impulse [kPa·ms].  Z in m/kg^(1/3), W in kg."""
    Z = np.maximum(Z, 0.05)
    scaled = P_ATM * 0.067 * np.sqrt(1.0 + (Z / 0.23) ** 4) / (
        Z ** 2 * np.sqrt(1.0 + (Z / 1.55) ** 3))
    return scaled * W ** (1.0 / 3.0)


def _solve_waveform_b(Pso: float, td: float, Is: float) -> float:
    """Friedlander waveform parameter b from integral of impulse."""
    def residual(b):
        if b < 1e-9:
            return Pso * td / 2.0 - Is
        return Pso * td * (1.0 / b - (1.0 - np.exp(-b)) / b ** 2) - Is
    try:
        return brentq(residual, 0.05, 30.0, xtol=1e-6)
    except ValueError:
        return 1.0


def friedlander(t_ms: np.ndarray, Pso: float, td_ms: float, b: float) -> np.ndarray:
    """Modified Friedlander waveform.  t_ms: shifted time [ms], returns kPa."""
    t = np.asarray(t_ms)
    return np.where(
        (t >= 0.0) & (t <= td_ms),
        Pso * (1.0 - t / td_ms) * np.exp(-b * t / td_ms),
        0.0,
    )


def rankine_hugoniot_reflected(Pso: float) -> float:
    """Peak normally reflected overpressure [kPa] via Rankine-Hugoniot."""
    return 2.0 * Pso * (7.0 * P_ATM + 4.0 * Pso) / (7.0 * P_ATM + Pso)


# ---------------------------------------------------------------------------
# Blast source
# ---------------------------------------------------------------------------

class BlastSource:
    """
    Point blast source using Kinney-Graham scaling and Friedlander waveform.

    burst_type:
        'surface'  – hemispherical (gas leak indoor/outdoor ground burst);
                     effective TNT doubled to account for ground reflection.
        'free_air' – spherical free-air burst.
    """

    def __init__(self, x: float, y: float, z: float,
                 tnt_kg: float, burst_type: str = 'surface'):
        self.pos = np.array([x, y, z], dtype=float)
        self.W_input = tnt_kg
        # Surface burst: energy confined to hemisphere → equivalent to 2× free-air
        self.W_eff = tnt_kg * 2.0 if burst_type == 'surface' else tnt_kg
        self.burst_type = burst_type

    def panel_loading(self,
                      center: np.ndarray,
                      normal: np.ndarray,
                      t_arr: np.ndarray):
        """
        Compute reflected pressure time-history on a panel.

        Parameters
        ----------
        center : array-like (3,)  – panel centroid [m]
        normal : array-like (3,)  – outward unit normal of the panel
        t_arr  : ndarray          – time array [s]

        Returns
        -------
        pressure : ndarray  [kPa]  same shape as t_arr
        meta     : dict     diagnostic scalars
        """
        r = np.asarray(center, dtype=float) - self.pos
        R = max(float(np.linalg.norm(r)), 0.10)

        Z = np.clip(scaled_distance(R, self.W_eff), 0.05, 40.0)

        Pso = float(peak_overpressure(Z))
        td  = float(positive_duration(Z, self.W_eff))   # ms
        Is  = float(specific_impulse(Z, self.W_eff))    # kPa·ms
        b   = _solve_waveform_b(Pso, td, Is)

        # Angle of incidence: 0 = normal incidence (max reflection),
        #                     90 = grazing (no enhancement)
        n = np.asarray(normal, dtype=float)
        n = n / max(np.linalg.norm(n), 1e-12)
        r_unit = r / R
        cos_alpha = float(np.clip(abs(np.dot(-r_unit, n)), 0.0, 1.0))

        Pr_normal = rankine_hugoniot_reflected(Pso)
        # Linear interpolation: normal → side-on as cos_alpha goes 1 → 0
        Pr = Pso + (Pr_normal - Pso) * cos_alpha

        # Time of arrival (crude: ambient speed of sound)
        toa_s = R / C_SOUND

        t_ms = (np.asarray(t_arr, dtype=float) - toa_s) * 1000.0
        pressure = friedlander(t_ms, Pr, td, b)

        alpha_deg = float(np.degrees(np.arccos(cos_alpha)))
        meta = {
            'Pso': Pso, 'Pr': Pr, 'td_ms': td, 'Is': Is, 'b': b,
            'Z': Z, 'R': R, 'alpha_deg': alpha_deg, 'toa_s': toa_s,
        }
        return pressure, meta
