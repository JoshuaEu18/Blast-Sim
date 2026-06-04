"""
Blast casualty model based on peak overpressure.

Uses the Baker (1983) Probit function for lung-hemorrhage fatality
and threshold-based fractions for severe / minor injury.

Reference: Baker et al. (1983) Explosion Hazards and Evaluation.
"""

import numpy as np
from scipy.special import ndtr   # standard normal CDF


def _probit_fatality(Pr_kPa: float) -> float:
    """
    Probability of fatality from direct blast overpressure (lung hemorrhage).
    Baker (1983): Pr = -77.1 + 6.91·ln(Pso)  where Pso is in Pa.
    P(fatality) = Φ(Y − 5) where Φ is the standard normal CDF.
    """
    if Pr_kPa <= 0.0:
        return 0.0
    Pso_Pa = Pr_kPa * 1e3
    Y = -77.1 + 6.91 * np.log(Pso_Pa)
    return float(ndtr(Y - 5.0))


def _prob_severe(Pr_kPa: float) -> float:
    """
    Probability of severe injury (significant lung / ear damage).
    Threshold at ~35 kPa, 100% at ~200 kPa.
    """
    if Pr_kPa < 35.0:
        return 0.0
    if Pr_kPa >= 200.0:
        return 1.0
    return (Pr_kPa - 35.0) / 165.0


def _prob_minor(Pr_kPa: float) -> float:
    """
    Probability of minor injury (eardrum rupture, lacerations from glass).
    Threshold at ~7 kPa (1 psi), 100% at ~70 kPa.
    """
    if Pr_kPa < 7.0:
        return 0.0
    if Pr_kPa >= 70.0:
        return 1.0
    return (Pr_kPa - 7.0) / 63.0


def compute_room_casualties(peak_pressure_kPa: float, n_occupants: int) -> dict:
    """
    Returns expected number of casualties in a room.

    Categories are mutually exclusive from worst to least severe:
        fatal  : expected fatalities
        severe : severe injuries (not fatal)
        minor  : minor injuries (not severe or fatal)
    """
    p_fatal  = _probit_fatality(peak_pressure_kPa)
    p_severe = max(0.0, _prob_severe(peak_pressure_kPa) - p_fatal)
    p_minor  = max(0.0, _prob_minor(peak_pressure_kPa) - p_severe - p_fatal)

    p_fatal  = min(p_fatal,  1.0)
    p_severe = min(p_severe, 1.0 - p_fatal)
    p_minor  = min(p_minor,  1.0 - p_fatal - p_severe)

    return {
        'fatal':  p_fatal  * n_occupants,
        'severe': p_severe * n_occupants,
        'minor':  p_minor  * n_occupants,
    }
