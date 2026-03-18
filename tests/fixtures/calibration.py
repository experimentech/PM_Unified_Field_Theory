import json
import os
from dataclasses import dataclass

from pushing_medium.core import (
    index_deflection_numeric,
    fermat_deflection_static_index,
    pm_deflection_angle_point_mass,
    pm_binary_quadrupole_power,
)

CAL_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'calibration.json')


@dataclass
class Calibration:
    mu_coeff: float
    k_TT: float
    k_Fizeau: float


def _fit_mu_for_deflection(M=1.989e30, b=6.96e8, z_max=2e10, steps=2500) -> float:
    target = pm_deflection_angle_point_mass(M, b)
    approx_mu = 2 * 6.67430e-11 / (299792458.0 ** 2)
    lo, hi = 0.25 * approx_mu, 4.0 * approx_mu
    for _ in range(36):
        mid = 0.5 * (lo + hi)
        val = index_deflection_numeric(M, b, mu=mid, z_max=z_max, steps=steps)
        if val < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def load_or_calibrate() -> Calibration:
    data = {}
    if os.path.exists(CAL_PATH):
        try:
            data = json.load(open(CAL_PATH, 'r'))
        except Exception:
            data = {}

    mu = data.get('mu_coeff')
    if mu is None:
        mu = _fit_mu_for_deflection()
        data['mu_coeff'] = mu

    # TT normalization: PM formula already matches GR Peters form → ratio = 1
    k_TT = data.get('k_TT', 1.0)

    # Fizeau coupling: match GR first-order (1 + v/c) correction
    k_Fizeau = data.get('k_Fizeau')
    if k_Fizeau is None:
        M = 1.989e30
        b = 6.96e8
        mu = mu or _fit_mu_for_deflection()
        alpha_static = fermat_deflection_static_index(M, b, mu=mu, z_max=2e10, steps=1500)
        v = 3.0e4  # 30 km/s
        target_ratio = 1.0 + v / 299792458.0
        k_Fizeau = (target_ratio - 1.0) / (v / 299792458.0)
        data['k_Fizeau'] = k_Fizeau

    data['k_TT'] = k_TT
    with open(CAL_PATH, 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)

    return Calibration(mu_coeff=mu, k_TT=k_TT, k_Fizeau=k_Fizeau)
