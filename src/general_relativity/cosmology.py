"""ΛCDM cosmological functions: Hubble function, luminosity distance, growth factor, fσ8.

All distances in Mpc.  H0 in km/s/Mpc.  σ8 dimensionless.

Growth ODE (flat ΛCDM, d/da form):
    D'' + [(3 + d ln H / d ln a) / a] D' - [3 Ω_m0 / (2 a^5 (H/H0)^2)] D = 0
    growing-mode IC: D = a, dD/da = 1  at  a_ini = 1e-4 (deep matter domination)

Growth rate:  f(z) = d ln D / d ln a
fσ8 observable: fσ8(z) = f(z) × σ8,0 × D(z)/D(0)
"""

import math
import numpy as np
from scipy.integrate import solve_ivp, quad

# ── Planck 2018 defaults ──────────────────────────────────────────────────────
H0_PLANCK    = 67.4    # km/s/Mpc
OM0_PLANCK   = 0.315
OL0_PLANCK   = 0.685
SIGMA8_PLANCK = 0.811
C_KM_S       = 299792.458   # km/s

# ── Background ────────────────────────────────────────────────────────────────

def lcdm_hubble(z, H0=H0_PLANCK, Om0=OM0_PLANCK, Ol0=OL0_PLANCK):
    """H(z) = H0 sqrt(Ω_m (1+z)³ + Ω_Λ) [same units as H0]."""
    return H0 * math.sqrt(Om0 * (1.0 + z)**3 + Ol0)


def lcdm_luminosity_distance(z, H0=H0_PLANCK, Om0=OM0_PLANCK, Ol0=OL0_PLANCK):
    """Luminosity distance D_L(z) [Mpc], flat universe."""
    def integrand(zp):
        return 1.0 / lcdm_hubble(zp, H0, Om0, Ol0)
    comoving, _ = quad(integrand, 0.0, z, limit=500)
    return (1.0 + z) * C_KM_S * comoving


def lcdm_om_z(z, Om0=OM0_PLANCK, Ol0=OL0_PLANCK):
    """Effective Ω_m(z) = Ω_m0 (1+z)^3 / E^2(z)."""
    E2 = Om0 * (1.0 + z)**3 + Ol0
    return Om0 * (1.0 + z)**3 / E2

# ── Growth ODE ────────────────────────────────────────────────────────────────

def _lcdm_growth_rhs(a, y, Om0):
    """RHS of growth ODE in terms of scale factor a.

    State vector y = [D, D'] where ' = d/da.
    ODE:  D'' + [(3 + d ln H / d ln a) / a] D' − [3 Ω_m0 / (2 a^5 E^2)] D = 0
    """
    D, Dp = y
    Ol0 = 1.0 - Om0
    E2      = Om0 / a**3 + Ol0           # (H/H0)²
    dE2da   = -3.0 * Om0 / a**4          # d(E²)/da
    dlnH_dlna = a * dE2da / (2.0 * E2)   # d ln H / d ln a

    p_coeff = (3.0 + dlnH_dlna) / a
    q_coeff =  1.5 * Om0 / (a**5 * E2)

    return [Dp, -p_coeff * Dp + q_coeff * D]


def lcdm_growth_factor_rate(z_values, Om0=OM0_PLANCK):
    """Solve ΛCDM growth ODE; return (D(z)/D(0), f(z)) arrays.

    f(z) = d ln D / d ln a = (a / D) × dD/da  (computed from ODE solution).

    Parameters
    ----------
    z_values : array-like   redshifts at which to evaluate.
    Om0      : float        present-day matter density parameter.

    Returns
    -------
    D_ratio : ndarray   D(z)/D(0)
    f_arr   : ndarray   growth rate f(z)
    """
    z_arr = np.atleast_1d(np.asarray(z_values, dtype=float))
    a_arr = np.clip(1.0 / (1.0 + z_arr), 1e-4, 1.0 - 1e-12)

    a_ini = 1e-4
    a_fin = 1.0
    y0    = [a_ini, 1.0]          # D=a, dD/da=1 (MD growing mode)

    a_eval = np.unique(np.concatenate([[a_ini], np.sort(a_arr), [a_fin]]))

    sol = solve_ivp(
        _lcdm_growth_rhs,
        [a_ini, a_fin],
        y0,
        args=(Om0,),
        method='DOP853',
        t_eval=a_eval,
        dense_output=True,
        rtol=1e-9, atol=1e-11,
    )
    if not sol.success:
        raise RuntimeError(f"ΛCDM growth ODE failed: {sol.message}")

    D0  = sol.sol(1.0)[0]    # D at a = 1

    D_ratio = np.empty(len(z_arr))
    f_arr   = np.empty(len(z_arr))
    for i, a in enumerate(a_arr):
        y    = sol.sol(a)
        D_z  = y[0]
        Dp_z = y[1]           # dD/da
        D_ratio[i] = D_z / D0
        f_arr[i]   = a * Dp_z / D_z    # d ln D / d ln a

    return D_ratio, f_arr


def lcdm_fsigma8(z_values, Om0=OM0_PLANCK, sigma8_0=SIGMA8_PLANCK):
    """fσ8(z) = f(z) × σ8,0 × D(z)/D(0) using numerically solved growth ODE.

    Parameters
    ----------
    z_values  : float or array-like
    Om0       : float   Ω_m0
    sigma8_0  : float   σ8 at z = 0

    Returns
    -------
    fsig8 : ndarray  same shape as z_values
    """
    D_ratio, f_arr = lcdm_growth_factor_rate(z_values, Om0)
    return f_arr * sigma8_0 * D_ratio


def lcdm_growth_rate_approx(z, Om0=OM0_PLANCK, gamma=0.55):
    """Linder (2005) approximation: f(z) ≈ Ω_m(z)^γ, γ = 6/11 ≈ 0.55."""
    return lcdm_om_z(z, Om0, 1.0 - Om0) ** gamma
