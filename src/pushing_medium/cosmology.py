"""PM cosmology: refractive redshift, luminosity distance, structure growth.

Core hypothesis (from docs/notes/COSMOLOGY_FINDINGS.md):
  1 + z = n(t_emit) / n(t_obs)         [refractive redshift]
  H_eff(z) = H0 × (1+z)^β             [density-dependent decay rate]

β = 0 → "Empty Universe" (Milne): D_L = (c/H0) z  — too bright at high-z
β ≈ 0.8 → matches ΛCDM D_L to ~3% over z=0.1–1.5  [calibrated]

Structure growth in PM  (two scenarios):
─────────────────────────────────────────────────────────────────────────────
Scenario A — PM-Static (no kinematic friction from index decay):
  δ̈ = 4πGρ_m δ    →  δ ∝ cosh/sinh(Γt),  Γ = H0 √(3Ω_m/2)
  With ρ_m constant (no dilution), growth is EXPONENTIAL.
  This produces far too much structure → hard falsification.

Scenario B — PM-Drag (index-decay rate acts as Hubble friction):
  Same ODE as ΛCDM but with H replaced by H_PM(z) = H0(1+z)^β
  Growing mode in matter domination: D ∝ a^p where
    p² + (2−β)p − (3/2)Ω_m^eff = 0
  With β=0.8, Ω_m^eff=1: p ≈ 0.76  (vs p=1 in ΛCDM matter domination)
  → SLOWER growth rate, different fσ8(z) shape than ΛCDM
─────────────────────────────────────────────────────────────────────────────
"""

import math
import numpy as np
from scipy.integrate import solve_ivp, quad

# ── Constants ─────────────────────────────────────────────────────────────────
H0_PLANCK    = 67.4        # km/s/Mpc  (used as default)
OM0_DEFAULT  = 0.315       # same as ΛCDM for comparison
BETA_CALIB   = 0.8         # calibrated to match ΛCDM D_L within ~3%
SIGMA8_PLANCK = 0.811
C_KM_S       = 299792.458  # km/s

# ── Background ────────────────────────────────────────────────────────────────

def pm_hubble_effective(z, H0=H0_PLANCK, beta=BETA_CALIB):
    """H_PM(z) = H0 (1+z)^β [same units as H0].

    β = 0.0 → constant decay (empty/Milne universe, D_L linear)
    β ≈ 0.8 → calibrated to match ΛCDM luminosity distances at z ≤ 1.5
    """
    return H0 * (1.0 + z) ** beta


def pm_luminosity_distance(z, H0=H0_PLANCK, beta=BETA_CALIB):
    """Luminosity distance in PM refractive cosmology [Mpc].

    D_L(z) = (1+z) × c × ∫₀ᶻ dz' / H_PM(z')
    """
    def integrand(zp):
        return 1.0 / pm_hubble_effective(zp, H0, beta)
    comoving, _ = quad(integrand, 0.0, z, limit=500)
    return (1.0 + z) * C_KM_S * comoving

# ── Scenario A: static universe, no Hubble friction ──────────────────────────

def pm_static_growth_exponent(Om0=OM0_DEFAULT, H0=H0_PLANCK):
    """Growth exponent Γ / H0 = √(3 Ω_m / 2) for static PM universe.

    In the static case there is no kinematic friction:
        δ̈ = 4πG ρ_m δ  →  δ ∝ exp(Γ t),  Γ = H0 √(3Ω_m/2)

    fσ8 in this scenario is not well-defined in the usual RSD sense
    because there is no background expansion.  What can be quoted is the
    effective f = Γ / H_eff, which is constant with redshift.
    """
    return math.sqrt(1.5 * Om0)   # in units of H0


def pm_static_fsigma8(z_values, Om0=OM0_DEFAULT, sigma8_0=SIGMA8_PLANCK,
                      H0=H0_PLANCK, beta=BETA_CALIB):
    """fσ8 estimate for PM-Static scenario.

    In a static universe, the growth rate Γ = H0√(3Ω_m/2) is constant.
    The effective f (ratio to H_eff) is Γ/H_PM(z).
    σ8 is NOT suppressed by expansion (no dilution), so structure amplitude
    grows exponentially.  We use an exponential D(z) to parametrise this.

    CAUTION: This scenario produces catastrophic over-growth at early times
    and is principally used to demonstrate the falsification.
    """
    z_arr = np.atleast_1d(np.asarray(z_values, dtype=float))
    Gamma = pm_static_growth_exponent(Om0, H0) * H0  # in km/s/Mpc

    fs8 = np.empty(len(z_arr))
    for i, z in enumerate(z_arr):
        H_eff  = pm_hubble_effective(z, H0, beta)
        f_eff  = Gamma / H_eff   # dimensionless effective growth rate
        # σ8(z) = σ8,0 (no suppression — structure accumulates in static BG)
        fs8[i] = f_eff * sigma8_0
    return fs8

# ── Scenario B: PM-Drag growth ODE ───────────────────────────────────────────

def _pm_drag_growth_rhs(a, y, Om_eff, beta):
    """RHS of PM-Drag growth ODE in terms of a.

    H_PM(a) = H0 a^{-β}  →  d ln H / d ln a = -β

    ODE:  D'' + [(3−β)/a] D' − [3 Ω_eff / (2 a^{5−2β})] D = 0

    Y.C.:  D = a^p_grow,  dD/da = p_grow a^{p_grow-1}  (matter-dominated)
    where p_grow is the positive root of p² + (2−β)p − (3/2)Ω_eff = 0
    """
    D, Dp = y
    p_coeff = (3.0 - beta) / a
    q_coeff = 1.5 * Om_eff / a**(5.0 - 2.0*beta)
    return [Dp, -p_coeff * Dp + q_coeff * D]


def pm_drag_growing_power(Om_eff, beta):
    """Positive root of p² + (2−β)p − (3/2) Ω_eff = 0 (matter-dom IC)."""
    b = 2.0 - beta
    c = -1.5 * Om_eff
    return (-b + math.sqrt(b*b - 4.0*c)) / 2.0


def pm_growth_factor_rate(z_values, Om_eff=OM0_DEFAULT, beta=BETA_CALIB):
    """Solve PM-Drag growth ODE; return (D(z)/D(0), f(z)) arrays.

    Parameters
    ----------
    z_values : array-like
    Om_eff   : float   effective matter density parameter in PM (often 0.315 or 1.0)
    beta     : float   PM Hubble exponent (default 0.8, calibrated to D_L)

    Returns
    -------
    D_ratio : ndarray   D(z)/D(0)
    f_arr   : ndarray   growth rate f(z) = d ln D / d ln a
    """
    z_arr = np.atleast_1d(np.asarray(z_values, dtype=float))
    a_arr = np.clip(1.0 / (1.0 + z_arr), 1e-4, 1.0 - 1e-12)

    a_ini  = 1e-4
    p_grow = pm_drag_growing_power(Om_eff, beta)
    y0     = [a_ini**p_grow, p_grow * a_ini**(p_grow - 1.0)]   # MD growing mode

    a_eval = np.unique(np.concatenate([[a_ini], np.sort(a_arr), [1.0]]))

    sol = solve_ivp(
        _pm_drag_growth_rhs,
        [a_ini, 1.0],
        y0,
        args=(Om_eff, beta),
        method='DOP853',
        t_eval=a_eval,
        dense_output=True,
        rtol=1e-9, atol=1e-11,
    )
    if not sol.success:
        raise RuntimeError(f"PM-Drag growth ODE failed: {sol.message}")

    D0 = sol.sol(1.0)[0]

    D_ratio = np.empty(len(z_arr))
    f_arr   = np.empty(len(z_arr))
    for i, a in enumerate(a_arr):
        y    = sol.sol(a)
        D_z  = y[0]
        Dp_z = y[1]
        D_ratio[i] = D_z / D0
        f_arr[i]   = a * Dp_z / D_z

    return D_ratio, f_arr


def pm_drag_fsigma8(z_values, Om_eff=OM0_DEFAULT, beta=BETA_CALIB,
                    sigma8_0=SIGMA8_PLANCK):
    """fσ8(z) for PM-Drag (refractive Hubble friction) scenario.

    Parameters
    ----------
    z_values  : float or array-like
    Om_eff    : float   effective Ω_m (0.315 = same as ΛCDM; 1.0 = matter only)
    beta      : float   PM Hubble exponent (0.8 calibrated)
    sigma8_0  : float   amplitude today

    Returns
    -------
    fsig8 : ndarray
    """
    D_ratio, f_arr = pm_growth_factor_rate(z_values, Om_eff, beta)
    return f_arr * sigma8_0 * D_ratio
