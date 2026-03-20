"""Option-A Corrected Poisson Equation — Numerical Experiment.

Tests whether promoting U(φ) into the PM Lagrangian as a potential term
-V(φ) = -U(φ)/c_φ² changes the M-R curve in a way that addresses the
strong-field gaps identified in the formula sheet.

Three scenarios are compared:
  baseline  (α = 0)  : exact current PM, ∇²φ = -(8πG/c²)ρ
  Gap-1     (α = +1) : U' reduces the effective source (formula-sheet direction)
  Gap-3     (α = -1) : U' adds to the gravitating source (field energy gravitates)

Output: tabulated R_1.4, M_max, surface redshift z_s, compactness C at M_max,
        and photon sphere ratio r_ps/R for each scenario.  Also prints the
        φ(r) profile for a 1.4 M_sun star to show how the field is modified.

Usage:
    python scripts/option_a_corrected_poisson.py

Results are used to update pm-formula-sheet.tex §Option-A Corrected Poisson.
"""

import sys
import math
import numpy as np

sys.path.insert(0, 'src')

from pushing_medium.stellar_structure import (
    solve_pm_star,
    solve_pm_star_option_a,
    compute_mr_curve,
    compute_mr_curve_option_a,
    pm_surface_redshift,
    RHO_NUC, RHO_CRIT, M_SUN, G, c, MU_G,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_r14(M_arr, R_arr):
    """Interpolate R at M = 1.4 M_sun from the M-R array."""
    valid = np.isfinite(M_arr) & np.isfinite(R_arr)
    M_v = M_arr[valid]
    R_v = R_arr[valid]
    # Sort by mass
    idx = np.argsort(M_v)
    M_v, R_v = M_v[idx], R_v[idx]
    if M_v[0] > 1.4 or M_v[-1] < 1.4:
        return float('nan')
    return float(np.interp(1.4, M_v, R_v))


def compactness(M_msun, R_km):
    return G * M_msun * M_SUN / (c * c * R_km * 1e3)


def photon_sphere_ratio(M_msun, R_km):
    """r_ps/R_star where r_ps = 2GM/c² (PM photon sphere)."""
    r_ps = MU_G * M_msun * M_SUN   # = 2GM/c²
    return r_ps / (R_km * 1e3)


# ---------------------------------------------------------------------------
# Compute M-R curves for three scenarios
# ---------------------------------------------------------------------------
N = 50   # points per M-R curve (moderate for speed)

print("Computing M-R curves … (3 × ~50 ODE solutions)")
print()

_, M_base, R_base = compute_mr_curve(n_points=N)
_, M_gp1, R_gp1   = compute_mr_curve_option_a(alpha=+1.0, n_points=N)
_, M_gp3, R_gp3   = compute_mr_curve_option_a(alpha=-1.0, n_points=N)

# ---------------------------------------------------------------------------
# Key summary quantities
# ---------------------------------------------------------------------------

def summarise(name, M_arr, R_arr):
    valid = np.isfinite(M_arr) & np.isfinite(R_arr)
    if not valid.any():
        print(f"  {name}: NO CONVERGED MODELS"); return

    M_v   = M_arr[valid]
    R_v   = R_arr[valid]
    M_max = float(np.nanmax(M_v))

    # R at M_max
    idx_max = int(np.nanargmax(M_v))
    R_at_Mmax = R_v[idx_max]

    # R at 1.4 M_sun
    R_1p4 = find_r14(M_arr, R_arr)

    # Compactness and photon sphere at M_max
    C_max    = compactness(M_max, R_at_Mmax)
    rps_R    = photon_sphere_ratio(M_max, R_at_Mmax)

    # Surface redshift at 1.4 M_sun
    if not math.isnan(R_1p4):
        z_s14 = pm_surface_redshift(1.4 * M_SUN, R_1p4 * 1e3)
    else:
        z_s14 = float('nan')

    # Surface redshift at M_max
    z_s_max = pm_surface_redshift(M_max * M_SUN, R_at_Mmax * 1e3)

    print(f"  Scenario: {name}")
    print(f"    M_max           = {M_max:.2f} M_sun")
    print(f"    R at M_max      = {R_at_Mmax:.2f} km")
    print(f"    C at M_max      = {C_max:.4f}  (Buchdahl limit 4/9 ≈ 0.444)")
    print(f"    r_ps/R at M_max = {rps_R:.4f}  (>1 → photon sphere outside star)")
    print(f"    z_s at M_max    = {z_s_max:.4f}")
    print(f"    R_1.4           = {R_1p4:.2f} km  (GW170817 limit ≤ 13.3 km)")
    print(f"    z_s at 1.4 M☉   = {z_s14:.4f}")
    print()

print("=" * 60)
print("M-R CURVE COMPARISON: Option-A Corrected Poisson Equation")
print("=" * 60)
print()
summarise("Baseline  (α =  0.0)", M_base, R_base)
summarise("Gap-1     (α = +1.0)", M_gp1,  R_gp1)
summarise("Gap-3     (α = -1.0)", M_gp3,  R_gp3)

# ---------------------------------------------------------------------------
# φ(r) profile comparison for a star near 1.4 M_sun
# ---------------------------------------------------------------------------
print("=" * 60)
print("φ(r) PROFILE — star near 1.4 M_sun (ρ_c = 1.4 × ρ_nuc)")
print("=" * 60)

rho_test = 1.4 * RHO_NUC

star_base = solve_pm_star(rho_test)
star_gp1  = solve_pm_star_option_a(rho_test, alpha=+1.0)
star_gp3  = solve_pm_star_option_a(rho_test, alpha=-1.0)

print(f"\n  Baseline: M={star_base.M_star/M_SUN:.3f} M_sun, "
      f"R={star_base.R_star/1e3:.2f} km, φ_c={star_base.phi_central:.4f}")
print(f"  Gap-1:    M={star_gp1.M_star/M_SUN:.3f} M_sun,  "
      f"R={star_gp1.R_star/1e3:.2f} km, φ_c={star_gp1.phi_central:.4f}")
print(f"  Gap-3:    M={star_gp3.M_star/M_SUN:.3f} M_sun,  "
      f"R={star_gp3.R_star/1e3:.2f} km, φ_c={star_gp3.phi_central:.4f}")

# Print φ at several radial fractions for the comparison star
print()
print("  φ(r) at fractional radii (baseline vs Gap-1 vs Gap-3):")
print(f"  {'r/R':>8}   {'φ_base':>9}   {'φ_Gap1':>9}   {'φ_Gap3':>9}")
R_base_m = star_base.R_star
R_gp1_m  = star_gp1.R_star
R_gp3_m  = star_gp3.R_star

for frac in [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]:
    # Interpolate φ at r = frac × R for each scenario
    r_b = frac * R_base_m; phi_b = float(np.interp(max(r_b,1.0), star_base.r, star_base.phi))
    r_g = frac * R_gp1_m;  phi_g = float(np.interp(max(r_g,1.0), star_gp1.r,  star_gp1.phi))
    r_3 = frac * R_gp3_m;  phi_3 = float(np.interp(max(r_3,1.0), star_gp3.r,  star_gp3.phi))
    print(f"  {frac:>8.1f}   {phi_b:>9.5f}   {phi_g:>9.5f}   {phi_3:>9.5f}")

# ---------------------------------------------------------------------------
# Observational comparison table
# ---------------------------------------------------------------------------
print()
print("=" * 60)
print("OBSERVATIONAL BOUNDS COMPARISON")
print("=" * 60)
print(f"  GW170817 radius bound: R_1.4 ≤ 13.3 km (90% CI)")
print(f"  PM baseline R_1.4    : {find_r14(M_base, R_base):.2f} km")
print(f"  Gap-1     R_1.4      : {find_r14(M_gp1,  R_gp1):.2f} km")
print(f"  Gap-3     R_1.4      : {find_r14(M_gp3,  R_gp3):.2f} km")
print()
print(f"  PM M_max baseline: {float(np.nanmax(M_base[np.isfinite(M_base)])):.2f} M_sun")
print(f"  PM M_max Gap-1:    {float(np.nanmax(M_gp1[np.isfinite(M_gp1)])):.2f} M_sun")
print(f"  PM M_max Gap-3:    {float(np.nanmax(M_gp3[np.isfinite(M_gp3)])):.2f} M_sun")
print()
print("  GWTC mass-gap context: nearest LIGO heavy object = 22.2 M_sun (GW190814 M1)")
print("  Compact objects above PM M_max cannot be formed under baseline EOS.")
print()

# ---------------------------------------------------------------------------
# Summary verdict
# ---------------------------------------------------------------------------
R_obs_limit = 13.3   # km
R_base_14 = find_r14(M_base, R_base)
R_gp1_14  = find_r14(M_gp1,  R_gp1)
R_gp3_14  = find_r14(M_gp3,  R_gp3)

print("=" * 60)
print("VERDICT")
print("=" * 60)

def verdict(name, R_14, M_max):
    r_ok  = not math.isnan(R_14) and R_14 <= R_obs_limit
    m_gap = M_max < 22.2
    print(f"  {name}:")
    print(f"    R_1.4 ≤ 13.3 km? {'YES ✓' if r_ok else f'NO  ({R_14:.2f} km)'}")
    print(f"    Mass gap kept?   {'YES (M_max < 22.2)' if m_gap else f'NO  (M_max={M_max:.1f})'}")

verdict("Baseline (α=0.0)", R_base_14, float(np.nanmax(M_base[np.isfinite(M_base)])))
verdict("Gap-1    (α=+1.0)", R_gp1_14, float(np.nanmax(M_gp1[np.isfinite(M_gp1)])))
verdict("Gap-3    (α=-1.0)", R_gp3_14, float(np.nanmax(M_gp3[np.isfinite(M_gp3)])))
