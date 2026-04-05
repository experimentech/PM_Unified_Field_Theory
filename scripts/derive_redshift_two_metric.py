"""Two-metric derivation of PM surface gravitational redshift.

Key result
----------
The factor-of-2 discrepancy between the old PM formula (z = e^φ − 1) and
GR (z = 1/sqrt(1−2GM/c²R) − 1) is not a genuine physical prediction of PM.
It arises from applying the *optical* metric — appropriate for photon clocks —
to the emission by a *massive atomic* oscillator, which obeys the *particle*
metric.

PM has two distinct effective metrics
--------------------------------------
  Optical metric:   g_tt^opt  = −c²/n² = −c²e^{−2φ}    (governs photons)
  Particle metric:  g_tt^mass = −(1 − φ)                (governs massive matter)

(see §Effective Metric in pm-formula-sheet.tex)

Old (wrong-for-atoms) derivation — optical-metric clock ratio
--------------------------------------------------------------
  1 + z_old = sqrt(−g_tt^opt(∞) / −g_tt^opt(R))
            = sqrt(c² / (c²/n²(R)))
            = n(R) = e^{φ_s}

  z_old = e^{μ_G M/R} − 1             ← the previously boxed PM formula

At leading order: z_old ≈ 2GM/c²R ≈ 2 × z_GR.

This formula correctly gives the clock-rate ratio for a *photon* clock (a
resonant cavity whose tick rate is set by the round-trip light time).
But X-ray spectral lines from neutron stars are emitted by atomic transitions —
electrons in massive atoms — whose frequencies are governed by the particle
metric, not the optical metric.

Corrected derivation — particle-metric atom clock + conserved coord-ω
----------------------------------------------------------------------
Step 1 — Coordinate frequency is conserved in a static PM medium.
  The optical metric g_tt^opt = −c²/n²(r) is static (no t-dependence).
  The Killing vector ξ = ∂_t exists.  For a photon with 4-wavevector k^μ:
    E = −k_t = −g_tt^opt k^t = (c²/n²)(dt/dλ) = const along ray.
  For a monochromatic wave Ψ = A(r) e^{−iωt}: k_t = ∂_t(phase) = −ω.
  Therefore E = −k_t = ω = const.  Coordinate frequency is conserved. ✓

  (Alternatively, from the wave equation n²(r)/c² ∂_t²Ψ = ∇²Ψ: any
  monochromatic solution has time-dependence e^{−iωt} with fixed ω regardless
  of the spatial n(r) profile, confirming ω_coord is conserved.)

Step 2 — Atom proper time from particle metric.
  The emitting atom is at rest at r = R (stellar surface).
  Its proper time element from the particle metric g_tt^mass = −(1 − φ):
    dτ_atom = sqrt(1 − φ_s) dt,   φ_s = μ_G M/R = 2GM/c²R.

Step 3 — Convert proper emission frequency to coordinate frequency.
  The atom's transition has proper frequency ω₀ (in its rest frame,
  measured in proper time τ).  In coordinate time:
    ω_coord = ω₀ × (dτ_atom/dt) = ω₀ sqrt(1 − φ_s).

Step 4 — Receiver at infinity.
  At r → ∞: φ = 0, so dτ_recv = dt.
  The receiver's proper frequency equals the coordinate frequency:
    ω_recv = ω_coord = ω₀ sqrt(1 − φ_s).

Step 5 — Redshift.
    1 + z = ω₀ / ω_recv = 1 / sqrt(1 − φ_s)

  With φ_s = μ_G M/R = 2GM/c²R:
    z_PM_corrected = 1/sqrt(1 − 2GM/c²R) − 1  =  z_GR            ★

The corrected PM surface redshift is EXACTLY the GR Schwarzschild formula.

Why μ_G = 2G/c² is preserved (no conflict with solar-system tests)
-------------------------------------------------------------------
The factor-of-2 in μ_G = 2G/c² (relative to the Newtonian Φ_N = GM/r):
  - Correctly explains 4GM/c²b light deflection ✓  (optical metric, photon)
  - Correctly explains 2GM/c³ ln(…) Shapiro delay ✓ (optical metric, photon)
  - Was INCORRECTLY applied to atomic spectral-line redshift ✗ (now fixed)

The fix uses the particle metric g_tt^mass = −(1−φ) = −(1−2U) for the
emitting atom's clock rate.  Since φ_s = 2U, the particle metric gives
g_tt^mass = −(1−2U), numerically identical to the Schwarzschild g_tt at this
order.  All solar-system predictions remain unchanged.

In GR (Schwarzschild), there is only ONE metric: g_tt = −(1−2GM/c²R).
The distinction between an "atom clock" and a "photon clock" does not arise
because both use the same g_tt.  PM has two metrics, so the distinction matters.
"""

import math
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pushing_medium.stellar_structure import (
    MU_G, G, c, M_SUN,
    pm_surface_redshift,
    pm_surface_redshift_corrected,
)
from general_relativity.classical import gr_surface_redshift


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _compare(M_Msun: float, R_km: float) -> dict:
    """Return a row comparing old PM, corrected PM, and GR redshift."""
    M = M_Msun * M_SUN
    R = R_km * 1e3
    phi_s = MU_G * M / R

    z_old  = pm_surface_redshift(M, R)
    z_corr = pm_surface_redshift_corrected(M, R)
    z_gr   = gr_surface_redshift(M, R)

    ratio_old  = z_old  / z_gr if z_gr > 0 else float('nan')
    ratio_corr = z_corr / z_gr if (z_gr > 0 and math.isfinite(z_corr)) else float('nan')
    equals_gr  = math.isfinite(z_corr) and abs(z_corr - z_gr) < 1e-12 * (1 + abs(z_gr))

    return {
        'M': M_Msun, 'R': R_km, 'phi_s': phi_s,
        'z_old': z_old, 'z_corr': z_corr, 'z_gr': z_gr,
        'ratio_old': ratio_old, 'ratio_corr': ratio_corr,
        'equals_gr': equals_gr,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 76)
    print("Two-Metric Derivation — PM Surface Gravitational Redshift")
    print("=" * 76)
    print()

    # ── 1. Key analytic result ────────────────────────────────────────────────
    print("KEY RESULT")
    print("----------")
    print("Old formula  (optical-metric clock ratio, not applicable to atoms):")
    print("  z_PM_old       = exp(μ_G M/R) − 1  ≈  2GM/c²R  at leading order")
    print()
    print("Corrected formula (particle-metric atom clock + conserved coord-ω):")
    print("  z_PM_corrected = 1/sqrt(1 − μ_G M/R) − 1  =  z_GR   (EXACTLY)")
    print()
    print("Resolution: PM has two metrics. Atomic spectral lines are emitted by")
    print("massive particles (particle metric g_tt^mass = −(1−φ)), not by photon")
    print("clocks (optical metric g_tt^opt = −c²/n²). Coordinate photon frequency")
    print("is conserved regardless (static medium, Killing vector ∂_t).")
    print()

    # ── 2. Numerical verification ─────────────────────────────────────────────
    print("NUMERICAL VERIFICATION")
    hdr = (
        f"{'M/M☉':>7} {'R/km':>6} {'φ_s':>7} "
        f"{'z_old':>8} {'z_corr':>8} {'z_GR':>8} "
        f"{'old/GR':>8} {'corr/GR':>9} {'=GR?':>6}"
    )
    print(hdr)
    print("-" * 76)

    test_cases = [
        # (M/M_sun, R/km) — spanning the observational range
        (1.34, 12.71),   # NICER J0030+0451 (Riley 2019)
        (2.08, 12.39),   # NICER J0740+6620 (Riley 2021)
        (1.0,  14.0),
        (1.4,  13.2),    # joint-feasibility PM (f=1.18)
        (1.4,  10.0),    # compact end
        (2.0,  11.0),
        (3.0,  15.0),
        (0.001, 100.0),  # near-Newtonian limit
    ]

    for M, R in test_cases:
        r = _compare(M, R)
        flag = "YES" if r['equals_gr'] else (" NaN" if math.isnan(r['z_corr']) else "  NO")
        rc = f"{r['ratio_corr']:>9.6f}" if math.isfinite(r['ratio_corr']) else f"{'NaN':>9}"
        print(
            f"{r['M']:>7.3f} {r['R']:>6.2f} {r['phi_s']:>7.4f} "
            f"{r['z_old']:>8.5f} "
            f"{r['z_corr']:>8.5f} "
            f"{r['z_gr']:>8.5f} "
            f"{r['ratio_old']:>8.4f} "
            f"{rc} "
            f"{flag:>6}"
        )

    print("-" * 76)
    print()

    # ── 3. Domain of validity ─────────────────────────────────────────────────
    print("DOMAIN OF VALIDITY")
    print("  Formula requires φ_s = 2GM/c²R < 1  (compactness C < 0.5).")
    print("  NICER neutron stars: C ≈ 0.15–0.25  →  valid ✓")
    print("  PM max-mass stars can reach C ≈ 0.74 (baseline), where the particle")
    print("  metric g_tt^mass → 0 at φ_s = 1.  A higher-order g_tt^mass would be")
    print("  needed there; the function returns nan for those inputs.")
    print()

    # ── 4. Weak-field check ───────────────────────────────────────────────────
    print("WEAK-FIELD LIMIT CHECK")
    M = 0.001 * M_SUN
    R = 100e3
    x = MU_G * M / R
    z_c = pm_surface_redshift_corrected(M, R)
    z_g = gr_surface_redshift(M, R)
    print(f"  M = 0.001 M☉, R = 100 km  →  φ_s = {x:.4e}")
    print(f"  z_corr = {z_c:.6e}   (expect ≈ {x:.6e}  =  GM/c²R)")
    print(f"  z_GR   = {z_g:.6e}")
    print(f"  ratio  = {z_c/z_g:.8f}  (expect 1.000000)")
    print()

    # ── 5. Summary ────────────────────────────────────────────────────────────
    print("SUMMARY")
    print("-------")
    print("  z_PM_corrected = z_GR  to machine precision for all C < 0.5.")
    print("  The 1.6–1.8× discrepancy was a derivation error, not a prediction.")
    print("  Solar-system tests (light bending, Shapiro delay) are unaffected:")
    print("  those are purely optical-metric (photon) phenomena and μ_G is correct.")


if __name__ == "__main__":
    main()
