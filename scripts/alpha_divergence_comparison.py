"""Alpha-measure divergence comparison: where PM predictions split.

This script computes M-R curves for all three action-measure variants of PM:

  α=0  Standard PM      : ∇²φ = −κρ n          (coordinate measure d³x)
  α=2  N-field / area   : ∇²n = −κρ n²         (area measure n² d³x)
  α=3  Physical-measure : ∇²w = −κρ w^{4/3}    (volume measure n³ d³x,  w=n^{3/2})

All three are identical at low density (n→1, φ→0).  They diverge at high φ
because the action weight e^{αφ} departs significantly from 1.

Outputs printed:
  1. M_max and R(M_max) for each variant
  2. The M and φ_central where α=2 diverges >5% from α=0
  3. The M and φ_central where α=3 diverges >5% from α=0
  4. A text table of (M, R, φ_central) for each variant at shared ρ_c values
  5. GW event masses overlaid in the table for context
"""

import sys
import os
import numpy as np

# ── path setup ──────────────────────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.pushing_medium.stellar_structure import (
    G, c, M_SUN, RHO_NUC, RHO_CRIT,
    compute_mr_curve,
    compute_mr_curve_nfield,
    compute_mr_curve_physical_measure,
    solve_pm_star,
    solve_pm_star_nfield,
    solve_pm_star_physical_measure,
)

# ── GW events for reference ──────────────────────────────────────────────────
GW_EVENTS = {
    "GW170817 (BNS, primary)":  1.46,
    "GW170817 (BNS, secondary)": 1.27,
    "GW190814 (secondary)": 2.6,
    "GW190814 (primary)":  22.2,
    "GW150914 (secondary)": 30.6,
    "GW150914 (primary)":  35.6,
}

N_POINTS = 80


def run():
    print("=" * 72)
    print("PM α-Measure Divergence Analysis")
    print("=" * 72)

    # ── compute M-R curves ────────────────────────────────────────────────────
    print("\nComputing M-R curves (this may take ~30 s)...")

    rho_c0, M0, R0 = compute_mr_curve(n_points=N_POINTS,
                                       rho_max_factor=0.999)
    rho_c2, M2, R2 = compute_mr_curve_nfield(n_points=N_POINTS,
                                              rho_max_factor=0.999)
    rho_c3, M3, R3 = compute_mr_curve_physical_measure(n_points=N_POINTS,
                                                        rho_max_factor=0.999)

    # φ_central for each model
    phi_c0 = np.log(rho_c0 / RHO_NUC)
    phi_c2 = np.log(rho_c2 / RHO_NUC)
    phi_c3 = np.log(rho_c3 / RHO_NUC)

    # ── M_max summary ─────────────────────────────────────────────────────────
    def peak(M, R):
        valid = np.isfinite(M) & np.isfinite(R)
        if not valid.any():
            return np.nan, np.nan
        idx = np.nanargmax(M[valid])
        return np.nanmax(M[valid]), R[valid][idx]

    Mmax0, Rmax0 = peak(M0, R0)
    Mmax2, Rmax2 = peak(M2, R2)
    Mmax3, Rmax3 = peak(M3, R3)

    print("\n── M_max comparison ─────────────────────────────────────────────")
    print(f"  α=0  standard PM      : M_max = {Mmax0:6.3f} M☉,  R = {Rmax0:5.2f} km")
    print(f"  α=2  n-field/area     : M_max = {Mmax2:6.3f} M☉,  R = {Rmax2:5.2f} km")
    print(f"  α=3  physical-measure : M_max = {Mmax3:6.3f} M☉,  R = {Rmax3:5.2f} km")

    print("\n  GW event reference masses:")
    for name, mass in sorted(GW_EVENTS.items(), key=lambda x: x[1]):
        ok  = [f"α={a}" for a, mx in ((0, Mmax0), (2, Mmax2), (3, Mmax3)) if mass <= mx]
        exc = [f"α={a}" for a, mx in ((0, Mmax0), (2, Mmax2), (3, Mmax3)) if mass > mx]
        if not ok:
            status = "exceeds ALL PM variants"
        elif not exc:
            status = "within ALL PM variants"
        else:
            status = f"within {', '.join(ok)}; exceeds {', '.join(exc)}"
        print(f"    {name:35s}: {mass:5.2f} M☉  ({status})")

    # ── divergence onset ──────────────────────────────────────────────────────
    print("\n── Divergence onset (α=2 vs α=0, α=3 vs α=0) ───────────────────")
    print("  Measured at shared φ_central grid (interpolating M onto φ grid)\n")

    # Build φ grid from α=0 curve (reference)
    # Use only valid points
    valid0 = np.isfinite(M0) & np.isfinite(R0)
    phi0_v = phi_c0[valid0]
    M0_v   = M0[valid0]
    R0_v   = R0[valid0]

    valid2 = np.isfinite(M2) & np.isfinite(R2)
    phi2_v = phi_c2[valid2]
    M2_v   = M2[valid2]

    valid3 = np.isfinite(M3) & np.isfinite(R3)
    phi3_v = phi_c3[valid3]
    M3_v   = M3[valid3]

    # Interpolate α=2 and α=3 onto α=0's φ grid
    phi_common = phi0_v[phi0_v <= min(phi2_v.max(), phi3_v.max())]
    M0_interp  = np.interp(phi_common, phi0_v, M0_v)
    M2_interp  = np.interp(phi_common, phi2_v, M2_v)
    M3_interp  = np.interp(phi_common, phi3_v, M3_v)

    rel_diff_2 = np.abs(M2_interp - M0_interp) / (M0_interp + 1e-30)
    rel_diff_3 = np.abs(M3_interp - M0_interp) / (M0_interp + 1e-30)

    THRESHOLD = 0.05  # 5%
    idx_div2 = np.where(rel_diff_2 > THRESHOLD)[0]
    idx_div3 = np.where(rel_diff_3 > THRESHOLD)[0]

    if idx_div2.size:
        i = idx_div2[0]
        print(f"  α=2 diverges >5% from α=0 at:")
        print(f"    φ_central = {phi_common[i]:.3f}  "
              f"(n_c = e^φ = {np.exp(phi_common[i]):.3f},  "
              f"ρ_c/ρ_nuc = {np.exp(phi_common[i]):.3f})")
        print(f"    M(α=0) = {M0_interp[i]:.3f} M☉,  "
              f"M(α=2) = {M2_interp[i]:.3f} M☉  "
              f"(diff = {rel_diff_2[i]*100:.1f}%)")
    else:
        print("  α=2 vs α=0: never diverges >5% in the valid range")

    if idx_div3.size:
        i = idx_div3[0]
        print(f"\n  α=3 diverges >5% from α=0 at:")
        print(f"    φ_central = {phi_common[i]:.3f}  "
              f"(n_c = e^φ = {np.exp(phi_common[i]):.3f},  "
              f"ρ_c/ρ_nuc = {np.exp(phi_common[i]):.3f})")
        print(f"    M(α=0) = {M0_interp[i]:.3f} M☉,  "
              f"M(α=3) = {M3_interp[i]:.3f} M☉  "
              f"(diff = {rel_diff_3[i]*100:.1f}%)")
    else:
        print("  α=3 vs α=0: never diverges >5% in the valid range")

    # ── detailed table at key φ values ───────────────────────────────────────
    print("\n── Comparison table at selected φ_central values ────────────────")
    print(f"  {'φ_c':>6}  {'n_c':>6}  {'M(α=0)':>9}  {'R(α=0)':>8}  "
          f"{'M(α=2)':>9}  {'R(α=2)':>8}  {'M(α=3)':>9}  {'R(α=3)':>8}  "
          f"{'Δ(2-0)%':>8}  {'Δ(3-0)%':>8}")
    print("  " + "-" * 85)

    phi_table = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99])

    for phi_t in phi_table:
        rho_t = RHO_NUC * np.exp(phi_t)
        if rho_t > RHO_CRIT * 0.999:
            continue
        try:
            s0 = solve_pm_star(rho_t)
            s2 = solve_pm_star_nfield(rho_t)
            s3 = solve_pm_star_physical_measure(rho_t)
            if not (s0.converged and s2.converged and s3.converged):
                continue
            m0 = s0.M_star / M_SUN;  r0_km = s0.R_star / 1e3
            m2 = s2.M_star / M_SUN;  r2_km = s2.R_star / 1e3
            m3 = s3.M_star / M_SUN;  r3_km = s3.R_star / 1e3
            d20 = (m2 - m0) / (m0 + 1e-30) * 100
            d30 = (m3 - m0) / (m0 + 1e-30) * 100
            print(f"  {phi_t:6.2f}  {np.exp(phi_t):6.3f}  "
                  f"{m0:9.4f}  {r0_km:8.3f}  "
                  f"{m2:9.4f}  {r2_km:8.3f}  "
                  f"{m3:9.4f}  {r3_km:8.3f}  "
                  f"{d20:+8.2f}  {d30:+8.2f}")
        except Exception as e:
            print(f"  {phi_t:6.2f}  [error: {e}]")

    # ── measure weight at M_max for each variant ──────────────────────────────
    print("\n── Action-weight factor e^{αφ} at M_max central compression ──────")
    print("  (this is the ratio by which the two actions differ at M_max)")
    for label, alpha, rho_c_arr, M_arr in [
        ("α=0 standard", 0, rho_c0, M0),
        ("α=2 n-field",  2, rho_c2, M2),
        ("α=3 phys-meas",3, rho_c3, M3),
    ]:
        valid = np.isfinite(M_arr)
        if not valid.any():
            continue
        idx_max = np.nanargmax(M_arr[valid])
        phi_at_max = np.log(rho_c_arr[valid][idx_max] / RHO_NUC)
        weight = np.exp(alpha * phi_at_max)
        print(f"  {label:20s}: φ_central = {phi_at_max:.3f},  "
              f"e^{{α·φ}} = e^{{{alpha}×{phi_at_max:.3f}}} = {weight:.2f}")

    print("\n" + "=" * 72)
    print("Interpretation:")
    print("  The action weight at M_max is e^{αφ_max}.  When this factor is")
    print("  large, the three variants weight the stellar interior differently")
    print("  and predict different M_max.  The divergence onset (5% threshold)")
    print("  above marks the mass scale where the choice of α becomes the physics.")
    print("=" * 72)


if __name__ == "__main__":
    run()
