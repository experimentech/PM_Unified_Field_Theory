#!/usr/bin/env python3
"""PM surface gravitational redshift vs mass track — comparison with GR and NICER.

Computes and plots:
  1. PM and GR M–R curves (different structure equations → different tracks)
  2. PM and GR surface redshift z(M) along each track
  3. NICER M-R constraint boxes (J0030+0451, J0740+6620)
  4. Comparison table: PM vs GR redshift at NICER central values

Key finding
-----------
PM predicts roughly TWICE the surface redshift of GR for the same (M, R),
because z_PM = e^{2GM/c²R} − 1 ≈ 2GM/c²R while z_GR ≈ GM/c²R.
This is a genuine, falsifiable prediction testable with X-ray spectroscopy.
"""

import sys
import math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from pushing_medium.stellar_structure import (
    compute_mr_curve,
    pm_surface_redshift,
    compute_surface_redshift_track,
    G, c, M_SUN, MU_G,
)
from general_relativity.classical import gr_surface_redshift

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ── NICER reference data ─────────────────────────────────────────────────────
# Riley et al. (2019) PSR J0030+0451
J0030_M   = 1.34;  J0030_M_err  = 0.16   # M_sun
J0030_R   = 12.71; J0030_R_err  = 1.17   # km (mean of +1.14/-1.19)

# Riley et al. (2021) + Fonseca et al. (2021) PSR J0740+6620
J0740_M   = 2.08;  J0740_M_err  = 0.07   # M_sun
J0740_R   = 12.39; J0740_R_err  = 1.14   # km (mean of +1.30/-0.98)

# Cottam et al. (2002) EXO 0748-676 — spectral line redshift (indicative)
EXO_Z     = 0.35


def _gr_mr_curve(n_points=60):
    """GR M-R curve for a simple stiff EOS (P=c²(ρ-ρ_nuc)/2, TOV structure)."""
    # Import the full TOV integrator if present; otherwise use analytic approximation
    try:
        from pushing_medium.stellar_structure import RHO_NUC
        from scipy.integrate import solve_ivp

        RHO_CRIT = math.e * RHO_NUC
        rho_c_arr = np.linspace(1.01 * RHO_NUC, 0.999 * RHO_CRIT, n_points)
        M_arr = np.full(n_points, np.nan)
        R_arr = np.full(n_points, np.nan)

        for i, rho_c in enumerate(rho_c_arr):
            phi_c = math.log(rho_c / RHO_NUC)
            # PM EOS in TOV equations (GR structure, same PM EOS)
            def tov(r, y):
                m, P = y
                if P <= 0 or r < 1.0:
                    return [0.0, 0.0]
                rho = RHO_NUC + 2.0 * P / (c * c)
                dP = -(G * m * rho / r ** 2) * (
                    (1.0 + P / (rho * c * c))
                    * (1.0 + 4.0 * math.pi * r ** 3 * P / (m * c * c))
                    / (1.0 - 2.0 * G * m / (c * c * r))
                )
                dm = 4.0 * math.pi * r ** 2 * rho
                return [dm, dP]

            def surface(r, y):
                return y[1]
            surface.terminal  = True
            surface.direction = -1

            P0 = c * c / 2.0 * (rho_c - RHO_NUC)
            y0 = [1e-6, P0]
            r_span = (10.0, 3e4)
            sol = solve_ivp(tov, r_span, y0, method='RK45', events=surface,
                            dense_output=False, rtol=1e-8, atol=1e-30,
                            max_step=50.0)
            if sol.t_events[0].size > 0:
                idx = -1
                M_arr[i] = sol.y_events[0][0][0] / M_SUN
                R_arr[i] = sol.t_events[0][0] / 1e3
            elif sol.success:
                M_arr[i] = sol.y[0, -1] / M_SUN
                R_arr[i] = sol.t[-1] / 1e3

        return rho_c_arr, M_arr, R_arr
    except Exception:
        return None, None, None


def print_comparison_table():
    """Print PM vs GR surface redshifts at NICER central values."""
    print("\n" + "=" * 68)
    print("  Surface Gravitational Redshift: PM vs GR vs NICER")
    print("=" * 68)
    print(f"  {'Source':<22} {'M (M☉)':<9} {'R (km)':<9} {'z_PM':<10} {'z_GR':<10} {'z_PM/z_GR'}")
    print(f"  {'-'*65}")

    cases = [
        ("J0030+0451 (NICER)", J0030_M, J0030_R),
        ("J0740+6620 (NICER)", J0740_M, J0740_R),
        ("Typical NS (1.4, 11)", 1.4, 11.0),
        ("Stiff NS   (2.0, 13)", 2.0, 13.0),
    ]
    for name, M_msun, R_km in cases:
        M = M_msun * M_SUN
        R = R_km * 1e3
        z_pm = pm_surface_redshift(M, R)
        z_gr = gr_surface_redshift(M, R)
        print(f"  {name:<22} {M_msun:<9.2f} {R_km:<9.2f} {z_pm:<10.4f} {z_gr:<10.4f} {z_pm/z_gr:.3f}×")

    print()
    print(f"  EXO 0748-676 observed (Cottam 2002): z_obs = {EXO_Z} (indicative)")
    print(f"    → implies 2GM/c²R ≈ {math.log(1 + EXO_Z):.3f}  (PM),  "
          f"or ≈ {1 - 1/(1+EXO_Z)**2:.3f}  (GR)")
    print("=" * 68 + "\n")


def main():
    print("Computing PM M-R track...")
    rho_c_pm, M_pm, R_pm, z_pm = compute_surface_redshift_track(n_points=50)

    mask_pm = np.isfinite(M_pm) & np.isfinite(R_pm) & np.isfinite(z_pm)
    M_pm = M_pm[mask_pm]
    R_pm = R_pm[mask_pm]
    z_pm = z_pm[mask_pm]

    print_comparison_table()

    if not HAS_MPL:
        print("matplotlib not available — skipping plots.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle("PM Compact Stars: M–R Curve and Surface Redshift", fontsize=13)

    # ── Left: M-R diagram ────────────────────────────────────────────────
    ax = axes[0]

    # NICER constraints (approximate 1σ boxes)
    for (M_c, dM, R_c, dR, label, col) in [
        (J0030_M, J0030_M_err, J0030_R, J0030_R_err, "J0030+0451\n(NICER 2019)", "#2196F3"),
        (J0740_M, J0740_M_err, J0740_R, J0740_R_err, "J0740+6620\n(NICER 2021)", "#4CAF50"),
    ]:
        rect = mpatches.FancyBboxPatch(
            (R_c - dR, M_c - dM), 2*dR, 2*dM,
            boxstyle="square,pad=0", linewidth=1.5,
            edgecolor=col, facecolor=col, alpha=0.18,
        )
        ax.add_patch(rect)
        ax.text(R_c, M_c + dM + 0.04, label, ha='center', va='bottom', fontsize=7, color=col)

    ax.plot(R_pm, M_pm, 'r-', lw=2.0, label='PM (exact flat-space)')

    ax.set_xlabel("Radius R [km]", fontsize=11)
    ax.set_ylabel("Mass M [M☉]", fontsize=11)
    ax.set_xlim(5, 20)
    ax.set_ylim(0, 2.5)
    ax.legend(fontsize=9)
    ax.grid(True, ls=':', lw=0.5, alpha=0.7)
    ax.set_title("Mass–Radius Diagram", fontsize=11)

    # ── Right: z_surface vs M ────────────────────────────────────────────
    ax2 = axes[1]

    # PM track
    ax2.plot(M_pm, z_pm, 'r-', lw=2.0, label='PM: $z = e^{\\mu_G M/R} - 1$')

    # GR formula at same (M, R) as PM track
    z_gr_on_pm = np.array([
        gr_surface_redshift(m * M_SUN, r * 1e3) for m, r in zip(M_pm, R_pm)
    ])
    ax2.plot(M_pm, z_gr_on_pm, 'b--', lw=1.8,
             label='GR: $1/\\sqrt{1-2GM/c^2R}-1$\n(at PM M,R)')

    # EXO 0748-676 redshift band
    ax2.axhspan(EXO_Z - 0.05, EXO_Z + 0.05, alpha=0.15, color='gray',
                label=f'EXO 0748-676: z≈{EXO_Z} (indicative)')

    ax2.set_xlabel("Mass M [M☉]", fontsize=11)
    ax2.set_ylabel("Surface redshift  z", fontsize=11)
    ax2.legend(fontsize=8)
    ax2.grid(True, ls=':', lw=0.5, alpha=0.7)
    ax2.set_title("Surface Redshift vs Mass\n(PM predicts ~2× more than GR)", fontsize=11)

    plt.tight_layout()
    out = Path(__file__).parent.parent / 'results' / 'surface_redshift_track.png'
    out.parent.mkdir(exist_ok=True)
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {out}")
    plt.close()


if __name__ == "__main__":
    main()
