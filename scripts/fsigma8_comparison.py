#!/usr/bin/env python3
"""PM vs ΛCDM: fσ8(z) structure-growth comparison against RSD measurements.

This is one of the most constraining cosmological tests:
  - fσ8(z) = f(z) × σ8(z)  is directly measured by redshift-space distortions
  - ΛCDM predicts a nearly flat curve ~0.40–0.50 over z=0–1.5
  - PM Scenario A (static, no friction): predicts constant f ≈ 0.69 with a
    β-dependent z-slope — shape is DIFFERENT from ΛCDM
  - PM Scenario B (drag, β=0.8): predicts slower growth rate in MD;
    fσ8 shape deviates by >10% from ΛCDM at z>0.5

Key finding
──────────────────────────────────────────────────────────────────────────────
PM cannot simultaneously:
  1. Match D_L(z) data  [β=0.8 calibrated here]
  2. Match fσ8(z) data  [requires different β or Ω_eff — degeneracy]

This is NOT an immediate falsification (PM could tune Ω_eff independently),
but it shows PM has no free-parameter-free prediction for structure growth.
A measurement with σ(fσ8) ~ 0.01 over multiple bins would tightly constrain
the [β, Ω_eff] plane and cross-check the D_L calibration.
──────────────────────────────────────────────────────────────────────────────
"""

import sys
import math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from general_relativity.cosmology import (
    lcdm_hubble,
    lcdm_luminosity_distance,
    lcdm_fsigma8,
    H0_PLANCK, OM0_PLANCK, SIGMA8_PLANCK,
)
from pushing_medium.cosmology import (
    pm_hubble_effective,
    pm_luminosity_distance,
    pm_static_fsigma8,
    pm_drag_fsigma8,
    pm_drag_growing_power,
    BETA_CALIB,
)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ── Published RSD fσ8 data ────────────────────────────────────────────────────
# Compilation from Nesseris et al. 2017 (PRD 96, 023542) and originals.
# Format: (z_eff, fσ8, σ_{fσ8}, survey, reference)
RSD_DATA = [
    (0.02,   0.428,  0.047, "2MTF",        "Howlett+2017"),
    (0.067,  0.423,  0.055, "6dFGS",       "Beutler+2012"),
    (0.10,   0.370,  0.130, "SDSS veloc.", "Hudson+2012"),
    (0.15,   0.490,  0.145, "SDSS MGS",   "Howlett+2015"),
    (0.22,   0.420,  0.070, "WiggleZ",    "Blake+2012"),
    (0.25,   0.351,  0.058, "SDSS LRG",   "Samushia+2012"),
    (0.32,   0.427,  0.056, "BOSS LOWZ",  "Gil-Marin+2016"),
    (0.37,   0.460,  0.038, "SDSS DR7",   "Samushia+2012"),
    (0.44,   0.413,  0.080, "WiggleZ",    "Blake+2012"),
    (0.57,   0.427,  0.066, "BOSS CMASS", "Beutler+2014"),
    (0.60,   0.390,  0.063, "WiggleZ",    "Blake+2012"),
    (0.78,   0.380,  0.044, "WiggleZ",    "Blake+2012"),
    (0.80,   0.470,  0.080, "VIPERS",     "de la Torre+2013"),
    (1.40,   0.482,  0.116, "FastSound",  "Okumura+2016"),
]

z_data   = np.array([d[0] for d in RSD_DATA])
fs8_data = np.array([d[1] for d in RSD_DATA])
err_data = np.array([d[2] for d in RSD_DATA])


def chi2(fs8_pred_at_data, fs8_obs, err_obs):
    return np.sum(((fs8_pred_at_data - fs8_obs) / err_obs)**2)


def build_predictions(z_points):
    """Build fσ8 predictions for all models at the given z values.

    PM-Drag with Ω_eff=1.0 and β=0.8 raises RuntimeError because there is
    no valid matter-dominated era (explosive early growth) — this is itself
    a falsification and is returned as NaN.
    """
    # ΛCDM (Planck 2018)
    fs8_lcdm = lcdm_fsigma8(z_points)

    # PM-Static: β=0.8, Ω_m=0.315
    fs8_static = pm_static_fsigma8(z_points, Om0=OM0_PLANCK, beta=BETA_CALIB)

    # PM-Drag: β=0.8, Ω_eff=0.315 (same matter density as ΛCDM)
    fs8_drag_315 = pm_drag_fsigma8(z_points, Om_eff=0.315, beta=BETA_CALIB)

    # PM-Drag: β=0.8, Ω_eff=1.0 (no dark energy component)
    # NOTE: This scenario has no valid MD era — ODE fails, meaning PM with
    # full Ω_eff=1 and β=0.8 predicts catastrophic over-growth at high z.
    try:
        fs8_drag_1 = pm_drag_fsigma8(z_points, Om_eff=1.0, beta=BETA_CALIB)
    except RuntimeError:
        fs8_drag_1 = np.full(len(np.atleast_1d(z_points)), np.nan)

    return fs8_lcdm, fs8_static, fs8_drag_315, fs8_drag_1


def print_fsigma8_table():
    print("\n" + "=" * 100)
    print("  fσ8(z) Predictions: ΛCDM vs PM Scenarios")
    print("=" * 100)
    print(f"  {'z':<7} {'ΛCDM':<10} {'PM-Static':<12} {'PM-Drag Ω=0.315':<18} {'PM-Drag Ω=1':<14} {'Observed'}")
    print("  " + "-" * 97)

    z_plot = np.array([0.067, 0.15, 0.25, 0.32, 0.44, 0.57, 0.78, 1.40])
    lcdm, static, drag315, drag1 = build_predictions(z_plot)

    for i, z in enumerate(z_plot):
        # Find matching data point if available
        obs_str = ""
        for d in RSD_DATA:
            if abs(d[0] - z) < 0.01:
                obs_str = f"{d[1]:.3f}±{d[2]:.3f} ({d[3]})"
                break
        print(f"  {z:<7.3f} {lcdm[i]:<10.4f} {static[i]:<12.4f} {drag315[i]:<18.4f} "
              f"{'N/A':<14} {obs_str}")
    print()
    print("  NOTE: PM-Drag Ω_eff=1.0 omitted — ODE fails (runaway growth, no MD era)")

    print()
    print(f"  σ8,0 = {SIGMA8_PLANCK},  H0 = {H0_PLANCK},  β_calib = {BETA_CALIB}")
    print("=" * 100 + "\n")


def print_chi2_table():
    """Evaluate χ² for each model relative to observed RSD data."""
    lcdm_pred    = lcdm_fsigma8(z_data)
    static_pred  = pm_static_fsigma8(z_data, Om0=OM0_PLANCK, beta=BETA_CALIB)
    drag315_pred = pm_drag_fsigma8(z_data, Om_eff=0.315, beta=BETA_CALIB)
    try:
        drag1_pred = pm_drag_fsigma8(z_data, Om_eff=1.0, beta=BETA_CALIB)
        drag1_c2   = chi2(drag1_pred, fs8_data, err_data)
        drag1_str  = f"{drag1_c2:.2f}"
    except RuntimeError:
        drag1_c2   = np.inf
        drag1_str  = "N/A (runaway growth — no MD era)"

    n = len(z_data)

    rows = [
        ("ΛCDM (Planck fixed)",        chi2(lcdm_pred,    fs8_data, err_data)),
        ("PM-Static β=0.8",            chi2(static_pred,  fs8_data, err_data)),
        ("PM-Drag β=0.8, Ω_eff=0.315", chi2(drag315_pred, fs8_data, err_data)),
        ("PM-Drag β=0.8, Ω_eff=1.0",   drag1_c2),
    ]

    print("=" * 65)
    print("  χ² Goodness of Fit  (N_data = {})".format(n))
    print("=" * 65)
    print(f"  {'Model':<35} {'χ²':<22} {'χ²/N':<10}")
    print("  " + "-" * 62)
    for name, c2 in rows:
        if np.isinf(c2):
            print(f"  {name:<35} {drag1_str:<22} {'—':10}")
        else:
            print(f"  {name:<35} {c2:<22.2f} {c2/n:<10.3f}")
    print("=" * 65 + "\n")

    return rows


def print_dl_match_table():
    """Confirm the β=0.8 D_L calibration from COSMOLOGY_FINDINGS.md."""
    print("=" * 75)
    print("  Luminosity Distance Calibration: PM β=0.8 vs ΛCDM")
    print("=" * 75)
    print(f"  {'z':<8} {'D_L ΛCDM [Mpc]':<20} {'D_L PM β=0.8 [Mpc]':<22} {'Diff %'}")
    print("  " + "-" * 72)
    for z in [0.1, 0.5, 1.0, 1.5]:
        dl_lcdm = lcdm_luminosity_distance(z)
        dl_pm   = pm_luminosity_distance(z,  beta=BETA_CALIB)
        diff    = (dl_pm / dl_lcdm - 1.0) * 100
        print(f"  {z:<8.1f} {dl_lcdm:<20.0f} {dl_pm:<22.0f} {diff:+.2f}%")
    print("  (β=0.8 calibrated to match D_L within ~3%, but fσ8 shape differs)")
    print("=" * 75 + "\n")


def print_pm_falsification_summary(chi2_rows):
    lcdm_c2   = chi2_rows[0][1]
    static_c2 = chi2_rows[1][1]
    drag315   = chi2_rows[2][1]
    drag1     = chi2_rows[3][1]   # may be inf
    n = len(z_data)

    print("COSMOLOGY FALSIFICATION SUMMARY")
    print("=" * 70)
    print(f"  Data: {n} published RSD fσ8 measurements, z = 0.02–1.40")
    print()
    print(f"  ΛCDM (fiducial):           χ²/N = {lcdm_c2/n:.3f}  ✓ REFERENCE")
    print(f"  PM-Static (no friction):   χ²/N = {static_c2/n:.3f}", end="")
    print("  ✗ POOR FIT" if static_c2/n > 3*lcdm_c2/n else "  ≈ marginal")
    print(f"  PM-Drag Ω_eff=0.315:       χ²/N = {drag315/n:.3f}", end="")
    print("  ✗ POOR FIT" if drag315/n > 2*lcdm_c2/n else "  ≈ marginal")
    if np.isinf(drag1):
        print(f"  PM-Drag Ω_eff=1.0:         χ²/N = N/A  ✗ NO VALID MD ERA (ODE failure = runaway growth)")
    else:
        print(f"  PM-Drag Ω_eff=1.0:         χ²/N = {drag1/n:.3f}", end="")
        print("  ✗ POOR FIT" if drag1/n > 2*lcdm_c2/n else "  ≈ marginal")
    print()
    print("  Key conclusion:")
    print("  PM has NO parameter-free prediction for structure growth.")
    print("  β calibrated on D_L (luminosity distance) does NOT")
    print("  simultaneously reproduce fσ8 without separately tuning Ω_eff.")
    print()
    print("  Best-fit PM might be found at ~Ω_eff ≈ 0.1–0.2 for β=0.8,")
    print("  but this would imply a very matter-underdense universe.")
    print("  A full joint fit [β, Ω_eff, σ8,0] with D_L + fσ8 data")
    print("  would tightly constrain or rule out the refractive model.")
    print("=" * 70 + "\n")


def make_plot(chi2_rows):
    if not HAS_MPL:
        return

    z_fine = np.linspace(0.01, 1.55, 150)
    lcdm_f, static_f, drag315_f, drag1_f = build_predictions(z_fine)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("PM vs ΛCDM: fσ8(z) Structure Growth Comparison", fontsize=13)

    # ── Left: fσ8 curves + data ───────────────────────────────────────────
    ax = axes[0]
    ax.errorbar(z_data, fs8_data, yerr=err_data, fmt='ko', ms=4, capsize=3,
                label='RSD data (14 surveys)', zorder=5, alpha=0.9)
    ax.plot(z_fine, lcdm_f,    '-',  color='steelblue',  lw=2.0, label='ΛCDM (Planck 2018)')
    ax.plot(z_fine, static_f,  '--', color='darkorange',  lw=1.8, label=f'PM-Static β={BETA_CALIB}')
    ax.plot(z_fine, drag315_f, '-.',  color='crimson',    lw=1.8, label=f'PM-Drag β={BETA_CALIB}, Ω_eff=0.315')
    ax.plot(z_fine, drag1_f,   ':',  color='purple',     lw=1.8, label=f'PM-Drag β={BETA_CALIB}, Ω_eff=1.0 (ODE fails)')
    ax.set_xlabel("Redshift z", fontsize=11)
    ax.set_ylabel("fσ8(z)", fontsize=11)
    ax.set_xlim(-0.05, 1.6)
    ax.set_ylim(0.0, 0.95)
    ax.legend(fontsize=8.5, loc='upper right')
    ax.grid(True, ls=':', lw=0.5, alpha=0.6)
    ax.set_title("fσ8 vs redshift\n(PM curves use β=0.8 calibrated on D_L)", fontsize=10)

    # ── Right: fractional residual vs ΛCDM ────────────────────────────────
    ax2 = axes[1]
    ax2.axhline(0, color='steelblue', lw=1.5, label='ΛCDM (reference)')
    ax2.fill_between(z_fine, -0.1, 0.1, alpha=0.08, color='steelblue')
    ax2.plot(z_fine, (static_f  / lcdm_f - 1)*100, '--', color='darkorange', lw=1.8,
             label=f'PM-Static')
    ax2.plot(z_fine, (drag315_f / lcdm_f - 1)*100, '-.', color='crimson',    lw=1.8,
             label=f'PM-Drag Ω_eff=0.315')
    ax2.plot(z_fine, (drag1_f   / lcdm_f - 1)*100, ':',  color='purple',     lw=1.8,
             label=f'PM-Drag Ω_eff=1.0')

    # Data residuals
    lcdm_at_data = lcdm_fsigma8(z_data)
    data_res = (fs8_data / lcdm_at_data - 1.0) * 100
    data_err_pct = err_data / lcdm_at_data * 100
    ax2.errorbar(z_data, data_res, yerr=data_err_pct, fmt='ko', ms=4,
                 capsize=3, label='Data residuals', zorder=5, alpha=0.9)

    ax2.axhline(-10, color='gray', lw=0.8, ls=':')
    ax2.axhline(+10, color='gray', lw=0.8, ls=':')
    ax2.set_xlabel("Redshift z", fontsize=11)
    ax2.set_ylabel("(model − ΛCDM) / ΛCDM  [%]", fontsize=11)
    ax2.set_xlim(-0.05, 1.6)
    ax2.set_ylim(-60, 120)
    ax2.legend(fontsize=8.5)
    ax2.grid(True, ls=':', lw=0.5, alpha=0.6)
    ax2.set_title("Fractional deviation from ΛCDM\n(>10%: potentially measurable)", fontsize=10)

    plt.tight_layout()
    out = Path(__file__).parent.parent / 'results' / 'fsigma8_comparison.png'
    out.parent.mkdir(exist_ok=True)
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {out}")
    plt.close()


def main():
    print_dl_match_table()
    print_fsigma8_table()
    chi2_rows = print_chi2_table()
    make_plot(chi2_rows)
    print_pm_falsification_summary(chi2_rows)


if __name__ == "__main__":
    main()
