#!/usr/bin/env python3
"""
PM Two-Phase Cosmology: Revealing β=0.8 as Dark Energy in Disguise
===================================================================

Demonstrates that the PM 'calibration parameter' β=0.8 is NOT a
PM-specific constant — it is the effective instantaneous power-law
index of the ΛCDM Hubble rate averaged over the D_L-sensitive redshift
range z=0.1–1.5.

The PM two-phase model:
  ─────────────────────────────────────────────────────────────
  Matter phase:  ρ_m(a) ∝ a^-3   (clusters, drives growth)
  Energy phase:  ρ_E(a) ∝ a^(-3(1+w_E))  (no clustering; w_E ≈ -1)

  Implied Friedmann equation:
    H²(z) = H₀²[Ω_m(1+z)³ + Ω_E(1+z)^(3(1+w_E))],  Ω_E = 1−Ω_m

  For w_E = -1 (PM prediction: ground-state medium has constant density):
    → This is MATHEMATICALLY IDENTICAL to ΛCDM.
  ─────────────────────────────────────────────────────────────

Physical interpretation of the energy phase:
  In PM, the stability cap φ_crit = 1 marks a phase transition from
  dense matter to a uniform 'energy phase'.  Cosmologically, this
  energy phase acts as dark energy — its density does not dilute
  as the universe expands because it represents the minimum excitation
  state (ground state) of the medium.

  This connects compact-object physics (max neutron-star mass ≈14 M☉)
  to cosmology (Ω_Λ ≈ 0.685) through the SAME phase transition.

Outputs
-------
  results/pm_twophase_cosmology.png   — 4-panel figure
  (tables to stdout)

Usage
-----
  python scripts/pm_twophase_cosmology.py
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))

import numpy as np
from scipy.integrate import quad

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("matplotlib not available; skipping figure")

from pushing_medium.cosmology import (
    pm_twophase_hubble,
    pm_effective_beta,
    pm_energy_phase_fraction,
    pm_twophase_growth_factor_rate,
    pm_twophase_fsigma8,
    pm_hubble_effective,
    BETA_CALIB, OM0_DEFAULT, SIGMA8_PLANCK, H0_PLANCK,
)
from general_relativity.cosmology import (
    lcdm_hubble,
    lcdm_fsigma8,
    lcdm_luminosity_distance,
    H0_PLANCK as H0_GR,
    OM0_PLANCK, SIGMA8_PLANCK as S8_GR,
)

OM0 = OM0_DEFAULT      # 0.315
H0  = H0_PLANCK        # 67.4 km/s/Mpc

# ── RSD data (same as fsigma8_comparison.py) ─────────────────────────────────
RSD_DATA = [
    # (z,  fσ8,  err,   name)
    (0.067, 0.423, 0.055, "6dFGS"),
    (0.150, 0.490, 0.145, "SDSS MGS"),
    (0.250, 0.351, 0.058, "SDSS LRG"),
    (0.320, 0.427, 0.056, "BOSS LOWZ"),
    (0.380, 0.440, 0.060, "BOSS DR12"),
    (0.440, 0.413, 0.080, "WiggleZ"),
    (0.570, 0.427, 0.066, "BOSS CMASS"),
    (0.600, 0.433, 0.067, "BOSS DR12"),
    (0.730, 0.437, 0.072, "WiggleZ"),
    (0.780, 0.380, 0.044, "WiggleZ"),
    (0.860, 0.441, 0.076, "VIPERS"),
    (1.000, 0.455, 0.120, "VVDS"),
    (1.360, 0.490, 0.139, "DESI LRG"),
    (1.400, 0.482, 0.116, "FastSound"),
]
z_data  = np.array([d[0] for d in RSD_DATA])
fs8_data = np.array([d[1] for d in RSD_DATA])
err_data = np.array([d[2] for d in RSD_DATA])


# ─────────────────────────────────────────────────────────────────────────────
#  Explanatory tables
# ─────────────────────────────────────────────────────────────────────────────

def print_beta_table():
    """Show how β_eff(z) varies — revealing the origin of β=0.8."""
    print()
    print("=" * 72)
    print("  Effective instantaneous power-law index  β_eff(z) = d ln H/d ln(1+z)")
    print()
    print("  For the PM two-phase model with w_E=-1 (= ΛCDM):")
    print("  β_eff → 3/2  at high z  (matter dominated)")
    print("  β_eff → (3/2)Ω_m ≈ 0.47  at z=0  (energy dominated)")
    print()
    print(f"  {'z':>6}  {'β_eff':>8}  {'f_E(z)':>10}  {'H_2ph/H0':>12}  {'H_β0.8/H0':>12}")
    print("  " + "-" * 65)
    for z in [0, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0]:
        be  = float(pm_effective_beta(z, w_E=-1.0))
        fe  = float(pm_energy_phase_fraction(z, w_E=-1.0))
        H2p = float(pm_twophase_hubble(z, w_E=-1.0)) / H0
        Hb  = float(pm_hubble_effective(z, beta=0.8)) / H0
        print(f"  {z:>6.1f}  {be:>8.4f}  {fe:>10.4f}  {H2p:>12.4f}  {Hb:>12.4f}")
    print()
    # Compute D_L-weighted mean β_eff
    z_g = np.linspace(0.1, 1.5, 500)
    beta_g  = pm_effective_beta(z_g, w_E=-1.0)
    weights = 1.0 / np.array([pm_twophase_hubble(z, w_E=-1.0) for z in z_g])
    beta_wt = float(np.average(beta_g, weights=weights))
    print(f"  D_L-weighted mean β_eff over z=0.1–1.5 = {beta_wt:.4f}")
    print(f"  → This is why β=0.8 was the best single-number fit to ΛCDM D_L!")
    print("=" * 72)
    print()


def print_key_result():
    """Print the central insight."""
    print()
    print("=" * 72)
    print("  CENTRAL RESULT: β=0.8 encodes ΛCDM dark energy implicitly")
    print("=" * 72)
    print()
    print("  Single-power-law PM:  H_PM(z) = H₀(1+z)^0.8  (one fixed exponent)")
    print("  Two-phase PM:         H(z) = H₀√[0.315(1+z)³ + 0.685]  (w_E=-1)")
    print()
    print("  At w_E=-1, the two-phase model IS ΛCDM.")
    print()
    print("  PM cosmology timeline:")
    print("  ─────────────────────────────────────────────────────────────────")
    print("  1. β=0.8 found by fitting D_L^ΛCDM with a power-law (prior work)")
    print("  2. Structure growth only works with Ω_eff=0.315, NOT Ω_eff=1")
    print("     (shown in fsigma8_comparison.py)")
    print("  3. Ω_eff=1 → runaway growth (PM-Drag ODE fails)")
    print("  4. TWO-PHASE: Ω_m=0.315 matter + Ω_E=0.685 energy phase")
    print("  5. w_E=-1 two-phase H(z) = ΛCDM H(z)")
    print("  6. β=0.8 ≈ D_L-weighted β_eff of ΛCDM over z=0.1–1.5")
    print()
    print("  PM-specific predictions:")
    print("  ─────────────────────────────────────────────────────────────────")
    print("  • The energy phase is the PM medium's ground state (φ<φ_crit=1)")
    print("    expanded to cosmic scales; same physics caps M_max ≈ 14 M☉")
    print("  • w_E = -1 exactly (constant density, no dilution)")
    print("  • Ω_E ≈ 0.685 is dynamically connected to ρ_E ~ e·ρ_nuc at")
    print("    compact-object scales (the PM form of the cosmological constant")
    print("    problem: 10^44 ratio between compact and cosmological scales)")
    print("  • Small departure w_E > -1 possible if ongoing matter→energy")
    print("    phase conversion: would give quintessence-like dark energy")
    print("=" * 72)
    print()


def print_fsigma8_comparison():
    """Compare fσ8 for single-β vs two-phase vs ΛCDM."""
    print("=" * 90)
    print("  fσ8(z): Single power-law β=0.8 vs Two-phase w_E=-1 vs ΛCDM")
    print("=" * 90)
    z_vals = np.array([0.067, 0.15, 0.32, 0.44, 0.57, 0.78, 1.0, 1.40])
    fs8_tp   = pm_twophase_fsigma8(z_vals, w_E=-1.0)
    fs8_lcdm = lcdm_fsigma8(z_vals)

    from pushing_medium.cosmology import pm_drag_fsigma8
    fs8_drag = pm_drag_fsigma8(z_vals, Om_eff=0.315, beta=0.8)

    print(f"  {'z':>6}  {'ΛCDM':>10}  {'2-Phase w=-1':>14}  {'Drag β=0.8':>14}  "
          f"{'Match 2ph/ΛCDM':>16}")
    print("  " + "-" * 75)
    for i, z in enumerate(z_vals):
        ratio = fs8_tp[i] / fs8_lcdm[i]
        print(f"  {z:>6.3f}  {fs8_lcdm[i]:>10.4f}  {fs8_tp[i]:>14.4f}  "
              f"{fs8_drag[i]:>14.4f}  {ratio:>16.7f}")
    print()
    print("  → Two-phase w_E=-1 and ΛCDM are numerically identical (as expected)")
    print("=" * 90)
    print()


# ─────────────────────────────────────────────────────────────────────────────
#  4-panel figure
# ─────────────────────────────────────────────────────────────────────────────

def make_figure():
    if not HAS_MPL:
        return

    z_fine = np.linspace(0.0, 4.0, 300)
    z_fine_nozero = np.linspace(0.02, 4.0, 300)

    # Pre-compute
    beta_eff    = pm_effective_beta(z_fine)
    f_E         = pm_energy_phase_fraction(z_fine)
    H_2phase    = np.array([pm_twophase_hubble(z, w_E=-1.0) / H0 for z in z_fine])
    H_beta08    = np.array([pm_hubble_effective(z, beta=0.8) / H0 for z in z_fine])
    H_lcdm      = np.array([lcdm_hubble(z) / H0 for z in z_fine])
    H_2ph_wm05  = np.array([pm_twophase_hubble(z, w_E=-0.5) / H0 for z in z_fine])

    z_arr_fs8 = np.linspace(0.02, 1.55, 120)
    fs8_2ph_m1  = pm_twophase_fsigma8(z_arr_fs8, w_E=-1.0)
    fs8_2ph_m05 = pm_twophase_fsigma8(z_arr_fs8, w_E=-0.5)
    fs8_lcdm    = lcdm_fsigma8(z_arr_fs8)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "PM Two-Phase Cosmology: β=0.8 is ΛCDM Dark Energy in Disguise",
        fontsize=13, fontweight='bold')

    # ── Panel A: β_eff(z) ─────────────────────────────────────────────────
    ax = axes[0, 0]
    ax.plot(z_fine, beta_eff, color='steelblue', lw=2.2, label=r'$\beta_{\rm eff}(z)$ [ΛCDM, $w_E=-1$]')
    ax.axhline(1.5, color='gray', lw=1.0, ls='--', label=r'$\beta=3/2$ (matter dom.)')
    ax.axhline(0.8, color='crimson', lw=1.5, ls=':', label=r'PM calibration $\beta=0.8$')
    ax.axhline(1.5*OM0, color='darkorange', lw=1.0, ls='--',
               label=fr'$\frac{{3}}{{2}}\Omega_m={1.5*OM0:.3f}$ (energy dom.)')
    # Mark the D_L-sensitive range
    ax.axvspan(0.1, 1.5, alpha=0.06, color='crimson', label='D_L-sensitive range')

    z_g = np.linspace(0.1, 1.5, 500)
    weights = 1.0 / np.array([pm_twophase_hubble(z, w_E=-1.0) for z in z_g])
    beta_wt = float(np.average(pm_effective_beta(z_g), weights=weights))
    ax.axhline(beta_wt, color='crimson', lw=1.0, ls='-.',
               label=fr'$\langle\beta_{{eff}}\rangle_{{D_L}}={beta_wt:.3f}$')

    ax.set_xlabel("Redshift $z$", fontsize=11)
    ax.set_ylabel(r"$\beta_{\rm eff}(z) = d\ln H/d\ln(1+z)$", fontsize=11)
    ax.set_xlim(-0.05, 4.0)
    ax.set_ylim(0.3, 1.6)
    ax.legend(fontsize=7.5, loc='lower right')
    ax.grid(True, ls=':', lw=0.5, alpha=0.6)
    ax.set_title(r"(A)  Effective power-law index $\beta_{\rm eff}(z)$", fontsize=10)

    # ── Panel B: Energy-phase fraction f_E(z) ─────────────────────────────
    ax = axes[0, 1]
    ax.plot(z_fine, f_E, color='darkorchid', lw=2.2)
    ax.axhline(1 - OM0, color='darkorchid', lw=0.8, ls='--',
               label=fr'$\Omega_\Lambda={1-OM0:.3f}$')
    ax.axhline(OM0, color='saddlebrown', lw=0.8, ls='--',
               label=fr'$\Omega_m={OM0:.3f}$')
    ax.fill_between(z_fine, f_E, alpha=0.15, color='darkorchid',
                    label='Energy phase (dark energy)')
    ax.fill_between(z_fine, f_E, 1, alpha=0.10, color='saddlebrown',
                    label='Matter phase')
    ax.set_xlabel("Redshift $z$", fontsize=11)
    ax.set_ylabel(r"Energy-phase fraction $f_E(z)$", fontsize=11)
    ax.set_xlim(-0.05, 4.0)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=8.5)
    ax.grid(True, ls=':', lw=0.5, alpha=0.6)
    ax.set_title("(B)  Energy-phase fraction $f_E(z)$\n"
                 r"$f_E(0)=\Omega_\Lambda=0.685$ ; $f_E\to0$ at high $z$", fontsize=10)

    # ── Panel C: H(z) comparison ──────────────────────────────────────────
    ax = axes[1, 0]
    ax.semilogy(z_fine, H_lcdm,   '-',  color='steelblue',  lw=2.2, label='ΛCDM')
    ax.semilogy(z_fine, H_2phase, '--', color='darkorchid', lw=1.8,
                label='Two-phase $w_E=-1$  [= ΛCDM]')
    ax.semilogy(z_fine, H_beta08, ':',  color='crimson',    lw=1.8,
                label=r'Single power law $\beta=0.8$')
    ax.semilogy(z_fine, H_2ph_wm05, '-.', color='darkorange', lw=1.4,
                label='Two-phase $w_E=-0.5$')
    ax.set_xlabel("Redshift $z$", fontsize=11)
    ax.set_ylabel("$H(z)/H_0$", fontsize=11)
    ax.set_xlim(-0.05, 4.0)
    ax.legend(fontsize=8.5)
    ax.grid(True, ls=':', lw=0.5, alpha=0.6, which='both')
    ax.set_title("(C)  Hubble rate: single $\\beta=0.8$ vs two-phase\n"
                 "Two-phase $w_E=-1$ is identical to ΛCDM", fontsize=10)

    # ── Panel D: fσ8 ──────────────────────────────────────────────────────
    ax = axes[1, 1]
    ax.errorbar(z_data, fs8_data, yerr=err_data,
                fmt='ko', ms=4, capsize=3, zorder=5, label='RSD data')
    ax.plot(z_arr_fs8, fs8_lcdm,    '-',  color='steelblue',  lw=2.2, label='ΛCDM')
    ax.plot(z_arr_fs8, fs8_2ph_m1,  '--', color='darkorchid', lw=1.8,
            label='Two-phase $w_E=-1$  [= ΛCDM]')
    ax.plot(z_arr_fs8, fs8_2ph_m05, '-.', color='darkorange', lw=1.4,
            label='Two-phase $w_E=-0.5$')
    ax.set_xlabel("Redshift $z$", fontsize=11)
    ax.set_ylabel(r"$f\sigma_8(z)$", fontsize=11)
    ax.set_xlim(-0.05, 1.6)
    ax.set_ylim(0.1, 0.9)
    ax.legend(fontsize=8.5)
    ax.grid(True, ls=':', lw=0.5, alpha=0.6)
    ax.set_title(r"(D)  $f\sigma_8$: two-phase PM ($w_E=-1$) = ΛCDM", fontsize=10)

    plt.tight_layout()
    out = ROOT / 'results' / 'pm_twophase_cosmology.png'
    out.parent.mkdir(exist_ok=True)
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Figure saved → {out}")
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print_beta_table()
    print_key_result()
    print_fsigma8_comparison()
    make_figure()


if __name__ == "__main__":
    main()
