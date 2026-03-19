#!/usr/bin/env python3
"""
EHT Shadow Comparison: PM vs GR
================================
Compares PM compact-star shadow predictions with GR Schwarzschild BH shadows
and the EHT observations of M87* and Sgr A*.

Key results printed:
  1. EHT data vs GR vs PM (maximum-compactness case) angular diameters
  2. PM/GR shadow ratio as a function of stellar compactness k = R_star/R_s
  3. Mass-gap falsification: EHT source masses >> PM max compact-object mass
  4. Shadow-size comparison for PM-allowed compact stars vs GR BH equivalent

Figure saved to results/eht_shadow_comparison.png (4 panels).
"""

import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.general_relativity.classical import (
    gr_shadow_radius,
    gr_compact_star_shadow_radius,
)
from src.pushing_medium.core import pm_shadow_radius, G, c

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
M_SUN = 1.989e30
RAD_TO_UAS = 180.0 / math.pi * 3600.0 * 1.0e6   # rad → μas

# EHT sources
EHT_SOURCES = {
    "M87*":   dict(M_sun=6.5e9,  D_Mpc=16.8,  theta_uas=42.0,  err_uas=3.0),
    "Sgr A*": dict(M_sun=4.15e6, D_kpc=8.178,  theta_uas=51.8,  err_uas=2.3),
}

PC  = 3.086e16   # m
KPC = 1e3 * PC
MPC = 1e6 * PC

def source_distance(name, s):
    if "D_Mpc" in s:
        return s["D_Mpc"] * MPC
    return s["D_kpc"] * KPC

def shadow_uas(b, D):
    return 2.0 * b / D * RAD_TO_UAS

PM_M_MAX = 14.0   # M_sun, conservative upper bound on PM compact-object mass


# ---------------------------------------------------------------------------
# 1. Print EHT summary table
# ---------------------------------------------------------------------------
def print_eht_table():
    print("=" * 72)
    print("EHT Shadow: Observed vs GR Prediction vs PM (max-compact) Prediction")
    print("=" * 72)
    hdr = f"{'Source':<10} {'M/M_sun':>12} {'D':>10} {'θ_obs (μas)':>13} {'θ_GR':>8} {'θ_PM':>8} {'PM/GR':>7}"
    print(hdr)
    print("-" * 72)
    for name, s in EHT_SOURCES.items():
        M   = s["M_sun"] * M_SUN
        D   = source_distance(name, s)
        R_s = 2.0 * G * M / (c * c)
        b_gr   = gr_shadow_radius(M)
        b_pm   = pm_shadow_radius(M, R_s)   # most compact PM object
        theta_gr  = shadow_uas(b_gr, D)
        theta_pm  = shadow_uas(b_pm, D)
        dist_str = f"{s.get('D_Mpc', s.get('D_kpc'))} {'Mpc' if 'D_Mpc' in s else 'kpc'}"
        obs_str  = f"{s['theta_uas']:.1f} ± {s['err_uas']:.1f}"
        print(f"{name:<10} {s['M_sun']:>12.3e} {dist_str:>10} {obs_str:>13} {theta_gr:>8.1f} {theta_pm:>8.1f} {theta_pm/theta_gr:>7.3f}")
    print("-" * 72)
    print("Note: θ_PM is for the MOST COMPACT possible PM object (R_star = R_s).")
    print(f"      PM/GR ratio = 2e/(3√3) ≈ {2*math.e/(3*math.sqrt(3)):.4f} at max compactness.")
    print()
    print("MASS-GAP FALSIFICATION:")
    pm_max_kg = PM_M_MAX * M_SUN
    for name, s in EHT_SOURCES.items():
        ratio = s["M_sun"] / PM_M_MAX
        print(f"  {name}: {s['M_sun']:.2e} M_sun  =  {ratio:.1e} × PM M_max ({PM_M_MAX:.0f} M_sun)")
    print()


# ---------------------------------------------------------------------------
# 2. Shadow ratio b_PM / b_GR_BH vs compactness k = R_star / R_s
# ---------------------------------------------------------------------------
def print_compactness_table():
    print("=" * 65)
    print("PM Shadow vs GR BH Shadow as a Function of Compactness")
    print(f"  (b_GR_BH = 3√3/2 × R_s ≈ {1.5*math.sqrt(3):.4f} R_s)")
    print("=" * 65)
    print(f"  {'k':>6}  {'b_PM/R_s':>10}  {'b_GR_star/R_s':>14}  {'b_PM/b_GR_BH':>14}  {'b_PM/b_GR_star':>15}")
    print("-" * 65)
    M = 1.4 * M_SUN
    R_s = 2.0 * G * M / (c * c)
    b_gr_bh = gr_shadow_radius(M)
    ks = [1.01, 1.1, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]
    for k in ks:
        R_star = k * R_s
        b_pm   = pm_shadow_radius(M, R_star) / R_s
        b_grst = gr_compact_star_shadow_radius(M, R_star) / R_s
        ratio_bh   = b_pm * R_s / b_gr_bh
        ratio_star = b_pm / b_grst
        print(f"  {k:>6.2f}  {b_pm:>10.4f}  {b_grst:>14.4f}  {ratio_bh:>14.4f}  {ratio_star:>15.4f}")
    print("-" * 65)
    print(f"  k→1  {math.e:>10.4f}  {'→ ∞':>14}  {2*math.e/(3*math.sqrt(3)):>14.4f}  {'→ 0':>15}")
    print()
    print("Key: b_PM/b_GR_BH ≥ 1.046 for ALL R_star ≥ R_s")
    print("     PM compact stars always cast larger shadows than GR BHs of the same mass.")
    print()


# ---------------------------------------------------------------------------
# 3. PM vs GR shadow for PM-allowed masses (≤ 14 M_sun)
# ---------------------------------------------------------------------------
def print_pm_allowed_table():
    print("=" * 60)
    print("PM-Allowed Compact Stars: Shadow vs GR BH (typical NS)")
    print(f"  R_star = 3 R_s (k=3, typical PM neutron-star compactness)")
    print("=" * 60)
    print(f"  {'M/M_sun':>10}  {'b_PM (km)':>12}  {'b_GR_BH (km)':>14}  {'ratio':>8}")
    print("-" * 60)
    for M_sun in [1.0, 1.4, 2.0, 5.0, 10.0, 14.0]:
        M = M_sun * M_SUN
        R_s = 2.0 * G * M / (c * c)
        R_star = 3.0 * R_s   # typical PM NS compactness
        b_pm = pm_shadow_radius(M, R_star)
        b_gr = gr_shadow_radius(M)
        print(f"  {M_sun:>10.1f}  {b_pm/1e3:>12.2f}  {b_gr/1e3:>14.2f}  {b_pm/b_gr:>8.4f}")
    print("-" * 60)
    print()


# ---------------------------------------------------------------------------
# 4. Figure
# ---------------------------------------------------------------------------
def make_figure():
    M = 1.4 * M_SUN
    R_s = 2.0 * G * M / (c * c)

    # Panel A: b_PM/R_s and b_GR_star/R_s vs k (linear scale, k ∈ [1.05,6])
    k_vals = np.linspace(1.05, 6.0, 400)
    R_star_vals = k_vals * R_s
    b_pm_vals   = np.array([pm_shadow_radius(M, R) / R_s for R in R_star_vals])
    b_gr_vals   = np.array([gr_compact_star_shadow_radius(M, R) / R_s for R in R_star_vals])
    b_gr_bh     = gr_shadow_radius(M) / R_s   # horizontal line ≈ 2.598
    b_pm_min    = math.e                        # = 2.718

    # Panel B: PM/GR_BH ratio vs k
    ratio_vals = b_pm_vals / b_gr_bh

    # Panel C: Angular diameter vs mass for M ∈ [1,14] M_sun at k=3 and GR BH
    M_arr = np.linspace(1.0, 14.0, 200) * M_SUN
    D_fixed = 1.0e19   # 1 kpc in m (arbitrary reference distance)
    theta_pm_k3 = np.array([
        shadow_uas(pm_shadow_radius(M_, 3.0 * 2.0 * G * M_ / (c*c)), D_fixed)
        for M_ in M_arr
    ])
    theta_gr_bh = np.array([shadow_uas(gr_shadow_radius(M_), D_fixed) for M_ in M_arr])

    # Panel D: EHT angular-diameter comparison bar chart
    sources = list(EHT_SOURCES.keys())
    theta_obs  = [EHT_SOURCES[s]["theta_uas"] for s in sources]
    theta_errs = [EHT_SOURCES[s]["err_uas"]   for s in sources]
    theta_gr_eht  = []
    theta_pm_eht  = []
    for s in sources:
        src = EHT_SOURCES[s]
        M_  = src["M_sun"] * M_SUN
        D_  = source_distance(s, src)
        R_s_ = 2.0 * G * M_ / (c * c)
        theta_gr_eht.append(shadow_uas(gr_shadow_radius(M_), D_))
        theta_pm_eht.append(shadow_uas(pm_shadow_radius(M_, R_s_), D_))   # max-compact

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("PM vs GR: Black-Hole Shadow Comparison", fontsize=14, fontweight="bold")

    # --- Panel A: shadow impact parameter vs compactness ---
    ax = axes[0, 0]
    ax.plot(k_vals, b_pm_vals, "b-", lw=2, label=r"PM  $b/R_s = e^{R_s/R_\star} \cdot k$")
    ax.plot(k_vals, b_gr_vals, "r--", lw=2, label=r"GR star  $b/R_s = k/\sqrt{1-1/k}$")
    ax.axhline(b_gr_bh, color="red", lw=1.5, ls=":", label=rf"GR BH  $b/R_s = 3\sqrt{{3}}/2 \approx {b_gr_bh:.3f}$")
    ax.axhline(b_pm_min, color="blue", lw=1.5, ls=":", label=rf"PM min  $b/R_s = e \approx {b_pm_min:.3f}$")
    ax.set_xlabel(r"Compactness $k = R_\star / R_s$")
    ax.set_ylabel(r"Shadow impact parameter  $b / R_s$")
    ax.set_xlim(1.05, 6)
    ax.set_ylim(0, 12)
    ax.legend(fontsize=8)
    ax.set_title("(A)  Shadow size vs compactness")
    ax.grid(alpha=0.3)

    # --- Panel B: PM/GR_BH ratio ---
    ax = axes[0, 1]
    ax.plot(k_vals, ratio_vals, "b-", lw=2, label=r"$b_\mathrm{PM} / b_\mathrm{GR\,BH}$")
    ax.axhline(1.0, color="red", lw=1.2, ls="--", label="GR BH = 1")
    ax.axhline(2.0 * math.e / (3.0 * math.sqrt(3.0)), color="blue", lw=1, ls=":",
               label=rf"$2e/(3\sqrt{{3}}) \approx {2*math.e/(3*math.sqrt(3)):.3f}$  (k→1 limit)")
    ax.set_xlabel(r"Compactness $k = R_\star / R_s$")
    ax.set_ylabel(r"$b_\mathrm{PM} / b_\mathrm{GR\,BH}$")
    ax.set_xlim(1.05, 6)
    ax.set_ylim(0.9, 2.5)
    ax.legend(fontsize=9)
    ax.set_title("(B)  PM shadow / GR BH shadow ratio")
    ax.grid(alpha=0.3)

    # --- Panel C: Angular diameter vs mass ---
    ax = axes[1, 0]
    M_sun_arr = M_arr / M_SUN
    ax.plot(M_sun_arr, theta_pm_k3, "b-", lw=2, label=r"PM star ($k=3$, 1 kpc)")
    ax.plot(M_sun_arr, theta_gr_bh, "r--", lw=2, label=r"GR BH (1 kpc)")
    ax.fill_between(M_sun_arr, theta_gr_bh, theta_pm_k3, alpha=0.15, color="blue",
                    label="PM excess")
    ax.set_xlabel(r"Mass  $[M_\odot]$")
    ax.set_ylabel(r"Shadow angular diameter  $[\mu\mathrm{as}]$ at 1 kpc")
    ax.legend(fontsize=9)
    ax.set_title("(C)  Predicted shadow for PM-allowed masses")
    ax.grid(alpha=0.3)

    # --- Panel D: EHT comparison ---
    ax = axes[1, 1]
    x = np.arange(len(sources))
    width = 0.25
    ax.bar(x - width, theta_obs, width, label="EHT observed", color="black", alpha=0.7)
    ax.bar(x,         theta_gr_eht, width, label="GR prediction", color="red", alpha=0.7)
    ax.bar(x + width, theta_pm_eht, width, label="PM max-compact", color="blue", alpha=0.7)
    ax.errorbar(x - width, theta_obs, yerr=theta_errs, fmt="none", color="white",
                capsize=4, lw=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(sources)
    ax.set_ylabel(r"Angular diameter  $[\mu\mathrm{as}]$")
    ax.legend(fontsize=9)
    ax.set_title("(D)  EHT angular diameters")
    # Annotate mass ratio
    for i, s in enumerate(sources):
        ratio = EHT_SOURCES[s]["M_sun"] / PM_M_MAX
        ax.text(i, max(theta_obs[i], theta_gr_eht[i], theta_pm_eht[i]) + 2,
                f"M/M$_{{max}}$\n= {ratio:.0e}", ha="center", fontsize=7, color="purple")
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "results", "eht_shadow_comparison.png")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Figure saved → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print_eht_table()
    print_compactness_table()
    print_pm_allowed_table()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    ratio_min = 2.0 * math.e / (3.0 * math.sqrt(3.0))
    print(f"  PM minimum shadow / GR BH shadow = 2e/(3√3) = {ratio_min:.6f}")
    print(f"  → PM always predicts shadows ≥ {(ratio_min-1)*100:.1f}% LARGER than GR BH")
    print()
    print("  BUT: EHT sources (M87*, Sgr A*) have masses:")
    for name, s in EHT_SOURCES.items():
        print(f"    {name}: {s['M_sun']:.2e} M_sun = {s['M_sun']/PM_M_MAX:.1e} × PM M_max")
    print("  → PM cannot form these objects at all (mass gap falsification)")
    print("  → Shadow-shape comparison is moot; MASS ITSELF falsifies PM here.")
    print()
    print("  For PM-accessible masses (≤ 14 M_sun) at k=3:")
    M = 1.4 * M_SUN
    R_s = 2.0 * G * M / (c * c)
    b_pm_k3 = pm_shadow_radius(M, 3.0 * R_s)
    b_gr_bh  = gr_shadow_radius(M)
    print(f"    1.4 M_sun NS: b_PM = {b_pm_k3/1e3:.1f} km, b_GR_BH = {b_gr_bh/1e3:.1f} km "
          f"(ratio {b_pm_k3/b_gr_bh:.3f})")
    make_figure()


if __name__ == "__main__":
    main()
