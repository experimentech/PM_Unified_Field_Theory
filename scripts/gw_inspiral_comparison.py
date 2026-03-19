"""
GW Inspiral Comparison: Pushing Medium vs General Relativity.

At leading (0PN) order, PM and GR produce identical chirp waveforms
because:
  (1) Quadrupole power P_GW is identical (proven by Hulse-Taylor agreement)
  (2) Orbital energy E = -GM/(2a) is identical (Newtonian limit)

The dominant LIGO constraint on PM is the MASS GAP:
  • All confident LIGO BBH detections in GWTC have component masses
    far exceeding PM M_max ≈ 13–14 M☉
  • Only BNS events (e.g. GW170817) have components within PM mass range,
    but the PM neutron-star radius is ~30% too large, predicting
    Λ_tidal ~ 4× larger than GW170817 observational upper bound (~900).

Usage:
    python scripts/gw_inspiral_comparison.py
"""

import math
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pushing_medium import (
    pm_chirp_mass,
    pm_gw_frequency_deriv,
    pm_time_to_coalescence,
    pm_gw_strain_amplitude,
    pm_gw_chirp_waveform,
)
from general_relativity import gr_chirp_mass

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
G = 6.674e-11        # m³ kg⁻¹ s⁻²
c = 2.998e8          # m s⁻¹
M_sun = 1.989e30     # kg
Mpc = 3.086e22       # m

PM_M_MAX = 14.0 * M_sun   # conservative PM compact-object mass limit

# ---------------------------------------------------------------------------
# GWTC catalog events (component masses from public GWTC-3 release)
# ---------------------------------------------------------------------------
GWTC_EVENTS = [
    # name,       M1/M☉,  M2/M☉,  D_L/Mpc,   notes
    ("GW150914",  35.6,   30.6,   410,  "BBH — both >> PM M_max"),
    ("GW151226",  13.7,    7.7,   440,  "BBH — primary ≈ PM M_max"),
    ("GW170814",  30.5,   25.3,   540,  "BBH — both >> PM M_max"),
    ("GW170817",   1.46,   1.27,   40,  "BNS — both within PM range"),
    ("GW190814",  22.2,    2.6,   241,  "NSBH — primary >> PM M_max"),
    ("GW190521",  85.0,   66.0,  3931,  "IMBH — far above PM M_max"),
]


def print_table():
    print("\n" + "="*80)
    print("GWTC EVENT TABLE — Pushing Medium vs GR")
    print("="*80)
    header = (
        f"{'Name':<12} {'M1/M☉':>7} {'M2/M☉':>7} {'Mc/M☉':>7} "
        f"{'PM ok?':>8} {'PM/GR 0PN':>12} {'Notes'}"
    )
    print(header)
    print("-"*80)
    for name, m1, m2, dL, notes in GWTC_EVENTS:
        M1 = m1 * M_sun
        M2 = m2 * M_sun
        M_c_pm = pm_chirp_mass(M1, M2)
        M_c_gr = gr_chirp_mass(M1, M2)
        pm_ok = (M1 < PM_M_MAX) and (M2 < PM_M_MAX)
        pm_ok_str = "YES" if pm_ok else "NO (gap)"
        ratio = M_c_pm / M_c_gr   # should be 1.000 at 0PN
        print(
            f"{name:<12} {m1:>7.1f} {m2:>7.1f} "
            f"{M_c_pm/M_sun:>7.3f} {pm_ok_str:>8} {ratio:>12.6f} {notes}"
        )
    print("="*80)
    print("PM/GR 0PN column is exactly 1.000000 (waveforms indistinguishable at 0PN).\n")


def print_tidal_note():
    print("GW170817 tidal deformability note:")
    print("  GW170817 observed: Λ̃ ≤ 900  (90% credible, GW measurement)")
    print("  GR canonical NS (R≈11 km): Λ ≈ 300–600")
    print("  PM NS (R≈30% larger):       Λ_PM ≈ (1.3)^5 × Λ_GR ≈ 3.7× Λ_GR")
    print("  PM prediction: Λ_PM ~ 1100–2200  →  potentially inconsistent")
    print("  (requires full PM stellar-structure integration for definitive result)\n")


def make_figure():
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle(
        "GW Inspiral: Pushing Medium vs General Relativity\n"
        "(0PN — PM and GR are identical at this order)",
        fontsize=13, fontweight='bold', y=0.98
    )

    # ------------------------------------------------------------------
    # Panel A: Frequency chirp f(t) during last 30 s before merger
    # ------------------------------------------------------------------
    ax = axes[0, 0]
    systems = [
        ("GW170817-like", 1.4 * M_sun, 1.4 * M_sun, "steelblue",  True),
        ("GW151226-like",13.7 * M_sun, 7.7 * M_sun,  "darkorange", False),
        ("GW150914-like",35.6 * M_sun,30.6 * M_sun,  "firebrick",  False),
    ]
    tau_max = 30.0

    for label, M1, M2, color, pm_ok in systems:
        M_c = pm_chirp_mass(M1, M2)
        tau_stop = pm_time_to_coalescence(1000.0, M_c)
        tau_arr = np.linspace(tau_max, tau_stop, 4096)
        t_arr = tau_max - tau_arr
        f_arr = (5.0**(3.0/8.0) / (8.0*math.pi)) * (c**3/(G*M_c))**(5.0/8.0) * tau_arr**(-3.0/8.0)
        ls = '-' if pm_ok else '--'
        lw = 2.0 if pm_ok else 1.5
        ax.plot(t_arr, f_arr, color=color, lw=lw, ls=ls, label=label)

    ax.set_xlabel("Time before merger [s]", fontsize=10)
    ax.set_ylabel("GW frequency [Hz]", fontsize=10)
    ax.set_title("(A) Frequency Chirp  f (t)", fontsize=11)
    ax.set_xlim(0, tau_max)
    ax.set_ylim(0, 1050)
    ax.axhspan(10, 1000, alpha=0.06, color='green', label="LIGO band (10–1000 Hz)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    pm_note = mpatches.Patch(color='none', label="Solid = PM accessible\nDashed = PM mass-gap")
    ax.legend(fontsize=8, handles=ax.get_legend_handles_labels()[0] + [pm_note])

    # ------------------------------------------------------------------
    # Panel B: h+(t) waveform for GW170817-like BNS (last 2 s)
    # ------------------------------------------------------------------
    ax = axes[0, 1]
    M1_bns = M2_bns = 1.4 * M_sun
    D_bns = 40.0 * Mpc
    tau_bns = 2.0
    t_bns, f_bns, hp_bns, hx_bns = pm_gw_chirp_waveform(
        M1_bns, M2_bns, D_bns, tau_bns, N=8192, iota=0.0
    )
    ax.plot(t_bns, hp_bns * 1e21, color='steelblue', lw=0.8, label=r"$h_+$ (PM = GR at 0PN)")
    ax.set_xlabel("Time before merger [s]", fontsize=10)
    ax.set_ylabel(r"Strain $h_+$ [$\times 10^{-21}$]", fontsize=10)
    ax.set_title("(B) GW170817-like BNS Waveform (D = 40 Mpc)", fontsize=11)
    ax.set_xlim(0, tau_bns)
    ax.legend(fontsize=9)
    ax.text(
        0.02, 0.97,
        "PM accessible (M₁, M₂ < PM M_max)\nWaveform identical to GR at 0PN\nDifference: tidal deformability Λ",
        transform=ax.transAxes, fontsize=8, va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8)
    )
    ax.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # Panel C: Time to coalescence vs frequency
    # ------------------------------------------------------------------
    ax = axes[1, 0]
    f_arr_log = np.logspace(1.0, 3.0, 300)   # 10 Hz to 1000 Hz
    for label, M1, M2, color, pm_ok in systems:
        M_c = pm_chirp_mass(M1, M2)
        tc = np.array([pm_time_to_coalescence(f, M_c) for f in f_arr_log])
        ls = '-' if pm_ok else '--'
        lw = 2.0 if pm_ok else 1.5
        ax.loglog(f_arr_log, tc, color=color, lw=lw, ls=ls, label=label)
    ax.axvline(10.0, color='gray', ls=':', lw=1.2, label="LIGO entry ~10 Hz")
    ax.set_xlabel("GW frequency [Hz]", fontsize=10)
    ax.set_ylabel("Time to coalescence [s]", fontsize=10)
    ax.set_title(r"(C) Merger Time  $t_c(f)$", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, which='both', alpha=0.3)

    # ------------------------------------------------------------------
    # Panel D: GWTC component masses vs PM M_max (mass-gap visualisation)
    # ------------------------------------------------------------------
    ax = axes[1, 1]
    event_names = [e[0] for e in GWTC_EVENTS]
    m1_vals = [e[1] for e in GWTC_EVENTS]
    m2_vals = [e[2] for e in GWTC_EVENTS]
    x = np.arange(len(event_names))
    width = 0.35
    pm_mmax_solar = PM_M_MAX / M_sun

    for i, (m1, m2) in enumerate(zip(m1_vals, m2_vals)):
        c1 = 'firebrick' if m1 > pm_mmax_solar else 'steelblue'
        c2 = 'firebrick' if m2 > pm_mmax_solar else 'steelblue'
        ax.bar(x[i] - width/2, m1, width, color=c1, alpha=0.8)
        ax.bar(x[i] + width/2, m2, width, color=c2, alpha=0.8)

    ax.axhline(pm_mmax_solar, color='black', ls='--', lw=2.0, label=f"PM $M_{{\\rm max}}$ ≈ {pm_mmax_solar:.0f} M☉")
    ax.set_xticks(x)
    ax.set_xticklabels([n.replace("GW", "GW\n") for n in event_names], fontsize=8)
    ax.set_ylabel("Component mass [$M_☉$]", fontsize=10)
    ax.set_title("(D) GWTC Masses vs PM Mass Limit", fontsize=11)
    ax.legend(fontsize=9)
    red_patch = mpatches.Patch(color='firebrick', alpha=0.8, label='Above PM $M_{\\rm max}$ (mass gap)')
    blue_patch = mpatches.Patch(color='steelblue', alpha=0.8, label='Within PM range')
    ax.legend(handles=[red_patch, blue_patch], fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = os.path.join(os.path.dirname(__file__), '..', 'results', 'gw_inspiral_comparison.png')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {os.path.abspath(out_path)}")
    plt.close()


def main():
    print_table()
    print_tidal_note()

    # Numerical summary
    print("Numerical checks")
    print("-"*50)
    M_c_170817 = pm_chirp_mass(1.46 * M_sun, 1.27 * M_sun)
    tc_20hz = pm_time_to_coalescence(20.0, M_c_170817)
    tc_100hz = pm_time_to_coalescence(100.0, M_c_170817)
    print(f"GW170817 chirp mass:          {M_c_170817/M_sun:.4f} M☉")
    print(f"GW170817 t_coal at f=20 Hz:   {tc_20hz:.1f} s")
    print(f"GW170817 t_coal at f=100 Hz:  {tc_100hz:.2f} s")
    dfdt_170817_50hz = pm_gw_frequency_deriv(50.0, M_c_170817)
    print(f"GW170817 df/dt at f=50 Hz:    {dfdt_170817_50hz:.5f} Hz/s")

    M_c_150914 = pm_chirp_mass(35.6 * M_sun, 30.6 * M_sun)
    dfdt_150914_50hz = pm_gw_frequency_deriv(50.0, M_c_150914)
    h_150914 = pm_gw_strain_amplitude(M_c_150914, 410.0 * Mpc, 150.0)
    print(f"\nGW150914 chirp mass:          {M_c_150914/M_sun:.2f} M☉")
    print(f"GW150914 df/dt at f=50 Hz:    {dfdt_150914_50hz:.2f} Hz/s")
    print(f"GW150914 strain at 150 Hz:    {h_150914:.2e}  (LIGO peak ~10⁻²¹)")

    print(f"\nPM M_max:                     {PM_M_MAX/M_sun:.1f} M☉")
    print(f"GW150914 M1/PM_Mmax:          {35.6*M_sun/PM_M_MAX:.1f}×  → MASS GAP")
    print(f"GW190814 M1/PM_Mmax:          {22.2*M_sun/PM_M_MAX:.1f}×  → MASS GAP")
    print(f"GW170817 M1/PM_Mmax:          {1.46*M_sun/PM_M_MAX:.2f}×  → accessible")
    print()

    make_figure()


if __name__ == "__main__":
    main()
