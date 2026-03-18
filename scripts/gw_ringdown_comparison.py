#!/usr/bin/env python3
"""PM vs GR: GW ringdown, QNMs, and surface echo comparison.

Summary of predictions
-----------------------
GR (black holes):
  After a merger, the remnant BH emits GWs at the quasi-normal mode (QNM)
  frequency set by its horizon geometry:
    f_QNM = 0.37367 c³ / (2π GM)  ≈ 12.1 kHz × (M_sun/M)
    Q_QNM ≈ 2.10  (critically damped — energy absorbed by horizon)

PM (compact objects, no horizon):
  - No horizon →  no GR-like QNMs from trapped photon orbits
  - PM photon sphere r_ps = 2GM/c² < R_star for all realistic stars
  - Physical surface → GW reflection → echoes at τ_echo = 2R/c
  - Structural breathing modes at f_surf = c/(2√2 R)  ≈ 9–12 kHz
  - Q > GR (no energy sink comparable to a horizon)
  - MASS GAP: PM compact objects max out at ~2.4 M_sun →
    PM cannot form remnants like GW150914 (62 M_sun)

LIGO events
-----------
  GW150914: M_final ≈ 62 M_sun → GR f_QNM ≈ 195 Hz (in LIGO band)
                                  PM: NO such object possible → NO ringdown
                                  *** PM FALSIFICATION CHALLENGE ***

  GW170817: M_final ≈ 2.7 M_sun → GR f_QNM ≈ 4,500 Hz (near LIGO edge)
                                   PM: f_surf ≈ 8,200 Hz (if star forms)
                                   PM max mass ≈ 2.4 M_sun → marginally over limit
"""

import sys
import math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from general_relativity.classical import (
    gr_qnm_frequency, gr_qnm_damping_time, gr_photon_sphere_radius, G, c,
)
from pushing_medium.core import (
    pm_photon_sphere_radius, pm_surface_echo_delay,
    pm_surface_mode_frequency, pm_surface_mode_damping_time,
)
from pushing_medium.stellar_structure import compute_mr_curve, M_SUN

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ── Selected LIGO events ─────────────────────────────────────────────────────
LIGO_EVENTS = [
    # name,           M_final [M_sun],  notes
    ("GW150914",      62.0,             "BNS/BBH — far above PM max mass"),
    ("GW170104",      49.1,             "BBH — far above PM max mass"),
    ("GW170814",      53.4,             "BBH — far above PM max mass"),
    ("GW170817",       2.74,            "BNS merger — marginally over PM max"),
    ("GW190521",     150.0,             "IMR BBH — extreme, PM has no description"),
]

# PM max compact object mass: stability cap ρ_c ≤ ρ_crit = e×ρ_nuc.
# PM uses Newtonian structure (no TOV), so mass grows monotonically with ρ_c;
# there is no GR-style turnaround. The max mass at ρ_c = ρ_crit is ~13-14 M_sun.
_PM_MAX_MASS_MSUN = 13.5   # approximate; from compute_mr_curve at ρ_c = ρ_crit

# GR BH QNM quality factor
_Q_GR = 0.37367 / (2 * 0.08896)  # ≈ 2.10


def print_ligo_table():
    print("\n" + "=" * 90)
    print("  LIGO Merger Remnants: GR QNM vs PM Prediction")
    print("=" * 90)
    print(f"  {'Event':<14} {'M_fin [M☉]':<13} {'f_QNM_GR [Hz]':<16} "
          f"{'τ_damp_GR [ms]':<17} {'PM status'}")
    print("  " + "-" * 86)
    for name, M_msun, note in LIGO_EVENTS:
        f = gr_qnm_frequency(M_msun * M_SUN)
        tau = gr_qnm_damping_time(M_msun * M_SUN) * 1e3  # ms
        pm_status = "NO COMPACT OBJECT (progenitor > PM max mass ~13 M☉)" if M_msun > _PM_MAX_MASS_MSUN else "within PM range"
        print(f"  {name:<14} {M_msun:<13.1f} {f:<16.1f} {tau:<17.3f} {pm_status}")
    print()
    print(f"  GR BH QNM quality factor Q ≈ {_Q_GR:.2f}  (critically damped from horizon)")
    print("=" * 90 + "\n")


def print_pm_star_modes():
    """Print PM surface mode frequencies along the M-R track."""
    print("=" * 72)
    print("  PM Compact Stars: Surface Modes and Echo Delays")
    print("=" * 72)
    print(f"  {'M [M☉]':<10} {'R [km]':<10} {'f_surf [kHz]':<15} "
          f"{'τ_echo [μs]':<14} {'f_GR_QNM [kHz]':<15} {'ratio f_PM/f_GR'}")
    print("  " + "-" * 68)

    _, M_arr, R_arr = compute_mr_curve(n_points=20)
    mask = np.isfinite(M_arr) & np.isfinite(R_arr)
    M_arr, R_arr = M_arr[mask], R_arr[mask]

    # Pick a representative subset
    indices = np.linspace(0, len(M_arr) - 1, 8, dtype=int)
    for i in indices:
        M_msun, R_km = M_arr[i], R_arr[i]
        R_m = R_km * 1e3
        M_kg = M_msun * M_SUN
        f_surf = pm_surface_mode_frequency(R_m) / 1e3        # kHz
        tau_echo = pm_surface_echo_delay(R_m) * 1e6          # μs
        f_gr = gr_qnm_frequency(M_kg) / 1e3                   # kHz
        ratio = f_surf / f_gr if f_gr > 0 else float('nan')
        print(f"  {M_msun:<10.3f} {R_km:<10.2f} {f_surf:<15.2f} "
              f"{tau_echo:<14.1f} {f_gr:<15.2f} {ratio:.2f}×")
    print()
    print(f"  PM sound speed: c_s = c/√2 ≈ {c/math.sqrt(2)/1e3:.0f} km/s")
    print(f"  PM surface mode Q_factor (no horizon) >> GR QNM Q ≈ {_Q_GR:.2f}")
    print("=" * 72 + "\n")


def print_photon_sphere_comparison():
    print("=" * 60)
    print("  Photon Sphere: PM vs GR")
    print("=" * 60)
    print(f"  {'M [M☉]':<12} {'r_ps_PM [km]':<16} {'r_ps_GR [km]':<16} {'ratio'}")
    print("  " + "-" * 56)
    for M_msun in [1.0, 1.4, 2.0, 2.4, 10.0, 60.0]:
        M_kg = M_msun * M_SUN
        r_pm = pm_photon_sphere_radius(M_kg) / 1e3
        r_gr = gr_photon_sphere_radius(M_kg) / 1e3
        print(f"  {M_msun:<12.1f} {r_pm:<16.3f} {r_gr:<16.3f} {r_pm/r_gr:.3f}×")
    print()
    print("  PM r_ps = 2GM/c²  vs  GR r_ps = 3GM/c²")
    print("  For typical PM stars (R > 8 km), r_ps is inside the star →")
    print("  no exterior photon trapping → no QNMs in PM")
    print("=" * 60 + "\n")


def make_plots():
    if not HAS_MPL:
        print("matplotlib not available — skipping plots.")
        return

    _, M_arr, R_arr = compute_mr_curve(n_points=50)
    mask = np.isfinite(M_arr) & np.isfinite(R_arr)
    M_arr, R_arr = M_arr[mask], R_arr[mask]
    M_kg = M_arr * M_SUN
    R_m  = R_arr * 1e3

    f_surf_pm = np.array([pm_surface_mode_frequency(r) for r in R_m]) / 1e3   # kHz
    tau_echo  = np.array([pm_surface_echo_delay(r) * 1e6 for r in R_m])       # μs
    f_qnm_gr  = np.array([gr_qnm_frequency(m) for m in M_kg]) / 1e3           # kHz

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle("PM vs GR: GW Ringdown and Surface Modes", fontsize=13)

    # ── Left: Frequency comparison ──────────────────────────────────────────
    ax = axes[0]
    ax.plot(M_arr, f_surf_pm, 'r-',  lw=2.0, label='PM surface mode  $f = c_s/(2R)$')
    ax.plot(M_arr, f_qnm_gr,  'b--', lw=2.0, label='GR QNM  $f = 0.37367\\,c^3/(2\\pi GM)$')

    # LIGO events GR QNM
    for name, M_msun, _ in LIGO_EVENTS:
        if M_msun < 5:
            f_ev = gr_qnm_frequency(M_msun * M_SUN) / 1e3
            ax.axhline(f_ev, ls=':', lw=0.8, color='blue', alpha=0.5)
            ax.text(0.1, f_ev + 0.15, f'{name} (GR)', fontsize=7, color='blue')

    ax.set_xlabel("Mass M [M☉]",    fontsize=11)
    ax.set_ylabel("Frequency [kHz]", fontsize=11)
    ax.set_ylim(0, 20)
    ax.set_xlim(0, max(M_arr) * 1.05)
    ax.legend(fontsize=9)
    ax.grid(True, ls=':', lw=0.5, alpha=0.7)
    ax.set_title("Ring Frequency: PM surface mode vs GR QNM\n"
                 "(PM ~2× higher frequency; qualitatively different)", fontsize=10)

    # ── Right: Echo delay and damping ───────────────────────────────────────
    ax2 = axes[1]
    ax2.plot(M_arr, tau_echo, 'r-', lw=2.0, label='PM echo delay  $2R/c$ [μs]')

    # GR QNM damping time (μs)
    tau_gr_us = np.array([gr_qnm_damping_time(m) * 1e6 for m in M_kg])
    ax2.plot(M_arr, tau_gr_us, 'b--', lw=2.0,
             label='GR QNM damping time $\\tau_{\\rm QNM}$ [μs]')

    ax2.set_xlabel("Mass M [M☉]", fontsize=11)
    ax2.set_ylabel("Time [μs]",    fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, ls=':', lw=0.5, alpha=0.7)
    ax2.set_title("PM Echo Delay vs GR Damping Time\n"
                  "(PM echo short; GR damps fast from horizon energy loss)", fontsize=10)

    plt.tight_layout()
    out = Path(__file__).parent.parent / 'results' / 'gw_ringdown_comparison.png'
    out.parent.mkdir(exist_ok=True)
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {out}")
    plt.close()


def main():
    print_photon_sphere_comparison()
    print_ligo_table()
    print_pm_star_modes()
    make_plots()

    print("KEY FALSIFICATION CONCLUSIONS")
    print("-" * 60)
    print("1. GW150914 component masses (~29 + 36 M_sun): both exceed the")
    print("   PM stability cap max mass of ~13-14 M_sun. PM cannot form")
    print("   these objects. Yet LIGO observes them as compact (from chirp")
    print("   mass scaling). This is a direct PM falsification challenge.")
    print()
    print("2. GW170817 (1.46 + 1.27 M_sun components, 2.74 M_sun remnant):")
    print("   Both components AND remnant are within PM's mass range.")
    print("   HOWEVER: PM predicts R ~ 15-17 km at these masses while NICER")
    print("   measures R ~ 11-13 km. PM over-predicts radii by ~30%.")
    print()
    print("3. For PM stars within its mass range (< 13 M_sun), the ringdown")
    print("   spectrum is surface modes at ~5-12 kHz (f_surf = c/(2\u221a2 R))")
    print("   vs GR QNMs at ~1-12 kHz (f_QNM \u221d 1/M). Near M ~ 1.5 M_sun these")
    print("   happen to cross; at M > 2 M_sun PM gives consistently HIGHER")
    print("   frequencies. Distinguishable with Einstein Telescope.")


if __name__ == "__main__":
    main()
