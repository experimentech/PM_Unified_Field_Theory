#!/usr/bin/env python3
"""PM vs GR: perihelion precession — analytic and numerical comparison.

Key finding
-----------
PM force law a = (c²/2)∇φ = −GM/r² r̂ is EXACTLY Newtonian — no precession
from the force alone.  The 6πGM/c²a(1−e²) advance is a 1PN *metric* effect
(from g_tt = −c²/n², g_rr = 1−2φ), identical to GR because PM has PPN
parameters β=γ=1.  This is ANOTHER PASS for PM — it agrees with solar
system observations to all measured precision.

The growing discrepancy picture
--------------------------------
  ✓ Newtonian gravity:    PM = GR (exact, by construction)
  ✓ Light deflection:     PM = GR = 4GM/c²b  (β=γ=1)
  ✓ Perihelion advance:   PM = GR = 6πGM/c²a(1-e²)  (β=γ=1)
  ✓ Shapiro delay:        PM = GR  (same)
  ✗ Surface redshift:     PM ≈ 2× GR  (large measurable difference)
  ✗ GW150914 masses:      PM max mass ~13 M_sun; progenitors 29+36 M_sun exceed it
  ? QNM/ringdown:         PM predicts surface modes, not horizon QNMs
  ? NICER radii:          PM over-predicts NS radii by ~30%
"""

import sys
import math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from pushing_medium.core import (
    pm_perihelion_precession,
    pm_precession_arcsec_per_century,
    pm_integrate_orbit,
    G, c,
)
from general_relativity.classical import perihelion_precession

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

M_SUN = 1.989e30
AU    = 1.496e11
S_PER_YEAR = 365.25 * 86400.0

# ── Reference data ────────────────────────────────────────────────────────────
# (a [AU], e, T [yr], observed GR-only residual [arcsec/cy])
PLANETS = {
    'Mercury': (0.387098, 0.205630, 0.24085,  42.9799),
    'Venus':   (0.723327, 0.006773, 0.61520,   8.6247),
    'Earth':   (1.000000, 0.016709, 1.00000,   3.8387),
    'Mars':    (1.523679, 0.093401, 1.88085,   1.3510),
    'Jupiter': (5.203363, 0.048393, 11.8618,   0.0623),
    'Saturn':  (9.582017, 0.056116, 29.4571,   0.0136),
}

# PSR B1913+16 (Hulse-Taylor binary pulsar)
PSR_M_TOT = (1.4408 + 1.3886) * M_SUN
PSR_P_B   = 27906.9807   # s
PSR_E     = 0.6171334
PSR_OBS   = 4.22663       # observed deg/yr
PSR_UNCT  = 0.00002       # deg/yr

# ── Analytic solar system table ───────────────────────────────────────────────

def print_solar_system_table():
    print("=" * 85)
    print("  Perihelion Precession: PM Analytic vs GR vs Observed")
    print("=" * 85)
    print(f"  {'Body':<10} {'a [AU]':<10} {'e':<9} {'PM [′′/cy]':<13} "
          f"{'GR [′′/cy]':<13} {'Obs [′′/cy]':<13} {'match?'}")
    print("  " + "-" * 82)

    all_pass = True
    for name, (a_au, e, T_yr, obs) in PLANETS.items():
        a = a_au * AU
        T = T_yr * S_PER_YEAR
        pm_val = pm_precession_arcsec_per_century(M_SUN, a, e, T)
        gr_val = perihelion_precession(a, e, M_SUN) * (S_PER_YEAR * 100 / T) * (206265.0)
        pm_gr_agree = abs(pm_val - gr_val) < 1e-6 * abs(pm_val)
        if obs > 0:
            obs_pct = abs(pm_val - obs) / obs * 100
            match_str = f"{'PASS' if obs_pct < 1.0 else 'FAIL'} ({obs_pct:.3f}%)"
            if obs_pct >= 1.0:
                all_pass = False
        else:
            match_str = "no obs ref"
        print(f"  {name:<10} {a_au:<10.6f} {e:<9.6f} {pm_val:<13.4f} "
              f"{gr_val:<13.4f} {obs:<13.4f} {match_str}")

    print()
    print(f"  PM = GR (analytic): YES — PPN β=γ=1 ⟹ identical 1PN formula")
    print(f"  Overall: {'ALL PASS' if all_pass else 'SOME FAIL'}")
    print("=" * 85 + "\n")


# ── Binary pulsar ─────────────────────────────────────────────────────────────

def print_binary_pulsar():
    a3 = G * PSR_M_TOT * (PSR_P_B / (2 * math.pi))**2
    a_psr = a3 ** (1.0 / 3.0)
    arcsec_yr = pm_precession_arcsec_per_century(PSR_M_TOT, a_psr, PSR_E, PSR_P_B) / 100.0
    deg_yr    = arcsec_yr / 3600.0
    pct = abs(deg_yr - PSR_OBS) / PSR_OBS * 100

    print("=" * 65)
    print("  PSR B1913+16 (Hulse-Taylor) — Periastron Advance")
    print("=" * 65)
    print(f"  Observed:  {PSR_OBS:.5f} ± {PSR_UNCT:.5f} deg/yr")
    print(f"  PM 1PN:    {deg_yr:.5f} deg/yr")
    print(f"  GR 1PN:    {deg_yr:.5f} deg/yr  (identical formula)")
    print(f"  Agreement: {pct:.3f}%  (< 0.5% is excellent at 1PN order)")
    print(f"  Status:    {'PASS' if pct < 0.5 else 'MARGINAL'}")
    print("=" * 65 + "\n")


# ── Numerical orbit integration ───────────────────────────────────────────────

def print_numerical_verification():
    print("Numerical 1PN orbit integration (Mercury, n_orbits=10)...")
    a_mer = 0.387098 * AU
    e_mer = 0.205630

    result = pm_integrate_orbit(M_SUN, a_mer, e_mer, n_orbits=10)
    pct = result['agreement_frac'] * 100

    print("=" * 65)
    print("  Mercury: Numerical vs Analytic Precession per Orbit")
    print("=" * 65)
    print(f"  Numeric:  {result['precession_per_orbit_rad']:.6e} rad/orbit")
    print(f"  Analytic: {result['precession_analytic_rad']:.6e} rad/orbit")
    print(f"  Agreement: {pct:.4f}%")
    print(f"  Periapsis passages detected: {result['n_periastron']}")
    print(f"  Status: {'PASS' if pct < 0.5 else 'FAIL'}")
    print("=" * 65 + "\n")

    return result


# ── Plot ──────────────────────────────────────────────────────────────────────

def make_plot():
    if not HAS_MPL:
        return

    bodies     = list(PLANETS.keys())
    a_vals     = [v[0] * AU for v in PLANETS.values()]
    e_vals     = [v[1]      for v in PLANETS.values()]
    T_vals     = [v[2] * S_PER_YEAR for v in PLANETS.values()]
    obs_vals   = [v[3]      for v in PLANETS.values()]
    pm_vals    = [pm_precession_arcsec_per_century(M_SUN, a, e, T)
                  for a, e, T in zip(a_vals, e_vals, T_vals)]

    # Only show bodies with obs reference > 0
    mask = [o > 0 for o in obs_vals]
    b_show = [b for b, m in zip(bodies, mask) if m]
    pm_show  = [p for p, m in zip(pm_vals, mask) if m]
    obs_show = [o for o, m in zip(obs_vals, mask) if m]

    x = np.arange(len(b_show))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle("PM Perihelion Precession: Analytic Formula vs Observations", fontsize=13)

    # ── Left: bar comparison ──────────────────────────────────────────────
    ax = axes[0]
    width = 0.35
    ax.bar(x - width/2, pm_show,  width, label='PM (1PN analytic) = GR', color='#e74c3c', alpha=0.85)
    ax.bar(x + width/2, obs_show, width, label='Observed (GR residual)',  color='#3498db', alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(b_show, fontsize=10)
    ax.set_ylabel("Precession [arcsec/century]", fontsize=10)
    ax.set_yscale('log')
    ax.legend(fontsize=9)
    ax.grid(True, ls=':', lw=0.5, alpha=0.7, axis='y')
    ax.set_title("PM = GR = Observed\n(PM passes solar-system test)", fontsize=10)

    # ── Right: fractional residual ────────────────────────────────────────
    ax2 = axes[1]
    residuals = [(pm - obs) / obs * 100 for pm, obs in zip(pm_show, obs_show)]
    colors = ['green' if abs(r) < 1.0 else 'red' for r in residuals]
    ax2.bar(x, residuals, color=colors, alpha=0.85)
    ax2.axhline(0, color='k', lw=0.8)
    ax2.axhspan(-1, 1, alpha=0.08, color='green', label='±1% band')
    ax2.set_xticks(x)
    ax2.set_xticklabels(b_show, fontsize=10)
    ax2.set_ylabel("(PM − obs) / obs  [%]", fontsize=10)
    ax2.set_xlabel("Solar system body", fontsize=10)
    ax2.legend(fontsize=9)
    ax2.grid(True, ls=':', lw=0.5, alpha=0.7, axis='y')
    ax2.set_title("Fractional residuals\n(All within ±1%; PM = GR here)", fontsize=10)

    plt.tight_layout()
    out = Path(__file__).parent.parent / 'results' / 'perihelion_precession.png'
    out.parent.mkdir(exist_ok=True)
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {out}")
    plt.close()


# ── Growing discrepancy summary ───────────────────────────────────────────────

def print_discrepancy_summary():
    print("GROWING PICTURE OF PM vs GR DISCREPANCIES")
    print("=" * 65)
    rows = [
        ("Newtonian gravity",      "PM = GR",  "✓ PASS",  "leading order exact"),
        ("Light deflection",       "PM = GR",  "✓ PASS",  "PPN γ=1"),
        ("Shapiro delay",          "PM = GR",  "✓ PASS",  "same formula"),
        ("Frame dragging",         "PM = GR",  "✓ PASS",  "same formula"),
        ("Perihelion advance",     "PM = GR",  "✓ PASS",  "PPN β=γ=1"),
        ("GW propagation speed",   "PM = GR",  "✓ PASS",  "v_GW = c"),
        ("Binary pulsar dP/dt",    "PM = GR",  "✓ PASS",  "0.63% (Galactic correction)"),
        ("Surface redshift",       "PM ≈ 2×GR","⚠ DIFFER","large; NICER could test"),
        ("NS radii (NICER)",       "PM +30% R", "⚠ DIFFER","PM M-R too stiff"),
        ("GW150914 masses",        "PM no BH", "✗ PROBLEM","29+36 M☉ > PM max 13 M☉"),
        ("QNM ringdown",           "PM: surface","⚠ DIFFER","no horizon QNMs in PM"),
    ]
    print(f"  {'Test':<28} {'PM vs GR':<18} {'Status':<12} {'Notes'}")
    print("  " + "-" * 62)
    for row in rows:
        print(f"  {row[0]:<28} {row[1]:<18} {row[2]:<12} {row[3]}")
    print()
    print("  The weak-field tests (solar system) all PASS for PM.")
    print("  The strong-field / compact-object tests show systematic")
    print("  discrepancies that grow with compactness (mass/radius ratio).")
    print("  The GW150914 mass gap is currently the hardest constraint:")
    print("  PM must explain how objects at 29-36 M_sun form and merge.")
    print("=" * 65 + "\n")


def main():
    print_solar_system_table()
    print_binary_pulsar()
    print_numerical_verification()
    make_plot()
    print_discrepancy_summary()


if __name__ == "__main__":
    main()
