#!/usr/bin/env python3
"""
Strong-field light deflection: PM (eikonal/Bouguer) vs GR (Schwarzschild).

PM index: n(r) = 1 + R_s/r,  R_s = 2GM/c^2 (Schwarzschild radius)

Key analytic results derived here:
  - PM turning point:  n(r_min)*r_min = b  =>  r_min = b - R_s
  - PM critical b:     b_crit_PM = R_s   (light captured when b < R_s; no photon sphere)
  - GR critical b:     b_crit_GR = 3*sqrt(3)*GM/c^2 = (3*sqrt(3)/2)*R_s ≈ 2.598 R_s
  - GR photon sphere:  r_ps = 3GM/c^2 = 1.5 R_s

Both deflection integrals use the substitution p = b/r (avoids integrand clustering),
then scipy.integrate.quad for accuracy.

PM deflection (Bouguer eikonal):
  alpha = 2 * integral_0^{p_max} dp / sqrt(1 - p^2*(1-eps^2) + 2*eps*p) - pi
  where eps = R_s/b,  p_max = b/r_min = b/(b-R_s)

GR Schwarzschild deflection:
  alpha = 2 * integral_0^{p_max} dp / sqrt(1 - p^2 + p^3/brs) - pi
  where brs = b/R_s,  p_max from 1 - p_max^2 + p_max^3/brs = 0

Usage:
  source .venv/bin/activate
  python scripts/strong_field_deflection.py

Outputs to docs/figures/.
"""
import math
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))

import numpy as np
from scipy import integrate

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("matplotlib not available; skipping figures")

G = 6.67430e-11
c = 299792458.0
M_sun = 1.989e30
R_s = 2 * G * M_sun / (c * c)

GR_B_CRIT = 1.5 * math.sqrt(3)   # in units of R_s


def gr_p_max(brs: float) -> float:
    """Largest root of f(p) = 1 - p^2 + p^3/brs = 0 (turning point in p=b/r)."""
    f = lambda p: 1.0 - p * p + p ** 3 / brs
    lo, hi = 1.0, 1.0
    step = 1e-4
    while f(hi) > 0 and hi < brs:
        hi += step
    if f(hi) > 0:
        return float('nan')
    lo = hi - step
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def gr_deflection(brs: float) -> float:
    """GR Schwarzschild light deflection for b = brs * R_s (radians)."""
    if brs <= GR_B_CRIT:
        return math.inf
    pm = gr_p_max(brs)
    if math.isnan(pm):
        return math.nan

    def ig(p):
        v = 1.0 - p * p + p ** 3 / brs
        return 1.0 / math.sqrt(max(v, 1e-30))

    res, _ = integrate.quad(ig, 0, pm, points=[0.99 * pm, 0.999 * pm], limit=200)
    return 2.0 * res - math.pi


def pm_deflection(brs: float) -> float:
    """PM eikonal (Bouguer) light deflection for b = brs * R_s (radians).
    Uses n(r) = 1 + R_s/r; turning point r_min = b - R_s."""
    if brs <= 1.0:
        return math.inf
    eps = 1.0 / brs
    p_max = brs / (brs - 1.0)

    def ig(p):
        v = 1.0 - p * p * (1.0 - eps * eps) + 2.0 * eps * p
        return 1.0 / math.sqrt(max(v, 1e-30))

    res, _ = integrate.quad(ig, 0, p_max, points=[0.99 * p_max, 0.999 * p_max], limit=200)
    return 2.0 * res - math.pi


def weak_field(brs: float) -> float:
    """Weak-field analytic: 4GM/c^2 b = 2 R_s / b  (in radians)."""
    return 2.0 / brs


def print_table():
    brs_vals = [1e4, 1e3, 100, 50, 20, 10, 5, 4, 3.0, 2.7, 2.6,
                2.55, 2.0, 1.8, 1.5, 1.3, 1.1, 1.05, 1.01]
    print(f"\nR_s = {R_s:.6e} m  (solar mass)")
    print(f"GR photon sphere: r = 1.5 R_s,  b_crit = {GR_B_CRIT:.4f} R_s")
    print(f"PM no photon sphere;  b_crit = 1.0 R_s\n")
    print(f"{'b/R_s':>8}  {'alpha_PM':>14}  {'alpha_GR':>14}  {'alpha_WF':>14}  {'PM/GR':>8}")
    print("-" * 70)
    for brs in brs_vals:
        pm = pm_deflection(brs)
        gr = gr_deflection(brs)
        wf = weak_field(brs)
        ratio = pm / gr if (math.isfinite(gr) and gr > 0) else float('nan')
        pm_s = f"{pm:.5e}" if math.isfinite(pm) else "  CAPTURE"
        gr_s = f"{gr:.5e}" if math.isfinite(gr) else "  DIVERGE"
        rt_s = f"{ratio:.4f}" if math.isfinite(ratio) else "     ---"
        print(f"{brs:>8.2f}  {pm_s:>14}  {gr_s:>14}  {wf:>14.5e}  {rt_s:>8}")


def make_figures():
    if not HAS_MATPLOTLIB:
        return
    out = ROOT / 'docs' / 'figures'
    out.mkdir(parents=True, exist_ok=True)

    brs_gr = np.concatenate([np.linspace(2.65, 4, 60), np.linspace(4, 10, 40), np.linspace(10, 100, 30)])
    brs_pm = np.concatenate([np.linspace(1.02, 1.5, 50), np.linspace(1.5, 4, 60), np.linspace(4, 100, 40)])
    a_gr = [gr_deflection(b) for b in brs_gr]
    a_pm = [pm_deflection(b) for b in brs_pm]
    a_wf = [weak_field(b) for b in brs_gr]

    # Fig 1: deflection angle
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(brs_gr, a_gr, color='#e31a1c', lw=2, label='GR (Schwarzschild)')
    ax.semilogy(brs_pm, a_pm, color='#1f78b4', lw=2, label='PM (eikonal)')
    ax.semilogy(brs_gr, a_wf, 'k--', lw=1, label='Weak-field')
    ax.axvline(GR_B_CRIT, color='#e31a1c', lw=0.8, ls=':', alpha=0.7,
               label=r'GR $b_{\rm crit}=2.598\,R_s$')
    ax.axvline(1.0, color='#1f78b4', lw=0.8, ls=':', alpha=0.7,
               label=r'PM $b_{\rm crit}=R_s$')
    ax.set_xlabel(r'$b / R_s$', fontsize=12)
    ax.set_ylabel(r'Deflection angle $\alpha$ (rad)', fontsize=12)
    ax.set_title('Strong-field light deflection: PM vs GR', fontsize=12)
    ax.legend(fontsize=9)
    ax.set_xlim(1.0, 20); ax.set_ylim(1e-3, 30)
    ax.grid(True, which='both', ls=':', lw=0.4, alpha=0.7)
    fig.tight_layout()
    fig.savefig(out / 'strong_field_deflection_pm_vs_gr.png', dpi=300)
    plt.close(fig)
    print(f"Saved: {out / 'strong_field_deflection_pm_vs_gr.png'}")

    # Fig 2: PM/GR ratio
    brs_c = [2.65, 2.7, 2.8, 3.0, 3.5, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100]
    ratios = [(b, pm_deflection(b) / gr_deflection(b))
              for b in brs_c if math.isfinite(gr_deflection(b))]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([r[0] for r in ratios], [r[1] for r in ratios], 'o-', color='#2ca02c', lw=2)
    ax.axhline(1.0, color='k', ls='--', lw=1, label='PM = GR')
    ax.set_xlabel(r'$b / R_s$', fontsize=12)
    ax.set_ylabel(r'$\alpha_{\rm PM} / \alpha_{\rm GR}$', fontsize=12)
    ax.set_title('PM vs GR deflection ratio', fontsize=12)
    ax.legend(fontsize=10); ax.grid(True, ls=':', lw=0.4); ax.set_xlim(2.5, 50)
    fig.tight_layout()
    fig.savefig(out / 'strong_field_pm_gr_ratio.png', dpi=300)
    plt.close(fig)
    print(f"Saved: {out / 'strong_field_pm_gr_ratio.png'}")

    # Fig 3: r_min vs b
    brs_r = np.linspace(1.05, 20, 200)
    rmin_pm_vals = [b - 1.0 for b in brs_r]
    rmin_gr_vals = [brs / gr_p_max(brs) if not math.isnan(gr_p_max(brs)) else math.nan for brs in brs_r]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(brs_r, rmin_pm_vals, color='#1f78b4', lw=2, label='PM')
    ax.plot(brs_r, rmin_gr_vals, color='#e31a1c', lw=2, label='GR Schwarzschild')
    ax.plot([GR_B_CRIT], [1.5], 'r*', ms=12, label=r'GR photon sphere ($1.5\,R_s$)')
    ax.axhline(0, color='k', lw=0.5, ls=':')
    ax.set_xlabel(r'$b / R_s$', fontsize=12)
    ax.set_ylabel(r'$r_{\rm min} / R_s$', fontsize=12)
    ax.set_title('Closest approach vs impact parameter', fontsize=12)
    ax.legend(fontsize=9); ax.set_xlim(1, 10); ax.set_ylim(-0.2, 8)
    ax.grid(True, ls=':', lw=0.4)
    fig.tight_layout()
    fig.savefig(out / 'strong_field_rmin_pm_vs_gr.png', dpi=300)
    plt.close(fig)
    print(f"Saved: {out / 'strong_field_rmin_pm_vs_gr.png'}")


if __name__ == '__main__':
    print_table()
    make_figures()
