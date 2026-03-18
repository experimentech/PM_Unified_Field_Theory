"""
pm_compactness_sweep.py
=======================
Sweep the compactness k = R_star / R_s from 5 down to 1.01 and track:

  * phi_c = phi(0)           central PM potential
  * n(0) = exp(phi_c)        central refractive index
  * n(R_star) = exp(phi_s)   surface refractive index
  * b_crit = F(R_star)       PM shadow radius (star-surface set)
  * F'(r) sign               photon-sphere check

Physics
-------
G = c = 1,  M = 1,  R_s = 2 GM = 2.
R_star = k * R_s,  rho_0 = 3M / (4 pi R_star^3).

Exact analytic field (Section 11 of PM paper):
  Interior:  phi(r) = phi_c - A r^2,   A = GM / R_star^3
  Exterior:  phi(r) = GM / r

  phi_c = 3 GM / (2 R_star)   [= phi_surf + (1/2)*phi_surf, i.e. 3/2 times surf]

F'(r) = 0 in interior at r_ext = R_star^(3/2) / sqrt(2).
  r_ext < R_star  <=>  R_star < R_s  (= 2 in these units).
  So for any R_star >= R_s: F'(r) > 0 everywhere  -> no photon sphere (analytic).

F'(r) = 0 in exterior at r = GM = 1, but 1 < R_star for all k >= 1.
  Never reached by exterior photons.

Conclusion: PM compact objects have NO photon sphere for any R_star >= R_s.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Global constants
# ---------------------------------------------------------------------------
G   = 1.0
M   = 1.0
R_s = 2.0 * G * M   # = 2

FIGURES_DIR = Path(__file__).parent.parent / "docs" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Per-compactness helpers (pure analytic, no spline)
# ---------------------------------------------------------------------------

def make_star(k):
    """Return dict of analytic quantities for compactness k = R_star/R_s.

    Interior phi(r) = phi_c - A*r^2, where A is the SCALAR FIELD coefficient.
    From the ODE solution:  A = GM / (2 R_star^3).
    Note: dphi/dr = -2*A*r = -GM/R_star^3 * r, which at r=R_star gives
    -GM/R_star^2 = (d/dr)(GM/r)|_{R_star}  — matching the exterior.  ✓
    phi_c = phi_s + A*R_star^2 = GM/R_star + GM/(2*R_star) = 3GM/(2*R_star). ✓
    """
    R_star = k * R_s
    rho_0  = 3.0 * M / (4.0 * np.pi * R_star**3)
    A      = G * M / (2.0 * R_star**3)    # phi coeff: (2pi G rho_0)/3 = GM/(2R*^3)
    phi_s  = G * M / R_star               # phi(R_star) = surface potential
    phi_c  = 3.0 * G * M / (2.0 * R_star) # continuity: phi_s + A*R_star^2
    return dict(k=k, R_star=R_star, rho_0=rho_0, A=A, phi_s=phi_s, phi_c=phi_c)


def phi_fn(r, s):
    """phi(r) for star s (analytic)."""
    r = np.atleast_1d(np.asarray(r, float))
    return np.where(r <= s['R_star'],
                    s['phi_c'] - s['A'] * r**2,
                    G * M / r)


def n_fn(r, s):
    return np.exp(phi_fn(r, s))


def F_fn(r, s):
    r = np.atleast_1d(np.asarray(r, float))
    return n_fn(r, s) * r


def Fprime_fn(r, s):
    """Analytic F'(r) = d/dr[exp(phi)*r]."""
    r  = np.atleast_1d(np.asarray(r, float))
    n  = n_fn(r, s)
    dp = np.where(r <= s['R_star'],
                  -2.0 * s['A'] * r,
                  -G * M / r**2)
    return n * (1.0 + r * dp)   # n*(1 + r * dphi/dr)


def b_crit_fn(s):
    """PM shadow radius = F(R_star)."""
    return F_fn(s['R_star'], s)[0]


def r_min_fn(b, s):
    """Turning point: smallest r >= R_star where F(r) = b."""
    R_star = s['R_star']
    r_lo   = R_star
    r_hi   = max(4.0 * b, 5.0 * R_star)
    return brentq(lambda r: F_fn(r, s)[0] - b, r_lo, r_hi, xtol=1e-13)


def deflection(b, s):
    """
    alpha(b) = 2 int_0^1 b/(u*sqrt(F(r_min/u)^2 - b^2)) du - pi
    Returns nan if b <= b_crit.
    """
    bc = b_crit_fn(s)
    if b <= bc * (1.0 - 1e-4):
        return np.nan
    r_m = r_min_fn(b, s)
    b2  = b * b

    def integrand(u):
        if u <= 0.0:
            return 0.0
        r   = r_m / u
        Fu  = F_fn(r, s)[0]
        arg = Fu**2 - b2
        if arg <= 0.0:
            return 0.0
        return b / (u * np.sqrt(arg))

    val, _ = quad(integrand, 0.0, 1.0,
                  points=[0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999],
                  limit=500, epsabs=1e-11, epsrel=1e-8)
    return 2.0 * val - np.pi


# ---------------------------------------------------------------------------
# Compactness sweep
# ---------------------------------------------------------------------------
k_vals = np.array([5.0, 4.0, 3.0, 2.5, 2.0, 1.75, 1.5, 1.25, 1.1, 1.05, 1.01])
stars  = [make_star(k) for k in k_vals]

# Analytic photon-sphere check:
#   Interior extremum at r_ext = R_star^(3/2)/sqrt(2);  inside star iff R_star < R_s
def r_ext_analytic(s):
    """Location of interior F'=0: r^2 = 1/(2A) = R_star^3/GM => r = R_star^(3/2).
    With A = GM/(2*R_star^3): 1/(2A) = R_star^3/GM.
    Enters the star (r_ext < R_star) iff R_star^(3/2) < R_star, i.e. R_star < 1 < R_s/2.
    So for R_star >= R_s=2: no interior photon sphere. (Analytic result.)
    """
    return s['R_star']**1.5 / np.sqrt(G * M)   # = R_star^(3/2) with GM=1

print(f"{'k=R*/Rs':>8}  {'R*':>6}  {'phi_c':>7}  {'n(0)':>7}  {'n(R*)':>7}"
      f"  {'b_crit':>8}  {'b_c/Rs':>7}  {'b_c/R*':>7}  {'F\' zero?':>10}")
rows = []
for s in stars:
    bc    = b_crit_fn(s)
    r_ext = r_ext_analytic(s)
    ps    = "YES (inside)" if r_ext < s['R_star'] else "no"
    row   = dict(k=s['k'], R_star=s['R_star'], phi_c=s['phi_c'],
                 n0=np.exp(s['phi_c']), nRs=np.exp(s['phi_s']),
                 bc=bc, bc_over_Rs=bc/R_s, bc_over_Rstar=bc/s['R_star'],
                 photon_sphere=ps, r_ext=r_ext)
    rows.append(row)
    print(f"{s['k']:8.2f}  {s['R_star']:6.2f}  {s['phi_c']:7.4f}  "
          f"{row['n0']:7.4f}  {row['nRs']:7.4f}  {bc:8.4f}  "
          f"{bc/R_s:7.4f}  {bc/s['R_star']:7.4f}  {ps:>10}")

# ---------------------------------------------------------------------------
# Deflection curves for a subset of compactnesses
# ---------------------------------------------------------------------------
k_plot = [5.0, 3.0, 2.0, 1.5, 1.1]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

print("\nComputing deflection curves …")
defl_curves = {}
for k, col in zip(k_plot, colors):
    s  = make_star(k)
    bc = b_crit_fn(s)
    b_arr = np.concatenate([
        np.geomspace(200 * R_s, 5.0 * bc, 20),
        np.linspace(4.8 * bc, 1.5 * bc, 15, endpoint=False),
        np.linspace(1.45 * bc, 1.01 * bc, 15),
    ])
    b_arr = np.unique(b_arr)[::-1]
    alphas = [deflection(b, s) for b in b_arr]
    defl_curves[k] = dict(b=b_arr, alpha=np.array(alphas), bc=bc, color=col, s=s)
    print(f"  k={k:.2f}  b_crit={bc:.4f}  computed {len(b_arr)} points")

# ---------------------------------------------------------------------------
# Figure 1: n(0) and n(R_star) vs k
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(k_vals, [r['n0']  for r in rows], 'b-o', ms=5, lw=1.5, label=r'$n(0)=e^{\phi_c}$')
ax.plot(k_vals, [r['nRs'] for r in rows], 'r-s', ms=5, lw=1.5, label=r'$n(R_\star)=e^{\phi(R_\star)}$')
ax.axhline(1.0, color='k', lw=0.5, ls='--')
ax.set_xlabel(r'Compactness $k = R_\star / R_s$')
ax.set_ylabel(r'Refractive index $n$')
ax.set_title('PM refractive index vs compactness')
ax.invert_xaxis()
ax.legend(); fig.tight_layout()
fig.savefig(FIGURES_DIR / 'pm_sweep_n.png', dpi=200); plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 2: b_crit / R_s  and  b_crit / R_star  vs k
# ---------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(7, 4))
ax2 = ax1.twinx()
ax1.plot(k_vals, [r['bc_over_Rs']    for r in rows], 'b-o', ms=5, lw=1.5, label=r'$b_{\rm crit}/R_s$')
ax2.plot(k_vals, [r['bc_over_Rstar'] for r in rows], 'r-s', ms=5, lw=1.5, label=r'$b_{\rm crit}/R_\star$')
ax1.set_xlabel(r'Compactness $k = R_\star / R_s$')
ax1.set_ylabel(r'$b_{\rm crit}/R_s$',  color='b')
ax2.set_ylabel(r'$b_{\rm crit}/R_\star$', color='r')
ax1.invert_xaxis()
ax1.set_title(r'PM shadow radius $b_{\rm crit}$ vs compactness')
lines1, labs1 = ax1.get_legend_handles_labels()
lines2, labs2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labs1 + labs2)
fig.tight_layout()
fig.savefig(FIGURES_DIR / 'pm_sweep_bcrit.png', dpi=200); plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 3: F(r) profiles for several compactnesses
# ---------------------------------------------------------------------------
r_plot = np.linspace(1e-3, 25.0, 3000)
fig, ax = plt.subplots(figsize=(8, 5))
for k, col in zip(k_plot, colors):
    s   = make_star(k)
    Fv  = F_fn(r_plot, s)
    bc  = b_crit_fn(s)
    ax.plot(r_plot / R_s, Fv / R_s, color=col, lw=1.5, label=fr'$k={k}$')
    ax.axvline(s['R_star'] / R_s, color=col, ls=':', lw=0.8)
ax.set_xlabel(r'$r / R_s$'); ax.set_ylabel(r'$F(r) / R_s$')
ax.set_title(r'$F(r)=n(r)\,r$ profiles (dashed verticals = $R_\star$)')
ax.set_xlim(0, 20); ax.set_ylim(0)
ax.legend(fontsize=9); fig.tight_layout()
fig.savefig(FIGURES_DIR / 'pm_sweep_F.png', dpi=200); plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 4: F'(r) profiles (confirm no sign change)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
for k, col in zip(k_plot, colors):
    s    = make_star(k)
    Fpv  = Fprime_fn(r_plot, s)
    ax.plot(r_plot / R_s, Fpv, color=col, lw=1.5, label=fr'$k={k}$')
    ax.axvline(s['R_star'] / R_s, color=col, ls=':', lw=0.8)
ax.axhline(0.0, color='k', lw=0.8, ls='--', label=r"$F'=0$  (photon-sphere condition)")
ax.set_xlabel(r'$r / R_s$')
ax.set_ylabel(r"$F'(r) = \mathrm{d}[n(r)\,r]/\mathrm{d}r$")
ax.set_title(r"$F'(r)$ profiles — no zero crossing for any $k \geq 1$")
ax.set_xlim(0, 20)
ax.legend(fontsize=9); fig.tight_layout()
fig.savefig(FIGURES_DIR / 'pm_sweep_Fprime.png', dpi=200); plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 5: alpha(b) curves
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5))
for k, col in zip(k_plot, colors):
    d     = defl_curves[k]
    valid = ~np.isnan(d['alpha'])
    ax.semilogy(d['b'][valid] / R_s, d['alpha'][valid], color=col, lw=1.5,
                label=fr'$k={k}$  ($b_{{\rm crit}}={d["bc"]/R_s:.2f}\,R_s$)')
ax.set_xlabel(r'$b / R_s$'); ax.set_ylabel(r'$\alpha$ (rad, log scale)')
ax.set_title(r'PM deflection angle $\alpha(b)$ for varying compactness')
ax.legend(fontsize=9); fig.tight_layout()
fig.savefig(FIGURES_DIR / 'pm_sweep_deflection.png', dpi=200); plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 6: b_crit / R_s vs k  (clean single-axis for paper)
# ---------------------------------------------------------------------------
k_dense = np.linspace(1.001, 5.5, 500)
bc_dense = np.array([b_crit_fn(make_star(k)) for k in k_dense])
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(k_dense, bc_dense / R_s, 'b-', lw=2)
ax.scatter(k_vals, [r['bc_over_Rs'] for r in rows], color='b', zorder=5)
ax.axhline(1.5 * np.sqrt(3), color='r', ls='--', lw=1,
           label=r'GR photon sphere $b_{\rm crit}^{\rm GR}=\frac{3\sqrt{3}}{2}R_s$')
ax.axvline(1.0, color='gray', ls=':', lw=1, label=r'$R_\star = R_s$ (Schwarzschild limit)')
ax.set_xlabel(r'Compactness $k = R_\star / R_s$')
ax.set_ylabel(r'$b_{\rm crit} / R_s$')
ax.set_title(r'PM shadow radius $b_{\rm crit}$ vs compactness')
ax.invert_xaxis(); ax.legend(); fig.tight_layout()
fig.savefig(FIGURES_DIR / 'pm_sweep_bcrit2.png', dpi=200); plt.close(fig)

print(f"\nAll figures saved to {FIGURES_DIR}")

# ---------------------------------------------------------------------------
# Summary table (paper-ready)
# ---------------------------------------------------------------------------
print("\n" + "="*90)
print("COMPACTNESS SWEEP SUMMARY")
print("="*90)
print(f"{'k':>5}  {'R*/Rs':>6}  {'phi_c':>7}  {'n(0)':>7}  {'n(R*)':>7}  "
      f"{'b_crit/Rs':>10}  {'b_crit/R*':>10}  {'Photon sphere':>14}")
print("-"*90)
for r in rows:
    print(f"{r['k']:5.2f}  {r['k']:6.2f}  {r['phi_c']:7.4f}  "
          f"{r['n0']:7.4f}  {r['nRs']:7.4f}  "
          f"{r['bc_over_Rs']:10.4f}  {r['bc_over_Rstar']:10.4f}  "
          f"{r['photon_sphere']:>14}")
print("="*90)
print(f"\nAnalytic result: interior F'(r)=0 at r_ext = R*^(3/2)/sqrt(GM) = R*^(3/2).")
print(f"r_ext < R*  iff  R* < 1.  So photon sphere requires R* < 1 (<<  R_s = {R_s}).")
print(f"For all k >= 1 (R* >= R_s): PM compact objects have NO photon sphere.")
