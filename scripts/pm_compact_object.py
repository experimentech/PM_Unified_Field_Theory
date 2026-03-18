"""
pm_compact_object.py
====================
Self-consistent PM compact object: solve the static scalar potential ODE,
compute the refractive index profile n(r) = exp(phi(r)), locate any
photon-sphere-analogue, and compute the Bouguer light-deflection integral.

Physics recap
-------------
PM static scalar equation (spherical symmetry):

    (1/r^2) d/dr [ r^2 dphi/dr ] = -4 pi G_eff rho(r)

with G_eff = G (recovers Newtonian gravity in the weak-field limit).
Refractive index:  n(r) = exp( phi(r) )

Boundary conditions:
    * dphi/dr |_{r=0} = 0   (regularity at centre)
    * phi(r) -> 0           as r -> inf (asymptotic flatness)

For a uniform-density star of mass M and radius R_star:
    rho(r) = rho_0 = 3M / (4 pi R_star^3)   r <= R_star
    rho(r) = 0                               r >  R_star

Analytic solution (PM: phi = +G_eff*M/r NOT Newtonian -GM/r):
    The PM Poisson equation nabla^2 phi = -4pi G_eff rho  has Green function
    phi = +G_eff M/r  (positive), giving n = exp(phi) > 1 (focusing medium).

    Interior (r <= R_star):
        dphi/dr = -(4 pi G_eff rho_0 / 3) r         (from regularity at centre)
        phi(r)  = phi_c - (2 pi G_eff rho_0 / 3) r^2
        phi_c   = +GM/R_star + (2 pi G_eff rho_0 / 3) R_star^2  (continuity at surface)

    Exterior (r > R_star):
        phi(r)  = +G_eff M / r = +GM / r  > 0
        dphi/dr = -GM / r^2                         (gradient of +GM/r)

Strategy: the analytic solution IS the exact ODE solution.  We build a dense
numerical interpolant of it (which also serves as the "self-consistent" solution)
and verify it satisfies the ODE everywhere.  For completeness we also run
solve_ivp to integrate the ODE inward-outward from the analytic IC at the
surface, confirming the analytic result.

We work in units G = c = 1, M = 1.  R_s = 2GM = 2.

Light-deflection integral (Bouguer / eikonal):
    alpha(b) = 2 int_{r_min}^{inf} [ b / (r sqrt( (n(r) r)^2 - b^2 )) ] dr - pi

where r_min satisfies n(r_min) * r_min = b  (turning point).

A photon-sphere analogue exists if d/dr[n(r)*r] = 0 has a solution.
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Physical / numerical parameters
# ---------------------------------------------------------------------------
G   = 1.0           # Newton constant (geometric units)
M   = 1.0           # Total mass
R_s = 2.0 * G * M  # Schwarzschild radius = 2 in these units

R_star = 5.0 * R_s  # Stellar radius = 10 (compactness R_s/R_star = 0.2)
rho_0  = 3.0 * M / (4.0 * np.pi * R_star**3)  # Uniform density

r_max = 200.0 * R_s  # Outer boundary for integration / plots

FIGURES_DIR = Path(__file__).parent.parent / "docs" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Analytic (exact) solution for phi(r)
# ---------------------------------------------------------------------------
# Interior constants
# PM exterior potential: phi = +G_eff*M/r > 0  (see Section 2.2 of PM theory docs)
phi_surface = G * M / R_star                                           # phi(R_star) = +GM/R_star > 0
phi_c       = phi_surface + (2.0 * np.pi * G * rho_0 / 3.0) * R_star**2  # phi(0) > phi_surface > 0


def phi_analytic(r):
    """Exact solution of the PM scalar equation for a uniform-density sphere."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    phi = np.empty_like(r)
    mask_in = r <= R_star
    phi[mask_in]  = phi_c - (2.0 * np.pi * G * rho_0 / 3.0) * r[mask_in]**2
    phi[~mask_in] = G * M / r[~mask_in]    # PM: +GM/r (positive, focusing)
    return phi


def dphi_analytic(r):
    """Exact first derivative (dphi/dr) of the analytic solution."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    dp = np.empty_like(r)
    mask_in = r <= R_star
    dp[mask_in]  = -(4.0 * np.pi * G * rho_0 / 3.0) * r[mask_in]
    dp[~mask_in] = -G * M / r[~mask_in]**2      # d(+GM/r)/dr = -GM/r^2
    return dp


# ---------------------------------------------------------------------------
# 2. Verify analytic solution by integrate.solve_ivp (shooting from surface)
# ---------------------------------------------------------------------------
# Integrate outward from R_star with IC from analytic solution; must match phi_out = -GM/r.
# Also integrate inward from R_star to r=0; must match phi_c at centre.

def ode_rhs(r, y, rho_val):
    """y = [phi, dphi/dr];  d/dr[r^2 dphi/dr]/r^2 = -4piG rho"""
    dphi = y[1]
    d2phi = -4.0 * np.pi * G * rho_val - (2.0 / r) * dphi
    return [dphi, d2phi]

# Outward integration (exterior, rho=0)
r_out = np.linspace(R_star, r_max, 5000)
ic_out = [phi_analytic(R_star)[0], dphi_analytic(R_star)[0]]
sol_out = solve_ivp(lambda r, y: ode_rhs(r, y, 0.0),
                    [R_star, r_max], ic_out, dense_output=True,
                    t_eval=r_out, rtol=1e-12, atol=1e-15, method='DOP853')

# Inward integration (interior, rho=rho_0); integrate from R_star to ~0
r_in = np.linspace(R_star, 1e-4, 6000)
ic_in = [phi_analytic(R_star)[0], dphi_analytic(R_star)[0]]
sol_in = solve_ivp(lambda r, y: ode_rhs(r, y, rho_0),
                   [R_star, 1e-4], ic_in, dense_output=True,
                   t_eval=r_in, rtol=1e-12, atol=1e-15, method='DOP853')

# Compare at boundary and centre
phi_out_check = sol_out.y[0, :]   # should match -GM/r
phi_analytic_out = G * M / r_out    # PM exterior: +GM/r
max_err_out = np.max(np.abs(phi_out_check - phi_analytic_out))

phi_in_check = sol_in.y[0, :]
phi_analytic_in = phi_c - (2.0 * np.pi * G * rho_0 / 3.0) * r_in**2
max_err_in = np.max(np.abs(phi_in_check - phi_analytic_in))

print(f"solve_ivp vs analytic: exterior max err = {max_err_out:.2e}")
print(f"solve_ivp vs analytic: interior max err = {max_err_in:.2e}")
print(f"phi(0)     = {phi_c:.6f} (analytic), {sol_in.y[0,-1]:.6f} (ivp)")
print(f"phi(R_max) = {G*M/r_max:.2e} (analytic +GM/r), {sol_out.y[0,-1]:.2e} (ivp)")

# ---------------------------------------------------------------------------
# 3. Build a dense cubic-spline interpolant of phi(r) over [0, r_max]
# ---------------------------------------------------------------------------
r_dense = np.concatenate([
    np.linspace(1e-6, R_star, 4000, endpoint=False),
    np.linspace(R_star, r_max, 8000),
])
phi_dense = phi_analytic(r_dense)
dphi_dense = dphi_analytic(r_dense)

# CubicSpline with known derivatives at both ends for maximum accuracy
phi_spline = CubicSpline(r_dense, phi_dense)

# ---------------------------------------------------------------------------
# 4. Refractive index helpers (use spline interpolant)
# ---------------------------------------------------------------------------
r_eps = 1e-4  # small lower bound; avoid r=0 exactly


def phi_of_r(r):
    """phi(r) from the analytic formula — valid for all r >= 0 with no domain cutoff.

    The CubicSpline phi_spline is only defined on [0, r_max] and would extrapolate
    badly at r >> r_max.  The analytic formula phi_analytic is exact for all r,
    so the integral routines (which evaluate r = r_min/u -> inf as u -> 0) work correctly.
    """
    return phi_analytic(np.atleast_1d(np.asarray(r, float)))


def n_of_r(r):
    """n(r) = exp(phi(r))."""
    return np.exp(phi_of_r(r))


def f_of_r(r):
    """f(r) = n(r)*r  (photon-orbit function)."""
    r = np.atleast_1d(np.asarray(r, float))
    return n_of_r(r) * r


# Fine output grid for plotting
r_plot  = np.linspace(r_eps, 30.0 * R_s, 10_000)
phi_plot = phi_of_r(r_plot)
n_plot   = np.exp(phi_plot)
f_plot   = n_plot * r_plot

# ---------------------------------------------------------------------------
# 5. Check for photon-sphere analogue: extremum of f(r) = n(r)*r
# ---------------------------------------------------------------------------
# f'(r) = n'(r)*r + n(r);  a sign change in f' -> extremum / photon sphere.

r_check  = np.linspace(r_eps, 30.0 * R_s, 100_000)
f_check  = f_of_r(r_check)
df_check = np.gradient(f_check, r_check)

sign_changes = np.where(np.diff(np.sign(df_check)))[0]

print(f"\nPhoton-sphere check: f(r)=n(r)*r extrema found: {len(sign_changes)}")
r_extrema = []
for idx in sign_changes:
    ra, rb = r_check[idx], r_check[idx + 1]
    try:
        # df/dr using finite differences of spline
        def df_scalar(r_val):
            dr = r_val * 1e-7
            return (f_of_r(r_val + dr)[0] - f_of_r(r_val - dr)[0]) / (2 * dr)
        r_ext = brentq(df_scalar, ra, rb, xtol=1e-12)
        r_extrema.append(r_ext)
        print(f"  extremum at r = {r_ext:.6f}  ({r_ext/R_s:.4f} R_s)")
    except ValueError:
        pass

if not r_extrema:
    print("  No extremum found – PM compact object has no photon sphere in this configuration.")

# ---------------------------------------------------------------------------
# 6. Light deflection integral
# ---------------------------------------------------------------------------
# alpha(b) = 2 * int_{r_min}^{r_max} [ b / (r * sqrt((n r)^2 - b^2)) ] dr - pi
#
# r_min(b): turning point where n(r_min)*r_min = b.
# Integrable singularity at r = r_min; scipy.quad handles it via 'points'.

# Global minimum of f(r) on (0, r_max]: for a monotone-increasing f this is f(r_eps).
f_min = f_of_r(r_eps)[0]

# Also note: for b < f(R_star) the light ray enters the star.
f_surface = f_of_r(R_star)[0]
print(f"\nf(R_star) = n(R_star)*R_star = {f_surface:.6f}")

# Determine the critical b: if f has no minimum > 0 the star surface sets the
# shadow edge (light rays that graze the surface or enter it).
# b_crit is the smallest b for which r_min >= R_star (grazing shot).
# For monotone f: b_crit = f(R_star).
b_crit = f_surface
print(f"Critical impact parameter b_crit = {b_crit:.6f}  ({b_crit/R_s:.4f} R_s)")


def r_min_for_b(b):
    """Turning point r_min where F(r_min) = n(r_min)*r_min = b.

    With phi = +GM/r, F(r) = exp(GM/r)*r has dF/dr = n(r)*(1 - GM/r).
    The minimum of F is at r = GM = 1 (inside R_star = 10).
    For r >= R_star = 10: F is monotone increasing, so for b > b_crit = F(R_star)
    there is a unique turning point in (R_star, inf).
    """
    r_lo = R_star
    r_hi = max(2.0 * b, 3.0 * R_star)   # F(r) ~ r for large r, so r_hi ~ b suffices
    return brentq(lambda r: f_of_r(r)[0] - b, r_lo, r_hi, xtol=1e-13, maxiter=300)


def deflection_angle(b):
    """
    Full PM deflection alpha(b) = 2 * integral - pi.

    Substitution u = r_min / r maps [r_min, inf) -> (0, 1]:

        I = int_r_min^inf  b / (r * sqrt(F(r)^2 - b^2))  dr
          = int_0^1        b / (u * sqrt(F(r_min/u)^2 - b^2))  du

    where F(r) = n(r)*r and u = r_min/r, dr = -(r_min/u^2) du.
    The integrand is b/(u*sqrt(F^2-b^2)); singularity at u=1 is integrable.

    Returns np.nan for b <= b_crit (light captured by the stellar surface).
    """
    if b <= b_crit * 0.9999:
        return np.nan

    r_m   = r_min_for_b(b)
    b_sq  = b**2
    r_m_  = r_m  # closure

    def integrand_u(u):
        """Integrand in u = r_min/r variable (u in (0,1])."""
        if u <= 0.0:
            return 0.0
        r    = r_m_ / u
        fu   = f_of_r(r)[0]   # n(r)*r
        arg  = fu**2 - b_sq
        if arg <= 0.0:
            return 0.0
        return b / (u * np.sqrt(arg))   # correct PM Bouguer integrand: b/(u*sqrt(F^2-b^2))

    # Singularity at u=1 (turning point), integrand decays as u -> 0 (r -> inf)
    # Split at u=0.999 to help adaptive integrator near the singular end
    try:
        val, _ = quad(integrand_u, 0.0, 1.0,
                      points=[0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999],
                      limit=400, epsabs=1e-10, epsrel=1e-7)
    except Exception as exc:
        print(f"  quad failed at b={b/R_s:.4f} R_s: {exc}")
        return np.nan

    # alpha = 2 * int_0^1  b/(u*sqrt(F^2-b^2)) du  -  pi
    return 2.0 * val - np.pi


# ---------------------------------------------------------------------------
# 7. Compute deflection over a range of impact parameters
# ---------------------------------------------------------------------------
# b_crit is the shadow edge: b < b_crit -> light captured by the star.
# We sweep from large b down to just above b_crit.

b_large  = np.geomspace(200.0 * R_s, 5.0 * b_crit, 40)
b_medium = np.linspace(4.8 * b_crit, 1.5 * b_crit, 25, endpoint=False)
b_near   = np.linspace(1.45 * b_crit, 1.005 * b_crit, 20)
b_all    = np.unique(np.concatenate([b_large, b_medium, b_near]))[::-1]

print(f"\nComputing deflection for {len(b_all)} impact parameters …", flush=True)
alphas = []
for b in b_all:
    a = deflection_angle(b)
    alphas.append(a)
    if not np.isnan(a):
        print(f"  b/R_s = {b/R_s:8.3f}   alpha = {a:.6f} rad")

b_arr     = np.array(b_all)
alpha_arr = np.array(alphas)
valid     = ~np.isnan(alpha_arr)

# Weak-field PM formula for comparison: alpha ~ 2 G_eff M / b = 2GM/b  (with G_eff = G)
alpha_wf = 2.0 * G * M / b_arr

# ---------------------------------------------------------------------------
# 8. Plots
# ---------------------------------------------------------------------------
fig_dir = FIGURES_DIR

# --- 8a. phi(r) ---
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(r_plot / R_s, phi_plot, 'b-', lw=1.5, label=r'$\phi(r)$  (analytic / IVP validated)')
ax.axvline(R_star / R_s, color='gray', ls='--', lw=1, label=r'$R_\star$')
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel(r'$r / R_s$')
ax.set_ylabel(r'$\phi(r)$')
ax.set_title('PM scalar potential – self-consistent solution')
ax.set_xlim(0, 30)
ax.legend()
fig.tight_layout()
fig.savefig(fig_dir / 'pm_co_phi.png', dpi=200)
plt.close(fig)

# --- 8b. n(r) ---
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(r_plot / R_s, n_plot, 'r-', lw=1.5, label=r'$n(r) = e^{\phi(r)}$')
ax.axvline(R_star / R_s, color='gray', ls='--', lw=1, label=r'$R_\star$')
ax.axhline(1.0, color='k', lw=0.5)
ax.set_xlabel(r'$r / R_s$')
ax.set_ylabel(r'$n(r)$')
ax.set_title('PM refractive index – self-consistent solution')
ax.set_xlim(0, 30)
ax.legend()
fig.tight_layout()
fig.savefig(fig_dir / 'pm_co_n.png', dpi=200)
plt.close(fig)

# --- 8c. f(r) = n(r)*r ---
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(r_plot / R_s, f_plot / R_s, 'g-', lw=1.5, label=r'$f(r) = n(r)\,r$')
for r_ext in r_extrema:
    f_ext = f_of_r(r_ext)[0]
    ax.axvline(r_ext / R_s, color='orange', ls=':', lw=1.5,
               label=fr'extremum $r={r_ext/R_s:.2f}\,R_s$ (photon sphere)')
    ax.plot(r_ext / R_s, f_ext / R_s, 'o', color='orange', ms=7)
ax.axvline(R_star / R_s, color='gray', ls='--', lw=1, label=r'$R_\star$')
ax.axhline(b_crit / R_s, color='purple', ls=':', lw=1,
           label=fr'$b_{{crit}}={b_crit/R_s:.3f}\,R_s$')
ax.set_xlabel(r'$r / R_s$')
ax.set_ylabel(r'$f(r)\;/\;R_s$')
ax.set_title(r'$f(r)=n(r)\,r$ — photon sphere condition $f^\prime=0$')
ax.set_xlim(0, 30)
ax.legend()
fig.tight_layout()
fig.savefig(fig_dir / 'pm_co_f.png', dpi=200)
plt.close(fig)

# --- 8d. F'(r) = d/dr[n(r)*r]  (photon-sphere condition: F'=0) ---
r_fp   = np.linspace(r_eps, 30.0 * R_s, 50_000)
f_fp   = f_of_r(r_fp)
df_fp  = np.gradient(f_fp, r_fp)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(r_fp / R_s, df_fp, 'c-', lw=1.5,
        label=r"$F'(r) = \mathrm{d}[n(r)\,r]/\mathrm{d}r$")
ax.axhline(0.0, color='k', lw=0.8, ls='--', label=r"$F'=0$ (photon-sphere condition)")
ax.axvline(R_star / R_s, color='gray', ls='--', lw=1, label=r'$R_\star$')
for r_ext in r_extrema:
    ax.axvline(r_ext / R_s, color='orange', ls=':', lw=1.5,
               label=fr'photon sphere $r={r_ext/R_s:.2f}\,R_s$')
ax.set_xlabel(r'$r / R_s$')
ax.set_ylabel(r"$F'(r)$")
ax.set_title(r"$F'(r)=\mathrm{d}[n(r)\,r]/\mathrm{d}r$ — sign determines focusing")
ax.set_xlim(0, 30)
ax.legend()
fig.tight_layout()
fig.savefig(fig_dir / 'pm_co_fprime.png', dpi=200)
plt.close(fig)

# --- 8e. deflection angle alpha(b) linear scale ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(b_arr[valid] / R_s, alpha_arr[valid], 'b-o', ms=3, lw=1.5,
        label=r'PM self-consistent $\alpha(b)$')
ax.plot(b_arr / R_s, alpha_wf, 'k--', lw=1, label=r'Weak field $2GM/b$')
ax.axvline(b_crit / R_s, color='purple', ls=':', lw=1.5,
           label=fr'$b_{{crit}} = {b_crit/R_s:.3f}\,R_s$ (shadow edge)')
ax.set_xlabel(r'$b / R_s$')
ax.set_ylabel(r'$\alpha$ (rad)')
ax.set_title('PM self-consistent deflection angle vs impact parameter')
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.legend()
fig.tight_layout()
fig.savefig(fig_dir / 'pm_co_deflection.png', dpi=200)
plt.close(fig)

# --- 8e. deflection angle log–log ---
fig, ax = plt.subplots(figsize=(8, 5))
b_pos  = b_arr[valid & (alpha_arr > 0)]
a_pos  = alpha_arr[valid & (alpha_arr > 0)]
aw_pos = 4.0 * G * M / b_pos
ax.loglog(b_pos / R_s, a_pos, 'b-o', ms=3, lw=1.5,
          label=r'PM self-consistent $\alpha(b)$')
ax.loglog(b_pos / R_s, aw_pos, 'k--', lw=1, label=r'Weak field $2GM/b$')
ax.axvline(b_crit / R_s, color='purple', ls=':', lw=1.5,
           label=fr'$b_{{crit}} = {b_crit/R_s:.3f}\,R_s$')
ax.set_xlabel(r'$b / R_s$')
ax.set_ylabel(r'$\alpha$ (rad)')
ax.set_title('PM deflection angle (log–log)')
ax.legend()
fig.tight_layout()
fig.savefig(fig_dir / 'pm_co_deflection_log.png', dpi=200)
plt.close(fig)

# --- 8f. r_min(b) ---
r_mins = np.array([r_min_for_b(b) for b in b_arr[valid]])
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(b_arr[valid] / R_s, r_mins / R_s, 'm-', lw=1.5,
        label=r'$r_{\min}(b)$ PM self-consistent')
ax.plot(b_arr[valid] / R_s, b_arr[valid] / R_s, 'k:', lw=1,
        label=r'$r_{\min}=b$ (vacuum / no refraction)')
ax.axvline(R_star / R_s, color='gray', ls='--', lw=1, label=r'$R_\star$')
ax.axhline(R_star / R_s, color='gray', ls='--', lw=1)
ax.set_xlabel(r'$b / R_s$')
ax.set_ylabel(r'$r_{\min} / R_s$')
ax.set_title(r'PM closest approach $r_{\min}(b)$')
ax.legend()
fig.tight_layout()
fig.savefig(fig_dir / 'pm_co_rmin.png', dpi=200)
plt.close(fig)

print("\nAll figures saved to", fig_dir)
print(f"\nSummary: b_crit = {b_crit:.4f} ({b_crit/R_s:.4f} R_s)")
print(f"         phi_c  = {phi_c:.6f}  (centre potential)")
print(f"         n(0)   = {np.exp(phi_c):.6f}")
print(f"         n(R*)  = {n_of_r(R_star)[0]:.6f}")
print()
print(f"{'b/R_s':>10}  {'alpha_PM (rad)':>16}  {'alpha_WF (rad)':>16}  {'PM/WF':>8}")
for b, a, aw in zip(b_arr[valid], alpha_arr[valid], alpha_wf[valid]):
    print(f"{b/R_s:10.3f}  {a:16.6e}  {aw:16.6e}  {a/aw:8.5f}")
