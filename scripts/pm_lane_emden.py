"""PM Lane-Emden: dimensionless stellar structure ODE.

The PM stellar ODE in dimensional form:
  dm/dr = 4π r² ρ
  dP/dr = -(G m ρ)/r²  (Newtonian; PM EOS is sub-relativistic for ρ ≈ ρ_nuc)

with PM EOS:
  P = (c²/2)(ρ - ρ_nuc)   →   ρ = ρ_nuc + 2P/c²
  dP/dρ = c²/2            →   c_s = c/√2

Natural scales at the stability boundary:
  λ = c_s / √(G ρ_crit)  where ρ_crit = e ρ_nuc, c_s = c/√2
  M_* = ρ_nuc λ³

Dimensionless variables:
  x = r / λ
  θ = (ρ - ρ_nuc) / ρ_nuc    so ρ = ρ_nuc (1 + θ)
  m̃ = m / M_*

Substituting into the structure equations yields:
  dm̃/dx = 4π x² (1 + θ)
  dθ/dx = -(1/e) m̃ / x²

Boundary conditions:
  m̃(0) = 0
  θ(0) = θ_c  (central overdensity, 0 ≤ θ_c ≤ e−1 for stability)

The stability boundary is θ_c = e−1  (ρ_c = ρ_crit = e ρ_nuc).

Key result: M_max = Ξ M_* where Ξ = m̃(surface, θ_c = e−1).
The surface is where ρ → ρ_nuc, i.e. θ → 0  (P → 0).

This script:
  1. Numerically solves the dimensionless ODE
  2. Scans θ_c from 0 to e−1 to produce the M-R curve
  3. Reports Ξ = M_max / M_* and checks for a closed-form candidate
  4. Tests the exact identity μ_G ρ_nuc λ² = 1/e
"""

import sys
import math
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, "src")

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
G = 6.67430e-11        # m³ kg⁻¹ s⁻²
c = 299792458.0        # m s⁻¹
RHO_NUC = 2.3e17       # kg m⁻³  (nuclear saturation density)
MSUN = 1.989e30        # kg

# Derived PM constants
RHO_CRIT = math.e * RHO_NUC                     # stability threshold
C_S = c / math.sqrt(2.0)                         # PM sound speed
LAM = C_S / math.sqrt(G * RHO_CRIT)             # natural length scale [m]
M_STAR = RHO_NUC * LAM**3                        # natural mass scale [kg]
MU_G = 2.0 * G / c**2                           # gravitational index coefficient [m/kg]

print("=" * 60)
print("PM Lane-Emden analysis")
print("=" * 60)
print(f"  ρ_nuc   = {RHO_NUC:.3e} kg/m³")
print(f"  ρ_crit  = {RHO_CRIT:.3e} kg/m³  (= e × ρ_nuc)")
print(f"  c_s     = c / √2 = {C_S:.6e} m/s")
print(f"  λ       = {LAM/1e3:.4f} km")
print(f"  M_*     = {M_STAR/MSUN:.4f} M☉")

# ---------------------------------------------------------------------------
# Exact identity verification
# ---------------------------------------------------------------------------
identity = MU_G * RHO_NUC * LAM**2
print(f"\nExact identity: μ_G ρ_nuc λ² = {identity:.8f}  (should be 1/e = {1/math.e:.8f})")
print(f"  Relative error: {abs(identity - 1/math.e) / (1/math.e):.2e}")

# ---------------------------------------------------------------------------
# Dimensionless ODE
# ---------------------------------------------------------------------------
# The coupling factor appearing in dθ/dx:
# dP/dr = -(Gm ρ)/r²
# With P = c²ρ_nuc θ/2, dP/dr = (c²ρ_nuc/2) dθ/dr = (c²ρ_nuc/2) (dθ/dx)/λ
# And -(Gm ρ)/r² = -(G M_* m̃ ρ_nuc(1+θ))/(λ² x²)
# So: dθ/dx = -(2G M_* ρ_nuc) / (c²ρ_nuc λ²) × m̃(1+θ)/x²
#           = -(2G M_* / c² λ²) × m̃(1+θ)/x²
#           = -(μ_G M_* / λ²) × m̃(1+θ)/x²   [since μ_G = 2G/c²]
#           = -(μ_G ρ_nuc λ) × m̃(1+θ)/x²    [M_* = ρ_nuc λ³]
#           = -(μ_G ρ_nuc λ²) × m̃(1+θ)/(x²) × (1/λ)... let me redo.
#
# Full substitution:
# dm̃/dx = (1/(ρ_nuc λ³)) × 4π (λx)² ρ_nuc(1+θ) × λ = 4π x²(1+θ)
# dθ/dx = (2/(c²ρ_nuc)) × [-(G M_* ρ_nuc m̃ ρ_nuc(1+θ))/(λ² x²)] × λ/(1)
#       = -(2G M_* ρ_nuc m̃ (1+θ)) / (c² ρ_nuc λ x²)
#       = -(2G ρ_nuc λ² m̃ (1+θ)) / (c² x²)   [M_* = ρ_nuc λ³ → M_*/λ = ρ_nuc λ²]
#       = -(μ_G ρ_nuc λ²) × m̃(1+θ)/x²
#       = -(1/e) × m̃(1+θ)/x²                  [exact identity!]
#
# The coupling constant is EXACTLY 1/e — no free parameters.
KAPPA = 1.0 / math.e   # = μ_G ρ_nuc λ² (exact, verified above)


def rhs(x, state):
    """ODE right-hand side for [m̃, θ]."""
    m_tilde, theta = state
    if x < 1e-20:
        return [0.0, 0.0]
    dm = 4.0 * math.pi * x**2 * (1.0 + theta)
    dtheta = -KAPPA * m_tilde * (1.0 + theta) / x**2
    return [dm, dtheta]


def surface_event(x, state):
    """Terminate when θ → 0 (surface: ρ → ρ_nuc, P → 0)."""
    return state[1]  # θ = 0 at surface


surface_event.terminal = True
surface_event.direction = -1


def solve_star(theta_c: float, x_max: float = 40.0):
    """Integrate from near-centre to surface.

    Parameters
    ----------
    theta_c : float
        Dimensionless central overdensity  (ρ_c - ρ_nuc)/ρ_nuc.
    x_max : float
        Maximum dimensionless radius for integration.

    Returns
    -------
    x_surface : float
        Dimensionless surface radius.
    m_surface : float
        Dimensionless total mass m̃.
    """
    # Power-series near x = 0 to avoid division by zero:
    #   m̃ ≈ (4π/3)(1 + θ_c) x³
    #   θ ≈ θ_c - (2π κ/3)(1 + θ_c) x²   [from Taylor expansion of θ]
    x0 = 1e-4
    m0 = (4.0 * math.pi / 3.0) * (1.0 + theta_c) * x0**3
    # Second-order correction for θ off-centre:
    theta0 = theta_c - (2.0 * math.pi * KAPPA / 3.0) * (1.0 + theta_c)**2 * x0**2

    sol = solve_ivp(
        rhs,
        [x0, x_max],
        [m0, theta0],
        events=surface_event,
        rtol=1e-9,
        atol=1e-12,
        dense_output=False,
    )

    if sol.t_events[0].size > 0:
        x_s = sol.t_events[0][0]
        m_s = sol.y_events[0][0][0]
    else:
        # Didn't reach surface — return last point
        x_s = sol.t[-1]
        m_s = sol.y[0, -1]

    return x_s, m_s


# ---------------------------------------------------------------------------
# Scan: M-R curve
# ---------------------------------------------------------------------------
print("\nM-R curve scan (dimensionless units):")
print(f"{'θ_c':>10}  {'ρ_c/ρ_nuc':>12}  {'x_s (λ)':>12}  {'m̃ (M_*)':>12}")
print("-" * 55)

theta_vals = np.linspace(0.02, math.e - 1.0, 60)
results = []
for theta_c in theta_vals:
    x_s, m_s = solve_star(theta_c)
    results.append((theta_c, 1.0 + theta_c, x_s, m_s))
    if abs(theta_c - 0.5) < 0.03 or abs(theta_c - 1.0) < 0.03 or abs(theta_c - (math.e - 1)) < 0.03:
        print(f"{theta_c:10.4f}  {1+theta_c:12.4f}  {x_s:12.4f}  {m_s:12.4f}")

results = np.array(results)

# Find maximum mass
idx_max = np.argmax(results[:, 3])
theta_c_max = results[idx_max, 0]
xi_max = results[idx_max, 3]  # M_max / M_*
rho_c_max = results[idx_max, 1] * RHO_NUC

print("\n" + "=" * 60)
print("MAXIMUM MASS RESULT")
print("=" * 60)
print(f"  θ_c at M_max   = {theta_c_max:.6f}  (stability cap θ_max = e-1 = {math.e-1:.6f})")
print(f"  ρ_c at M_max   = {results[idx_max,1]:.4f} × ρ_nuc  = {rho_c_max:.4e} kg/m³")
print(f"  Ξ = M_max/M_*  = {xi_max:.6f}")
print(f"  M_max          = {xi_max * M_STAR / MSUN:.4f} M☉")
print(f"  R at M_max     = {results[idx_max, 2] * LAM / 1e3:.2f} km")

# Specifically solve at the exact stability cap θ_c = e - 1
print("\nExact cap θ_c = e - 1:")
x_cap, m_cap = solve_star(math.e - 1.0)
print(f"  Ξ (cap)        = {m_cap:.8f}")
print(f"  M_max (cap)    = {m_cap * M_STAR / MSUN:.4f} M☉")
print(f"  R at cap       = {x_cap * LAM / 1e3:.3f} km")

# ---------------------------------------------------------------------------
# Closed-form candidates for Ξ = 3.28...
# ---------------------------------------------------------------------------
xi = m_cap
print(f"\nΞ = {xi:.8f}")
print("\nClosed-form candidates:")
candidates = {
    "π/√(e)":         math.pi / math.sqrt(math.e),
    "4π/(3√(2e))":    4*math.pi / (3*math.sqrt(2*math.e)),
    "e + 1/e":        math.e + 1/math.e,
    "√(4π²/e)":       math.sqrt(4*math.pi**2 / math.e),
    "π(e-1)":         math.pi * (math.e - 1),
    "4e/√(3π)":       4*math.e / math.sqrt(3*math.pi),
    "3.28":           3.28,
}
for name, val in candidates.items():
    diff = abs(val - xi)
    print(f"  {name:<20} = {val:.8f}  (diff = {diff:.2e})")

print(f"\n  Best match: the prefactor Ξ ≈ {xi:.4f} has no obvious closed form.")
print(f"  It is determined numerically from the PM stellar ODE.")
