"""PM Compact Star Structure — exact within the PM framework.

Conceptual foundation
---------------------
PM uses **flat space**.  The optical metric n(r) = e^φ describes wave
propagation in a variable-density medium — it is NOT spacetime curvature.
There is no Einstein equation, no event horizon, no Schwarzschild radius,
and no GR compactness bound GM/(c²R) < 1/2.  These are GR concepts that
do not apply here.

PM force law for massive particles
-----------------------------------
    a_PM = +(c²/2) ∇φ

Light uses both the time-dilation and spatial-compression components of
the optical metric (factor 2), giving α = 4GM/(c²b).  Massive particles
couple only to the time-dilation component (factor 1), giving the 1/2.

Why the PM structure equations are EXACT — not a Newtonian approximation
-------------------------------------------------------------------------
Step 1 — PM Poisson equation:
    ∇²φ = −(8πG/c²) ρ(φ),   ρ(φ) = ρ_nuc exp(φ)

Step 2 — Apply Gauss's theorem in flat space (spherical symmetry):
    ∮ ∇φ · dA = ∫ ∇²φ dV  ⟹  4πr² (dφ/dr) = −(8πG/c²) m(r)

    dφ/dr = −(2G/c²) m(r)/r²  =  −μ_G m(r)/r²       [exact]

Step 3 — Substitute into the force law:
    a_PM = +(c²/2) × (−2G/c²) × m(r)/r²  =  −G m(r)/r²   [exact]

The μ_G = 2G/c² and c²/2 factors cancel algebraically.  The resulting
acceleration is exactly Newtonian at all compactnesses.  This is not a
weak-field limit — it is a consequence of flat-space Gauss + PM Poisson.

Absent corrections (GR artifacts that do NOT appear in PM):
  • metric factor 1/(1 − 2Gm/c²r):          absent — flat space
  • pressure-as-active-mass +4πr³P/c²:      absent — source = ρ only
  • pressure-as-inertial-mass +(P/c²):       absent

Step 4 — Hydrostatic equation (exact):
    dP/dr = ρ(P) × (c²/2) × dφ/dr = −G m(r) ρ(r) / r²

PM self-consistent EOS
----------------------
Integrating dP = ρ(c²/2)dφ with ρ = ρ_nuc exp(φ), P → 0 at ρ → ρ_nuc:

    P = c²(ρ − ρ_nuc)/2

Properties:
  • P(ρ_nuc) = 0:   correct surface condition
  • c_s = c/√2:     causal (sub-luminal) sound speed
  • φ ↔ P are equivalent: φ = ln(1 + 2P/(ρ_nuc c²))

PM maximum density
------------------
From the deformation-energy stability criterion φ_crit = 1:

    ρ_crit = e·ρ_nuc ≈ 6.25×10¹⁷ kg/m³

Above ρ_crit the medium becomes radiation-dominated.  No GR horizon forms —
the medium saturates instead.  A PM star with GM/(c²R) ~ 0.5–0.6 is a
compact horizonless object, not an unphysical singularity.

Structure system (exact within PM)
------------------------------------
    dm/dr  = 4π r² ρ(r)
    dP/dr  = −G m(r) ρ(r) / r²
    dφ/dr  = −μ_G m(r) / r²     [PM Poisson via Gauss; consistent with EOS]

Boundary conditions:
    m(0) = 0,   P(0) = P(ρ_c),   φ(0) = ln(ρ_c/ρ_nuc),   ρ_c ≤ ρ_crit
    Surface at r = R where P(R) = 0;  M_star = m(R)

References
----------
  docs/latex/pm-stellar-structure.tex  — full derivation and predictions
  docs/latex/hamiltonian-formulation-v1.tex §7  — deformation energy
  src/pushing_medium/critical_state.py  — compute_critical_state()
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp

from .critical_state import (
    RHO_NUC,
    pm_density_from_phi,
    compute_critical_state,
)

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
G = 6.67430e-11        # N m² kg⁻²
c = 299792458.0        # m s⁻¹
M_SUN = 1.98892e30     # kg

# Derived from the critical-state computation (exact result: φ_crit = 1)
_CRITICAL = compute_critical_state()
RHO_CRIT: float = _CRITICAL.rho_crit   # ≈ 6.25×10¹⁷ kg m⁻³
MU_G: float = 2.0 * G / (c * c)        # Gravitational index coefficient (m² kg⁻¹)


# ---------------------------------------------------------------------------
# Equation of state
# ---------------------------------------------------------------------------

def pm_eos_pressure(rho: float | np.ndarray) -> float | np.ndarray:
    """PM self-consistent EOS: P = c²(ρ − ρ_nuc)/2.

    Derived by integrating the PM hydrostatic equation dP/dr = ρ(c²/2)∂_r φ
    with ρ(φ) = ρ_nuc exp(φ) and boundary condition P → 0 at ρ → ρ_nuc.

    Parameters
    ----------
    rho : float or ndarray
        Mass density [kg m⁻³].  Values ≤ ρ_nuc return P = 0 (no tensile
        pressure outside matter configuration).

    Returns
    -------
    float or ndarray
        Pressure [Pa].
    """
    rho = np.asarray(rho, dtype=float) if not isinstance(rho, float) else rho
    P = c * c / 2.0 * (rho - RHO_NUC)
    return np.maximum(P, 0.0) if isinstance(P, np.ndarray) else max(P, 0.0)


def pm_eos_density(P: float | np.ndarray) -> float | np.ndarray:
    """Inverse PM EOS: ρ = ρ_nuc + 2P/c².

    Parameters
    ----------
    P : float or ndarray
        Pressure [Pa].  Values ≤ 0 return ρ_nuc.

    Returns
    -------
    float or ndarray
        Mass density [kg m⁻³].
    """
    P = np.asarray(P, dtype=float) if not isinstance(P, float) else P
    rho = RHO_NUC + 2.0 * P / (c * c)
    return np.maximum(rho, RHO_NUC) if isinstance(rho, np.ndarray) else max(rho, RHO_NUC)


def pm_eos_sound_speed() -> float:
    """PM EOS sound speed c_s = c/√2  [m s⁻¹].

    Since P = c²(ρ − ρ_nuc)/2, we have dP/dρ = c²/2, so c_s = c/√2.
    This is causal (c_s < c) and is a fixed property of the EOS, independent
    of density.
    """
    return c / math.sqrt(2.0)


# ---------------------------------------------------------------------------
# Compact-star structure solver
# ---------------------------------------------------------------------------

@dataclass
class StarSolution:
    """Result of a PM compact-star structure integration.

    Attributes
    ----------
    M_star : float
        Total gravitational mass [kg].
    R_star : float
        Stellar radius [m] (where P → 0).
    rho_central : float
        Central density [kg m⁻³].
    phi_central : float
        Central PM compression field φ_c.
    converged : bool
        True if the integrator found the stellar surface (P = 0) within
        the integration domain.
    r : ndarray
        Radial coordinate array [m].
    m : ndarray
        Enclosed mass  m(r)  [kg].
    rho : ndarray
        Mass density  ρ(r)  [kg m⁻³].
    P : ndarray
        Pressure  P(r)  [Pa].
    phi : ndarray
        PM scalar field φ(r) integrated directly from the PM Poisson equation,
        dφ/dr = −μ_G m(r)/r² (exact from Gauss's theorem in flat space).
        Consistent with φ = ln(ρ/ρ_nuc) at every point by construction.
        [dimensionless]
    """
    M_star:       float
    R_star:       float
    rho_central:  float
    phi_central:  float
    converged:    bool
    r:     np.ndarray
    m:     np.ndarray
    rho:   np.ndarray
    P:     np.ndarray
    phi:   np.ndarray


def solve_pm_star(
    rho_central: float,
    r_max: float = 5.0e4,     # 50 km outer boundary [m]
    n_eval: int = 5000,
    rtol: float = 1e-9,
    atol: float = 1e-6,
) -> StarSolution:
    """Integrate the exact PM compact-star structure equations from centre to surface.

    The PM structure system (exact within PM — see module docstring):

        dm/dr  = 4π r² ρ(P)
        dP/dr  = −G m(r) ρ(P) / r²
        dφ/dr  = −μ_G m(r) / r²

    with PM EOS  P = c²(ρ − ρ_nuc)/2  and central density ρ_c ≤ ρ_crit.

    Parameters
    ----------
    rho_central : float
        Central density [kg m⁻³].  Must satisfy 0 < ρ_c ≤ ρ_crit.
    r_max : float
        Maximum integration radius [m].  Integration stops earlier at the
        stellar surface (P = 0).  Default 50 km.
    n_eval : int
        Number of evaluation points for the output arrays.
    rtol, atol : float
        Tolerances for the ODE integrator.

    Returns
    -------
    StarSolution
        Filled data class with the full radial profile and star parameters.

    Raises
    ------
    ValueError
        If rho_central > ρ_crit (PM stability limit exceeded).
    """
    if rho_central > RHO_CRIT * (1.0 + 1e-6):
        raise ValueError(
            f"rho_central = {rho_central:.3e} kg/m³ exceeds PM critical density "
            f"ρ_crit = {RHO_CRIT:.3e} kg/m³.  Stable PM matter configurations "
            f"require ρ_c ≤ ρ_crit."
        )

    P_central = pm_eos_pressure(rho_central)
    phi_central = math.log(rho_central / RHO_NUC)  # φ_c = ln(ρ_c/ρ_nuc)

    # ODE right-hand side: y = [m, P, φ]
    # dφ/dr = −μ_G m/r²  comes from PM Poisson + Gauss (exact in flat space)
    def rhs(r, y):
        m, P, phi = y
        rho = pm_eos_density(P)
        if r < 1.0:
            # Near the centre m → 0; dφ/dr and dP/dr → 0.
            dm_dr = 4.0 * math.pi * r * r * rho
            return [dm_dr, 0.0, 0.0]
        dm_dr   = 4.0 * math.pi * r * r * rho
        dP_dr   = -G * m * rho / (r * r)
        dphi_dr = -MU_G * m / (r * r)   # exact: PM Poisson via Gauss
        return [dm_dr, dP_dr, dphi_dr]

    # Event: surface where P = 0 (pressure hits zero → stellar surface)
    def surface_event(r, y):
        return y[1]  # P = 0

    surface_event.terminal = True
    surface_event.direction = -1    # P decreasing toward zero

    # Start integration at a small but finite radius (avoid r=0 singularity).
    # Near-centre series: m ≈ (4π/3)r³ρ_c, P ≈ P_c, φ ≈ φ_c
    r_start   = 1.0   # 1 m from centre (negligible vs stellar radii ~10 km)
    m_start   = (4.0 / 3.0) * math.pi * r_start**3 * rho_central
    P_start   = P_central     # essentially unchanged at r = 1 m
    phi_start = phi_central   # φ barely changes over 1 m

    r_eval = np.linspace(r_start, r_max, n_eval)

    sol = solve_ivp(
        rhs,
        [r_start, r_max],
        [m_start, P_start, phi_start],

        method='DOP853',
        t_eval=r_eval,
        events=surface_event,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    # Trim output to the surface
    r_arr   = sol.t
    m_arr   = sol.y[0]
    P_arr   = np.maximum(sol.y[1], 0.0)   # clamp P ≥ 0
    phi_arr = sol.y[2]                     # from PM Poisson integration
    rho_arr = pm_eos_density(P_arr)

    # Locate surface: last point where P > 0
    surface_mask = P_arr > 0
    if surface_mask.any():
        i_surf = np.where(surface_mask)[0][-1]
        R_star = r_arr[i_surf]
        M_star = m_arr[i_surf]
        converged = True
    else:
        # Surface was not reached (star extends beyond r_max)
        R_star = r_max
        M_star = m_arr[-1]
        converged = False

    return StarSolution(
        M_star=M_star,
        R_star=R_star,
        rho_central=rho_central,
        phi_central=phi_central,
        converged=converged,
        r=r_arr,
        m=m_arr,
        rho=rho_arr,
        P=P_arr,
        phi=phi_arr,
    )


# ---------------------------------------------------------------------------
# M–R curve sweeper
# ---------------------------------------------------------------------------

def compute_mr_curve(
    n_points: int = 60,
    rho_min_factor: float = 1.01,     # just above ρ_nuc (× factor)
    rho_max_factor: float = 1.0,      # fraction of ρ_crit (≤ 1)
    **solve_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the PM mass–radius (M–R) curve by sweeping central density.

    Sweeps central density from ρ_min = rho_min_factor × ρ_nuc to
    ρ_max = rho_max_factor × ρ_crit, integrating the structure equations at
    each central density.

    Parameters
    ----------
    n_points : int
        Number of central-density values (star models) to integrate.
    rho_min_factor : float
        Lower end as a multiple of ρ_nuc.
    rho_max_factor : float
        Upper end as a fraction of ρ_crit (default 1.0 = exactly the critical
        density, the PM maximum).
    **solve_kwargs
        Passed directly to :func:`solve_pm_star`.

    Returns
    -------
    rho_c : ndarray  [kg m⁻³]
        Sequence of central densities used.
    M : ndarray  [M_☉]
        Corresponding total masses in solar masses.
    R : ndarray  [km]
        Corresponding stellar radii in kilometres.
    """
    rho_c_arr = np.linspace(
        rho_min_factor * RHO_NUC,
        rho_max_factor * RHO_CRIT,
        n_points,
    )

    M_arr = np.full(n_points, np.nan)
    R_arr = np.full(n_points, np.nan)

    for i, rho_c in enumerate(rho_c_arr):
        try:
            star = solve_pm_star(rho_c, **solve_kwargs)
            if star.converged:
                M_arr[i] = star.M_star / M_SUN
                R_arr[i] = star.R_star / 1e3     # m → km
        except (ValueError, Exception):
            pass   # keep as NaN

    return rho_c_arr, M_arr, R_arr
