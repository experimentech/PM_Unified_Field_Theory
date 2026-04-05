"""PM Compact Star Structure — exact within the PM framework.

Conceptual foundation
---------------------
PM uses a **flat Minkowski background** — but the medium field n(r) = e^φ
creates a **curved effective optical metric** for all wave propagation and
dynamics.  These are distinct levels:

  • Background geometry:   flat (Euclidean/Minkowski).  No Einstein equations,
                           no Riemann curvature of spacetime, no horizons.
  • Effective optical metric: ds² = −(c²/n²)dt² + dr²  — geometrically curved.
                           Null geodesics bend, time dilation occurs in this metric.

The phrase "flat space" alone is imprecise.  The correct statement is:
"flat background, but curved effective geometry from the medium density field."
There is no GR compactness bound GM/(c²R) < 1/2 and no event horizon because
these are properties of background curvature, which PM does not have.

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

Step 2 — Apply Gauss's theorem on the flat Minkowski background (spherical symmetry):
    ∮ ∇φ · dA = ∫ ∇²φ dV  ⟹  4πr² (dφ/dr) = −(8πG/c²) m(r)

    dφ/dr = −(2G/c²) m(r)/r²  =  −μ_G m(r)/r²       [exact]

Step 3 — Substitute into the force law:
    a_PM = +(c²/2) × (−2G/c²) × m(r)/r²  =  −G m(r)/r²   [exact]

The μ_G = 2G/c² and c²/2 factors cancel algebraically.  The resulting
acceleration is exactly Newtonian at all compactnesses.  This is not a
    weak-field limit — it is a consequence of Gauss's theorem on the flat
    Minkowski background + PM Poisson equation.

Absent corrections (GR artifacts that do NOT appear in PM):
  • metric factor 1/(1 − 2Gm/c²r):          absent — no background curvature
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

# ---------------------------------------------------------------------------
# Tidal Love number and deformability
# ---------------------------------------------------------------------------

def pm_clairaut_k2(y_R: float) -> float:
    """Tidal Love number k₂ from the Clairaut-Radau equation (Newtonian/PM).

    Parameters
    ----------
    y_R : float
        y(R) at the stellar surface, from integrating the Clairaut-Radau ODE:
            dy/dr = (−y² − 5y + 6ν(1+y)) / r
        where ν = ρ/ρ̄ = ρ·4πr³/(3m(r)).  Initial condition y(0) = 2.

    Returns
    -------
    float
        k₂ (Newtonian tidal Love number, > 0 for y_R > 2).

    Notes
    -----
    The Clairaut-Radau ODE is the correct Newtonian formulation for tidal
    deformability in flat-space (PM) gravity.

    For a uniform-density star: ν = 1 throughout → y_R = 3 → k₂ = 1/18 ≈ 0.056.
    For compressible PM stars: y_R ≈ 2.7–3.0 → k₂ ≈ 0.04–0.07.

    PM stars have near-uniform density (ρ ≈ ρ_nuc throughout), so y_R ≈ 2.7–2.9
    and k₂ ≈ 0.04–0.07, giving Λ₁.₄ ≈ 300–500.  This is within but at the
    upper end of the GW170817 90% CI (Λ̃ < 800), in ~1–2σ tension with the
    preferred value (Λ̃ ~ 70–300 for favoured GR EOSs).

    References
    ----------
    Poisson & Will, "Gravity" (2014), eq. 14.34 (Clairaut-Radau equation).
    Abbott et al. 2018, PRL 121, 161101 (GW170817 tidal constraint).
    """
    return (y_R - 2.0) / (2.0 * (2.0 * y_R + 3.0))


def pm_tidal_lambda(k2: float, C: float) -> float:
    """Dimensionless tidal deformability Λ = (2/3) k₂ / C⁵.

    This is the quantity measured by LIGO/Virgo from the GW phase.
    GW170817 constraint: Λ̃ ≤ 800 (combined binary parameter, 90% CL).

    PM stars with large radii (R ~ 25–35 km vs GR ~ 12 km) give Λ ∝ R⁵/M⁵;
    the typical PM prediction is Λ ~ 3000–40000 for a 1.4 M☉ star,
    roughly 4–50× the GW170817 upper bound.
    """
    return (2.0 / 3.0) * k2 / C**5


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
        dφ/dr = −μ_G m(r)/r² (exact from Gauss's theorem on flat Minkowski background).
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
    k2_love:      Optional[float] = None    # Tidal Love number k₂ (Hinderer 2008)
    lambda_tidal: Optional[float] = None    # Dimensionless tidal deformability Λ = (2/3)k₂/C⁵


def solve_pm_star(
    rho_central: float,
    r_max: float = 5.0e4,     # 50 km outer boundary [m]
    n_eval: int = 5000,
    rtol: float = 1e-9,
    atol: float = 1e-6,
    compute_tidal: bool = False,
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
    compute_tidal : bool
        If True, also integrate the Newtonian tidal perturbation ODE alongside
        the stellar structure to compute the tidal Love number k₂ and the
        dimensionless tidal deformability Λ = (2/3)k₂/C⁵.  Stored in
        StarSolution.k2_love and .lambda_tidal.

    Returns
    -------
    StarSolution
        Filled data class with the full radial profile and star parameters.
        If compute_tidal=True, also sets k2_love and lambda_tidal.

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

    if not compute_tidal:
        # ── 3-variable ODE: y = [m, P, φ] ──────────────────────────────────
        def rhs(r, y):
            m, P, phi = y
            rho = pm_eos_density(P)
            if r < 1.0:
                dm_dr = 4.0 * math.pi * r * r * rho
                return [dm_dr, 0.0, 0.0]
            dm_dr   = 4.0 * math.pi * r * r * rho
            dP_dr   = -G * m * rho / (r * r)
            dphi_dr = -MU_G * m / (r * r)
            return [dm_dr, dP_dr, dphi_dr]

        def surface_event(r, y):
            return y[1]

        surface_event.terminal  = True
        surface_event.direction = -1

        r_start   = 1.0
        m_start   = (4.0 / 3.0) * math.pi * r_start**3 * rho_central
        P_start   = P_central
        phi_start = phi_central
        y0 = [m_start, P_start, phi_start]

    else:
        # ── 4-variable ODE: y = [m, P, φ, y2] ──────────────────────────────
        # y2 is the Clairaut-Radau variable for the ℓ=2 tidal perturbation.
        #
        # Clairaut-Radau ODE (Newtonian, correct for PM flat-space gravity):
        #   dy2/dr = (−y2² − 5y2 + 6ν(1+y2)) / r
        # where ν = ρ/ρ̄ = ρ·4πr³/(3m(r)) is the local/mean-density ratio.
        #
        # Initial condition: y2 → 2 at r → 0 (from l=2 regular solution).
        # Near r=0 with uniform density: ν → 1, dy2/dr → (−4−10+12)/r = −2/r
        # → stable fixed-point structure starting from y=2.
        #
        # Physical interpretation:
        # For a near-uniform-density PM star, ν ≈ 1 throughout, so y2 evolves
        # toward 3 (the incompressible fixed point) → k₂ ≈ 1/18 ≈ 0.056.
        # This gives Λ₁.₄ ≈ 300–500, within the GW170817 90% CI (< 800) but
        # ~1–2σ above the preferred value (~70–300 for GR EOSs).
        #
        # Reference: Poisson & Will, "Gravity" (2014), eq. 14.34.

        def rhs(r, y):
            m, P, phi, y2 = y
            rho = pm_eos_density(P)
            dm_dr   = 4.0 * math.pi * r * r * rho
            if r < 1.0 or m < 1.0:
                return [dm_dr, 0.0, 0.0, 0.0]
            dP_dr    = -G * m * rho / (r * r)
            dphi_dr  = -MU_G * m / (r * r)
            # Clairaut-Radau: ν = ρ / ρ̄,  ρ̄ = 3m/(4πr³)
            rho_bar  = 3.0 * m / (4.0 * math.pi * r * r * r)
            nu       = rho / rho_bar if rho_bar > 0.0 else 1.0
            dy2_dr   = (-y2 * y2 - 5.0 * y2 + 6.0 * nu * (1.0 + y2)) / r
            return [dm_dr, dP_dr, dphi_dr, dy2_dr]

        def surface_event(r, y):
            return y[1]

        surface_event.terminal  = True
        surface_event.direction = -1

        r_start   = 1.0
        m_start   = (4.0 / 3.0) * math.pi * r_start**3 * rho_central
        P_start   = P_central
        phi_start = phi_central
        y2_start  = 2.0   # exact to better than 1e-6 at r = 1 m
        y0 = [m_start, P_start, phi_start, y2_start]

    r_eval = np.linspace(r_start, r_max, n_eval)

    sol = solve_ivp(
        rhs,
        [r_start, r_max],
        y0,
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

    # Compute tidal Love number and deformability if requested
    k2_love      = None
    lambda_tidal = None
    if compute_tidal and converged:
        y2_arr = sol.y[3]
        y_R    = float(y2_arr[i_surf])
        C      = G * M_star / (c * c * R_star)   # compactness
        if y_R > 2.0:   # CR gives y_R > 2 for physical stars
            k2_love      = pm_clairaut_k2(y_R)
            lambda_tidal = pm_tidal_lambda(k2_love, C)

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
        k2_love=k2_love,
        lambda_tidal=lambda_tidal,
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


def compute_lambda_M_curve(
    n_points: int = 40,
    rho_min_factor: float = 1.01,
    rho_max_factor: float = 1.0,
    **solve_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the PM tidal deformability–mass curve Λ(M).

    Sweeps central density and returns M, R, and Λ for each converged star.
    Uses the Clairaut-Radau tidal ODE (see pm_clairaut_k2).

    Parameters
    ----------
    n_points : int
        Number of central-density values to integrate.
    rho_min_factor, rho_max_factor : float
        Range of central densities as multiples of ρ_nuc (min) and ρ_crit (max).

    Returns
    -------
    M : ndarray  [M_☉]
    R : ndarray  [km]
    Lambda : ndarray  [dimensionless]
    k2 : ndarray  [dimensionless]
    """
    rho_c_arr = np.linspace(
        rho_min_factor * RHO_NUC,
        rho_max_factor * RHO_CRIT,
        n_points,
    )
    M_arr  = np.full(n_points, np.nan)
    R_arr  = np.full(n_points, np.nan)
    Lam_arr = np.full(n_points, np.nan)
    k2_arr  = np.full(n_points, np.nan)

    for i, rho_c in enumerate(rho_c_arr):
        try:
            star = solve_pm_star(rho_c, compute_tidal=True, **solve_kwargs)
            if star.converged and star.lambda_tidal is not None:
                M_arr[i]   = star.M_star / M_SUN
                R_arr[i]   = star.R_star / 1e3
                Lam_arr[i] = star.lambda_tidal
                k2_arr[i]  = star.k2_love
        except (ValueError, Exception):
            pass

    return M_arr, R_arr, Lam_arr, k2_arr


# ---------------------------------------------------------------------------
# Two-zone solver: matter-phase shell + energy-phase core
# ---------------------------------------------------------------------------
#
# Physical basis
# --------------
# PM's EOS P = c²/2(ρ − ρ_nuc) is only valid in the matter phase: φ < 1, i.e.
# ρ < ρ_crit = e·ρ_nuc.  Stars with ρ_central > ρ_crit contain a core where
# the medium has undergone the second-order phase transition (m²_eff → 0 at
# φ = 1, negative above).  The matter-phase EOS cannot be extrapolated there.
#
# Two-zone structure:
#   r  < r_core :  energy-phase core  — φ ≥ 1 throughout
#   r_core ≤ r  :  matter-phase shell — PM EOS valid
#
# The phase boundary is at φ = 1  ⟺  ρ = ρ_crit.  Inside the boundary, the
# integration switches to the energy-phase EOS.
#
# Energy-phase EOS options
# ------------------------
# The unambiguous PM prediction for the energy-phase pressure comes from the
# deformation potential at the transition:
#
#   U(φ_crit) = ε₀(φ² − φ³/3)|_{φ=1} = ε₀ · 2/3  where ε₀ = ρ_nuc c²/2
#
# This gives  P_core = (2/3) ε₀ = (c²/3) ρ_nuc
#
# Three natural choices are implemented:
#   'constant'  : P_core = U(1) = c² ρ_nuc / 3  (deformation energy at transition)
#   'radiation' : P = ρ c² / 3                   (ultra-relativistic / radiation EOS)
#   'vacuum'    : P_core = 0                     (collapsed core, no pressure support)
#
# The 'constant' EOS is the PM-theory-derived choice.  'radiation' is what
# critical_state.py already used for P_crit.  'vacuum' is the lower bound.
#
# Gravitational structure in the core uses the same ODE (mass accumulation
# continues), because the energy-phase density still sources gravity.  The
# energy density of the core is treated as ρ_core = ρ_crit (constant at the
# transition value), consistent with a stalled phase transition.

def _energy_phase_eos(rho: float, mode: str, eps0: float) -> float:
    """Pressure in the energy-phase core.

    Parameters
    ----------
    rho : float
        Local density [kg m⁻³].  Ignored for 'constant' and 'vacuum' modes.
    mode : str
        'constant'  — P = U(φ_crit) = (2/3)ε₀ = 2ρ_nuc c²/3
        'radiation' — P = ρ c²/3
        'vacuum'    — P = 0
    eps0 : float
        ε₀ = ρ_nuc c²  [J m⁻³]  (same convention as critical_state.py)
    """
    if mode == 'constant':
        return 2.0 * eps0 / 3.0          # U(1) = ε₀(1 − 1/3) = 2ε₀/3
    elif mode == 'radiation':
        return rho * c * c / 3.0
    elif mode == 'vacuum':
        return 0.0
    else:
        raise ValueError(f"Unknown energy-phase EOS mode: {mode!r}")


def solve_pm_star_two_zone(
    rho_central: float,
    core_eos: str = 'constant',
    r_max: float = 5.0e4,
    n_eval: int = 5000,
    rtol: float = 1e-9,
    atol: float = 1e-6,
) -> StarSolution:
    """Integrate PM stellar structure with a two-zone (matter + energy phase) model.

    For central densities ρ_c ≤ ρ_crit, the result is identical to
    :func:`solve_pm_star` (no energy-phase core exists).

    For ρ_c > ρ_crit, the core (r < r_core, where φ reaches 1) uses the
    energy-phase EOS and the shell uses the standard PM matter EOS.

    Parameters
    ----------
    rho_central : float
        Central density [kg m⁻³].  May exceed ρ_crit.
    core_eos : str
        Energy-phase EOS: 'constant' (PM-derived, default), 'radiation', 'vacuum'.
    r_max, n_eval, rtol, atol
        Passed to the ODE integrator.

    Returns
    -------
    StarSolution
        Same structure as :func:`solve_pm_star`, with an additional attribute
        ``r_core`` accessible via ``sol.r_core`` (the core radius in metres,
        NaN if no core exists).
    """
    eps0 = RHO_NUC * c * c          # ε₀ = ρ_nuc c²  (consistent with critical_state.py)
    has_core = rho_central > RHO_CRIT

    # ------------------------------------------------------------------ #
    # Special case: vacuum core + φ > 1.                                  #
    # P = 0 throughout.  At r_core (φ = 1) the matter shell would start   #
    # at P = 0 which IS the surface condition, so the star has no shell.   #
    # A standard 3-variable ODE is numerically stiff here because the full  #
    # ODE tries to integrate into the matter phase with P = 0 and ρ = ρ_nuc, #
    # causing dP/dr = −Gmρ_nuc/r² → step-size → 0 immediately.            #
    # Instead, solve a simpler 2-variable (m, φ) ODE in the core only.    #
    # ------------------------------------------------------------------ #
    if core_eos == 'vacuum' and has_core:
        phi_central = 1.0 + math.log(rho_central / RHO_CRIT)

        def _rhs_vac(r, y2):
            m2, phi2 = y2
            rho = RHO_CRIT
            return [4.0 * math.pi * r * r * rho, -MU_G * m2 / (r * r)]

        def _phase_ev(r, y2):
            return y2[1] - 1.0   # φ − 1 = 0

        _phase_ev.terminal  = True
        _phase_ev.direction = -1

        r_start = 1.0
        sol2 = solve_ivp(
            _rhs_vac,
            [r_start, r_max],
            [(4.0 / 3.0) * math.pi * r_start**3 * rho_central, phi_central],
            events=[_phase_ev],
            rtol=rtol, atol=atol, dense_output=False,
        )

        if sol2.t_events[0].size > 0:
            r_core_val = float(sol2.t_events[0][0])
            M_vac      = float(sol2.y_events[0][0][0])       # m at r_core
            r_arr_vac  = np.append(sol2.t, r_core_val)
            m_arr_vac  = np.append(sol2.y[0], M_vac)
            phi_arr_vac = np.append(sol2.y[1], 1.0)
            P_arr_vac  = np.zeros_like(r_arr_vac)
            rho_arr_vac = np.full_like(r_arr_vac, RHO_CRIT)
            result = StarSolution(
                M_star=M_vac,
                R_star=r_core_val,
                rho_central=rho_central,
                phi_central=phi_central,
                converged=True,
                r=r_arr_vac,
                m=m_arr_vac,
                rho=rho_arr_vac,
                P=P_arr_vac,
                phi=phi_arr_vac,
            )
        else:
            # Core fills r_max — never reached phase boundary
            result = StarSolution(
                M_star=float(sol2.y[0, -1]),
                R_star=r_max,
                rho_central=rho_central,
                phi_central=phi_central,
                converged=False,
                r=sol2.t, m=sol2.y[0],
                rho=np.full(sol2.t.size, RHO_CRIT),
                P=np.zeros(sol2.t.size),
                phi=sol2.y[1],
            )
        result.r_core = r_core_val if sol2.t_events[0].size > 0 else math.nan  # type: ignore[attr-defined]
        return result

    if has_core:
        phi_central = 1.0 + math.log(rho_central / RHO_CRIT)   # φ_c > 1
        # For radiation: P varies freely; start at P(ρ_central) = ρc²/3.
        # For constant: P is fixed at U(1) = 2ε₀/3 throughout. dP/dr=0, P = P_central throughout.
        #   At r_core the matter-phase shell inherits this pressure and integrates down to P=0.
        P_central = _energy_phase_eos(rho_central, core_eos, eps0)
    else:
        phi_central = math.log(rho_central / RHO_NUC)
        P_central = pm_eos_pressure(rho_central)

    # The ODE state is y = [m, P, φ].
    # The phase boundary is tracked via φ: when φ falls below 1 (outward),
    # we switch from energy-phase to matter-phase EOS.
    # Because the ODE is integrated outward from the centre, we enter the
    # matter phase once φ drops below 1.

    r_core_found = [math.nan]   # mutable container for the core radius

    def rhs(r, y):
        m, P, phi = y
        if r < 1.0:
            # Near-centre: avoid singularity
            rho = rho_central
            dm_dr = 4.0 * math.pi * r * r * rho
            return [dm_dr, 0.0, 0.0]

        # Determine which phase we are in based on φ
        in_core = phi >= 1.0
        if in_core:
            # Energy-phase core.
            # Density: for radiation, ρ = 3P/c² (normal fluid EOS);
            #          for constant/vacuum, ρ is pinned at ρ_crit (phase boundary).
            if core_eos == 'radiation':
                rho = max(3.0 * P / (c * c), RHO_CRIT)   # radiation fluid in core
                dP_dr = -G * m * rho / (r * r)
            else:
                # 'constant' or 'vacuum': pressure is fixed by the phase state,
                # not by hydrostatic weight.  dP/dr = 0 in the core.
                # Mass still accumulates at ρ_crit.
                rho = RHO_CRIT
                dP_dr = 0.0
        else:
            rho = pm_eos_density(max(P, 0.0))
            dP_dr = -G * m * rho / (r * r)

        dm_dr   = 4.0 * math.pi * r * r * rho
        dphi_dr = -MU_G * m / (r * r)
        return [dm_dr, dP_dr, dphi_dr]

    # Event 1: stellar surface (P = 0 in matter phase, φ < 1)
    # N.B. we return (P + 1) in the core so the function is continuous
    # across the phase boundary (avoids integrator stalling on a
    # sudden 1→0 jump when P_core = 0 in vacuum mode).
    def surface_event(r, y):
        m, P, phi = y
        if phi >= 1.0:
            return P + 1.0      # always > 0 inside core; continuous at φ=1
        return P                # matter shell: fires when P → 0

    surface_event.terminal = True
    surface_event.direction = -1

    # Event 2: phase boundary (φ = 1, crossing outward into matter phase)
    def phase_boundary_event(r, y):
        return y[2] - 1.0       # φ − 1 = 0

    # For vacuum mode the matter shell has P = 0 immediately at r_core,
    # so the phase boundary IS the stellar surface.  Stop there.
    _vacuum = (core_eos == 'vacuum' and has_core)
    phase_boundary_event.terminal = _vacuum
    phase_boundary_event.direction = -1       # φ decreasing through 1

    r_start   = 1.0
    m_start   = (4.0 / 3.0) * math.pi * r_start**3 * rho_central
    P_start   = P_central
    phi_start = phi_central

    r_eval = np.linspace(r_start, r_max, n_eval)

    sol = solve_ivp(
        rhs,
        [r_start, r_max],
        [m_start, P_start, phi_start],
        method='DOP853',
        t_eval=r_eval,
        events=[surface_event, phase_boundary_event],
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    # Extract core radius from phase-boundary event (event index 1)
    if sol.t_events[1].size > 0:
        r_core_found[0] = float(sol.t_events[1][0])

    r_arr   = sol.t
    m_arr   = sol.y[0]
    P_arr   = np.maximum(sol.y[1], 0.0)
    phi_arr = sol.y[2]
    # Density: use matter-phase EOS in shell, ρ_crit in core
    rho_arr = np.where(phi_arr >= 1.0, RHO_CRIT, pm_eos_density(P_arr))

    surface_mask = (P_arr > 0) & (phi_arr < 1.0)
    if surface_mask.any():
        i_surf    = np.where(surface_mask)[0][-1]
        R_star    = r_arr[i_surf]
        M_star    = m_arr[i_surf]
        converged = True
    elif _vacuum and not math.isnan(r_core_found[0]):
        # Vacuum star: no matter shell, surface = core boundary
        R_star    = r_core_found[0]
        M_star    = m_arr[-1]
        converged = True
    else:
        R_star    = r_max
        M_star    = m_arr[-1]
        converged = False

    result = StarSolution(
        M_star=rho_central,
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
    # Patch M_star (dataclass is not frozen, can assign)
    result.M_star = M_star
    # Attach core radius as extra attribute
    result.r_core = r_core_found[0]    # type: ignore[attr-defined]
    return result


def compute_mr_curve_two_zone(
    n_points: int = 60,
    rho_min_factor: float = 1.01,
    rho_max_factor: float = 3.0,
    core_eos: str = 'constant',
    **solve_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """M–R curve sweeping above and below ρ_crit with the two-zone model.

    Parameters
    ----------
    n_points : int
        Number of central density models to integrate.
    rho_min_factor : float
        Lower bound as multiple of ρ_nuc.
    rho_max_factor : float
        Upper bound as multiple of ρ_crit (default 3.0 — well into energy-phase
        core territory).
    core_eos : str
        Energy-phase EOS for cores: 'constant', 'radiation', or 'vacuum'.

    Returns
    -------
    rho_c, M, R, r_core : ndarrays
        Central densities [kg/m³], masses [M☉], radii [km],
        core radii [km] (NaN for stars without cores).
    """
    rho_c_arr = np.geomspace(
        rho_min_factor * RHO_NUC,
        rho_max_factor * RHO_CRIT,
        n_points,
    )

    M_arr      = np.full(n_points, np.nan)
    R_arr      = np.full(n_points, np.nan)
    r_core_arr = np.full(n_points, np.nan)

    for i, rho_c in enumerate(rho_c_arr):
        try:
            star = solve_pm_star_two_zone(rho_c, core_eos=core_eos, **solve_kwargs)
            if star.converged:
                M_arr[i]      = star.M_star / M_SUN
                R_arr[i]      = star.R_star / 1e3
                r_core_val    = getattr(star, 'r_core', math.nan)
                r_core_arr[i] = r_core_val / 1e3 if math.isfinite(r_core_val) else math.nan
        except Exception:
            pass

    return rho_c_arr, M_arr, R_arr, r_core_arr


# ---------------------------------------------------------------------------
# Option-A: corrected Poisson equation with U'(φ) self-coupling
# ---------------------------------------------------------------------------

def solve_pm_star_option_a(
    rho_central: float,
    alpha: float = 1.0,
    r_max: float = 5.0e4,
    n_eval: int = 5000,
    rtol: float = 1e-9,
    atol: float = 1e-6,
) -> StarSolution:
    """PM compact-star structure with Option-A corrected Poisson self-coupling.

    Physical motivation
    -------------------
    The standard PM Poisson equation (Gap 1 of the formula-sheet diagnostic) is:

        ∇²φ = −(8πG/c²) ρ          [baseline]

    If the deformation energy U(φ) is promoted into the Lagrangian as a
    potential term −V(φ) = −U(φ)/c_φ², the static field equation picks up
    a U'(φ) self-coupling.  Written in terms of an effective "self-coupling
    density"  ρ_U = ρ_nuc (2φ − φ²) = U'(φ)/c², this becomes:

        ∇²φ = −(8πG/c²) [ρ  −  α · ρ_U]             [corrected, Gap-1 direction]

    where α controls the coupling strength:
      •  α = 0   → exact baseline (identical to solve_pm_star)
      •  α = +1  → full Gap-1 correction as written in the formula sheet;
                   U' acts as a *repulsive* self-coupling that reduces the
                   effective gravitating source.  Expect LARGER R and HIGHER
                   M_max compared to baseline.
      •  α = −1  → Gap-3 style: U-field energy ADDS to the gravitating source
                   (as in GR where all energy gravitates).  Expect SMALLER R
                   and reduced M_max.

    Dimensional derivation
    ----------------------
    U'(φ)/c² = ρ_nuc c²(2φ−φ²)/c² = ρ_nuc(2φ−φ²) [kg m⁻³]

    Applying Gauss's theorem to the corrected ∇²φ source in spherical symmetry:

        4πr² (dφ/dr) = −(8πG/c²) m(r) + α (8πG/c²) · 4π ∫₀ʳ ρ_nuc(2φ−φ²) r'² dr'
                     = −(8πG/c²) [m(r) − α m_U(r)]

    where  m_U(r) = 4π ∫₀ʳ ρ_nuc(2φ−φ²) r'² dr'  has the same dimensions as m [kg].
    Hence:

        dφ/dr = −μ_G (m − α m_U) / r²

    The PM force law  a = (c²/2)∇φ  then consistently gives:

        dP/dr = ρ (c²/2) (dφ/dr) = −G m ρ / r²  +  α G m_U ρ / r²

    ODE system (4 variables: m, P, φ, m_U)
    ----------------------------------------
        dm/dr   = 4π r² ρ(P)
        dP/dr   = −G m ρ / r²  +  α G m_U ρ / r²
        dφ/dr   = −μ_G (m − α m_U) / r²
        dm_U/dr = 4π r² ρ_nuc (2φ − φ²)

    Parameters
    ----------
    rho_central : float
        Central density [kg m⁻³].  Must satisfy 0 < ρ_c ≤ ρ_crit.
    alpha : float
        Self-coupling strength.  0 = baseline; +1 = Gap-1 (formula sheet);
        −1 = Gap-3 style (field energy gravitates).  Default +1.
    r_max, n_eval, rtol, atol :
        ODE integration parameters (same meaning as solve_pm_star).

    Returns
    -------
    StarSolution
        Same structure as solve_pm_star.  The phi array reflects the corrected
        field profile.
    """
    if rho_central > RHO_CRIT * (1.0 + 1e-6):
        raise ValueError(
            f"rho_central = {rho_central:.3e} kg/m³ exceeds PM critical density "
            f"ρ_crit = {RHO_CRIT:.3e} kg/m³."
        )

    P_central   = pm_eos_pressure(rho_central)
    phi_central = math.log(rho_central / RHO_NUC)

    def u_prime_normalised(phi: float) -> float:
        """U'(φ)/c² / ρ_nuc = (2φ − φ²)  [dimensionless; multiply by ρ_nuc → kg/m³]."""
        phi_c = max(phi, 0.0)   # avoid negative values outside the star
        return 2.0 * phi_c - phi_c * phi_c

    def rhs(r, y):
        m, P, phi, m_U = y
        rho = pm_eos_density(P)
        up_norm = u_prime_normalised(phi)   # dimensionless ∈ [0, 1] for φ ∈ [0, 1]

        if r < 1.0:
            # Near-centre: m ≈ (4π/3)r³ρ_c and m_U ≈ (4π/3)r³ρ_nuc·up
            # dφ/dr → 0 and dP/dr → 0 by symmetry
            dm_dr   = 4.0 * math.pi * r * r * rho
            dm_U_dr = 4.0 * math.pi * r * r * RHO_NUC * up_norm
            return [dm_dr, 0.0, 0.0, dm_U_dr]

        dm_dr    = 4.0 * math.pi * r * r * rho
        dm_U_dr  = 4.0 * math.pi * r * r * RHO_NUC * up_norm
        m_eff    = m - alpha * m_U
        dphi_dr  = -MU_G * m_eff / (r * r)
        dP_dr    = rho * (c * c / 2.0) * dphi_dr
        return [dm_dr, dP_dr, dphi_dr, dm_U_dr]

    def surface_event(r, y):
        return y[1]   # P = 0

    surface_event.terminal  = True
    surface_event.direction = -1

    r_start    = 1.0
    m_start    = (4.0 / 3.0) * math.pi * r_start**3 * rho_central
    P_start    = P_central
    phi_start  = phi_central
    m_U_start  = (4.0 / 3.0) * math.pi * r_start**3 * RHO_NUC * u_prime_normalised(phi_central)

    r_eval = np.linspace(r_start, r_max, n_eval)

    sol = solve_ivp(
        rhs,
        [r_start, r_max],
        [m_start, P_start, phi_start, m_U_start],
        method='DOP853',
        t_eval=r_eval,
        events=surface_event,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    r_arr   = sol.t
    m_arr   = sol.y[0]
    P_arr   = np.maximum(sol.y[1], 0.0)
    phi_arr = sol.y[2]
    rho_arr = pm_eos_density(P_arr)

    surface_mask = P_arr > 0
    if surface_mask.any():
        i_surf  = np.where(surface_mask)[0][-1]
        R_star  = r_arr[i_surf]
        M_star  = m_arr[i_surf]
        converged = True
    else:
        R_star  = r_max
        M_star  = m_arr[-1]
        converged = False

    return StarSolution(
        M_star=M_star, R_star=R_star,
        rho_central=rho_central, phi_central=phi_central,
        converged=converged,
        r=r_arr, m=m_arr, rho=rho_arr, P=P_arr, phi=phi_arr,
    )


def compute_mr_curve_option_a(
    alpha: float = 1.0,
    n_points: int = 60,
    rho_min_factor: float = 1.01,
    rho_max_factor: float = 1.0,
    **solve_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """M–R curve with the Option-A corrected Poisson self-coupling.

    Parameters
    ----------
    alpha : float
        Self-coupling strength passed to solve_pm_star_option_a.
        0 = baseline; +1 = Gap-1 correction; −1 = Gap-3 style.
    n_points : int
        Number of central-density models.
    rho_min_factor : float
        Lower central density as a multiple of ρ_nuc.
    rho_max_factor : float
        Upper central density as a fraction of ρ_crit.
    **solve_kwargs
        Passed to solve_pm_star_option_a.

    Returns
    -------
    rho_c : ndarray  [kg m⁻³]
    M : ndarray      [M_☉]
    R : ndarray      [km]
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
            star = solve_pm_star_option_a(rho_c, alpha=alpha, **solve_kwargs)
            if star.converged:
                M_arr[i] = star.M_star / M_SUN
                R_arr[i] = star.R_star / 1e3   # m → km
        except (ValueError, Exception):
            pass

    return rho_c_arr, M_arr, R_arr


# ---------------------------------------------------------------------------
# N-field (area-measure) compact-star structure solver
# ---------------------------------------------------------------------------

def solve_pm_star_nfield(
    rho_central: float,
    r_max: float = 5.0e4,
    n_eval: int = 5000,
    rtol: float = 1e-9,
    atol: float = 1e-6,
) -> StarSolution:
    """PM compact-star structure using the n-field (area-measure) Poisson equation.

    Physical motivation
    -------------------
    The standard PM Poisson equation is derived from the coordinate-measure
    action (α=0).  The n-field equation arises from the area-measure action:

        S₂ = ∫ ½|∇φ|² n² d³x

    whose Euler-Lagrange vacuum equation is ∇²φ + |∇φ|² = 0, linearising to
    ∇²n = 0  via  n = e^φ.  With the standard matter source the sourced equation
    is:

        ∇²φ + |∇φ|² = −(8πG/c²) ρ_nuc n

    which in the n-variable reads:

        ∇²n = −(8πG ρ_nuc / c²) n²

    The n² source replaces the standard n source — a stronger coupling at high
    density (n > 1) that is invisible in the weak field (n → 1).

    Structure system
    ----------------
    Define  m̃_n(r) = ∫₀ʳ 4π r'² n(r')² dr'  (n²-accumulator, units m³).

    Applying Gauss's theorem to ∇²n = −(8πGρ_nuc/c²)n²:

        dn/dr = −(2G ρ_nuc / c²) m̃_n / r²

    The PM force law a = (c²/2)∇φ = (c²/2n)∇n gives the hydrostatic equation:

        dP/dr = ρ a = ρ_nuc n · (c²/2n) · (dn/dr) = (c²ρ_nuc/2) dn/dr

    This is already the exact differential of the PM EOS  P = (c²ρ_nuc/2)(n−1).
    The n-field EOS and the hydrostatic equation are thus **automatically consistent**
    — the inconsistency that existed in the standard formulation does not appear here.

    The 3-variable ODE state is  y = [m_phys, m̃_n, n]:

        dm_phys/dr = 4π r² ρ_nuc n       (physical mass, for output)
        dm̃_n/dr   = 4π r² n²             (Poisson accumulator)
        dn/dr      = −(2G ρ_nuc/c²) m̃_n / r²

    Surface: n = 1  (P = 0, ρ = ρ_nuc).

    Auto-consistency check
    ----------------------
    Because P = (c²ρ_nuc/2)(n−1) and dP/dr = (c²ρ_nuc/2)(dn/dr), the φ profile
    computed as φ = ln n is exactly consistent with the EOS φ = ln(ρ/ρ_nuc) at
    every grid point — no accumulation of EOS/Poisson drift.

    Parameters
    ----------
    rho_central : float
        Central density [kg m⁻³].  Must satisfy 0 < ρ_c ≤ ρ_crit.
    r_max : float
        Maximum integration radius [m].  Default 50 km.
    n_eval : int
        Number of output evaluation points.
    rtol, atol : float
        ODE integration tolerances.

    Returns
    -------
    StarSolution
        Same structure as solve_pm_star.  phi = ln n is exact by construction.
    """
    if rho_central > RHO_CRIT * (1.0 + 1e-6):
        raise ValueError(
            f"rho_central = {rho_central:.3e} kg/m³ exceeds PM critical density "
            f"ρ_crit = {RHO_CRIT:.3e} kg/m³.  Stable PM matter configurations "
            f"require ρ_c ≤ ρ_crit."
        )

    n_central   = rho_central / RHO_NUC      # n_c = ρ_c / ρ_nuc = e^{φ_c}
    phi_central = math.log(n_central)         # φ_c = ln n_c

    # Coefficient: 2G ρ_nuc / c²  [m⁻²]  (same dimensional role as μ_G = 2G/c²
    # in the standard solver, but weighted by ρ_nuc because m̃_n has units m³
    # rather than kg)
    kappa_n = 2.0 * G * RHO_NUC / (c * c)

    def rhs(r, y):
        m_phys, m_n, n = y
        if r < 1.0:
            # Near centre: all gradients → 0 by spherical symmetry
            dm_phys = 4.0 * math.pi * r * r * RHO_NUC * n
            dm_n    = 4.0 * math.pi * r * r * n * n
            return [dm_phys, dm_n, 0.0]
        dm_phys = 4.0 * math.pi * r * r * RHO_NUC * n
        dm_n    = 4.0 * math.pi * r * r * n * n
        dn_dr   = -kappa_n * m_n / (r * r)
        return [dm_phys, dm_n, dn_dr]

    # Surface event: n = 1  (equivalent to P = 0)
    def surface_event(r, y):
        return y[2] - 1.0   # n − 1 = 0

    surface_event.terminal  = True
    surface_event.direction = -1   # n decreasing toward 1

    # Near-centre series: m ≈ (4π/3)r³ρ_nuc n_c,  m̃_n ≈ (4π/3)r³ n_c²
    r_start    = 1.0
    m_phys_0   = (4.0 / 3.0) * math.pi * r_start**3 * RHO_NUC * n_central
    m_n_0      = (4.0 / 3.0) * math.pi * r_start**3 * n_central**2
    n_0        = n_central

    r_eval = np.linspace(r_start, r_max, n_eval)

    sol = solve_ivp(
        rhs,
        [r_start, r_max],
        [m_phys_0, m_n_0, n_0],
        method='DOP853',
        t_eval=r_eval,
        events=surface_event,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    r_arr      = sol.t
    m_arr      = sol.y[0]
    n_arr      = np.maximum(sol.y[2], 1.0)    # clamp n ≥ 1 (no vacuum below nuclear density)
    phi_arr    = np.log(n_arr)                 # exact by construction
    P_arr      = np.maximum(
        c * c * RHO_NUC / 2.0 * (n_arr - 1.0), 0.0
    )
    rho_arr    = RHO_NUC * n_arr

    surface_mask = n_arr > 1.0 + 1e-12
    if surface_mask.any():
        i_surf    = np.where(surface_mask)[0][-1]
        R_star    = r_arr[i_surf]
        M_star    = m_arr[i_surf]
        converged = True
    else:
        R_star    = r_max
        M_star    = m_arr[-1]
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


def compute_mr_curve_nfield(
    n_points: int = 60,
    rho_min_factor: float = 1.01,
    rho_max_factor: float = 1.0,
    **solve_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """M–R curve using the n-field (area-measure) stellar structure solver.

    Parameters
    ----------
    n_points : int
        Number of central-density models.
    rho_min_factor : float
        Lower central density as a multiple of ρ_nuc.
    rho_max_factor : float
        Upper central density as a fraction of ρ_crit.
    **solve_kwargs
        Passed to solve_pm_star_nfield.

    Returns
    -------
    rho_c : ndarray  [kg m⁻³]
    M : ndarray      [M_☉]
    R : ndarray      [km]
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
            star = solve_pm_star_nfield(rho_c, **solve_kwargs)
            if star.converged:
                M_arr[i] = star.M_star / M_SUN
                R_arr[i] = star.R_star / 1e3    # m → km
        except (ValueError, Exception):
            pass

    return rho_c_arr, M_arr, R_arr


# ---------------------------------------------------------------------------
# N-field + stiffened potential stellar structure solver  (action-derived)
# ---------------------------------------------------------------------------

def solve_pm_star_nfield_stiffened(
    rho_central: float,
    r_max: float = 5.0e4,
    n_eval: int = 5000,
    rtol: float = 1e-9,
    atol: float = 1e-6,
) -> StarSolution:
    """PM compact-star with the action-derived self-stiffening potential.

    Action and vacuum-subtracted potential
    ----------------------------------------
    The full n-field action with the medium self-energy potential is:

        S_full = ∫ [½|∇φ|² n²  −  A · V̂(n)] d³x

    where  V̂(n) = ½n² + ln n  and  A = (8πG/c²) ρ_nuc = κ ρ_nuc.

    V̂'(n) = n + 1/n, but V̂'(1) = 2 ≠ 0, meaning the bare potential provides
    a non-zero source even at the vacuum boundary n = 1.  For stellar structure
    we require the surface (n = 1) to be a stable fixed point of the field
    equation, so we use the vacuum-subtracted potential:

        V̂_vac(n) = V̂(n) − V̂(1) − V̂'(1)(n − 1)
                  = ½n² + ln n − 3/2 − 2(n − 1)

    which satisfies  V̂_vac(1) = 0  and  V̂_vac'(1) = 0.

        V̂_vac'(n) = n + 1/n − 2 = (n − 1)² / n  ≥ 0  (for n ≥ 1)

    This is the "above-vacuum" stiffening: it vanishes exactly at the phase
    boundary (n = 1), is strictly positive inside the medium (n > 1), and is
    maximally small near the surface, growing only quadratically in (n − 1).

    Euler-Lagrange field equation (exact, vacuum-subtracted)
    --------------------------------------------------------
        ∇²n = κ ρ_nuc (n − 1)²/n  −  κ ρ_nuc n²

    The first term (self-stiffening) reduces the effective inward pull; the
    second is the full n-field gravity source.  At n = 1 (surface):
        ∇²n = 0 − κ ρ_nuc = −κ ρ_nuc  [same as bare n-field],
    so the surface physics is unchanged.  Inside the medium the stiffening
    provides a positive correction to n at each r.

    At n = 2: stiffening = (1/2) = 0.5 / gravity = 4  →  12 % reduction.
    At n = 2.72: stiffening ≈ 1.09 / gravity ≈ 7.40  →  15 % reduction.

    ODE state  y = [m_phys, m̃_n, m̃_stiff, n]
    ------------------------------------------
        dm_phys/dr    = 4π r² ρ_nuc n
        dm̃_n/dr      = 4π r² n²                   (gravity accumulator)
        dm̃_stiff/dr  = 4π r² (n − 1)²/n           (vac-subtracted stiffening)
        dn/dr         = −κ_n (m̃_n − m̃_stiff) / r²

    where  κ_n = 2G ρ_nuc / c²  (same coefficient for both terms).

    Surface: n = 1  (P = 0).

    Parameters
    ----------
    rho_central : float
        Central density [kg m⁻³].  Must satisfy 0 < ρ_c ≤ ρ_crit.

    Returns
    -------
    StarSolution
    """
    if rho_central > RHO_CRIT * (1.0 + 1e-6):
        raise ValueError(
            f"rho_central = {rho_central:.3e} kg/m³ exceeds PM critical density "
            f"ρ_crit = {RHO_CRIT:.3e} kg/m³."
        )

    n_central   = rho_central / RHO_NUC
    phi_central = math.log(n_central)

    # Both n-field and stiffening terms share the same coefficient
    #   κ_n = 2G ρ_nuc / c²  (this is the Gauss-theorem form: κ/4π × ρ_nuc)
    kappa_n = 2.0 * G * RHO_NUC / (c * c)

    def rhs(r, y):
        m_phys, m_n, m_stiff, n_val = y
        n_s = max(n_val, 1.0)   # clamp at surface

        # Vacuum-subtracted stiffening source: (n-1)²/n = n + 1/n - 2
        stiff_src = (n_s - 1.0) ** 2 / n_s

        if r < 1.0:
            dm_phys   = 4.0 * math.pi * r * r * RHO_NUC * n_s
            dm_n      = 4.0 * math.pi * r * r * n_s * n_s
            dm_stiff  = 4.0 * math.pi * r * r * stiff_src
            return [dm_phys, dm_n, dm_stiff, 0.0]

        dm_phys  = 4.0 * math.pi * r * r * RHO_NUC * n_s
        dm_n     = 4.0 * math.pi * r * r * n_s * n_s
        dm_stiff = 4.0 * math.pi * r * r * stiff_src
        dn_dr    = -kappa_n * (m_n - m_stiff) / (r * r)
        return [dm_phys, dm_n, dm_stiff, dn_dr]

    def surface_event(r, y):
        return y[3] - 1.0   # n = 1

    surface_event.terminal  = True
    surface_event.direction = -1

    r_start      = 1.0
    m_phys_0     = (4.0 / 3.0) * math.pi * r_start**3 * RHO_NUC * n_central
    m_n_0        = (4.0 / 3.0) * math.pi * r_start**3 * n_central**2
    stiff_src_0  = (n_central - 1.0) ** 2 / n_central   # vacuum-subtracted
    m_stiff_0    = (4.0 / 3.0) * math.pi * r_start**3 * stiff_src_0
    n_0          = n_central

    r_eval = np.linspace(r_start, r_max, n_eval)

    sol = solve_ivp(
        rhs,
        [r_start, r_max],
        [m_phys_0, m_n_0, m_stiff_0, n_0],
        method='DOP853',
        t_eval=r_eval,
        events=surface_event,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    r_arr   = sol.t
    m_arr   = sol.y[0]
    n_arr   = np.maximum(sol.y[3], 1.0)
    phi_arr = np.log(n_arr)
    P_arr   = np.maximum(c * c * RHO_NUC / 2.0 * (n_arr - 1.0), 0.0)
    rho_arr = RHO_NUC * n_arr

    surface_mask = n_arr > 1.0 + 1e-12
    if surface_mask.any():
        i_surf    = np.where(surface_mask)[0][-1]
        R_star    = r_arr[i_surf]
        M_star    = m_arr[i_surf]
        converged = True
    else:
        R_star    = r_max
        M_star    = m_arr[-1]
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


def compute_mr_curve_nfield_stiffened(
    n_points: int = 60,
    rho_min_factor: float = 1.01,
    rho_max_factor: float = 1.0,
    **solve_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """M–R curve for the action-derived stiffened n-field solver.

    Returns
    -------
    rho_c : ndarray  [kg m⁻³]
    M : ndarray      [M_☉]
    R : ndarray      [km]
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
            star = solve_pm_star_nfield_stiffened(rho_c, **solve_kwargs)
            if star.converged:
                M_arr[i] = star.M_star / M_SUN
                R_arr[i] = star.R_star / 1e3
        except (ValueError, Exception):
            pass

    return rho_c_arr, M_arr, R_arr


# ---------------------------------------------------------------------------
# α=3 physical-measure stellar structure solver
# ---------------------------------------------------------------------------

def solve_pm_star_physical_measure(
    rho_central: float,
    r_max: float = 5.0e4,
    n_eval: int = 5000,
    rtol: float = 1e-9,
    atol: float = 1e-6,
) -> StarSolution:
    """PM compact-star structure using the physical-measure (α=3) Poisson equation.

    Physical motivation
    -------------------
    The physical-measure action uses the volume element n³ d³x:

        S₃ = ∫ ½|∇φ|² n³ d³x

    Euler-Lagrange vacuum equation:

        ∇²φ + (3/2)|∇φ|² = 0

    Linearisation via w = n^{3/2} = e^{3φ/2}:

        ∇²w = 0    (linear in w, exact superposition)

    Exact point-mass vacuum solution:

        w(r) = 1 + 3GM/(c²r)  →  φ(r) = (2/3) ln(1 + 3GM/c²r)

    Force law (from a = (c²/3) ∇w / w):

        a(r) = −GM/r² / (1 + 3GM/c²r)

    Phase boundary (φ=1):

        r_phase = 3GM / [c²(e^{3/2} − 1)] ≈ 0.431 r_s

    sourced equation with matter:

        ∇²w = −K_w w^{1/3}    where K_w = 12πGρ_nuc/c²

    In the weak field (w≈1, n≈1) this reduces to ∇²w ≈ −K_w  which matches
    the Newton–Poisson source at the same rate as the α=0 and α=2 solvers.

    Structure system (w-variable ODE)
    ----------------------------------
    Define  m̃_w(r) = ∫₀ʳ 4π r'² w(r')^{1/3} dr'  (w^{1/3} accumulator).

    Gauss's theorem on ∇²w = −κ_w w^{4/3}  with  κ_w = (3/2)(8πGρ_nuc/c²):

        dw/dr = −(κ_w / (4πr²)) × 4πr² w^{4/3}  [can't factor simply]

    More directly: the sourced equation in w-form is:

        ∇²w = −κ_w w^{4/3}     κ_w = 12πGρ_nuc / c²

    Gauss integral: dw/dr = −K_w/(4πr²) ∫₀ʳ 4πr'² w^{1/3} dr'
                          = −κ_w m̃_w / r²    where κ_w = K_w/4π = 3Gρ_nuc/c²

    ODE state y = [m_phys, m̃_w, w]:

        dm_phys/dr = 4π r² ρ_nuc n = 4π r² ρ_nuc w^{2/3}
        dm̃_w/dr   = 4π r² w^{1/3}
        dw/dr      = −κ_w m̃_w / r²

    Surface: n = 1 → w = n^{3/2} = 1  (P = 0, ρ = ρ_nuc).

    Parameters
    ----------
    rho_central : float
        Central density [kg m⁻³].  Must satisfy 0 < ρ_c ≤ ρ_crit.
    r_max : float
        Maximum integration radius [m].
    n_eval : int
        Number of output evaluation points.
    rtol, atol : float
        ODE tolerances.

    Returns
    -------
    StarSolution
        Same structure as solve_pm_star.  phi = (2/3) ln w by construction.
    """
    if rho_central > RHO_CRIT * (1.0 + 1e-6):
        raise ValueError(
            f"rho_central = {rho_central:.3e} kg/m³ exceeds PM critical density "
            f"ρ_crit = {RHO_CRIT:.3e} kg/m³."
        )

    n_central   = rho_central / RHO_NUC
    phi_central = math.log(n_central)
    w_central   = n_central ** 1.5           # w = n^{3/2}

    # κ_w = 3Gρ_nuc / c²  (= K_w / 4π where K_w = 12πGρ_nuc/c² is the
    # coefficient in ∇²w = −K_w w^{1/3}; same Gauss-theorem convention as
    # kappa_n = 2Gρ_nuc/c² for the n-field solver)
    kappa_w = 3.0 * G * RHO_NUC / (c * c)

    def rhs(r, y):
        m_phys, m_w, w = y
        # n = w^{2/3}, clamp w ≥ 1 to avoid sub-nuclear densities
        w_safe = max(w, 1.0)
        n      = w_safe ** (2.0 / 3.0)
        if r < 1.0:
            dm_phys = 4.0 * math.pi * r * r * RHO_NUC * n
            dm_w    = 4.0 * math.pi * r * r * w_safe ** (1.0 / 3.0)
            return [dm_phys, dm_w, 0.0]
        dm_phys = 4.0 * math.pi * r * r * RHO_NUC * n
        dm_w    = 4.0 * math.pi * r * r * w_safe ** (1.0 / 3.0)
        dw_dr   = -kappa_w * m_w / (r * r)
        return [dm_phys, dm_w, dw_dr]

    def surface_event(r, y):
        return y[2] - 1.0   # w = n^{3/2} = 1  ↔  n = 1

    surface_event.terminal  = True
    surface_event.direction = -1

    r_start  = 1.0
    m_phys_0 = (4.0 / 3.0) * math.pi * r_start**3 * RHO_NUC * n_central
    m_w_0    = (4.0 / 3.0) * math.pi * r_start**3 * w_central ** (1.0 / 3.0)
    w_0      = w_central

    r_eval = np.linspace(r_start, r_max, n_eval)

    sol = solve_ivp(
        rhs,
        [r_start, r_max],
        [m_phys_0, m_w_0, w_0],
        method='DOP853',
        t_eval=r_eval,
        events=surface_event,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    r_arr   = sol.t
    m_arr   = sol.y[0]
    w_arr   = np.maximum(sol.y[2], 1.0)
    n_arr   = w_arr ** (2.0 / 3.0)
    phi_arr = (2.0 / 3.0) * np.log(w_arr)   # φ = (2/3) ln w, exact
    P_arr   = np.maximum(c * c * RHO_NUC / 2.0 * (n_arr - 1.0), 0.0)
    rho_arr = RHO_NUC * n_arr

    surface_mask = n_arr > 1.0 + 1e-12
    if surface_mask.any():
        i_surf    = np.where(surface_mask)[0][-1]
        R_star    = r_arr[i_surf]
        M_star    = m_arr[i_surf]
        converged = True
    else:
        R_star    = r_max
        M_star    = m_arr[-1]
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


def compute_mr_curve_physical_measure(
    n_points: int = 60,
    rho_min_factor: float = 1.01,
    rho_max_factor: float = 1.0,
    **solve_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """M–R curve using the physical-measure (α=3) stellar structure solver.

    Returns
    -------
    rho_c : ndarray  [kg m⁻³]
    M : ndarray      [M_☉]
    R : ndarray      [km]
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
            star = solve_pm_star_physical_measure(rho_c, **solve_kwargs)
            if star.converged:
                M_arr[i] = star.M_star / M_SUN
                R_arr[i] = star.R_star / 1e3
        except (ValueError, Exception):
            pass

    return rho_c_arr, M_arr, R_arr


# ---------------------------------------------------------------------------
# Surface gravitational redshift
# ---------------------------------------------------------------------------

def pm_surface_redshift(M_star: float, R_star: float) -> float:
    """PM surface gravitational redshift (exact in PM, all compactness).

    Derivation
    ----------
    The PM optical metric has g_tt = −c²/n² = −c² e^{−2φ}.
    For a static source at radius R_star and a receiver at infinity:

        1 + z = ν_emit/ν_obs = sqrt(−g_tt(∞) / −g_tt(R))
                             = sqrt(c² / c² e^{−2φ_R}) = e^{φ_R}

    The relevant φ_R is the EXTERIOR field at the stellar surface:

        φ_surface = μ_G M_star / R_star   (exact exterior Poisson solution)

    This is different from the interior ODE value (which equals ln(ρ_nuc/ρ_nuc)=0
    at the topmost cell, reflecting the EOS boundary condition).

    Result:
        z_PM = e^{μ_G M_star / R_star} − 1 = e^{2GM/(c²R)} − 1

    PM vs GR comparison
    --------------------
    At leading order in x = 2GM/c²R:
      z_PM ≈ x + x²/2 + ...          → x leading term
      z_GR = 1/√(1−x) − 1 ≈ x/2 + ...  → x/2 leading term

    PM predicts roughly TWICE the GR surface redshift for the same (M, R).
    This is a genuine, observationally testable difference: NICER spectral
    observations can constrain z_surface independently of mass.

    Parameters
    ----------
    M_star : float
        Total stellar mass [kg].
    R_star : float
        Stellar radius [m].

    Returns
    -------
    float
        Dimensionless surface gravitational redshift z = Δν/ν.
    """
    phi_surface = MU_G * M_star / R_star
    return math.exp(phi_surface) - 1.0


def compute_surface_redshift_track(
    n_points: int = 60,
    rho_min_factor: float = 1.01,
    rho_max_factor: float = 1.0,
    **solve_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute PM surface redshift along the M–R track.

    Returns
    -------
    rho_c : ndarray  [kg m⁻³]
    M : ndarray  [M_☉]
    R : ndarray  [km]
    z_surf : ndarray  PM surface gravitational redshift z = e^{μ_G M/R} − 1
    """
    rho_c_arr, M_arr, R_arr = compute_mr_curve(
        n_points=n_points,
        rho_min_factor=rho_min_factor,
        rho_max_factor=rho_max_factor,
        **solve_kwargs,
    )

    z_arr = np.full(n_points, np.nan)
    for i in range(n_points):
        if np.isfinite(M_arr[i]) and np.isfinite(R_arr[i]):
            M_si = M_arr[i] * M_SUN
            R_si = R_arr[i] * 1e3   # km → m
            z_arr[i] = pm_surface_redshift(M_si, R_si)

    return rho_c_arr, M_arr, R_arr, z_arr
