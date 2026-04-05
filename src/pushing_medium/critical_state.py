"""PM Critical State: matter/energy transition point.

In the Pushing-Medium framework matter is a stable localised configuration of
the PM medium held together by pressure balance around a compression gradient.
The configuration remains stable as long as the curvature of the deformation-
energy landscape U(φ) is positive (∂²U/∂φ² > 0).  The critical compression
φ_crit is the smallest positive root of U''(φ) = 0.

----------------------------------------------------------------------
Deformation energy (Landau-Ginzburg form for the PM scalar field)
----------------------------------------------------------------------

    U(φ) = ε₀ [ φ² − φ³/3 ]

    U'(φ) = ε₀ [ 2φ − φ² ]

    U''(φ) = 2 ε₀ ( 1 − φ )

where  ε₀ = ρ_nuc c²  with  ρ_nuc  the nuclear saturation density.

Physical origin of the cubic term:
  The PM Lagrangian contains the coupling φ ∇·(φ u).  Adiabatically
  eliminating the flow field u (driven by ∇φ in the static limit) yields an
  effective cubic self-interaction of the scalar field, the leading
  non-linear correction to the quadratic elastic energy φ².  The result is
  the classic Landau free energy for a second-order stability transition.

Stability boundaries:
  U''(φ) > 0  ↔  φ < 1  – stable matter configuration
  U''(φ) = 0  ↔  φ = 1  – critical point (onset of instability)
  U''(φ) < 0  ↔  φ > 1  – unstable; configuration dissolves into excitations

Reference:
  docs/LLM_notes_and_conversations/matter_energy_transition.md

----------------------------------------------------------------------
Density and pressure at the critical state
----------------------------------------------------------------------

PM density-from-φ mapping:

    ρ(φ) = ρ_nuc · exp(φ) = ρ_nuc · n

The PM refractive index n = exp(φ) directly amplifies the medium density;
ρ_nuc is the reference (nuclear saturation) density at which the PM medium
is in its baseline compressed state for stellar matter.

PM equation of state (linear compressible EOS):

    P(φ) = (c²/2) [ρ(φ) − ρ_nuc]

This is the PM-native EOS P = (c²/2)(ρ − ρ_nuc), derived from the PM
field action.  At the critical compression φ_crit = 1:
  n_crit  = e          ≈ 2.718
  ρ_crit  = ρ_nuc e   ≈ 6.25 × 10¹⁷ kg m⁻³  (above nuclear density ✓)
  P_crit  = (c²/2) ρ_nuc (e − 1)  ≈ 1.78 × 10³⁴ Pa   (finite, positive ✓)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Physical constants (SI) – matching src/pushing_medium/core.py
# ---------------------------------------------------------------------------
G = 6.67430e-11       # Gravitational constant (N m² kg⁻²)
c = 299792458.0       # Speed of light (m s⁻¹)

# Nuclear saturation density: reference density for PM matter configurations.
# This is the density at which the PM medium is in its baseline compressed
# state for nuclear matter (the onset of quark-gluon interplay).
RHO_NUC: float = 2.3e17   # kg m⁻³


# ---------------------------------------------------------------------------
# Deformation energy helpers
# ---------------------------------------------------------------------------

def pm_deformation_energy(phi: float, rho_ref: float = RHO_NUC) -> float:
    """PM deformation energy density  U(φ) = ε₀ [φ² − φ³/3].

    Parameters
    ----------
    phi : float
        PM compression field  φ = ln n  (dimensionless, ≥ 0 for compressed
        medium).
    rho_ref : float
        Reference energy density scale  ε₀ = rho_ref × c²  [J m⁻³].
        Defaults to nuclear saturation density.

    Returns
    -------
    float
        Deformation energy density  U(φ)  [J m⁻³].
    """
    eps0 = rho_ref * c * c
    return eps0 * (phi * phi - phi * phi * phi / 3.0)


def pm_deformation_energy_deriv1(phi: float, rho_ref: float = RHO_NUC) -> float:
    """First derivative  U'(φ) = ε₀ [2φ − φ²].

    Equals the PM medium pressure (up to a sign convention for the compression
    variable): P ~ U'(φ).  Positive for 0 < φ < 2.

    Parameters
    ----------
    phi : float
        PM compression field.
    rho_ref : float
        Reference energy density scale.

    Returns
    -------
    float
        U'(φ)  [J m⁻³].
    """
    eps0 = rho_ref * c * c
    return eps0 * (2.0 * phi - phi * phi)


def pm_deformation_energy_deriv2(phi: float, rho_ref: float = RHO_NUC) -> float:
    """Second derivative  U''(φ) = 2 ε₀ (1 − φ).

    Stability:
      U'' > 0  ↔  φ < 1  – stable matter configuration.
      U'' = 0  ↔  φ = 1  – critical point.
      U'' < 0  ↔  φ > 1  – unstable (matter dissolves).

    Parameters
    ----------
    phi : float
        PM compression field.
    rho_ref : float
        Reference energy density scale.

    Returns
    -------
    float
        U''(φ)  [J m⁻³].
    """
    eps0 = rho_ref * c * c
    return 2.0 * eps0 * (1.0 - phi)


# ---------------------------------------------------------------------------
# PM equation-of-state helpers
# ---------------------------------------------------------------------------

def pm_density_from_phi(phi: float, rho_ref: float = RHO_NUC) -> float:
    """PM density-from-φ mapping:  ρ(φ) = ρ_ref · exp(φ).

    The PM refractive index n = exp(φ) directly amplifies the medium density.
    At φ = 0 (vacuum compression) the density equals the reference density
    ρ_ref; at the critical compression φ = 1 it equals ρ_ref · e.

    Parameters
    ----------
    phi : float
        PM compression field φ = ln n.
    rho_ref : float
        Reference density [kg m⁻³].  Defaults to nuclear saturation density.

    Returns
    -------
    float
        Local matter density  ρ(φ)  [kg m⁻³].
    """
    return rho_ref * math.exp(phi)


def pm_pressure_from_phi(phi: float, rho_ref: float = RHO_NUC) -> float:
    """PM pressure law:  P(φ) = (c²/2) [ρ(φ) − ρ_ref].

    PM-native linear compressible EOS derived from the PM field action.
    The sound speed c_s = c/√2 follows directly from dP/dρ = c²/2.

    Parameters
    ----------
    phi : float
        PM compression field φ = ln n.
    rho_ref : float
        Reference (nuclear saturation) density [kg m⁻³].

    Returns
    -------
    float
        Pressure  P(φ)  [Pa].
    """
    rho = pm_density_from_phi(phi, rho_ref=rho_ref)
    return 0.5 * c * c * (rho - rho_ref)


# ---------------------------------------------------------------------------
# CriticalState dataclass
# ---------------------------------------------------------------------------

@dataclass
class CriticalState:
    """Critical state of a PM matter configuration at the matter/energy transition.

    Attributes
    ----------
    phi_crit : float
        Critical compression  φ_crit  where  U''(φ_crit) = 0  (dimensionless).
    rho_crit : float
        Critical density  ρ_crit = ρ_nuc exp(φ_crit)  [kg m⁻³].
    p_crit : float
        Critical pressure  P_crit = (c²/2)(ρ_crit − ρ_nuc)  [Pa].
    """
    phi_crit: float
    rho_crit: float
    p_crit: float


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

def compute_critical_state(
    phi_min: float = 0.0,
    phi_max: float = 20.0,
    n_scan: int = 1000,
    rho_ref: float = RHO_NUC,
) -> CriticalState:
    """Find the PM critical compression where matter configurations become unstable.

    Algorithm
    ---------
    1. Scan  U''(φ)  on a uniform grid over [phi_min, phi_max].
    2. Locate the first sub-interval  [φ_a, φ_b]  where  U''  changes sign
       (U''(φ_a) · U''(φ_b) ≤ 0).
    3. Refine the root with 64 bisection steps (≈ machine-precision accuracy
       for smooth functions).
    4. Compute  ρ_crit  and  P_crit  via the PM density-from-φ mapping and
       the PM pressure law.

    Parameters
    ----------
    phi_min : float
        Lower end of the φ scan range.  Default 0.0 (vacuum compression).
    phi_max : float
        Upper end of the φ scan range.  Default 20.0 (well beyond any
        physically realised compression in known compact objects).
    n_scan : int
        Number of uniform scan steps.  Default 1000.
    rho_ref : float
        Reference density for the deformation energy and EOS helpers.
        Defaults to nuclear saturation density RHO_NUC.

    Returns
    -------
    CriticalState
        Dataclass containing  phi_crit, rho_crit, p_crit.

    Raises
    ------
    RuntimeError
        If no root of  U''(φ)  is found in  [phi_min, phi_max].
    """
    dphi = (phi_max - phi_min) / n_scan
    phi_root: float | None = None

    phi_a = phi_min
    U_pp_a = pm_deformation_energy_deriv2(phi_a, rho_ref=rho_ref)

    for i in range(n_scan):
        phi_b = phi_min + (i + 1) * dphi
        U_pp_b = pm_deformation_energy_deriv2(phi_b, rho_ref=rho_ref)

        if U_pp_a * U_pp_b <= 0.0:   # sign change → root bracketed
            # ---- bisection refinement ----
            lo, hi = phi_a, phi_b
            for _ in range(64):
                mid = 0.5 * (lo + hi)
                if (pm_deformation_energy_deriv2(lo, rho_ref=rho_ref)
                        * pm_deformation_energy_deriv2(mid, rho_ref=rho_ref)) <= 0.0:
                    hi = mid
                else:
                    lo = mid
            phi_root = 0.5 * (lo + hi)
            break

        phi_a = phi_b
        U_pp_a = U_pp_b

    if phi_root is None:
        raise RuntimeError(
            f"No root of U''(φ) found in [{phi_min}, {phi_max}]. "
            "Try extending phi_max or increasing n_scan."
        )

    phi_crit = phi_root
    rho_crit = pm_density_from_phi(phi_crit, rho_ref=rho_ref)
    p_crit = pm_pressure_from_phi(phi_crit, rho_ref=rho_ref)

    return CriticalState(
        phi_crit=phi_crit,
        rho_crit=rho_crit,
        p_crit=p_crit,
    )


# ---------------------------------------------------------------------------
# Perturbation stability mass-squared (dynamical stability condition)
# ---------------------------------------------------------------------------

def pm_phase_stability_mass_sq(
    phi: float,
    c_phi: float = c,
    rho_ref: float = RHO_NUC,
) -> float:
    """Perturbation mass-squared for the PM medium phase stability.

    Derivation
    ----------
    The area-measure action with deformation potential:

        S = ∫ [½|∇φ|² n² − U(φ)/c_φ²] d³x

    gives the Euler-Lagrange field equation (vacuum + source term):

        ∇²φ + |∇φ|² = −U'(φ) / (c_φ² n²) − (8πGρ_nuc/c²) n

    Linearising around a uniform background φ₀ (no spatial gradients):

        ∂²(δφ)/∂t² = c_φ² [∇²(δφ) − m²_eff(φ₀) δφ]

    where the effective perturbation mass-squared is:

        m²_eff(φ₀) = U''(φ₀) / (c_φ² n₀²) = 2ε₀(1−φ₀) / (c_φ² n₀²)

    Sign interpretation:
      m²_eff > 0  (φ < 1):  massive mode — perturbations oscillate and damp
                            (stable matter phase).
      m²_eff = 0  (φ = 1):  massless mode — perturbations are long-range,
                            marginal stability (second-order phase transition).
      m²_eff < 0  (φ > 1):  tachyonic — long-wavelength perturbations grow
                            exponentially (matter phase dissolves to energy phase).

    Consequence
    -----------
    The stability cap φ < 1 is not an external postulate — it is the condition
    m²_eff(φ) ≥ 0, which follows from U''(φ) ≥ 0.  The deformation potential
    U(φ) encodes the phase transition: its Hessian goes soft (m² → 0) at the
    critical compression.  Above that compression the matter configuration is
    dynamically unstable and the medium flows to the energy phase.

    The factor U''(φ)/n² reflects the same softening seen in the sound-speed
    comparison table (numerical probe: c²_U/c²_EOS = 4(1−φ)/n, decreasing
    from 4 at φ=0 to 0 at φ=1).

    Parameters
    ----------
    phi : float
        PM compression field.
    c_phi : float
        Wave speed for φ-mode perturbations [m s⁻¹].  Default = c (the speed
        of light).  Use c/√2 to match the EOS sound speed instead.
    rho_ref : float
        Reference density for ε₀ = rho_ref × c².  Default = RHO_NUC.

    Returns
    -------
    float
        m²_eff(φ)  [kg m⁻³].  Positive → stable; zero → marginal; negative → unstable.

        Dimensional analysis: [U''] = J/m³ = kg/(m·s²),
        [c_phi²] = m²/s², [n²] = dimensionless → [m²_eff] = kg/m³.
    """
    n = math.exp(phi)
    u_double_prime = pm_deformation_energy_deriv2(phi, rho_ref=rho_ref)   # = 2ε₀(1−φ)
    return u_double_prime / (c_phi * c_phi * n * n)

