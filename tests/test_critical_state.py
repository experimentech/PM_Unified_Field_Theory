"""Test suite for the PM critical state (matter/energy transition).

Tests
-----
1. φ_crit is positive and finite.
2. ρ_crit is above nuclear density scales (> 1e17 kg m⁻³).
3. P_crit is finite and positive.
4. U''(φ_crit) ≈ 0 to numerical precision.
5. ρ_crit and P_crit are self-consistent with the PM EOS helpers.
6. U'' changes sign across φ_crit (confirms instability threshold).
7. Print helper for manual inspection (captured by pytest -s).
"""

import math
import pytest

from pushing_medium.critical_state import (
    CriticalState,
    compute_critical_state,
    pm_deformation_energy,
    pm_deformation_energy_deriv1,
    pm_deformation_energy_deriv2,
    pm_density_from_phi,
    pm_pressure_from_phi,
    RHO_NUC,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def critical_state() -> CriticalState:
    """Compute the PM critical state once per test module."""
    return compute_critical_state()


# ---------------------------------------------------------------------------
# Core assertions
# ---------------------------------------------------------------------------

def test_phi_crit_positive(critical_state: CriticalState) -> None:
    """φ_crit must be strictly positive (vacuum φ = 0 is always stable)."""
    assert critical_state.phi_crit > 0, (
        f"Expected φ_crit > 0, got φ_crit = {critical_state.phi_crit}"
    )


def test_phi_crit_finite(critical_state: CriticalState) -> None:
    """φ_crit must be a finite number (not inf or nan)."""
    assert math.isfinite(critical_state.phi_crit), (
        f"φ_crit is not finite: {critical_state.phi_crit}"
    )


def test_rho_crit_above_nuclear_density(critical_state: CriticalState) -> None:
    """ρ_crit must exceed the nuclear density scale (1e17 kg m⁻³).

    At nuclear density the PM medium enters the strongly compressed regime.
    The critical compression must lie at or above this scale.
    """
    assert critical_state.rho_crit > 1e17, (
        f"ρ_crit = {critical_state.rho_crit:.3e} kg/m³ must exceed 1e17 kg/m³ "
        f"(nuclear density scale)"
    )


def test_p_crit_positive(critical_state: CriticalState) -> None:
    """P_crit must be strictly positive."""
    assert critical_state.p_crit > 0.0, (
        f"Expected P_crit > 0, got {critical_state.p_crit}"
    )


def test_p_crit_finite(critical_state: CriticalState) -> None:
    """P_crit must be finite."""
    assert math.isfinite(critical_state.p_crit), (
        f"P_crit is not finite: {critical_state.p_crit}"
    )


def test_p_crit_numerical_value(critical_state: CriticalState) -> None:
    """P_crit must match the PM-EOS analytic value (c²/2)ρ_nuc(e−1) ≈ 1.776×10³⁴ Pa.

    This pins the EOS to the PM-native linear compressible formula, not the
    radiation formula ρ_crit c²/3 (which gives 1.873×10³⁴ Pa, 5.5% too high).
    Tolerance 0.1% covers floating-point constants across platforms.
    """
    c = 299792458.0
    RHO_NUC = 2.3e17
    p_analytic = 0.5 * c * c * RHO_NUC * (math.e - 1.0)
    rel_err = abs(critical_state.p_crit - p_analytic) / p_analytic
    assert rel_err < 1e-3, (
        f"P_crit = {critical_state.p_crit:.6e} Pa differs from PM-EOS analytic "
        f"{p_analytic:.6e} Pa by {rel_err:.2e} (should be < 0.1%)"
    )


# ---------------------------------------------------------------------------
# Deformation-energy curvature vanishes at φ_crit
# ---------------------------------------------------------------------------

def test_U_pp_vanishes_at_phi_crit(critical_state: CriticalState) -> None:
    """U''(φ_crit) should be numerically zero.

    Tolerance: relative to the scale  U''(0) = 2 ε₀  with ε₀ = ρ_nuc c².
    After 64 bisection steps from an initial interval of width 0.02, the
    residual is at the level of floating-point rounding (~10⁻¹⁵ relative).
    """
    U_pp = pm_deformation_energy_deriv2(critical_state.phi_crit)
    U_pp_scale = abs(pm_deformation_energy_deriv2(0.0))   # = 2 ε₀ ≈ 4.1e34
    relative_residual = abs(U_pp) / U_pp_scale

    assert relative_residual < 1e-10, (
        f"|U''(φ_crit)| / U''(0) = {relative_residual:.2e}; "
        f"expected < 1e-10 (numerical root of U'')"
    )


# ---------------------------------------------------------------------------
# Sign-change: stability below, instability above φ_crit
# ---------------------------------------------------------------------------

def test_U_pp_sign_change_across_phi_crit(critical_state: CriticalState) -> None:
    """U'' must be positive just below φ_crit and negative just above."""
    eps = 1e-4
    U_pp_below = pm_deformation_energy_deriv2(critical_state.phi_crit - eps)
    U_pp_above = pm_deformation_energy_deriv2(critical_state.phi_crit + eps)
    assert U_pp_below > 0.0, (
        f"U''(φ_crit − ε) = {U_pp_below:.3e} should be positive (stable side)"
    )
    assert U_pp_above < 0.0, (
        f"U''(φ_crit + ε) = {U_pp_above:.3e} should be negative (unstable side)"
    )


# ---------------------------------------------------------------------------
# Self-consistency: CriticalState matches the PM EOS helpers
# ---------------------------------------------------------------------------

def test_rho_crit_consistent_with_pm_density_mapping(critical_state: CriticalState) -> None:
    """ρ_crit must match pm_density_from_phi(φ_crit) to machine precision."""
    rho_check = pm_density_from_phi(critical_state.phi_crit)
    rel_err = abs(rho_check - critical_state.rho_crit) / critical_state.rho_crit
    assert rel_err < 1e-12, (
        f"ρ_crit inconsistency: stored {critical_state.rho_crit:.6e}, "
        f"recomputed {rho_check:.6e}, relative error {rel_err:.2e}"
    )


def test_p_crit_consistent_with_pm_pressure_law(critical_state: CriticalState) -> None:
    """P_crit must match pm_pressure_from_phi(φ_crit) to machine precision."""
    p_check = pm_pressure_from_phi(critical_state.phi_crit)
    rel_err = abs(p_check - critical_state.p_crit) / critical_state.p_crit
    assert rel_err < 1e-12, (
        f"P_crit inconsistency: stored {critical_state.p_crit:.6e}, "
        f"recomputed {p_check:.6e}, relative error {rel_err:.2e}"
    )


# ---------------------------------------------------------------------------
# Regression: expected approximate value of φ_crit
# ---------------------------------------------------------------------------

def test_phi_crit_approximate_value(critical_state: CriticalState) -> None:
    """For the PM Landau cubic U''(φ) = 2ε₀(1 − φ), the exact root is φ_crit = 1."""
    assert abs(critical_state.phi_crit - 1.0) < 1e-10, (
        f"Expected φ_crit ≈ 1.0 (exact root of 2ε₀(1-φ)), "
        f"got {critical_state.phi_crit:.12f}"
    )


# ---------------------------------------------------------------------------
# Manual inspection helper
# ---------------------------------------------------------------------------

def test_print_critical_state(critical_state: CriticalState) -> None:
    """Print the critical state values for manual inspection (pytest -s)."""
    c = 299792458.0
    print()
    print("=" * 60)
    print("  PM Critical State: Matter / Energy Transition")
    print("=" * 60)
    print(f"  φ_crit   = {critical_state.phi_crit:.8f}    (dimensionless)")
    print(f"  n_crit   = exp(φ_crit) = {math.exp(critical_state.phi_crit):.8f}")
    print(f"  ρ_crit   = {critical_state.rho_crit:.6e} kg m⁻³")
    print(f"  P_crit   = {critical_state.p_crit:.6e} Pa")
    print(f"  ρ_nuc    = {RHO_NUC:.3e} kg m⁻³  (reference)")
    print(f"  ρ_crit / ρ_nuc = {critical_state.rho_crit / RHO_NUC:.6f}")
    print(f"  U''(φ_crit) = {pm_deformation_energy_deriv2(critical_state.phi_crit):.3e} J m⁻³")
    print("=" * 60)
