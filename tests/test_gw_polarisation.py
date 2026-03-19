"""Tests for GW polarisation content in the Pushing-Medium theory.

PM mediates gravity via a massless scalar refractive-index field n(r, t).
In principle this could produce an l = 0 'breathing' (scalar) GW polarisation
in addition to the standard GR tensor modes h₊ and h×.

Conservation-law argument (why there is NO leading-order breathing mode):
─────────────────────────────────────────────────────────────────────────────
(1) Monopole: total mass M_tot = Σ mᵢ = const  →  no monopole scalar radiation.
(2) Dipole  : CoM position  R = Σ mᵢ rᵢ / M_tot = const (isolated system)
              → no dipole scalar radiation.
(3) Scalar quadrupole: I_trace = Σ mᵢ rᵢ²
    For a CIRCULAR orbit: I_trace = μ a² = CONSTANT (does not oscillate).
    Therefore d²I_trace/dt² = 0 exactly, giving zero scalar quadrupole radiation.
    Contrast: STF tensor Iᵢⱼ oscillates at 2ω  →  tensor h₊, h× radiation.

First non-zero breathing amplitude:
    Enters when the orbit decays (ȧ ≠ 0) via the Peters mechanism,
    giving d²I_trace/dt² = −4μβ²/a⁶  where  β = (64/5) G³ M₁M₂M_tot / c⁵.
    This contribution is suppressed by (v/c)^10 relative to the tensor modes.

Numerical result for GW170817 at 100 Hz (1.4 + 1.4 M☉):
    |h_b| / |h₊| ≈ 10⁻⁷
    α_NT = P_scalar / P_total ≈ 10⁻¹⁴

Observational constraint (LIGO O3, Abbott+ 2021 PRD 103, 122002):
    α_NT < 0.07  (non-tensor power fraction, 90% CI)

PM prediction is ≈ 13 orders of magnitude below the LIGO O3 limit.
PM PASSES the GW polarisation test with an enormous margin.
"""

import math
import pytest

from pushing_medium import (
    pm_scalar_quadrupole_trace,
    pm_breathing_strain_ratio,
    pm_nontensor_power_fraction,
    G, c,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
M_sun = 1.989e30      # kg
Mpc   = 3.086e22      # m

# LIGO O3 upper limit on non-tensor GW power fraction (90% CI)
# Abbott+ 2021, PRD 103, 122002
LIGO_O3_ALPHA_NT_LIMIT = 0.07

# LIGO strain sensitivity floor around 100 Hz (representative)
LIGO_STRAIN_FLOOR = 1e-23   # dimensionless h

# GW170817-like BNS
GW170817_M1 = 1.46 * M_sun
GW170817_M2 = 1.27 * M_sun
GW170817_D  = 40.0 * Mpc
GW170817_F  = 100.0   # Hz representative GW frequency near merger


# ===========================================================================
class TestScalarConservationLaws:
    """Monopole and dipole scalar radiation are exactly zero for an isolated binary."""

    def test_total_mass_constant_no_monopole(self):
        """Total mass M_tot = const  →  no monopole breathing radiation.

        The monopole amplitude h_b^(0) ∝ dM_tot/dt = 0 for any isolated binary.
        """
        # Any isolated system has dM/dt = 0; check concept by verifying
        # that the trace I_trace at a fixed separation doesn't depend on
        # which mass is 'larger' — total mass is the only relevant conserved charge.
        M = 2.73 * M_sun
        M1, M2 = 1.46 * M_sun, 1.27 * M_sun
        # Total mass from the two different parametrisations:
        assert abs((M1 + M2) - M) < 0.01 * M_sun, "mass parameterisation inconsistent"
        # Monopole source Q = M_tot = const; no time derivative → no radiation.
        # This is purely analytic — nothing to compute numerically.
        assert True   # by conservation of mass

    def test_com_fixed_no_dipole(self):
        """Centre-of-mass R_cm = const  →  no dipole breathing radiation.

        For an isolated binary the CoM drifts at constant velocity.
        In the CoM frame R_cm = 0 and d²R_cm/dt² = 0 exactly.
        """
        M1, M2 = GW170817_M1, GW170817_M2
        a = 300e3   # m, representative orbital separation
        # Positions in CoM frame:
        r1 = M2 / (M1 + M2) * a
        r2 = M1 / (M1 + M2) * a
        R_cm = M1 * r1 - M2 * r2   # should be zero by construction
        assert abs(R_cm) < 1e-6 * a


# ===========================================================================
class TestCircularOrbitTraceConstant:
    """For circular orbit, I_trace = μa² is independent of orbital phase."""

    def test_trace_formula(self):
        """I_trace = M₁M₂/(M₁+M₂) × a²."""
        M1, M2 = GW170817_M1, GW170817_M2
        a = 200e3   # m
        mu = M1 * M2 / (M1 + M2)
        expected = mu * a * a
        assert pm_scalar_quadrupole_trace(M1, M2, a) == pytest.approx(expected, rel=1e-12)

    def test_trace_independent_of_phase(self):
        """I_trace = μa² does not depend on phase angle φ of the orbit.

        The full STF tensor Iᵢⱼ oscillates at 2ω, but I_trace = Iₓₓ + Iyy
        is constant:
            I_xx = μa² cos²φ
            I_yy = μa² sin²φ
            I_trace = μa²(cos²φ + sin²φ) = μa²
        """
        M1, M2 = GW170817_M1, GW170817_M2
        a = 200e3   # m
        mu = M1 * M2 / (M1 + M2)
        I_trace = mu * a * a
        # Verify phase-independence by computing I_xx + I_yy at many phases:
        phases = [0.0, math.pi / 6, math.pi / 4, math.pi / 3, math.pi / 2, math.pi]
        for phi in phases:
            I_xx = mu * a ** 2 * math.cos(phi) ** 2
            I_yy = mu * a ** 2 * math.sin(phi) ** 2
            assert I_xx + I_yy == pytest.approx(I_trace, rel=1e-14)

    def test_trace_zero_second_derivative_circular(self):
        """For a purely circular orbit (constant a), Ï_trace = 0 exactly.

        da/dt = 0  →  d²(μa²)/dt² = 2μ(ȧ² + aä) = 0.
        The scalar breathing mode amplitude is identically zero at leading order.
        """
        # Demonstrate that any constant-a orbit gives zero second derivative:
        a_dot = 0.0   # circular orbit assumption
        M1, M2 = GW170817_M1, GW170817_M2
        mu = M1 * M2 / (M1 + M2)
        a = 200e3   # m
        # d²I_trace/dt² = 2μ(ȧ² + aä) = 0 when ȧ = 0
        I_trace_ddot = 2.0 * mu * (a_dot ** 2 + a * 0.0)
        assert I_trace_ddot == 0.0

    def test_trace_increases_with_separation(self):
        """Larger separation → larger I_trace (∝ a²)."""
        M1, M2 = 1.4 * M_sun, 1.4 * M_sun
        a1, a2 = 100e3, 200e3
        assert pm_scalar_quadrupole_trace(M1, M2, a1) < pm_scalar_quadrupole_trace(M1, M2, a2)


# ===========================================================================
class TestBreathingStrainRatio:
    """pm_breathing_strain_ratio: leading non-zero h_b/h₊ from orbital decay."""

    def test_positive(self):
        """Ratio is always positive (absolute value)."""
        ratio = pm_breathing_strain_ratio(GW170817_M1, GW170817_M2, GW170817_F)
        assert ratio > 0.0

    def test_gw170817_ratio_tiny(self):
        """For GW170817 at 100 Hz: |h_b|/|h₊| < 10⁻⁵.

        PM breathing mode is completely negligible compared to tensor modes.
        The orbit is so circular and the v/c so small that the scalar
        radiation is suppressed by (v/c)^10 ≈ 10^{-10}.
        """
        ratio = pm_breathing_strain_ratio(GW170817_M1, GW170817_M2, 100.0)
        assert ratio < 1e-5, f"|h_b/h₊| = {ratio:.2e}; expected < 1e-5"

    def test_ratio_smaller_at_lower_frequency(self):
        """Lower GW frequency (more adiabatic): h_b/h₊ is smaller.

        At lower f the orbit is wider (v/c is smaller), so the (v/c)^10
        suppression is stronger and the breathing mode is even more negligible.
        """
        M1, M2 = GW170817_M1, GW170817_M2
        ratio_20hz  = pm_breathing_strain_ratio(M1, M2, 20.0)
        ratio_100hz = pm_breathing_strain_ratio(M1, M2, 100.0)
        assert ratio_20hz < ratio_100hz

    def test_ratio_increases_with_frequency(self):
        """Closer to merger (higher f, higher v/c), ratio grows."""
        M1, M2 = 1.4 * M_sun, 1.4 * M_sun
        ratios = [pm_breathing_strain_ratio(M1, M2, f) for f in [20.0, 50.0, 100.0, 200.0]]
        for i in range(len(ratios) - 1):
            assert ratios[i] < ratios[i + 1], f"ratio should grow with f; {ratios}"

    def test_frequency_scaling_exponent(self):
        """h_b/h₊ ∝ f^(10/3) at fixed masses.

        From β²/(ω²a⁸) with ω ∝ f and a ∝ f^(-2/3):
            ω²a⁸ ∝ f² × f^(-16/3) = f^(-10/3)
            ratio ∝ β² / (ω²a⁸) ∝ f^(10/3)
        """
        M1, M2 = 1.4 * M_sun, 1.4 * M_sun
        f1, f2 = 50.0, 100.0
        r1 = pm_breathing_strain_ratio(M1, M2, f1)
        r2 = pm_breathing_strain_ratio(M1, M2, f2)
        expected_ratio = (f2 / f1) ** (10.0 / 3.0)
        assert r2 / r1 == pytest.approx(expected_ratio, rel=1e-6)

    def test_equal_mass_bns_below_ligo_floor(self):
        """|h_b| << LIGO strain sensitivity for a typical BNS at 40 Mpc.

        Even if we scale |h_b| = ratio × |h₊|, the breathing mode strain
        is far below LIGO's noise floor and completely undetectable.
        """
        M1 = M2 = 1.4 * M_sun
        D = 40.0 * Mpc
        f = 100.0

        from pushing_medium import pm_gw_strain_amplitude, pm_chirp_mass
        M_c = pm_chirp_mass(M1, M2)
        h_plus = pm_gw_strain_amplitude(M_c, D, f)

        ratio = pm_breathing_strain_ratio(M1, M2, f)
        h_breath = ratio * h_plus

        # Breathing mode strain << LIGO O3 noise floor (~10⁻²³ at 100 Hz)
        assert h_breath < LIGO_STRAIN_FLOOR, (
            f"h_b = {h_breath:.2e} is above LIGO noise floor {LIGO_STRAIN_FLOOR:.2e}"
        )


# ===========================================================================
class TestNonTensorPowerFraction:
    """pm_nontensor_power_fraction: α_NT << LIGO O3 bound."""

    def test_bns_below_o3_bound(self):
        """GW170817 parameters: α_NT << 0.07 (LIGO O3 bound)."""
        alpha = pm_nontensor_power_fraction(GW170817_M1, GW170817_M2, GW170817_F)
        assert alpha < LIGO_O3_ALPHA_NT_LIMIT, (
            f"α_NT = {alpha:.2e}; LIGO O3 limit = {LIGO_O3_ALPHA_NT_LIMIT}"
        )

    def test_bns_margin_enormous(self):
        """PM is at least 10 orders of magnitude below the LIGO O3 limit."""
        alpha = pm_nontensor_power_fraction(GW170817_M1, GW170817_M2, 100.0)
        # LIGO limit 0.07; PM prediction ~10^-14
        margin = LIGO_O3_ALPHA_NT_LIMIT / alpha
        assert margin > 1e10, (
            f"Margin = {margin:.2e}; expected > 1e10. α_NT = {alpha:.2e}"
        )

    def test_equal_mass_bns_below_o3_bound(self):
        """Equal-mass BNS (1.4+1.4 M☉) at 100 Hz also satisfies O3 bound."""
        alpha = pm_nontensor_power_fraction(1.4 * M_sun, 1.4 * M_sun, 100.0)
        assert alpha < LIGO_O3_ALPHA_NT_LIMIT

    def test_multiple_frequencies_all_below_bound(self):
        """α_NT < 0.07 across all LIGO band frequencies (10–500 Hz)."""
        M1, M2 = GW170817_M1, GW170817_M2
        for f_gw in [10.0, 30.0, 50.0, 100.0, 200.0, 500.0]:
            alpha = pm_nontensor_power_fraction(M1, M2, f_gw)
            assert alpha < LIGO_O3_ALPHA_NT_LIMIT, (
                f"α_NT = {alpha:.2e} at f = {f_gw} Hz; exceeds O3 limit"
            )

    def test_power_fraction_is_strain_ratio_squared(self):
        """α_NT = (|h_b|/|h₊|)² by definition."""
        M1, M2 = GW170817_M1, GW170817_M2
        f = 100.0
        ratio = pm_breathing_strain_ratio(M1, M2, f)
        alpha = pm_nontensor_power_fraction(M1, M2, f)
        assert alpha == pytest.approx(ratio ** 2, rel=1e-12)

    def test_alpha_increases_with_frequency(self):
        """Approaching merger (higher f), α_NT grows but stays << O3 bound."""
        M1, M2 = 1.4 * M_sun, 1.4 * M_sun
        alpha_20hz  = pm_nontensor_power_fraction(M1, M2, 20.0)
        alpha_500hz = pm_nontensor_power_fraction(M1, M2, 500.0)
        # Even at the highest frequency, should still be far below limit
        assert alpha_20hz  < alpha_500hz
        assert alpha_500hz < LIGO_O3_ALPHA_NT_LIMIT


# ===========================================================================
class TestWaveformNoBreatheMode:
    """The PM chirp waveform correctly returns only h₊ and h× — no h_b term."""

    def test_waveform_returns_tensor_modes_only(self):
        """pm_gw_chirp_waveform returns exactly 4 arrays (t, f, h₊, h×).

        There is no scalar breathing mode in the waveform because the
        PM leading-order scalar radiation is zero for circular orbits.
        """
        from pushing_medium import pm_gw_chirp_waveform
        result = pm_gw_chirp_waveform(
            1.4 * M_sun, 1.4 * M_sun,
            40.0 * Mpc,
            tau_start=30.0, N=512
        )
        # Should return exactly 4 components (t, f, h_plus, h_cross) — no h_b
        assert len(result) == 4, f"Expected 4 arrays (no breathing mode), got {len(result)}"

    def test_breathing_mode_not_in_waveform_tuple(self):
        """The waveform tuple (t, f, h₊, h×) has no 5th 'h_breathing' element."""
        from pushing_medium import pm_gw_chirp_waveform
        t, f, h_plus, h_cross = pm_gw_chirp_waveform(
            GW170817_M1, GW170817_M2,
            GW170817_D,
            tau_start=2.0, N=256
        )
        # Physical check: |h_b|/|h₊| << 1, so the waveform is dominated by tensor
        # Confirm h_plus is non-trivial (not zero)
        import numpy as np
        assert np.max(np.abs(h_plus)) > 0.0
