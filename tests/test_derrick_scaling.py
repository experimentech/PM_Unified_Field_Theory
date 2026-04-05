"""
Tests for the Derrick scaling argument applied to PM's vacuum field.

These tests codify the analytic results from scripts/derrick_soliton_analysis.py
and ensure they remain sound as the PM codebase evolves.

Physical context
----------------
PM's vacuum action (α=2 measure) gives static energy:

    E[n] = I_grad + I_pot
    I_grad = 4π ∫ ½(dn/dr)² r² dr  ≥ 0
    I_pot  = 4π ∫ V̂_vac(n) n² r² dr

Under Derrick rescaling n_λ(r) = n(λr):

    E(λ) = λ⁻¹ I_grad + λ⁻³ I_pot

Stationary point: dE/dλ|_{λ=1} = 0  ⟹  I_grad + 3 I_pot = 0
Stability:  d²E/dλ²|_{λ*} = −2 I_grad  <  0  (always a saddle)

Consequence: PM's classical field admits no stable static solitons.
ρ_nuc cannot be derived from the classical action alone.
"""

import math
import numpy as np
import pytest

from scripts.derrick_soliton_analysis import (
    v_vac,
    v_vac_prime,
    energy_integrals,
    derrick_energy,
    derrick_stationary,
    derrick_d2E,
    gaussian_soliton,
    lorentzian_soliton,
    analyse_soliton,
)


# ---------------------------------------------------------------------------
# 1. Vacuum potential properties
# ---------------------------------------------------------------------------

class TestVacuumPotential:
    """V̂_vac(n) fixed point and sign structure."""

    def test_fixed_point_value(self):
        """V̂_vac(1) = 0  (vacuum is the fixed point)."""
        assert abs(v_vac(1.0)) < 1e-12

    def test_fixed_point_derivative(self):
        """V̂_vac'(1) = 0  (n=1 is a critical point)."""
        assert abs(v_vac_prime(1.0)) < 1e-12

    def test_v_prime_nonnegative_everywhere(self):
        """V̂_vac'(n) = (n−1)²/n ≥ 0 for all n > 0."""
        n_vals = np.concatenate([
            np.linspace(0.01, 0.99, 100),
            np.linspace(1.01, 5.0,  100),
        ])
        assert np.all(v_vac_prime(n_vals) >= -1e-14)

    def test_v_prime_equals_analytic_form(self):
        """V̂_vac'(n) = (n−1)²/n  (AM-GM identity)."""
        n_vals = np.array([0.2, 0.5, 0.8, 1.2, 2.0, 5.0])
        analytic = (n_vals - 1.0)**2 / n_vals
        computed = v_vac_prime(n_vals)
        np.testing.assert_allclose(computed, analytic, rtol=1e-12)

    def test_v_vac_negative_for_n_less_than_1(self):
        """V̂_vac < 0 for n ∈ (0, 1) — rarefied medium is sub-vacuum."""
        n_vals = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 0.99])
        assert np.all(v_vac(n_vals) < 0)

    def test_v_vac_positive_for_n_greater_than_1(self):
        """V̂_vac > 0 for n > 1 — compressed medium is super-vacuum."""
        n_vals = np.array([1.01, 1.1, 1.5, 2.0, 5.0])
        assert np.all(v_vac(n_vals) > 0)

    def test_v_vac_formula(self):
        """V̂_vac(n) = ½n² + ln(n) − ½ − 2(n−1)."""
        n_vals = np.array([0.3, 0.7, 1.0, 1.5, 3.0])
        expected = 0.5*n_vals**2 + np.log(n_vals) - 0.5 - 2.0*(n_vals - 1.0)
        np.testing.assert_allclose(v_vac(n_vals), expected, rtol=1e-12)


# ---------------------------------------------------------------------------
# 2. Energy scaling law
# ---------------------------------------------------------------------------

class TestDerrickScalingLaw:
    """E(λ) = λ⁻¹ I_grad + λ⁻³ I_pot is exact for any trial function."""

    @pytest.fixture
    def rarefaction(self):
        """A rarefaction soliton that has I_pot < 0."""
        def n(r): return gaussian_soliton(r, A=-0.3, width=1.0)
        I_grad, I_pot, E = energy_integrals(n, r_max=20.0, n_pts=5000)
        return I_grad, I_pot

    def test_lambda_1_recovers_total_energy(self, rarefaction):
        """E(λ=1) = I_grad + I_pot."""
        I_grad, I_pot = rarefaction
        assert abs(derrick_energy(1.0, I_grad, I_pot) - (I_grad + I_pot)) < 1e-10

    def test_scaling_law_at_lambda_2(self, rarefaction):
        """E(2) = I_grad/2 + I_pot/8."""
        I_grad, I_pot = rarefaction
        expected = I_grad / 2.0 + I_pot / 8.0
        assert abs(derrick_energy(2.0, I_grad, I_pot) - expected) < 1e-10

    def test_scaling_law_at_lambda_half(self, rarefaction):
        """E(0.5) = 2·I_grad + 8·I_pot."""
        I_grad, I_pot = rarefaction
        expected = 2.0 * I_grad + 8.0 * I_pot
        assert abs(derrick_energy(0.5, I_grad, I_pot) - expected) < 1e-10

    def test_I_grad_always_nonneg(self):
        """Gradient energy I_grad ≥ 0 for any compression profile."""
        for A in [+0.1, +0.3, -0.3, -0.5]:
            def n(r, _A=A): return gaussian_soliton(r, A=_A, width=1.0)
            I_grad, _, _ = energy_integrals(n, r_max=15.0, n_pts=3000)
            assert I_grad >= 0.0, f"I_grad < 0 for A={A}"


# ---------------------------------------------------------------------------
# 3. Derrick first-order condition
# ---------------------------------------------------------------------------

class TestDerrickFirstOrderCondition:
    """Stationary point at λ* requires I_grad + 3·I_pot = 0."""

    def test_compression_no_stationary_point(self):
        """Compression soliton has I_pot > 0  →  no Derrick stationary point."""
        def n(r): return gaussian_soliton(r, A=+0.3, width=1.0)
        I_grad, I_pot, _ = energy_integrals(n)
        assert I_pot > 0, "I_pot should be positive for compression"
        assert derrick_stationary(I_grad, I_pot) is None

    def test_rarefaction_has_stationary_point(self):
        """Rarefaction soliton has I_pot < 0  →  stationary point exists."""
        def n(r): return gaussian_soliton(r, A=-0.5, width=1.0)
        I_grad, I_pot, _ = energy_integrals(n)
        assert I_pot < 0, "I_pot should be negative for rarefaction"
        lam_star = derrick_stationary(I_grad, I_pot)
        assert lam_star is not None
        assert lam_star > 0

    def test_stationary_point_formula(self):
        """λ* = sqrt(−3·I_pot / I_grad)."""
        I_grad, I_pot = 1.0, -0.2   # synthetic values
        lam_star = derrick_stationary(I_grad, I_pot)
        expected = math.sqrt(3.0 * 0.2 / 1.0)
        assert abs(lam_star - expected) < 1e-12

    def test_derrick_condition_at_lambda_star(self):
        """At λ*, dE/dλ = 0:  −I_grad/λ² − 3·I_pot/λ⁴ = 0."""
        def n(r): return lorentzian_soliton(r, A=-0.4, width=1.5)
        I_grad, I_pot, _ = energy_integrals(n, r_max=20.0, n_pts=5000)
        lam_star = derrick_stationary(I_grad, I_pot)
        if lam_star is None:
            pytest.skip("No stationary point for this profile")
        # Numerical derivative of E(λ) at λ*
        eps = 1e-5
        dE = (derrick_energy(lam_star + eps, I_grad, I_pot)
              - derrick_energy(lam_star - eps, I_grad, I_pot)) / (2 * eps)
        assert abs(dE) < 1e-4, f"dE/dλ at λ* = {dE:.2e}, expected ≈ 0"


# ---------------------------------------------------------------------------
# 4. Derrick second-order instability — the central result
# ---------------------------------------------------------------------------

class TestDerrickInstability:
    """
    The key theorem: d²E/dλ² at any Derrick stationary point equals
    −2·I_grad < 0.  All solitons are scale-unstable.
    """

    def _get_rarefaction(self, A: float = -0.4, width: float = 1.5):
        def n(r): return gaussian_soliton(r, A=A, width=width)
        I_grad, I_pot, _ = energy_integrals(n, r_max=20.0, n_pts=5000)
        lam_star = derrick_stationary(I_grad, I_pot)
        return I_grad, I_pot, lam_star

    def test_second_derivative_formula(self):
        """d²E/dλ²|_{λ=1} = 2·I_grad + 12·I_pot (general formula)."""
        I_grad, I_pot = 1.0, -0.1
        eps = 1e-4
        lam = 1.0
        numerical = (
            derrick_energy(lam + eps, I_grad, I_pot)
            - 2 * derrick_energy(lam, I_grad, I_pot)
            + derrick_energy(lam - eps, I_grad, I_pot)
        ) / eps**2
        analytic = 2.0 * I_grad + 12.0 * I_pot
        assert abs(numerical - analytic) < 1e-3

    def test_d2E_at_stationary_point_equals_minus2_I_grad(self):
        """At λ*:  d²E/dλ² = −2·I_grad/λ*³  (derived analytically).

        From the stationary condition:  I_pot = −I_grad·λ*²/3
        Substituting into  d²E/dλ²|_{λ*} = 2·I_grad/λ*³ + 12·I_pot/λ*⁵:
            = 2·I_grad/λ*³ + 12·(−I_grad·λ*²/3)/λ*⁵
            = 2·I_grad/λ*³ − 4·I_grad/λ*³
            = −2·I_grad/λ*³
        """
        I_grad, I_pot, lam_star = self._get_rarefaction()
        if lam_star is None:
            pytest.skip("No stationary point")
        d2E = derrick_d2E(lam_star, I_grad, I_pot)
        # Stationary condition: I_pot = −I_grad·λ*²/3
        I_pot_at_stat = -I_grad * lam_star**2 / 3.0
        # Analytic result: 2/λ³·I_grad + 12/λ⁵·I_pot
        expected_exact = 2.0 * I_grad / lam_star**3 + 12.0 * I_pot_at_stat / lam_star**5
        assert abs(d2E - expected_exact) < 1e-6

    def test_d2E_negative_at_stationary_point(self):
        """d²E/dλ² < 0 at the Derrick stationary point — always a saddle."""
        I_grad, I_pot, lam_star = self._get_rarefaction()
        if lam_star is None:
            pytest.skip("No stationary point")
        d2E = derrick_d2E(lam_star, I_grad, I_pot)
        assert d2E < 0, (
            f"d²E/dλ² = {d2E:.6f} at λ* = {lam_star:.4f};  "
            f"expected < 0 (scale instability)"
        )

    def test_instability_multiple_profiles(self):
        """Scale instability holds for all tested rarefaction profiles."""
        configs = [
            (-0.2, 1.0), (-0.4, 1.0), (-0.4, 2.0), (-0.6, 0.5),
        ]
        for A, width in configs:
            I_grad, I_pot, lam_star = self._get_rarefaction(A=A, width=width)
            if lam_star is None:
                continue  # no stationary point → trivially no stable soliton
            d2E = derrick_d2E(lam_star, I_grad, I_pot)
            assert d2E < 0, (
                f"Stable soliton found for A={A}, width={width}!  "
                f"d²E/dλ² = {d2E:.6f} (Derrick theorem violated)"
            )

    def test_no_stable_classical_soliton_in_3d(self):
        """
        Summary test: no combination of amplitude and width yields a stable
        static soliton in PM's 3D classical field.
        """
        # Try a wide sweep of parameters
        for A in [-0.1, -0.3, -0.5, -0.7]:
            for width in [0.5, 1.0, 2.0, 4.0]:
                def n(r, _A=A, _w=width):
                    return gaussian_soliton(r, A=_A, width=_w)
                I_grad, I_pot, _ = energy_integrals(n, r_max=20.0, n_pts=2000)
                lam_star = derrick_stationary(I_grad, I_pot)
                if lam_star is not None:
                    d2E = derrick_d2E(lam_star, I_grad, I_pot)
                    assert d2E < 0, (
                        f"Stable soliton at A={A}, width={width}: d²E={d2E:.6f}"
                    )


# ---------------------------------------------------------------------------
# 5. Physical consequence: ρ_nuc
# ---------------------------------------------------------------------------

class TestRhoNucConsequence:
    """
    The Derrick result implies ρ_nuc cannot come from the classical action.
    These tests verify the logical chain.
    """

    def test_compression_I_pot_positive(self):
        """For n > 1 (compressed), V̂_vac > 0, so I_pot > 0."""
        def n(r): return gaussian_soliton(r, A=+0.3, width=1.0)
        _, I_pot, _ = energy_integrals(n)
        assert I_pot > 0

    def test_rarefaction_I_pot_negative(self):
        """For n < 1 (rarefied), V̂_vac < 0, so I_pot < 0."""
        def n(r): return gaussian_soliton(r, A=-0.3, width=1.0)
        _, I_pot, _ = energy_integrals(n)
        assert I_pot < 0

    def test_derrick_condition_forces_rarefaction(self):
        """
        The Derrick stationary-point condition I_grad + 3·I_pot = 0
        requires I_pot = −I_grad/3 < 0, which requires the field to go
        sub-vacuum (n < 1) somewhere.
        A pure compression field cannot satisfy the condition.
        """
        I_grad = 2.5   # positive by construction
        # For a compression: I_pot > 0
        I_pot_compression = 0.8
        lam_comp = derrick_stationary(I_grad, I_pot_compression)
        assert lam_comp is None, "No stationary point for pure compression"

        # For tuned rarefaction satisfying Derrick exactly:
        I_pot_rarefaction = -I_grad / 3.0
        lam_rare = derrick_stationary(I_grad, I_pot_rarefaction)
        assert lam_rare is not None
        assert abs(lam_rare - 1.0) < 1e-10, "λ* = 1 when already at stationary point"

    def test_analytical_instability_proof(self):
        """
        Prove analytically: at any Derrick stationary point,
        d²E/dλ²|_{λ*} = −2·I_grad/λ*³ < 0.

        From the Derrick condition: I_pot = −λ*² · I_grad / 3
        Substituting into d²E/dλ² = 2/λ*³ · I_grad + 12/λ*⁵ · I_pot:
            = 2/λ*³ · I_grad + 12/λ*⁵ · (−λ*² I_grad / 3)
            = 2/λ*³ · I_grad − 4/λ*³ · I_grad
            = −2/λ*³ · I_grad  <  0
        """
        for I_grad in [0.1, 1.0, 5.0]:
            for lam_star in [0.5, 1.0, 2.0]:
                I_pot_star = -I_grad * lam_star**2 / 3.0
                d2E = derrick_d2E(lam_star, I_grad, I_pot_star)
                expected = -2.0 * I_grad / lam_star**3
                assert abs(d2E - expected) < 1e-10, (
                    f"I_grad={I_grad}, λ*={lam_star}: "
                    f"d²E={d2E:.6f}, expected={expected:.6f}"
                )
                assert d2E < 0, "Must be negative — scale instability"
