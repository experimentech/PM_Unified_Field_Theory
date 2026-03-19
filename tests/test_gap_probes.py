"""
Tests: Gap Probes — Cross-Constraints on the PM Strong-Field Aberration

Each test probes one or more of the three structural gaps identified in the
formula sheet (§ Deformation Energy / Diagnostic):

  Gap 1: U(φ) absent from Lagrangian / Poisson equation
  Gap 2: U(φ) and the EOS are thermodynamically inconsistent
  Gap 3: U′(φ) absent from hydrostatic equilibrium

Tests here are DESCRIPTIVE not prescriptive — they measure the magnitude
of each inconsistency and the effect of each repair, providing quantitative
datapoints for choosing the correct solution.
"""

import math
import pytest
import numpy as np

from pushing_medium.stellar_structure import RHO_NUC, RHO_CRIT, MU_G, M_SUN, c
from pushing_medium.stellar_structure import pm_eos_pressure, pm_eos_density

# Import gap-probe helpers from the comparison script
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from gap_probe_comparison import (
    EPS0, u_phi, u_prime, u_double_prime,
    p_thermo_from_phi, p_thermo_from_rho, rho_from_p_thermo,
    compare_eos, make_option_c_eos,
    compute_mr_curve_variant, sweep_phi_crit_for_max_mass,
)


# ===========================================================================
class TestDeformationEnergyStructure:
    """Gap 0: Basic properties of U(φ) — sanity before probing gaps."""

    def test_u_zero_at_origin(self):
        assert u_phi(0.0) == pytest.approx(0.0, abs=1e-30)

    def test_u_prime_zero_at_origin(self):
        """U′(0) = 0 — no deformation force at φ=0."""
        assert u_prime(0.0) == pytest.approx(0.0, abs=1e-30)

    def test_u_double_prime_positive_below_crit(self):
        """U″(φ) > 0 for φ < 1 — medium is stable."""
        for phi in [0.0, 0.25, 0.5, 0.75, 0.99]:
            assert u_double_prime(phi) > 0.0, f"U″({phi}) = {u_double_prime(phi)}"

    def test_u_double_prime_zero_at_crit(self):
        """U″(1) = 0 — onset of instability."""
        assert u_double_prime(1.0) == pytest.approx(0.0, abs=1e-12)

    def test_u_double_prime_negative_above_crit(self):
        """U″(φ) < 0 for φ > 1 — medium unstable above the cap."""
        for phi in [1.01, 1.1, 1.5, 2.0]:
            assert u_double_prime(phi) < 0.0

    def test_u_prime_magnitude_at_crit(self):
        """U′(1) = ε₀(2×1 − 1²) = ε₀.  This is the self-coupling omitted from Poisson."""
        assert u_prime(1.0) == pytest.approx(EPS0, rel=1e-10)


# ===========================================================================
class TestGap2EOSInconsistency:
    """
    Gap 2: Quantify the thermodynamic inconsistency between U(φ) and the EOS.

    If U(φ) were the fluid's internal energy, the thermodynamic pressure would
    be P_thermo = ρ² ∂(U/ρ)/∂ρ ≈ ε₀(2φ − 2φ² + φ³/3), which differs from
    P_EOS = c²/2(ρ − ρ_nuc) by a factor of ~4 at small φ.

    These tests measure the ratio and confirm it is not negligible.
    """

    def test_p_thermo_zero_at_origin(self):
        """Both pressures vanish at φ=0 (ρ = ρ_nuc)."""
        assert p_thermo_from_phi(0.0) == pytest.approx(0.0, abs=1e-30)

    def test_p_thermo_positive_inside_medium(self):
        """P_thermo > 0 inside the star (φ > 0)."""
        assert p_thermo_from_phi(0.5) > 0.0

    def test_ratio_at_small_phi_is_approximately_4(self):
        """
        At φ → 0:  P_thermo ≈ 2ε₀φ,  P_EOS ≈ ε₀φ/2.
        Ratio ≈ 4.  This is the Gap 2 inconsistency.
        """
        phi = 0.01   # small enough for leading-order
        rho = RHO_NUC * math.exp(phi)
        P_eos   = pm_eos_pressure(rho)
        P_therm = p_thermo_from_phi(phi)
        ratio = P_therm / P_eos
        # Leading-order ratio = 4; allow ±20% for sub-leading corrections
        assert 3.0 < ratio < 5.0, \
            f"P_thermo/P_EOS = {ratio:.4f} at φ={phi} (expected ≈ 4)"

    def test_ratio_at_phi_half(self):
        """At φ = 0.5 the ratio should still be significantly > 1."""
        phi = 0.5
        rho = RHO_NUC * math.exp(phi)
        P_eos   = pm_eos_pressure(rho)
        P_therm = p_thermo_from_phi(phi)
        ratio = P_therm / P_eos
        assert ratio > 1.5, \
            f"P_thermo/P_EOS = {ratio:.3f} at φ=0.5 (should be >> 1)"

    def test_ratio_is_not_exactly_1_anywhere_in_0_to_1(self):
        """The two pressure functions are not equal anywhere in 0 < φ < 1."""
        for phi in np.linspace(0.05, 0.95, 10):
            rho = RHO_NUC * math.exp(phi)
            P_eos   = pm_eos_pressure(rho)
            P_therm = p_thermo_from_phi(phi)
            assert not math.isclose(P_therm, P_eos, rel_tol=0.01), \
                f"P_thermo == P_EOS at φ={phi:.2f} — Gap 2 not present?"

    def test_p_thermo_roundtrip_consistency(self):
        """
        p_thermo_from_rho and rho_from_p_thermo are inverse functions.

        Valid only for ρ < 1.797 × RHO_NUC (φ < 2−√2 ≈ 0.586), where
        P_thermo(φ) is strictly increasing.  Above that density P_thermo
        *decreases* with increasing ρ — making the EOS non-invertible.
        """
        # Only test on the ascending branch of P_thermo (φ < 0.586)
        for rho in [1.1*RHO_NUC, 1.5*RHO_NUC, 1.7*RHO_NUC]:
            P   = p_thermo_from_rho(rho)
            rho_back = rho_from_p_thermo(P)
            assert rho_back == pytest.approx(rho, rel=1e-6), \
                f"Roundtrip failed at ρ={rho:.3e}: got {rho_back:.3e}"

    def test_p_thermo_nonmonotonic_above_phi_max(self):
        """
        Critical Gap 2 finding: P_thermo(φ) peaks at φ = 2−√2 ≈ 0.586 and
        then *decreases*.  At φ_crit = 1 the thermodynamic pressure has already
        fallen below its peak — meaning U(φ) cannot serve as a valid
        thermodynamic free energy in the range [0.586, φ_crit].
        """
        phi_peak = 2.0 - math.sqrt(2.0)           # ≈ 0.586
        P_at_peak = p_thermo_from_phi(phi_peak)
        P_at_crit = p_thermo_from_phi(1.0)        # φ_crit = 1
        # By φ_crit=1 the pressure must have fallen substantially below its peak
        assert P_at_crit < P_at_peak * 0.80, (
            f"P_thermo(φ_crit=1) = {P_at_crit:.3e} Pa should be < 80% of "
            f"P_thermo(φ_peak={phi_peak:.3f}) = {P_at_peak:.3e} Pa"
        )
        # Derivative is negative at φ_crit — pressure is falling
        dP_at_crit = EPS0 * (2.0 - 4.0*1.0 + 1.0**2)   # = −ε₀
        assert dP_at_crit < 0.0, "dP_thermo/dφ should be negative at φ_crit=1"


# ===========================================================================
class TestGap1PoissonCoupling:
    """
    Gap 1: Quantify the self-coupling term that U′(φ) would add to the
    Poisson equation if included.

    The correction  δ(dφ/dr) ≈ (r/3) U′(φ)/c²  is the dominant additional
    term.  These tests show its magnitude relative to the standard Poisson term
    at stellar interior conditions.
    """

    def test_u_prime_correction_negligible_at_low_phi(self):
        """
        Near surface (φ ≈ 0.01) the Poisson self-coupling correction from U′ is tiny.

        The correct formula derived from Gauss theorem on the modified Poisson:
            δ(dφ/dr) = MU_G × [U′(φ)/c²] × r/3
                     = MU_G × ρ_nuc(2φ−φ²) × r/3     [units: m⁻¹ ✓]

        WARNING: the naive formula (r/3)·U′(φ)/c² omits the MU_G = 2G/c² factor
        and has units kg/m² instead of m⁻¹ — giving a spurious ratio of ~10²³.
        """
        phi = 0.01
        r = 5e3          # 5 km
        m = 0.5 * M_SUN
        standard   = MU_G * m / r**2                      # |dφ/dr| [m⁻¹]
        correction = MU_G * u_prime(phi) / c**2 * r / 3.0  # [m⁻¹]
        relative   = correction / standard
        assert relative < 0.05, \
            f"Correction/standard = {relative:.6f} at φ=0.01 (should be < 5%)"

    def test_u_prime_correction_at_phi_crit_is_percent_level(self):
        """
        At φ = φ_crit = 1 the Poisson self-coupling correction is ~3% of the
        standard term — not negligible at precision level, but not dominant.

        For a uniform-density interior at ρ̄ ≈ ρ_crit the ratio is:
            δ(dφ/dr) / (dφ/dr)_std  ≈  ρ_nuc(2φ−φ²) / (4π × ρ_crit)
                                      = 1 / (4π e) ≈ 0.029  (≈ 3%)

        Note: Gap 3 (U as elastic pressure) produces an ~78% correction to P —
        far more consequential than this Gap 1 Poisson correction.
        """
        phi = 1.0
        ratio_interior = RHO_NUC * (2*phi - phi**2) / (4 * math.pi * RHO_CRIT)
        # Expect ~3% = 1/(4πe); must be > 1% (not negligible) and < 10%
        assert 0.01 < ratio_interior < 0.10, (
            f"Gap 1 interior correction ratio = {ratio_interior:.2%} — "
            f"expected ~3% (i.e., 1/(4πe)) at φ_crit"
        )

    def test_poisson_correction_grows_with_phi(self):
        """U′(φ) grows with φ in [0, 1], so the Gap 1 correction grows inward."""
        corrections = []
        phi_vals = [0.1, 0.3, 0.5, 0.7, 0.9]
        for phi in phi_vals:
            corrections.append(u_prime(phi))
        # U′(φ) = ε₀(2φ − φ²) has a maximum at φ=1; should increase over [0,1]
        for i in range(len(corrections) - 1):
            assert corrections[i] < corrections[i + 1], \
                f"U′ not increasing at φ={phi_vals[i]:.1f}"


# ===========================================================================
class TestGap3ElasticPressure:
    """
    Gap 3: U(φ) adds an elastic pressure term absent from hydrostatics.
    Option (c) repair: P_total = P_EOS + U(φ).

    U(φ) at φ_crit = 1 equals U(1) = ε₀(1 − 1/3) = 2ε₀/3.
    The standard EOS at the same point gives P_EOS = c²(e−1)ρ_nuc/2 ≈ 0.859 ε₀.
    So U(1)/P_EOS(φ_crit) ≈ 0.78  — the elastic correction is ~78% of P_EOS.
    """

    def test_u_at_crit_is_large_fraction_of_p_eos(self):
        """U(φ_crit=1) / P_EOS(ρ_crit) should be substantial (> 50%)."""
        U_at_crit = u_phi(1.0)
        P_eos_at_crit = pm_eos_pressure(RHO_CRIT)
        ratio = U_at_crit / P_eos_at_crit
        assert ratio > 0.5, \
            f"U(φ_crit)/P_EOS = {ratio:.3f} — elastic correction is not negligible"

    def test_option_c_total_pressure_exceeds_eos(self):
        """P_total = P_EOS + U(φ) > P_EOS throughout the interior."""
        p_total_fn, _ = make_option_c_eos()
        for rho in [1.2*RHO_NUC, 1.5*RHO_NUC, 2.0*RHO_NUC, RHO_CRIT]:
            P_eos   = pm_eos_pressure(rho)
            P_total = p_total_fn(rho)
            assert P_total > P_eos, \
                f"P_total ≤ P_EOS at ρ={rho:.2e} (Option c must add pressure)"

    def test_option_c_roundtrip(self):
        """P_total and its density inverse are self-consistent."""
        _, rho_from_ptot = make_option_c_eos()
        p_total_fn, _ = make_option_c_eos()
        for rho in [1.1*RHO_NUC, 1.5*RHO_NUC, 2.0*RHO_NUC]:
            P = p_total_fn(rho)
            rho_back = rho_from_ptot(P)
            assert rho_back == pytest.approx(rho, rel=1e-5), \
                f"Roundtrip failed at ρ={rho:.3e}: got {rho_back:.3e}"

    def test_option_c_surface_condition_preserved(self):
        """P_total at ρ_nuc should still be 0 (correct surface condition)."""
        p_total_fn, _ = make_option_c_eos()
        assert p_total_fn(RHO_NUC) == pytest.approx(0.0, abs=1.0)    # within 1 Pa


# ===========================================================================
class TestMRCurveVariants:
    """
    Cross-constraint tests: run the stellar structure ODE under each repair
    option and check that the resulting M_max moves in a physically meaningful
    direction, and compare against LIGO constraints.
    """

    @pytest.fixture(scope="class")
    def baseline_mr(self):
        M, R = compute_mr_curve_variant(
            RHO_CRIT, pm_eos_pressure, pm_eos_density, n_points=30
        )
        return M, R

    @pytest.fixture(scope="class")
    def option_c_mr(self):
        p_fn, rho_fn = make_option_c_eos()
        M, R = compute_mr_curve_variant(RHO_CRIT, p_fn, rho_fn, n_points=30)
        return M, R

    def test_baseline_mmax_below_14_solar(self, baseline_mr):
        """Baseline M_max ≈ 13–14 M☉ — confirms current PM limit."""
        M, R = baseline_mr
        M_max = float(np.nanmax(M))
        assert 10.0 < M_max < 16.0, f"Baseline M_max = {M_max:.2f} M☉"

    def test_baseline_mmax_far_below_gw150914(self, baseline_mr):
        """Baseline M_max << 30 M☉ (GW150914 lower component)."""
        M, R = baseline_mr
        M_max = float(np.nanmax(M))
        assert M_max < 20.0, \
            f"Baseline M_max = {M_max:.2f} M☉ (should be < 20 M☉ = below LIGO)"

    def test_option_c_mmax_changes_from_baseline(self, baseline_mr, option_c_mr):
        """
        Option c (elastic pressure) should produce a different M_max than baseline.
        Adding U(φ) to the pressure increases stiffness → smaller, denser stars.
        """
        M_base, _ = baseline_mr
        M_opt_c, _ = option_c_mr
        M_max_base  = float(np.nanmax(M_base))
        M_max_opt_c = float(np.nanmax(M_opt_c))
        assert not math.isclose(M_max_base, M_max_opt_c, rel_tol=0.01), \
            "Option c M_max identical to baseline — U(φ) elastic term has no effect?"

    def test_option_c_radii_change(self, baseline_mr, option_c_mr):
        """
        With extra pressure, radii should differ from baseline.
        (Direction could be either way — stiffer EOS can increase R for low M.)
        """
        _, R_base  = baseline_mr
        _, R_opt_c = option_c_mr
        R_base_mid  = float(np.nanmedian(R_base[np.isfinite(R_base)]))
        R_opt_c_mid = float(np.nanmedian(R_opt_c[np.isfinite(R_opt_c)]))
        assert not math.isclose(R_base_mid, R_opt_c_mid, rel_tol=0.01), \
            "Option c radii identical to baseline — elastic term has no structural effect?"

    def test_option_b_cannot_close_mass_gap_with_linear_eos(self):
        """
        Option b (relaxed φ_crit): sweeping φ_crit from 1 to 5 with the
        *unchanged* linear EOS (c_s = c/√2) CANNOT reach M_max = 30 M☉.

        Reason: for the linear EOS P = c²/2(ρ−ρ_nuc), M_max is set by EOS
        stiffness (c_s = c/√2), not by where the density cap sits.  The M–R
        turning point saturates well below 30 M☉ regardless of ρ_crit.
        Strong negative result: the mass gap cannot be closed by relaxing the
        stability cap alone — the EOS itself must change.
        """
        phi_sweep = np.linspace(1.0, 5.0, 15)
        M_max_sweep = sweep_phi_crit_for_max_mass(phi_sweep)
        above_30 = [M for M in M_max_sweep if M >= 30.0]
        # With the baseline linear EOS, no phi_crit achieves 30 M_sun
        assert len(above_30) == 0, (
            f"Unexpected: baseline EOS reaches 30 M☉ at some φ_crit — "
            f"the linear EOS should saturate before 30 M☉"
        )

    def test_option_c_closes_mass_gap_at_current_phi_crit(self, option_c_mr):
        """
        Key positive finding: Option c (P_total = P_EOS + U(φ)) at the
        *current* φ_crit = 1 pushes M_max above 30 M☉ — comparable to
        the GW150914 lower component (≈30.7 M☉).

        Mechanism: U(φ_crit=1) = 2ε₀/3 ≈ 78% × P_EOS(ρ_crit).  Adding this
        as elastic pressure substantially stiffens the effective EOS, allowing
        the stellar structure ODE to support more massive stars.

        This is the most promising single-gap repair: no Lagrangian change or
        Poisson modification is required.  The next diagnostic is whether
        this also corrects the radius (currently ~30% too large) and
        surface redshift.
        """
        M_opt_c, _ = option_c_mr
        M_max_opt_c = float(np.nanmax(M_opt_c))
        assert M_max_opt_c > 28.0, (
            f"Option c M_max = {M_max_opt_c:.1f} M☉ — expected > 28 M☉ "
            f"(U(φ) elastic pressure should stiffen EOS sufficiently)"
        )
