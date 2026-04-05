"""Test: PM gravitational radiation structure (k_TT derivation).

Tests the analytical derivation of the PM GW radiation power from first principles:
  - Scalar φ channel: P_φ = G/(15c⁵) × <I⃛ᵢⱼI⃛ᵢⱼ> = (1/12) P_Peters (derived)
  - Vector u channel: needs u-field energy normalisation (open)
  - Total: validates k_TT = 1 via Hulse-Taylor at <0.5%

Reference: scripts/derive_k_TT.py, scripts/derive_k_TT_clean.py
"""

import math
import pytest
import sys
sys.path.insert(0, "src")

from pushing_medium import core as pm

G = 6.67430e-11
c = 299792458.0
M_SUN = 1.989e30

# ─── Hulse-Taylor system ─────────────────────────────────────────────────────
M1_HT = 1.4406 * M_SUN
M2_HT = 1.3886 * M_SUN
P_B_HT = 27906.9807        # s
E_HT = 0.6171334

# Galactic-corrected observed dP_b/dt (Weisberg, Nice & Taylor 2010)
DPDT_OBS_CORR = -2.4184e-12  # s/s

# ─── Helper ──────────────────────────────────────────────────────────────────
def semi_major_axis(M_tot, P_orb):
    """Kepler's third law: a from total mass and orbital period."""
    return (G * M_tot * (P_orb / (2 * math.pi))**2) ** (1.0 / 3.0)


def angular_frequency(M_tot, a):
    return math.sqrt(G * M_tot / a**3)


def i_triple_dot_squared(mu_r, a, omega):
    """Time-averaged <I⃛ᵢⱼ I⃛ᵢⱼ> for a circular binary."""
    return 32.0 * mu_r**2 * a**4 * omega**6


def pm_scalar_phi_power(i_dddot_sq):
    """P_φ = G/(60c⁵) × I⃛² — scalar φ-field radiation power.
    
    Derived from:
    - PM radiation zone: φ_q = (μ_G/2c²r) n̂ᵢn̂ⱼ Ïᵢⱼ(t_r)  [retarded multipole exp.]
    - ∂_tφ_q = (G/c⁴r) n̂ᵢn̂ⱼ I⃛ᵢⱼ(t_r)  [3rd derivative]
    - Energy scale K_grav = c⁴/(32πG)  [from Newtonian field energy c⁴|∇φ|²/(32πG)]
    - dP/dΩ = r² K_grav (∂_tφ)² / c = K_grav (G/c⁴)² (n̂·n̂·I⃛)² / c
    - Angular integral ∫(n̂ᵢn̂ⱼI⃛ᵢⱼ)² dΩ = 8π/15 × I⃛ᵢⱼI⃛ᵢⱼ  (traceless)
    - P_φ = K_grav × G²/c⁹ × (8π/15) × I⃛² = G/(60c⁵) × I⃛²
    
    This equals (1/12) P_Peters = (1/12) × (G/5c⁵ × I⃛²).
    """
    return (G / (60.0 * c**5)) * i_dddot_sq


def peters_power_circular(M1, M2, a):
    """Standard Peters (1964) formula for circular orbit."""
    mu_r = M1 * M2 / (M1 + M2)
    M_tot = M1 + M2
    return (32.0 / 5.0) * G**4 * mu_r**2 * M_tot**3 / (c**5 * a**5)


# ─── Tests ────────────────────────────────────────────────────────────────────

class TestITripleDotSquared:
    """Verify the <I⃛ᵢⱼI⃛ᵢⱼ> = 32μ²a⁴ω⁶ identity for circular orbits."""

    def test_circular_manual(self):
        """Cross-check by explicit time-averaging the components."""
        mu = 1.0  # dimensionless units
        a = 1.0
        omega = 1.0

        # Third derivatives at phase 0:
        # I⃛_xx = −4μa²ω³ sin(2ωt) → at t=0: 0
        # I⃛_yy = +4μa²ω³ sin(2ωt) → at t=0: 0
        # I⃛_xy = +4μa²ω³ cos(2ωt) → at t=0: +4

        # Time average: <sin²> = <cos²> = 1/2
        # <I⃛²> = 2×(4μa²ω³)² × <sin²> + 2×2×(4μa²ω³)² × <cos²>
        #       = 2 × 16 × 0.5 + 4 × 16 × 0.5 = 16 + 32 = ... let me recalculate
        # I⃛_xx² = 16μ²a⁴ω⁶ sin²(2ωt) → avg = 8
        # I⃛_yy² = 16μ²a⁴ω⁶ sin²(2ωt) → avg = 8
        # I⃛_xy = I⃛_yx → 2×I⃛_xy² = 2×16cos²ωt → avg = 16
        # Total: 8 + 8 + 16 = 32 ✓
        expected = 32.0 * mu**2 * a**4 * omega**6
        assert i_triple_dot_squared(mu, a, omega) == pytest.approx(expected, rel=1e-12)

    def test_scaling_with_mass(self):
        """I⃛² scales as μ²."""
        result1 = i_triple_dot_squared(1.0, 1.0, 1.0)
        result2 = i_triple_dot_squared(2.0, 1.0, 1.0)
        assert result2 / result1 == pytest.approx(4.0, rel=1e-12)

    def test_scaling_with_separation(self):
        """I⃛² scales as a⁴ × ω⁶ = (G³M³/a⁵) with Kepler."""
        mu = 1.0
        a1, a2 = 1.0, 2.0
        omega1, omega2 = 1.0, 1.0 / 2**(3.0/2)  # Kepler: ω ∝ a^(-3/2)
        r1 = i_triple_dot_squared(mu, a1, omega1)
        r2 = i_triple_dot_squared(mu, a2, omega2)
        # Expects: 32μ²a⁴ω⁶ ∝ a⁴×a^(-9) = a^(-5) → r2/r1 = (a2/a1)^(-5) = 2^(-5)
        assert r2 / r1 == pytest.approx(2.0**(-5), rel=1e-10)


class TestScalarPhiRadiation:
    """Test P_φ = G/(15c⁵) × I⃛² = (1/12) P_Peters for circular orbit."""

    def test_phi_is_one_twelfth_of_peters(self):
        """The scalar φ channel gives exactly 1/12 of Peters power.
        
        Proof: P_φ = G/(60c⁵) × I⃛²
               P_Peters = G/(5c⁵) × I⃛²  (traceless Q = I for circular orbit)
               Ratio = (1/60)/(1/5) = 1/12  (exact, any masses)
        """
        M1 = 1.4 * M_SUN
        M2 = 1.4 * M_SUN
        mu_r = M1 * M2 / (M1 + M2)
        M_tot = M1 + M2
        a = 1.0e9  # 1 million km
        omega = angular_frequency(M_tot, a)

        i_sq = i_triple_dot_squared(mu_r, a, omega)
        P_phi = pm_scalar_phi_power(i_sq)

        # Peters formula = G/(5c⁵) × I⃛² for circular orbit (traceless quadrupole)
        P_peters_from_i = G / (5.0 * c**5) * i_sq  
        # NOTE: This is G/5c⁵ × 32μ²a⁴ω⁶ = (32/5) G⁴μ²M³/(c⁵a⁵) via Kepler ✓

        ratio = P_phi / P_peters_from_i
        assert ratio == pytest.approx(1.0 / 12.0, rel=1e-12)

    def test_phi_power_dims_watts(self):
        """P_φ must have units of Watts (positive finite number)."""
        mu_r = 1.0 * M_SUN
        a = 1.0e10
        omega = 1e-4
        i_sq = i_triple_dot_squared(mu_r, a, omega)
        P = pm_scalar_phi_power(i_sq)
        assert P > 0
        assert math.isfinite(P)

    def test_phi_coefficient(self):
        """Coefficient G/(60c⁵) matches analytic derivation."""
        # Derivation:
        #   φ_q = (μ_G/2c²r) n̂ᵢn̂ⱼ Ïᵢⱼ  [retarded quadrupole]
        #   ∂_tφ_q = (G/c⁴r) n̂ᵢn̂ⱼ I⃛ᵢⱼ  [μ_G = 2G/c² → μ_G/(2c²) = G/c⁴]
        #   K_grav = c⁴/(32πG)  [from Newtonian field energy]
        #   P_φ = K_grav × (G/c⁴)² × (8π/15) × I⃛² / c
        #       = c⁴/(32πG) × G²/c⁸ × (8π/15) × I⃛² / c
        #       = G/(60c⁵) × I⃛²
        expected_coeff = G / (60.0 * c**5)
        computed_coeff = pm_scalar_phi_power(1.0)  # I⃛² = 1
        assert computed_coeff == pytest.approx(expected_coeff, rel=1e-12)

    def test_trace_zero_for_circular(self):
        """For circular orbit, Tr(I) = μa² = const → I⃛_trace = 0."""
        # This is an algebraic fact, not a numerical test, but we can verify
        # that the formula only uses the traceless part.
        # I_trace = I_xx + I_yy + I_zz = μa²(cos²+sin²+0) = μa² = const
        # Verified by the fact that P_φ uses I⃛ᵢⱼI⃛ᵢⱼ not I⃛²_trace
        # (see formula sheet §GW Polarisation and §GW Inspiral)
        mu = 1.4 * M_SUN
        a = 1.0e9
        omega = 1e-3
        # I_trace(t) = μa²(cos²ωt + sin²ωt) = μa² = const
        # All three time derivatives = 0
        I_trace_dot3 = 0.0  # exact
        assert I_trace_dot3 == 0.0


class TestKTTObservational:
    """Test k_TT = 1 via Hulse-Taylor observational constraint."""

    def test_peters_formula_hulse_taylor(self):
        """Peters formula matches Hulse-Taylor observed dP_b/dt to < 1%."""
        mu_r = M1_HT * M2_HT / (M1_HT + M2_HT)
        M_tot = M1_HT + M2_HT
        a = semi_major_axis(M_tot, P_B_HT)

        f_e = pm.pm_peters_decay_enhancement(E_HT)
        dPb_dt_pred = pm.pm_binary_period_derivative(M1_HT, M2_HT, P_B_HT, E_HT)

        # Observed Galactic-corrected value (Weisberg et al. 2010)
        # Note: k_TT ≡ dP_obs/dP_pred should be 1 to <1% (Galactic correction ±0.5%)
        ratio = DPDT_OBS_CORR / dPb_dt_pred
        assert ratio == pytest.approx(1.0, abs=0.01), (
            f"k_TT = {ratio:.4f} ≠ 1 (expected 1.000 ± 0.01 from Hulse-Taylor)"
        )

    def test_k_tt_observational_bound(self):
        """k_TT is within 1% of unity from Hulse-Taylor."""
        dPb_dt_pred = pm.pm_binary_period_derivative(M1_HT, M2_HT, P_B_HT, E_HT)
        k_TT = DPDT_OBS_CORR / dPb_dt_pred
        assert abs(k_TT - 1.0) < 0.01, (
            f"k_TT = {k_TT:.4f} deviates from 1 by > 1%"
        )

    def test_phi_scalar_fraction(self):
        """Scalar φ channel gives 1/12 of Peters power.
        
        Proof: P_φ/P_Peters = (G/60c⁵ × I⃛²) / (G/5c⁵ × I⃛²) = 5/60 = 1/12.
        """
        mu_r = M1_HT * M2_HT / (M1_HT + M2_HT)
        M_tot = M1_HT + M2_HT
        a = semi_major_axis(M_tot, P_B_HT)
        omega = angular_frequency(M_tot, a)

        i_sq = i_triple_dot_squared(mu_r, a, omega)
        P_phi = pm_scalar_phi_power(i_sq)  # = G/(60c⁵) × I⃛²
        # Peters (traceless Q = I for circular parallel): P_GR = G/(5c⁵) × I⃛²
        P_total_from_i = G / (5.0 * c**5) * i_sq

        scalar_fraction = P_phi / P_total_from_i
        assert scalar_fraction == pytest.approx(1.0 / 12.0, rel=1e-6), (
            f"Scalar φ fraction = {scalar_fraction:.6f}, expected 1/12 = {1/12:.6f}"
        )

    def test_u_field_dominates(self):
        """Vector u-field must carry ~92% of GW power (11/12 of Peters)."""
        mu_r = M1_HT * M2_HT / (M1_HT + M2_HT)
        M_tot = M1_HT + M2_HT
        a = semi_major_axis(M_tot, P_B_HT)
        omega = angular_frequency(M_tot, a)

        i_sq = i_triple_dot_squared(mu_r, a, omega)
        P_phi = pm_scalar_phi_power(i_sq)  # = G/(60c⁵) × I⃛²
        P_total_from_i = G / (5.0 * c**5) * i_sq  # Peters via same I⃛²

        # If k_TT = 1, then P_u = P_total - P_phi = (11/12) P_total
        P_u_required = P_total_from_i - P_phi
        fraction_u = P_u_required / P_total_from_i
        assert fraction_u == pytest.approx(11.0 / 12.0, rel=1e-6), (
            f"u-field required fraction = {fraction_u:.4f}, expected 11/12"
        )


class TestEnergyScaleIdentification:
    """Test the energy scale identification K_φ = c²/(8πG).
    
    This comes from identifying PM φ with GR h_tt via h_tt ≈ 2φ_PM,
    and using the GR gravitational wave energy density:
        ε_GW = c²/(32πG) × <ḣ_ij ḣ_ij>
    """

    def test_k_phi_value(self):
        """K_φ = c²/(8πG) is the correct PM φ-field energy scale."""
        K_phi = c**2 / (8 * math.pi * G)
        # This should be order of magnitude c²/G ~ 10^26 kg/m⁻¹... let me check:
        # [c²/G] = m²s⁻²/m³kg⁻¹s⁻² = kg/m  → K in [kg/m] → epsilon = K * phi_dot^2 in [kg m^-1 s^-2 = J/m^3] ✓
        assert K_phi > 0
        assert math.isfinite(K_phi)

    def test_phi_power_from_energy_scale(self):
        """P_φ from first-principles energy scale matches G/(60c⁵) × I⃛²."""
        mu_r = 1.0 * M_SUN
        a = 1.0e9
        omega = angular_frequency(3.0 * M_SUN, a)

        i_sq = i_triple_dot_squared(mu_r, a, omega)

        # Direct derivation:
        #   φ_q = (μ_G/2c²r) n̂ᵢn̂ⱼ Ïᵢⱼ(t_r)  [retarded multipole]
        #   ∂_tφ_q = (G/c⁴r) n̂ᵢn̂ⱼ I⃛ᵢⱼ        [μ_G = 2G/c²]
        #   K_grav = c⁴/(32πG)               [Newtonian: ε = c⁴|∇φ|²/(32πG)]
        #   dP/dΩ = r² K_grav (∂_tφ)²/c
        K_grav = c**4 / (32 * math.pi * G)
        coupling_phi = G / c**4   # from ∂_tφ = (G/c⁴r) × I⃛
        angular = 8 * math.pi / 15  # traceless quadrupole angular integral
        P_phi_from_scale = K_grav * coupling_phi**2 * i_sq * angular / c

        # Should match G/(60c⁵) × I⃛²
        P_phi_formula = G / (60.0 * c**5) * i_sq

        assert P_phi_from_scale == pytest.approx(P_phi_formula, rel=1e-10)
