"""
Tests for PM vs GR shadow radius predictions, cross-checked against EHT data
for M87* and Sgr A*.

Key physical results tested:
  - GR BH shadow: b = 3√3 GM/c² ≈ 2.598 R_s  (photon-sphere formula)
  - PM compact-star shadow: b = exp(2GM/c²R_star) × R_star (Bouguer/eikonal)
  - PM minimum shadow (R_star → R_s): b_min = e R_s ≈ 2.718 R_s  > GR shadow
  - PM always predicts shadows ≥ 4.6 % larger than GR BH for same mass
  - EHT observations: M87* θ=42±3 μas, Sgr A* θ=51.8±2.3 μas
  - EHT sources have masses >> PM max compact-object mass (~14 M_sun)
"""

import math
import pytest

from src.general_relativity.classical import (
    gr_photon_sphere_radius,
    gr_shadow_radius,
    gr_compact_star_shadow_radius,
)
from src.pushing_medium.core import pm_shadow_radius, G, c

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
M_SUN = 1.989e30          # kg
RAD_TO_MICRO_AS = 180.0 / math.pi * 3600.0 * 1.0e6   # 1 rad → μas

# EHT sources
M_M87  = 6.5e9  * M_SUN   # M87* mass [kg]
D_M87  = 16.8 * 3.086e22    # 16.8 Mpc in metres (1 Mpc = 3.086e22 m)
THETA_M87_OBS   = 42.0    # μas observed
THETA_M87_ERR   = 3.0     # μas 1σ

M_SGRA = 4.15e6 * M_SUN   # Sgr A* mass [kg]
D_SGRA = 8.178 * 3.086e19   # 8.178 kpc in metres (1 kpc = 3.086e19 m)
THETA_SGRA_OBS  = 51.8    # μas observed
THETA_SGRA_ERR  = 2.3     # μas 1σ

PM_M_MAX_SUN = 14.0        # conservative PM max compact-object mass [M_sun]


def shadow_angular_diameter_uas(b_shadow: float, distance: float) -> float:
    """Angular diameter θ = 2 b_shadow / D, converted to μas."""
    return 2.0 * b_shadow / distance * RAD_TO_MICRO_AS


# ===========================================================================
# GR shadow radius
# ===========================================================================
class TestGRShadowRadius:
    """gr_shadow_radius = 3√3 GM/c²"""

    def test_formula_numerical(self):
        M = M_SUN
        expected = 3.0 * math.sqrt(3.0) * G * M / (c * c)
        assert math.isclose(gr_shadow_radius(M), expected, rel_tol=1e-12)

    def test_scales_linearly_with_mass(self):
        assert math.isclose(
            gr_shadow_radius(2 * M_SUN) / gr_shadow_radius(M_SUN), 2.0, rel_tol=1e-12
        )

    def test_ratio_to_photon_sphere(self):
        """b_shadow = sqrt(3) × r_ps = sqrt(3) × 3GM/c²."""
        M = 10 * M_SUN
        ratio = gr_shadow_radius(M) / gr_photon_sphere_radius(M)
        assert math.isclose(ratio, math.sqrt(3.0), rel_tol=1e-12)

    def test_ratio_to_schwarzschild_radius(self):
        """b_shadow = 3√3/2 × R_s.  3√3/2 ≈ 2.598."""
        M = 10 * M_SUN
        R_s = 2.0 * G * M / (c * c)
        ratio = gr_shadow_radius(M) / R_s
        assert math.isclose(ratio, 1.5 * math.sqrt(3.0), rel_tol=1e-12)

    def test_m87_angular_diameter_within_eht_bounds(self):
        """GR shadow angular diameter for M87* must lie within the measured EHT value."""
        theta = shadow_angular_diameter_uas(gr_shadow_radius(M_M87), D_M87)
        assert abs(theta - THETA_M87_OBS) < 3 * THETA_M87_ERR, (
            f"GR θ_M87 = {theta:.1f} μas, expected {THETA_M87_OBS}±{3*THETA_M87_ERR} μas"
        )

    def test_sgra_angular_diameter_within_eht_bounds(self):
        """GR shadow angular diameter for Sgr A* must lie within the measured EHT value."""
        theta = shadow_angular_diameter_uas(gr_shadow_radius(M_SGRA), D_SGRA)
        assert abs(theta - THETA_SGRA_OBS) < 3 * THETA_SGRA_ERR, (
            f"GR θ_SgrA* = {theta:.1f} μas, expected {THETA_SGRA_OBS}±{3*THETA_SGRA_ERR} μas"
        )


# ===========================================================================
# GR compact-star shadow
# ===========================================================================
class TestGRCompactStarShadow:
    """gr_compact_star_shadow_radius = R_star / √(1 − R_s/R_star)"""

    def test_formula_numerical(self):
        M = 1.4 * M_SUN
        R_star = 12e3          # 12 km
        R_s = 2.0 * G * M / (c * c)
        expected = R_star / math.sqrt(1.0 - R_s / R_star)
        assert math.isclose(gr_compact_star_shadow_radius(M, R_star), expected, rel_tol=1e-12)

    def test_at_photon_sphere_equals_bh_shadow(self):
        """When R_star = 3GM/c² = r_ps, the Schwarzschild grazing-ray formula
        gives exactly the same b as the photon-sphere shadow formula."""
        M = 1.4 * M_SUN
        r_ps = 3.0 * G * M / (c * c)
        # grazing b at R_star = r_ps
        b_star = gr_compact_star_shadow_radius(M, r_ps)
        b_bh = gr_shadow_radius(M)
        assert math.isclose(b_star, b_bh, rel_tol=1e-10)

    def test_larger_star_has_larger_shadow(self):
        """Larger star (smaller compactness) gives larger GR shadow.
        b = R_star/sqrt(1-R_s/R_star) is monotonically increasing with R_star."""
        M = 1.4 * M_SUN
        b_large_star = gr_compact_star_shadow_radius(M, 15e3)   # 15 km
        b_small_star = gr_compact_star_shadow_radius(M, 10e3)   # 10 km
        assert b_large_star > b_small_star

    def test_typical_neutron_star(self):
        """Typical NS: M=1.4 M_sun, R=12 km → shadow impact parameter check."""
        M = 1.4 * M_SUN
        R_star = 12e3
        b = gr_compact_star_shadow_radius(M, R_star)
        R_s = 2.0 * G * M / (c * c)
        # Should be only modestly above R_star (low compactness)
        assert R_star < b < 5 * R_star


# ===========================================================================
# PM shadow radius
# ===========================================================================
class TestPMShadowRadius:
    """pm_shadow_radius = exp(2GM/c²R_star) × R_star"""

    def test_formula_numerical(self):
        M = 1.4 * M_SUN
        R_star = 12e3
        mu_G = 2.0 * G / (c * c)
        expected = math.exp(mu_G * M / R_star) * R_star
        assert math.isclose(pm_shadow_radius(M, R_star), expected, rel_tol=1e-12)

    def test_at_schwarzschild_radius_is_e_times_Rs(self):
        """At R_star = R_s = 2GM/c² (maximum compactness), b_PM = e × R_s."""
        M = 10.0 * M_SUN
        R_s = 2.0 * G * M / (c * c)
        b = pm_shadow_radius(M, R_s)
        assert math.isclose(b, math.e * R_s, rel_tol=1e-10)

    def test_minimum_shadow_exceeds_gr_bh_shadow(self):
        """PM minimum shadow (e R_s) is always larger than GR BH shadow (2.598 R_s)."""
        M = 10.0 * M_SUN
        R_s = 2.0 * G * M / (c * c)
        b_pm_min = pm_shadow_radius(M, R_s)     # most compact case
        b_gr_bh  = gr_shadow_radius(M)
        ratio = b_pm_min / b_gr_bh
        # ratio = 2e / (3√3) ≈ 1.046
        assert ratio > 1.0, "PM shadow must exceed GR BH shadow at max compactness"
        assert math.isclose(ratio, 2.0 * math.e / (3.0 * math.sqrt(3.0)), rel_tol=1e-8)

    def test_shadow_ratio_at_max_compactness(self):
        """PM/GR ratio = 2e/(3√3) ≈ 1.046."""
        M = 5.0 * M_SUN
        R_s = 2.0 * G * M / (c * c)
        ratio = pm_shadow_radius(M, R_s) / gr_shadow_radius(M)
        assert math.isclose(ratio, 2.0 * math.e / (3.0 * math.sqrt(3.0)), rel_tol=1e-8)

    def test_larger_pm_star_has_larger_shadow(self):
        """For R_star > R_s, larger star gives larger PM shadow.
        d/dR_star [exp(μ_G M/R_star) R_star] = exp(μ_G M/R_star)(1 - μ_G M/R_star) > 0
        for R_star > R_s (where μ_G M/R_star = 1), so b is increasing."""
        M = 1.4 * M_SUN
        b_large_star = pm_shadow_radius(M, 20e3)   # 20 km (less compact)
        b_small_star = pm_shadow_radius(M, 10e3)   # 10 km (more compact)
        assert b_large_star > b_small_star

    def test_shadow_approaches_star_radius_when_dilute(self):
        """For very extended star (R_star >> R_s), n(R_star) → 1 and b → R_star."""
        M = M_SUN
        R_very_large = 1e12   # 1 Gm — nearly flat spacetime
        b = pm_shadow_radius(M, R_very_large)
        assert math.isclose(b / R_very_large, 1.0, rel_tol=1e-3)

    def test_pm_shadow_larger_than_gr_star_shadow_for_typical_ns(self):
        """PM shadow > GR star shadow for typical neutron-star compactness."""
        M = 1.4 * M_SUN
        R_star = 12e3   # 12 km; R_s ≈ 4.1 km → k ≈ 2.9
        b_pm = pm_shadow_radius(M, R_star)
        b_gr = gr_compact_star_shadow_radius(M, R_star)
        assert b_pm > b_gr


# ===========================================================================
# EHT falsification via mass gap
# ===========================================================================
class TestEHTMassGapFalsification:
    """EHT supermassive BHs have masses >> PM maximum compact-object mass."""

    def test_m87_exceeds_pm_max_by_large_factor(self):
        """M87* mass is hundreds of millions of times PM M_max."""
        ratio = M_M87 / (PM_M_MAX_SUN * M_SUN)
        assert ratio > 1e8, f"M87*/PM_Mmax ratio = {ratio:.2e}, expected > 1e8"

    def test_sgra_exceeds_pm_max_by_large_factor(self):
        """Sgr A* mass is hundreds of thousands of times PM M_max."""
        ratio = M_SGRA / (PM_M_MAX_SUN * M_SUN)
        assert ratio > 1e5, f"Sgr A*/PM_Mmax ratio = {ratio:.2e}, expected > 1e5"

    def test_gr_m87_shadow_consistent_with_eht(self):
        """GR prediction for M87* shadow is within 2σ of EHT measurement."""
        theta = shadow_angular_diameter_uas(gr_shadow_radius(M_M87), D_M87)
        assert abs(theta - THETA_M87_OBS) < 2 * THETA_M87_ERR

    def test_gr_sgra_shadow_consistent_with_eht(self):
        """GR prediction for Sgr A* shadow is within 2σ of EHT measurement."""
        theta = shadow_angular_diameter_uas(gr_shadow_radius(M_SGRA), D_SGRA)
        assert abs(theta - THETA_SGRA_OBS) < 2 * THETA_SGRA_ERR

    def test_pm_max_compact_shadow_would_be_larger_than_eht(self):
        """Even the most compact PM object at EHT-source mass would cast a larger
        shadow than GR, but PM cannot form such objects (mass > M_max).
        Tests that the PM formula gives b > b_GR for equal mass."""
        for M in [M_SGRA, M_M87]:
            R_s = 2.0 * G * M / (c * c)
            b_pm  = pm_shadow_radius(M, R_s)      # most compact PM case
            b_gr  = gr_shadow_radius(M)
            assert b_pm > b_gr

    def test_pm_shadow_ratio_independent_of_mass(self):
        """The PM/GR shadow ratio at maximum compactness is a pure number
        2e/(3√3) ≈ 1.046, independent of mass."""
        expected_ratio = 2.0 * math.e / (3.0 * math.sqrt(3.0))
        for M_fac in [1.0, 5.0, 14.0, 1e6, 1e9]:
            M = M_fac * M_SUN
            R_s = 2.0 * G * M / (c * c)
            ratio = pm_shadow_radius(M, R_s) / gr_shadow_radius(M)
            assert math.isclose(ratio, expected_ratio, rel_tol=1e-8), (
                f"ratio={ratio} ≠ {expected_ratio} at M={M_fac} M_sun"
            )
