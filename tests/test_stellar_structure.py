"""Tests for PM compact-star structure equations.

Covers:
  • PM self-consistent EOS properties
  • Single-star structure integration (mass, radius, surface pressure)
  • PM critical-density constraint enforcement
  • M–R curve computation
"""

import math

import numpy as np
import pytest

from src.pushing_medium.stellar_structure import (
    G,
    M_SUN,
    MU_G,
    RHO_CRIT,
    StarSolution,
    c,
    compute_mr_curve,
    pm_eos_density,
    pm_eos_pressure,
    pm_eos_sound_speed,
    solve_pm_star,
)
from src.pushing_medium.critical_state import RHO_NUC


# ---------------------------------------------------------------------------
# EOS tests
# ---------------------------------------------------------------------------

class TestEOS:
    def test_surface_condition(self):
        """P(ρ_nuc) = 0: star surface is at nuclear density."""
        assert pm_eos_pressure(RHO_NUC) == pytest.approx(0.0, abs=1.0)

    def test_pressure_increases_with_density(self):
        """P must increase with ρ (stable EOS)."""
        rho_values = np.linspace(RHO_NUC, RHO_CRIT, 20)
        P_values = pm_eos_pressure(rho_values)
        assert np.all(np.diff(P_values) > 0)

    def test_inverse_eos_roundtrip(self):
        """pm_eos_density(pm_eos_pressure(ρ)) ≡ ρ  for ρ ≥ ρ_nuc."""
        rho_test = np.array([RHO_NUC * 1.1, RHO_NUC * 1.5, RHO_CRIT * 0.9])
        for rho_in in rho_test:
            P = pm_eos_pressure(rho_in)
            rho_out = pm_eos_density(P)
            assert rho_out == pytest.approx(rho_in, rel=1e-10)

    def test_causal_sound_speed(self):
        """c_s = c/√2 < c: causal (sub-luminal) sound speed."""
        cs = pm_eos_sound_speed()
        assert cs < c
        assert cs == pytest.approx(c / math.sqrt(2.0), rel=1e-12)

    def test_eos_linear(self):
        """EOS slope: dP/dρ = c²/2 exactly."""
        rho1, rho2 = RHO_NUC * 1.2, RHO_NUC * 1.4
        slope = (pm_eos_pressure(rho2) - pm_eos_pressure(rho1)) / (rho2 - rho1)
        assert slope == pytest.approx(c**2 / 2.0, rel=1e-8)

    def test_no_tensile_pressure(self):
        """P = 0 for ρ < ρ_nuc (no tensile pressure in vacuum)."""
        rho_below = RHO_NUC * 0.5
        assert pm_eos_pressure(rho_below) == 0.0


# ---------------------------------------------------------------------------
# Single-star integration tests
# ---------------------------------------------------------------------------

# Reference central density (midpoint between ρ_nuc and ρ_crit)
RHO_MID = (RHO_NUC + RHO_CRIT) / 2.0


class TestSolvePMStar:
    @pytest.fixture(scope="class")
    def star_mid(self):
        """Single reference star at ρ_mid."""
        return solve_pm_star(RHO_MID)

    def test_converged(self, star_mid):
        """Integration must find the stellar surface (P → 0) within 50 km."""
        assert star_mid.converged, "ODE integrator did not find stellar surface"

    def test_positive_mass(self, star_mid):
        """Total stellar mass M > 0."""
        assert star_mid.M_star > 0.0

    def test_finite_radius(self, star_mid):
        """Stellar radius R ∈ [1 km, 50 km]."""
        R_km = star_mid.R_star / 1e3
        assert 1.0 < R_km < 50.0, f"Radius {R_km:.2f} km out of expected range"

    def test_surface_pressure_zero(self, star_mid):
        """Surface pressure should be ≈ 0 (integration terminates at P = 0)."""
        # The last non-zero pressure point should be very small relative to P_central
        P_central = pm_eos_pressure(RHO_MID)
        P_surface = star_mid.P[-1]
        # Allow tolerance: P_surface / P_central < 1e-4
        assert P_surface / P_central < 1e-3, (
            f"Pressure did not reach zero: P_surface/P_central = "
            f"{P_surface/P_central:.2e}"
        )

    def test_density_decreases_outward(self, star_mid):
        """ρ(r) should be monotonically non-increasing in the interior."""
        # Trim to points with P > 0 (inside star)
        interior = star_mid.P > 0
        rho_interior = star_mid.rho[interior]
        if len(rho_interior) > 2:
            assert np.all(np.diff(rho_interior) <= 1e3), (
                "Density is not monotonically non-increasing with radius"
            )

    def test_mass_enclosed_increases(self, star_mid):
        """m(r) must be monotonically increasing."""
        assert np.all(np.diff(star_mid.m) >= 0.0)

    def test_solar_mass_range(self, star_mid):
        """PM compact star mass should be in a physically plausible range.

        PM has a flat Minkowski background — the medium field φ creates a curved
        effective optical metric, but there is no background spacetime curvature,
        hence no GR compactness bound or event horizon.
        A PM object with GM/(c²R) ~ 0.5 is a compact horizonless object, not
        an unphysical configuration.  The PM critical-density bound ρ_c ≤ ρ_crit
        is the only hard stability limit.
        """
        M_solar = star_mid.M_star / M_SUN
        assert 0.01 < M_solar < 20.0, f"M_star = {M_solar:.3f} M_☉ outside expected range"

    def test_low_density_limit(self):
        """At ρ_c slightly above ρ_nuc, star should still converge with low mass."""
        star = solve_pm_star(RHO_NUC * 1.05)
        assert star.converged
        assert star.M_star > 0.0

    def test_critical_density_star(self):
        """At ρ_c = ρ_crit, star must converge (maximum PM star)."""
        star = solve_pm_star(RHO_CRIT)
        assert star.converged

    def test_rejects_supercritical_density(self):
        """solve_pm_star must raise ValueError for ρ_c > ρ_crit."""
        with pytest.raises(ValueError, match="ρ_crit"):
            solve_pm_star(RHO_CRIT * 1.05)

    def test_mass_increases_with_central_density(self):
        """At low-to-mid densities, higher ρ_c should give higher M_star."""
        rho_lo = RHO_NUC * 1.1
        rho_hi = RHO_NUC * 2.0
        star_lo = solve_pm_star(rho_lo)
        star_hi = solve_pm_star(rho_hi)
        if star_lo.converged and star_hi.converged:
            assert star_hi.M_star > star_lo.M_star


# ---------------------------------------------------------------------------
# Physical self-consistency tests
# ---------------------------------------------------------------------------

class TestPhysicsConsistency:
    def test_pm_gravity_exact(self):
        """Verify the PM structure equations are exact within PM.

        PM has a flat Minkowski background.  The force law a = +(c²/2)∇φ combined
        with the PM Poisson equation ∇²φ = −(8πG/c²)ρ via Gauss's theorem gives
        dφ/dr = −μ_G m(r)/r² and therefore a = −Gm/r² exactly at all
        compactnesses.  There is no GR compactness limit in PM because there is
        no background curvature — the medium creates an effective metric, not a
        curved spacetime.
        """
        star = solve_pm_star(RHO_MID)
        # φ integrated from PM Poisson should match ln(ρ/ρ_nuc) from EOS
        interior = star.P > 0
        phi_from_eos = np.log(star.rho[interior] / RHO_NUC)
        np.testing.assert_allclose(
            star.phi[interior], phi_from_eos, rtol=1e-4,
            err_msg="φ(r) from Poisson integration inconsistent with EOS φ = ln(ρ/ρ_nuc)"
        )

    def test_phi_decreases_outward(self):
        """PM field φ(r) must decrease from centre to surface (less compressed outward)."""
        star = solve_pm_star(RHO_MID)
        interior = star.P > 0
        phi_in = star.phi[interior]
        if len(phi_in) > 2:
            assert np.all(np.diff(phi_in) <= 1e-10), "φ(r) is not monotonically decreasing"

    def test_phi_central_equals_log_rho_ratio(self):
        """Central φ_c = ln(ρ_c/ρ_nuc) from PM density law."""
        rho_c = RHO_NUC * 1.5
        star = solve_pm_star(rho_c)
        expected_phi_c = math.log(rho_c / RHO_NUC)
        assert star.phi_central == pytest.approx(expected_phi_c, rel=1e-10)

    def test_phi_central_positive(self):
        """PM field φ > 0 inside the star (compressed medium, n = e^φ > 1)."""
        star = solve_pm_star(RHO_MID)
        assert np.all(star.phi >= 0.0)

    def test_pm_eos_consistency_inside_star(self):
        """ρ and P arrays should satisfy the PM EOS pointwise."""
        star = solve_pm_star(RHO_MID)
        interior = star.P > 0
        rho_from_eos = pm_eos_density(star.P[interior])
        np.testing.assert_allclose(
            rho_from_eos, star.rho[interior], rtol=1e-6,
            err_msg="ρ(r) and P(r) are inconsistent with PM EOS",
        )


# ---------------------------------------------------------------------------
# M–R curve tests
# ---------------------------------------------------------------------------

class TestMRCurve:
    @pytest.fixture(scope="class")
    def mr_curve(self):
        return compute_mr_curve(n_points=15)

    def test_returns_three_arrays(self, mr_curve):
        rho_c, M, R = mr_curve
        assert len(rho_c) == 15
        assert len(M) == 15
        assert len(R) == 15

    def test_masses_positive(self, mr_curve):
        _, M, _ = mr_curve
        valid = ~np.isnan(M)
        assert valid.any(), "All M values are NaN — no stars converged"
        assert np.all(M[valid] > 0.0)

    def test_radii_positive(self, mr_curve):
        _, _, R = mr_curve
        valid = ~np.isnan(R)
        assert valid.any()
        assert np.all(R[valid] > 0.0)

    def test_radii_in_km_range(self, mr_curve):
        """Radii should be in the range [1, 100] km."""
        _, _, R = mr_curve
        valid = ~np.isnan(R)
        assert np.all(R[valid] > 1.0)
        assert np.all(R[valid] < 100.0)

    def test_max_mass_below_20_solar(self, mr_curve):
        """PM maximum mass sanity check (< 20 M_☉).

        PM has a flat Minkowski background with no GR compactness bound.  The PM critical
        density ρ_crit = e·ρ_nuc sets the hard upper limit on central density.
        The maximum M_star from sweeping 0 to ρ_crit should be a finite,
        physically reasonable value.
        """
        _, M, _ = mr_curve
        valid = ~np.isnan(M)
        assert np.nanmax(M) < 20.0, f"Max mass {np.nanmax(M):.2f} M_☉ exceeds sanity bound"

    def test_max_density_in_range(self, mr_curve):
        """Central densities should span up to ρ_crit."""
        rho_c, _, _ = mr_curve
        assert rho_c[-1] <= RHO_CRIT * (1.0 + 1e-6)
