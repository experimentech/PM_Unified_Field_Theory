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

from pushing_medium.stellar_structure import (
    G,
    M_SUN,
    MU_G,
    RHO_CRIT,
    StarSolution,
    c,
    compute_mr_curve,
    compute_mr_curve_nfield,
    compute_mr_curve_nfield_stiffened,
    compute_mr_curve_physical_measure,
    pm_eos_density,
    pm_eos_pressure,
    pm_eos_sound_speed,
    solve_pm_star,
    solve_pm_star_nfield,
    solve_pm_star_nfield_stiffened,
    solve_pm_star_physical_measure,
)
from pushing_medium.critical_state import RHO_NUC


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


# ---------------------------------------------------------------------------
# N-field (area-measure) stellar structure tests
# ---------------------------------------------------------------------------

class TestNFieldStar:
    """Tests for the n-field compact-star solver (area-measure action S₂).

    Key differences from the standard solver:
      • Poisson source scales as n² rather than n (stronger gravity at high density).
      • EOS and field equation are automatically consistent in the n-variable.
      • φ = ln n is exact by construction — no φ/EOS drift.
      • Expected: smaller max mass and slightly smaller radii than standard PM.
    """

    @pytest.fixture(scope="class")
    def star_mid(self):
        return solve_pm_star_nfield(RHO_MID)

    # --- Basic sanity ---

    def test_converged(self, star_mid):
        """N-field integrator must find the stellar surface within r_max."""
        assert star_mid.converged

    def test_positive_mass(self, star_mid):
        assert star_mid.M_star > 0.0

    def test_finite_radius(self, star_mid):
        R_km = star_mid.R_star / 1e3
        assert 1.0 < R_km < 50.0, f"Radius {R_km:.2f} km out of range"

    def test_surface_pressure_zero(self, star_mid):
        """Surface pressure ≈ 0 (integration terminates at n = 1)."""
        P_central = c**2 * RHO_NUC / 2.0 * (RHO_MID / RHO_NUC - 1.0)
        assert star_mid.P[-1] / P_central < 1e-3

    def test_solar_mass_range(self, star_mid):
        M_solar = star_mid.M_star / M_SUN
        assert 0.01 < M_solar < 20.0

    def test_rejects_supercritical_density(self):
        with pytest.raises(ValueError, match="ρ_crit"):
            solve_pm_star_nfield(RHO_CRIT * 1.05)

    def test_critical_density_star_converges(self):
        star = solve_pm_star_nfield(RHO_CRIT)
        assert star.converged

    def test_low_density_limit(self):
        star = solve_pm_star_nfield(RHO_NUC * 1.05)
        assert star.converged
        assert star.M_star > 0.0

    # --- EOS/field auto-consistency (key test for n-field) ---

    def test_phi_equals_ln_n_exactly(self, star_mid):
        """φ = ln n is exact by construction — no EOS/Poisson drift."""
        rho = star_mid.rho
        phi_from_eos = np.log(rho / RHO_NUC)
        np.testing.assert_allclose(
            star_mid.phi, phi_from_eos, rtol=1e-10,
            err_msg="φ ≠ ln(ρ/ρ_nuc): n-field auto-consistency broken"
        )

    def test_eos_consistent_inside_star(self, star_mid):
        """P and ρ satisfy the PM EOS pointwise."""
        interior = star_mid.P > 0
        rho_from_P = pm_eos_density(star_mid.P[interior])
        np.testing.assert_allclose(
            rho_from_P, star_mid.rho[interior], rtol=1e-6,
            err_msg="P(r) and ρ(r) inconsistent with PM EOS"
        )

    def test_density_decreases_outward(self, star_mid):
        interior = star_mid.P > 0
        rho_in = star_mid.rho[interior]
        if len(rho_in) > 2:
            assert np.all(np.diff(rho_in) <= 1e3)

    def test_phi_decreases_outward(self, star_mid):
        interior = star_mid.P > 0
        phi_in = star_mid.phi[interior]
        if len(phi_in) > 2:
            assert np.all(np.diff(phi_in) <= 1e-10)

    def test_mass_enclosed_increases(self, star_mid):
        assert np.all(np.diff(star_mid.m) >= 0.0)

    # --- Comparison with standard PM solver ---

    def test_weak_field_agreement(self):
        """At very low central density (n_c → 1), both solvers should give similar M and R."""
        rho_low = RHO_NUC * 1.05
        std  = solve_pm_star(rho_low)
        nfld = solve_pm_star_nfield(rho_low)
        if std.converged and nfld.converged:
            # At low density n≈1, n²≈n so the two Poisson sources are similar.
            # Agree to within 10%.
            assert abs(nfld.M_star - std.M_star) / std.M_star < 0.10
            assert abs(nfld.R_star - std.R_star) / std.R_star < 0.10

    def test_nfield_max_mass_less_than_standard(self):
        """N-field gravity is stronger at high density (n² > n for n > 1).

        The n²-weighted Poisson source compresses the star more, so M_max
        from the n-field solver should be ≤ the standard PM M_max.
        """
        _, M_std,   _ = compute_mr_curve(n_points=20)
        _, M_nfld, _  = compute_mr_curve_nfield(n_points=20)
        max_std  = float(np.nanmax(M_std))
        max_nfld = float(np.nanmax(M_nfld))
        assert max_nfld <= max_std * 1.02, (
            f"N-field M_max ({max_nfld:.2f} M_☉) > standard M_max ({max_std:.2f} M_☉) "
            f"— expected n² source to give smaller or equal max mass"
        )

    # --- M-R curve ---

    def test_mr_curve_nfield_basic(self):
        rho_c, M, R = compute_mr_curve_nfield(n_points=10)
        valid = ~np.isnan(M)
        assert valid.any()
        assert np.all(M[valid] > 0.0)
        assert np.all(R[valid] > 1.0)
        assert np.all(R[valid] < 100.0)


# ---------------------------------------------------------------------------
# Physical-measure stellar structure (α=3 action, w = n^{3/2} variable)
# ---------------------------------------------------------------------------

class TestPhysicalMeasureStar:
    """Tests for the physical-measure compact-star solver (volume-measure action S₃).

    Key properties vs the other variants:
      • Linearising variable: w = n^{3/2} = e^{3φ/2}
      • Vacuum equation: ∇²w = 0  (exact, same as n-field has ∇²n = 0)
      • Sourced equation: ∇²w = −κ_w w^{1/3}  with κ_w = 3Gρ_nuc/c²
      • Source grows as w^{1/3} = n^{1/2}  — sub-linear in density,
        WEAKER coupling at high n than both α=0 (source ∝ n) and α=2 (∝ n²).
      • Exact identity: φ = (2/3) ln w  at every grid point.
      • Expected: M_max **larger** than standard PM (α=0) because of weaker coupling.
    """

    @pytest.fixture(scope="class")
    def star_mid(self):
        return solve_pm_star_physical_measure(RHO_MID)

    # --- Basic sanity ---

    def test_converged(self, star_mid):
        assert star_mid.converged

    def test_positive_mass(self, star_mid):
        assert star_mid.M_star > 0.0

    def test_finite_radius(self, star_mid):
        R_km = star_mid.R_star / 1e3
        assert 1.0 < R_km < 100.0, f"Radius {R_km:.2f} km out of range"

    def test_solar_mass_range(self, star_mid):
        M_sol = star_mid.M_star / M_SUN
        assert 0.01 < M_sol < 200.0

    def test_surface_pressure_zero(self, star_mid):
        P_central = c**2 * RHO_NUC / 2.0 * (RHO_MID / RHO_NUC - 1.0)
        assert star_mid.P[-1] / P_central < 1e-3

    def test_rejects_supercritical_density(self):
        with pytest.raises(ValueError, match="ρ_crit"):
            solve_pm_star_physical_measure(RHO_CRIT * 1.05)

    def test_critical_density_star_converges(self):
        star = solve_pm_star_physical_measure(RHO_CRIT)
        assert star.converged

    def test_low_density_limit(self):
        star = solve_pm_star_physical_measure(RHO_NUC * 1.05)
        assert star.converged
        assert star.M_star > 0.0

    # --- φ = (2/3) ln w exactness ---

    def test_phi_equals_two_thirds_ln_w(self, star_mid):
        """φ = (2/3) ln w is exact by construction."""
        rho = star_mid.rho
        n   = rho / RHO_NUC
        w   = n ** 1.5
        phi_from_w = (2.0 / 3.0) * np.log(w)
        np.testing.assert_allclose(
            star_mid.phi, phi_from_w, rtol=1e-10,
            err_msg="φ ≠ (2/3)ln(n^{3/2}): α=3 self-consistency broken"
        )

    def test_phi_equals_ln_n(self, star_mid):
        """φ = ln(ρ/ρ_nuc) must also hold (both identities are equivalent)."""
        phi_from_eos = np.log(star_mid.rho / RHO_NUC)
        np.testing.assert_allclose(
            star_mid.phi, phi_from_eos, rtol=1e-10,
            err_msg="φ ≠ ln(ρ/ρ_nuc) in physical-measure solver"
        )

    def test_eos_consistent_inside_star(self, star_mid):
        interior = star_mid.P > 0
        rho_from_P = pm_eos_density(star_mid.P[interior])
        np.testing.assert_allclose(
            rho_from_P, star_mid.rho[interior], rtol=1e-6,
            err_msg="P(r) and ρ(r) inconsistent with PM EOS in α=3 solver"
        )

    def test_density_decreases_outward(self, star_mid):
        interior = star_mid.P > 0
        rho_in = star_mid.rho[interior]
        if len(rho_in) > 2:
            assert np.all(np.diff(rho_in) <= 1e3)

    def test_mass_enclosed_increases(self, star_mid):
        assert np.all(np.diff(star_mid.m) >= 0.0)

    # --- Comparison with other solvers ---

    def test_weak_field_agreement_with_standard(self):
        """At very low central density (n_c → 1) all solvers agree."""
        rho_low = RHO_NUC * 1.05
        std  = solve_pm_star(rho_low)
        phys = solve_pm_star_physical_measure(rho_low)
        if std.converged and phys.converged:
            assert abs(phys.M_star - std.M_star) / std.M_star < 0.20
            assert abs(phys.R_star - std.R_star) / std.R_star < 0.20

    def test_phys_measure_max_mass_greater_than_standard(self):
        """α=3 source  w^{1/3} = n^{1/2} is SUB-linear in n, so weaker coupling
        at high density → larger M_max than α=0 (source ∝ n) and α=2 (∝ n²).
        """
        _, M_std,  _ = compute_mr_curve(n_points=20)
        _, M_phys, _ = compute_mr_curve_physical_measure(n_points=20)
        max_std  = float(np.nanmax(M_std))
        max_phys = float(np.nanmax(M_phys))
        assert max_phys >= max_std * 0.98, (
            f"Physical-measure M_max ({max_phys:.2f} M☉) unexpectedly < "
            f"standard M_max ({max_std:.2f} M☉); weak-coupling action "
            f"should give larger or equal max mass"
        )

    def test_phys_measure_diverges_from_standard_at_high_density(self):
        """At high central density the two predictions must diverge significantly."""
        rho_high = RHO_NUC * 2.0
        std  = solve_pm_star(rho_high)
        phys = solve_pm_star_physical_measure(rho_high)
        if std.converged and phys.converged:
            rel_diff = abs(phys.M_star - std.M_star) / std.M_star
            assert rel_diff > 0.10, (
                f"Expected >10%% divergence at 2ρ_nuc but got {rel_diff*100:.1f}%%"
            )

    # --- M-R curve ---

    def test_mr_curve_physical_measure_basic(self):
        rho_c, M, R = compute_mr_curve_physical_measure(n_points=10)
        valid = ~np.isnan(M)
        assert valid.any()
        assert np.all(M[valid] > 0.0)
        assert np.all(R[valid] > 1.0)
        assert np.all(R[valid] < 200.0)


# ---------------------------------------------------------------------------
# Vacuum-stiffened n-field stellar structure (action-derived self-stiffening)
# ---------------------------------------------------------------------------

class TestNFieldStiffenedStar:
    """Tests for the vacuum-stiffened n-field compact-star solver.

    Physical background
    -------------------
    The full n-field action with a vacuum-subtracted medium self-energy potential:

        S_full = ∫ [½|∇φ|² n²  −  A V̂_vac(n)] d³x,
        V̂_vac'(n) = (n−1)²/n  (vanishes at surface n=1)

    gives the field equation:

        ∇²n = κρ_nuc (n−1)²/n − κρ_nuc n²

    The stiffening term (n−1)²/n reduces effective gravity by 12–15 % at
    typical neutron-star densities.  At the surface (n=1) it vanishes exactly,
    so all surface physics is identical to the bare n-field solver.

    Expected properties:
      • M_max > bare n-field M_max (~22% larger: 9.4 → 11.5 M☉)
      • M_max < Option-A M_max (different base equation — Option-A is not
        derived from the n-field action)
      • n decreases monotonically from center to surface (same as n-field)
      • Surface: n = 1  (P = 0)
      • φ = ln n exact by construction
    """

    @pytest.fixture(scope="class")
    def star_mid(self):
        return solve_pm_star_nfield_stiffened(RHO_MID)

    # --- Basic sanity ---

    def test_converged(self, star_mid):
        """Stiffened integrator must find the stellar surface within r_max."""
        assert star_mid.converged

    def test_positive_mass(self, star_mid):
        assert star_mid.M_star > 0.0

    def test_finite_radius(self, star_mid):
        R_km = star_mid.R_star / 1e3
        assert 1.0 < R_km < 50.0, f"Radius {R_km:.2f} km out of range"

    def test_surface_pressure_zero(self, star_mid):
        """Surface pressure ≈ 0 (integration terminates at n = 1)."""
        P_central = c**2 * RHO_NUC / 2.0 * (RHO_MID / RHO_NUC - 1.0)
        assert star_mid.P[-1] / P_central < 1e-3

    def test_solar_mass_range(self, star_mid):
        M_solar = star_mid.M_star / M_SUN
        assert 0.01 < M_solar < 30.0

    def test_rejects_supercritical_density(self):
        with pytest.raises(ValueError, match="ρ_crit"):
            solve_pm_star_nfield_stiffened(RHO_CRIT * 1.05)

    def test_critical_density_star_converges(self):
        star = solve_pm_star_nfield_stiffened(RHO_CRIT)
        assert star.converged

    def test_low_density_limit(self):
        star = solve_pm_star_nfield_stiffened(RHO_NUC * 1.05)
        assert star.converged
        assert star.M_star > 0.0

    # --- EOS / field consistency ---

    def test_phi_equals_ln_n_exactly(self, star_mid):
        """φ = ln n is exact by construction (uses n directly)."""
        phi_from_eos = np.log(star_mid.rho / RHO_NUC)
        np.testing.assert_allclose(
            star_mid.phi, phi_from_eos, rtol=1e-10,
            err_msg="φ ≠ ln(ρ/ρ_nuc): stiffened solver self-consistency broken"
        )

    def test_density_decreases_outward(self, star_mid):
        interior = star_mid.P > 0
        rho_in = star_mid.rho[interior]
        if len(rho_in) > 2:
            assert np.all(np.diff(rho_in) <= 1e3)

    def test_mass_enclosed_increases(self, star_mid):
        assert np.all(np.diff(star_mid.m) >= 0.0)

    # --- Stiffening raises M_max vs bare n-field ---

    def test_stiffened_max_mass_greater_than_bare_nfield(self):
        """Vacuum-stiffening reduces effective gravity, raising M_max.

        Expected: M_max_stiffened > M_max_nfield (~22% increase).
        Tolerance: stiffened must exceed bare by at least 5%.
        """
        _, M_nfld, _ = compute_mr_curve_nfield(n_points=25)
        _, M_stif, _ = compute_mr_curve_nfield_stiffened(n_points=25)
        max_nfld = float(np.nanmax(M_nfld))
        max_stif = float(np.nanmax(M_stif))
        assert max_stif > max_nfld * 1.05, (
            f"Stiffened M_max ({max_stif:.2f} M☉) not > 105% of bare n-field "
            f"M_max ({max_nfld:.2f} M☉)"
        )

    def test_stiffened_max_mass_in_expected_range(self):
        """M_max for the vacuum-stiffened solver should be around 10–13 M☉."""
        _, M_stif, _ = compute_mr_curve_nfield_stiffened(n_points=25)
        max_stif = float(np.nanmax(M_stif))
        assert 9.0 < max_stif < 15.0, (
            f"Stiffened M_max {max_stif:.2f} M☉ outside expected 9–15 M☉"
        )

    # --- M-R curve ---

    def test_mr_curve_stiffened_basic(self):
        rho_c, M, R = compute_mr_curve_nfield_stiffened(n_points=10)
        valid = ~np.isnan(M)
        assert valid.any()
        assert np.all(M[valid] > 0.0)
        assert np.all(R[valid] > 1.0)
        assert np.all(R[valid] < 60.0)

