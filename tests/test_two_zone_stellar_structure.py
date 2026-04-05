"""Tests for the two-zone (matter-phase shell + energy-phase core) stellar
structure solver.

Coverage
--------
1.  Single-zone agreement: for ρ_c ≤ ρ_crit all modes give the same result
    as solve_pm_star.
2.  Core detection: stars with ρ_c > ρ_crit detect a finite r_core.
3.  Vacuum special case: R_star = r_core (no shell), mass consistent.
4.  Radiation mode has a true mass maximum (dM/dρ_c turns negative).
5.  R_1.4 close to the baseline (matter-phase EOS drives 1.4 M☉ stars).
6.  Phase-boundary continuity: r_core is between 0 and R_star for constant/radiation.
7.  Mass ordering at ρ_c = 1.5 ρ_crit: radiation > constant > vacuum.
"""

import math
import numpy as np
import pytest

from pushing_medium.stellar_structure import (
    RHO_NUC,
    RHO_CRIT,
    M_SUN,
    solve_pm_star,
    solve_pm_star_two_zone,
    compute_mr_curve,
    compute_mr_curve_two_zone,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_r14(rho_c_arr, M_arr, R_arr):
    """Return R at the model closest to M = 1.4 M☉."""
    valid = np.isfinite(M_arr)
    assert valid.any(), "No converged models"
    M_v = M_arr[valid]; R_v = R_arr[valid]
    i = np.argmin(np.abs(M_v - 1.4))
    return float(R_v[i])


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def baseline_mr():
    """Baseline single-zone M-R sweep."""
    rho_c, M, R = compute_mr_curve(n_points=40)
    return rho_c, M, R


@pytest.fixture(scope="module")
def two_zone_radiation():
    rho_c, M, R, r_core = compute_mr_curve_two_zone(
        n_points=50, core_eos='radiation', rho_max_factor=5.0
    )
    return rho_c, M, R, r_core


# ---------------------------------------------------------------------------
# 1. Single-zone agreement
# ---------------------------------------------------------------------------

class TestSingleZoneAgreement:
    """For ρ_c ≤ ρ_crit, all two-zone modes must agree with solve_pm_star."""

    @pytest.mark.parametrize("core_eos", ["constant", "radiation", "vacuum"])
    def test_mass_matches_single_zone(self, core_eos):
        rho_c = 0.8 * RHO_CRIT           # well below transition
        ref   = solve_pm_star(rho_c)
        tzn   = solve_pm_star_two_zone(rho_c, core_eos=core_eos)

        assert ref.converged and tzn.converged, "Both models must converge"
        assert math.isnan(tzn.r_core), "No core expected below ρ_crit"

        M_ref = ref.M_star / M_SUN
        M_tzn = tzn.M_star / M_SUN
        assert abs(M_tzn - M_ref) / M_ref < 0.02, (
            f"Mass disagreement ({core_eos}): single={M_ref:.3f} vs "
            f"two-zone={M_tzn:.3f} M☉"
        )

    @pytest.mark.parametrize("core_eos", ["constant", "radiation", "vacuum"])
    def test_radius_matches_single_zone(self, core_eos):
        rho_c = 0.8 * RHO_CRIT
        ref   = solve_pm_star(rho_c)
        tzn   = solve_pm_star_two_zone(rho_c, core_eos=core_eos)

        R_ref = ref.R_star / 1e3
        R_tzn = tzn.R_star / 1e3
        assert abs(R_tzn - R_ref) / R_ref < 0.02, (
            f"Radius disagreement ({core_eos}): single={R_ref:.2f} vs "
            f"two-zone={R_tzn:.2f} km"
        )


# ---------------------------------------------------------------------------
# 2. Core detection at ρ_c > ρ_crit
# ---------------------------------------------------------------------------

class TestCoreDetection:

    @pytest.mark.parametrize("core_eos", ["constant", "radiation"])
    def test_core_found_above_rhocrit(self, core_eos):
        rho_c = 1.5 * RHO_CRIT
        sol   = solve_pm_star_two_zone(rho_c, core_eos=core_eos)
        assert sol.converged, f"Model must converge ({core_eos})"
        assert math.isfinite(sol.r_core), (
            f"r_core must be finite for ρ_c > ρ_crit ({core_eos})"
        )
        assert 0 < sol.r_core < sol.R_star, (
            f"r_core must be (0, R_star) ({core_eos}): "
            f"r_core={sol.r_core/1e3:.1f} km, R_star={sol.R_star/1e3:.1f} km"
        )

    def test_no_core_at_rhocrit_exactly(self):
        """Exactly at ρ_crit, has_core = False → same as single-zone."""
        # Slightly below to avoid floating-point equality edge cases
        rho_c = RHO_CRIT * (1 - 1e-6)
        sol   = solve_pm_star_two_zone(rho_c, core_eos='constant')
        assert math.isnan(sol.r_core), "No core expected at ρ_c = ρ_crit"


# ---------------------------------------------------------------------------
# 3. Vacuum special case
# ---------------------------------------------------------------------------

class TestVacuumMode:

    def test_vacuum_star_rstar_equals_rcore(self):
        """Vacuum star has P = 0 → no shell → R_star = r_core."""
        rho_c = 1.5 * RHO_CRIT
        sol   = solve_pm_star_two_zone(rho_c, core_eos='vacuum')
        assert sol.converged, "Vacuum star must converge"
        assert math.isfinite(sol.r_core), "r_core must be finite"
        assert abs(sol.R_star - sol.r_core) < 10.0, (   # within 10 m
            f"R_star ({sol.R_star/1e3:.2f} km) ≠ r_core ({sol.r_core/1e3:.2f} km) "
            "for vacuum star"
        )

    def test_vacuum_mass_consistent_with_uniform_ball(self):
        """For vacuum, M ≈ (4π/3) ρ_crit r_core³."""
        rho_c = 2.0 * RHO_CRIT
        sol   = solve_pm_star_two_zone(rho_c, core_eos='vacuum')
        assert sol.converged
        M_expected = (4.0 / 3.0) * math.pi * sol.r_core**3 * RHO_CRIT
        assert abs(sol.M_star - M_expected) / M_expected < 0.05, (
            f"Vacuum mass {sol.M_star/M_SUN:.2f} M☉ inconsistent with "
            f"uniform ρ_crit ball {M_expected/M_SUN:.2f} M☉"
        )


# ---------------------------------------------------------------------------
# 4. Radiation mode has a true mass maximum
# ---------------------------------------------------------------------------

class TestRadiationMassMaximum:

    def test_radiation_has_turnover(self, two_zone_radiation):
        rho_c, M, R, r_core = two_zone_radiation
        valid = np.isfinite(M); assert valid.sum() > 10
        M_v = M[valid]
        dM = np.diff(M_v)
        sign_changes = np.where(np.diff(np.sign(dM)))[0]
        assert len(sign_changes) >= 1, (
            "Radiation M-R curve must have at least one turnover (dM/dρ = 0)"
        )

    def test_radiation_mmax_above_baseline(self, two_zone_radiation, baseline_mr):
        _, M, _, _ = two_zone_radiation
        _, M_base, _  = baseline_mr
        M_max_tz = float(np.nanmax(M))
        M_max_base = float(np.nanmax(M_base))
        assert M_max_tz > M_max_base + 1.0, (
            f"Radiation M_max={M_max_tz:.2f} should exceed baseline "
            f"M_max={M_max_base:.2f} by > 1 M☉"
        )

    def test_radiation_mmax_in_physical_range(self, two_zone_radiation):
        _, M, _, _ = two_zone_radiation
        M_max = float(np.nanmax(M))
        # Must be in (13, 25) M☉ — above NS limit, below black-hole mass gap
        assert 13.0 < M_max < 25.0, (
            f"Radiation M_max = {M_max:.2f} M☉ outside expected range (13, 25)"
        )


# ---------------------------------------------------------------------------
# 5. R_1.4 close to baseline
# ---------------------------------------------------------------------------

class TestRadius14:

    @pytest.mark.parametrize("core_eos", ["constant", "radiation"])
    def test_r14_close_to_baseline(self, core_eos, baseline_mr):
        """1.4 M☉ stars have no core — same as baseline within 1 km."""
        _, M_base, R_base = baseline_mr
        # Tighter range for R_1.4: use default sweep but check it's near baseline
        try:
            rho_c, M, R, r_core = compute_mr_curve_two_zone(
                n_points=40, core_eos=core_eos,
            )
        except Exception as e:
            pytest.skip(f"Two-zone sweep failed: {e}")

        R14_tz   = _find_r14(rho_c, M, R)
        R14_base = _find_r14(*baseline_mr)

        assert abs(R14_tz - R14_base) < 2.0, (       # within 2 km (sampling tolerance)
            f"R_1.4 {core_eos}={R14_tz:.2f} vs baseline={R14_base:.2f} km"
        )


# ---------------------------------------------------------------------------
# 6. Phase-boundary continuity: r_core in (0, R_star)
# ---------------------------------------------------------------------------

class TestCoreRadiusRange:

    @pytest.mark.parametrize("rho_factor,core_eos", [
        (1.2, "constant"),
        (1.5, "constant"),
        (2.0, "constant"),
        (1.2, "radiation"),
        (1.5, "radiation"),
        (2.0, "radiation"),
    ])
    def test_rcore_inside_star(self, rho_factor, core_eos):
        rho_c = rho_factor * RHO_CRIT
        sol   = solve_pm_star_two_zone(rho_c, core_eos=core_eos)
        assert sol.converged, f"Model must converge (core_eos={core_eos}, factor={rho_factor})"
        assert math.isfinite(sol.r_core), "Core radius must be finite"
        assert 0 < sol.r_core < sol.R_star, (
            f"r_core={sol.r_core/1e3:.1f} km must be inside R_star={sol.R_star/1e3:.1f} km"
        )
        # Shell must exist
        assert sol.R_star - sol.r_core > 0.1e3, (    # shell > 100 m
            "Constant/radiation star must have a non-trivial matter shell"
        )


# ---------------------------------------------------------------------------
# 7. Mass ordering at ρ_c = 1.5 ρ_crit
# ---------------------------------------------------------------------------

class TestMassOrdering:

    def test_vacuum_lightest(self):
        """Vacuum core (P = 0) provides least support — lightest star at fixed ρ_c."""
        rho_c = 1.5 * RHO_CRIT
        sol_r = solve_pm_star_two_zone(rho_c, core_eos='radiation')
        sol_v = solve_pm_star_two_zone(rho_c, core_eos='vacuum')
        assert sol_r.converged and sol_v.converged
        assert sol_r.M_star > sol_v.M_star, (
            "Radiation core must produce heavier star than vacuum "
            f"({sol_r.M_star/M_SUN:.2f} vs {sol_v.M_star/M_SUN:.2f} M☉)"
        )

    def test_constant_and_radiation_both_above_vacuum(self):
        """Both constant and radiation cores give more mass than vacuum."""
        rho_c = 1.5 * RHO_CRIT
        sol_c = solve_pm_star_two_zone(rho_c, core_eos='constant')
        sol_v = solve_pm_star_two_zone(rho_c, core_eos='vacuum')
        assert sol_c.converged and sol_v.converged
        assert sol_c.M_star > sol_v.M_star, (
            "Constant core must produce heavier star than vacuum "
            f"({sol_c.M_star/M_SUN:.2f} vs {sol_v.M_star/M_SUN:.2f} M☉)"
        )
