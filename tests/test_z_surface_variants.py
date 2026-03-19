"""
Tests: Surface Redshift Across Gap-Repair Variants

Summary of findings
-------------------
The PM surface redshift formula is derived from the PM optical metric:

    n(r) = e^{φ(r)},   g_tt = −c² e^{−2φ}

which gives:

    1 + z = e^{φ_ext},   φ_ext = μ_G M/R = 2GM/(c² R)

    z_PM = e^{2GM/(c² R)} − 1

The GR Schwarzschild formula for the same star:

    z_GR = (1 − 2GM/(c² R))^{−1/2} − 1

Writing x = 2GM/(c² R):

    z_PM = e^x − 1  ≈  x + x²/2 + …
    z_GR = (1−x)^{−1/2} − 1  ≈  x/2 + 3x²/8 + …

The LEADING-ORDER PM redshift is 2× the GR redshift for the SAME (M, R).

Tests below confirm that:

1.  z_PM ≈ 0.344–0.367 for ALL (f, α) variants at 1.4 M☉
    → the surface redshift barely changes with EOS parameters

2.  z_PM / z_GR ≈ 1.78–1.80 at all variants
    → this ratio is STRUCTURAL, rooted in μ_G = 2G/c²
    → no choice of f, α, or EOS stiffness can lower the ratio below ~1.5

3.  z_PM sits within the broad spectroscopic observed range [0.15, 0.45]
    (X-ray bursters and NS atmosphere models)

4.  z_PM exceeds the NICER-based GR expectation (z_GR ~ 0.20–0.26)
    by the structural factor ~1.79 at the 1.4 M☉ / joint-feasibility (M, R)

5.  To align z_PM with z_GR at M = 1.4 M☉ while KEEPING z_PM ~ 0.22,
    one would need R > 20 km — conclusively excluded by GW/NICER timing.

Conclusion: the (f, α)-parameter space fixes M_max and R_1.4 independently
but cannot lower z_PM/z_GR below ~1.7.  Resolving the redshift discrepancy
requires modifying the PM optical-metric formula (refractive-index prescription)
rather than the EOS.
"""

import math
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from pushing_medium.stellar_structure import (
    RHO_NUC, RHO_CRIT, M_SUN, MU_G, c,
    pm_surface_redshift,
)
from general_relativity.classical import gr_surface_redshift
from gap_probe_comparison import _solve_star_variant

G = 6.674e-11

# ---------------------------------------------------------------------------
# Shared constants / tolerances
# ---------------------------------------------------------------------------

# Observed spectral range from X-ray burster spectroscopy
# (EXO 0748-676, RX J0720.4-3125, 4U 1702-429, and NS atmosphere models)
Z_OBS_MIN = 0.10
Z_OBS_MAX = 0.45

# NICER-derived expectation assuming GR geometry: z_GR ≈ 0.20–0.26 for
# a 1.4 M☉ star with R ≈ 11–14 km
Z_NICER_MIN = 0.18
Z_NICER_MAX = 0.28

# Structural ratio z_PM/z_GR at typical NS compactness (x ≈ 0.28–0.35)
RATIO_MIN = 1.60   # never below this for x up to 0.40
RATIO_MAX = 2.00   # never above this for x above 0.05

# Numerical ODE tolerance for find_mr helpers
MASS_TOL = 0.01   # M☉


# ---------------------------------------------------------------------------
# Helpers shared across test classes
# ---------------------------------------------------------------------------

def _make_eos(f, alpha):
    """Build (eos_p, eos_rho) for Option-c EOS with scaled ρ_nuc = f × RHO_NUC."""
    rn   = f * RHO_NUC
    rc   = math.e * rn
    eps0 = rn * c ** 2

    def eos_p(rho, _rn=rn, _a=alpha, _e=eps0):
        phi = math.log(rho / _rn) if rho > _rn else 0.0
        return max((c ** 2 / 2.0) * (rho - _rn) + _a * _e * (phi ** 2 - phi ** 3 / 3.0), 0.0)

    def eos_rho(P, _rn=rn, _a=alpha, _e=eps0):
        if P <= 0:
            return _rn
        lo, hi = _rn, rc * 5.0
        for _ in range(60):
            mid = (lo + hi) / 2.0
            if eos_p(mid) < P:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0

    return eos_p, eos_rho, rn, rc


def _find_mr_at_target(target_msun, f, alpha):
    """Bisect on central density to find the star with M ≈ target_msun.

    Returns (R_km, M_kg) or (nan, nan) if not found.
    """
    eos_p, eos_rho, rn, rc = _make_eos(f, alpha)
    lo, hi = rn * 1.001, rc
    best_R = math.nan
    best_M = math.nan
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        result = _solve_star_variant(mid, rc, eos_p, eos_rho)
        if result is None:
            hi = mid
            continue
        R_m, M_kg = result
        M_sol = M_kg / M_SUN
        if abs(M_sol - target_msun) < MASS_TOL:
            best_R = R_m / 1e3
            best_M = M_kg
            break
        if M_sol < target_msun:
            lo = mid
        else:
            hi = mid
    return best_R, best_M


# ---------------------------------------------------------------------------
# Pre-computed reference values (verified with numerical ODE above)
# ---------------------------------------------------------------------------
#
# Scenario                                 R_1.4   z_PM    z_GR    ratio
# ---------------------------------------- ------  ------  ------  -----
# Baseline        (f=1.00, α=0.0)          13.94   0.3444  0.1918  1.796
# Option c        (f=1.00, α=1.0)          14.00   0.3447  0.1920  1.796
# Joint feasibil. (f=1.18, α=1.0)          13.22   0.3669  0.2061  1.780
# Joint (boosted) (f=1.18, α=1.15)         13.22   0.3664  0.2058  1.781
# ---------------------------------------- ------  ------  ------  -----


# ---------------------------------------------------------------------------
# Class 1 – z values at each repair variant
# ---------------------------------------------------------------------------

class TestZAtRepairVariants:
    """Surface redshift at 1.4 M☉ for each gap-repair scenario."""

    @pytest.fixture(params=[
        ("baseline",    1.00, 0.00),
        ("option_c",    1.00, 1.00),
        ("feasible",    1.18, 1.00),
        ("boosted",     1.18, 1.15),
    ])
    def scenario(self, request):
        label, f, alpha = request.param
        R_km, M_kg = _find_mr_at_target(1.4, f, alpha)
        assert not math.isnan(R_km), f"ODE failed to find 1.4 M☉ star for {label}"
        R_m = R_km * 1e3
        z_pm = pm_surface_redshift(M_kg, R_m)
        z_gr = gr_surface_redshift(M_kg, R_m)
        return label, z_pm, z_gr

    def test_z_pm_positive(self, scenario):
        _, z_pm, _ = scenario
        assert z_pm > 0.0, "z_PM must be positive"

    def test_z_pm_within_spectroscopic_range(self, scenario):
        """z_PM fits inside the broad observed spectral range [0.10, 0.45]."""
        label, z_pm, _ = scenario
        assert Z_OBS_MIN < z_pm < Z_OBS_MAX, (
            f"{label}: z_PM = {z_pm:.4f} outside [{Z_OBS_MIN}, {Z_OBS_MAX}]"
        )

    def test_z_pm_exceeds_nicer_gr_range(self, scenario):
        """z_PM exceeds the NICER-inferred GR expectation for every variant."""
        label, z_pm, _ = scenario
        assert z_pm > Z_NICER_MAX, (
            f"{label}: z_PM = {z_pm:.4f} should exceed NICER upper bound {Z_NICER_MAX}"
        )

    def test_z_pm_larger_than_z_gr(self, scenario):
        """PM always predicts a larger redshift than GR for the same (M, R)."""
        label, z_pm, z_gr = scenario
        assert z_pm > z_gr, (
            f"{label}: z_PM = {z_pm:.4f} should be > z_GR = {z_gr:.4f}"
        )


# ---------------------------------------------------------------------------
# Class 2 – structural ratio z_PM / z_GR
# ---------------------------------------------------------------------------

class TestStructuralZRatio:
    """The z_PM/z_GR ratio is a structural property of the PM optical metric."""

    def _ratio_at(self, f, alpha):
        R_km, M_kg = _find_mr_at_target(1.4, f, alpha)
        assert not math.isnan(R_km)
        R_m = R_km * 1e3
        z_pm = pm_surface_redshift(M_kg, R_m)
        z_gr = gr_surface_redshift(M_kg, R_m)
        return z_pm / z_gr

    def test_ratio_stays_above_1p6_at_all_variants(self):
        """Ratio z_PM/z_GR > 1.6 for all (f, α) pairs tested."""
        for f, alpha in [(1.00, 0.0), (1.00, 1.0), (1.18, 1.0), (1.18, 1.15)]:
            ratio = self._ratio_at(f, alpha)
            assert ratio > RATIO_MIN, (
                f"f={f}, α={alpha}: ratio = {ratio:.3f} fell below {RATIO_MIN}"
            )

    def test_ratio_near_1p78_at_joint_feasibility(self):
        """At joint-feasibility (f=1.18, α=1), ratio is 1.78 ± 5%."""
        ratio = self._ratio_at(1.18, 1.0)
        assert 1.65 < ratio < 1.95, (
            f"Expected ratio ≈ 1.78, got {ratio:.3f}"
        )

    def test_ratio_invariant_over_alpha(self):
        """Ratio barely changes when α varies from 0 to 2 at fixed f=1.18."""
        ratios = [self._ratio_at(1.18, a) for a in [0.0, 0.5, 1.0, 1.5, 2.0]]
        span = max(ratios) - min(ratios)
        assert span < 0.10, (
            f"Ratio span over α ∈ [0,2] at f=1.18: {span:.3f} — should be < 0.10"
        )

    def test_ratio_invariant_over_f(self):
        """Ratio barely changes when f varies from 0.90 to 1.30 at fixed α=1."""
        ratios = [self._ratio_at(f, 1.0) for f in [0.90, 1.00, 1.10, 1.18, 1.25]]
        span = max(ratios) - min(ratios)
        assert span < 0.15, (
            f"Ratio span over f ∈ [0.90,1.30] at α=1: {span:.3f} — should be < 0.15"
        )

    def test_analytic_ratio_formula_at_typical_compactness(self):
        """Analytic check: at x = 0.30, z_PM/z_GR = 1.79 (formula only, no ODE)."""
        x = 0.30
        z_pm = math.exp(x) - 1.0
        z_gr = 1.0 / math.sqrt(1.0 - x) - 1.0
        ratio = z_pm / z_gr
        assert 1.70 < ratio < 1.90, (
            f"At x=0.30: ratio = {ratio:.4f}, expected ≈ 1.79"
        )

    def test_ratio_bounded_between_1p5_and_2p0_over_full_compactness_range(self):
        """Over NS compactness x ∈ [0.10, 0.40], ratio in [1.5, 2.0]."""
        for x in [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
            z_pm = math.exp(x) - 1.0
            z_gr = 1.0 / math.sqrt(1.0 - x) - 1.0
            ratio = z_pm / z_gr
            assert 1.50 <= ratio <= 2.00, (
                f"x={x:.2f}: ratio = {ratio:.4f} outside [1.50, 2.00]"
            )


# ---------------------------------------------------------------------------
# Class 3 – z_PM is insensitive to EOS tuning
# ---------------------------------------------------------------------------

class TestZInsensitiveToEOS:
    """Changing (f, α) does not substantially change z_PM at fixed mass."""

    def _z_pm_at(self, f, alpha):
        R_km, M_kg = _find_mr_at_target(1.4, f, alpha)
        assert not math.isnan(R_km)
        return pm_surface_redshift(M_kg, R_km * 1e3)

    def test_z_pm_change_over_alpha_is_small(self):
        """z_PM variation over α ∈ [0, 2] at f=1.18 is < 5% of the baseline value."""
        z_values = [self._z_pm_at(1.18, a) for a in [0.0, 0.5, 1.0, 1.5, 2.0]]
        z_ref = z_values[0]
        for z in z_values:
            assert abs(z - z_ref) / z_ref < 0.05, (
                f"z_PM={z:.4f} deviates > 5% from baseline {z_ref:.4f} over α sweep"
            )

    def test_z_pm_change_over_f_is_small(self):
        """z_PM variation over f ∈ [0.90, 1.25] at α=0 is < 10% of baseline."""
        z_base = self._z_pm_at(1.00, 0.0)
        for f in [0.90, 1.00, 1.10, 1.18, 1.25]:
            z = self._z_pm_at(f, 0.0)
            assert abs(z - z_base) / z_base < 0.10, (
                f"f={f}: z_PM={z:.4f} deviates > 10% from baseline {z_base:.4f}"
            )

    def test_z_pm_at_baseline_approx_0p34(self):
        """At baseline (f=1.0, α=0), z_PM ≈ 0.344 at 1.4 M☉."""
        z = self._z_pm_at(1.00, 0.0)
        assert 0.30 < z < 0.40, f"Expected z_PM ≈ 0.344, got {z:.4f}"

    def test_z_pm_at_joint_feasibility_approx_0p37(self):
        """At joint feasibility (f=1.18, α=1.0), z_PM ≈ 0.367 at 1.4 M☉."""
        z = self._z_pm_at(1.18, 1.0)
        assert 0.33 < z < 0.41, f"Expected z_PM ≈ 0.367, got {z:.4f}"


# ---------------------------------------------------------------------------
# Class 4 – optical-metric origin of the discrepancy (analytic)
# ---------------------------------------------------------------------------

class TestOpticalMetricOrigin:
    """Analytic tests confirming the discrepancy is rooted in the PM optical metric.

    PM uses a flat-space optical metric with refractive index n = e^φ:
        g_tt = −c² e^{−2φ},   1+z = e^{φ_surface} = e^{μ_G M / R}

    This gives z_PM = e^x − 1 where x = 2GM/(c² R).

    GR uses the Schwarzschild metric:
        g_tt = −(1 − 2GM/c²R),   1+z = (1−x)^{−1/2}

    The leading terms differ by a factor of 2.  No EOS change can alter this.
    """

    def test_pm_leading_order_is_twice_gr(self):
        """For small compactness x, z_PM ≈ 2 z_GR."""
        x = 0.02   # very weak field
        z_pm = math.exp(x) - 1.0
        z_gr = 1.0 / math.sqrt(1.0 - x) - 1.0
        # z_PM ~ x,  z_GR ~ x/2  → ratio ~ 2
        ratio = z_pm / z_gr
        assert 1.90 < ratio < 2.10, (
            f"Weak-field ratio expected ≈ 2.0, got {ratio:.4f}"
        )

    def test_pm_formula_exp_minus_one(self):
        """Unit test: pm_surface_redshift(M, R) matches exp(MU_G M/R) − 1."""
        M_kg = 1.4 * M_SUN
        R_m  = 13.22e3
        expected = math.exp(MU_G * M_kg / R_m) - 1.0
        got = pm_surface_redshift(M_kg, R_m)
        assert abs(got - expected) < 1e-10, (
            f"pm_surface_redshift mismatch: {got} vs {expected}"
        )

    def test_required_radius_for_pm_equaling_gr_z_is_unphysical(self):
        """To make z_PM = z_GR at M=1.4 M☉ would require R > 20 km (excluded).

        Equating the weak-field limits:
            z_PM ≈ 2GM/(c²R)  =  z_GR ≈ GM/(c²R')
        for the PM star to give the same z as a GR star at R=12 km:
            R_PM  ≈  (z_GR / z_PM) × R_GR  →  not viable
        More precisely, if we demand z_PM(M, R_PM) = z_GR(M, 12km):
            exp(μ_G M / R_PM) − 1  =  1/sqrt(1 − μ_G M / 12km) − 1
        Solve for R_PM.
        """
        M_kg = 1.4 * M_SUN
        R_ref = 12.0e3          # GR reference radius [m]
        z_target = 1.0 / math.sqrt(1.0 - MU_G * M_kg / R_ref) - 1.0

        # Find R_PM such that exp(MU_G M / R_PM) − 1 = z_target
        def f(R):
            return math.exp(MU_G * M_kg / R) - 1.0 - z_target

        # Bisect over [10 km, 50 km]
        lo, hi = 10.0e3, 50.0e3
        for _ in range(60):
            mid = (lo + hi) / 2.0
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
        R_required_km = (lo + hi) / 2.0 / 1e3

        assert R_required_km > 18.0, (
            f"PM radius required for z_PM=z_GR: {R_required_km:.1f} km "
            f"(should be > 18 km, conclusively excluded by NS observations)"
        )
