"""Surface gravitational redshift tests for PM compact stars.

Physical context
----------------
A photon emitted from the surface of a compact star and received at infinity
is gravitationally redshifted. The surface redshift z encodes the stellar
compactness.

PM formula (exact within PM, all compactness):
    z_PM  = e^{μ_G M/R} − 1   where μ_G = 2G/c²

GR formula (exact Schwarzschild):
    z_GR  = 1/√(1 − 2GM/c²R) − 1

Key difference: both equal ~GM/c²R at leading order... wait:
    z_PM  ≈ 2GM/c²R  (first-order)
    z_GR  ≈  GM/c²R  (first-order)

→ PM predicts roughly TWICE the surface redshift of GR for the same (M, R).
This is a large, observable difference testable with X-ray spectroscopy (NICER).

NICER observations used for reference bounds
--------------------------------------------
Riley et al. (2019)  PSR J0030+0451:  M = 1.34 ± 0.16 M_sun,  R = 12.71 km
Riley et al. (2021)  PSR J0740+6620:  M = 2.08 ± 0.07 M_sun,  R = 12.39 km

Spectral redshift measurements (X-ray bursts):
Cottam et al. (2002) EXO 0748-676: z_obs = 0.35  (disputed; used as indicative)
"""

import math
import pytest
import numpy as np

from pushing_medium.stellar_structure import (
    pm_surface_redshift,
    compute_surface_redshift_track,
    RHO_NUC,
    MU_G,
    G,
    c,
    M_SUN,
)
from general_relativity.classical import gr_surface_redshift

# ── Helpers ──────────────────────────────────────────────────────────────────
def _mu(M, R):
    """Dimensionless PM compactness: μ_G M/R = 2GM/c²R."""
    return MU_G * M / R


# ── Analytic formula tests ────────────────────────────────────────────────────
class TestPMSurfaceRedshiftFormula:
    def test_zero_compactness_gives_zero(self):
        """Massless star: no compactness → no redshift."""
        z = pm_surface_redshift(0.0, 1e4)
        assert z == pytest.approx(0.0, abs=1e-15)

    def test_positive(self):
        """Redshift must be positive for any positive M, R."""
        M = 1.4 * M_SUN
        R = 10e3
        assert pm_surface_redshift(M, R) > 0.0

    def test_weak_field_linear(self):
        """For small compactness x = μ_G M/R: z_PM ≈ x  (first order)."""
        M = 0.001 * M_SUN   # very low mass
        R = 100e3            # large radius → small compactness
        x = MU_G * M / R
        z = pm_surface_redshift(M, R)
        assert z == pytest.approx(x, rel=1e-3)

    def test_formula_value(self):
        """Spot check: M=1.4 M_sun, R=10 km → z = e^(μ_G M/R) - 1."""
        M = 1.4 * M_SUN
        R = 10e3
        expected = math.exp(MU_G * M / R) - 1.0
        assert pm_surface_redshift(M, R) == pytest.approx(expected, rel=1e-12)

    def test_increases_with_mass(self):
        """Higher mass at fixed radius → more compactness → larger redshift."""
        R = 10e3
        z1 = pm_surface_redshift(1.0 * M_SUN, R)
        z2 = pm_surface_redshift(1.4 * M_SUN, R)
        z3 = pm_surface_redshift(2.0 * M_SUN, R)
        assert z1 < z2 < z3

    def test_decreases_with_radius(self):
        """Larger radius at fixed mass → less compactness → smaller redshift."""
        M = 1.4 * M_SUN
        z1 = pm_surface_redshift(M, 8e3)
        z2 = pm_surface_redshift(M, 10e3)
        z3 = pm_surface_redshift(M, 14e3)
        assert z1 > z2 > z3


# ── GR formula tests ─────────────────────────────────────────────────────────
class TestGRSurfaceRedshiftFormula:
    def test_positive(self):
        M = 1.4 * M_SUN
        R = 10e3
        assert gr_surface_redshift(M, R) > 0.0

    def test_weak_field_linear(self):
        """For small x = 2GM/c²R: z_GR ≈ x/2 = GM/c²R."""
        M = 0.001 * M_SUN
        R = 100e3
        x = 2.0 * G * M / (c * c * R)
        z = gr_surface_redshift(M, R)
        assert z == pytest.approx(x / 2.0, rel=1e-3)

    def test_formula_value(self):
        M = 1.4 * M_SUN
        R = 10e3
        compactness = 2.0 * G * M / (c * c * R)
        expected = 1.0 / math.sqrt(1.0 - compactness) - 1.0
        assert gr_surface_redshift(M, R) == pytest.approx(expected, rel=1e-12)


# ── PM vs GR comparison ───────────────────────────────────────────────────────
class TestPMvsGRRedshift:
    def test_pm_larger_than_gr(self):
        """PM predicts larger surface redshift than GR for any (M, R).

        This is a genuine, testable prediction: z_PM = e^{2GM/c²R} − 1
        while z_GR = 1/sqrt(1−2GM/c²R) − 1.  Since e^x > 1/sqrt(1−x) for
        x > 0, PM always predicts more redshift for the same compactness.
        """
        test_cases = [
            (1.0 * M_SUN, 12e3),
            (1.4 * M_SUN, 10e3),
            (2.0 * M_SUN, 11e3),
        ]
        for M, R in test_cases:
            z_pm = pm_surface_redshift(M, R)
            z_gr = gr_surface_redshift(M, R)
            assert z_pm > z_gr, f"M={M/M_SUN:.1f} M_sun R={R/1e3:.0f} km: z_PM={z_pm:.4f} z_GR={z_gr:.4f}"

    def test_pm_approximately_twice_gr_weak_field(self):
        """In weak field, z_PM ≈ 2 × z_GR for same (M, R)."""
        M = 0.01 * M_SUN   # small compactness
        R = 50e3
        z_pm = pm_surface_redshift(M, R)
        z_gr = gr_surface_redshift(M, R)
        # At leading order: z_PM ≈ 2GM/c²R, z_GR ≈ GM/c²R → ratio ≈ 2
        assert z_pm / z_gr == pytest.approx(2.0, rel=0.01)

    def test_nicer_j0030_pm_vs_gr(self):
        """PSR J0030+0451 NICER central values: PM vs GR redshift.

        Riley et al. (2019): M ≈ 1.34 M_sun, R ≈ 12.71 km.
        PM should give a larger redshift than GR — both are computed here
        for comparison, not as a pass/fail bound (NICER constrains M-R, not z).
        """
        M = 1.34 * M_SUN
        R = 12.71e3
        z_pm = pm_surface_redshift(M, R)
        z_gr = gr_surface_redshift(M, R)
        # Both should be in a physically sensible range
        assert 0.05 < z_pm < 1.0
        assert 0.05 < z_gr < 1.0
        assert z_pm > z_gr

    def test_nicer_j0740_pm_vs_gr(self):
        """PSR J0740+6620 NICER central values: PM vs GR redshift.

        Riley et al. (2021): M ≈ 2.08 M_sun, R ≈ 12.39 km.
        """
        M = 2.08 * M_SUN
        R = 12.39e3
        z_pm = pm_surface_redshift(M, R)
        z_gr = gr_surface_redshift(M, R)
        assert 0.1 < z_pm < 2.0
        assert 0.1 < z_gr < 2.0
        assert z_pm > z_gr


# ── Redshift along M–R track (integration smoke test) ──────────────────────
class TestSurfaceRedshiftTrack:
    @pytest.fixture(scope="class")
    def track(self):
        rho_c, M, R, z = compute_surface_redshift_track(n_points=15)
        return rho_c, M, R, z

    def test_returns_finite_values(self, track):
        _, M, R, z = track
        mask = np.isfinite(M) & np.isfinite(R) & np.isfinite(z)
        assert mask.sum() >= 5, "Too few converged star models"

    def test_redshift_positive_on_track(self, track):
        """Every converged star on the PM M-R track has z > 0."""
        _, M, R, z = track
        mask = np.isfinite(z)
        assert np.all(z[mask] > 0.0)

    def test_redshift_increases_with_compactness(self, track):
        """More compact PM stars (higher M/R) have higher surface redshift."""
        _, M, R, z = track
        mask = np.isfinite(M) & np.isfinite(R) & np.isfinite(z)
        compactness = M[mask] / R[mask]   # proportional to μ_G M/R
        z_sel = z[mask]
        # Check global trend: correlation should be positive
        corr = np.corrcoef(compactness, z_sel)[0, 1]
        assert corr > 0.8, f"z_surf not well correlated with compactness: r={corr:.3f}"
