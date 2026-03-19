"""Tests for the PM photon-sphere guarantee across the full M-R curve.

The PM optical metric has a photon sphere (circular light orbit) at:
    r_ps = 2 GM/c²   (= the Schwarzschild radius, 2/3 of GR's 3GM/c²)

The question is: where does this sit relative to the stellar surface?

The answer depends on the compactness  C = GM/(c²R):
  • C < 0.5  →  R_star > r_ps  →  photon sphere INSIDE the star
                               →  no exterior trapped photon orbits
  • C > 0.5  →  R_star < r_ps  →  photon sphere EXTERIOR to the star
                               →  exterior photon ring exists
                               →  light-ring features possible in GW ringdown

PM-specific result (from integrating the stellar structure equations):
──────────────────────────────────────────────────────────────────────
  Low-mass regime  (M ≲ 8 M☉, typical NS):
    C ≈ 0.10–0.49 ; photon sphere safely interior.

  High-mass regime (M ≳ 8 M☉, near M_max ≈ 13.4 M☉):
    C > 0.5 ; photon sphere exterior.
    At maximum mass: C ≈ 0.74, r_ps ≈ 1.49 × R_star.

  All models: C < 1.0 (no GR-type horizon; PM has a flat Minkowski background).

Why PM can exceed C = 0.5 (unlike GR):
  GR imposes the Buchdahl limit: C < 4/9 ≈ 0.44 for uniform stars.
  PM uses a flat Minkowski background — there is no background curvature,
  hence no Buchdahl limit and no event horizon.  A PM star with C > 0.5
  is simply an
  ultracompact horizonless object; its photon ring can be observed
  in principle by photon-orbit-sensitive instruments (EHT, future).

LIGO / EM constraint status:
  No existing observation constrains C > 0.5 for PM-like stars.
  The EHT M87* shadow measures the shadow size (∝ M/D), not C.
  This is a structural prediction, not a violation.
"""

import math
import sys
import pytest
import numpy as np

sys.path.insert(0, 'src')

from pushing_medium import pm_photon_sphere_radius, G, c
from pushing_medium.stellar_structure import (
    solve_pm_star, compute_mr_curve, RHO_NUC, RHO_CRIT, M_SUN,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def compactness(M_kg: float, R_m: float) -> float:
    """C = GM/(c²R)."""
    return G * M_kg / (c * c * R_m)


def photon_sphere_outside(M_kg: float, R_m: float) -> bool:
    """True when the PM photon sphere lies outside the stellar surface."""
    r_ps = pm_photon_sphere_radius(M_kg)
    return r_ps > R_m


# Typical astrophysical parameters
M_14 = 1.4 * M_SUN        # 1.4 M_sun in kg
M_14_RHO_C = 2.688e17     # rho_c [kg/m³] giving M ≈ 1.4 M_sun (from scan)


# ===========================================================================
class TestPhotonSphereFormulas:
    """Basic formula checks — pm_photon_sphere_radius == 2GM/c²."""

    def test_formula(self):
        """r_ps = 2GM/c² exactly."""
        M = 1.0 * M_SUN
        assert pm_photon_sphere_radius(M) == pytest.approx(2.0 * G * M / c**2, rel=1e-12)

    def test_pm_smaller_than_gr(self):
        """PM r_ps = 2GM/c² < GR r_ps = 3GM/c²."""
        from general_relativity import gr_photon_sphere_radius
        M = 5.0 * M_SUN
        r_pm = pm_photon_sphere_radius(M)
        r_gr = gr_photon_sphere_radius(M)
        assert r_pm == pytest.approx(r_gr * 2.0 / 3.0, rel=1e-12)

    def test_photon_sphere_scales_linearly(self):
        """r_ps ∝ M: doubling M doubles r_ps."""
        r1 = pm_photon_sphere_radius(1.0 * M_SUN)
        r2 = pm_photon_sphere_radius(2.0 * M_SUN)
        assert r2 == pytest.approx(2.0 * r1, rel=1e-12)

    def test_numerical_1msun(self):
        """r_ps(1 M☉) ≈ 2953 m ≈ 3.0 km."""
        r = pm_photon_sphere_radius(M_SUN)
        assert 2900.0 < r < 3100.0


# ===========================================================================
class TestLowMassPhotonSphereInside:
    """For low-mass PM stars (< 3 M☉), the photon sphere is well inside."""

    @pytest.fixture(scope="class")
    def star_14msun(self):
        """Integrate the 1.4 M_sun PM star."""
        # Use a central density near the known 1.4 M_sun solution
        rho_c = 2.688e17   # ≈ 1.164 × ρ_nuc
        return solve_pm_star(rho_c, r_max=3e4, n_eval=2000)

    def test_1p4msun_photon_sphere_inside(self, star_14msun):
        """1.4 M☉ star: R_star >> r_ps (photon sphere deep inside)."""
        star = star_14msun
        assert star.converged
        r_ps = pm_photon_sphere_radius(star.M_star)
        assert star.R_star > r_ps, (
            f"R_star={star.R_star/1e3:.2f} km, r_ps={r_ps/1e3:.2f} km"
        )

    def test_1p4msun_compactness_low(self, star_14msun):
        """1.4 M☉: C ≈ 0.148, well below the C=0.5 exterior-photon-sphere threshold."""
        star = star_14msun
        C = compactness(star.M_star, star.R_star)
        assert C < 0.50, f"C = {C:.4f}"
        assert C > 0.05, f"Unexpectedly low C = {C:.4f}"

    def test_1p4msun_radius_to_photon_sphere_ratio(self, star_14msun):
        """R_star / r_ps > 1 for the 1.4 M☉ star (margin at least 2×)."""
        star = star_14msun
        r_ps = pm_photon_sphere_radius(star.M_star)
        ratio = star.R_star / r_ps
        assert ratio > 2.0, f"R/r_ps = {ratio:.2f}; expected > 2"

    def test_low_mass_sector_all_inside(self):
        """All models with M < 3 M☉ have photon sphere inside: max C < 0.5."""
        rho_c, M, R = compute_mr_curve(n_points=30, rho_max_factor=0.999)
        mask = np.isfinite(M) & (M < 3.0)
        if not mask.any():
            pytest.skip("No M < 3 M_sun models in sample")
        C = G * M[mask] * M_SUN / (c**2 * R[mask] * 1e3)
        assert C.max() < 0.50, (
            f"Found M < 3 M_sun model with C = {C.max():.4f} (expected < 0.5)"
        )


# ===========================================================================
class TestHighMassPhotonSphereExterior:
    """Near maximum mass, PM stars are ultracompact: C > 0.5, r_ps exterior."""

    @pytest.fixture(scope="class")
    def mr_curve(self):
        rho_c, M, R = compute_mr_curve(n_points=50, rho_max_factor=0.998)
        mask = np.isfinite(M) & np.isfinite(R)
        return M[mask], R[mask]

    def test_max_mass_compactness_above_half(self, mr_curve):
        """Maximum-mass PM star has C > 0.5 (photon sphere exterior)."""
        M, R = mr_curve
        idx_max = M.argmax()
        C = compactness(M[idx_max] * M_SUN, R[idx_max] * 1e3)
        assert C > 0.50, f"C_max = {C:.4f}; expected > 0.5"

    def test_max_mass_compactness_below_unity(self, mr_curve):
        """C < 1.0 for all models — no event-horizon-scale compactness."""
        M, R = mr_curve
        C = G * M * M_SUN / (c**2 * R * 1e3)
        assert C.max() < 1.0, f"Max C = {C.max():.4f}; expected < 1.0"

    def test_max_mass_photon_sphere_outside(self, mr_curve):
        """At maximum mass, r_ps > R_star (photon sphere exterior)."""
        M, R = mr_curve
        idx_max = M.argmax()
        M_kg = M[idx_max] * M_SUN
        R_m  = R[idx_max] * 1e3
        r_ps = pm_photon_sphere_radius(M_kg)
        assert r_ps > R_m, (
            f"r_ps = {r_ps/1e3:.2f} km, R_star = {R_m/1e3:.2f} km; "
            f"expected photon sphere exterior for maximum-mass star"
        )

    def test_max_mass_r_ps_to_R_ratio(self, mr_curve):
        """At max mass, r_ps / R_star ≈ 1.4–1.6 (ultracompact regime)."""
        M, R = mr_curve
        idx_max = M.argmax()
        M_kg = M[idx_max] * M_SUN
        R_m  = R[idx_max] * 1e3
        ratio = pm_photon_sphere_radius(M_kg) / R_m
        assert 1.3 < ratio < 1.7, (
            f"r_ps/R = {ratio:.3f}; expected 1.3–1.7 for ultracompact PM star"
        )

    def test_maximum_mass_value(self, mr_curve):
        """PM maximum mass ≈ 13–14 M☉ (baseline EOS)."""
        M, R = mr_curve
        assert 11.0 < M.max() < 16.0, f"M_max = {M.max():.2f} M_sun"


# ===========================================================================
class TestPhotonSphereTransition:
    """The C = 0.5 transition occurs at an intermediate mass on the M-R curve."""

    @pytest.fixture(scope="class")
    def mr_curve(self):
        rho_c, M, R = compute_mr_curve(n_points=60, rho_max_factor=0.999)
        mask = np.isfinite(M) & np.isfinite(R)
        return M[mask], R[mask]

    def test_transition_mass_in_range(self, mr_curve):
        """The C = 0.5 transition (R = r_ps) occurs at M ≈ 5–12 M☉."""
        M, R = mr_curve
        C = G * M * M_SUN / (c**2 * R * 1e3)
        # There should be models both below and above C = 0.5
        assert (C < 0.5).any(), "No models with C < 0.5 found"
        assert (C > 0.5).any(), "No models with C > 0.5 found"
        # Transition mass
        trans_idx = np.where(C > 0.5)[0][0]
        M_trans = M[trans_idx]
        assert 5.0 < M_trans < 12.0, (
            f"Transition mass = {M_trans:.2f} M_sun; expected 5–12 M_sun"
        )

    def test_compactness_monotonically_increases_with_mass(self, mr_curve):
        """C increases monotonically along the M-R curve as M increases."""
        M, R = mr_curve
        C = G * M * M_SUN / (c**2 * R * 1e3)
        # Sort by M and check C is non-decreasing (with small tolerance for numerical noise)
        idx_sort = np.argsort(M)
        C_sorted = C[idx_sort]
        diff = np.diff(C_sorted)
        # Allow tiny numerical dips (< 0.01 in C); overall trend must increase
        n_violations = (diff < -0.01).sum()
        assert n_violations == 0, (
            f"Compactness is non-monotonic in {n_violations} places"
        )

    def test_three_solar_mass_photon_sphere_inside(self, mr_curve):
        """3 M☉ representative: photon sphere still inside the star."""
        M, R = mr_curve
        # Find closest model to 3 M_sun
        idx = np.argmin(np.abs(M - 3.0))
        M_kg = M[idx] * M_SUN
        R_m  = R[idx] * 1e3
        r_ps = pm_photon_sphere_radius(M_kg)
        # At 3 M_sun: C ≈ 0.24 → R_star >> r_ps
        assert R_m > r_ps, (
            f"3 M_sun model: R={R_m/1e3:.2f} km, r_ps={r_ps/1e3:.2f} km"
        )


# ===========================================================================
class TestNoBuchdahlLimit:
    """PM has no Buchdahl limit — stars can exceed C = 4/9 (GR Buchdahl)."""

    def test_maximum_compactness_exceeds_buchdahl(self):
        """PM max C > 4/9 ≈ 0.444 (Buchdahl limit for uniform-density GR stars)."""
        BUCHDAHL = 4.0 / 9.0   # ≈ 0.444
        rho_c, M, R = compute_mr_curve(n_points=50, rho_max_factor=0.998)
        mask = np.isfinite(M) & np.isfinite(R)
        C = G * M[mask] * M_SUN / (c**2 * R[mask] * 1e3)
        assert C.max() > BUCHDAHL, (
            f"Max C = {C.max():.4f}; expected > Buchdahl limit {BUCHDAHL:.4f}"
        )

    def test_maximum_compactness_no_horizon(self):
        """C < 0.5 is the GR horizon threshold (for comparison only).
        PM CAN exceed this — no event horizon forms; the background is flat Minkowski.
        """
        rho_c, M, R = compute_mr_curve(n_points=50, rho_max_factor=0.998)
        mask = np.isfinite(M) & np.isfinite(R)
        C = G * M[mask] * M_SUN / (c**2 * R[mask] * 1e3)
        # PM allows C > 0.5 — this is a feature, not a bug.
        # All models must still be sub-luminal (C < 1.0):
        assert C.max() < 1.0, f"Max C = {C.max():.4f}; PM should not reach C = 1"
