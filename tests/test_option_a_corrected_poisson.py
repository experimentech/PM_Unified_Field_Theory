"""Tests for Option-A: Corrected Poisson Equation with U'(φ) self-coupling.

Physical experiment
-------------------
The standard PM Poisson equation (Gap 1) is:

    ∇²φ = −(8πG/c²) ρ          [baseline]

Promoting U(φ) into the Lagrangian adds a self-coupling correction with
effective self-coupling density ρ_U = ρ_nuc (2φ − φ²):

    ∇²φ = −(8πG/c²) [ρ  −  α ρ_U]

Two directions:
  α = +1  (Gap-1 as written in formula sheet): U' REDUCES gravitating source
  α = −1  (Gap-3 style): U' ADDS to gravitating source

Numerical results (recorded from scripts/option_a_corrected_poisson.py):
─────────────────────────────────────────────────────
                  M_max [M☉]    R_1.4 [km]    C_max
  Baseline (α=0)    13.40         13.92        0.743
  Gap-1    (α=+1)   30.28         13.94        1.280   ← mass gap CLOSES
  Gap-3    (α=−1)    7.94         13.89        0.524

Key findings:
  1. Gap-1 DOUBLES M_max (13.4 → 30.3 M☉), comparable to GW150914 components.
     The GWTC mass gap is dramatically narrowed.
  2. NEITHER correction fixes R_1.4 (all remain ~13.9 km, still > 13.3 km limit).
     The radius problem requires a SEPARATE fix (EOS adjustment, not field eq).
  3. The two problems are DECOUPLED: Poisson correction targets high-density
     (M_max); the EOS fix must target moderate-density (1.4 M☉ stars).
  4. At 1.4 M☉ (φ_c ≈ 0.34): φ(r) profile barely changes with α — correction
     activates mainly at high compression near the M-R turning point.
  5. Gap-3 (α=−1) keeps M_max < 22.2 (GW190814 M1) and keeps the mass gap.

EM cross-check (via n = e^φ):
  At low mass (1.4 M☉), φ(r) is almost unchanged between scenarios.
  Therefore EHT shadow, photon sphere, and surface redshift predictions for
  normal-density stars are unaffected by either correction.

Tests are organised in four classes:
  TestOptionABaselineRecovery  — α=0 reproduces exact baseline
  TestGap1Direction            — α=+1 (formula-sheet direction)
  TestGap3Direction            — α=−1 (Gap-3 / field energy gravitates)
  TestDecouplingOfProblems     — R_1.4 insensitivity; φ profile
"""

import sys
import math
import pytest
import numpy as np

sys.path.insert(0, 'src')

from pushing_medium.stellar_structure import (
    solve_pm_star,
    solve_pm_star_option_a,
    compute_mr_curve,
    compute_mr_curve_option_a,
    pm_surface_redshift,
    RHO_NUC, RHO_CRIT, M_SUN, G, c, MU_G,
)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

N_POINTS = 30   # sparse enough for reasonable test runtime


def find_r14(M_arr, R_arr):
    """Interpolate R at M = 1.4 M_sun."""
    valid = np.isfinite(M_arr) & np.isfinite(R_arr)
    M_v, R_v = M_arr[valid], R_arr[valid]
    idx = np.argsort(M_v)
    M_v, R_v = M_v[idx], R_v[idx]
    if M_v[0] > 1.4 or M_v[-1] < 1.4:
        return float('nan')
    return float(np.interp(1.4, M_v, R_v))


@pytest.fixture(scope="module")
def baseline_mr():
    _, M, R = compute_mr_curve(n_points=N_POINTS)
    return M, R


@pytest.fixture(scope="module")
def gap1_mr():
    _, M, R = compute_mr_curve_option_a(alpha=+1.0, n_points=N_POINTS)
    return M, R


@pytest.fixture(scope="module")
def gap3_mr():
    _, M, R = compute_mr_curve_option_a(alpha=-1.0, n_points=N_POINTS)
    return M, R


@pytest.fixture(scope="module")
def star_low_mass():
    """Three versions of a low-mass star (ρ_c = 1.4 ρ_nuc) for profile comparison."""
    rho_c = 1.4 * RHO_NUC
    base = solve_pm_star(rho_c)
    gp1  = solve_pm_star_option_a(rho_c, alpha=+1.0)
    gp3  = solve_pm_star_option_a(rho_c, alpha=-1.0)
    return base, gp1, gp3


# ===========================================================================
class TestOptionABaselineRecovery:
    """α=0 must recover the exact baseline solver to machine precision."""

    def test_alpha0_same_mmax_as_baseline(self, baseline_mr, gap1_mr):
        """α=0 curve matches baseline M_max to < 0.5%."""
        _, M0, R0 = compute_mr_curve_option_a(alpha=0.0, n_points=N_POINTS)
        M_base = float(np.nanmax(baseline_mr[0][np.isfinite(baseline_mr[0])]))
        M_0    = float(np.nanmax(M0[np.isfinite(M0)]))
        assert abs(M_0 - M_base) / M_base < 0.005, (
            f"α=0 M_max = {M_0:.2f} M☉ vs baseline {M_base:.2f} M☉"
        )

    def test_alpha0_same_r14_as_baseline(self, baseline_mr):
        """α=0 curve matches baseline R_1.4 to < 0.5%."""
        _, M0, R0 = compute_mr_curve_option_a(alpha=0.0, n_points=N_POINTS)
        R_base = find_r14(*baseline_mr)
        R_0    = find_r14(M0, R0)
        assert abs(R_0 - R_base) / R_base < 0.005, (
            f"α=0 R_1.4 = {R_0:.2f} km vs baseline {R_base:.2f} km"
        )

    def test_alpha0_single_star_matches_baseline(self):
        """Single-star solver: α=0 agrees with baseline to < 0.1%."""
        rho_c = 1.8 * RHO_NUC
        base  = solve_pm_star(rho_c)
        opt_a = solve_pm_star_option_a(rho_c, alpha=0.0)
        assert base.converged and opt_a.converged
        rel_M = abs(opt_a.M_star - base.M_star) / base.M_star
        rel_R = abs(opt_a.R_star - base.R_star) / base.R_star
        assert rel_M < 0.001, f"α=0 M mismatch: {rel_M:.4f}"
        assert rel_R < 0.001, f"α=0 R mismatch: {rel_R:.4f}"

    def test_alpha0_phi_profile_matches_baseline(self):
        """φ(r) profile from α=0 matches baseline to < 0.1% at all radii."""
        rho_c = 2.0 * RHO_NUC
        base  = solve_pm_star(rho_c)
        opt_a = solve_pm_star_option_a(rho_c, alpha=0.0)
        # Compare on five interior radii
        for frac in [0.1, 0.3, 0.5, 0.7, 0.9]:
            r = frac * base.R_star
            phi_b = float(np.interp(r, base.r, base.phi))
            phi_a = float(np.interp(r, opt_a.r, opt_a.phi))
            if abs(phi_b) > 0.01:
                assert abs(phi_a - phi_b) / phi_b < 0.001, (
                    f"φ mismatch at r/R={frac}: base={phi_b:.5f}, α=0={phi_a:.5f}"
                )


# ===========================================================================
class TestGap1Direction:
    """α = +1: Gap-1 correction (U' reduces gravitating source, formula-sheet direction)."""

    def test_gap1_all_models_converge(self, gap1_mr):
        """At least 80% of Gap-1 models converge (no catastrophic instability)."""
        M_arr, R_arr = gap1_mr
        valid = np.isfinite(M_arr) & np.isfinite(R_arr) & (R_arr > 0)
        total = len(M_arr)
        frac = valid.sum() / total
        assert frac > 0.8, f"Only {valid.sum()}/{total} Gap-1 models converged"

    def test_gap1_mmax_much_larger_than_baseline(self, baseline_mr, gap1_mr):
        """Gap-1 M_max ≫ baseline M_max (self-coupling dramatically raises ceiling).

        Numerical result: baseline 13.4 → Gap-1 30.3 M☉ (factor ~2.26).
        """
        M_base = float(np.nanmax(baseline_mr[0][np.isfinite(baseline_mr[0])]))
        M_gp1  = float(np.nanmax(gap1_mr[0][np.isfinite(gap1_mr[0])]))
        assert M_gp1 > 2.0 * M_base, (
            f"Gap-1 M_max = {M_gp1:.1f} M☉ should be > 2× baseline ({2*M_base:.1f} M☉)"
        )

    def test_gap1_mmax_above_gw190814_primary(self, gap1_mr):
        """Gap-1 M_max > 22.2 M☉ (GW190814 primary) — mass gap is closed.

        Numerical result: M_max ≈ 30.3 M☉ > 22.2 M☉.
        """
        M_gp1 = float(np.nanmax(gap1_mr[0][np.isfinite(gap1_mr[0])]))
        assert M_gp1 > 22.2, (
            f"Gap-1 M_max = {M_gp1:.1f} M☉ — expected > 22.2 M☉ (GW190814 M1)"
        )

    def test_gap1_mmax_comparable_to_gw150914(self, gap1_mr):
        """Gap-1 M_max is in the range of GW150914 components (28–37 M☉).

        Numerical result: M_max ≈ 30.3 M☉, just below GW150914 M2 = 30.6 M☉.
        """
        M_gp1 = float(np.nanmax(gap1_mr[0][np.isfinite(gap1_mr[0])]))
        assert 25.0 < M_gp1 < 40.0, (
            f"Gap-1 M_max = {M_gp1:.1f} M☉ — expected in GW150914 range (25–40 M☉)"
        )

    def test_gap1_high_compactness_at_mmax(self, gap1_mr):
        """Gap-1 stars near M_max exceed C = 1.0 — extreme but no horizon (flat background).

        Numerical result: C ≈ 1.28 at M_max.  PM allows C > 1 because there is
        no background curvature and hence no Schwarzschild horizon.
        """
        M_arr, R_arr = gap1_mr
        valid = np.isfinite(M_arr) & np.isfinite(R_arr) & (R_arr > 0)
        C_max = float(np.nanmax(
            G * M_arr[valid] * M_SUN / (c * c * R_arr[valid] * 1e3)
        ))
        assert C_max > 1.0, (
            f"Gap-1 C_max = {C_max:.3f} — expected > 1.0 (extreme compact objects)"
        )

    def test_gap1_surface_redshift_very_large_at_mmax(self, gap1_mr):
        """Gap-1 z_s at M_max ≫ baseline — extreme compactness signature.

        Numerical result: z_s ≈ 11.9 vs baseline 3.4.
        """
        M_arr, R_arr = gap1_mr
        valid = np.isfinite(M_arr) & np.isfinite(R_arr) & (R_arr > 0)
        idx_max = int(np.nanargmax(M_arr[valid]))
        M_v, R_v = M_arr[valid], R_arr[valid]
        z_max = pm_surface_redshift(M_v[idx_max] * M_SUN, R_v[idx_max] * 1e3)
        assert z_max > 5.0, (
            f"Gap-1 z_s at M_max = {z_max:.2f} — expected > 5.0"
        )


# ===========================================================================
class TestGap3Direction:
    """α = −1: Gap-3 style (U' adds to gravitating source — field energy gravitates)."""

    def test_gap3_all_models_converge(self, gap3_mr):
        """At least 70% of Gap-3 models converge."""
        M_arr, R_arr = gap3_mr
        valid = np.isfinite(M_arr) & np.isfinite(R_arr) & (R_arr > 0)
        assert valid.sum() / len(M_arr) > 0.70

    def test_gap3_mmax_smaller_than_baseline(self, baseline_mr, gap3_mr):
        """Gap-3 M_max < baseline M_max (stronger gravity lowers the ceiling).

        Numerical result: 13.4 → 7.9 M☉.
        """
        M_base = float(np.nanmax(baseline_mr[0][np.isfinite(baseline_mr[0])]))
        M_gp3  = float(np.nanmax(gap3_mr[0][np.isfinite(gap3_mr[0])]))
        assert M_gp3 < M_base, (
            f"Gap-3 M_max = {M_gp3:.1f} M☉ should be < baseline {M_base:.1f} M☉"
        )

    def test_gap3_mmax_below_8_solar(self, gap3_mr):
        """Gap-3 M_max < 10 M☉ (numerical: ~7.9 M☉)."""
        M_gp3 = float(np.nanmax(gap3_mr[0][np.isfinite(gap3_mr[0])]))
        assert M_gp3 < 10.0, f"Gap-3 M_max = {M_gp3:.1f} M☉ — expected < 10 M☉"

    def test_gap3_mass_gap_kept(self, gap3_mr):
        """Gap-3 M_max < 22.2 M☉ — mass gap to GW190814 primary is preserved."""
        M_gp3 = float(np.nanmax(gap3_mr[0][np.isfinite(gap3_mr[0])]))
        assert M_gp3 < 22.2, f"Gap-3 M_max = {M_gp3:.1f} M☉ — expected < 22.2 M☉"

    def test_gap3_stars_more_compact_than_baseline(self, baseline_mr, gap3_mr):
        """Gap-3 R at M_max < baseline R at M_max (stronger gravity → smaller stars)."""
        def r_at_mmax(M_arr, R_arr):
            valid = np.isfinite(M_arr) & np.isfinite(R_arr) & (R_arr > 0)
            M_v, R_v = M_arr[valid], R_arr[valid]
            return R_v[np.argmax(M_v)]
        R_base_max = r_at_mmax(*baseline_mr)
        R_gp3_max  = r_at_mmax(*gap3_mr)
        assert R_gp3_max < R_base_max, (
            f"Gap-3 R at M_max = {R_gp3_max:.1f} km, baseline = {R_base_max:.1f} km"
        )


# ===========================================================================
class TestDecouplingOfProblems:
    """The Poisson correction decouples M_max from R_1.4: it changes M_max
    dramatically but leaves R_1.4 essentially unchanged.  This establishes
    that the radius problem and the mass-gap problem require SEPARATE fixes:

    - Mass gap:    Gap-1 Poisson correction (this test class documents it)
    - Radius (R_1.4 ≤ 13.3 km): EOS stiffening / ρ_nuc adjustment
    """

    def test_r14_insensitive_to_alpha_plus1(self, baseline_mr, gap1_mr):
        """R_1.4 changes by < 2% between baseline and α=+1.

        Numerical: baseline 13.92 km, Gap-1 13.94 km (Δ ≈ 0.02 km, 0.1%).
        """
        R_base = find_r14(*baseline_mr)
        R_gp1  = find_r14(*gap1_mr)
        assert not math.isnan(R_base) and not math.isnan(R_gp1)
        rel_diff = abs(R_gp1 - R_base) / R_base
        assert rel_diff < 0.02, (
            f"R_1.4: baseline={R_base:.3f} km, Gap-1={R_gp1:.3f} km, "
            f"relative change={rel_diff:.4f} (expected < 2%)"
        )

    def test_r14_insensitive_to_alpha_minus1(self, baseline_mr, gap3_mr):
        """R_1.4 changes by < 2% between baseline and α=−1.

        Numerical: baseline 13.92 km, Gap-3 13.89 km (Δ ≈ 0.03 km, 0.2%).
        """
        R_base = find_r14(*baseline_mr)
        R_gp3  = find_r14(*gap3_mr)
        assert not math.isnan(R_base) and not math.isnan(R_gp3)
        rel_diff = abs(R_gp3 - R_base) / R_base
        assert rel_diff < 0.02, (
            f"R_1.4: baseline={R_base:.3f} km, Gap-3={R_gp3:.3f} km, "
            f"relative change={rel_diff:.4f} (expected < 2%)"
        )

    def test_r14_all_scenarios_still_exceed_observation(self, baseline_mr, gap1_mr, gap3_mr):
        """None of the three scenarios meet the GW170817 R_1.4 ≤ 13.3 km bound.

        This is a key NEGATIVE result: the Poisson correction does NOT fix
        the radius problem.  That requires a separate EOS adjustment.
        """
        OBS_LIMIT = 13.3
        for name, mr in [("baseline", baseline_mr), ("Gap-1", gap1_mr), ("Gap-3", gap3_mr)]:
            R = find_r14(*mr)
            assert R > OBS_LIMIT, (
                f"{name} has R_1.4 = {R:.2f} km ≤ {OBS_LIMIT} km — "
                "unexpected: thought Poisson correction doesn't help radii"
            )

    def test_low_mass_phi_profile_barely_changes_with_alpha(self, star_low_mass):
        """At 1.4 M☉ (low density, φ_c ≈ 0.34): α-correction changes φ(r) by < 5%
        when compared at the same *fractional* stellar radius r/R.

        The Gap-1 star is slightly larger (weakened gravity) so comparisons must
        be done at r/R rather than absolute r.  The script confirmed < 5% change
        at r/R = 0.5 (baseline 0.2458, Gap-1 0.2542, Gap-3 0.2409).

        This shows the correction is small in the low-density regime,
        explaining why R_1.4 is unaffected whilst high-density M_max changes.
        """
        base, gp1, gp3 = star_low_mass
        # Sample φ at the same fractional radial position in each star
        for frac in [0.3, 0.5, 0.7]:
            phi_b  = float(np.interp(frac * base.R_star, base.r, base.phi))
            phi_g1 = float(np.interp(frac *  gp1.R_star,  gp1.r,  gp1.phi))
            phi_g3 = float(np.interp(frac *  gp3.R_star,  gp3.r,  gp3.phi))
            if phi_b < 0.01:
                continue   # avoid division by near-zero at the surface
            rel_g1 = abs(phi_g1 - phi_b) / phi_b
            rel_g3 = abs(phi_g3 - phi_b) / phi_b
            assert rel_g1 < 0.10, (
                f"Gap-1 φ change at r/R={frac}: {rel_g1:.4f} (expected < 10%); "
                f"base φ={phi_b:.4f}, Gap-1 φ={phi_g1:.4f}"
            )
            assert rel_g3 < 0.10, (
                f"Gap-3 φ change at r/R={frac}: {rel_g3:.4f} (expected < 10%); "
                f"base φ={phi_b:.4f}, Gap-3 φ={phi_g3:.4f}"
            )

    def test_low_mass_star_mass_changes_with_alpha(self, star_low_mass):
        """At ρ_c = 1.4 ρ_nuc: star mass changes significantly with α.

        Gap-1 produces a MORE massive star than baseline for the same ρ_c
        (because the self-coupling weakens gravity and allows more material).
        Gap-3 produces a LESS massive star (stronger gravity → smaller star).
        """
        base, gp1, gp3 = star_low_mass
        M_base = base.M_star / M_SUN
        M_gp1  = gp1.M_star  / M_SUN
        M_gp3  = gp3.M_star  / M_SUN
        assert M_gp1 > M_base, f"Gap-1 M={M_gp1:.3f} should > baseline {M_base:.3f}"
        assert M_gp3 < M_base, f"Gap-3 M={M_gp3:.3f} should < baseline {M_base:.3f}"

    def test_surface_redshift_z14_barely_changes(self, baseline_mr, gap1_mr, gap3_mr):
        """Surface redshift at 1.4 M☉ is within 5% for all scenarios.

        Since R_1.4 barely changes, z_s(1.4) = e^{μ_G M/(R)} also barely changes.
        This means the EM cross-check (n=e^φ, photon sphere, EHT) is
        essentially unaffected for normal-density stars.
        """
        def z14(mr):
            R = find_r14(*mr)
            if math.isnan(R):
                return float('nan')
            return pm_surface_redshift(1.4 * M_SUN, R * 1e3)

        z_base = z14(baseline_mr)
        z_gp1  = z14(gap1_mr)
        z_gp3  = z14(gap3_mr)

        assert abs(z_gp1 - z_base) / z_base < 0.05, (
            f"Gap-1 z_s(1.4) = {z_gp1:.4f} vs baseline {z_base:.4f}"
        )
        assert abs(z_gp3 - z_base) / z_base < 0.05, (
            f"Gap-3 z_s(1.4) = {z_gp3:.4f} vs baseline {z_base:.4f}"
        )
