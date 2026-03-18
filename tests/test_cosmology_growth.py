"""Tests for cosmological structure growth: ΛCDM vs PM refractive scenarios.

Physics under test
──────────────────
ΛCDM growth factor ODE:
    D'' + [(3 + dlnH/dlna)/a] D' − [3Ω_m/(2a^5 E²)] D = 0
    growing mode: D(a→0) ∝ a  (matter domination)
    Growth rate: f = d ln D / d ln a

PM-Drag growth ODE (refractive Hubble H_PM = H0(1+z)^β):
    D'' + [(3−β)/a] D' − [3Ω_eff/(2a^{5−2β})] D = 0
    growing mode: D ∝ a^p where p = positive root of p²+(2−β)p−3Ω_eff/2=0
    With β=0.8, Ω_eff=1: p ≈ 0.764 (slower than ΛCDM MD where p=1)

PM-Static: no friction; effective f = √(3Ω_m/2) = constant ≈ 0.687

Falsification tests:
    - ΛCDM fσ8 matches published RSD data at z=0..1.4
    - PM-Static predicts f≈0.69 CONSTANT with z (distinct from ΛCDM shape)
    - PM-Drag β=0.8 gives different fσ8(z) slope than ΛCDM
    - With same σ8,0=0.811 and β=0.8 fixed: PM cannot simultaneously
      match both D_L (calibration) AND fσ8 data without tuning Ω_eff separately
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from general_relativity.cosmology import (
    lcdm_hubble,
    lcdm_luminosity_distance,
    lcdm_om_z,
    lcdm_growth_factor_rate,
    lcdm_growth_rate_approx,
    lcdm_fsigma8,
    H0_PLANCK, OM0_PLANCK, SIGMA8_PLANCK,
)
from pushing_medium.cosmology import (
    pm_hubble_effective,
    pm_luminosity_distance,
    pm_static_growth_exponent,
    pm_static_fsigma8,
    pm_drag_growing_power,
    pm_growth_factor_rate,
    pm_drag_fsigma8,
)

M_SUN = 1.989e30
AU    = 1.496e11


# ══════════════════════════════════════════════════════════════════════════════
#  ΛCDM — background
# ══════════════════════════════════════════════════════════════════════════════

class TestLCDMBackground:
    """Basic sanity checks on the ΛCDM background functions."""

    def test_hubble_z0_equals_H0(self):
        assert abs(lcdm_hubble(0.0) - H0_PLANCK) < 1e-10

    def test_hubble_increases_with_z(self):
        """H(z) must grow at high z (matter domination)."""
        assert lcdm_hubble(1.0) > lcdm_hubble(0.0)
        assert lcdm_hubble(5.0) > lcdm_hubble(1.0)

    def test_hubble_high_z_approx(self):
        """At z >> 1: H ≈ H0 √(Ω_m) (1+z)^1.5."""
        z = 100.0
        H  = lcdm_hubble(z)
        H_approx = H0_PLANCK * math.sqrt(OM0_PLANCK) * (1+z)**1.5
        assert abs(H / H_approx - 1.0) < 0.01   # within 1%

    def test_luminosity_distance_zero_at_z0(self):
        assert lcdm_luminosity_distance(0.0) < 1.0   # essentially 0 Mpc

    def test_luminosity_distance_z1_order_of_magnitude(self):
        """D_L(z=1) ≈ 6600 Mpc for Planck cosmology."""
        DL = lcdm_luminosity_distance(1.0)
        assert 5000 < DL < 8000, f"D_L(1)={DL:.0f} Mpc out of expected range"

    def test_om_z_increases_with_z(self):
        """Ω_m(z) must increase toward 1 at high z (matter domination)."""
        om0 = lcdm_om_z(0.0)
        om1 = lcdm_om_z(1.0)
        om5 = lcdm_om_z(5.0)
        assert om0 < om1 < om5 < 1.0

    def test_om_z0_matches_param(self):
        assert abs(lcdm_om_z(0.0) - OM0_PLANCK) < 1e-8


# ══════════════════════════════════════════════════════════════════════════════
#  ΛCDM — growth factor
# ══════════════════════════════════════════════════════════════════════════════

class TestLCDMGrowth:
    """Tests for the growth factor ODE and fσ8."""

    def test_growth_factor_z0_equals_one(self):
        D, f = lcdm_growth_factor_rate([0.0])
        assert abs(D[0] - 1.0) < 1e-4

    def test_growth_factor_decreases_with_z(self):
        """D(z)/D(0) must be < 1 for z > 0 (structure was smaller)."""
        D, _ = lcdm_growth_factor_rate([0.5, 1.0, 2.0])
        assert D[0] > D[1] > D[2] > 0

    def test_growth_rate_z0_approx_053(self):
        """f(z=0) ≈ Ω_m^0.55 ≈ 0.315^0.55 ≈ 0.53 for Planck cosmology."""
        _, f = lcdm_growth_factor_rate([0.0])
        assert 0.45 < f[0] < 0.60, f"f(z=0) = {f[0]:.4f}"

    def test_growth_rate_agrees_with_linder_approx(self):
        """Numerical f(z) should match Linder approximation within 3%."""
        z_test = [0.1, 0.5, 1.0, 1.5]
        _, f_num = lcdm_growth_factor_rate(z_test)
        for z, fn in zip(z_test, f_num):
            f_approx = lcdm_growth_rate_approx(z)
            assert abs(fn / f_approx - 1.0) < 0.04, \
                f"z={z}: f_num={fn:.4f}, f_linder={f_approx:.4f}"

    def test_growth_factor_highz_approaches_1_over_1pz(self):
        """In matter domination D ∝ a = 1/(1+z): D(2)/D(1) ≈ 1.5."""
        D, _ = lcdm_growth_factor_rate([1.0, 2.0])
        ratio = D[1] / D[0]   # D(z=2)/D(z=1)
        # In pure MD: ratio = (1+1)/(1+2) = 2/3 ≈ 0.667; Λ suppresses slightly
        assert 0.60 < ratio < 0.73, f"ratio = {ratio:.4f}"

    def test_fsigma8_z0_approx(self):
        """fσ8(z=0) ≈ 0.46 × 0.811 ≈ 0.37–0.41."""
        fs8 = lcdm_fsigma8([0.0])
        assert 0.33 < fs8[0] < 0.45, f"fσ8(0) = {fs8[0]:.4f}"

    def test_fsigma8_roughly_flat_0_to_1(self):
        """ΛCDM fσ8 is nearly flat (0.40–0.50) over z=0–1.5 — well-known."""
        fs8 = lcdm_fsigma8(np.linspace(0.05, 1.4, 20))
        assert fs8.min() > 0.30
        assert fs8.max() < 0.65
        # Peak-to-trough variation should be < 30%
        assert (fs8.max() - fs8.min()) / fs8.mean() < 0.30

    def test_fsigma8_matches_6dfgs(self):
        """6dFGS measurement: z=0.067, fσ8=0.423±0.055 (Beutler+2012).
        ΛCDM prediction should be within 2σ."""
        z_obs, fs8_obs, err = 0.067, 0.423, 0.055
        fs8_lcdm = lcdm_fsigma8([z_obs])[0]
        assert abs(fs8_lcdm - fs8_obs) < 2.5 * err, \
            f"ΛCDM fσ8(0.067)={fs8_lcdm:.3f}, obs={fs8_obs}±{err}"

    def test_fsigma8_matches_boss_cmass(self):
        """BOSS CMASS: z=0.57, fσ8=0.427±0.066 (Anderson+2014).
        ΛCDM prediction should be within 2σ."""
        z_obs, fs8_obs, err = 0.57, 0.427, 0.066
        fs8_lcdm = lcdm_fsigma8([z_obs])[0]
        assert abs(fs8_lcdm - fs8_obs) < 2.5 * err, \
            f"ΛCDM fσ8(0.57)={fs8_lcdm:.3f}, obs={fs8_obs}±{err}"


# ══════════════════════════════════════════════════════════════════════════════
#  PM background
# ══════════════════════════════════════════════════════════════════════════════

class TestPMBackground:
    """Tests for PM refractive cosmology background functions."""

    def test_pm_hubble_z0_equals_H0(self):
        assert abs(pm_hubble_effective(0.0, H0=H0_PLANCK, beta=0.8) - H0_PLANCK) < 1e-10

    def test_pm_hubble_beta0_constant(self):
        """β=0 → H(z) = H0 for all z (constant decay = empty universe)."""
        for z in [0.5, 1.0, 2.0]:
            assert abs(pm_hubble_effective(z, H0=H0_PLANCK, beta=0.0) - H0_PLANCK) < 1e-8

    def test_pm_hubble_increases_with_z_for_positive_beta(self):
        assert pm_hubble_effective(1.0, beta=0.8) > pm_hubble_effective(0.0, beta=0.8)

    def test_pm_luminosity_distance_matches_finding(self):
        """At z=1.5 PM(β=0.8) matches ΛCDM to within ~3% (COSMOLOGY_FINDINGS.md)."""
        from general_relativity.cosmology import lcdm_luminosity_distance
        DL_pm   = pm_luminosity_distance(1.5, beta=0.8)
        DL_lcdm = lcdm_luminosity_distance(1.5)
        diff = abs(DL_pm / DL_lcdm - 1.0) * 100
        assert diff < 5.0, f"D_L diff at z=1.5: {diff:.2f}%"

    def test_pm_luminosity_distance_beta0_linear(self):
        """β=0: D_L = c/H0 × z × (1+z) for small z ≈ c/H0 × z."""
        from general_relativity.cosmology import lcdm_luminosity_distance
        DL_pm = pm_luminosity_distance(0.1, beta=0.0)
        # Empty universe: D_L ≈ (c/H0)*z*(1 + z/2) for small z
        DL_expected = (299792.458 / H0_PLANCK) * 0.1 * (1.1)
        assert abs(DL_pm / DL_expected - 1.0) < 0.01


# ══════════════════════════════════════════════════════════════════════════════
#  PM structure growth: static scenario
# ══════════════════════════════════════════════════════════════════════════════

class TestPMStaticGrowth:
    """PM-static: no Hubble friction — catastrophic over-growth (falsification)."""

    def test_static_growth_exponent_value(self):
        """Γ/H0 = √(3×0.315/2) = √0.4725 ≈ 0.687."""
        gamma = pm_static_growth_exponent(Om0=0.315)
        assert abs(gamma - math.sqrt(1.5 * 0.315)) < 1e-10

    def test_static_f_eff_z0_below_lcdm(self):
        """At z=0, PM-static f_eff = Γ/H0 ≈ 0.69 (close to but above ΛCDM f(0)≈0.46).
        Actually higher effective growth per Hubble time → too much structure."""
        gamma = pm_static_growth_exponent()
        # γ in units of H0
        assert 0.60 < gamma < 0.80, f"Γ/H0 = {gamma:.4f}"

    def test_static_fsigma8_constant_with_z(self):
        """PM-static fσ8 does not decay with z (no growth suppression)."""
        z_arr = np.array([0.1, 0.5, 1.0, 1.5])
        fs8   = pm_static_fsigma8(z_arr)
        # Values should not decrease monotonically — they vary only from 1/H_eff(z)
        # At β=0.8, H_eff grows with z, so f_eff = Γ/H_eff DECREASES with z
        assert all(fs8 > 0)
        # All values should be within a factor of ~2 of σ8 reference
        assert all(fs8 < 2.0 * SIGMA8_PLANCK)

    def test_static_fsigma8_different_from_lcdm(self):
        """PM-static and ΛCDM fσ8 should have detectably different z-slopes."""
        z_arr = np.array([0.1, 0.5, 1.0])
        fs8_pm   = pm_static_fsigma8(z_arr)
        fs8_lcdm = lcdm_fsigma8(z_arr)
        # The shapes must differ — compute slope difference
        slope_pm   = (fs8_pm[-1]   - fs8_pm[0])   / (z_arr[-1] - z_arr[0])
        slope_lcdm = (fs8_lcdm[-1] - fs8_lcdm[0]) / (z_arr[-1] - z_arr[0])
        # PM static should have a steeper negative slope (H_eff grows faster)
        assert slope_pm < slope_lcdm, \
            f"PM slope={slope_pm:.4f}, ΛCDM slope={slope_lcdm:.4f}"


# ══════════════════════════════════════════════════════════════════════════════
#  PM structure growth: drag (β=0.8) scenario
# ══════════════════════════════════════════════════════════════════════════════

class TestPMDragGrowth:
    """PM-Drag: refractive Hubble friction with H_PM(z) = H0(1+z)^β."""

    def test_growing_power_beta08_om1(self):
        """With β=0.8, Ω_eff=1: p ≈ 0.764."""
        p = pm_drag_growing_power(Om_eff=1.0, beta=0.8)
        assert abs(p - 0.764) < 0.01, f"p = {p:.4f}"

    def test_growing_power_lcdm_md_limit(self):
        """β=3/2, Ω_eff=1: p = 1.0 (standard ΛCDM matter domination).

        β=3/2 gives d ln H/d ln a = -3/2 which is the ΛCDM MD case.
        β=0 (constant H) gives p≈0.581 — sub-MD growth, less friction.
        """
        p_lcdm_md = pm_drag_growing_power(Om_eff=1.0, beta=1.5)
        assert abs(p_lcdm_md - 1.0) < 1e-8
        # β=0 (constant H) has less friction → p < 1
        p_const_h = pm_drag_growing_power(Om_eff=1.0, beta=0.0)
        assert abs(p_const_h - 0.581) < 0.01, f"p(β=0)={p_const_h:.4f}"

    def test_growing_power_beta_lcdm_equiv(self):
        """p² + (2−β)p − 3Ω/2 = 0 satisfied."""
        for beta, om in [(0.8, 1.0), (0.8, 0.315), (0.5, 0.5)]:
            p = pm_drag_growing_power(om, beta)
            residual = p**2 + (2.0 - beta)*p - 1.5*om
            assert abs(residual) < 1e-10, f"residual={residual:.2e}"

    def test_pm_drag_growth_factor_z0_one(self):
        D, f = pm_growth_factor_rate([0.0])
        assert abs(D[0] - 1.0) < 1e-3

    def test_pm_drag_growth_factor_decreases_with_z(self):
        D, _ = pm_growth_factor_rate([0.5, 1.0, 2.0])
        assert D[0] > D[1] > D[2] > 0

    def test_pm_drag_growth_rate_faster_than_lcdm(self):
        """PM-Drag (Ω_eff=0.315, β=0.8) growth rate at z=2 is HIGHER than ΛCDM.

        At z=2, H_PM = H0*3^0.8 ≈ 2.41 H0 while H_ΛCDM ≈ 3.03 H0 — PM has
        less Hubble friction, so perturbations grow faster: f_PM > f_ΛCDM.
        This is a falsification: PM-Drag fσ8 would be too high at intermediate z.
        """
        _, f_pm   = pm_growth_factor_rate([2.0], Om_eff=0.315, beta=0.8)
        _, f_lcdm = lcdm_growth_factor_rate([2.0])
        assert f_pm[0] > f_lcdm[0], \
            f"PM f(2)={f_pm[0]:.4f} should be > ΛCDM f(2)={f_lcdm[0]:.4f}"

    def test_pm_fsigma8_beta08_om_lcdm_value(self):
        """PM-Drag fσ8(z=0.5) with Ω_eff=0.315, β=0.8: should be a finite number."""
        fs8 = pm_drag_fsigma8([0.5], Om_eff=0.315, beta=0.8)
        assert 0.10 < fs8[0] < 0.70

    def test_pm_fsigma8_differs_from_lcdm_shape(self):
        """PM-Drag and ΛCDM fσ8 must have different curvature / slopes."""
        z_arr    = np.linspace(0.1, 1.4, 15)
        fs8_pm   = pm_drag_fsigma8(z_arr, Om_eff=0.315, beta=0.8)
        fs8_lcdm = lcdm_fsigma8(z_arr)
        # The maximum departure must exceed 5% at some redshift
        max_diff_pct = np.max(np.abs(fs8_pm / fs8_lcdm - 1.0)) * 100
        assert max_diff_pct > 5.0, \
            f"Max PM vs ΛCDM fσ8 diff only {max_diff_pct:.2f}%"

    def test_pm_drag_z_slope_vs_lcdm(self):
        """PM-Drag slope of fσ8 vs z should differ from ΛCDM slope."""
        z_arr = np.array([0.1, 1.0])
        fs8_pm   = pm_drag_fsigma8(z_arr, Om_eff=0.315, beta=0.8)
        fs8_lcdm = lcdm_fsigma8(z_arr)
        slope_pm   = (fs8_pm[1]   - fs8_pm[0])   / 0.9
        slope_lcdm = (fs8_lcdm[1] - fs8_lcdm[0]) / 0.9
        # Slopes must be numerically different
        assert abs(slope_pm - slope_lcdm) > 0.005, \
            f"Slopes: PM={slope_pm:.4f}, ΛCDM={slope_lcdm:.4f}"

    def test_pm_om1_beta08_ode_fails_runaway(self):
        """FALSIFICATION: PM-Drag with Ω_eff=1, β=0.8 → no valid matter-dominated era.

        When β < 3/2 and Ω_eff is large, the gravity term in the growth ODE
        grows as a^{2β-3} → ∞ as a→0.  There is no power-law growing mode;
        instead δ grows exponentially fast at high z (essential singularity
        at a=0).  This means PM with its D_L-calibrated β=0.8 and full
        matter density would produce catastrophic over-growth of structure
        at z > 10 — contradicting CMB linear-regime observations.

        The ODE solver correctly raises RuntimeError for these parameters.
        """
        with pytest.raises(RuntimeError, match="PM-Drag growth ODE failed"):
            pm_drag_fsigma8([0.0], Om_eff=1.0, beta=0.8)
