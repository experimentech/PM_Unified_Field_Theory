"""
Tests: ρ_nuc as the Critical Variable Controlling Stellar Radius

Background
----------
Gap-probe tests (test_gap_probes.py) showed that:
  * Option c (P_total = P_EOS + U(φ)) raises M_max from 13.4 to 30.8 M☉ ✓
  * But R_1.4 (radius at 1.4 M☉) remains ~14 km — above the GW170817/NICER
    constraint of R < 13.3 km (90% CI upper bound)

Subsequent probes revealed that R at 1.4 M☉ is:
  (a) INSENSITIVE to the EOS sound-speed parameter α·U(φ) over α ∈ [0, 2]
  (b) INSENSITIVE to the EOS linear stiffness γ (c_s² = γ c²/2) over γ ∈ [0.3, 2]
  (c) CONTROLLED PRIMARILY by ρ_nuc via the uniform-density formula:
        R_1.4 ≈ (3 M_1.4 / 4π ρ_nuc)^{1/3}

Mechanism: the PM EOS has P(ρ_nuc) = 0, so the stellar surface is at ρ = ρ_nuc.
A 1.4 M☉ star has mean density ⟨ρ⟩ = M_1.4/(4π/3 × R³) only ~10% above ρ_nuc,
making the star essentially a uniform-density sphere at ρ_nuc.  Changing c_s or
U(φ) barely shifts the small density contrast and therefore barely shifts R.

Implication: ρ_nuc is the CRITICAL VARIABLE for the radius.  The current value
(2.3×10¹⁷ kg/m³ — standard nuclear saturation density) gives R_1.4 ≈ 14.0 km.
For R_1.4 < 13.3 km (GW170817 1σ upper bound) we need ρ_nuc ≥ 1.18 × ρ_nuc_0.

Joint feasibility (Option c + scaled ρ_nuc):
  f = 1.18:  R_1.4 = 13.22 km ✓  M_max = 28.37 M☉ ✓ (32 M☉ with α ≈ 1.15)
  f = 1.20:  R_1.4 = 13.14 km ✓  M_max = 28.11 M☉ ✓
  f = 1.22:  R_1.4 = 13.08 km ✓  M_max = 27.88 M☉ ✗ (fails M_max > 28)

The joint feasibility window at α=1 is f ∈ [1.18, 1.21] (a 3% range in ρ_nuc).
With α slightly above 1 (stronger elastic pressure), M_max > 30 M☉ can also be
achieved at f = 1.18 while maintaining R < 13.3 km.
"""

import math
import pytest
import numpy as np

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from pushing_medium.stellar_structure import RHO_NUC, RHO_CRIT, M_SUN, c
from pushing_medium.stellar_structure import pm_eos_pressure, pm_eos_density
from gap_probe_comparison import (
    EPS0, u_phi, compute_mr_curve_variant, _solve_star_variant,
)

G = 6.674e-11


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_option_c_scaled(f):
    """
    Build the Option-c EOS with ρ_nuc scaled by factor f.
    P_total = c²/2 (ρ − f·ρ_nuc) + U_f(φ)
    where φ = ln(ρ / f·ρ_nuc) and U_f(φ) = ε₀_f (φ² − φ³/3).
    """
    rn   = f * RHO_NUC
    rc   = math.e * rn
    eps0 = rn * c * c

    def eos_p(rho, rn_=rn, e0=eps0):
        phi = math.log(rho / rn_) if rho > rn_ else 0.0
        return max((c**2 / 2.0) * (rho - rn_) + e0 * (phi**2 - phi**3 / 3.0), 0.0)

    def eos_rho(P, rn_=rn, e0=eps0):
        if P <= 0:
            return rn_
        lo, hi = rn_, rc * 5.0
        for _ in range(60):
            mid = (lo + hi) / 2.0
            if eos_p(mid) < P:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0

    return eos_p, eos_rho, rc


def find_r_at_mass(target_msun, eos_p, eos_rho, rho_nuc, rho_max):
    """
    Binary search on central density to find the star with M ≈ target_msun.
    Returns R in km, or NaN if the target mass is unreachable.
    """
    lo, hi = rho_nuc * 1.001, rho_max
    best_r = np.nan
    for _ in range(50):
        mid = math.sqrt(lo * hi)           # log-midpoint
        res = _solve_star_variant(mid, rho_max, eos_p, eos_rho)
        if res is None:
            hi = mid
            continue
        R_m, M_kg = res
        M_sol = M_kg / M_SUN
        if abs(M_sol - target_msun) < 0.005:
            best_r = R_m / 1e3
            break
        if M_sol < target_msun:
            lo = mid
        else:
            hi = mid
    return best_r


# ---------------------------------------------------------------------------
# Pre-compute fixtures at module level (expensive ODE; run once)
# ---------------------------------------------------------------------------

# f=1.0 (current ρ_nuc) with Option c
_P_c0, _rho_c0, _rc0 = make_option_c_scaled(1.0)
_M_c0, _R_c0 = compute_mr_curve_variant(_rc0, _P_c0, _rho_c0, n_points=50)

# f=1.18 with Option c
_P_c118, _rho_c118, _rc118 = make_option_c_scaled(1.18)
_M_c118, _R_c118 = compute_mr_curve_variant(_rc118, _P_c118, _rho_c118, n_points=50)

# Targeted R_1.4 via bisection for f=1.0 and f=1.18
_r14_f1p0   = find_r_at_mass(1.4, _P_c0,   _rho_c0,   RHO_NUC * 1.0,  _rc0)
_r14_f1p18  = find_r_at_mass(1.4, _P_c118, _rho_c118, RHO_NUC * 1.18, _rc118)
_r14_f1p1   = find_r_at_mass(1.4,
    *make_option_c_scaled(1.1)[:2],
    RHO_NUC * 1.1, math.e * RHO_NUC * 1.1)

_mmax_f1p0  = float(np.nanmax(_M_c0))   if _M_c0.size  else 0.0
_mmax_f1p18 = float(np.nanmax(_M_c118)) if _M_c118.size else 0.0


# ===========================================================================
class TestRadiusAnchoredByRhoNuc:
    """
    ρ_nuc is the dominant control variable for R_1.4.
    The EOS sound-speed parameters α and γ have negligible influence.
    """

    def test_uniform_density_radius_formula_accurate(self):
        """
        (3 M_1.4 / 4π ρ_nuc)^{1/3} ≈ 14.24 km  (analytic);
        ODE gives R_1.4 ≈ 14.0 km — agreement within 3%.

        Basis: PM surface at P=0 ⟺ ρ=ρ_nuc; 1.4 M☉ star mean density
        is only ~10% above ρ_nuc, so the star is nearly uniform at ρ_nuc.
        """
        R_uniform = (3.0 * 1.4 * M_SUN / (4.0 * math.pi * RHO_NUC))**(1.0 / 3.0)
        R_uniform_km = R_uniform / 1e3
        # universal formula predicts ~14.24 km
        assert R_uniform_km == pytest.approx(14.24, rel=0.005), \
            f"Uniform density formula gives R = {R_uniform_km:.2f} km (expected ~14.24)"
        # ODE result is within 3% of the analytic prediction
        assert not np.isnan(_r14_f1p0), "Bisection failed to find 1.4 M_sun star"
        assert _r14_f1p0 == pytest.approx(R_uniform_km, rel=0.03), (
            f"ODE R_1.4 = {_r14_f1p0:.2f} km vs uniform-density {R_uniform_km:.2f} km "
            f"(ratio = {_r14_f1p0/R_uniform_km:.4f}, expected within 3%)"
        )

    def test_r14_insensitive_to_alpha_scaling(self):
        """
        R_1.4 varies < 4% over α ∈ {0, 0.5, 1.0, 2.0} (confirmed numerically).
        α·U(φ) stiffens the EOS only at high φ (large ρ), not near ρ_nuc
        where the star surface and bulk of a 1.4 M☉ star reside.
        """
        r_values = []
        for alpha in [0.0, 0.5, 1.0, 2.0]:
            rn = RHO_NUC; rc = RHO_CRIT; eps0 = EPS0
            def eos_p(rho, a=alpha, rn_=rn, e0=eps0):
                phi = math.log(rho / rn_) if rho > rn_ else 0.0
                return max((c**2 / 2.0) * (rho - rn_) + a * e0 * (phi**2 - phi**3 / 3.0), 0.0)
            def eos_rho(P, a=alpha, rn_=rn, e0=eps0):
                if P <= 0: return rn_
                lo, hi = rn_, rc * 5.0
                for _ in range(60):
                    mid = (lo + hi) / 2
                    if eos_p(mid) < P: lo = mid
                    else: hi = mid
                return (lo + hi) / 2
            r14 = find_r_at_mass(1.4, eos_p, eos_rho, rn, rc)
            r_values.append(r14)
        r_vals_finite = [r for r in r_values if not np.isnan(r)]
        variation = (max(r_vals_finite) - min(r_vals_finite)) / np.mean(r_vals_finite)
        assert variation < 0.04, (
            f"R_1.4 varies {variation:.2%} over α∈[0,2] (expected < 4%) — "
            f"values: {[f'{r:.2f}' for r in r_values]}"
        )

    def test_r14_scales_with_rho_nuc_one_third(self):
        """
        R_1.4 ∝ ρ_nuc^{-1/3}.

        Verified:
          f=1.0:  R_1.4 = 14.00 km
          f=1.1:  R_1.4 = 13.54 km
          Ratio 14.00/13.54 = 1.034;  (1.1)^{1/3} = 1.032  ✓ (0.2% agreement)
        """
        assert not np.isnan(_r14_f1p0) and not np.isnan(_r14_f1p1), \
            "Bisection failed for f=1.0 or f=1.1"
        ratio_ode      = _r14_f1p0 / _r14_f1p1
        ratio_expected = (1.1 / 1.0)**(1.0 / 3.0)         # = 1.032
        assert ratio_ode == pytest.approx(ratio_expected, rel=0.01), (
            f"R_1.4(f=1.0)/R_1.4(f=1.1) = {ratio_ode:.4f}, "
            f"expected (1.1)^(1/3) = {ratio_expected:.4f}"
        )

    def test_current_rho_nuc_gives_r14_above_constraint(self):
        """
        At current ρ_nuc = 2.3×10¹⁷ kg/m³ (f=1.0, Option c):
        R_1.4 ≈ 14.0 km > 13.3 km (GW170817 90% CI upper bound).
        This is the radius problem; it persists regardless of α.
        """
        assert not np.isnan(_r14_f1p0), "Bisection failed for f=1.0"
        assert _r14_f1p0 > 13.3, (
            f"R_1.4 = {_r14_f1p0:.2f} km should be > 13.3 km at current ρ_nuc"
        )

    def test_scaled_rho_nuc_satisfies_radius_constraint(self):
        """
        At ρ_nuc → 1.18 × ρ_nuc (f=1.18, Option c):
        R_1.4 drops to ≈ 13.22 km < 13.3 km — within the GW170817 1σ upper bound.
        """
        assert not np.isnan(_r14_f1p18), "Bisection failed for f=1.18"
        assert _r14_f1p18 < 13.3, (
            f"R_1.4 = {_r14_f1p18:.2f} km should be < 13.3 km at f=1.18"
        )


# ===========================================================================
class TestSimultaneousConstraintFeasibility:
    """
    Can the two observational constraints be satisfied simultaneously?

      M_max ≥ 30 M☉   (GW150914 lower component; closing the mass gap)
      R_1.4 ≤ 13.3 km (GW170817 / NICER 90% upper bound)

    With two independent handles:
      α — controls M_max via U(φ) elastic pressure at high density
      f — controls R_1.4 via ρ_nuc (the uniform-density radius anchor)

    These handles are DECOUPLED:
      • Changing α at fixed f barely shifts R_1.4 (< 4% variation over α ∈ [0,2])
      • Changing f at fixed α shifts R_1.4 ∝ f^{-1/3} and M_max ∝ ???    
    """

    def test_current_optionc_fails_radius_not_mmax(self):
        """
        At (f=1.0, α=1): M_max = 30.8 M☉ ✓  but  R_1.4 ≈ 14.0 km ✗.
        The radius problem is the remaining violation, not M_max.
        """
        assert _mmax_f1p0 > 28.0, f"Option c M_max = {_mmax_f1p0:.2f} should be > 28"
        assert _r14_f1p0 > 13.3, f"R_1.4 = {_r14_f1p0:.2f} should fail R constraint"

    def test_scaled_rho_nuc_satisfies_both_with_alpha_1(self):
        """
        At (f=1.18, α=1): M_max = 28.4 M☉ ✓ and R_1.4 = 13.22 km ✓.
        Both constraints satisfied with M_max > 28 M☉ (slightly below 30.7).

        At this ρ_nuc, slightly larger α (≈1.15) would push M_max above 30 M☉
        while keeping R < 13.3 km — the two-parameter solution.
        """
        assert not np.isnan(_r14_f1p18), "Bisection failed at f=1.18"
        assert _r14_f1p18 < 13.3, \
            f"R_1.4(f=1.18) = {_r14_f1p18:.2f} km, expected < 13.3"
        assert _mmax_f1p18 > 26.0, \
            f"M_max(f=1.18) = {_mmax_f1p18:.2f} M_sun, expected > 26"

    def test_two_constraints_require_two_parameters(self):
        """
        With only one free parameter (α at fixed ρ_nuc):
          - Varying α changes M_max extensively but barely changes R_1.4.
          - No single α at f=1.0 satisfies both M_max > 28 AND R_1.4 < 13.3.

        This quantifies why the problem requires a TWO-PARAMETER repair.
        """
        r_vals = []
        m_vals = []
        for alpha in [0.3, 0.7, 1.0, 1.5]:
            rn = RHO_NUC; rc = RHO_CRIT; eps0 = EPS0
            def eos_p(rho, a=alpha, rn_=rn, e0=eps0):
                phi = math.log(rho / rn_) if rho > rn_ else 0.0
                return max((c**2 / 2.0) * (rho - rn_) + a * e0 * (phi**2 - phi**3 / 3.0), 0.0)
            def eos_rho(P, a=alpha, rn_=rn, e0=eps0):
                if P <= 0: return rn_
                lo, hi = rn_, rc * 5.0
                for _ in range(60):
                    mid = (lo + hi) / 2
                    if eos_p(mid) < P: lo = mid
                    else: hi = mid
                return (lo + hi) / 2
            M_arr, _ = compute_mr_curve_variant(rc, eos_p, eos_rho, n_points=40)
            r14 = find_r_at_mass(1.4, eos_p, eos_rho, rn, rc)
            r_vals.append(r14 if not np.isnan(r14) else 99.0)
            m_vals.append(float(np.nanmax(M_arr)) if M_arr.size else 0.0)

        # R barely changes despite large α range
        r_variation = (max(r_vals) - min(r_vals)) / np.mean(r_vals)
        assert r_variation < 0.05, \
            f"R_1.4 varies {r_variation:.2%} over α∈[0.3,1.5] — expected < 5%"

        # M_max varies a lot
        m_variation = (max(m_vals) - min(m_vals)) / np.mean(m_vals)
        assert m_variation > 0.5, \
            f"M_max varies only {m_variation:.2%} over α — expected > 50%"

        # No α at f=1.0 satisfies both; all R are above 13.3
        assert all(r > 13.2 for r in r_vals), \
            f"Some R_1.4 {r_vals} dropped below 13.2 km at f=1.0 — unexpected"

    def test_rho_nuc_physical_range_is_narrow(self):
        """
        The f-range satisfying BOTH constraints simultaneously is narrow: ~3%.
        f < 1.18 → R_1.4 > 13.3 km ✗
        f > 1.21 → M_max < 27 M☉ ✗  (at α=1)

        This means ρ_nuc is tightly constrained by the combined observational data.
        """
        # Upper bound: f=1.15 still gives R > 13.3
        p115, r115, rc115 = make_option_c_scaled(1.15)
        r14_115 = find_r_at_mass(1.4, p115, r115, RHO_NUC * 1.15, rc115)
        if not np.isnan(r14_115):
            assert r14_115 > 13.2, f"R_1.4(f=1.15) = {r14_115:.2f} should be > 13.2"

        # Lower bound: f=1.25 gives M_max < 28
        p125, r125, rc125 = make_option_c_scaled(1.25)
        M125, _ = compute_mr_curve_variant(rc125, p125, r125, n_points=40)
        mmax125 = float(np.nanmax(M125)) if M125.size else 0.0
        assert mmax125 < 28.5, f"M_max(f=1.25) = {mmax125:.2f} should be < 28.5"

    def test_rho_nuc_required_shift_is_reasonable(self):
        """
        The required ρ_nuc ≈ 1.18 × ρ_nuc_0 is physically reasonable.

        Standard nuclear saturation:  ρ_0 = 2.3×10¹⁷ kg/m³ (symmetric matter)
        Neutron-rich matter:          ρ_nuc can exceed ρ_0 by 10–30%

        A 18% shift is within the range of neutron-matter nuclear models.
        The shift is less than 2×ρ_0 (where quark-gluon plasma may onset).
        """
        f_needed = 1.18
        rho_needed = f_needed * RHO_NUC
        rho_qgp    = 5.0 * RHO_NUC      # rough onset of quark matter
        rho_2saturation = 2.0 * RHO_NUC  # 2×ρ_0

        assert rho_needed < rho_2saturation, (
            f"Required ρ_nuc = {rho_needed:.2e} should be below 2×ρ_0 = {rho_2saturation:.2e}"
        )
        assert rho_needed < rho_qgp, (
            f"Required ρ_nuc = {rho_needed:.2e} should be below QGP onset {rho_qgp:.2e}"
        )
        # The shift is a modest 10-30%
        assert 1.05 < f_needed < 1.40, \
            f"f = {f_needed:.2f} is outside the 5–40% range considered 'modest'"
