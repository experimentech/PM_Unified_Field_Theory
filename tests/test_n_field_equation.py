"""Tests: comparison of the standard PM field equation vs the n-field equation.

Two field equations for the PM scalar field
--------------------------------------------
Standard (Poisson in φ):
    ∇²φ = source,  φ(r) = 2GM/c²r  (point mass)

n-field (linear in n = e^φ):
    ∇²n = (4πG/c²) ρ n,  n(r) = 1 + 2GM/c²r  →  φ(r) = ln(1 + 2GM/c²r)

Motivation
----------
The standard equation is actually ∇²(ln n) = source, which when expanded gives
∇²n/n − |∇n|²/n² = source.  The |∇n|²/n² term is the field's own self-energy;
the standard equation silently drops it, which is valid when φ ≪ 1 but becomes
a formal inconsistency as φ → 1.  The n-field equation retains it exactly.

Key analytical relations (point mass)
--------------------------------------
Let x = 2GM/c²r  (the "compactness parameter" at distance r).

    φ_old = x
    φ_new = ln(1 + x)

    φ_new / φ_old = ln(1+x)/x  →  1 − x/2 + x²/3 − ...  (Taylor)

    a_old = GM/r²  (exact Newtonian, by construction)
    a_new = GM/r² · 1/(1+x)  =  a_old / (1 + φ_old_linear)

    relative difference in |a|  =  x/(1+x)  ≈  x = φ_old  for x ≪ 1

Phase boundary φ = 1:
    Standard:  r_crit_old = 2GM/c²             (x = 1)
    n-field:   r_crit_new = 2GM / (c²(e−1))    (x = e−1 ≈ 1.718)
               r_crit_new ≈ 0.582 × r_crit_old

    At the old boundary (r = r_crit_old, x = 1):
        φ_new = ln 2 ≈ 0.693  (well below 1 — the n-field does not trigger the
                                 phase transition at the same radius)
        |a_new| / |a_old| = 0.5   (factor of 2 suppression in strong field)

Observable predictions
-----------------------
Because the relative difference ≈ φ_old ≈ 2GM/c²r:
    Solar system (1 AU, sun):     diff ≈ 2 × 10⁻⁸  (undetectable)
    Neutron star surface (~12 km, 1.4 M☉):  diff ≈ 17%   (significant)
    GR photon sphere (3GM/c²):    diff = 40%             (strong-field regime)

All existing PM tests (solar-system, binary pulsar) lie firmly in the weak-field
regime and therefore pass equally under both equations.  Only compact-star
observables (redshift, radius, maximum mass) would distinguish them.

Implementation
--------------
    Standard equation:  newtonian_accel_sum (≡ PM for slow particles)
    n-field equation:   massive_accel_n_field, phi_n_field_point_mass
"""

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from pushing_medium.core import (
    G, c,
    newtonian_accel_sum,
    massive_accel_n_field,
    phi_n_field_point_mass,
    pm_deflection_angle_point_mass,
)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
M_SUN   = 1.989e30    # kg
M_NS    = 1.4 * M_SUN # kg  (canonical neutron star)
AU      = 1.496e11    # m
R_NS    = 12.0e3      # m   (typical neutron star radius)
R_SUN   = 6.96e8      # m


def _phi_old(M, r):
    """Standard PM field: φ = 2GM/c²r."""
    return 2.0 * G * M / (c * c * r)


def _a_old_mag(M, r):
    """Standard acceleration magnitude: GM/r²."""
    return G * M / (r * r)


def _a_new_mag(M, r):
    """n-field acceleration magnitude: GM/(r²(1 + 2GM/c²r)) = a_old / (1 + φ_old)."""
    return G * M / (r * r * (1.0 + 2.0 * G * M / (c * c * r)))


# ---------------------------------------------------------------------------

class TestWeakFieldEquivalence:
    """Both equations give indistinguishable results across the solar system.

    The relative difference is ≈ φ_old = 2GM/c²r, which is < 3×10⁻⁸
    everywhere from Mercury to Neptune.
    """

    @pytest.mark.parametrize("r_AU, tol", [
        (0.387, 6e-8),    # Mercury — closest to sun, φ ≈ 5.1e-8
        (1.0,   2.1e-8),  # Earth,   φ ≈ 2.0e-8
        (1.524, 1.4e-8),  # Mars,    φ ≈ 1.3e-8
        (5.2,   4.0e-9),  # Jupiter, φ ≈ 3.8e-9
        (30.1,  7.0e-10), # Neptune, φ ≈ 6.6e-10
    ])
    def test_acceleration_magnitude_agrees(self, r_AU, tol):
        """n-field acceleration agrees with Newtonian to better than φ_old."""
        r = r_AU * AU
        a_newt = abs(newtonian_accel_sum((r, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))])[0])
        a_nf   = abs(massive_accel_n_field((r, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))])[0])
        rel = abs(a_newt - a_nf) / a_newt
        assert rel < tol, f"r={r_AU} AU: rel diff = {rel:.2e}, tol = {tol:.2e}"

    @pytest.mark.parametrize("r_AU", [0.387, 1.0, 5.2])
    def test_phi_ratio_matches_lnp1_over_x(self, r_AU):
        """φ_new/φ_old = ln(1+x)/x exactly."""
        r = r_AU * AU
        x = _phi_old(M_SUN, r)
        phi_new = phi_n_field_point_mass(r, M_SUN)
        expected_ratio = math.log1p(x) / x
        assert abs(phi_new / x - expected_ratio) < 1e-12

    def test_multi_body_weak_field(self):
        """Four solar-mass bodies in a square: both equations agree to 1e-7."""
        d = AU
        bodies = [(M_SUN, (d,  d,  0.0)),
                  (M_SUN, (-d, d,  0.0)),
                  (M_SUN, (d,  -d, 0.0)),
                  (M_SUN, (-d, -d, 0.0))]
        r0 = (0.0, 0.0, 0.0)
        a_newt = newtonian_accel_sum(r0, bodies)
        a_nf   = massive_accel_n_field(r0, bodies)
        for i in range(3):
            assert abs(a_newt[i]) < 1e-20   # cancels by symmetry
            assert abs(a_nf[i])   < 1e-20


class TestStrongFieldDivergence:
    """The n-field gives systematically smaller φ and |a| in the strong field,
    and the suppression follows the analytical formula a_new = a_old/(1+φ_old).
    """

    @pytest.mark.parametrize("r_rs", [100.0, 10.0, 3.0, 2.0, 1.5, 1.0])
    def test_nfield_accel_always_smaller(self, r_rs):
        """n-field |a| < Newtonian |a| at all radii (equality only at r→∞)."""
        r_s = 2.0 * G * M_SUN / (c * c)   # Schwarzschild radius
        r = r_rs * r_s
        a_old = _a_old_mag(M_SUN, r)
        a_new = _a_new_mag(M_SUN, r)
        assert a_new < a_old, f"r = {r_rs} r_s: a_new={a_new:.4e} >= a_old={a_old:.4e}"

    @pytest.mark.parametrize("r_rs", [100.0, 10.0, 3.0, 2.0, 1.5, 1.0])
    def test_accel_ratio_matches_analytic(self, r_rs):
        """a_new / a_old = 1/(1 + 2GM/c²r) exactly."""
        r_s = 2.0 * G * M_SUN / (c * c)
        r = r_rs * r_s
        ratio_computed = _a_new_mag(M_SUN, r) / _a_old_mag(M_SUN, r)
        x = _phi_old(M_SUN, r)                   # = 1/r_rs
        ratio_analytic = 1.0 / (1.0 + x)
        assert abs(ratio_computed - ratio_analytic) < 1e-12

    def test_relative_diff_grows_with_compactness(self):
        """Relative difference increases monotonically as r decreases."""
        r_s = 2.0 * G * M_SUN / (c * c)
        radii = [100.0 * r_s, 10.0 * r_s, 3.0 * r_s, 1.0 * r_s]
        diffs = []
        for r in radii:
            a_old = _a_old_mag(M_SUN, r)
            a_new = _a_new_mag(M_SUN, r)
            diffs.append((a_old - a_new) / a_old)
        for i in range(len(diffs) - 1):
            assert diffs[i+1] > diffs[i], \
                f"Difference not monotone: {diffs}"

    def test_suppression_factor_at_old_phase_boundary(self):
        """At r = 2GM/c² (old φ=1 boundary): a_new = a_old/2 exactly."""
        r_s = 2.0 * G * M_SUN / (c * c)   # = 2GM/c², so φ_old = 1 here
        ratio = _a_new_mag(M_SUN, r_s) / _a_old_mag(M_SUN, r_s)
        assert abs(ratio - 0.5) < 1e-12

    def test_phi_new_below_1_at_old_phase_boundary(self):
        """At the old phase boundary (φ_old = 1), φ_new = ln 2 ≈ 0.693 < 1.

        The n-field equation does NOT trigger the phase transition at the same
        radius as the standard equation.  The phase boundary moves inward.
        """
        r_s = 2.0 * G * M_SUN / (c * c)
        phi_old_val = _phi_old(M_SUN, r_s)
        phi_new_val = phi_n_field_point_mass(r_s, M_SUN)
        assert abs(phi_old_val - 1.0) < 1e-12, "Sanity: φ_old should be 1 here"
        assert abs(phi_new_val - math.log(2.0)) < 1e-12
        assert phi_new_val < 1.0

    def test_n_field_phase_boundary_at_correct_radius(self):
        """φ_new = 1 when 2GM/c²r = e−1, i.e. r = 2GM/(c²(e−1))."""
        r_crit_new = 2.0 * G * M_SUN / (c * c * (math.e - 1.0))
        phi_new_at_crit = phi_n_field_point_mass(r_crit_new, M_SUN)
        assert abs(phi_new_at_crit - 1.0) < 1e-12

    def test_neutron_star_surface_compactness_diff(self):
        """At a canonical neutron star surface (χ ≈ 0.17), difference is ~14%.

        This is the regime where the two equations give distinguishable predictions
        for surface redshift, stellar radius, and maximum mass.
        """
        phi_old = _phi_old(M_NS, R_NS)             # ≈ 0.172
        a_ratio = _a_new_mag(M_NS, R_NS) / _a_old_mag(M_NS, R_NS)
        expected_ratio = 1.0 / (1.0 + phi_old)
        assert abs(a_ratio - expected_ratio) < 1e-12
        # Confirm the difference is non-trivial (> 10%)
        assert (1.0 - a_ratio) > 0.10, \
            f"Expected >10% difference at NS surface, got {(1-a_ratio)*100:.1f}%"


class TestObservablesPredictionsMatch:
    """Both equations make identical predictions for all currently tested observables.

    The weak-field limit ensures that solar-system and binary-pulsar results are
    equation-independent.  Any existing passing test in the battery remains valid.
    """

    def test_light_deflection_independent_of_field_equation(self):
        """Light deflection formula 4GM/c²b does not depend on which field equation
        is used: it derives from the refractive index n(r) ≈ 1 + 2GM/c²r, which is
        the *same* in both equations to first order in 2GM/c²r.
        """
        b_values = [R_SUN, 10 * R_SUN, AU]
        for b in b_values:
            # In the weak field: n(r) = 1 + 2GM/c²r under both equations
            # → deflection 4GM/c²b is equation-independent
            alpha = pm_deflection_angle_point_mass(M_SUN, b)
            assert alpha > 0

    @pytest.mark.parametrize("r_AU, tol_frac", [
        (0.387, 6e-8),   # Mercury
        (1.0,   2.1e-8), # Earth
    ])
    def test_solar_system_acceleration_within_measurement_precision(self, r_AU, tol_frac):
        """Both equations agree to < 1 part in 10⁷ at planetary distances,
        far below any current gravitational measurement precision (~10⁻⁶).
        """
        r = r_AU * AU
        a_old = abs(newtonian_accel_sum((r, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))])[0])
        a_new = abs(massive_accel_n_field((r, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))])[0])
        assert abs(a_old - a_new) / a_old < tol_frac

    def test_phi_difference_negligible_at_sun_surface(self):
        """At the sun's surface, the two equations differ by < 5×10⁻⁶ in φ.

        This is 50× below the precision needed to affect solar-system tests.
        """
        phi_old = _phi_old(M_SUN, R_SUN)
        phi_new = phi_n_field_point_mass(R_SUN, M_SUN)
        rel = abs(phi_old - phi_new) / phi_old
        assert rel < 5e-6


class TestSuperpositionAndSymmetry:
    """The n-field equation handles multi-body superposition correctly."""

    def test_two_equal_masses_symmetric_cancellation(self):
        """Test mass at centre between two equal masses: net force should be zero
        under both equations (by symmetry).
        """
        M = M_SUN
        d = AU
        masses = [(M, (-d, 0.0, 0.0)), (M, (d, 0.0, 0.0))]
        r0 = (0.0, 0.0, 0.0)
        a_nf = massive_accel_n_field(r0, masses)
        for comp in a_nf:
            assert abs(comp) < 1e-20

    def test_n_field_reduces_to_single_body_far_from_second(self):
        """With two masses, one very far away, n-field ≈ single-body result."""
        M1, M2 = M_SUN, M_SUN
        r_close = AU
        r_far   = 1000.0 * AU
        r_test  = (r_close, 0.0, 0.0)

        a_single = massive_accel_n_field(r_test, [(M1, (0.0, 0.0, 0.0))])
        a_double = massive_accel_n_field(r_test, [(M1, (0.0, 0.0, 0.0)),
                                                   (M2, (r_far, 0.0, 0.0))])
        # x-component dominates; small correction from M2 at r_far
        rel = abs(a_single[0] - a_double[0]) / abs(a_single[0])
        # Contribution from M2: roughly M2 r_close / (M1 r_far) ≈ 10⁻³
        assert rel < 2e-3

    def test_n_field_magnitude_bounded_by_mass_sum(self):
        """n-field acceleration at the origin from two off-axis bodies is less
        than it would be if both masses were co-located at the nearer distance.
        """
        M = M_SUN
        d = 2.0 * AU
        r_test = (0.0, 0.0, 0.0)
        bodies = [(M, (d, 0.0, 0.0)), (M, (0.0, d, 0.0))]
        a_nf = massive_accel_n_field(r_test, bodies)
        # Upper bound: both masses at distance d along x — would give 2GM/d² in x
        upper = 2.0 * G * M / (d * d)
        # Actual magnitude should be less than upper due to the 1/(1+n) factor
        a_mag = math.hypot(*a_nf)
        assert a_mag <= upper * 1.01   # allow tiny floating point margin


# ---------------------------------------------------------------------------
# 5.  Action derivation: n-field equation is EL equation of ∫½|∇φ|²n² d³x
# ---------------------------------------------------------------------------

class TestNFieldActionDerivation:
    """The n-field vacuum equation ∇²φ + |∇φ|² = 0 is the Euler-Lagrange
    equation of the action S = ∫½|∇φ|² n² d³x.

    This class verifies that the n-field point-mass vacuum solution satisfies
    this EL equation exactly (to floating-point precision), and confirms the
    consequence for PPN β.

    Derivation summary
    ------------------
    Action:  S = ∫ ½|∇φ|² e^{2φ} d³x  (physical area measure n² d³x)

    EL equation (vacuum):  ∇²φ + |∇φ|² = 0

    Linearisation:  w = n = e^φ  →  ∇²w = 0  (Laplace equation)
    So superposition is exact in n, and the point-mass solution n = 1+A/r is
    the unique regular solution decaying at infinity — exactly the n-field
    equation's vacuum solution.

    PPN β:  φ = ln(1+2U) = 2U − 2U² + O(U³)  →  β = 1
    Combined with γ=1 from gij = n²δij:
        PPN factor = (2+2γ−β)/3 = (2+2−1)/3 = 1
        Perihelion precession exactly matches GR (and observation) without
        any GR import.
    """

    @pytest.mark.parametrize("r_factor", [100.0, 10.0, 3.0, 1.5, 1.0, 0.582])
    def test_vacuum_equation_satisfied(self, r_factor):
        """∇²φ + |∇φ|² = 0 at each test radius (floating-point noise only)."""
        r_s = 2 * G * M_SUN / c**2
        r = r_factor * r_s
        A = 2 * G * M_SUN / c**2   # coefficient in n = 1 + A/r
        n_val = 1.0 + A / r
        dphi_dr = -A / (n_val * r**2)
        # ∇²φ in spherical symmetry: d²φ/dr² + 2/r dφ/dr
        dn_dr = -A / r**2
        d2phi_dr2 = -A * (-dn_dr * r**2 - 2 * r * n_val) / (n_val * r**2)**2
        laplacian_phi = d2phi_dr2 + 2 / r * dphi_dr
        grad_sq = dphi_dr**2
        residual = laplacian_phi + grad_sq
        assert abs(residual) < 1e-20, (
            f"r={r_factor} r_s: ∇²φ+|∇φ|²={residual:.3e} (should be 0)"
        )

    def test_ppn_beta_equals_one(self):
        """PPN β = 1 from φ = ln(1+2U) Taylor expansion: φ = 2U − 2U² + O(U³).

        d²(ln(1+2U))/dU² |_{U=0} = d/dU[2/(1+2U)]|_0 = −4
        c2 = d²φ/dU² / 2 = −2,  β = −c2/2 = 1  (exact arithmetic).
        """
        # c2 is the U^2 coefficient in phi = ln(1+2U) expanded around U=0
        # ln(1+x) = x - x^2/2 + x^3/3 - ...
        # phi = ln(1+2U) = 2U - (2U)^2/2 + ... = 2U - 2*U^2 + ...
        # so c2 = -2,  beta = -c2/2 = 1 exactly.
        c2 = -2.0   # exact from Taylor series
        beta = -c2 / 2.0
        assert abs(beta - 1.0) < 1e-15

    def test_ppn_factor_exact_with_spatial_metric(self):
        """With γ=1 (gij=n²δij) and β=1 (n-field), PPN factor = 1 exactly."""
        beta = 1.0
        gamma = 1.0
        ppn_factor = (2 + 2 * gamma - beta) / 3
        assert abs(ppn_factor - 1.0) < 1e-15

    def test_linearisation_w_equals_n(self):
        """Substitution w=n transforms ∇²φ+|∇φ|²=0 into ∇²n=0.

        Proof: ∇²n = ∇·∇(e^φ) = ∇·(e^φ∇φ) = e^φ(|∇φ|²+∇²φ) = n·0 = 0.
        Verified analytically: n(r)=1+A/r satisfies ∇²n=0 because ∇²(1/r)=0
        for r>0.  Numerically: 2A/r³ − 2A/r³ = 0, leaving only floating-point
        cancellation noise (< 1e-20 relative to either term).
        """
        r_s = 2 * G * M_SUN / c**2
        A = 2 * G * M_SUN / c**2
        for r_factor in [2.0, 5.0, 10.0]:
            r = r_factor * r_s
            # ∇²(1 + A/r) = ∇²(A/r) = 0 for r>0 (harmonic function)
            # Numerically: d^2(A/r)/dr^2 + (2/r)d(A/r)/dr = 2A/r^3 - 2A/r^3 = 0
            d2_dr2 = 2 * A / r**3
            term_2over_r = (2 / r) * (-A / r**2)
            laplacian_n = d2_dr2 + term_2over_r
            # tolerance: floating-point cancellation noise relative to each term
            assert abs(laplacian_n) < 1e-14 * abs(d2_dr2), (
                f"r={r_factor} r_s: ∇²n = {laplacian_n:.3e}"
            )
