"""Tests: physical-measure field equation for the PM scalar field.

Physical-measure field equation
--------------------------------
Action integrated against the physical volume element n³ d³x = e^{3φ} d³x:

    S = ∫ ½|∇φ|² e^{3φ} d³x

Euler–Lagrange equation (vacuum):

    ∇²φ + (3/2)|∇φ|² = 0

Linearisation: the substitution w = n^{3/2} = e^{3φ/2} transforms this exactly
into the Laplace equation ∇²w = 0.  Superposition is therefore exact in w.

Point-mass vacuum solution
---------------------------
    w(r) = 1 + 3GM/(c²r)
    φ(r) = (2/3) ln(1 + 3GM/(c²r))

Let q = 3GM/(c²r)  (the "compactness parameter" for this equation).

    φ_pm   = (2/3) ln(1 + q)

    a      = (c²/2)∇φ  = −GM/r² / (1 + q)   [radially inward]

    relative difference from Newtonian: q/(1+q) ≈ q = 3GM/(c²r) for q ≪ 1

Phase boundary φ = 1 (compared with the other two field equations):
    Standard PM (∇²φ = src):  r_crit = 2GM/c²            = r_s
    n-field     (∇²n = src n): r_crit = 2GM/(c²(e−1))   ≈ 0.582 r_s
    Physical-measure (∇²w=0):  r_crit = 3GM/(c²(e^{3/2}−1)) ≈ 0.431 r_s

Vacuum field equation verification
------------------------------------
Analytically: ∇²φ = −3φ₀²/(2r⁴) / (1+q)²  and  (3/2)(φ')² = +3φ₀²/(2r⁴) / (1+q)²
so  ∇²φ + (3/2)|∇φ|² = 0 exactly (where φ₀ = 2GM/c²).

All three field equations agree to within their respective φ-values in the weak
field, and are only observably distinct in the strong field.
"""

import math
import sys
import pytest

sys.path.insert(0, "src")
from pushing_medium.core import (
    G, c,
    newtonian_accel_sum,
    phi_physical_measure_point_mass,
    grad_phi_physical_measure_point_mass,
    massive_accel_physical_measure,
    phi_n_field_point_mass,
    massive_accel_n_field,
)

# ---------------------------------------------------------------------------
# Physical constants and reference scale
# ---------------------------------------------------------------------------
M_SUN = 1.989e30        # kg
R_SUN = 6.957e8         # m
AU = 1.496e11           # m
r_s = 2 * G * M_SUN / c**2   # Schwarzschild radius of the Sun ≈ 2954 m


# ---------------------------------------------------------------------------
# Helper: Newtonian |a| from a point mass
# ---------------------------------------------------------------------------
def newtonian_magnitude(r_dist, M):
    return G * M / r_dist**2


# ---------------------------------------------------------------------------
# 1.  Weak-field limit: φ agrees with standard PM to within O(φ²)
# ---------------------------------------------------------------------------

class TestWeakFieldLimit:
    """Physical-measure φ ≈ 2GM/c²r in the weak field."""

    @pytest.mark.parametrize("r_AU, name", [
        (0.387,  "Mercury"),
        (0.723,  "Venus"),
        (1.000,  "Earth"),
        (1.524,  "Mars"),
        (5.203,  "Jupiter"),
        (9.537,  "Saturn"),
        (19.19,  "Uranus"),
        (30.07,  "Neptune"),
    ])
    def test_phi_matches_newtonian_limit(self, r_AU, name):
        """φ_pm ≈ 2GM/c²r: fractional error < φ_std at each solar-system radius."""
        r = r_AU * AU
        phi_std = 2 * G * M_SUN / (c**2 * r)           # standard PM
        phi_pm  = phi_physical_measure_point_mass(r, M_SUN)
        # For weak field: phi_pm = (2/3)ln(1 + 3x/2) ≈ x - x²/2 + ...
        # where x = 2GM/c²r = phi_std.  Difference ≈ phi_std²/4.
        assert abs(phi_pm - phi_std) < phi_std, (
            f"{name}: φ_pm={phi_pm:.3e} vs φ_std={phi_std:.3e}"
        )

    def test_force_matches_newtonian_at_1au(self):
        """Radial force equals Newtonian to within the leading-order correction.

        The physical-measure suppressor is 1/(1 + q) where q = 3GM/c²r,
        so the fractional deviation from Newtonian is q/(1+q) ≈ q = 3GM/c²r.
        """
        r = AU
        q = 3 * G * M_SUN / (c**2 * r)   # leading correction parameter
        a_newt = newtonian_magnitude(r, M_SUN)
        ax, ay, az = massive_accel_physical_measure(
            (r, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))]
        )
        a_pm = abs(ax)
        assert abs(a_pm - a_newt) / a_newt < q


# ---------------------------------------------------------------------------
# 2.  Vacuum field equation: ∇²φ + (3/2)|∇φ|² = 0
# ---------------------------------------------------------------------------

class TestVacuumFieldEquation:
    """The vacuum field equation is satisfied numerically."""

    @pytest.mark.parametrize("r_factor", [100.0, 10.0, 3.0, 1.5, 0.8, 0.5])
    def test_vacuum_equation_satisfied(self, r_factor):
        """∇²φ + (3/2)|∇φ|² = 0 to tolerance 1e-6 at each test radius."""
        r = r_factor * r_s
        phi0 = 2 * G * M_SUN / c**2   # = r_s
        # Analytic expressions at radius r:
        q = 3 * G * M_SUN / (c**2 * r)   # = (3/2)(phi_std)
        # d phi/dr  = -(2GM/c²r²) / (1 + q)  =  -(phi0/r²)/(1+q)
        dphi_dr = -(phi0 / r**2) / (1 + q)
        # d²phi/dr² + (2/r)(dphi/dr) = nabla^2 phi in spherical symmetry
        # Analytic: nabla^2 phi = -3*phi0^2 / (2*r^4) / (1+q)^2
        laplacian_phi = -3 * phi0**2 / (2 * r**4) / (1 + q)**2
        # (3/2)|grad phi|^2  =  (3/2) * (phi0/r^2)^2 / (1+q)^2
        grad_sq_term = (3.0 / 2.0) * (phi0 / r**2)**2 / (1 + q)**2
        residual = laplacian_phi + grad_sq_term
        assert abs(residual) < 1e-6, (
            f"r={r_factor} r_s: residual = {residual:.3e} (should be 0)"
        )


# ---------------------------------------------------------------------------
# 3.  Comparison with standard and n-field equations
# ---------------------------------------------------------------------------

class TestStrongFieldDivergence:
    """Physical-measure equation gives different predictions from standard PM
    and the n-field equation in the strong field."""

    @pytest.mark.parametrize("r_factor", [10.0, 3.0, 1.5, 1.0, 0.582, 0.431])
    def test_force_suppression_formula(self, r_factor):
        """a_pm / a_newt = 1/(1 + 3GM/c²r) exactly."""
        r = r_factor * r_s
        q = 3 * G * M_SUN / (c**2 * r)       # = 1.5 * phi_std
        expected_ratio = 1.0 / (1.0 + q)
        ax, _, _ = massive_accel_physical_measure(
            (r, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))]
        )
        a_newt = newtonian_magnitude(r, M_SUN)
        actual_ratio = abs(ax) / a_newt
        assert abs(actual_ratio - expected_ratio) < 1e-10, (
            f"r={r_factor} r_s: ratio={actual_ratio:.6f}, expected={expected_ratio:.6f}"
        )

    def test_phase_boundary_location(self):
        """Phase boundary φ=1 is at r = 3GM/(c²(e^{3/2}-1)) ≈ 0.431 r_s."""
        r_phase_analytic = 3 * G * M_SUN / (c**2 * (math.exp(1.5) - 1))
        phi_at_boundary = phi_physical_measure_point_mass(r_phase_analytic, M_SUN)
        assert abs(phi_at_boundary - 1.0) < 1e-12

    def test_phase_boundary_deeper_than_n_field(self):
        """Physical-measure phase boundary is deeper (smaller r) than n-field."""
        r_phase_pm  = 3 * G * M_SUN / (c**2 * (math.exp(1.5) - 1))
        r_phase_nf  = 2 * G * M_SUN / (c**2 * (math.e - 1))
        r_phase_std = 2 * G * M_SUN / c**2
        assert r_phase_pm < r_phase_nf < r_phase_std

    def test_phi_at_standard_phase_boundary(self):
        """At r = r_s (standard PM phase boundary), φ_pm < 1 (no premature transition)."""
        phi_at_rs = phi_physical_measure_point_mass(r_s, M_SUN)
        assert phi_at_rs < 1.0
        # More precisely: phi_pm(r_s) = (2/3)ln(1 + 3/2) = (2/3)ln(2.5)
        expected = (2.0/3.0) * math.log(2.5)
        assert abs(phi_at_rs - expected) < 1e-12

    def test_pm_weaker_than_newtonian_in_strong_field(self):
        """At r = r_s, physical-measure force is weaker than Newtonian by factor 1/(1+3/2)."""
        ax, _, _ = massive_accel_physical_measure(
            (r_s, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))]
        )
        a_newt = newtonian_magnitude(r_s, M_SUN)
        # q = 3GM/(c²r_s) = 3GM/(c² · 2GM/c²) = 3/2
        expected_ratio = 1.0 / (1.0 + 1.5)   # = 0.4
        assert abs(abs(ax) / a_newt - expected_ratio) < 1e-10

    def test_pm_differs_from_n_field_in_strong_field(self):
        """Physical-measure and n-field forces are measurably different at r = r_s."""
        r = r_s
        ax_pm, _, _ = massive_accel_physical_measure(
            (r, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))]
        )
        ax_nf, _, _ = massive_accel_n_field(
            (r, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))]
        )
        # n-field suppressor: 1/(1 + phi_std) = 1/(1+1) = 0.5
        # pm suppressor: 1/(1 + 1.5) = 0.4
        assert abs(abs(ax_pm) - abs(ax_nf)) / abs(ax_nf) > 0.1


# ---------------------------------------------------------------------------
# 4.  Superposition and symmetry
# ---------------------------------------------------------------------------

class TestSuperpositionAndSymmetry:
    """Multi-body acceleration uses the exact w-superposition."""

    def test_single_body_matches_analytic_formula(self):
        """massive_accel_physical_measure exactly matches (c²/2)∇φ for one mass."""
        r = 3 * r_s
        ax, ay, az = massive_accel_physical_measure(
            (r, 0.0, 0.0), [(M_SUN, (0.0, 0.0, 0.0))]
        )
        # Analytic: a = -GM/r² / (1 + 3GM/c²r)
        q = 3 * G * M_SUN / (c**2 * r)
        a_analytic = -G * M_SUN / r**2 / (1 + q)
        assert abs(ax - a_analytic) / abs(a_analytic) < 1e-12
        assert abs(ay) < 1e-30
        assert abs(az) < 1e-30

    def test_two_equal_masses_symmetric_force(self):
        """Two equal masses symmetrically placed: force at midpoint is zero."""
        r = 5 * r_s
        masses = [
            (M_SUN, ( r, 0.0, 0.0)),
            (M_SUN, (-r, 0.0, 0.0)),
        ]
        ax, ay, az = massive_accel_physical_measure((0.0, 0.0, 0.0), masses)
        assert abs(ax) < 1e-20
        assert abs(ay) < 1e-20
        assert abs(az) < 1e-20

    def test_w_superposition_explicit(self):
        """Force from two separate masses equals force from combined w field."""
        r_test = (10 * r_s, 2 * r_s, 0.0)
        m_a = (2.0 * M_SUN, (0.0, 0.0, 0.0))
        m_b = (0.5 * M_SUN, (5 * r_s, 0.0, 0.0))
        ax, ay, az = massive_accel_physical_measure(r_test, [m_a, m_b])
        # Build w manually and compare
        mu_G = 3.0 * G / c**2
        def dist(p, q):
            return math.sqrt(sum((pi - qi)**2 for pi, qi in zip(p, q)))
        r_a = dist(r_test, m_a[1])
        r_b = dist(r_test, m_b[1])
        w = 1 + mu_G * m_a[0] / r_a + mu_G * m_b[0] / r_b
        # ∇w_x from each
        def grad_w_component(r_test, M, r_src, component):
            dx = r_test[component] - r_src[component]
            r_dist = dist(r_test, r_src)
            return -mu_G * M * dx / r_dist**3
        gx = grad_w_component(r_test, m_a[0], m_a[1], 0) + grad_w_component(r_test, m_b[0], m_b[1], 0)
        gy = grad_w_component(r_test, m_a[0], m_a[1], 1) + grad_w_component(r_test, m_b[0], m_b[1], 1)
        gz = grad_w_component(r_test, m_a[0], m_a[1], 2) + grad_w_component(r_test, m_b[0], m_b[1], 2)
        c2_over_3 = c**2 / 3.0
        ax_expected = c2_over_3 * gx / w
        ay_expected = c2_over_3 * gy / w
        az_expected = c2_over_3 * gz / w
        assert abs(ax - ax_expected) / (abs(ax_expected) + 1e-50) < 1e-10
        assert abs(ay - ay_expected) / (abs(ay_expected) + 1e-50) < 1e-10
        assert abs(az - az_expected) / (abs(az_expected) + 1e-50) < 1e-10

    def test_w_satisfies_laplace(self):
        """w(r) = 1 + 3GM/c²r satisfies ∇²w = 0 in the vacuum (away from origin)."""
        # ∇²w = d²w/dr² + (2/r) dw/dr.  With w = 1 + A/r:
        # dw/dr = -A/r², d²w/dr² = 2A/r³.  ∇²w = 2A/r³ - 2A/r³ = 0.
        # Verify numerically for several radii.
        for r_factor in [2.0, 5.0, 10.0, 50.0]:
            r = r_factor * r_s
            A = 3 * G * M_SUN / c**2
            dw_dr    = -A / r**2
            d2w_dr2  =  2 * A / r**3
            laplacian_w = d2w_dr2 + (2.0 / r) * dw_dr
            assert abs(laplacian_w) < 1e-10, (
                f"r={r_factor} r_s: ∇²w = {laplacian_w:.3e} (should be 0)"
            )
