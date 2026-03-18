from pushing_medium import (
    index_point_masses,
    flow_rotational,
    ray_advection,
    massive_accel_medium,
    newtonian_accel_sum,
)


def test_index_point_masses_increases_with_mass():
    mu_scale = lambda M: 1.0 * M  # dimensionless scaling for test visibility
    n1 = index_point_masses((1, 0, 0), [(1.0, (0, 0, 0))], mu_scale=mu_scale)
    n2 = index_point_masses((1, 0, 0), [(2.0, (0, 0, 0))], mu_scale=mu_scale)
    assert n2 > n1


def test_flow_rotational_zero_far_from_zero_spin():
    u = flow_rotational((1, 2, 3), [])
    assert u == (0.0, 0.0, 0.0)


def test_ray_advection_dimensions():
    v = ray_advection((0, 0, 0), (1, 0, 0), n_total=1.0, u_g=(0, 0, 0))
    assert len(v) == 3 and v[0] > 0


def test_massive_accel_medium_direction():
    # Physical scenario: mass at origin, test particle at +x.
    # ∇φ at +x points TOWARD the mass = negative-x direction.
    # Corrected PM law: a = +(c²/2)∇φ → also negative-x (attractive). ✓
    a = massive_accel_medium((-1e-10, 0, 0))
    assert a[0] < 0


def test_massive_accel_newtonian_magnitude():
    # For a solar-mass point source, the PM field gives the Newtonian magnitude.
    from pushing_medium import G, c
    M_sun = 1.989e30  # kg
    r = 1.496e11       # 1 AU in metres
    mu_G = 2 * G / (c * c)
    # ∇φ at distance r: dφ/dr = -μ_G M/r²  (inward = negative radial)
    grad_phi_r = -mu_G * M_sun / r**2
    a_pm = massive_accel_medium((grad_phi_r, 0, 0))
    a_newton = -G * M_sun / r**2
    assert abs(a_pm[0] - a_newton) / abs(a_newton) < 1e-6, (
        f"PM massive force {a_pm[0]:.4e} should match Newtonian {a_newton:.4e}"
    )


def test_newtonian_accel_sum_matches_sign():
    a = newtonian_accel_sum((1, 0, 0), [(1.0, (0, 0, 0))])
    assert a[0] < 0
