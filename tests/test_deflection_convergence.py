import math

from pushing_medium.core import index_deflection_numeric, pm_deflection_angle_point_mass


def test_deflection_converges_steps(pm_calibration):
    M = 1.989e30
    b = 6.96e8
    mu = pm_calibration.mu_coeff
    z_max = 2e10
    coarse = index_deflection_numeric(M, b, mu, z_max=z_max, steps=500)
    mid = index_deflection_numeric(M, b, mu, z_max=z_max, steps=2000)
    fine = index_deflection_numeric(M, b, mu, z_max=z_max, steps=6000)
    diff_coarse_mid = abs(coarse - mid)
    diff_mid_fine = abs(mid - fine)
    assert diff_mid_fine < diff_coarse_mid
    analytic = pm_deflection_angle_point_mass(M, b)
    rel_err = abs(fine - analytic) / analytic
    assert rel_err < 0.03


def test_deflection_converges_zmax(pm_calibration):
    M = 1.989e30
    b = 6.96e8
    mu = pm_calibration.mu_coeff
    configs = [
        (5e9, 1500),
        (1e10, 3000),
        (2e10, 6000),
    ]
    values = [index_deflection_numeric(M, b, mu, z_max=Z, steps=st) for Z, st in configs]
    d1 = abs(values[1] - values[0])
    d2 = abs(values[2] - values[1])
    assert d2 < d1
    analytic = pm_deflection_angle_point_mass(M, b)
    rel_err_final = abs(values[-1] - analytic) / analytic
    assert rel_err_final < 0.03
