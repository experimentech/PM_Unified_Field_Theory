from pushing_medium.core import (
    estimate_ppn_gamma_beta,
    effective_isotropic_metric,
    index_point_masses,
)


def test_ppn_gamma_beta_near_unity(pm_calibration):
    M = 1.989e30
    r = 1.5e11  # ~1 AU
    mu = pm_calibration.mu_coeff
    gamma, beta, U = estimate_ppn_gamma_beta(M, r, mu)
    assert abs(gamma - 1.0) < 1e-6
    assert abs(beta - 1.0) < 1e-6
    assert U < 1e-6


def test_effective_metric_signatures(pm_calibration):
    mu_fn = lambda M: pm_calibration.mu_coeff * M
    r_point = (1.5e11, 0.0, 0.0)
    n_val = index_point_masses(r_point, [(1.989e30, (0.0, 0.0, 0.0))], mu_scale=mu_fn)
    g_tt, g_rr = effective_isotropic_metric(n_val)
    assert g_tt < 0 and g_rr > 0
    assert abs(g_tt + 1.0) < 0.02
    assert abs(g_rr - 1.0) < 0.02
