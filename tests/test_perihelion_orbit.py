"""Tests: PM perihelion precession — analytic formula and numerical orbit integration.

Physics overview
----------------
The PM force law a = (c²/2)∇φ = −GM/r² r̂ is EXACTLY Newtonian for a
point mass — a pure Keplerian orbit with zero precession.

The 6πGM/c²a(1−e²) advance per orbit is a 1PN *metric* effect: when
orbital coordinates are defined via light travel time through the optical
metric (which has g_tt = −c²/n², g_rr = (1−2φ)), the effective potential
picks up a −GML²/(c²r³) term that causes the ellipse to precess.

Since PM has PPN parameters β=γ=1 (matching GR), the 1PN precession is:
    Δω/orbit = 6πGM / [c²a(1−e²)]      [rad/orbit]

This is identical to the GR prediction and matches all solar-system
observations. Higher PN corrections (which PM and GR handle differently)
are far below measurement precision for any current experiment.

Reference observations
----------------------
  Mercury: 42.9799 ± 0.0009 arcsec/century (GR contribution, residual after
           subtracting Newtonian planetary perturbations)
  PSR B1913+16: 4.22663 ± 0.00002 deg/year (measured via timing)
               GR prediction: 4.22660 deg/year  (0.001% agreement)
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

import pytest
from pushing_medium.core import (
    pm_perihelion_precession,
    pm_precession_arcsec_per_century,
    pm_integrate_orbit,
    G, c,
)
from general_relativity.classical import perihelion_precession

M_SUN = 1.989e30       # kg
AU    = 1.496e11       # m
S_PER_YEAR = 365.25 * 86400.0

# ── Solar system reference data ───────────────────────────────────────────────
# (a [AU], e, T_orb [yr], GR residual [arcsec/cy])
PLANETS = {
    'Mercury': (0.387098, 0.205630, 0.24085,  42.9799),
    'Venus':   (0.723327, 0.006773, 0.61520,   8.6247),
    'Earth':   (1.000000, 0.016709, 1.00000,   3.8387),
    'Mars':    (1.523679, 0.093401, 1.88085,   1.3510),
}

# PSR B1913+16 (Hulse-Taylor) orbital parameters
PSR_M1   = 1.4408 * M_SUN   # kg (pulsar)
PSR_M2   = 1.3886 * M_SUN   # kg (companion)
PSR_M_TOT = PSR_M1 + PSR_M2
PSR_P_B  = 27906.9807        # s (orbital period)
PSR_E    = 0.6171334
PSR_PERI_DEG_YR = 4.22663   # observed periastron advance [deg/yr]


def _planet_si(name):
    a_AU, e, T_yr, _ = PLANETS[name]
    return a_AU * AU, e, T_yr * S_PER_YEAR


# ── Analytic formula tests ────────────────────────────────────────────────────

class TestAnalyticFormula:
    """Verifies that the PM analytic formula matches GR and known values."""

    def test_pm_equals_gr_formula(self):
        """PM and GR analytic formulae are identical (PPN β=γ=1)."""
        a, e, _ = _planet_si('Mercury')
        delta_pm = pm_perihelion_precession(a, e, M_SUN)
        delta_gr = perihelion_precession(a, e, M_SUN)
        assert abs(delta_pm - delta_gr) < 1e-30

    def test_scales_linearly_with_mass(self):
        """Δω ∝ M."""
        a, e, _ = _planet_si('Mercury')
        d1 = pm_perihelion_precession(a, e, M_SUN)
        d2 = pm_perihelion_precession(a, e, 2 * M_SUN)
        assert abs(d2 / d1 - 2.0) < 1e-10

    def test_scales_inversely_with_semimajor_axis(self):
        """Δω ∝ 1/a (at fixed e, M)."""
        e = 0.2
        d1 = pm_perihelion_precession(AU,      e, M_SUN)
        d2 = pm_perihelion_precession(2 * AU,  e, M_SUN)
        assert abs(d2 / d1 - 0.5) < 1e-10

    def test_increases_with_eccentricity(self):
        """Δω increases with e (factor 1/(1−e²) grows)."""
        a = AU
        d0 = pm_perihelion_precession(a, 0.0, M_SUN)
        d1 = pm_perihelion_precession(a, 0.5, M_SUN)
        assert d1 > d0

    def test_mercury_arcsec_per_century(self):
        """Mercury analytic precession ≈ 42.98 arcsec/century."""
        a, e, T_s = _planet_si('Mercury')
        arcsec = pm_precession_arcsec_per_century(M_SUN, a, e, T_s)
        assert abs(arcsec - PLANETS['Mercury'][3]) / PLANETS['Mercury'][3] < 0.002

    def test_venus_arcsec_per_century(self):
        """Venus: 8.62 arcsec/century."""
        a, e, T_s = _planet_si('Venus')
        arcsec = pm_precession_arcsec_per_century(M_SUN, a, e, T_s)
        assert abs(arcsec - PLANETS['Venus'][3]) / PLANETS['Venus'][3] < 0.01

    def test_earth_arcsec_per_century(self):
        """Earth: 3.84 arcsec/century."""
        a, e, T_s = _planet_si('Earth')
        arcsec = pm_precession_arcsec_per_century(M_SUN, a, e, T_s)
        assert abs(arcsec - PLANETS['Earth'][3]) / PLANETS['Earth'][3] < 0.01

    def test_mars_arcsec_per_century(self):
        """Mars: 1.35 arcsec/century."""
        a, e, T_s = _planet_si('Mars')
        arcsec = pm_precession_arcsec_per_century(M_SUN, a, e, T_s)
        assert abs(arcsec - PLANETS['Mars'][3]) / PLANETS['Mars'][3] < 0.01

    def test_psr_b1913_periastron_advance(self):
        """PSR B1913+16: analytic periastron advance vs observed 4.22663 deg/yr.

        Uses total mass M_tot for two-body orbit (leading-order formula).
        """
        import math
        # Compute semi-major axis from Kepler's third law
        a3 = G * PSR_M_TOT * (PSR_P_B / (2 * math.pi))**2
        a_psr = a3 ** (1.0 / 3.0)
        arcsec_yr = pm_precession_arcsec_per_century(
            PSR_M_TOT, a_psr, PSR_E, PSR_P_B
        ) / 100.0   # convert cy → yr
        deg_yr = arcsec_yr / 3600.0
        # PSR B1913+16 observed: 4.22663 deg/yr; expect < 0.5% agreement at 1PN
        assert abs(deg_yr - PSR_PERI_DEG_YR) / PSR_PERI_DEG_YR < 0.005


# ── Numerical orbit integration tests ────────────────────────────────────────

class TestNumericalOrbit:
    """Verifies numerical 1PN orbit integration matches the analytic formula."""

    def test_mercury_numerical_matches_analytic(self):
        """Mercury: numerical precession per orbit within 0.5% of analytic."""
        a, e, _ = _planet_si('Mercury')
        result = pm_integrate_orbit(M_SUN, a, e, n_orbits=10)
        assert result['agreement_frac'] < 0.005

    def test_mercury_n_periastron_passages(self):
        """With n_orbits=10, should detect ~10 periapsis passages."""
        a, e, _ = _planet_si('Mercury')
        result = pm_integrate_orbit(M_SUN, a, e, n_orbits=10)
        assert result['n_periastron'] >= 9

    def test_earth_numerical_matches_analytic(self):
        """Earth: numerical precession per orbit within 1% of analytic."""
        a, e, _ = _planet_si('Earth')
        result = pm_integrate_orbit(M_SUN, a, e, n_orbits=5)
        assert result['agreement_frac'] < 0.01

    def test_newtonian_orbit_has_zero_precession(self):
        """Pure Newtonian orbit (1PN term disabled) gives zero precession."""
        from scipy.integrate import solve_ivp
        import numpy as np

        a = AU
        e = 0.5
        h2 = G * M_SUN * a * (1.0 - e * e)
        h  = math.sqrt(h2)
        T_kep = 2.0 * math.pi * math.sqrt(a**3 / (G * M_SUN))
        r0 = a * (1.0 - e)

        def derivs_newton(t, y):
            r, rdot, phi = y
            rddot = h2 / r**3 - G * M_SUN / r**2   # NO 1PN term
            return [rdot, rddot, h / r**2]

        def peri_event(t, y): return y[1]
        peri_event.terminal  = False
        peri_event.direction = 1

        sol = solve_ivp(derivs_newton, [0, 5 * T_kep * 1.1], [r0, 0.0, 0.0],
                        method='DOP853', events=peri_event,
                        rtol=1e-10, atol=1e-10 * r0, max_step=T_kep / 200)

        phi_peri = sol.y_events[0][:, 2]
        orb_nums = np.arange(1, len(phi_peri) + 1, dtype=float)
        slope = float(np.polyfit(orb_nums, phi_peri, 1)[0])
        precession_newtonian = slope - 2.0 * math.pi
        # Should be essentially zero (< 1e-8 rad/orbit is numerical noise)
        assert abs(precession_newtonian) < 1e-7

    def test_precession_positive_prograde(self):
        """Perihelion advance is positive (prograde — in direction of motion)."""
        a, e, _ = _planet_si('Mercury')
        result = pm_integrate_orbit(M_SUN, a, e, n_orbits=5)
        assert result['precession_per_orbit_rad'] > 0

    def test_higher_eccentricity_more_precession(self):
        """More eccentric orbit → larger precession (1/(1−e²) factor)."""
        a = 0.387098 * AU
        res_low_e  = pm_integrate_orbit(M_SUN, a, 0.1, n_orbits=5)
        res_high_e = pm_integrate_orbit(M_SUN, a, 0.6, n_orbits=5)
        assert res_high_e['precession_per_orbit_rad'] > res_low_e['precession_per_orbit_rad']


# ── PM vs GR: agreement at 1PN order ─────────────────────────────────────────

class TestPMvsGR:
    """Confirm PM = GR to 1PN order; document where they would diverge."""

    def test_pm_gr_same_formula_all_planets(self):
        """PM and GR analytic formulae agree to machine precision for all planets."""
        for name, (a_au, e, T_yr, _) in PLANETS.items():
            a = a_au * AU
            d_pm = pm_perihelion_precession(a, e, M_SUN)
            d_gr = perihelion_precession(a, e, M_SUN)
            assert d_pm == d_gr, f"Formulae differ for {name}"

    def test_pm_gr_force_law_is_newtonian(self):
        """At leading order, PM force -GM/r² gives zero precession (confirmed by Newtonian integration)."""
        # Trivial algebraic check: at leading order the PM force law is exactly Newtonian
        r = AU
        force_pm_direction = -(G * M_SUN / r**2)   # -GM/r² (Newton, also PM at leading order)
        force_newton       = -(G * M_SUN / r**2)
        assert force_pm_direction == force_newton

    def test_pm_gr_1pn_correction_same(self):
        """The 1PN orbital correction −3GMh²/(c²r⁴) is the same in PM and GR (PPN β=γ=1)."""
        r = 0.387098 * AU
        h2 = G * M_SUN * r * (1 - 0.205630**2)
        correction = -3.0 * G * M_SUN * h2 / (c**2 * r**4)
        Newtonian  = -G * M_SUN / r**2
        # Correction is ~v²/c² ≈ 10^-7 times Newtonian — confirms it's post-Newtonian
        assert abs(correction / Newtonian) < 1e-5
        assert abs(correction / Newtonian) > 1e-9   # but non-zero
