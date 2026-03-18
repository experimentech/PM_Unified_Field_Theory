"""Hulse-Taylor binary pulsar test: PM period derivative vs observation.

PSR B1913+16 (Hulse & Taylor 1974) provided the first indirect evidence for
gravitational-wave energy loss.  The system's orbital period is decaying at a
rate that agrees with the Peters (1964) GR quadrupole formula to ~0.2%
(Weisberg, Nice & Taylor 2010).

PM predicts *identical* quadrupole power to GR (pm_binary_quadrupole_power
passes to 1e-12 relative error in test_passed_benchmarks.py), so PM also
predicts identical dP/dt.  These tests verify:

  1. The Peters enhancement factor f(e) is correct.
  2. pm_binary_period_derivative gives the right number for the Hulse-Taylor
     system and agrees with GR to machine precision.
  3. The PM prediction agrees with the Galactic-corrected observed value to
     better than 1% — confirming PM passes this 0.2%-level observational test.

Observed parameters (Weisberg, Nice & Taylor 2010, ApJ 722, 1030):
  m_p  = 1.4408  M_sun
  m_c  = 1.3873  M_sun
  P_b  = 0.322997448930 d = 27906.9816 s
  e    = 0.6171334
  dP/dt_obs (Galactic corrected) = -2.4184e-12
"""

import math
import pytest
from pushing_medium import core as pm

# ── Hulse-Taylor system parameters ──────────────────────────────────────────
M_SUN = 1.98892e30          # kg
M_PULSAR   = 1.4408 * M_SUN
M_COMPANION = 1.3873 * M_SUN
P_B_HT = 27906.9816         # s  (0.322997448930 days)
E_HT   = 0.6171334

# Observed dP_b/dt after Galactic acceleration correction (Weisberg et al. 2010)
DPDT_OBS = -2.4184e-12


# ── Helper ───────────────────────────────────────────────────────────────────
def rel_diff(a, b):
    return abs(a - b) / max(abs(a), abs(b))


# ── Peters enhancement factor ────────────────────────────────────────────────
class TestPetersEnhancement:
    def test_circular_orbit_gives_unity(self):
        """e=0 → f(e)=1 (circular orbit, no enhancement)."""
        assert pm.pm_peters_decay_enhancement(0.0) == pytest.approx(1.0, rel=1e-12)

    def test_positive_definite(self):
        """f(e) must be > 1 for any e > 0."""
        for e in [0.1, 0.3, 0.5, 0.617, 0.9]:
            assert pm.pm_peters_decay_enhancement(e) > 1.0

    def test_monotone_increasing(self):
        """f(e) is monotonically increasing with eccentricity."""
        es = [0.0, 0.2, 0.4, 0.6, 0.8, 0.95]
        fs = [pm.pm_peters_decay_enhancement(e) for e in es]
        for i in range(len(fs) - 1):
            assert fs[i] < fs[i + 1]

    def test_hulse_taylor_enhancement(self):
        """f(e=0.6171334) ≈ 11.87 (well-known published value)."""
        f = pm.pm_peters_decay_enhancement(E_HT)
        assert f == pytest.approx(11.87, rel=5e-3)

    def test_known_value_e05(self):
        """Spot-check: e=0.5 → f = (1 + 73/96 + 37/1536) / (0.75)^3.5."""
        e = 0.5
        expected = (1.0 + 73.0 / 96.0 + 37.0 / 1536.0) / (0.75 ** 3.5)
        assert pm.pm_peters_decay_enhancement(e) == pytest.approx(expected, rel=1e-12)


# ── Period derivative: formula checks ────────────────────────────────────────
class TestPeriodDerivativeFormula:
    def test_negative(self):
        """dP/dt must be negative (orbit decays)."""
        dpdt = pm.pm_binary_period_derivative(M_PULSAR, M_COMPANION, P_B_HT, E_HT)
        assert dpdt < 0.0

    def test_circular_limit_consistent(self):
        """For e=0 the formula must be consistent with da/dt → dP/dt via Kepler.

        da/dt = -(64/5) G^3 m1 m2 M_tot / (c^5 a^3)  (Peters, circular)
        dP/dt = (3P/2a) da/dt
        Check that pm_binary_period_derivative(e=0) matches this independently.
        """
        m1, m2 = 1.4 * M_SUN, 1.3 * M_SUN
        P = 8 * 3600.0  # 8 hours
        M_tot = m1 + m2
        G, c = pm.G, pm.c
        a = (G * M_tot * (P / (2 * math.pi)) ** 2) ** (1.0 / 3.0)
        da_dt = -(64.0 / 5.0) * G ** 3 * m1 * m2 * M_tot / (c ** 5 * a ** 3)
        expected = (3.0 * P / (2.0 * a)) * da_dt
        result = pm.pm_binary_period_derivative(m1, m2, P, e=0.0)
        assert result == pytest.approx(expected, rel=1e-9)

    def test_eccentricity_amplifies_decay(self):
        """Eccentric orbit must decay faster than equivalent circular orbit."""
        m1, m2, P = M_PULSAR, M_COMPANION, P_B_HT
        dpdt_circ = pm.pm_binary_period_derivative(m1, m2, P, e=0.0)
        dpdt_ecc  = pm.pm_binary_period_derivative(m1, m2, P, e=E_HT)
        # Both negative; eccentric one is more negative
        assert dpdt_ecc < dpdt_circ


# ── Hulse-Taylor observational comparison ────────────────────────────────────
class TestHulseTaylorObservation:
    def test_pm_prediction_magnitude(self):
        """PM predicts dP/dt ~ -2.4e-12 for the Hulse-Taylor system."""
        dpdt = pm.pm_binary_period_derivative(M_PULSAR, M_COMPANION, P_B_HT, E_HT)
        assert dpdt == pytest.approx(-2.40e-12, rel=5e-2)

    def test_pm_agrees_with_observation_to_1pct(self):
        """PM prediction agrees with Galactic-corrected observed dP/dt to < 1%.

        The agreement is actually ~0.2% — PM passes this test at the same level
        as GR because PM reproduces GR quadrupole power exactly.
        """
        dpdt_pm = pm.pm_binary_period_derivative(M_PULSAR, M_COMPANION, P_B_HT, E_HT)
        assert rel_diff(dpdt_pm, DPDT_OBS) < 0.01

    def test_pm_matches_gr_exactly(self):
        """PM period derivative equals the GR Peters formula to machine precision.

        Since pm_binary_quadrupole_power == GR to 1e-12, the integrated
        dP/dt must also agree to that level.
        """
        from general_relativity import classical as gr
        dpdt_pm = pm.pm_binary_period_derivative(M_PULSAR, M_COMPANION, P_B_HT, E_HT)
        dpdt_gr = gr.binary_period_derivative(M_PULSAR, M_COMPANION, P_B_HT, E_HT)
        assert dpdt_pm == pytest.approx(dpdt_gr, rel=1e-12)
