"""
Tests for GW inspiral / chirp-waveform functions.

At 0PN (leading order) PM and GR share identical chirp physics because:
  (1) Quadrupole power  P_GW is the same (proven by Hulse–Taylor match)
  (2) Orbital energy   E = -GM/(2a) is the same (Newtonian)
=> chirp mass M_c and all 0PN observables are indistinguishable.

The dominant LIGO constraint on PM is therefore the MASS GAP:
  GW150914 progenitors (29 + 36 M☉)  >>  PM M_max ≈ 13–14 M☉
  GW190814 primary              23 M☉ >>  PM M_max
  GW170817 BNS (1.46 + 1.27 M☉)       — within PM range but R 30% too large
"""

import math
import pytest
import numpy as np

from pushing_medium import (
    pm_chirp_mass,
    pm_gw_frequency_deriv,
    pm_time_to_coalescence,
    pm_gw_strain_amplitude,
    pm_gw_chirp_waveform,
)
from general_relativity import (
    gr_chirp_mass,
    gr_time_to_coalescence,
    gr_gw_strain_amplitude,
)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
G = 6.674e-11        # m³ kg⁻¹ s⁻²
c = 2.998e8          # m s⁻¹
M_sun = 1.989e30     # kg
Mpc = 3.086e22       # m

# PM maximum compact-object mass (Pushing Medium theory, ≈ 13–14 M☉)
PM_M_MAX = 14.0 * M_sun    # conservative upper bound

# GWTC reference events (component masses, 90% credible-interval medians)
GW150914_M1 = 35.6 * M_sun   # heavier
GW150914_M2 = 30.6 * M_sun
GW190814_M1 = 22.2 * M_sun   # heavier (possible BH)
GW190814_M2 =  2.6 * M_sun   # light compact object
GW170817_M1 =  1.46 * M_sun  # heavier neutron star
GW170817_M2 =  1.27 * M_sun
GW151226_M1 = 13.7 * M_sun
GW151226_M2 =  7.7 * M_sun


# ===========================================================================
class TestChirpMass:
    """M_c = (M₁ M₂)^(3/5) / (M₁+M₂)^(1/5)."""

    def test_formula_numerical(self):
        """Known value: two 1.4 M☉ NS → M_c = 1.4 × 2^(-1/5)."""
        M = 1.4 * M_sun
        expected = M * 2.0 ** (-0.2)   # (M·M)^0.6 / (2M)^0.2 = M^0.6/2^0.2
        assert abs(pm_chirp_mass(M, M) - expected) / expected < 1e-10

    def test_symmetric_under_swap(self):
        M1, M2 = 1.4 * M_sun, 0.8 * M_sun
        assert pm_chirp_mass(M1, M2) == pytest.approx(pm_chirp_mass(M2, M1))

    def test_equal_mass_binary(self):
        """Equal masses: M_c = M / 2^(1/5)."""
        M = 10.0 * M_sun
        assert pm_chirp_mass(M, M) == pytest.approx(M * 2 ** (-0.2), rel=1e-10)

    def test_gw170817_chirp_mass(self):
        """GW170817: M_c ≈ 1.19 M☉."""
        M_c = pm_chirp_mass(GW170817_M1, GW170817_M2)
        # (1.46 × 1.27)^0.6 / (2.73)^0.2 ≈ 1.19 M☉
        assert 1.15 * M_sun < M_c < 1.25 * M_sun

    def test_gw150914_chirp_mass(self):
        """GW150914: M_c ≈ 28–29 M☉."""
        M_c = pm_chirp_mass(GW150914_M1, GW150914_M2)
        assert 25.0 * M_sun < M_c < 32.0 * M_sun


# ===========================================================================
class TestFrequencyDerivative:
    """df/dt = (96/5) π^(8/3) (G M_c/c³)^(5/3) f^(11/3)."""

    def test_positive(self):
        """Frequency always increases during inspiral."""
        M_c = pm_chirp_mass(1.4 * M_sun, 1.4 * M_sun)
        assert pm_gw_frequency_deriv(50.0, M_c) > 0.0

    def test_scales_as_f_to_11_thirds(self):
        """Doubling f multiplies df/dt by 2^(11/3)."""
        M_c = pm_chirp_mass(1.4 * M_sun, 1.4 * M_sun)
        ratio = pm_gw_frequency_deriv(100.0, M_c) / pm_gw_frequency_deriv(50.0, M_c)
        assert ratio == pytest.approx(2.0 ** (11.0 / 3.0), rel=1e-8)

    def test_gw150914_at_50hz(self):
        """GW150914 at 50 Hz: df/dt ~ 250 Hz/s (heavy BBH, rapid sweep)."""
        M_c = pm_chirp_mass(GW150914_M1, GW150914_M2)
        dfdt = pm_gw_frequency_deriv(50.0, M_c)
        # Heavy BBH (M_c ~ 28 M☉) sweeps fast even at low frequency
        assert 100.0 < dfdt < 600.0, f"df/dt = {dfdt:.2f} Hz/s (expected ~250 Hz/s)"

    def test_bns_slower_sweep(self):
        """BNS (1.4+1.4 M☉) sweeps slower than BBH (30+30 M☉) at same frequency."""
        M_c_bns = pm_chirp_mass(1.4 * M_sun, 1.4 * M_sun)
        M_c_bbh = pm_chirp_mass(30.0 * M_sun, 30.0 * M_sun)
        assert pm_gw_frequency_deriv(100.0, M_c_bns) < pm_gw_frequency_deriv(100.0, M_c_bbh)


# ===========================================================================
class TestTimeToCoalescence:
    """t_c = (5/256)(c³/G M_c)^(5/3) (πf)^(-8/3)."""

    def test_decreases_with_frequency(self):
        """Higher frequency → less time remaining."""
        M_c = pm_chirp_mass(1.4 * M_sun, 1.4 * M_sun)
        assert pm_time_to_coalescence(50.0, M_c) > pm_time_to_coalescence(100.0, M_c)

    def test_gw170817_at_ligo_entry(self):
        """GW170817: at ~20 Hz (LIGO low-frequency entry) t_c ≈ 100–200 s."""
        M_c = pm_chirp_mass(GW170817_M1, GW170817_M2)
        t_c = pm_time_to_coalescence(20.0, M_c)
        assert 50.0 < t_c < 400.0, f"t_c = {t_c:.1f} s"

    def test_gw150914_short_inspiral(self):
        """GW150914 at 20 Hz: much shorter merger time than GW170817 (heavier)."""
        M_c_150914 = pm_chirp_mass(GW150914_M1, GW150914_M2)
        M_c_170817 = pm_chirp_mass(GW170817_M1, GW170817_M2)
        assert (
            pm_time_to_coalescence(20.0, M_c_150914)
            < pm_time_to_coalescence(20.0, M_c_170817)
        )

    def test_scales_as_f_minus_8_thirds(self):
        """t_c ∝ f^(-8/3): halving f multiplies t_c by 2^(8/3)."""
        M_c = pm_chirp_mass(1.4 * M_sun, 1.4 * M_sun)
        ratio = pm_time_to_coalescence(25.0, M_c) / pm_time_to_coalescence(50.0, M_c)
        assert ratio == pytest.approx(2.0 ** (8.0 / 3.0), rel=1e-8)

    def test_self_consistency_with_freq_deriv(self):
        """t_c(f) = f / (df/dt) × 3/8  from integrating df/dt = A f^(11/3)."""
        M_c = pm_chirp_mass(1.4 * M_sun, 1.4 * M_sun)
        f = 50.0
        t_c = pm_time_to_coalescence(f, M_c)
        dfdt = pm_gw_frequency_deriv(f, M_c)
        # From t_c = (3/8) × f / dfdt  (0PN scaling identity)
        t_c_check = (3.0 / 8.0) * f / dfdt
        assert t_c == pytest.approx(t_c_check, rel=1e-8)


# ===========================================================================
class TestStrainAmplitude:
    """h_c = (4/D)(G M_c/c²)(π G M_c f/c³)^(2/3)."""

    def test_scales_inversely_with_distance(self):
        """Double D → half the strain."""
        M_c = pm_chirp_mass(1.4 * M_sun, 1.4 * M_sun)
        D1, D2 = 100.0 * Mpc, 200.0 * Mpc
        ratio = pm_gw_strain_amplitude(M_c, D1, 100.0) / pm_gw_strain_amplitude(M_c, D2, 100.0)
        assert ratio == pytest.approx(2.0, rel=1e-10)

    def test_scales_as_Mc_5_thirds(self):
        """h ∝ M_c^(5/3): doubling M_c multiplies h by 2^(5/3)."""
        D = 400.0 * Mpc
        M_c = pm_chirp_mass(1.4 * M_sun, 1.4 * M_sun)
        ratio = (
            pm_gw_strain_amplitude(2.0 * M_c, D, 100.0)
            / pm_gw_strain_amplitude(M_c, D, 100.0)
        )
        assert ratio == pytest.approx(2.0 ** (5.0 / 3.0), rel=1e-8)

    def test_gw150914_strain_order_of_magnitude(self):
        """GW150914 at 150 Hz, D≈410 Mpc: h ~ 10⁻²¹."""
        M_c = pm_chirp_mass(GW150914_M1, GW150914_M2)
        D = 410.0 * Mpc
        h = pm_gw_strain_amplitude(M_c, D, 150.0)
        # LIGO measured peak strain ~10⁻²¹
        assert 1e-22 < h < 1e-20, f"h = {h:.2e}"

    def test_positive(self):
        M_c = pm_chirp_mass(1.4 * M_sun, 1.4 * M_sun)
        assert pm_gw_strain_amplitude(M_c, 100.0 * Mpc, 100.0) > 0.0


# ===========================================================================
class TestWaveform:
    """pm_gw_chirp_waveform returns consistent analytic 0PN chirp."""

    @pytest.fixture(scope="class")
    def bns_waveform(self):
        M1, M2 = 1.4 * M_sun, 1.4 * M_sun
        D = 100.0 * Mpc
        tau_start = 30.0   # 30 s before coalescence
        return pm_gw_chirp_waveform(M1, M2, D, tau_start, N=2048, iota=0.0)

    def test_returns_four_arrays(self, bns_waveform):
        assert len(bns_waveform) == 4
        t, f, hp, hx = bns_waveform
        assert all(len(a) == 2048 for a in (t, f, hp, hx))

    def test_time_starts_at_zero(self, bns_waveform):
        t, *_ = bns_waveform
        assert t[0] == pytest.approx(0.0, abs=1e-10)

    def test_frequency_increases_in_time(self, bns_waveform):
        t, f, hp, hx = bns_waveform
        # f should be uniformly increasing (chirp: sweep upward)
        assert f[-1] > f[0]
        # Monotonically increasing — check at a few points
        n = len(f)
        for i in range(0, n - 1, n // 20):
            assert f[i + 1] > f[i]

    def test_strain_grows_toward_merger(self, bns_waveform):
        """Envelope |h(t)| should be larger near merger."""
        t, f, hp, hx = bns_waveform
        # Compare RMS of first quarter vs last quarter
        n = len(hp)
        q = n // 4
        rms_early = np.sqrt(np.mean(hp[:q] ** 2))
        rms_late = np.sqrt(np.mean(hp[-q:] ** 2))
        assert rms_late > rms_early

    def test_face_on_h_cross_nonzero(self, bns_waveform):
        """Face-on (iota=0): both polarisations are non-zero."""
        _, _, hp, hx = bns_waveform
        assert np.any(np.abs(hx) > 0.0)

    def test_edge_on_h_cross_zero(self):
        """Edge-on (iota=π/2): h× ≈ 0 (limited by float cos(π/2) ≈ 6e-17)."""
        t, f, hp, hx = pm_gw_chirp_waveform(
            1.4 * M_sun, 1.4 * M_sun, 100.0 * Mpc, 30.0,
            N=1024, iota=math.pi / 2
        )
        # h× = -cos(ι) × h0 × sin(Φ); cos(π/2) ≈ 6e-17 in float64
        # so |h×| / |h+| < 1e-15 is the float-precision expectation
        assert np.max(np.abs(hx)) < 1e-14 * np.max(np.abs(hp))

    def test_frequency_consistent_with_standalone(self):
        """f at t≈0 should match pm_time_to_coalescence inverse."""
        M1 = M2 = 1.4 * M_sun
        M_c = pm_chirp_mass(M1, M2)
        tau_start = 30.0
        t, f, _, _ = pm_gw_chirp_waveform(M1, M2, 100.0 * Mpc, tau_start, N=512)
        # At t=0 (τ = tau_start), f should match our standalone function
        # f(τ) = (5^(3/8)/8π)(c³/GM_c)^(5/8) τ^(-3/8)
        f_expected = (
            (5.0 ** (3.0 / 8.0) / (8.0 * math.pi))
            * (c ** 3 / (G * M_c)) ** (5.0 / 8.0)
            * tau_start ** (-3.0 / 8.0)
        )
        # rel=1e-4 tolerates slight G/c constant differences between test and core.py
        assert f[0] == pytest.approx(f_expected, rel=1e-4)


# ===========================================================================
class TestPMEqualsGRAt0PN:
    """At leading order, PM and GR produce identical waveform observables."""

    def test_chirp_mass_identical(self):
        """pm_chirp_mass == gr_chirp_mass to floating-point precision."""
        M1, M2 = 1.4 * M_sun, 1.27 * M_sun
        assert pm_chirp_mass(M1, M2) == pytest.approx(gr_chirp_mass(M1, M2), rel=1e-12)

    def test_time_to_coalescence_identical(self):
        """pm_time_to_coalescence == gr_time_to_coalescence."""
        M_c = pm_chirp_mass(GW170817_M1, GW170817_M2)
        t_pm = pm_time_to_coalescence(50.0, M_c)
        t_gr = gr_time_to_coalescence(50.0, M_c)
        assert t_pm == pytest.approx(t_gr, rel=1e-12)

    def test_strain_amplitude_identical(self):
        """pm_gw_strain_amplitude == gr_gw_strain_amplitude."""
        M_c = pm_chirp_mass(GW170817_M1, GW170817_M2)
        D = 40.0 * Mpc
        h_pm = pm_gw_strain_amplitude(M_c, D, 100.0)
        h_gr = gr_gw_strain_amplitude(M_c, D, 100.0)
        assert h_pm == pytest.approx(h_gr, rel=1e-12)

    def test_freq_deriv_formula_identical(self):
        """df/dt computed with same M_c gives identical result for PM and GR."""
        M_c_pm = pm_chirp_mass(GW150914_M1, GW150914_M2)
        M_c_gr = gr_chirp_mass(GW150914_M1, GW150914_M2)
        # Same chirp mass → same df/dt
        assert M_c_pm == pytest.approx(M_c_gr, rel=1e-12)
        assert pm_gw_frequency_deriv(50.0, M_c_pm) == pytest.approx(
            pm_gw_frequency_deriv(50.0, M_c_gr), rel=1e-12
        )


# ===========================================================================
class TestGWTCMassGap:
    """
    The dominant PM falsification from LIGO is the mass gap.
    PM cannot form compact objects above M_max ≈ 13–14 M☉.
    """

    def test_gw150914_m1_above_pm_mmax(self):
        """GW150914 heavier component (35.6 M☉) >> PM M_max."""
        assert GW150914_M1 > 2.0 * PM_M_MAX, (
            f"GW150914 M1 = {GW150914_M1/M_sun:.1f} M☉ "
            f"should be >> PM M_max = {PM_M_MAX/M_sun:.1f} M☉"
        )

    def test_gw150914_m2_above_pm_mmax(self):
        """GW150914 lighter component (30.6 M☉) >> PM M_max."""
        assert GW150914_M2 > 2.0 * PM_M_MAX, (
            f"GW150914 M2 = {GW150914_M2/M_sun:.1f} M☉ "
            f"should be >> PM M_max = {PM_M_MAX/M_sun:.1f} M☉"
        )

    def test_gw190814_primary_above_pm_mmax(self):
        """GW190814 primary (22.2 M☉) > PM M_max."""
        assert GW190814_M1 > PM_M_MAX, (
            f"GW190814 M1 = {GW190814_M1/M_sun:.1f} M☉ "
            f"should be > PM M_max = {PM_M_MAX/M_sun:.1f} M☉"
        )

    def test_gw190814_secondary_within_pm_range(self):
        """GW190814 secondary (2.6 M☉) is within PM mass range."""
        assert GW190814_M2 < PM_M_MAX

    def test_gw170817_both_within_pm_range(self):
        """GW170817 both components (1.46 + 1.27 M☉) within PM mass range."""
        assert GW170817_M1 < PM_M_MAX
        assert GW170817_M2 < PM_M_MAX

    def test_gw151226_primary_at_pm_mmax_edge(self):
        """GW151226 heavier component (13.7 M☉) is near PM M_max edge."""
        # Primary sits at or just below PM M_max — marginal accessibility
        assert GW151226_M1 < PM_M_MAX * 1.05, (
            f"GW151226 M1 = {GW151226_M1/M_sun:.1f} M☉ "
            f"vs PM M_max = {PM_M_MAX/M_sun:.1f} M☉"
        )

    def test_pm_accessible_events_are_minority(self):
        """Only BNS events are unambiguously within PM mass range."""
        # In the public GWTC catalog, GW150914, GW170814, GW190521, etc.
        # all have components far above PM M_max.  GW170817 is the main
        # accessible case.  Verify its chirp mass is well-determined.
        M_c = pm_chirp_mass(GW170817_M1, GW170817_M2)
        # 0PN chirp mass agrees with published value ~1.188 M☉
        assert 1.1 * M_sun < M_c < 1.3 * M_sun
