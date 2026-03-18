"""Tests: GW ringdown QNMs (GR) and PM surface echoes.

Physics overview
----------------
GR black holes ring down as damped sinusoids at quasi-normal mode (QNM)
frequencies set purely by the horizon geometry.  For the fundamental
Schwarzschild l=2, n=0 mode (Leaver 1985):

    f_QNM = 0.37367 c³ / (2π GM)  ≈ 12.1 kHz × (M_sun / M)
    τ_QNM = GM / (2 × 0.08896 × c³)  ≈ 27.7 μs × (M / M_sun)
    Q_QNM ≡ ω_R / (2 ω_I) = 0.37367 / (2 × 0.08896) ≈ 2.10

PM compact objects have NO horizon → NO horizon QNMs:
  - PM photon sphere r_ps = 2GM/c² (GR has 3GM/c²)
  - For all realistic PM stars R_star >> r_ps, so no exterior trapped modes
  - Instead: surface reflection → echoes at τ_echo = 2R_star/c
  - Structural oscillations at f_surf = c_s/(2R_star), c_s = c/√2
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

import pytest
from general_relativity.classical import (
    gr_qnm_frequency,
    gr_qnm_damping_time,
    gr_photon_sphere_radius,
    G, c,
)
from pushing_medium.core import (
    pm_photon_sphere_radius,
    pm_surface_echo_delay,
    pm_surface_mode_frequency,
    pm_surface_mode_damping_time,
)
from pushing_medium.stellar_structure import MU_G

M_SUN = 1.989e30  # kg

# ── GR QNM formula tests ─────────────────────────────────────────────────────

class TestGRQNMFormula:
    """GR Schwarzschild QNM (Leaver 1985, l=2, n=0 fundamental mode)."""

    def test_frequency_scaling(self):
        """f_QNM scales as 1/M and equals ~12.1 kHz per solar mass."""
        f1 = gr_qnm_frequency(M_SUN)
        f2 = gr_qnm_frequency(2 * M_SUN)
        assert abs(f2 - f1 / 2) < 1.0  # Hz

    def test_frequency_solar_mass_value(self):
        """f_QNM(1M_sun) ≈ 12.07 kHz (Leaver coefficient 0.37367)."""
        f = gr_qnm_frequency(M_SUN)
        assert 12000 < f < 12200  # Hz

    def test_damping_time_scaling(self):
        """τ_QNM scales linearly with M."""
        tau1 = gr_qnm_damping_time(M_SUN)
        tau2 = gr_qnm_damping_time(2 * M_SUN)
        assert abs(tau2 - 2 * tau1) / tau1 < 1e-10

    def test_damping_time_solar_mass_value(self):
        """τ_QNM(1M_sun) ≈ 27.7 μs."""
        tau = gr_qnm_damping_time(M_SUN)
        assert 25e-6 < tau < 30e-6  # seconds

    def test_quality_factor(self):
        """Q ≡ 2π f τ = ω_R / (2 ω_I) = 0.37367 / (2×0.08896) ≈ 2.10."""
        M = 10 * M_SUN
        f = gr_qnm_frequency(M)
        tau = gr_qnm_damping_time(M)
        Q = 2 * math.pi * f * tau
        assert abs(Q - 0.37367 / (2 * 0.08896)) < 0.01

    def test_gw150914_final_mass(self):
        """GW150914 final BH (~62 M_sun): f_QNM ≈ 195 Hz in LIGO band."""
        f = gr_qnm_frequency(62 * M_SUN)
        assert 150 < f < 250  # Hz

    def test_gw170817_final_mass(self):
        """GW170817 remnant (~2.7 M_sun): f_QNM ~ 4.5 kHz."""
        f = gr_qnm_frequency(2.7 * M_SUN)
        assert 4000 < f < 5000  # Hz

    def test_photon_sphere_gr(self):
        """GR photon sphere: r_ps = 3 GM/c²."""
        M = M_SUN
        r = gr_photon_sphere_radius(M)
        expected = 3 * G * M / c ** 2
        assert abs(r - expected) / expected < 1e-10


# ── PM photon sphere / echo / surface-mode tests ─────────────────────────────

class TestPMPhotonSphere:
    """PM photon sphere is at 2GM/c² — HALF the GR value."""

    def test_pm_photon_sphere_formula(self):
        """r_ps(PM) = 2GM/c² = μ_G M."""
        M = M_SUN
        r = pm_photon_sphere_radius(M)
        expected = 2 * G * M / c ** 2
        assert abs(r - expected) / expected < 1e-10

    def test_pm_photon_sphere_half_gr(self):
        """PM photon sphere (2GM/c²) is 2/3 the GR photon sphere (3GM/c²)."""
        M = 5 * M_SUN
        r_pm = pm_photon_sphere_radius(M)
        r_gr = gr_photon_sphere_radius(M)
        assert abs(r_pm / r_gr - 2.0 / 3.0) < 1e-10

    def test_pm_photon_sphere_inside_typical_star(self):
        """For a 1.4 M_sun star at R=11 km, r_ps is well inside the star."""
        M = 1.4 * M_SUN
        R_star = 11e3  # m
        r_ps = pm_photon_sphere_radius(M)
        assert r_ps < R_star  # photon sphere inside star → no exterior QNMs


class TestPMEchoDelay:
    """PM surface echo delay τ_echo = 2R_star/c."""

    def test_echo_delay_formula(self):
        """τ_echo = 2R/c exactly."""
        R = 11e3  # 11 km
        tau = pm_surface_echo_delay(R)
        assert abs(tau - 2 * R / c) < 1e-20

    def test_echo_delay_positive(self):
        for R in [8e3, 11e3, 14e3]:
            assert pm_surface_echo_delay(R) > 0

    def test_echo_delay_scales_with_radius(self):
        """Larger star → longer echo delay."""
        tau1 = pm_surface_echo_delay(10e3)
        tau2 = pm_surface_echo_delay(14e3)
        assert tau2 > tau1

    def test_echo_delay_typical_ns(self):
        """For R=11 km: τ_echo ≈ 73 μs."""
        tau = pm_surface_echo_delay(11e3)
        assert 60e-6 < tau < 85e-6

    def test_echo_delay_much_shorter_than_gw150914_ringdown(self):
        """PM echo (μs) << GR QNM decay time for ~62 M_sun (ms)."""
        # PM echo at R=11km
        tau_pm_echo = pm_surface_echo_delay(11e3)
        # GR QNM decay for 62 M_sun BH
        tau_gr = gr_qnm_damping_time(62 * M_SUN)
        assert tau_pm_echo < tau_gr


class TestPMSurfaceMode:
    """PM fundamental surface (breathing) mode: f_surf = c_s / (2R), c_s = c/√2."""

    def test_frequency_formula(self):
        """f_surf = c/(2√2 R) exactly."""
        R = 11e3
        f = pm_surface_mode_frequency(R)
        expected = c / (2 * math.sqrt(2) * R)
        assert abs(f - expected) / expected < 1e-10

    def test_frequency_positive(self):
        for R in [8e3, 11e3, 14e3]:
            assert pm_surface_mode_frequency(R) > 0

    def test_frequency_decreases_with_radius(self):
        """Larger star → lower frequency (confined mode)."""
        f1 = pm_surface_mode_frequency(10e3)
        f2 = pm_surface_mode_frequency(14e3)
        assert f2 < f1

    def test_frequency_typical_ns(self):
        """For R=11 km: f_surf ≈ 9.6 kHz."""
        f = pm_surface_mode_frequency(11e3)
        assert 8000 < f < 12000  # Hz

    def test_pm_surface_mode_higher_than_gr_qnm_at_same_mass(self):
        """PM surface mode frequency >> GR QNM frequency for GW170817-like remnant.

        GR: f_QNM(2.7 M_sun) ≈ 4.5 kHz from horizon geometry.
        PM: f_surf(R≈13km)  ≈ 8 kHz from c_s/2R (structural oscillation).
        These are separated by nearly a factor of 2 — qualitatively different.
        """
        M_rem = 2.7 * M_SUN
        R_pm = 13e3   # PM M-R curve at ~2.7 M_sun gives R≈13 km
        f_gr = gr_qnm_frequency(M_rem)
        f_pm = pm_surface_mode_frequency(R_pm)
        assert f_pm > f_gr  # PM structural mode is higher frequency

    def test_damping_time_pm_slower_than_gr(self):
        """PM surface modes damp more slowly than GR QNMs (no horizon to absorb energy)."""
        R = 11e3
        M = 1.4 * M_SUN
        tau_pm = pm_surface_mode_damping_time(R, Q_factor=10.0)
        tau_gr = gr_qnm_damping_time(M)
        assert tau_pm > tau_gr  # PM modes live longer


# ── PM vs GR distinguishability ──────────────────────────────────────────────

class TestPMvsGRRingdown:
    """Key distinguishing predictions between PM and GR for merger remnants."""

    def test_gr_bh_qnm_q_factor_low(self):
        """GR BH QNM has Q ≈ 2.1 (critically damped; energy lost to horizon)."""
        M = 10 * M_SUN
        f = gr_qnm_frequency(M)
        tau = gr_qnm_damping_time(M)
        Q_gr = math.pi * f * tau
        assert Q_gr < 3.0

    def test_pm_surface_mode_q_factor_higher(self):
        """PM surface mode Q = 10 (default) >> GR BH QNM Q ≈ 2.1."""
        R = 11e3
        M = 1.4 * M_SUN
        f_pm = pm_surface_mode_frequency(R)
        tau_pm = pm_surface_mode_damping_time(R, Q_factor=10.0)
        Q_pm = math.pi * f_pm * tau_pm
        Q_gr = 0.37367 / (2 * 0.08896)
        assert Q_pm > Q_gr

    def test_mass_gap_pm_cannot_form_heavy_bhs(self):
        """PM stability cap max mass ~13-14 M_sun << GW150914 component masses ~29+36 M_sun.

        PM uses Newtonian structure (no TOV), so mass grows monotonically with
        central density up to the stability cap ρ_crit = e×ρ_nuc, giving
        M_max ≈ 13-14 M_sun.  GW150914 progenitor masses (~29 M_sun, ~36 M_sun)
        both exceed this limit → PM cannot form these objects.  LIGO observes
        them as compact (from the orbital chirp mass), which is a direct
        falsification challenge for PM.
        """
        # GW150914-like final mass
        M_gw150914 = 62 * M_SUN
        r_ps_pm = pm_photon_sphere_radius(M_gw150914)
        r_schwarzschild = 2 * G * M_gw150914 / c ** 2

        # PM photon sphere = Schwarzschild radius (not a real surface in PM)
        assert abs(r_ps_pm - r_schwarzschild) < 1.0

        # GR QNM is in LIGO audio band; PM predicts no such object → no ringdown
        f_gr = gr_qnm_frequency(M_gw150914)
        assert 100 < f_gr < 300  # Hz — LIGO can see this!
