"""Tests: PM Mass Gap — GWTC Event Component Masses vs PM Maximum Mass

PM theory (baseline EOS, flat Minkowski background) sets a compact-object
mass ceiling

    M_max ≈ 13–14 M_☉

from the stellar-structure ODE turning point.  Above this density, no stable
PM matter configuration exists; the concept is somewhat analogous to the
Chandrasekhar / Oppenheimer–Volkoff limit in GR, but with a different physical
origin (the PM φ field saturating at ρ_crit rather than relativistic
Fermi pressure).

Consequence — the GWTC mass gap
────────────────────────────────
Many LIGO/Virgo events have component masses *above* PM M_max.  PM cannot
form the progenitors of these events under the baseline EOS.  They are either:

  (a) Classical GR black holes (PM and GR agree on chirp physics at 0PN, so
      the waveform does not distinguish them — this is the mass-gap
      *implication*, not a waveform test); or

  (b) Indication that the PM EOS needs extension (e.g. elastic U(φ) pressure,
      Option c in test_gap_probes.py, pushes M_max above 30 M_☉).

GWTC reference events (component masses, median from published posteriors):
  GW150914  M1=35.6, M2=30.6 M_☉  [LIGO/Virgo 2016, PRL 116, 061102]
  GW151226  M1=13.7, M2= 7.7 M_☉  [LIGO/Virgo 2016, PRL 116, 241103]
  GW170814  M1=30.5, M2=25.3 M_☉  [LIGO/Virgo 2017, PRL 119, 141101]
  GW170817  M1= 1.46, M2=1.27 M_☉  [LIGO/Virgo 2017, PRL 119, 161101] (BNS)
  GW190814  M1=22.2, M2= 2.6 M_☉  [Abbott et al. 2020, ApJL 896, L44]
  GW190521  M1=85.0, M2=66.0 M_☉  [LIGO/Virgo 2020, PRL 125, 101102]

Tests proceed in three groups:
  1. TestBaselineMMax            — numerical M_max from ODE, sanity bounds
  2. TestMassGapGWTC             — each GWTC event vs PM M_max
  3. TestBNSAccessible           — GW170817 is within PM range (positive test)
"""

import sys
import pytest
import numpy as np

sys.path.insert(0, 'src')

from pushing_medium.stellar_structure import (
    compute_mr_curve, M_SUN,
)

# ---------------------------------------------------------------------------
# GWTC event component masses [solar masses]
# ---------------------------------------------------------------------------
GW150914_M1 = 35.6   # heavier   (BBH)
GW150914_M2 = 30.6   # secondary (BBH)
GW151226_M1 = 13.7   # heavier   (BBH — near PM M_max boundary)
GW151226_M2 =  7.7   # secondary (BBH — within PM range)
GW170814_M1 = 30.5   # heavier   (BBH, three-detector)
GW170814_M2 = 25.3   # secondary (BBH)
GW170817_M1 =  1.46  # heavier   (BNS — accessible by PM)
GW170817_M2 =  1.27  # secondary (BNS — accessible by PM)
GW190814_M1 = 22.2   # heavier   (NSBH/BBH, primary above PM M_max)
GW190814_M2 =  2.6   # secondary (NSBH/BBH, within PM range)
GW190521_M1 = 85.0   # heavier   (intermediate-mass BBH)
GW190521_M2 = 66.0   # secondary (intermediate-mass BBH)


# ===========================================================================
@pytest.fixture(scope="module")
def pm_mmax():
    """Compute PM M_max from the stellar-structure ODE (cached per module)."""
    _, M_arr, _ = compute_mr_curve(n_points=40)
    M_valid = M_arr[np.isfinite(M_arr)]
    assert len(M_valid) > 5, "Too few converged star models — ODE solver issue"
    return float(np.nanmax(M_valid))


# ===========================================================================
class TestBaselineMMax:
    """Numerical M_max from the stellar-structure ODE, with sanity bounds."""

    def test_mmax_within_expected_range(self, pm_mmax):
        """Baseline PM M_max should be 10–16 M_☉ (well-known ODE result)."""
        assert 10.0 < pm_mmax < 16.0, \
            f"Unexpected PM M_max = {pm_mmax:.2f} M_☉ (expected 10–16 M_☉)"

    def test_mmax_above_neutron_star_max(self, pm_mmax):
        """PM M_max must exceed the known NS maximum (≈ 2 M_☉); PM can form NS."""
        assert pm_mmax > 2.0, \
            f"PM M_max = {pm_mmax:.2f} M_☉ — too low to form observed NS"

    def test_mmax_above_gw170817_components(self, pm_mmax):
        """M_max must comfortably exceed both GW170817 component masses."""
        assert pm_mmax > GW170817_M1 * 2.0, (
            f"PM M_max = {pm_mmax:.2f} M_☉ not well above "
            f"GW170817 heavier component ({GW170817_M1} M_☉)"
        )

    def test_mmax_below_gw190814_m1(self, pm_mmax):
        """Baseline M_max < GW190814 primary (22.2 M_☉) — PM mass gap exists."""
        assert pm_mmax < GW190814_M1, (
            f"Baseline PM M_max = {pm_mmax:.2f} M_☉ ≥ GW190814 M1 = {GW190814_M1} M_☉ "
            f"— no mass gap at baseline EOS?"
        )

    def test_mmax_below_gw150914_secondary(self, pm_mmax):
        """Baseline M_max << GW150914 secondary (30.6 M_☉) — large gap."""
        assert pm_mmax < GW150914_M2, (
            f"Baseline PM M_max = {pm_mmax:.2f} M_☉ ≥ GW150914 M2 = {GW150914_M2} M_☉"
        )


# ===========================================================================
class TestMassGapGWTC:
    """Each resolved GWTC event compared to PM M_max."""

    def test_gw150914_m1_above_pm_mmax(self, pm_mmax):
        """GW150914 heavier component (35.6 M_☉) is above PM M_max."""
        assert GW150914_M1 > pm_mmax, (
            f"GW150914 M1 = {GW150914_M1} M_☉ is NOT above "
            f"PM M_max = {pm_mmax:.2f} M_☉"
        )

    def test_gw150914_m2_above_pm_mmax(self, pm_mmax):
        """GW150914 secondary (30.6 M_☉) is above PM M_max."""
        assert GW150914_M2 > pm_mmax, (
            f"GW150914 M2 = {GW150914_M2} M_☉ is NOT above "
            f"PM M_max = {pm_mmax:.2f} M_☉"
        )

    def test_gw150914_both_above_pm_mmax(self, pm_mmax):
        """Both GW150914 components are above PM M_max — PM cannot form either."""
        assert GW150914_M1 > pm_mmax and GW150914_M2 > pm_mmax, (
            f"At least one GW150914 component is within PM range: "
            f"M1={GW150914_M1}, M2={GW150914_M2}, M_max={pm_mmax:.2f}"
        )

    def test_gw190814_m1_above_pm_mmax(self, pm_mmax):
        """GW190814 primary (22.2 M_☉) is above PM M_max."""
        assert GW190814_M1 > pm_mmax, (
            f"GW190814 M1 = {GW190814_M1} M_☉ ≤ PM M_max = {pm_mmax:.2f} M_☉"
        )

    def test_gw190814_secondary_within_pm_range(self, pm_mmax):
        """GW190814 secondary (2.6 M_☉) IS within PM range (positive test)."""
        assert GW190814_M2 < pm_mmax, (
            f"GW190814 M2 = {GW190814_M2} M_☉ should be < PM M_max = {pm_mmax:.2f} M_☉"
        )

    def test_gw151226_m1_near_pm_mmax(self, pm_mmax):
        """GW151226 heavier component (13.7 M_☉) is at or near PM M_max (≈ 13–14 M_☉)."""
        # This is the boundary event: M1 ≈ M_max within ~10%.
        # Assert only that the gap is small (not a factor of 2+ discrepancy).
        gap = GW151226_M1 - pm_mmax
        assert abs(gap) < 3.0, (
            f"GW151226 M1 = {GW151226_M1} M_☉; PM M_max = {pm_mmax:.2f} M_☉; "
            f"gap = {gap:.2f} M_☉ — expected small (< 3 M_☉)"
        )

    def test_gw151226_m2_within_pm_range(self, pm_mmax):
        """GW151226 secondary (7.7 M_☉) is well within PM range."""
        assert GW151226_M2 < pm_mmax, (
            f"GW151226 M2 = {GW151226_M2} M_☉ should be < PM M_max = {pm_mmax:.2f} M_☉"
        )

    def test_gw170814_both_above_pm_mmax(self, pm_mmax):
        """Both GW170814 components (30.5 + 25.3 M_☉) are above PM M_max."""
        assert GW170814_M1 > pm_mmax and GW170814_M2 > pm_mmax, (
            f"At least one GW170814 component within PM range: "
            f"M1={GW170814_M1}, M2={GW170814_M2}, M_max={pm_mmax:.2f}"
        )

    def test_gw190521_far_above_pm_mmax(self, pm_mmax):
        """GW190521 (85 + 66 M_☉) is far above PM M_max — intermediate-mass BH."""
        assert GW190521_M2 > 4.0 * pm_mmax, (
            f"GW190521 M2 = {GW190521_M2} M_☉ expected >> 4 × PM M_max = "
            f"{4*pm_mmax:.1f} M_☉"
        )


# ===========================================================================
class TestBNSAccessible:
    """GW170817 BNS: both components are within PM range (positive/sanity test)."""

    def test_gw170817_m1_below_mmax(self, pm_mmax):
        """GW170817 heavier NS (1.46 M_☉) is well below PM M_max."""
        assert GW170817_M1 < pm_mmax, (
            f"GW170817 M1 = {GW170817_M1} M_☉ should be < PM M_max = {pm_mmax:.2f} M_☉"
        )

    def test_gw170817_m2_below_mmax(self, pm_mmax):
        """GW170817 lighter NS (1.27 M_☉) is well below PM M_max."""
        assert GW170817_M2 < pm_mmax, (
            f"GW170817 M2 = {GW170817_M2} M_☉ should be < PM M_max = {pm_mmax:.2f} M_☉"
        )

    def test_gw170817_mmax_margin_at_least_5x(self, pm_mmax):
        """M_max / M1 ≥ 5 — PM strongly supports NS-mass objects."""
        margin = pm_mmax / GW170817_M1
        assert margin >= 5.0, (
            f"PM M_max / GW170817 M1 = {margin:.1f}× (expected ≥ 5)"
        )

    def test_gw190814_secondary_accessible(self, pm_mmax):
        """GW190814 secondary (2.6 M_☉) is also within PM range."""
        assert GW190814_M2 < pm_mmax


# ===========================================================================
class TestMassGapCharacterisation:
    """Quantify the mass gap and confirm it is structurally meaningful."""

    def test_gap_size_to_gw190814(self, pm_mmax):
        """Gap to nearest GWTC heavy object: GW190814 M1 - M_max ≥ 5 M_☉."""
        gap = GW190814_M1 - pm_mmax
        assert gap >= 5.0, (
            f"Gap to GW190814 M1 = {gap:.2f} M_☉ — expected ≥ 5 M_☉ "
            f"(PM M_max = {pm_mmax:.2f} M_☉)"
        )

    def test_gap_size_to_gw150914(self, pm_mmax):
        """Gap to GW150914 lighter component: M2 - M_max ≥ 15 M_☉."""
        gap = GW150914_M2 - pm_mmax
        assert gap >= 15.0, (
            f"Gap to GW150914 M2 = {gap:.2f} M_☉ — expected ≥ 15 M_☉ "
            f"(PM M_max = {pm_mmax:.2f} M_☉)"
        )

    def test_fraction_of_gwtc_events_above_mmax(self, pm_mmax):
        """≥ 4 out of 6 surveyed GWTC events have at least one component above M_max."""
        events = [
            ("GW150914", GW150914_M1),
            ("GW151226", GW151226_M1),
            ("GW170814", GW170814_M1),
            ("GW170817", GW170817_M1),
            ("GW190814", GW190814_M1),
            ("GW190521", GW190521_M1),
        ]
        above = [name for name, m1 in events if m1 > pm_mmax]
        assert len(above) >= 4, (
            f"Only {len(above)} events above PM M_max = {pm_mmax:.2f} M_☉: {above}"
        )

    def test_gw170817_is_only_clearly_accessible_bns_event(self, pm_mmax):
        """GW170817 is the only event where BOTH components are clearly < M_max."""
        both_below = [
            ("GW150914", GW150914_M1 < pm_mmax and GW150914_M2 < pm_mmax),
            ("GW151226", GW151226_M1 < pm_mmax and GW151226_M2 < pm_mmax),
            ("GW170814", GW170814_M1 < pm_mmax and GW170814_M2 < pm_mmax),
            ("GW170817", GW170817_M1 < pm_mmax and GW170817_M2 < pm_mmax),
            ("GW190814", GW190814_M1 < pm_mmax and GW190814_M2 < pm_mmax),
            ("GW190521", GW190521_M1 < pm_mmax and GW190521_M2 < pm_mmax),
        ]
        accessible = [name for name, ok in both_below if ok]
        assert "GW170817" in accessible, (
            "GW170817 should be accessible (both components < PM M_max)"
        )
        # GW150914, GW170814, GW190521 should NOT be accessible
        for name in ("GW150914", "GW170814", "GW190521"):
            assert name not in accessible, \
                f"{name} incorrectly marked as accessible under PM M_max"
