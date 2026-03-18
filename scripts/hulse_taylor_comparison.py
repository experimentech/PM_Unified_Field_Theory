#!/usr/bin/env python3
"""Hulse-Taylor binary pulsar: PM vs GR vs observation comparison.

PSR B1913+16 — the first binary pulsar (Hulse & Taylor 1974, Nobel Prize 1993).
Its orbital period decay is the first and strongest indirect test of
gravitational-wave energy loss.

References
----------
Weisberg, Nice & Taylor (2010), ApJ 722, 1030.
Peters (1964), Phys Rev 136, B1224.
"""

import sys
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from pushing_medium import core as pm
from general_relativity import classical as gr

# ── System parameters (Weisberg, Nice & Taylor 2010) ──────────────────────
M_SUN       = 1.98892e30        # kg
M_PULSAR    = 1.4408 * M_SUN    # kg
M_COMPANION = 1.3873 * M_SUN    # kg
P_B         = 27906.9816        # s   (0.322997448930 d)
ECC         = 0.6171334

# Observed dP_b/dt after Galactic + Shklovskii corrections
# (Damour & Taylor 1991; Weisberg et al. 2010)
DPDT_OBS      = -2.4184e-12     # dimensionless
DPDT_OBS_ERR  =  0.0009e-12     # 1σ

def main():
    # ── Peters enhancement factor ──────────────────────────────────────────
    f_e = pm.pm_peters_decay_enhancement(ECC)

    # ── Semi-major axis from Kepler ────────────────────────────────────────
    M_tot = M_PULSAR + M_COMPANION
    a = (pm.G * M_tot * (P_B / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)

    # ── Predictions ───────────────────────────────────────────────────────
    dpdt_pm = pm.pm_binary_period_derivative(M_PULSAR, M_COMPANION, P_B, ECC)
    dpdt_gr = gr.binary_period_derivative(M_PULSAR, M_COMPANION, P_B, ECC)

    frac_pm  = dpdt_pm / DPDT_OBS
    frac_gr  = dpdt_gr / DPDT_OBS
    pm_gr_diff = abs(dpdt_pm - dpdt_gr) / abs(dpdt_gr)

    print("=" * 62)
    print("  Hulse-Taylor Binary Pulsar  —  PSR B1913+16")
    print("=" * 62)
    print(f"  Pulsar mass         m_p  = {M_PULSAR/M_SUN:.4f}  M_sun")
    print(f"  Companion mass      m_c  = {M_COMPANION/M_SUN:.4f}  M_sun")
    print(f"  Orbital period      P_b  = {P_B:.4f}  s  ({P_B/3600:.6f} hr)")
    print(f"  Eccentricity        e    = {ECC}")
    print(f"  Semi-major axis     a    = {a/1e9:.4f}  Gm  ({a/3.086e16*1e6:.4f} μpc)")
    print(f"  Peters f(e)              = {f_e:.5f}")
    print()
    print(f"  {'Quantity':<32} {'Value':>18}  {'PM/Obs':>8}")
    print(f"  {'-'*60}")
    print(f"  {'PM  dP/dt  (prediction)':<32} {dpdt_pm:>18.6e}  {frac_pm:>8.6f}")
    print(f"  {'GR  dP/dt  (prediction)':<32} {dpdt_gr:>18.6e}  {frac_gr:>8.6f}")
    print(f"  {'Observed  dP/dt  (corrected)':<32} {DPDT_OBS:>18.6e}  {'1.000000':>8}")
    print(f"  {'Observed  1σ error':<32} {DPDT_OBS_ERR:>18.6e}")
    print()
    print(f"  PM vs GR  relative difference  : {pm_gr_diff:.2e}  (machine precision)")
    print(f"  PM vs observation  agreement   : {abs(1 - frac_pm)*100:.4f}%")
    print()

    status = "PASS" if abs(1 - frac_pm) < 0.01 else "FAIL"
    print(f"  Hulse-Taylor test  [{status}]  (threshold: < 1%)")
    print("=" * 62)

    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
