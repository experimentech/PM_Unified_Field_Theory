"""
Gap Probe Comparison: Three-Way Structural Test of PM Strong-Field Aberration
=============================================================================

Background
----------
The PM formula sheet (§ Deformation Energy / Diagnostic) identifies three
structural gaps in the current strong-field formulation:

  Gap 1: U(φ) is absent from the Lagrangian / Poisson equation.
         The linear Poisson  ∇²φ = −(8πG/c²)ρ_M  has no self-coupling;
         φ_crit=1 is externally imposed rather than dynamically emergent.
         Repair: ∇²φ − U′(φ)/c² = −(8πG/c²)ρ_M  (Option a).

  Gap 2: U(φ) and the EOS are thermodynamically inconsistent.
         Thermodynamic pressure from U: P_thermo = ε₀(2φ − 2φ² + φ³/3)
         Stated EOS:                   P_EOS    = c²/2(ρ − ρ_nuc)
         Ratio at small φ: P_thermo/P_EOS ≈ 4.
         Repair: replace P_EOS with P_thermo (Option a, or Option c).

  Gap 3: U′(φ) does not appear in hydrostatic equilibrium.
         Adding it as an elastic pressure: P_total = P_EOS + U(φ)
         modifies the equilibrium without changing the Poisson equation
         (Option c).

This script integrates the PM stellar structure ODE under four variants:

  Variant 0 (baseline): current PM  — linear Poisson, P_EOS
  Variant 1 (Option a): nonlinear Poisson with U′(φ) self-coupling
                         + EOS replaced by thermodynamic P from U
  Variant 2 (Option b): linear Poisson + P_EOS, but φ_crit is a free
                         parameter;  sweep φ_crit to find what value
                         is required to reach 30 M☉
  Variant 3 (Option c): linear Poisson + P_total = P_EOS + U(φ)
                         (U as additional elastic pressure)

Outputs
-------
  results/gap_probe_comparison.png   — 4-panel figure
  results/gap_probe_summary.txt      — numerical summary table

Usage
-----
    python scripts/gap_probe_comparison.py
"""

import math
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from pushing_medium.stellar_structure import (
    RHO_NUC, RHO_CRIT, MU_G, M_SUN, pm_eos_pressure, pm_eos_density
)

G = 6.67430e-11
c = 299792458.0
EPS0 = RHO_NUC * c * c   # ε₀ = ρ_nuc c²

# ---------------------------------------------------------------------------
# EOS helpers for the variants
# ---------------------------------------------------------------------------

def u_phi(phi):
    """Deformation energy density: U(φ) = ε₀(φ² − φ³/3)."""
    return EPS0 * (phi**2 - phi**3 / 3.0)

def u_prime(phi):
    """U′(φ) = ε₀(2φ − φ²)."""
    return EPS0 * (2.0 * phi - phi**2)

def u_double_prime(phi):
    """U″(φ) = 2ε₀(1 − φ).  Zero at φ_crit = 1."""
    return 2.0 * EPS0 * (1.0 - phi)

def p_thermo_from_phi(phi):
    """
    Thermodynamic pressure derived from U(φ) as internal energy density:
        P_thermo = ρ² ∂(U/ρ)/∂ρ  with  ρ = ρ_nuc e^φ

    Expanding:
        U/ρ = (ε₀/ρ)(φ² − φ³/3) = (c²/e^φ)(φ² − φ³/3)
        ∂(U/ρ)/∂ρ = (1/ρ_nuc)e^{-φ} d/dφ [(c²)(φ² − φ³/3)e^{-φ}]

    Result (exact):
        P_thermo = ε₀(2φ − 2φ² + φ³/3) + O(φ⁴)

    Note at small φ:  P_thermo ≈ 2ε₀φ = 4 × P_EOS(φ)
    """
    return EPS0 * (2.0*phi - 2.0*phi**2 + phi**3/3.0)

def p_thermo_from_rho(rho):
    """P_thermo expressed in terms of ρ via φ = ln(ρ/ρ_nuc)."""
    phi = math.log(rho / RHO_NUC) if rho > RHO_NUC else 0.0
    return max(p_thermo_from_phi(phi), 0.0)

def rho_from_p_thermo(P):
    """
    Invert P_thermo(ρ) numerically via Newton's method.

    Valid only on the ascending branch of P_thermo, i.e.
    φ < 2−√2 ≈ 0.586  (ρ < 1.797 × RHO_NUC).  Above that density
    P_thermo decreases with ρ and the EOS is non-invertible.

    Uses a leading-order initial guess from P_thermo ≈ 2 ε₀ φ instead of
    the P_EOS inverse (which is off by a factor of ~4 due to Gap 2).
    """
    if P <= 0.0:
        return RHO_NUC
    # Better initial guess: from leading-order P_thermo ≈ 2 ε₀ φ
    # Prevents Newton iteration from jumping to the wrong branch
    phi_init = min(P / (2.0 * EPS0), 0.57)   # cap below P_thermo maximum
    rho = RHO_NUC * math.exp(phi_init)
    for _ in range(40):
        phi = math.log(rho / RHO_NUC)
        P_current = p_thermo_from_phi(phi)
        # dP_thermo/dρ = (dP_thermo/dφ)(dφ/dρ) = (1/ρ) × ε₀(2 − 4φ + φ²)
        dfdphi = EPS0 * (2.0 - 4.0*phi + phi**2)
        if abs(dfdphi) < EPS0 * 1e-8:   # too close to P_thermo maximum
            break
        dPdρ = dfdphi / rho
        drho = (P - P_current) / dPdρ
        rho += drho
        rho = max(rho, RHO_NUC)
        rho = min(rho, RHO_NUC * math.exp(0.57))   # stay on ascending branch
        if abs(drho) < rho * 1e-12:
            break
    return rho


# ---------------------------------------------------------------------------
# Generic ODE integrator — same structure as stellar_structure.py
# ---------------------------------------------------------------------------

def _solve_star_variant(rho_central, rho_crit_override, eos_pressure_fn,
                        eos_density_fn, poisson_extra_fn=None,
                        r_max=1e5, n_eval=5000, rtol=1e-9, atol=1e-6):
    """
    Integrate compact-star ODE with swappable EOS and optional Poisson coupling.

    Parameters
    ----------
    rho_central         : float   central density [kg/m³]
    rho_crit_override   : float   stability cap (replaced or relaxed)
    eos_pressure_fn     : callable  P(ρ) → Pa
    eos_density_fn      : callable  ρ(P) → kg/m³
    poisson_extra_fn    : callable or None
                          If not None: called as f(phi) and added to rhs of
                          Poisson as  d²φ/dr² correction:
                          dφ/dr_effective = −μ_G m/r² + correction.
                          NB: in the nonlinear-Poisson option (a), the correction
                          enters through a modified relationship between φ and ρ,
                          not literally as a separate ODE term at fixed m.
                          Here we implement it as a source correction:
                            dφ/dr += U′(φ) × (r/c²) × some geometric factor.
                          See notes below.
    """
    if rho_central > rho_crit_override * 1.0001:
        return None   # outside stability cap

    P_central   = eos_pressure_fn(rho_central)
    phi_central = math.log(rho_central / RHO_NUC)

    def rhs(r, y):
        m, P, phi = y
        rho = eos_density_fn(P)
        phi_eos = math.log(max(rho, RHO_NUC) / RHO_NUC)  # φ from EOS

        if r < 1.0:
            dm_dr = 4.0 * math.pi * r * r * rho
            return [dm_dr, 0.0, 0.0]

        dm_dr   = 4.0 * math.pi * r * r * rho
        dP_dr   = -G * m * rho / (r * r)
        dphi_dr = -MU_G * m / (r * r)

        # Option (a) — nonlinear Poisson correction:
        # The modified Poisson ∇²φ − U′(φ)/c² = source adds a restoring term.
        # In spherical symmetry, Gauss theorem gives:
        #   dφ/dr = −μ_G m/r² + μ_G × [U′(φ)/c²] × r/3
        #         = −μ_G m/r² + MU_G × ρ_nuc(2φ−φ²) × r/3
        # Note: the naive (r/3) U′(φ)/c² form is missing the MU_G = 2G/c²
        # factor and has wrong units (kg/m² instead of m⁻¹).
        if poisson_extra_fn is not None:
            # poisson_extra_fn returns U′(φ) in Pa; divide by c² → kg/m³
            correction = MU_G * poisson_extra_fn(phi_eos) / (c * c) * r / 3.0
            dphi_dr += correction

        return [dm_dr, dP_dr, dphi_dr]

    def surface_event(r, y):
        return y[1]
    surface_event.terminal = True
    surface_event.direction = -1

    r_start   = 1.0
    m_start   = (4.0/3.0) * math.pi * r_start**3 * rho_central
    r_eval    = np.linspace(r_start, r_max, n_eval)

    sol = solve_ivp(rhs, [r_start, r_max],
                    [m_start, P_central, phi_central],
                    method='DOP853', t_eval=r_eval,
                    events=surface_event, rtol=rtol, atol=atol)

    P_arr   = np.maximum(sol.y[1], 0.0)
    rho_arr = np.vectorize(eos_density_fn)(P_arr)
    surface_mask = P_arr > 0
    if surface_mask.any():
        i_surf = np.where(surface_mask)[0][-1]
        return sol.t[i_surf], sol.y[0][i_surf]   # R [m], M [kg]
    return None


def compute_mr_curve_variant(rho_crit_override, eos_pressure_fn, eos_density_fn,
                              poisson_extra_fn=None, n_points=50):
    """Return M [M☉] and R [km] arrays for a full M–R sweep."""
    rho_c_arr = np.linspace(RHO_NUC * 1.01, rho_crit_override, n_points)
    M_arr, R_arr = [], []
    for rho_c in rho_c_arr:
        result = _solve_star_variant(rho_c, rho_crit_override,
                                     eos_pressure_fn, eos_density_fn,
                                     poisson_extra_fn)
        if result is not None:
            R_m, M_kg = result
            M_arr.append(M_kg / M_SUN)
            R_arr.append(R_m / 1e3)
        else:
            M_arr.append(np.nan)
            R_arr.append(np.nan)
    return np.array(M_arr), np.array(R_arr)


# ---------------------------------------------------------------------------
# Variant definitions
# ---------------------------------------------------------------------------

def make_option_c_eos():
    """Option (c): P_total = P_EOS + U(φ)  as a callable of ρ."""
    def p_total(rho):
        phi = math.log(max(rho, RHO_NUC) / RHO_NUC)
        return max(pm_eos_pressure(rho) + u_phi(phi), 0.0)

    def rho_from_p_total(P_tot):
        """Invert P_total(ρ) numerically."""
        if P_tot <= 0.0:
            return RHO_NUC
        rho = pm_eos_density(P_tot)  # initial guess
        for _ in range(40):
            phi = math.log(max(rho, RHO_NUC) / RHO_NUC)
            P_cur = max(pm_eos_pressure(rho) + u_phi(phi), 0.0)
            # dP_total/dρ = c²/2 + U′(φ)/ρ
            dPdrho = c*c/2.0 + u_prime(phi) / rho
            if abs(dPdrho) < 1e-10:
                break
            drho = (P_tot - P_cur) / dPdrho
            rho += drho
            rho = max(rho, RHO_NUC)
            if abs(drho) < rho * 1e-12:
                break
        return rho

    return p_total, rho_from_p_total


# ---------------------------------------------------------------------------
# Equation-of-state comparison (Gap 2 quantification)
# ---------------------------------------------------------------------------

def compare_eos():
    """Return φ values and both pressures for plotting."""
    phi_arr = np.linspace(0.0, 1.0, 200)
    rho_arr = RHO_NUC * np.exp(phi_arr)
    P_eos   = np.array([pm_eos_pressure(rho) for rho in rho_arr])
    P_therm = np.array([p_thermo_from_phi(phi) for phi in phi_arr])
    return phi_arr, P_eos, P_therm


# ---------------------------------------------------------------------------
# Option (b): swept φ_crit
# ---------------------------------------------------------------------------

def sweep_phi_crit_for_max_mass(phi_crit_values):
    """For each φ_crit, compute maximum PM mass with baseline EOS."""
    results = []
    for phi_c in phi_crit_values:
        rho_crit_new = RHO_NUC * math.exp(phi_c)
        M_vals, R_vals = compute_mr_curve_variant(
            rho_crit_new, pm_eos_pressure, pm_eos_density,
            poisson_extra_fn=None, n_points=40
        )
        valid = M_vals[np.isfinite(M_vals)]
        M_max = float(np.nanmax(M_vals)) if len(valid) > 0 else np.nan
        results.append(M_max)
    return np.array(results)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Computing M–R curve variants (this may take ~30 s)...")
    print()

    # --- Variant 0: baseline ---
    M0, R0 = compute_mr_curve_variant(RHO_CRIT, pm_eos_pressure, pm_eos_density,
                                       n_points=50)
    M0_max = float(np.nanmax(M0))
    print(f"Variant 0 (baseline):           M_max = {M0_max:.2f} M☉")

    # --- Variant 1: Option (a) — nonlinear Poisson + thermodynamic EOS ---
    # Find φ_crit for the thermodynamic EOS: U″(φ) > 0 → φ < 1 (same)
    # but now the medium also has P_thermo which may cap at a different φ
    # (P_thermo has a maximum at φ = 1 + √(4/3) ≈ 2.15; P_crit is higher)
    rho_crit_thermo = RHO_NUC * math.exp(1.0)   # keep φ_crit = 1 for option (a)
    M1, R1 = compute_mr_curve_variant(rho_crit_thermo,
                                       p_thermo_from_rho, rho_from_p_thermo,
                                       poisson_extra_fn=u_prime, n_points=50)
    M1_max = float(np.nanmax(M1))
    print(f"Variant 1 (Option a, U in Poisson + P_thermo): M_max = {M1_max:.2f} M☉")

    # --- Variant 2: Option (b) — sweep φ_crit ---
    phi_crit_sweep = np.linspace(1.0, 4.5, 20)
    M_max_sweep = sweep_phi_crit_for_max_mass(phi_crit_sweep)
    # Find φ_crit required for M_max ≥ 30 M☉ (GW150914 lower component)
    above_30 = phi_crit_sweep[M_max_sweep >= 30.0]
    phi_crit_needed = float(above_30[0]) if len(above_30) > 0 else float('inf')
    print(f"Variant 2 (Option b, φ_crit sweep):  φ_crit needed for M_max ≥ 30 M☉ = "
          f"{phi_crit_needed:.2f}  (rho_crit = {RHO_NUC*math.exp(phi_crit_needed):.2e} kg/m³)")

    # M-R at the required φ_crit
    rho_crit_b = RHO_NUC * math.exp(phi_crit_needed) if math.isfinite(phi_crit_needed) else RHO_CRIT
    M2, R2 = compute_mr_curve_variant(rho_crit_b, pm_eos_pressure, pm_eos_density,
                                       n_points=50)
    M2_max = float(np.nanmax(M2))
    print(f"  M-R at that φ_crit:  M_max = {M2_max:.2f} M☉")

    # --- Variant 3: Option (c) — P_total = P_EOS + U(φ) ---
    p_c, rho_c = make_option_c_eos()
    # φ_crit is same as baseline (stability from U″ = 0 at φ=1)
    M3, R3 = compute_mr_curve_variant(RHO_CRIT, p_c, rho_c,
                                       poisson_extra_fn=None, n_points=50)
    M3_max = float(np.nanmax(M3))
    print(f"Variant 3 (Option c, P_EOS + U(φ)):    M_max = {M3_max:.2f} M☉")

    # --- EOS comparison (Gap 2 quantification) ---
    phi_arr, P_eos, P_therm = compare_eos()
    ratio_at_1 = P_therm[-1] / P_eos[-1] if P_eos[-1] > 0 else float('inf')
    ratio_at_half = P_therm[len(phi_arr)//2] / P_eos[len(phi_arr)//2]
    print()
    print("Gap 2 — EOS comparison:")
    print(f"  P_thermo / P_EOS at φ = 0.5:  {ratio_at_half:.3f}")
    print(f"  P_thermo / P_EOS at φ = 1.0:  {ratio_at_1:.3f}")

    print()
    print("="*60)
    print("SUMMARY TABLE")
    print("="*60)
    print(f"{'Variant':<40} {'M_max/M☉':>8} {'Direction'}")
    print("-"*60)
    print(f"{'Baseline (current PM)':<40} {M0_max:>8.2f}  (reference)")
    print(f"{'Option a (U in Poisson + P_thermo)':<40} {M1_max:>8.2f}  {'↑' if M1_max > M0_max else '↓'}")
    print(f"{'Option b (free φ_crit, needed for 30 M☉)':<40} {phi_crit_needed:>8.2f}  (φ_crit value)")
    print(f"{'Option c (P_EOS + U(φ) elastic)':<40} {M3_max:>8.2f}  {'↑' if M3_max > M0_max else '↓'}")
    print(f"{'LIGO GW150914 lower component':<40} {'30.6':>8}  (target)")
    print("="*60)

    # --- Save summary to file ---
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'gap_probe_summary.txt'), 'w') as f:
        f.write("Gap Probe Comparison — PM Strong-Field Variants\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Gap 2 (EOS inconsistency):\n")
        f.write(f"  P_thermo / P_EOS at φ=0.5:  {ratio_at_half:.4f}\n")
        f.write(f"  P_thermo / P_EOS at φ=1.0:  {ratio_at_1:.4f}\n\n")
        f.write(f"M_max results:\n")
        f.write(f"  Baseline:              {M0_max:.4f} M☉\n")
        f.write(f"  Option a:              {M1_max:.4f} M☉\n")
        f.write(f"  Option b φ_crit:       {phi_crit_needed:.4f}  (needed for M_max≥30 M☉)\n")
        f.write(f"  Option c:              {M3_max:.4f} M☉\n")
        f.write(f"\nLIGO constraint: GW150914 components 35.6 + 30.6 M☉\n")
        f.write(f"GW170817 components: 1.46 + 1.27 M☉ (within PM range)\n")

    # --- Figure ---
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle(
        "PM Strong-Field Gap Probes: Three Repair Options vs Baseline",
        fontsize=13, fontweight='bold', y=0.98
    )

    # Panel A: M-R curve comparison
    ax = axes[0, 0]
    ax.plot(R0[np.isfinite(R0)], M0[np.isfinite(M0)],
            'k-', lw=2.5, label=f"Baseline (current PM)  $M_{{\\rm max}}={M0_max:.1f}\\,M_\\odot$")
    ax.plot(R1[np.isfinite(R1)], M1[np.isfinite(M1)],
            'C0--', lw=2, label=f"Option a: $U$ in Poisson + $P_{{\\rm thermo}}$  $M_{{\\rm max}}={M1_max:.1f}$")
    ax.plot(R2[np.isfinite(R2)], M2[np.isfinite(M2)],
            'C1-.', lw=2, label=f"Option b: $\\phi_{{\\rm crit}}={phi_crit_needed:.1f}$ (for 30 M$_\\odot$)")
    ax.plot(R3[np.isfinite(R3)], M3[np.isfinite(M3)],
            'C3:', lw=2.5, label=f"Option c: $P_{{\\rm EOS}} + U(\\phi)$  $M_{{\\rm max}}={M3_max:.1f}$")
    ax.axhline(30.6, color='firebrick', ls=':', lw=1.5, alpha=0.7, label="GW150914 M2 = 30.6 M☉")
    ax.axhline(1.46, color='steelblue', ls=':', lw=1.5, alpha=0.7, label="GW170817 M1 = 1.46 M☉")
    ax.set_xlabel("Radius [km]", fontsize=10)
    ax.set_ylabel("Mass [$M_\\odot$]", fontsize=10)
    ax.set_title("(A) M–R Curves: Four Variants", fontsize=11)
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 50)
    ax.set_ylim(0, max(M2_max, M1_max, M0_max) * 1.2)

    # Panel B: Gap 2 — EOS comparison
    ax = axes[0, 1]
    ax.semilogy(phi_arr, P_eos / 1e34, 'k-', lw=2.5,
                label=r"$P_{\rm EOS} = c^2(\rho - \rho_{\rm nuc})/2$")
    ax.semilogy(phi_arr, P_therm / 1e34, 'C0--', lw=2,
                label=r"$P_{\rm thermo}$ from $U(\phi)$")
    ax.semilogy(phi_arr, (P_eos + np.array([u_phi(p) for p in phi_arr])) / 1e34,
                'C3:', lw=2, label=r"$P_{\rm EOS} + U(\phi)$ (Option c)")
    ax.axvline(1.0, color='gray', ls='--', lw=1.2, label=r"$\phi_{\rm crit}=1$")
    ax.set_xlabel(r"$\phi$", fontsize=10)
    ax.set_ylabel(r"Pressure [$\times 10^{34}$ Pa]", fontsize=10)
    ax.set_title("(B) Gap 2: EOS Inconsistency", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ratio_arr = np.where(P_eos > 0, P_therm / P_eos, np.nan)
    ax2 = ax.twinx()
    ax2.plot(phi_arr, ratio_arr, 'C2-', lw=1.5, alpha=0.6)
    ax2.set_ylabel(r"Ratio $P_{\rm thermo}/P_{\rm EOS}$", color='C2', fontsize=9)
    ax2.tick_params(axis='y', labelcolor='C2')
    ax2.axhline(1.0, color='C2', ls=':', lw=0.8, alpha=0.5)

    # Panel C: Option (b) - φ_crit sweep
    ax = axes[1, 0]
    ax.plot(phi_crit_sweep, M_max_sweep, 'C1-o', lw=2, ms=4)
    ax.axhline(30.6, color='firebrick', ls='--', lw=1.5, label="GW150914 M2 = 30.6 M☉")
    ax.axhline(14.0, color='gray', ls=':', lw=1.2, label="Current PM M_max ≈ 14 M☉")
    ax.axvline(1.0, color='black', ls=':', lw=1, label=r"Current $\phi_{\rm crit}=1$")
    if math.isfinite(phi_crit_needed):
        ax.axvline(phi_crit_needed, color='C1', ls='--', lw=1.5,
                   label=fr"Needed $\phi_{{crit}}={phi_crit_needed:.1f}$")
    ax.set_xlabel(r"$\phi_{\rm crit}$ (relaxed stability cap)", fontsize=10)
    ax.set_ylabel(r"$M_{\rm max}$ [$M_\odot$]", fontsize=10)
    ax.set_title(r"(C) Option b: $M_{\rm max}$ vs $\phi_{\rm crit}$", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel D: U(φ) and its derivatives
    ax = axes[1, 1]
    phi_full = np.linspace(0, 1.2, 200)
    U_vals   = np.array([u_phi(p) / EPS0 for p in phi_full])
    Up_vals  = np.array([u_prime(p) / EPS0 for p in phi_full])
    Upp_vals = np.array([u_double_prime(p) / EPS0 for p in phi_full])
    ax.plot(phi_full, U_vals,  'k-',  lw=2, label=r"$U(\phi)/\varepsilon_0$")
    ax.plot(phi_full, Up_vals, 'C0--', lw=2, label=r"$U'(\phi)/\varepsilon_0$")
    ax.plot(phi_full, Upp_vals,'C3:',  lw=2, label=r"$U''(\phi)/\varepsilon_0$ (stability)")
    ax.axhline(0.0, color='black', lw=0.8, ls=':')
    ax.axvline(1.0, color='gray', ls='--', lw=1.2, label=r"$\phi_{\rm crit}=1$ ($U''=0$)")
    ax.set_xlabel(r"$\phi$", fontsize=10)
    ax.set_ylabel(r"Normalised value ($/ \varepsilon_0$)", fontsize=10)
    ax.set_title(r"(D) Deformation Energy $U(\phi)$ Structure", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.text(0.62, 0.92,
            "Gap 1: $U'(\\phi)$ absent from\nPoisson equation\n"
            "Gap 2: $U$ ≠ thermodynamic\nfree energy of EOS",
            transform=ax.transAxes, fontsize=8, va='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.85))

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = os.path.join(out_dir, 'gap_probe_comparison.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to {os.path.abspath(out_path)}")
    plt.close()


if __name__ == "__main__":
    main()
