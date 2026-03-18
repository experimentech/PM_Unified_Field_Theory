#!/usr/bin/env python3
"""
Analysis of current-induced dragging in the Pushing-Medium unified field theory.

This script explores the novel prediction that electric currents create flow fields
that drag light and particles, analogous to (but distinct from) gravitational
frame-dragging from mass currents.
"""

import sys
import math
sys.path.insert(0, 'src')

from pushing_medium import (
    G, c, mu_0, epsilon_0,
    estimate_A_M_from_GR,
    current_dragging_ratio,
    flow_translational_mass_current,
    flow_charge_current,
    flow_from_vector_potential,
    magnetic_vector_potential_wire,
)


def analyze_current_dragging():
    """Compare gravitational and electromagnetic dragging mechanisms."""
    
    print("="*70)
    print("PUSHING-MEDIUM UNIFIED FIELD THEORY")
    print("Current-Induced Dragging Analysis")
    print("="*70)
    print()
    
    # Extract GR-constrained value for A_M
    A_M = estimate_A_M_from_GR()
    print(f"Gravitomagnetic coupling from GR analogy:")
    print(f"  A_M = G/c² = {A_M:.6e} m³/(kg·s)")
    print()
    
    # Scenario: Compare dragging near Earth vs a strong electromagnet
    print("-"*70)
    print("SCENARIO 1: Rotating Earth vs Laboratory Solenoid")
    print("-"*70)
    print()
    
    # Earth parameters
    M_earth = 5.972e24  # kg
    R_earth = 6.371e6  # m
    omega_earth = 7.292e-5  # rad/s (1 day)
    v_earth = omega_earth * R_earth  # ~465 m/s at equator
    
    # Solenoid parameters
    I_solenoid = 1e6  # 1 MA (superconducting coil)
    L_solenoid = 1.0  # 1 meter length
    
    # Observation distance
    r_obs = 10.0  # 10 meters from both
    
    print(f"Earth rotation at equator:")
    print(f"  Mass: {M_earth:.3e} kg")
    print(f"  Surface velocity: {v_earth:.1f} m/s")
    print(f"  Observation distance: {r_obs:.1f} m")
    print()
    
    print(f"Laboratory solenoid:")
    print(f"  Current: {I_solenoid:.3e} A")
    print(f"  Length: {L_solenoid:.1f} m")
    print(f"  Observation distance: {r_obs:.1f} m")
    print()
    
    # Calculate dragging ratio
    ratios = current_dragging_ratio(I_solenoid, L_solenoid, M_earth, v_earth, r_obs)
    
    print("Flow field magnitudes:")
    print(f"  |u_G| (mass current): {ratios['u_G']:.6e} m/s")
    print(f"  |u_q| (charge current, direct): {ratios['u_q']:.6e} m/s")
    print(f"  |u_B| (via magnetic potential): {ratios['u_B']:.6e} m/s")
    print()
    
    print("Dragging ratios:")
    print(f"  u_q/u_G = {ratios['u_q/u_G']:.6e}")
    print(f"  u_B/u_G = {ratios['u_B/u_G']:.6e}")
    print()
    
    print("Coupling coefficients used:")
    print(f"  A_M = {ratios['A_M']:.6e} m³/(kg·s)")
    print(f"  A_q = {ratios['A_q']:.6e} m³/(C·s)")
    print(f"  kappa_B = {ratios['kappa_B']:.6e} m/s per (T·m²)")
    print()
    
    # Angular deflection estimate
    path_length = 20.0  # meters through interaction region
    
    deflection_G = (ratios['u_G'] / c) * path_length / r_obs
    deflection_B = (ratios['u_B'] / c) * path_length / r_obs
    
    print("Estimated angular deflection (order of magnitude):")
    print(f"  From Earth's rotation: {deflection_G:.3e} radians = {deflection_G * 206265:.3e} arcsec")
    print(f"  From solenoid: {deflection_B:.3e} radians = {deflection_B * 206265:.3e} arcsec")
    print()
    
    print("-"*70)
    print("SCENARIO 2: Astrophysical Magnetar")
    print("-"*70)
    print()
    
    # Magnetar parameters
    M_ns = 2.8e30  # 1.4 solar masses
    R_ns = 1.2e4  # 12 km
    B_magnetar = 1e11  # 10¹¹ Tesla at surface
    
    # Estimate current from surface B field
    # For dipole: B ~ μ₀ μ / (4π r³) where μ = I A is magnetic moment
    # At surface: B ~ μ₀ I R / (4π R³) ~ μ₀ I / (4π R²)
    I_eff = B_magnetar * 4 * math.pi * R_ns * R_ns / mu_0
    
    print(f"Magnetar parameters:")
    print(f"  Mass: {M_ns/1.989e30:.1f} solar masses")
    print(f"  Radius: {R_ns/1e3:.1f} km")
    print(f"  Surface B field: {B_magnetar:.3e} T")
    print(f"  Effective current: {I_eff:.3e} A")
    print()
    
    # At distance of one radius above surface
    r_obs_ns = 2 * R_ns
    
    # Spinning at ~1 Hz (typical young pulsar)
    omega_ns = 2 * math.pi * 1.0  # rad/s
    v_ns = omega_ns * R_ns  # ~75 km/s at surface
    
    ratios_ns = current_dragging_ratio(I_eff, R_ns, M_ns, v_ns, r_obs_ns)
    
    print(f"At r = {r_obs_ns/1e3:.1f} km from center:")
    print(f"  |u_G| (rotation): {ratios_ns['u_G']:.6e} m/s")
    print(f"  |u_B| (magnetic): {ratios_ns['u_B']:.6e} m/s")
    print(f"  Ratio u_B/u_G = {ratios_ns['u_B/u_G']:.6e}")
    print()
    
    if ratios_ns['u_B/u_G'] > 0.01:
        print("*** SIGNIFICANT INTERACTION: Magnetic dragging comparable to rotational! ***")
        print("This would modify pulsar timing and light propagation near magnetars.")
    else:
        print("Magnetic dragging still subdominant to rotation for magnetars.")
    print()
    
    print("-"*70)
    print("KEY INSIGHTS")
    print("-"*70)
    print()
    print("1. STRUCTURAL INSIGHT:")
    print("   PM predicts electric currents create flow fields u_q just as")
    print("   mass currents create u_G. Same Poisson equation, different sources.")
    print()
    print("2. GR CONSTRAINT:")
    print("   Frame-dragging tests constrain A_M ≈ G/c².")
    print("   This sets the reference scale for all flow-field couplings.")
    print()
    print("3. TESTABLE PREDICTION:")
    print("   The ratio A_q/A_M determines how much stronger charge currents")
    print("   drag compared to mass currents at equal current densities.")
    print()
    print("4. OBSERVATIONAL VENUES:")
    print("   - Magnetars: extreme B fields → large u_B")
    print("   - Lab tests: superconducting coils near precision gyroscopes")
    print("   - Plasma astrophysics: current sheets in accretion disks")
    print()
    print("5. FUNDAMENTAL DIFFERENCE FROM GR:")
    print("   GR: Only rotating mass-energy drags frames (via metric perturbation)")
    print("   PM: Any current (mass or charge) drags medium (via flow field)")
    print("   This is a qualitatively different physics, not just a quantitative tweak.")
    print()


def estimate_laboratory_test():
    """Estimate detectability of current-dragging in lab."""
    
    print("="*70)
    print("LABORATORY TEST FEASIBILITY")
    print("="*70)
    print()
    
    # Best-case lab scenario
    I = 1e6  # 1 MA superconducting
    L = 1.0  # 1 m
    r = 0.1  # 10 cm from coil center
    path_length = 1.0  # 1 m beam path
    
    A_M = estimate_A_M_from_GR()
    kappa_B = c
    
    # Magnetic field flow
    A_mag = mu_0 * I * L / (4 * math.pi * r)
    u_B = kappa_B * A_mag
    
    # Deflection angle
    theta = (u_B / c) * (path_length / r)
    
    print(f"Setup:")
    print(f"  Solenoid current: {I/1e6:.1f} MA")
    print(f"  Coil length: {L:.1f} m")
    print(f"  Beam distance from coil: {r*100:.1f} cm")
    print(f"  Path length through region: {path_length:.1f} m")
    print()
    
    print(f"Predictions:")
    print(f"  |A| (vector potential): {A_mag:.6e} T·m")
    print(f"  |u_B| (flow field): {u_B:.6e} m/s")
    print(f"  |u_B|/c: {u_B/c:.6e}")
    print()
    
    print(f"  Deflection angle: {theta:.6e} radians")
    print(f"                  = {theta * 206265:.6e} arcseconds")
    print(f"                  = {theta * 206265 * 1e6:.6e} microarcseconds")
    print()
    
    # Compare to precision limits
    print(f"Detection feasibility:")
    print(f"  Atom interferometry: ~10⁻¹¹ rad sensitivity")
    print(f"  Optical interferometry: ~10⁻⁹ rad sensitivity")
    print(f"  Required for detection: {theta:.3e} rad")
    
    if theta > 1e-11:
        print("  → Potentially detectable with atom interferometry!")
    elif theta > 1e-9:
        print("  → Potentially detectable with optical interferometry!")
    else:
        print("  → Below current experimental sensitivity")
    print()
    
    # What would be needed
    if theta < 1e-11:
        factor_needed = 1e-11 / theta
        print(f"To reach atom interferometry sensitivity:")
        print(f"  Need {factor_needed:.1f}× improvement in:")
        print(f"    - Higher current (MA → GA range)")
        print(f"    - Closer proximity (cm → mm)")
        print(f"    - Longer path length (m → km)")
        print(f"    - Or kappa_B is larger than dimensional estimate")
    print()


if __name__ == '__main__':
    analyze_current_dragging()
    print()
    estimate_laboratory_test()
