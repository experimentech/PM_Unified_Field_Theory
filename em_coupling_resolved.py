#!/usr/bin/env python3
"""
Clean analysis of PM electromagnetic coupling after constraint resolution.

Shows that the unified current-coupling mechanism (S_u = A_M J_M + A_q J_q) is
structurally correct and observationally constrained, while the magnetic potential
coupling (u_B = kappa_B A) is ruled out.
"""

import sys
sys.path.insert(0, 'src')

from pushing_medium import (
    estimate_A_M_from_GR,
)
from pushing_medium.core import (
    derive_A_q_upper_bound_from_lab,
    current_coupling_comparison,
)


def main():
    print("="*70)
    print("PUSHING-MEDIUM ELECTROMAGNETIC COUPLING: RESOLVED STRUCTURE")
    print("="*70)
    print()
    
    print("The PM unified theory predicts flow fields are sourced by currents:")
    print()
    print("  ∇²u = -S_u = -A_M J_M - A_q J_q")
    print()
    print("where:")
    print("  J_M = ρ_M v  (mass current density)")
    print("  J_q          (electric current density)")
    print()
    print("-"*70)
    
    # Get comparison
    comparison = current_coupling_comparison()
    
    print("COUPLING MECHANISM STATUS:")
    print()
    
    print("1. MASS CURRENTS (Gravitational)")
    mc = comparison['mass_current']
    print(f"   Coefficient: {mc['coefficient']} = {mc['value']:.6e} m³/(kg·s)")
    print(f"   Formula: {mc['formula']}")
    print(f"   Source: {mc['source']}")
    print(f"   Status: {mc['status']}")
    print()
    
    print("2. CHARGE CURRENTS (Direct Coupling)")
    cc = comparison['charge_current_direct']
    print(f"   Coefficient: {cc['coefficient']} < {cc['upper_bound']:.6e} m³/(C·s)")
    print(f"   Formula: {cc['formula']}")
    print(f"   Source: {cc['source']}")
    print(f"   Status: {cc['status']}")
    print(f"   Ratio: A_q/A_M < {cc['ratio_to_mass']:.3e}")
    print()
    
    print("3. MAGNETIC POTENTIAL (Ruled Out)")
    mp = comparison['magnetic_potential']
    print(f"   Coefficient: {mp['coefficient']}")
    print(f"   Formula: {mp['formula']}")
    print(f"   Status: {mp['status']}")
    print(f"   Reason: {mp['reason']}")
    print(f"   Bound: {mp['bound']}")
    print()
    
    print("-"*70)
    print("PHYSICAL INTERPRETATION")
    print("-"*70)
    print()
    
    print("The PM medium couples to CURRENTS (sources), not POTENTIALS (derived fields).")
    print()
    print("This means:")
    print("  ✓ Mass currents J_M = ρ_M v  → flow field u_G")
    print("  ✓ Charge currents J_q        → flow field u_q")
    print("  ✗ Magnetic potential A       → NO direct coupling")
    print()
    print("Why this is correct:")
    print("  1. Avoids double-counting (A is already derived from J_q)")
    print("  2. Treats all currents on equal structural footing")
    print("  3. Respects 'no special cases' principle")
    print("  4. Consistent with all laboratory null results")
    print()
    
    print("-"*70)
    print("OBSERVATIONAL CONSTRAINTS")
    print("-"*70)
    print()
    
    # Detailed constraint from specific lab setup
    print("Laboratory constraint from superconducting solenoid:")
    lab = derive_A_q_upper_bound_from_lab(
        I=1e6,              # 1 MA current
        L=1.0,              # 1 m length
        r=0.1,              # 10 cm distance
        path_length=1.0,    # 1 m beam path
        detection_threshold=1e-10  # ~10⁻¹⁰ rad (optical interferometry limit)
    )
    
    print(f"  Setup: {lab['lab_config']['current']/1e6:.1f} MA, "
          f"{lab['lab_config']['length']:.1f} m, "
          f"r={lab['lab_config']['distance']*100:.1f} cm")
    print(f"  Detection threshold: {lab['lab_config']['threshold']:.3e} rad")
    print()
    print(f"  Upper bound: A_q < {lab['A_q_max']:.6e} m³/(C·s)")
    print(f"  Compared to A_M: A_q/A_M < {lab['A_q/A_M_max']:.3e}")
    print()
    print(f"  {lab['interpretation']}")
    print()
    
    # What this means for different scenarios
    print("Implications:")
    print(f"  - Charge current dragging is at least {lab['A_q/A_M_max']:.0e}× weaker than mass dragging")
    print(f"  - For practical purposes: A_q ≈ 0")
    print(f"  - But nonzero coupling not forbidden—just unmeasurably small")
    print()
    
    print("-"*70)
    print("KEY INSIGHTS FOR PM THEORY DEVELOPMENT")
    print("-"*70)
    print()
    
    print("1. STRUCTURAL CLEANLINESS")
    print("   The unified source structure S_u = A_M J_M + A_q J_q is correct.")
    print("   Both currents couple via Poisson equations at source level.")
    print()
    
    print("2. NO SPECIAL CASES")
    print("   Mass and charge both treated as currents, same mathematical form.")
    print("   No ad-hoc 'potential coupling' that would privilege EM over gravity.")
    print()
    
    print("3. OBSERVATIONAL HIERARCHY")
    print("   Nature tells us: A_M ≫ A_q (by at least 10¹⁶ in current/current basis).")
    print("   This isn't a fine-tuning—it's just how strongly the medium responds")
    print("   to different types of charge flow.")
    print()
    
    print("4. FUTURE TESTS")
    print("   While A_q ≈ 0 for now, precision experiments could probe:")
    print("   - High-precision atom interferometry near extreme currents")
    print("   - Pulsar timing near magnetar current sheets")
    print("   - Gravitational wave + EM correlations in neutron star mergers")
    print()
    
    print("5. THEORETICAL ROBUSTNESS")
    print("   By eliminating κ_B A coupling and constraining A_q, PM avoids")
    print("   the 'landmine' of predicting huge lab effects. The theory remains")
    print("   clean, minimal, and consistent with all observations.")
    print()
    

if __name__ == '__main__':
    main()
