
import math
import sys
import os

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from pushing_medium.core import (
    unified_index_field, 
    unified_flow_field, 
    G, c, e, epsilon_0, m_e, mu_0
)

def print_header(title):
    print("\n" + "="*60)
    print(f" {title}")
    print("="*60)

def demo_charged_lensing():
    print_header("DEMO 1: Unified Lensing (Mass + Charge)")
    
    # Parameters for a hypothetical charged object
    mu_E_dim = e * e / (epsilon_0 * c * c) # ~ 3e-44 m/C
    print(f"Dimensional mu_E estimate: {mu_E_dim:.2e} m/C")
    
    # Compare Electron vs Proton vs Macroscopic
    scenarios = [
        ("Electron", m_e, -e),
        ("Proton", 1.67e-27, e),
        ("Charged Sphere (1kg, 1mC)", 1.0, 1e-3)
    ]
    
    r_obs = (1e-10, 0, 0) # Observer at 1 Angstrom for particles
    
    print(f"\nComparing Refractive Index Perturbations (delta_n) at r = {r_obs[0]:.1e} m")
    print("-" * 105)
    print(f"{'Object':<25} | {'Mass (kg)':<10} | {'Charge (C)':<10} | {'dn_G':<12} | {'dn_E':<12} | {'Ratio E/G':<10}")
    print("-" * 105)
    
    for name, M, Q in scenarios:
        # If macroscopic, move observer away
        dist = r_obs
        if M > 1e-5:
            dist = (1.0, 0, 0) # 1 meter
        
        # Calculate Gravitational part directly to avoid precision loss
        # n = 1 + mu M / r
        mu_G = 2 * G / (c * c)
        r = math.sqrt(dist[0]**2 + dist[1]**2 + dist[2]**2)
        dn_G = mu_G * M / r
        
        # Calculate Electric part directly
        # n = 1 + mu_E Q / (4 pi eps0 r)
        phi_E = abs(Q) / (4 * math.pi * epsilon_0 * r)
        dn_E = mu_E_dim * phi_E
        
        ratio = dn_E / dn_G if abs(dn_G) > 0 else float('inf')
        
        print(f"{name:<25} | {M:<10.2e} | {Q:<10.2e} | {dn_G:<12.2e} | {dn_E:<12.2e} | {ratio:<10.2e}")

    print("\nKey Insight: For elementary particles, the electric contribution to the")
    print("refractive index (if mu_E ~ classical radius) is comparable in magnitude to")
    print("the gravitational contribution, potentially unifying mass and charge as")
    print("similar defects in the medium.")

def demo_current_dragging():
    print_header("DEMO 2: Current-Induced Frame Dragging")
    
    print("PM predicts that electric currents drag spacetime, just like mass currents.")
    print("Equation: ∇²u = -A_M J_M - A_q J_q")
    
    # Current loop parameters
    I = 1e6  # 1 Mega-Ampere (extreme lab current)
    R_loop = 1.0 # 1 meter radius
    circumference = 2 * math.pi * R_loop
    
    # Rotating sphere of water (1000 kg/m^3), R=1m
    M_sphere = (4/3) * math.pi * (R_loop**3) * 1000
    v_surface = 100.0 # 100 m/s rotation
    
    r_obs = (1.1, 0, 0) # Just outside
    
    # Gravitational flow
    u_G_vec = unified_flow_field(r_obs, moving_masses=[(M_sphere, (0,0,0), (0,v_surface,0))])
    u_G_mag = math.sqrt(sum(x*x for x in u_G_vec))
    
    # Natural A_q (tiny)
    A_q_natural = e * e / (epsilon_0 * c * c) # ~3e-44
    
    # Maximum allowed A_q by lab constraints (<1e-11)
    A_q_bound_max = 1e-11 
    
    # Loop calculation
    segments = [
        ((R_loop, 0, 0), (0, R_loop, 0), I),
        ((0, R_loop, 0), (-R_loop, 0, 0), I),
        ((-R_loop, 0, 0), (0, -R_loop, 0), I),
        ((0, -R_loop, 0), (R_loop, 0, 0), I)
    ]
    
    u_q_natural_vec = unified_flow_field(r_obs, current_elements=segments, A_q=A_q_natural)
    u_q_natural_mag = math.sqrt(sum(x*x for x in u_q_natural_vec))
    
    u_q_max_vec = unified_flow_field(r_obs, current_elements=segments, A_q=A_q_bound_max)
    u_q_max_mag = math.sqrt(sum(x*x for x in u_q_max_vec))
    
    print(f"\nComparison at r = 1.1m:")
    print(f"Gravity Source: {M_sphere:.0f} kg sphere spinning at {v_surface} m/s")
    print(f"EM Source:      {I:.0e} A current loop (1m radius)")
    print("-" * 80)
    print(f"{'Field Type':<35} | {'Flow Velocity u (m/s)':<25}")
    print("-" * 80)
    print(f"{'Gravitational Drag (u_G)':<35} | {u_G_mag:.4e}")
    print(f"{'EM Drag (Natural A_q ~ 10^-44)':<35} | {u_q_natural_mag:.4e}  <-- Theoretical min")
    print(f"{'EM Drag (Max Allowed A_q ~ 10^-11)':<35} | {u_q_max_mag:.4e}  <-- Observation limit")
    
    print("\nInterpretation:")
    print("If A_q were as large as allowed by current limits (10^-11), a 1 MA current")
    print("would produce drag (10^-5 m/s) comparable to a spinning planet.")
    print("However, the 'natural' coupling (10^-44) suggests the effect is truly")
    print("negligible unless new physics enhances it.")

def demo_plasma_index():
    print_header("DEMO 3: Plasma Lensing Dispersions")
    
    print("PM predicts frequency-dependent lensing in plasma.")
    print("Formula: n_plasma = 1 - ω_p²/(2ω²)")
    
    # Interstellar medium density
    n_e = 1e6 # m^-3 (approx ISM)
    
    # Frequencies
    freqs = [
        ("Radio (100 MHz)", 100e6),
        ("Visible (500 THz)", 500e12),
        ("X-Ray (1 EHz)", 1e18)
    ]
    
    print(f"\nRefractive Index Deviation (1-n) for ISM density n_e={n_e:.0e} m^-3")
    print("-" * 50)
    
    for name, f in freqs:
        omega = 2 * math.pi * f
        n = unified_index_field((0,0,0), electron_density=n_e, alpha_plasma=None) 
        # Note: core.py has simplistic linear plasma. 
        # Let's use the explicit frequency dependent logic if available or calculate it
        
        # Re-using logic from core.py manually since the unified_index_field 
        # might use the linearized alpha_plasma default
        omega_p_sq = n_e * e * e / (epsilon_0 * m_e)
        delta_n = omega_p_sq / (2 * omega * omega)
        
        print(f"{name:<20} | {delta_n:.4e}")
        
    print("\nEffect:")
    print("Radio waves are lensed LESS strongly effectively (or bent differently)")
    print("than optical light in a PM plasma environment, providing a testable prediction.")

if __name__ == "__main__":
    print("PUSHING-MEDIUM: UNIQUE ELECTROMAGNETIC CAPABILITIES EXPLORATION")
    demo_charged_lensing()
    demo_current_dragging()
    demo_plasma_index()
