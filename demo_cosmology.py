
import math

# Physical constants
c = 299792458.0
H0_Planck = 67.4  # km/s/Mpc (Planck 2018)
H0_Shoes = 74.0   # km/s/Mpc (SH0ES 2019)
Mpc_to_km = 3.08567758e19

def print_header(title):
    print("\n" + "="*60)
    print(f" {title}")
    print("="*60)

def flrw_luminosity_distance(z, H0, Omega_M=0.315, Omega_L=0.685):
    """
    Standard LCDM Luminosity Distance.
    """
    # H(z) = H0 * sqrt(Omega_M * (1+z)^3 + Omega_L)
    # D_L = (1+z) * c * integral(dz / H(z))
    
    # Numerical integration
    dz = 0.001
    # Simple loop for integration
    integral = 0.0
    curr_z = 0.0
    while curr_z < z:
        H_z = H0 * math.sqrt(Omega_M * (1+curr_z)**3 + Omega_L)
        integral += dz / H_z
        curr_z += dz
    
    return (1+z) * c * integral

def pm_luminosity_distance(z, n0, n_decay_rate):
    """
    PM 'Refractive' Luminosity Distance.
    Hypothesis: The universe is static but the refractive index n(t) decays.
    """
    # Simple model for exponential decay index: D_L = (c/H) * z
    return (c/n_decay_rate) * z 

def demo_hubble_tension():
    print_header("PM Refractive Cosmology & Hubble Tension")
    
    print("Hypothesis: Cosmic Redshift is due to a decaying refractive index n(t).")
    print("           1 + z = n(t_emit) / n(t_now)\n")
    
    # Redshifts to test (Supernova range)
    # 20 points from 0.01 to 1.5
    z_values = []
    val = 0.01
    step = (1.5 - 0.01) / 19
    for _ in range(20):
        z_values.append(val)
        val += step
    
    # 1. Planck LCDM
    dist_planck = [flrw_luminosity_distance(z, H0_Planck)/Mpc_to_km for z in z_values]
    
    # 2. SH0ES LCDM
    dist_shoes = [flrw_luminosity_distance(z, H0_Shoes)/Mpc_to_km for z in z_values]
    
    # 3. PM Static (Linear D_L = c/H * z) - The "Empty Universe" / Decaying Index Model
    # Try to fit Planck H0 scale
    dist_pm = [(c/H0_Planck) * z / Mpc_to_km for z in z_values]
    
    print(f"{'Redshift z':<10} | {'D_L (Mpc) Planck':<20} | {'D_L (Mpc) PM (Linear)':<20} | {'Diff %':<10}")
    print("-" * 70)
    
    for i, z in enumerate(z_values):
        if i % 4 == 0: # Print every 4th
            print(f"{z:<10.2f} | {dist_planck[i]:<20.1f} | {dist_pm[i]:<20.1f} | {100*(dist_pm[i]-dist_planck[i])/dist_planck[i]:<10.1f}")
            
    print("-" * 70)
    print("\nObservation:")
    print("The simple PM model (exponential index decay) yields D_L ~ z.")
    print("This matches 'Empty Universe' cosmology.")
    print("LCDM D_L curves UPWARDS (fainter SNe) due to Dark Energy.")
    print("PM's linear prediction is too bright at high z compared to SNe data.")
    
    print("\nResolving the Hubble Tension?")
    print("The Hubble Tension is a mismatch between H0 (local) and H0 (CMB).")
    print("Local (SH0ES): 74 km/s/Mpc")
    print("CMB (Planck): 67.4 km/s/Mpc")
    print("\nIf PM is correct, H0 is just the local decay rate of the refractive index.")
    print("To distinguish, we'd need to show that n(t) decay is NOT uniform,")
    print("or that local density (Laniakea supercluster) affects the index decay rate.")
    
    print("\nConclusion: PM provides a mechanism for redshift (n-decay) without expansion.")
    print("However, the simplest constant-decay model creates a 'linear' Hubble diagram")
    print("which conflicts with high-z Supernovae (Dark Energy) data.")
    print("To work, PM would need a decay rate H(n) that changes over time.")

if __name__ == "__main__":
    demo_hubble_tension()
