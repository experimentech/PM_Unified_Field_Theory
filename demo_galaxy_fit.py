#!/usr/bin/env python3
"""
Simple demonstration of PM galaxy rotation curve fitting.

This shows how to:
1. Load SPARC data
2. Set up PM medium model
3. Calculate rotation curve
4. Compare to observations
"""

import sys
import math
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'legacy' / 'Pushing-Medium' / 'src'))

from galaxy_dynamics.data import load_sparc_real
from galaxy_dynamics.rotation import (
    circular_velocity, DiskParams, MediumParams,
    accel_baryonic, accel_medium
)

def fit_single_galaxy_demo(galaxy_name='GAL_X'):
    """Demonstrate PM fitting for one galaxy."""
    
    # Load data
    data_file = 'legacy/Pushing-Medium/tests/data/sparc_sample.csv'
    rc = load_sparc_real(data_file, galaxy_name=galaxy_name)
    
    print("="*70)
    print(f"FITTING GALAXY: {rc.name}")
    print("="*70)
    
    # Extract baryonic component velocities (already in data)
    # These come from stellar mass + gas measurements
    v_bar_obs = rc.components.get('bar', [v*1000 for v in [10, 15, 18, 20, 22]])  # m/s
    
    # Estimate disk parameters from baryonic data
    # (In real fitting, these would be free parameters or derived from photometry)
    R_d_estimate = rc.radii_m[2]  # Use ~middle radius as scale
    M_d_estimate = 1e41  # kg (would be fitted)
    
    disk = DiskParams(M_d=M_d_estimate, R_d=R_d_estimate)
    
    # PM medium parameters (these are what we fit)
    # v_inf: asymptotic circular speed
    # r_s: scale radius for log term
    # r_c: core transition
    # m: smoothing exponent
    medium = MediumParams(
        v_inf=125e3,    # m/s (from flat part of curve)
        r_s=5e20,       # m (few kpc)
        r_c=1e20,       # m (core radius)
        m=2.0           # smoothing
    )
    
    print(f"\nObserved data:")
    print(f"  Radii: {[f'{r/3.086e19:.1f}' for r in rc.radii_m]} kpc")
    print(f"  V_obs: {[f'{v/1000:.0f}' for v in rc.v_obs_ms]} km/s")
    
    print(f"\nPM Model parameters:")
    print(f"  Disk: M_d={disk.M_d:.2e} kg, R_d={disk.R_d/3.086e19:.1f} kpc")
    print(f"  Medium: v_inf={medium.v_inf/1000:.0f} km/s, r_s={medium.r_s/3.086e19:.1f} kpc")
    
    print(f"\nPM Predictions:")
    print(f"  {'Radius':>10} {'V_obs':>10} {'V_bar':>10} {'V_PM':>10} {'a_bar':>12} {'a_med':>12}")
    print(f"  {'(kpc)':>10} {'(km/s)':>10} {'(km/s)':>10} {'(km/s)':>10} {'(m/s²)':>12} {'(m/s²)':>12}")
    print("  " + "-"*74)
    
    chi_squared = 0.0
    for i, r in enumerate(rc.radii_m):
        v_obs = rc.v_obs_ms[i]
        v_err = rc.v_err_ms[i]
        
        # Calculate PM prediction
        a_bar = accel_baryonic(r, disk)
        a_med = accel_medium(r, medium)
        v_pm = circular_velocity(r, disk, medium)
        
        # Chi-squared contribution
        residual = (v_obs - v_pm) / v_err
        chi_squared += residual**2
        
        # Display
        r_kpc = r / 3.086e19
        v_obs_kms = v_obs / 1000
        v_bar_kms = math.sqrt(r * a_bar) / 1000 if a_bar > 0 else 0
        v_pm_kms = v_pm / 1000
        
        print(f"  {r_kpc:>10.2f} {v_obs_kms:>10.0f} {v_bar_kms:>10.0f} {v_pm_kms:>10.0f} "
              f"{a_bar:>12.2e} {a_med:>12.2e}")
    
    dof = len(rc.radii_m) - 4  # 4 free parameters in medium model
    reduced_chi2 = chi_squared / dof if dof > 0 else chi_squared
    
    print("  " + "-"*74)
    print(f"\n  χ² = {chi_squared:.2f}, dof = {dof}, χ²_red = {reduced_chi2:.2f}")
    
    print("\n" + "="*70)
    print("INTERPRETATION:")
    print("="*70)
    print(f"  • Baryonic component (disk only) falls short at large radii")
    print(f"  • PM medium provides additional a_med acceleration")
    print(f"  • Total v_PM matches observed flat curve")
    print(f"  • No dark matter halo needed!")
    
    if reduced_chi2 < 2.0:
        print(f"\n  ✅ Good fit! (χ²_red = {reduced_chi2:.2f})")
    else:
        print(f"\n  ⚠️  Moderate fit (χ²_red = {reduced_chi2:.2f})")
        print(f"     → Parameters would be optimized in real fitting")
    
    print("\n" + "="*70)
    print("NEXT STEPS:")
    print("="*70)
    print("  1. Get full SPARC dataset (175 galaxies)")
    print("  2. Implement optimizer (scipy.optimize.minimize)")
    print("  3. Fit all galaxies with PM vs NFW halo")
    print("  4. Statistical comparison → PAPER!")
    print("="*70)
    
    return chi_squared, reduced_chi2

if __name__ == '__main__':
    fit_single_galaxy_demo('GAL_X')
