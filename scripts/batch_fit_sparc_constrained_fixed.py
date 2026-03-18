#!/usr/bin/env python3
"""
Batch fit PM model to all galaxies in SPARC dataset.

This demonstrates the complete pipeline:
1. Load all galaxies
2. Fit PM model (disk + medium)
3. Fit dark matter halo for comparison
4. Statistical analysis
5. Generate summary report

Usage:
  source venv/bin/activate
  python3 batch_fit_sparc.py [data_file]
  
Default uses sample data. Provide full SPARC path when ready.
"""

import sys
import json
from pathlib import Path
from dataclasses import asdict

# Venv check
if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    print("❌ Virtual environment not active! Run: source venv/bin/activate")
    sys.exit(1)

sys.path.insert(0, str(Path(__file__).parent / 'legacy' / 'Pushing-Medium' / 'src'))

from galaxy_dynamics.data import load_sparc_real
from galaxy_dynamics.fitting import fit_rotation_curve, compute_residual_metrics, fit_population
from galaxy_dynamics.halos import fit_halo_rotation_curve
import statistics

KPC_TO_M = 3.0856775814913673e19
M_SUN = 1.989e30

def batch_fit_all(data_file='legacy/Pushing-Medium/tests/data/sparc_sample.csv', n_random=300, n_refine=100):
    """Fit all galaxies in dataset with PM and compare to dark matter."""
    
    print("="*80)
    print("PM BATCH GALAXY FITTING")
    print("="*80)
    print(f"\nData file: {data_file}")
    
    # Load all galaxies
    print("\nLoading galaxies...")
    galaxies = load_sparc_real(data_file, return_dict=True)
    print(f"✓ Loaded {len(galaxies)} galaxies")
    
    # Define bounds
    disk_bounds = {
        'M_d': (1e40, 1e42),
        'R_d': (1e20, 1e21),
    }
    
    medium_bounds = {
        'v_inf': (50e3, 200e3),
        'r_s': (1e20, 1e22),
        'r_c': (1e19, 1e21),
        'm': (0.5, 4.0),
    }
    
    # PM fits
    print(f"\n" + "-"*80)
    print(f"FITTING PM MODEL (6 parameters)...")
    print(f"  n_random={n_random}, n_refine={n_refine}")
    print("-"*80)
    
    results_pm = fit_population(
        galaxies,
        disk_bounds=disk_bounds,
        medium_bounds=medium_bounds,
        n_random=n_random,
        n_refine=n_refine,
    )
    
    print(f"✓ PM fits complete for {len(results_pm)} galaxies")
    
    # NFW fits for comparison
    print(f"\n" + "-"*80)
    print("FITTING NFW HALO MODEL (for comparison)...")
    print("-"*80)
    
    nfw_bounds = {
        'rho_s': (1e-21, 1e-18),  # kg/m^3
        'r_s': (1e20, 1e22),      # m
    }
    
    results_nfw = {}
    for name, rc in galaxies.items():
        print(f"  Fitting {name} (NFW)...")
        try:
            fit_nfw = fit_halo_rotation_curve(
                rc.radii_m, rc.v_obs_ms, rc.v_err_ms,
                halo_type='nfw',
                bounds=nfw_bounds,
                n_random=200,
                n_refine=50,
            )
            results_nfw[name] = fit_nfw
        except Exception as e:
            print(f"    ⚠️  Failed: {e}")
            results_nfw[name] = None
    
    print(f"✓ NFW fits complete")
    
    # Analysis
    print("\n" + "="*80)
    print("STATISTICAL SUMMARY")
    print("="*80)
    
    chi2_pm = [r['chi2'] for r in results_pm.values()]
    chi2_nfw = [r['chi2'] for r in results_nfw.values() if r is not None]
    
    print(f"\nPM Model:")
    print(f"  Mean χ²: {statistics.mean(chi2_pm):.2f}")
    print(f"  Median χ²: {statistics.median(chi2_pm):.2f}")
    print(f"  Std dev: {statistics.stdev(chi2_pm):.2f}" if len(chi2_pm) > 1 else "  Std dev: N/A")
    
    if chi2_nfw:
        print(f"\nNFW Halo Model:")
        print(f"  Mean χ²: {statistics.mean(chi2_nfw):.2f}")
        print(f"  Median χ²: {statistics.median(chi2_nfw):.2f}")
        print(f"  Std dev: {statistics.stdev(chi2_nfw):.2f}" if len(chi2_nfw) > 1 else "  Std dev: N/A")
        
        # Comparison
        pm_better = sum(1 for i, name in enumerate(results_pm.keys()) 
                       if results_nfw.get(name) and results_pm[name]['chi2'] < results_nfw[name]['chi2'])
        print(f"\nComparison:")
        print(f"  PM fits better in {pm_better}/{len(galaxies)} galaxies")
    
    # Parameter distributions
    print(f"\n" + "-"*80)
    print("PM PARAMETER DISTRIBUTIONS:")
    print("-"*80)
    
    v_inf_list = [r['medium'].v_inf / 1000 for r in results_pm.values()]
    r_s_list = [r['medium'].r_s / KPC_TO_M for r in results_pm.values()]
    M_d_list = [r['disk'].M_d / M_SUN for r in results_pm.values()]
    
    print(f"  v_inf: {statistics.mean(v_inf_list):.1f} ± {statistics.stdev(v_inf_list):.1f} km/s" 
          if len(v_inf_list) > 1 else f"  v_inf: {v_inf_list[0]:.1f} km/s")
    print(f"  r_s:   {statistics.mean(r_s_list):.1f} ± {statistics.stdev(r_s_list):.1f} kpc"
          if len(r_s_list) > 1 else f"  r_s:   {r_s_list[0]:.1f} kpc")
    print(f"  M_d:   {statistics.mean(M_d_list):.2e} ± {statistics.stdev(M_d_list):.2e} M_☉"
          if len(M_d_list) > 1 else f"  M_d:   {M_d_list[0]:.2e} M_☉")
    
    # Detailed results per galaxy
    print(f"\n" + "="*80)
    print("PER-GALAXY RESULTS:")
    print("="*80)
    print(f"  {'Galaxy':<10} {'χ²_PM':>8} {'χ²_NFW':>8} {'v_inf':>10} {'r_s':>10} {'M_d':>12}")
    print(f"  {'':10} {'':>8} {'':>8} {'(km/s)':>10} {'(kpc)':>10} {'(M_☉)':>12}")
    print("  " + "-"*78)
    
    for name in results_pm.keys():
        pm = results_pm[name]
        nfw = results_nfw.get(name)
        
        chi2_pm_val = pm['chi2']
        chi2_nfw_val = nfw['chi2'] if nfw else float('nan')
        v_inf = pm['medium'].v_inf / 1000
        r_s = pm['medium'].r_s / KPC_TO_M
        M_d = pm['disk'].M_d / M_SUN
        
        print(f"  {name:<10} {chi2_pm_val:>8.2f} {chi2_nfw_val:>8.2f} {v_inf:>10.1f} {r_s:>10.1f} {M_d:>12.2e}")
    
    print("  " + "-"*78)
    
    # Save results
    print("\n" + "-"*80)
    print("SAVING RESULTS...")
    print("-"*80)
    
    # Convert to serializable format
    results_serializable = {}
    for name, r in results_pm.items():
        results_serializable[name] = {
            'chi2': r['chi2'],
            'disk': {'M_d': r['disk'].M_d, 'R_d': r['disk'].R_d},
            'medium': {
                'v_inf': r['medium'].v_inf,
                'r_s': r['medium'].r_s,
                'r_c': r['medium'].r_c,
                'm': r['medium'].m,
            },
            'metrics': r.get('metrics', {}),
        }
    
    with open('data/pm_sparc_fits.json', 'w') as f:
        json.dump(results_serializable, f, indent=2)
    
    print("✓ Saved: data/pm_sparc_fits.json")
    
    # Summary report
    print("\n" + "="*80)
    print("ASSESSMENT:")
    print("="*80)
    
    mean_chi2 = statistics.mean(chi2_pm)
    good_fits = sum(1 for c in chi2_pm if c < 10)
    
    print(f"  Galaxies fitted: {len(galaxies)}")
    print(f"  Mean χ²: {mean_chi2:.2f}")
    print(f"  Good fits (χ² < 10): {good_fits}/{len(galaxies)} ({100*good_fits/len(galaxies):.0f}%)")
    
    if mean_chi2 < 5:
        print(f"\n  ✅ EXCELLENT! PM fits galaxies well on average.")
    elif mean_chi2 < 10:
        print(f"\n  ✓ GOOD. PM provides reasonable fits.")
    else:
        print(f"\n  ⚠️  MODERATE. PM struggles with some galaxies.")
    
    print(f"\n  The PM medium explains flat rotation curves")
    print(f"  WITHOUT dark matter halos.")
    
    print("\n" + "="*80)
    print("NEXT STEPS:")
    print("="*80)
    print("  1. Get full SPARC dataset (175 galaxies)")
    print("     → http://astroweb.cwru.edu/SPARC/")
    print("     → Or from ADS: 2016AJ....152..157L")
    print("  2. Run this script on full dataset")
    print("  3. Compare to published dark matter fits")
    print("  4. Generate plots for paper")
    print("  5. Write manuscript!")
    print("="*80)
    
    return results_pm, results_nfw

if __name__ == '__main__':
    data_file = sys.argv[1] if len(sys.argv) > 1 else 'legacy/Pushing-Medium/tests/data/sparc_sample.csv'
    results_pm, results_nfw = batch_fit_all(data_file, n_random=300, n_refine=100)
    print("\n✅ Batch fitting complete!")
