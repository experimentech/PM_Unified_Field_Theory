#!/usr/bin/env python3
"""
Fit PM model to full SPARC dataset (175 galaxies).
This is the real test of PM vs dark matter.
"""
import sys
import json
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]

# Venv check
if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    print("❌ Virtual environment not active! Run: source .venv/bin/activate")
    sys.exit(1)

sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'legacy' / 'Pushing-Medium' / 'src'))

from load_sparc_full import load_all_sparc
from galaxy_dynamics.fitting import fit_population
import statistics

KPC_TO_M = 3.0856775814913673e19
M_SUN = 1.989e30

def main():
    print("="*80)
    print("PM vs DARK MATTER: FULL SPARC DATASET (175 GALAXIES)")
    print("="*80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Load all 175 galaxies
    print("\nLoading SPARC rotation curves...")
    galaxies = load_all_sparc()
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
    
    # Fit all galaxies
    print(f"\n" + "-"*80)
    print(f"FITTING PM MODEL...")
    print(f"  n_random=300, n_refine=100")
    print(f"  This will take a while (~10-20 min)...")
    print("-"*80 + "\n")
    
    results_pm = fit_population(
        galaxies,
        disk_bounds=disk_bounds,
        medium_bounds=medium_bounds,
        n_random=300,
        n_refine=100,
    )
    
    print(f"\n✓ Completed {len(results_pm)} fits")
    
    # Analysis
    print("\n" + "="*80)
    print("RESULTS SUMMARY")
    print("="*80)
    
    chi2_pm = [r['chi2'] for r in results_pm.values()]
    
    print(f"\nPM Model Statistics:")
    print(f"  Galaxies fitted: {len(chi2_pm)}")
    print(f"  Mean χ²: {statistics.mean(chi2_pm):.2f}")
    print(f"  Median χ²: {statistics.median(chi2_pm):.2f}")
    print(f"  Std dev: {statistics.stdev(chi2_pm):.2f}")
    print(f"  Min χ²: {min(chi2_pm):.2f}")
    print(f"  Max χ²: {max(chi2_pm):.2f}")
    
    # Quality breakdown
    excellent = sum(1 for c in chi2_pm if c < 2)
    good = sum(1 for c in chi2_pm if 2 <= c < 5)
    acceptable = sum(1 for c in chi2_pm if 5 <= c < 10)
    poor = sum(1 for c in chi2_pm if c >= 10)
    
    print(f"\nFit Quality:")
    print(f"  Excellent (χ² < 2):  {excellent:3d} ({100*excellent/len(chi2_pm):5.1f}%)")
    print(f"  Good (2 ≤ χ² < 5):   {good:3d} ({100*good/len(chi2_pm):5.1f}%)")
    print(f"  Acceptable (5 ≤ χ² < 10): {acceptable:3d} ({100*acceptable/len(chi2_pm):5.1f}%)")
    print(f"  Poor (χ² ≥ 10):      {poor:3d} ({100*poor/len(chi2_pm):5.1f}%)")
    
    # Parameter distributions
    print(f"\n" + "-"*80)
    print("PM PARAMETER DISTRIBUTIONS:")
    print("-"*80)
    
    v_inf_list = [r['medium'].v_inf / 1000 for r in results_pm.values()]
    r_s_list = [r['medium'].r_s / KPC_TO_M for r in results_pm.values()]
    M_d_list = [r['disk'].M_d / M_SUN for r in results_pm.values()]
    
    print(f"  v_inf: {statistics.mean(v_inf_list):6.1f} ± {statistics.stdev(v_inf_list):6.1f} km/s")
    print(f"  r_s:   {statistics.mean(r_s_list):6.1f} ± {statistics.stdev(r_s_list):6.1f} kpc")
    print(f"  M_d:   {statistics.mean(M_d_list):.2e} ± {statistics.stdev(M_d_list):.2e} M_☉")
    
    # Save results
    print("\n" + "-"*80)
    print("SAVING RESULTS...")
    print("-"*80)
    
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
    
    output_file = ROOT / 'results' / 'sparc_fit_results.json'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(results_serializable, f, indent=2)
    
    print(f"✓ Saved: {output_file}")
    
    # Final assessment
    print("\n" + "="*80)
    print("VERDICT:")
    print("="*80)
    
    mean_chi2 = statistics.mean(chi2_pm)
    success_rate = 100 * sum(1 for c in chi2_pm if c < 10) / len(chi2_pm)
    
    if mean_chi2 < 5 and success_rate > 80:
        print("  ✅ STRONG SUCCESS")
        print(f"     Mean χ² = {mean_chi2:.2f}")
        print(f"     Success rate = {success_rate:.0f}%")
        print("\n  PM explains galaxy rotation curves without dark matter!")
    elif mean_chi2 < 10 and success_rate > 60:
        print("  ✓ PROMISING")
        print(f"     Mean χ² = {mean_chi2:.2f}")
        print(f"     Success rate = {success_rate:.0f}%")
        print("\n  PM shows potential but needs refinement.")
    else:
        print("  ⚠️  MIXED RESULTS")
        print(f"     Mean χ² = {mean_chi2:.2f}")
        print(f"     Success rate = {success_rate:.0f}%")
        print("\n  PM struggles with many galaxies. Investigate.")
    
    print("\n" + "="*80)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)

if __name__ == '__main__':
    main()
