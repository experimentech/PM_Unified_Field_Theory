#!/usr/bin/env python3
"""
Batch fit PM model with CONSTRAINED r_c/r_s ratio.

Based on empirical observation that excellent fits cluster around r_c/r_s ≈ 0.55,
this script tests whether imposing that constraint:
1. Improves overall fit quality
2. Reduces parameter space degeneracy
3. Reveals a fundamental PM structural ratio

Usage:
  source venv/bin/activate
  python3 batch_fit_sparc_constrained.py [data_file] [--ratio RATIO]
"""

import sys
import json
from pathlib import Path
from dataclasses import asdict
import argparse

# Venv check
if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    print("❌ Virtual environment not active! Run: source venv/bin/activate")
    sys.exit(1)

sys.path.insert(0, str(Path(__file__).parent / 'legacy' / 'Pushing-Medium' / 'src'))

from galaxy_dynamics.data import load_sparc_real
from galaxy_dynamics.fitting import fit_rotation_curve, compute_residual_metrics
from galaxy_dynamics.halos import fit_halo_rotation_curve
import statistics
import numpy as np
from scipy.optimize import differential_evolution

KPC_TO_M = 3.0856775814913673e19
M_SUN = 1.989e30

def fit_with_ratio_constraint(rc, disk_bounds, medium_bounds, ratio=0.55, n_random=300, n_refine=100):
    """
    Fit PM model with r_c = ratio * r_s constraint.
    
    Reduces from 6 free parameters to 5:
    - M_d, R_d (disk)
    - v_inf, r_s, m (medium)
    - r_c is computed as r_c = ratio * r_s
    """
    
    # Combined bounds for constrained fit
    bounds = [
        disk_bounds['M_d'],
        disk_bounds['R_d'],
        medium_bounds['v_inf'],
        medium_bounds['r_s'],
        # r_c is constrained by ratio
        medium_bounds['m'],
    ]
    
    def objective(params):
        M_d, R_d, v_inf, r_s, m = params
        r_c = ratio * r_s  # CONSTRAINT
        
        try:
            from galaxy_dynamics.models import PMGalaxy
            from galaxy_dynamics.disk_models import ExponentialDisk
            from galaxy_dynamics.medium import PMParameters
            
            disk = ExponentialDisk(M_d=M_d, R_d=R_d)
            medium = PMParameters(v_inf=v_inf, r_s=r_s, r_c=r_c, m=m)
            galaxy = PMGalaxy(disk=disk, medium=medium)
            
            v_model = galaxy.rotation_curve(rc.radii_m)
            residuals = (rc.v_obs_ms - v_model) / rc.v_err_ms
            chi2 = np.sum(residuals**2)
            
            return chi2
        except Exception:
            return 1e10
    
    # Random search
    best_chi2 = float('inf')
    best_params = None
    
    for _ in range(n_random):
        params = [np.random.uniform(b[0], b[1]) for b in bounds]
        chi2 = objective(params)
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_params = params
    
    # Refine
    result = differential_evolution(
        objective,
        bounds=bounds,
        maxiter=n_refine,
        x0=best_params,
        seed=42,
        atol=0.01,
        tol=0.01,
    )
    
    M_d, R_d, v_inf, r_s, m = result.x
    r_c = ratio * r_s
    
    from galaxy_dynamics.models import PMGalaxy
    from galaxy_dynamics.disk_models import ExponentialDisk
    from galaxy_dynamics.medium import PMParameters
    
    disk = ExponentialDisk(M_d=M_d, R_d=R_d)
    medium = PMParameters(v_inf=v_inf, r_s=r_s, r_c=r_c, m=m)
    
    return {
        'chi2': result.fun,
        'disk': disk,
        'medium': medium,
        'ratio': ratio,
    }

def batch_fit_constrained(data_file, ratio=0.55, n_random=300, n_refine=100):
    """Fit all galaxies with r_c/r_s constraint."""
    
    print("="*80)
    print(f"PM BATCH FITTING WITH RATIO CONSTRAINT: r_c/r_s = {ratio:.2f}")
    print("="*80)
    print(f"\nData file: {data_file}")
    
    # Load galaxies
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
        'r_c': (1e19, 1e21),  # Not used directly
        'm': (0.5, 4.0),
    }
    
    # Constrained fits
    print(f"\n" + "-"*80)
    print(f"FITTING WITH CONSTRAINT (5 parameters)...")
    print(f"  r_c = {ratio:.2f} × r_s")
    print(f"  n_random={n_random}, n_refine={n_refine}")
    print("-"*80)
    
    results_constrained = {}
    for name, rc in galaxies.items():
        print(f"  Fitting {name}...")
        try:
            fit = fit_with_ratio_constraint(
                rc, disk_bounds, medium_bounds,
                ratio=ratio,
                n_random=n_random,
                n_refine=n_refine,
            )
            results_constrained[name] = fit
            print(f"    χ² = {fit['chi2']:.2f}")
        except Exception as e:
            print(f"    ⚠️  Failed: {e}")
            results_constrained[name] = None
    
    print(f"✓ Constrained fits complete for {len([r for r in results_constrained.values() if r])} galaxies")
    
    # Load unconstrained fits for comparison
    print("\nLoading unconstrained fits for comparison...")
    try:
        with open('data/pm_sparc_fits.json', 'r') as f:
            results_unconstrained = json.load(f)
        print(f"✓ Loaded {len(results_unconstrained)} unconstrained fits")
    except FileNotFoundError:
        print("⚠️  No unconstrained fits found. Run batch_fit_sparc.py first.")
        results_unconstrained = {}
    
    # Analysis
    print("\n" + "="*80)
    print("COMPARISON: CONSTRAINED vs UNCONSTRAINED")
    print("="*80)
    
    chi2_constrained = [r['chi2'] for r in results_constrained.values() if r is not None]
    
    if chi2_constrained:
        print(f"\nConstrained (r_c/r_s = {ratio:.2f}):")
        print(f"  Mean χ²: {statistics.mean(chi2_constrained):.2f}")
        print(f"  Median χ²: {statistics.median(chi2_constrained):.2f}")
        print(f"  Std dev: {statistics.stdev(chi2_constrained):.2f}" if len(chi2_constrained) > 1 else "  Std dev: N/A")
    
    if results_unconstrained:
        chi2_unconstrained = [r['chi2'] for r in results_unconstrained.values()]
        print(f"\nUnconstrained (free r_c and r_s):")
        print(f"  Mean χ²: {statistics.mean(chi2_unconstrained):.2f}")
        print(f"  Median χ²: {statistics.median(chi2_unconstrained):.2f}")
        print(f"  Std dev: {statistics.stdev(chi2_unconstrained):.2f}" if len(chi2_unconstrained) > 1 else "  Std dev: N/A")
        
        # Per-galaxy comparison
        print(f"\n" + "-"*80)
        print("PER-GALAXY IMPACT:")
        print("-"*80)
        print(f"  {'Galaxy':<10} {'χ²_free':>10} {'χ²_const':>10} {'Δχ²':>10} {'r_c/r_s':>10}")
        print(f"  {'':10} {'':>10} {'':>10} {'':>10} {'(free)':>10}")
        print("  " + "-"*68)
        
        better_count = 0
        worse_count = 0
        
        for name in sorted(results_constrained.keys()):
            if name not in results_unconstrained or results_constrained[name] is None:
                continue
            
            chi2_free = results_unconstrained[name]['chi2']
            chi2_const = results_constrained[name]['chi2']
            delta_chi2 = chi2_const - chi2_free
            
            # Original ratio from unconstrained fit
            r_c_orig = results_unconstrained[name]['medium']['r_c']
            r_s_orig = results_unconstrained[name]['medium']['r_s']
            ratio_orig = r_c_orig / r_s_orig
            
            marker = ""
            if delta_chi2 < 0:
                marker = "✓ Better"
                better_count += 1
            elif delta_chi2 < 2:
                marker = "≈ Same"
            else:
                marker = "✗ Worse"
                worse_count += 1
            
            print(f"  {name:<10} {chi2_free:>10.2f} {chi2_const:>10.2f} {delta_chi2:>+10.2f} {ratio_orig:>10.3f}  {marker}")
        
        print("  " + "-"*68)
        print(f"\n  Summary:")
        print(f"    Better with constraint: {better_count}")
        print(f"    Worse with constraint: {worse_count}")
        print(f"    Constraint cost: {statistics.mean([results_constrained[n]['chi2'] - results_unconstrained[n]['chi2'] for n in results_constrained.keys() if n in results_unconstrained and results_constrained[n]]):.2f} Δχ² on average")
    
    # Save constrained results
    results_serializable = {}
    for name, r in results_constrained.items():
        if r is None:
            continue
        results_serializable[name] = {
            'chi2': r['chi2'],
            'disk': {'M_d': r['disk'].M_d, 'R_d': r['disk'].R_d},
            'medium': {
                'v_inf': r['medium'].v_inf,
                'r_s': r['medium'].r_s,
                'r_c': r['medium'].r_c,
                'm': r['medium'].m,
            },
            'ratio_constraint': ratio,
        }
    
    with open(f'data/pm_sparc_fits_constrained_{ratio:.2f}.json', 'w') as f:
        json.dump(results_serializable, f, indent=2)
    
    print(f"\n✓ Saved: data/pm_sparc_fits_constrained_{ratio:.2f}.json")
    
    # Assessment
    print("\n" + "="*80)
    print("VERDICT:")
    print("="*80)
    
    if results_unconstrained and chi2_constrained:
        avg_cost = statistics.mean([results_constrained[n]['chi2'] - results_unconstrained[n]['chi2'] 
                                   for n in results_constrained.keys() 
                                   if n in results_unconstrained and results_constrained[n]])
        
        if avg_cost < 1:
            print(f"  ✅ MINIMAL COST: Constraint barely degrades fits (Δχ² ≈ {avg_cost:.2f})")
            print(f"  → This ratio appears to be a STRUCTURAL PREFERENCE of PM.")
            print(f"  → Consider making r_c/r_s ≈ {ratio:.2f} a theoretical prediction.")
        elif avg_cost < 5:
            print(f"  ✓ SMALL COST: Constraint degrades fits slightly (Δχ² ≈ {avg_cost:.2f})")
            print(f"  → r_c/r_s ≈ {ratio:.2f} is a good approximation for most galaxies.")
        else:
            print(f"  ⚠️  SIGNIFICANT COST: Constraint degrades fits (Δχ² ≈ {avg_cost:.2f})")
            print(f"  → r_c and r_s may need independent freedom.")
    
    print("="*80)
    
    return results_constrained

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fit SPARC with constrained r_c/r_s ratio')
    parser.add_argument('data_file', nargs='?', 
                       default='legacy/Pushing-Medium/tests/data/sparc_sample.csv',
                       help='Path to SPARC data file')
    parser.add_argument('--ratio', type=float, default=0.55,
                       help='r_c/r_s ratio constraint (default: 0.55)')
    parser.add_argument('--n-random', type=int, default=300,
                       help='Number of random search iterations')
    parser.add_argument('--n-refine', type=int, default=100,
                       help='Number of refinement iterations')
    
    args = parser.parse_args()
    
    results = batch_fit_constrained(
        args.data_file, 
        ratio=args.ratio,
        n_random=args.n_random,
        n_refine=args.n_refine
    )
    print("\n✅ Constrained batch fitting complete!")
