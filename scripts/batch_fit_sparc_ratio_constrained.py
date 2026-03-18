#!/usr/bin/env python3
"""
Batch fit PM model with CONSTRAINED r_c/r_s ratio.

Tests the hypothesis that excellent PM fits cluster around r_c/r_s ≈ 0.55.
Instead of fitting 6 free parameters, we constrain r_c = ratio × r_s,
leaving only 5 free parameters.

Usage:
  source venv/bin/activate
  python3 batch_fit_sparc_ratio_constrained.py [data_file] [--ratio RATIO]
"""

import sys
import json
from pathlib import Path
from dataclasses import asdict
import argparse

ROOT = Path(__file__).resolve().parents[1]

# Venv check
if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    print("❌ Virtual environment not active! Run: source .venv/bin/activate")
    sys.exit(1)

sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'legacy' / 'Pushing-Medium' / 'src'))

from galaxy_dynamics.data import load_sparc_real
from galaxy_dynamics.fitting import fit_rotation_curve, compute_residual_metrics
from galaxy_dynamics.halos import fit_halo_rotation_curve
import statistics
import random

KPC_TO_M = 3.0856775814913673e19
M_SUN = 1.989e30

def fit_with_ratio_constraint(rc, disk_bounds, medium_bounds, ratio=0.55, n_random=300, n_refine=100):
    """
    Fit with r_c = ratio × r_s constraint using iterative sampling.
    
    Since fit_rotation_curve doesn't support direct constraints, we:
    1. Sample r_s from bounds
    2. Set r_c = ratio × r_s 
    3. Call fit_rotation_curve with r_c fixed
    4. Repeat n_random times, keeping best result
    """
    from galaxy_dynamics.fitting import DiskParams, MediumParams, _generate_model, chi_square
    from dataclasses import replace
    
    rng = random.Random(42)
    
    def sample_param(bounds, name):
        lo, hi = bounds[name]
        return lo + rng.random() * (hi - lo)
    
    best = None
    
    # Random search with constraint
    for _ in range(n_random):
        r_s = sample_param(medium_bounds, 'r_s')
        r_c = ratio * r_s  # APPLY CONSTRAINT
        
        # Sample other parameters
        disk = DiskParams(
            M_d=sample_param(disk_bounds, 'M_d'),
            R_d=sample_param(disk_bounds, 'R_d'),
        )
        medium = MediumParams(
            v_inf=sample_param(medium_bounds, 'v_inf'),
            r_s=r_s,
            r_c=r_c,
            m=sample_param(medium_bounds, 'm'),
        )
        
        model = _generate_model(rc.radii_m, disk, medium)
        chi2 = chi_square(rc.radii_m, rc.v_obs_ms, rc.v_err_ms, model)
        
        if best is None or chi2 < best['chi2']:
            best = {'disk': disk, 'medium': medium, 'chi2': chi2, 'model': model}
    
    # Local refinement maintaining constraint
    def perturb(value, scale):
        return value * (1.0 + rng.uniform(-scale, scale)) if value > 0 else value
    
    for _ in range(n_refine):
        disk_new = replace(best['disk'])
        medium_new = replace(best['medium'])
        
        disk_new.M_d = perturb(disk_new.M_d, 0.3)
        disk_new.R_d = perturb(disk_new.R_d, 0.3)
        medium_new.v_inf = perturb(medium_new.v_inf, 0.15)
        
        r_s_new = perturb(medium_new.r_s, 0.3)
        medium_new.r_s = r_s_new
        medium_new.r_c = ratio * r_s_new  # MAINTAIN CONSTRAINT
        
        medium_new.m = perturb(medium_new.m, 0.3)
        
        model = _generate_model(rc.radii_m, disk_new, medium_new)
        chi2 = chi_square(rc.radii_m, rc.v_obs_ms, rc.v_err_ms, model)
        
        if chi2 < best['chi2']:
            best = {'disk': disk_new, 'medium': medium_new, 'chi2': chi2, 'model': model}
    
    best['radii_m'] = rc.radii_m
    best['ratio'] = ratio
    return best


def batch_fit_constrained(data_file, ratio=0.55, n_random=300, n_refine=100):
    """Fit all galaxies with r_c/r_s constraint and compare to unconstrained."""
    
    print("="*80)
    print(f"PM BATCH FITTING WITH RATIO CONSTRAINT: r_c/r_s = {ratio}")
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
        'r_c': (1e19, 1e21),  # Only for unconstrained
        'm': (0.5, 4.0),
    }
    
    # Fit with constraint
    print(f"\n" + "-"*80)
    print(f"FITTING WITH CONSTRAINT (5 parameters)...")
    print(f"  r_c = {ratio} × r_s")
    print(f"  n_random={n_random}, n_refine={n_refine}")
    print("-"*80)
    
    results_constrained = {}
    for name, rc in galaxies.items():
        print(f"  Fitting {name}...")
        try:
            result = fit_with_ratio_constraint(rc, disk_bounds, medium_bounds, ratio, n_random, n_refine)
            results_constrained[name] = result
            print(f"    ✓ χ² = {result['chi2']:.2f}, r_c/r_s = {ratio}")
        except Exception as e:
            print(f"    ⚠️  Failed: {e}")
    
    print(f"✓ Constrained fits complete for {len(results_constrained)} galaxies")
    
    # Load unconstrained results
    print("\nLoading unconstrained fits for comparison...")
    unconstrained_file = ROOT / 'results' / 'sparc_fit_results.json'
    if unconstrained_file.exists():
        with open(unconstrained_file) as f:
            data = json.load(f)
            results_unconstrained = data.get('pm_fits', {})
        print(f"✓ Loaded {len(results_unconstrained)} unconstrained fits")
    else:
        print("⚠️  No unconstrained results found. Run batch_fit_sparc.py first.")
        results_unconstrained = {}
    
    # Analysis
    print("\n" + "="*80)
    print("COMPARISON: CONSTRAINED vs UNCONSTRAINED")
    print("="*80)
    
    if results_unconstrained:
        chi2_free = [r['chi2'] for r in results_unconstrained.values() if r]
        print(f"\nUnconstrained (free r_c and r_s):")
        print(f"  Mean χ²: {statistics.mean(chi2_free):.2f}")
        print(f"  Median χ²: {statistics.median(chi2_free):.2f}")
        print(f"  Std dev: {statistics.stdev(chi2_free):.2f}")
    
    if results_constrained:
        chi2_const = [r['chi2'] for r in results_constrained.values()]
        print(f"\nConstrained (r_c = {ratio} × r_s):")
        print(f"  Mean χ²: {statistics.mean(chi2_const):.2f}")
        print(f"  Median χ²: {statistics.median(chi2_const):.2f}")
        print(f"  Std dev: {statistics.stdev(chi2_const):.2f}")
    
    # Per-galaxy comparison
    print(f"\n" + "-"*80)
    print("PER-GALAXY IMPACT:")
    print("-"*80)
    print(f"  {'Galaxy':<12}  {'χ²_free':>8}  {'χ²_const':>10}  {'Δχ²':>10}  {'r_c/r_s':>8}")
    print(f"  {'':>12}  {'':>8}  {'':>10}  {'':>10}  {'(free)':>8}")
    print("  " + "-"*68)
    
    better = 0
    worse = 0
    
    for name in sorted(results_constrained.keys()):
        if name in results_unconstrained and results_unconstrained[name]:
            chi2_free = results_unconstrained[name]['chi2']
            chi2_const = results_constrained[name]['chi2']
            delta = chi2_const - chi2_free
            
            # Compute original ratio
            med = results_unconstrained[name]['medium']
            orig_ratio = med['r_c'] / med['r_s']
            
            marker = "→" if abs(delta) < 1.0 else ("↑" if delta > 0 else "↓")
            print(f"  {name:<12}  {chi2_free:8.2f}  {chi2_const:10.2f}  {delta:+10.2f} {marker}  {orig_ratio:8.3f}")
            
            if delta < 0:
                better += 1
            elif delta > 1.0:
                worse += 1
    
    print("  " + "-"*68)
    print(f"\n  Summary:")
    print(f"    Better with constraint: {better}")
    print(f"    Worse with constraint: {worse}")
    
    if results_constrained and results_unconstrained:
        # Overall cost
        common = set(results_constrained.keys()) & set(results_unconstrained.keys())
        common = [n for n in common if results_unconstrained[n]]
        
        if common:
            deltas = [results_constrained[n]['chi2'] - results_unconstrained[n]['chi2'] for n in common]
            avg_cost = statistics.mean(deltas)
            print(f"    Constraint cost: {avg_cost:+.2f} Δχ² on average")
    
    # Save results
    output_file = ROOT / 'results' / 'sparc_constrained_results.json'
    output_data = {
        'ratio': ratio,
        'n_galaxies': len(results_constrained),
        'constrained_fits': {
            name: {
                'disk': asdict(r['disk']),
                'medium': asdict(r['medium']),
                'chi2': r['chi2'],
                'ratio': ratio,
            }
            for name, r in results_constrained.items()
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n✓ Results saved to {output_file}")
    
    return results_constrained


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Batch fit with r_c/r_s constraint')
    default_data = ROOT / 'legacy' / 'Pushing-Medium' / 'tests' / 'data' / 'sparc_sample.csv'
    parser.add_argument('data_file', nargs='?', default=str(default_data))
    parser.add_argument('--ratio', type=float, default=0.55, help='r_c/r_s ratio constraint')
    parser.add_argument('--n-random', type=int, default=300, help='Random samples')
    parser.add_argument('--n-refine', type=int, default=100, help='Refinement iterations')
    
    args = parser.parse_args()
    
    results = batch_fit_constrained(
        args.data_file,
        ratio=args.ratio,
        n_random=args.n_random,
        n_refine=args.n_refine
    )
