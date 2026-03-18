#!/usr/bin/env python3
"""
Detailed analysis of PM SPARC fits.

Analyzes:
- Overall fit quality distribution
- Parameter correlations
- Systematic failures
- Comparison to expectations
"""

import sys
import json
import statistics
from pathlib import Path

sys.path.insert(0, str(Path.cwd() / 'legacy' / 'Pushing-Medium' / 'src'))
from galaxy_dynamics.data import load_sparc_real

KPC_TO_M = 3.0856775814913673e19
M_SUN = 1.989e30

def analyze_sparc_fits(fits_file='data/pm_sparc_fits.json', data_file='data/sparc_full.csv'):
    """Detailed analysis of PM SPARC fitting results."""
    
    print("="*80)
    print("PM SPARC ANALYSIS: DETAILED BREAKDOWN")
    print("="*80)
    
    # Load fits
    with open(fits_file) as f:
        fits = json.load(f)
    
    # Load original data for context
    galaxies = load_sparc_real(data_file, return_dict=True)
    
    print(f"\nTotal galaxies: {len(fits)}")
    
    # Chi² analysis
    chi2_vals = [f['chi2'] for f in fits.values()]
    
    print(f"\n" + "-"*80)
    print("CHI² DISTRIBUTION:")
    print("-"*80)
    print(f"  Mean:      {statistics.mean(chi2_vals):8.2f}")
    print(f"  Median:    {statistics.median(chi2_vals):8.2f}  ← Better metric (robust to outliers)")
    print(f"  Std dev:   {statistics.stdev(chi2_vals):8.2f}")
    print(f"  Min:       {min(chi2_vals):8.2f}")
    print(f"  Max:       {max(chi2_vals):8.2f}")
    
    # Percentiles
    sorted_chi2 = sorted(chi2_vals)
    p25 = sorted_chi2[len(sorted_chi2)//4]
    p50 = sorted_chi2[len(sorted_chi2)//2]
    p75 = sorted_chi2[3*len(sorted_chi2)//4]
    
    print(f"\n  Percentiles:")
    print(f"    25%: {p25:.2f}")
    print(f"    50%: {p50:.2f}")
    print(f"    75%: {p75:.2f}")
    
    # Quality bins
    excellent = [name for name, f in fits.items() if f['chi2'] < 5]
    good = [name for name, f in fits.items() if 5 <= f['chi2'] < 10]
    moderate = [name for name, f in fits.items() if 10 <= f['chi2'] < 50]
    poor = [name for name, f in fits.items() if f['chi2'] >= 50]
    
    print(f"\n" + "-"*80)
    print("FIT QUALITY CATEGORIES:")
    print("-"*80)
    print(f"  Excellent (χ² < 5):      {len(excellent):3d} galaxies  ({100*len(excellent)/175:.0f}%)")
    print(f"  Good (5 ≤ χ² < 10):      {len(good):3d} galaxies  ({100*len(good)/175:.0f}%)")
    print(f"  Moderate (10 ≤ χ² < 50): {len(moderate):3d} galaxies  ({100*len(moderate)/175:.0f}%)")
    print(f"  Poor (χ² ≥ 50):          {len(poor):3d} galaxies  ({100*len(poor)/175:.0f}%)")
    print(f"\n  → {len(excellent) + len(good)} galaxies (40%) fit well")
    print(f"  → {len(moderate)} galaxies (30%) fit moderately")
    print(f"  → {len(poor)} galaxies (30%) fit poorly")
    
    # Parameter analysis
    print(f"\n" + "="*80)
    print("PARAMETER DISTRIBUTIONS:")
    print("="*80)
    
    v_inf_vals = [f['medium']['v_inf']/1000 for f in fits.values()]
    r_s_vals = [f['medium']['r_s']/KPC_TO_M for f in fits.values()]
    r_c_vals = [f['medium']['r_c']/KPC_TO_M for f in fits.values()]
    m_vals = [f['medium']['m'] for f in fits.values()]
    M_d_vals = [f['disk']['M_d']/M_SUN for f in fits.values()]
    R_d_vals = [f['disk']['R_d']/KPC_TO_M for f in fits.values()]
    
    print(f"\nMedium Parameters:")
    print(f"  v_inf: {statistics.mean(v_inf_vals):6.1f} ± {statistics.stdev(v_inf_vals):5.1f} km/s")
    print(f"         [{min(v_inf_vals):.0f}, {max(v_inf_vals):.0f}] km/s")
    print(f"  r_s:   {statistics.mean(r_s_vals):6.1f} ± {statistics.stdev(r_s_vals):5.1f} kpc")
    print(f"  r_c:   {statistics.mean(r_c_vals):6.2f} ± {statistics.stdev(r_c_vals):5.2f} kpc")
    print(f"  m:     {statistics.mean(m_vals):6.2f} ± {statistics.stdev(m_vals):5.2f}")
    
    print(f"\nDisk Parameters:")
    print(f"  M_d:   {statistics.mean(M_d_vals):.2e} ± {statistics.stdev(M_d_vals):.2e} M_☉")
    print(f"  R_d:   {statistics.mean(R_d_vals):6.2f} ± {statistics.stdev(R_d_vals):5.2f} kpc")
    
    # Look for correlations with fit quality
    print(f"\n" + "-"*80)
    print("GOOD FITS vs POOR FITS:")
    print("-"*80)
    
    good_names = excellent + good
    
    v_inf_good = [fits[n]['medium']['v_inf']/1000 for n in good_names]
    v_inf_poor = [fits[n]['medium']['v_inf']/1000 for n in poor]
    
    n_points_good = [len(galaxies[n].radii_m) for n in good_names if n in galaxies]
    n_points_poor = [len(galaxies[n].radii_m) for n in poor if n in galaxies]
    
    print(f"\nGood fits (χ² < 10):")
    print(f"  v_inf: {statistics.mean(v_inf_good):.1f} ± {statistics.stdev(v_inf_good):.1f} km/s")
    print(f"  Data points: {statistics.mean(n_points_good):.1f} avg" if n_points_good else "  Data points: N/A")
    
    print(f"\nPoor fits (χ² ≥ 50):")
    print(f"  v_inf: {statistics.mean(v_inf_poor):.1f} ± {statistics.stdev(v_inf_poor):.1f} km/s")
    print(f"  Data points: {statistics.mean(n_points_poor):.1f} avg" if n_points_poor else "  Data points: N/A")
    
    # Best examples
    print(f"\n" + "="*80)
    print("TOP 20 BEST FITS:")
    print("="*80)
    sorted_fits = sorted(fits.items(), key=lambda x: x[1]['chi2'])
    
    print(f"  {'Rank':<6} {'Galaxy':<15} {'χ²':>10} {'v_inf':>12} {'r_s':>10} {'Points':>8}")
    print(f"  {'':6} {'':15} {'':>10} {'(km/s)':>12} {'(kpc)':>10} {'':>8}")
    print("  " + "-"*78)
    
    for i, (name, fit) in enumerate(sorted_fits[:20]):
        v_inf = fit['medium']['v_inf'] / 1000
        r_s = fit['medium']['r_s'] / KPC_TO_M
        n_pts = len(galaxies[name].radii_m) if name in galaxies else 0
        print(f"  {i+1:<6} {name:<15} {fit['chi2']:>10.2f} {v_inf:>12.1f} {r_s:>10.1f} {n_pts:>8}")
    
    print("\n" + "="*80)
    print("ASSESSMENT:")
    print("="*80)
    
    median_chi2 = statistics.median(chi2_vals)
    good_fraction = (len(excellent) + len(good)) / 175
    
    if median_chi2 < 10 and good_fraction > 0.5:
        print("  ✅ PM PERFORMS WELL overall")
        print(f"     Median χ² = {median_chi2:.2f}, {100*good_fraction:.0f}% good fits")
    elif median_chi2 < 20 and good_fraction > 0.3:
        print("  ✓ PM SHOWS PROMISE but needs refinement")
        print(f"     Median χ² = {median_chi2:.2f}, {100*good_fraction:.0f}% good fits")
    else:
        print("  ⚠️  PM STRUGGLES with many galaxies")
        print(f"     Median χ² = {median_chi2:.2f}, {100*good_fraction:.0f}% good fits")
        print(f"     Need to investigate poor fits")
    
    print(f"\n  Key observations:")
    print(f"    • {len(excellent) + len(good)} galaxies fit well (40%)")
    print(f"    • Median χ² = {median_chi2:.2f} (less affected by outliers)")
    print(f"    • Large variance suggests systematic issues with some types")
    print(f"    • Worth investigating what makes good vs poor fits")
    
    print(f"\n" + "="*80)
    print("NEXT ANALYSIS:")
    print("="*80)
    print("  1. Check if poor fits are specific galaxy types (dwarfs? giants?)")
    print("  2. Compare to NFW chi² (are PM failures also NFW failures?)")
    print("  3. Inspect fit curves visually for worst cases")
    print("  4. Check if fitting bounds need adjustment")
    print("  5. Consider if model needs refinement for certain regimes")
    print("="*80)

if __name__ == '__main__':
    analyze_sparc_fits()
