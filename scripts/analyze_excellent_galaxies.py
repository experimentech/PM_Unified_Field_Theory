#!/usr/bin/env python3
"""
Analyze the characteristics of galaxies with excellent PM fits.
"""

import json
import numpy as np

def analyze_fits():
    with open('sparc_fit_results.json', 'r') as f:
        results = json.load(f)
    
    excellent = []
    good = []
    poor = []
    
    for name, data in results.items():
        if 'error' in data or 'chi2' not in data:
            continue
        
        chi2 = data['chi2']
        rc = data['medium']['r_c']
        rs = data['medium']['r_s']
        ratio = rc / rs
        
        entry = {
            'name': name,
            'chi2': chi2,
            'rc': rc / 3.086e21,  # Convert to kpc
            'rs': rs / 3.086e21,
            'ratio': ratio,
            'rms': data['metrics']['rms'],
            'frac_rms': data['metrics']['frac_rms']
        }
        
        if chi2 < 2:
            excellent.append(entry)
        elif chi2 < 5:
            good.append(entry)
        else:
            poor.append(entry)
    
    print(f"=== FIT QUALITY BREAKDOWN ===\n")
    print(f"Excellent (χ² < 2): {len(excellent)}")
    print(f"Good (2 ≤ χ² < 5): {len(good)}")
    print(f"Poor (χ² ≥ 5): {len(poor)}")
    print(f"Total: {len(excellent) + len(good) + len(poor)}\n")
    
    if excellent:
        print("=== EXCELLENT FITS ===\n")
        ratios = [g['ratio'] for g in excellent]
        print(f"r_c/r_s ratio: {np.median(ratios):.3f} ± {np.std(ratios):.3f}")
        print(f"Range: [{np.min(ratios):.3f}, {np.max(ratios):.3f}]\n")
        
        print(f"{'Galaxy':<15} {'χ²':<8} {'r_c/r_s':<8} {'r_c[kpc]':<10} {'r_s[kpc]':<10}")
        print("-" * 65)
        for g in sorted(excellent, key=lambda x: x['chi2'])[:15]:
            print(f"{g['name']:<15} {g['chi2']:<8.2f} {g['ratio']:<8.3f} {g['rc']:<10.2f} {g['rs']:<10.2f}")
    
    if poor:
        print("\n\n=== POOR FITS ===\n")
        ratios = [g['ratio'] for g in poor]
        print(f"r_c/r_s ratio: {np.median(ratios):.3f} ± {np.std(ratios):.3f}")
        print(f"Range: [{np.min(ratios):.3f}, {np.max(ratios):.3f}]\n")
        
        print(f"{'Galaxy':<15} {'χ²':<8} {'r_c/r_s':<8} {'r_c[kpc]':<10} {'r_s[kpc]':<10}")
        print("-" * 65)
        for g in sorted(poor, key=lambda x: x['chi2'], reverse=True)[:15]:
            print(f"{g['name']:<15} {g['chi2']:<8.1f} {g['ratio']:<8.3f} {g['rc']:<10.2f} {g['rs']:<10.2f}")
    
    if excellent and poor:
        print("\n\n=== STATISTICAL COMPARISON ===\n")
        exc_ratios = [g['ratio'] for g in excellent]
        poor_ratios = [g['ratio'] for g in poor]
        
        separation = abs(np.median(exc_ratios) - np.median(poor_ratios))
        combined_std = np.sqrt(np.std(exc_ratios)**2 + np.std(poor_ratios)**2)
        
        print(f"Median r_c/r_s:")
        print(f"  Excellent: {np.median(exc_ratios):.3f}")
        print(f"  Poor:      {np.median(poor_ratios):.3f}")
        print(f"  Separation: {separation:.3f} ({separation/combined_std:.2f}σ)")

if __name__ == '__main__':
    analyze_fits()
