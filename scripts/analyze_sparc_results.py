#!/usr/bin/env python3
"""
Analyze SPARC fitting results to understand what went wrong.
"""
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Venv check
if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    print("❌ Virtual environment not active! Run: source .venv/bin/activate")
    sys.exit(1)

sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'legacy' / 'Pushing-Medium' / 'src'))

from load_sparc_full import load_all_sparc
import statistics

KPC_TO_M = 3.0856775814913673e19
M_SUN = 1.989e30

# Load results
results_path = ROOT / 'results' / 'sparc_fit_results.json'
with open(results_path, 'r') as f:
    results = json.load(f)

galaxies = load_all_sparc()

print("="*80)
print("SPARC FITTING ANALYSIS")
print("="*80)

# Sort by chi2
sorted_results = sorted(results.items(), key=lambda x: x[1]['chi2'])

print("\n--- BEST FITS (χ² < 2) ---")
for name, r in sorted_results[:29]:
    if r['chi2'] < 2:
        gal = galaxies[name]
        print(f"{name:20s} χ²={r['chi2']:6.2f}  n_pts={len(gal.radii_m):3d}  "
              f"v_inf={r['medium']['v_inf']/1000:5.1f} km/s")

print("\n--- WORST FITS (χ² > 100) ---")
worst = [item for item in sorted_results if item[1]['chi2'] > 100]
print(f"Total: {len(worst)} galaxies")
for name, r in worst[-10:]:
    gal = galaxies[name]
    print(f"{name:20s} χ²={r['chi2']:8.2f}  n_pts={len(gal.radii_m):3d}  "
          f"v_inf={r['medium']['v_inf']/1000:5.1f} km/s")

# Check for patterns
print("\n" + "="*80)
print("DIAGNOSTIC: Data point count vs fit quality")
print("="*80)

chi2_by_npts = {}
for name, r in results.items():
    gal = galaxies[name]
    n_pts = len(gal.radii_m)
    if n_pts not in chi2_by_npts:
        chi2_by_npts[n_pts] = []
    chi2_by_npts[n_pts].append(r['chi2'])

for n_pts in sorted(chi2_by_npts.keys()):
    chi2_list = chi2_by_npts[n_pts]
    mean_chi2 = statistics.mean(chi2_list)
    print(f"  n_pts={n_pts:3d}: {len(chi2_list):3d} galaxies, mean χ²={mean_chi2:8.2f}")

print("\n--- Hypothesis: Are low data point galaxies problematic? ---")
low_npts = [r['chi2'] for name, r in results.items() if len(galaxies[name].radii_m) < 10]
high_npts = [r['chi2'] for name, r in results.items() if len(galaxies[name].radii_m) >= 10]

if low_npts:
    print(f"  Galaxies with <10 points: {len(low_npts)}, mean χ²={statistics.mean(low_npts):.2f}")
if high_npts:
    print(f"  Galaxies with ≥10 points: {len(high_npts)}, mean χ²={statistics.mean(high_npts):.2f}")
