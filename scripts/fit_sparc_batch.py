#!/usr/bin/env python3
"""Batch fit SPARC rotation curves using PM model.

Requires SPARC data in data/SPARC/ directory.
Downloads from https://astroweb.case.edu/SPARC/ if needed.
"""

import sys
import os
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Dict, Tuple
import csv

# Check venv
if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    print("ERROR: Virtual environment not activated!")
    print("Run: source venv/bin/activate")
    sys.exit(1)

# Constants
G = 6.67430e-11  # m^3/(kg·s^2)
c = 299792458.0  # m/s
KPC_TO_M = 3.0856775814913673e19
KM_TO_M = 1000.0

# Load calibration
_CAL_PATH = Path(__file__).resolve().parent.parent / 'tests' / 'calibration.json'
with open(_CAL_PATH) as f:
    CALIB = json.load(f)
    MU_COEFF = CALIB['mu_coeff']


@dataclass
class RotationCurve:
    """Simple rotation curve container."""
    name: str
    radii_m: np.ndarray
    v_obs_ms: np.ndarray
    v_err_ms: np.ndarray
    
    
def load_sparc_galaxy(path: Path) -> RotationCurve:
    """Load single SPARC galaxy file."""
    radii, vobs, verr = [], [], []
    
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # Try various column name variants
                r = float(row.get('R_kpc', row.get('Rad', 0))) * KPC_TO_M
                v = float(row.get('V_obs', row.get('Vobs', row.get('V_obs_kms', 0))))
                # V_obs might be in km/s or m/s - check magnitude
                if v < 1000:  # Likely km/s
                    v *= KM_TO_M
                e = float(row.get('Err_V', row.get('errV', row.get('eV_kms', v*0.1))))
                if e < 1000 and v > 1000:  # Error in km/s, velocity already converted
                    e *= KM_TO_M
                if r > 0 and v > 0:
                    radii.append(r)
                    vobs.append(v)
                    verr.append(e if e > 0 else v * 0.1)
            except (ValueError, KeyError, TypeError):
                continue
                
    return RotationCurve(
        name=path.stem,
        radii_m=np.array(radii),
        v_obs_ms=np.array(vobs),
        v_err_ms=np.array(verr)
    )


def pm_circular_velocity(r: float, M_d: float, R_d: float, 
                         v_inf: float, r_s: float, r_c: float, m: float) -> float:
    """PM model circular velocity at radius r.
    
    v_c^2 = r * (a_bar + a_med)
    where a_bar = GM(<r)/r^2 for exponential disk
    and a_med = -c^2 * d(delta_n)/dr with delta_n = mu_coeff * Phi
    """
    if r <= 0:
        return 0.0
        
    # Baryonic: exponential disk
    x = r / R_d
    M_enc = M_d * (1.0 - (1.0 + x) * np.exp(-x))
    a_bar = G * M_enc / (r * r) if M_enc > 0 else 0.0
    
    # Medium contribution via index gradient
    # Phi(r) for disk (approximation for thin exponential)
    Phi = -G * M_d / r * np.exp(-r/R_d)
    
    # delta_n = mu_coeff * Phi with smoothing
    S = (1.0 + (r_c / r)**m)**(-1.0/m) if r_c > 0 and m > 0 else 1.0
    factor = (v_inf**2 / c**2) * np.log(1.0 + r/r_s) * S if r_s > 0 else 0.0
    
    # Simplified gradient (numerical derivative would be more accurate)
    dr = r * 0.001
    r1, r2 = max(r - dr, 1e10), r + dr
    
    Phi1 = -G * M_d / r1 * np.exp(-r1/R_d)
    Phi2 = -G * M_d / r2 * np.exp(-r2/R_d)
    S1 = (1.0 + (r_c / r1)**m)**(-1.0/m) if r_c > 0 and m > 0 else 1.0
    S2 = (1.0 + (r_c / r2)**m)**(-1.0/m) if r_c > 0 and m > 0 else 1.0
    
    delta_n1 = MU_COEFF * Phi1 * ((v_inf**2 / c**2) * np.log(1.0 + r1/r_s) * S1 if r_s > 0 else 0.0)
    delta_n2 = MU_COEFF * Phi2 * ((v_inf**2 / c**2) * np.log(1.0 + r2/r_s) * S2 if r_s > 0 else 0.0)
    
    grad_n = (delta_n2 - delta_n1) / (2.0 * dr)
    a_med = -c**2 * grad_n
    
    a_total = a_bar + a_med
    return np.sqrt(r * a_total) if a_total > 0 else 0.0


def chi_square(v_obs, v_err, v_model) -> float:
    """Compute chi-square."""
    valid = np.isfinite(v_obs) & np.isfinite(v_model) & (v_err > 0)
    if not np.any(valid):
        return np.inf
    residual = (v_obs[valid] - v_model[valid]) / v_err[valid]
    return np.sum(residual**2)


def fit_galaxy(rc: RotationCurve, n_trials: int = 500) -> Dict:
    """Simple random search fit."""
    # Estimate total mass from outer velocity
    v_max = np.max(rc.v_obs_ms)
    r_max = np.max(rc.radii_m)
    M_est = v_max**2 * r_max / G
    
    # Parameter bounds
    M_d_bounds = (M_est * 0.1, M_est * 10.0)
    R_d_bounds = (r_max * 0.05, r_max * 0.5)
    v_inf_bounds = (v_max * 0.5, v_max * 2.0)
    r_s_bounds = (r_max * 0.1, r_max * 5.0)
    r_c_bounds = (r_max * 0.01, r_max * 0.5)
    m_bounds = (1.0, 4.0)
    
    best_chi2 = np.inf
    best_params = None
    
    for _ in range(n_trials):
        M_d = np.random.uniform(*M_d_bounds)
        R_d = np.random.uniform(*R_d_bounds)
        v_inf = np.random.uniform(*v_inf_bounds)
        r_s = np.random.uniform(*r_s_bounds)
        r_c = np.random.uniform(*r_c_bounds)
        m = np.random.uniform(*m_bounds)
        
        v_model = np.array([pm_circular_velocity(r, M_d, R_d, v_inf, r_s, r_c, m) 
                           for r in rc.radii_m])
        
        chi2 = chi_square(rc.v_obs_ms, rc.v_err_ms, v_model)
        
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_params = {
                'M_d': M_d, 'R_d': R_d, 'v_inf': v_inf,
                'r_s': r_s, 'r_c': r_c, 'm': m
            }
    
    dof = len(rc.radii_m) - 6
    return {
        'params': best_params,
        'chi2': best_chi2,
        'dof': dof,
        'chi2_dof': best_chi2 / dof if dof > 0 else np.inf
    }


def download_sparc_data():
    """Download SPARC data files."""
    print("Checking for SPARC data...")
    data_dir = Path('data/SPARC')
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if we already have data
    csv_files = list(data_dir.glob('*.csv'))
    if len(csv_files) > 0:
        print(f"Found {len(csv_files)} SPARC data files")
        return data_dir
        
    print("\nSPARC data not found locally.")
    print("Please download from: https://astroweb.case.edu/SPARC/")
    print("Extract rotation curve files to: data/SPARC/")
    return None


def main():
    """Batch fit all SPARC galaxies."""
    data_dir = download_sparc_data()
    if data_dir is None:
        return
        
    csv_files = sorted(data_dir.glob('*.csv'))
    if not csv_files:
        print("No CSV files found in data/SPARC/")
        return
        
    print(f"\nFitting {len(csv_files)} galaxies...")
    results = []
    
    for i, path in enumerate(csv_files, 1):
        try:
            print(f"\n[{i}/{len(csv_files)}] {path.stem}...", end=' ', flush=True)
            rc = load_sparc_galaxy(path)
            
            if len(rc.radii_m) < 7:
                print("SKIP (too few points)")
                continue
                
            fit_result = fit_galaxy(rc, n_trials=500)
            fit_result['name'] = rc.name
            results.append(fit_result)
            
            print(f"χ²/dof = {fit_result['chi2_dof']:.2f}")
            
        except Exception as e:
            print(f"ERROR: {e}")
            continue
    
    # Save results
    output = Path('sparc_fit_results.json')
    with open(output, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n\nResults saved to: {output}")
    
    # Summary statistics
    chi2_dof = [r['chi2_dof'] for r in results if np.isfinite(r['chi2_dof'])]
    if chi2_dof:
        print(f"\nSummary: {len(results)} successful fits")
        print(f"  Mean χ²/dof: {np.mean(chi2_dof):.2f}")
        print(f"  Median χ²/dof: {np.median(chi2_dof):.2f}")
        print(f"  Good fits (χ²/dof < 2): {sum(x < 2.0 for x in chi2_dof)}")


if __name__ == '__main__':
    main()
