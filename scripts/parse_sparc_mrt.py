#!/usr/bin/env python3
"""Parse SPARC MassModels catalog and run batch fitting."""

import sys
import json
from pathlib import Path
from collections import defaultdict
import numpy as np

# Activate venv check
if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    print("ERROR: Virtual environment not activated. Run: source venv/bin/activate")
    sys.exit(1)

from src.pushing_medium.core import PushingMedium
from src.galaxy_dynamics.models import NFWHalo, ExponentialDisk
from scipy.optimize import minimize


def parse_sparc_mrt(filepath):
    """Parse SPARC MRT file into galaxy rotation curves."""
    galaxies = defaultdict(lambda: {
        'R_kpc': [], 'Vobs': [], 'errV': [],
        'Vdisk': [], 'Vbul': [], 'Vgas': []
    })
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('-'):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            
            try:
                name = parts[0]
                dist_mpc = float(parts[1])
                r_kpc = float(parts[2])
                vobs = float(parts[3])
                errv = float(parts[4])
                vdisk = float(parts[5])
                vbul = float(parts[6])
                vgas = float(parts[7])
                
                galaxies[name]['R_kpc'].append(r_kpc)
                galaxies[name]['Vobs'].append(vobs)
                galaxies[name]['errV'].append(errv)
                galaxies[name]['Vdisk'].append(vdisk)
                galaxies[name]['Vbul'].append(vbul)
                galaxies[name]['Vgas'].append(vgas)
                galaxies[name]['distance_mpc'] = dist_mpc
            except (ValueError, IndexError):
                continue
    
    return dict(galaxies)


def fit_galaxy_pm(galaxy_data, pm):
    """Fit PM model to galaxy rotation curve."""
    R = np.array(galaxy_data['R_kpc']) * 3.086e19  # kpc to m
    V_obs = np.array(galaxy_data['Vobs']) * 1000  # km/s to m/s
    V_err = np.array(galaxy_data['errV']) * 1000
    
    # Baryonic velocity from components
    V_disk = np.array(galaxy_data['Vdisk']) * 1000
    V_bul = np.array(galaxy_data['Vbul']) * 1000
    V_gas = np.array(galaxy_data['Vgas']) * 1000
    V_bar = np.sqrt(V_disk**2 + V_bul**2 + V_gas**2)
    
    def chi2(params):
        rho0, r_s = params
        if rho0 <= 0 or r_s <= 0:
            return 1e10
        
        V_pred = []
        for r_val, v_bar in zip(R, V_bar):
            phi = pm.gravitational_potential(r_val, rho0, r_s)
            v_pm = np.sqrt(max(0, v_bar**2 - 2*phi))
            V_pred.append(v_pm)
        
        V_pred = np.array(V_pred)
        weights = 1.0 / (V_err**2 + 1.0)
        return np.sum(weights * (V_obs - V_pred)**2) / len(R)
    
    result = minimize(chi2, [1e7, 1e4], method='Nelder-Mead',
                     options={'maxiter': 1000})
    
    return {
        'rho0': result.x[0],
        'r_s': result.x[1],
        'chi2': result.fun,
        'success': result.success
    }


def fit_galaxy_nfw(galaxy_data):
    """Fit NFW halo model."""
    R = np.array(galaxy_data['R_kpc']) * 3.086e19
    V_obs = np.array(galaxy_data['Vobs']) * 1000
    V_err = np.array(galaxy_data['errV']) * 1000
    V_disk = np.array(galaxy_data['Vdisk']) * 1000
    V_bul = np.array(galaxy_data['Vbul']) * 1000
    V_gas = np.array(galaxy_data['Vgas']) * 1000
    V_bar = np.sqrt(V_disk**2 + V_bul**2 + V_gas**2)
    
    nfw = NFWHalo(rho_s=1e7, r_s=1e4)
    
    def chi2(params):
        rho_s, r_s = params
        if rho_s <= 0 or r_s <= 0:
            return 1e10
        nfw.rho_s = rho_s
        nfw.r_s = r_s
        V_pred = np.array([np.sqrt(v_bar**2 + nfw.rotation_velocity(r)**2)
                           for r, v_bar in zip(R, V_bar)])
        weights = 1.0 / (V_err**2 + 1.0)
        return np.sum(weights * (V_obs - V_pred)**2) / len(R)
    
    result = minimize(chi2, [1e7, 1e4], method='Nelder-Mead',
                     options={'maxiter': 1000})
    
    return {
        'rho_s': result.x[0],
        'r_s': result.x[1],
        'chi2': result.fun,
        'success': result.success
    }


if __name__ == '__main__':
    print("Parsing SPARC catalog...")
    galaxies = parse_sparc_mrt('data/sparc/MassModels_Lelli2016c.mrt')
    print(f"Found {len(galaxies)} galaxies")
    
    pm = PushingMedium()
    results = {}
    
    for i, (name, data) in enumerate(galaxies.items(), 1):
        if len(data['R_kpc']) < 3:
            continue
        
        print(f"[{i}/{len(galaxies)}] Fitting {name}...")
        
        try:
            pm_fit = fit_galaxy_pm(data, pm)
            nfw_fit = fit_galaxy_nfw(data)
            
            results[name] = {
                'pm': pm_fit,
                'nfw': nfw_fit,
                'n_points': len(data['R_kpc'])
            }
            
            print(f"  PM χ²={pm_fit['chi2']:.2f}, NFW χ²={nfw_fit['chi2']:.2f}")
        except Exception as e:
            print(f"  Failed: {e}")
    
    # Save results
    with open('sparc_full_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nCompleted {len(results)} galaxies. Results saved to sparc_full_results.json")
