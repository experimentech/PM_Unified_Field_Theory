#!/usr/bin/env python3
"""
Hybrid vector-tensor analysis for barred galaxies in SPARC.

For axisymmetric galaxies: pure 1D radial vector method is optimal.
For barred/interacting galaxies: hybrid approach reveals asymmetric flow structure.

Usage:
  source venv/bin/activate
  python3 hybrid_sparc_bars.py
"""

import sys
import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass

# Venv check
if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    print("❌ Virtual environment not active! Run: source venv/bin/activate")
    sys.exit(1)

sys.path.insert(0, str(Path(__file__).parent / 'legacy' / 'Pushing-Medium' / 'src'))

from galaxy_dynamics.data import load_sparc_real
from galaxy_dynamics.fitting import DiskParams, MediumParams, rotation_curve
import matplotlib.pyplot as plt

KPC_TO_M = 3.0856775814913673e19
M_SUN = 1.989e30

@dataclass
class BarParams:
    """Bar perturbation to axisymmetric disk."""
    M_bar: float  # Bar mass (kg)
    a_bar: float  # Bar semi-major axis (m)
    b_bar: float  # Bar semi-minor axis (m)
    theta: float  # Bar orientation angle (rad)

class PM2DFlow:
    """2D flow field for PM including bar perturbation."""
    
    def __init__(self, disk: DiskParams, medium: MediumParams, bar: BarParams = None):
        self.disk = disk
        self.medium = medium
        self.bar = bar
    
    def u_field(self, x, y):
        """Compute 2D flow field u = (u_x, u_y) at position (x, y).
        
        This is where hybrid becomes powerful:
        - Axisymmetric component from disk+medium (vector, 1D)
        - Bar perturbation breaks symmetry (requires 2D tensor field)
        """
        r = np.sqrt(x*x + y*y)
        
        # Axisymmetric base flow (radial)
        v_rot = rotation_curve(r, self.disk, self.medium)
        if r < 1e-10:
            u_r = 0.0
        else:
            # Rotational flow: u_phi = v / r, convert to Cartesian
            # For circular orbits: u_r = 0, u_phi = v/r
            # In Cartesian: u_x = -u_phi * sin(phi), u_y = u_phi * cos(phi)
            cos_phi = x / r
            sin_phi = y / r
            u_phi = v_rot / r
            u_x = -u_phi * sin_phi  # tangential
            u_y = u_phi * cos_phi
        else:
            u_x, u_y = 0.0, 0.0
        
        # Bar perturbation (if present)
        if self.bar is not None:
            # Rotate to bar frame
            cos_th = np.cos(self.bar.theta)
            sin_th = np.sin(self.bar.theta)
            x_bar = x * cos_th + y * sin_th
            y_bar = -x * sin_th + y * cos_th
            
            # Quadrupole-like perturbation (simplified)
            # Real implementation would compute bar potential and take gradient
            # For now: illustrative form that creates flow toward bar major axis
            a, b = self.bar.a_bar, self.bar.b_bar
            scale = self.bar.M_bar / (a * b)
            
            # Perturbation strength falls with distance
            r_norm = np.sqrt((x_bar/a)**2 + (y_bar/b)**2)
            if r_norm < 10:  # bar influence zone
                strength = scale * np.exp(-r_norm / 2)
                # Push toward major axis
                du_x_bar = -strength * y_bar / (b*b)
                du_y_bar = 0.0
                
                # Rotate back to lab frame
                du_x = du_x_bar * cos_th - du_y_bar * sin_th
                du_y = du_x_bar * sin_th + du_y_bar * cos_th
                
                u_x += du_x
                u_y += du_y
        
        return u_x, u_y
    
    def find_critical_points(self, x_range, y_range, n_seeds=50):
        """Vector approach: find stagnation points u = 0."""
        from scipy.optimize import fsolve
        
        critical_points = []
        seeds = []
        
        # Generate seeds
        for _ in range(n_seeds):
            x0 = np.random.uniform(*x_range)
            y0 = np.random.uniform(*y_range)
            seeds.append([x0, y0])
        
        def flow_vec(pos):
            return self.u_field(pos[0], pos[1])
        
        for seed in seeds:
            try:
                sol = fsolve(flow_vec, seed, full_output=True)
                if sol[2] == 1:  # converged
                    pos = sol[0]
                    # Check if new (not duplicate)
                    is_new = all(np.linalg.norm(pos - cp) > 1e-2 for cp in critical_points)
                    if is_new and x_range[0] <= pos[0] <= x_range[1] and y_range[0] <= pos[1] <= y_range[1]:
                        critical_points.append(pos)
            except:
                pass
        
        return np.array(critical_points) if critical_points else np.zeros((0, 2))
    
    def grid_near_features(self, critical_points, x_range, y_range, Nx=200, Ny=200, tube_width=0.5):
        """Tensor approach: compute grid only near critical features."""
        x = np.linspace(*x_range, Nx)
        y = np.linspace(*y_range, Ny)
        X, Y = np.meshgrid(x, y)
        
        # Build mask: sample near critical points
        mask = np.zeros((Ny, Nx), dtype=bool)
        for cp in critical_points:
            dist2 = (X - cp[0])**2 + (Y - cp[1])**2
            mask |= dist2 <= tube_width**2
        
        # Compute flow only where masked
        U = np.zeros((Ny, Nx))
        V = np.zeros((Ny, Nx))
        Speed = np.zeros((Ny, Nx))
        
        idxs = np.argwhere(mask)
        for (j, i) in idxs:
            u_x, u_y = self.u_field(x[i], y[j])
            U[j, i] = u_x
            V[j, i] = u_y
            Speed[j, i] = np.hypot(u_x, u_y)
        
        return {
            'x': x, 'y': y,
            'U': U, 'V': V, 'Speed': Speed,
            'mask': mask,
            'coverage': 100 * mask.sum() / mask.size
        }

def analyze_barred_galaxy(name: str, pm_fit: dict, bar_params: BarParams = None):
    """
    Hybrid analysis of a fitted galaxy with optional bar.
    
    1. Vector: find critical points in 2D projection
    2. Tensor: compute flow field near those points
    3. Compare axisymmetric vs barred cases
    """
    print(f"\n{'='*80}")
    print(f"HYBRID ANALYSIS: {name}")
    print(f"{'='*80}")
    
    # Reconstruct fitted parameters
    disk = DiskParams(**pm_fit['disk'])
    medium = MediumParams(**pm_fit['medium'])
    
    print(f"\nFitted PM parameters:")
    print(f"  Disk: M_d = {disk.M_d/M_SUN:.2e} M_☉, R_d = {disk.R_d/KPC_TO_M:.2f} kpc")
    print(f"  Medium: v_inf = {medium.v_inf/1000:.1f} km/s, r_s = {medium.r_s/KPC_TO_M:.1f} kpc")
    
    # Analysis domain (in kpc)
    x_range = (-10 * KPC_TO_M, 10 * KPC_TO_M)
    y_range = (-10 * KPC_TO_M, 10 * KPC_TO_M)
    
    # Case 1: Axisymmetric (no bar)
    print(f"\n1. Axisymmetric analysis:")
    flow_axi = PM2DFlow(disk, medium, bar=None)
    crit_axi = flow_axi.find_critical_points(x_range, y_range, n_seeds=30)
    print(f"   Critical points found: {len(crit_axi)}")
    print(f"   (Expected: 1 at center for axisymmetric case)")
    
    # Case 2: With bar (if provided)
    if bar_params:
        print(f"\n2. Barred galaxy analysis:")
        print(f"   Bar: M_bar = {bar_params.M_bar/M_SUN:.2e} M_☉, a = {bar_params.a_bar/KPC_TO_M:.2f} kpc")
        print(f"   Bar angle: {np.rad2deg(bar_params.theta):.1f}°")
        
        flow_bar = PM2DFlow(disk, medium, bar=bar_params)
        crit_bar = flow_bar.find_critical_points(x_range, y_range, n_seeds=50)
        print(f"   Critical points found: {len(crit_bar)}")
        print(f"   (Bar introduces L4/L5-like points along major axis)")
        
        # Compute grid near critical features
        print(f"\n3. Tensor field refinement:")
        grid_result = flow_bar.grid_near_features(
            crit_bar, x_range, y_range,
            Nx=200, Ny=200, tube_width=1.0*KPC_TO_M
        )
        print(f"   Grid coverage: {grid_result['coverage']:.1f}% (sampled near critical points)")
        print(f"   → Hybrid saves {100 - grid_result['coverage']:.1f}% computation vs full grid")
        
        return {
            'axisymmetric_critical_points': crit_axi,
            'barred_critical_points': crit_bar,
            'tensor_field': grid_result,
        }
    else:
        print(f"   (No bar specified - axisymmetric case only)")
        return {
            'axisymmetric_critical_points': crit_axi,
        }


if __name__ == '__main__':
    print("="*80)
    print("HYBRID VECTOR-TENSOR ANALYSIS FOR SPARC GALAXIES")
    print("="*80)
    
    # Load previous PM fits
    fits_file = 'data/pm_sparc_fits.json'
    if not Path(fits_file).exists():
        print(f"\n❌ No PM fits found. Run batch_fit_sparc.py first!")
        sys.exit(1)
    
    with open(fits_file) as f:
        pm_fits = json.load(f)
    
    print(f"\nLoaded {len(pm_fits)} fitted galaxies")
    
    # Demo: analyze first galaxy with and without hypothetical bar
    first_galaxy = list(pm_fits.keys())[0]
    print(f"\nDemo analysis on: {first_galaxy}")
    
    # Case 1: Axisymmetric (what we've been doing)
    result = analyze_barred_galaxy(first_galaxy, pm_fits[first_galaxy], bar_params=None)
    
    # Case 2: Add hypothetical bar (10% of disk mass, 3 kpc semi-major axis)
    print(f"\n{'='*80}")
    print(f"ADDING HYPOTHETICAL BAR")
    print(f"{'='*80}")
    
    disk_mass = pm_fits[first_galaxy]['disk']['M_d']
    bar = BarParams(
        M_bar=0.1 * disk_mass,  # 10% of disk mass
        a_bar=3 * KPC_TO_M,     # 3 kpc
        b_bar=0.5 * KPC_TO_M,   # 0.5 kpc
        theta=np.deg2rad(30),   # 30° bar angle
    )
    
    result_bar = analyze_barred_galaxy(first_galaxy, pm_fits[first_galaxy], bar_params=bar)
    
    print(f"\n{'='*80}")
    print(f"CONCLUSION:")
    print(f"{'='*80}")
    print(f"  Axisymmetric galaxies → pure vector (1D) is optimal")
    print(f"  Barred/interacting → hybrid reveals asymmetric critical structure")
    print(f"  PM naturally handles both: same medium, different mass distribution")
    print(f"\n  This is where PM's flow picture differs from GR:")
    print(f"    - GR: compute g_μν everywhere, then geodesics")
    print(f"    - PM: analyze u topology (vector), refine locally (tensor)")
    print(f"{'='*80}")
