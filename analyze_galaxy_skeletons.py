"""
Skeleton analysis for PM galaxy refractive index fields.

Purpose:
    Analyze the topological structure of n(r) fields for fitted galaxies
    to understand why some fits succeed and others fail.
    
Key Idea:
    - Simple galaxies → simple skeleton topology (few critical points)
    - Complex galaxies → complex skeleton topology (many ridges/saddles)
    
This may reveal that fit quality correlates with skeleton complexity,
suggesting we need multi-component models for complex topologies.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d
import json
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from pushing_medium.core import unified_index_field


def compute_gradients_2d(field, dx, dy):
    """Compute ∂n/∂x and ∂n/∂y using central differences."""
    nx = np.zeros_like(field)
    ny = np.zeros_like(field)
    
    # x-direction
    nx[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / (2*dx)
    nx[:, 0] = (field[:, 1] - field[:, 0]) / dx
    nx[:, -1] = (field[:, -1] - field[:, -2]) / dx
    
    # y-direction
    ny[1:-1, :] = (field[2:, :] - field[:-2, :]) / (2*dy)
    ny[0, :] = (field[1, :] - field[0, :]) / dy
    ny[-1, :] = (field[-1, :] - field[-2, :]) / dy
    
    return nx, ny


def compute_hessian_2d(field, dx, dy):
    """Compute Hessian components n_xx, n_xy, n_yy."""
    nxx = np.zeros_like(field)
    nyy = np.zeros_like(field)
    
    # n_xx
    nxx[:, 1:-1] = (field[:, 2:] - 2*field[:, 1:-1] + field[:, :-2]) / (dx*dx)
    nxx[:, 0] = nxx[:, 1]
    nxx[:, -1] = nxx[:, -2]
    
    # n_yy
    nyy[1:-1, :] = (field[2:, :] - 2*field[1:-1, :] + field[:-2, :]) / (dy*dy)
    nyy[0, :] = nyy[1, :]
    nyy[-1, :] = nyy[-2, :]
    
    # n_xy via successive differentiation
    nx, _ = compute_gradients_2d(field, dx, dy)
    _, nxy = compute_gradients_2d(nx, dx, dy)
    
    return nxx, nxy, nyy


def analyze_skeleton_topology(field, dx, dy):
    """
    Analyze topological complexity of a 2D refractive index field.
    
    Returns:
        dict with:
            - n_maxima: number of local maxima
            - n_minima: number of local minima
            - n_saddles: number of saddle points
            - ridge_length: total length of ridge structures (pixels)
            - valley_length: total length of valley structures (pixels)
            - complexity_score: combined metric
    """
    # Smooth to stabilize derivatives
    field_smooth = gaussian_filter(field, sigma=1.5)
    
    # Gradients and Hessian
    nx, ny = compute_gradients_2d(field_smooth, dx, dy)
    nxx, nxy, nyy = compute_hessian_2d(field_smooth, dx, dy)
    
    grad_mag = np.hypot(nx, ny)
    
    # Hessian eigenvalues
    trace = nxx + nyy
    diff = nxx - nyy
    disc = np.sqrt(0.25*diff*diff + nxy*nxy)
    lam1 = 0.5*trace + disc  # larger eigenvalue
    lam2 = 0.5*trace - disc  # smaller eigenvalue
    
    # Find stationary points (low gradient)
    g_thresh = np.percentile(grad_mag, 1.0)
    stationary = grad_mag < g_thresh
    
    # Classify by Hessian eigenvalues
    n_maxima = np.sum(stationary & (lam1 < 0) & (lam2 < 0))
    n_minima = np.sum(stationary & (lam1 > 0) & (lam2 > 0))
    n_saddles = np.sum(stationary & (lam1 * lam2 < 0))
    
    # Ridge/valley detection
    # Ridges: negative curvature along principal direction
    ridge_mask = (lam1 < -0.01)  # threshold to avoid noise
    ridge_length = np.sum(ridge_mask)
    
    # Valleys: positive curvature
    valley_mask = (lam2 > 0.01)
    valley_length = np.sum(valley_mask)
    
    # Complexity score: weighted combination
    complexity = (
        n_maxima * 10 +      # Each maximum adds structure
        n_minima * 5 +       # Minima less important usually
        n_saddles * 15 +     # Saddles indicate topological complexity
        ridge_length / 100 + # Ridge extent
        valley_length / 100
    )
    
    return {
        'n_maxima': int(n_maxima),
        'n_minima': int(n_minima),
        'n_saddles': int(n_saddles),
        'ridge_length': int(ridge_length),
        'valley_length': int(valley_length),
        'complexity_score': float(complexity)
    }


def build_2d_index_field(params, r_max=30.0, n_grid=200):
    """
    Build 2D axisymmetric refractive index field n(x,y) from PM parameters.
    
    Assumes circular symmetry: n(x,y) = n(r) where r = sqrt(x² + y²)
    """
    x = np.linspace(-r_max, r_max, n_grid)
    y = np.linspace(-r_max, r_max, n_grid)
    X, Y = np.meshgrid(x, y)
    R = np.hypot(X, Y)
    
    # Extract galaxy parameters (in SI units already from fits)
    M_d = params['M_d']
    R_d = params['R_d']
    
    # Flatten to compute n for each radius  
    r_flat = R.flatten()
    n_flat = np.zeros_like(r_flat)
    
    for i, ri_kpc in enumerate(r_flat):
        if ri_kpc < 1e-3:
            ri_kpc = 1e-3  # avoid singularity
        # Convert kpc to meters
        ri_m = ri_kpc * 3.086e19  # kpc to m
        # Compute M(<r) for exponential disk
        x_norm = ri_m / R_d
        M_enclosed = M_d * (1 - (1 + x_norm) * np.exp(-x_norm))
        # Single point mass approximation for index calculation
        masses = [(M_enclosed, (0, 0, 0))]
        n_flat[i] = unified_index_field((ri_m, 0, 0), masses=masses)
    
    n_field = n_flat.reshape(R.shape)
    
    return X, Y, n_field, x, y


def analyze_galaxy_skeleton(galaxy_name, params):
    """Analyze skeleton topology for a specific galaxy."""
    print(f"\nAnalyzing skeleton for {galaxy_name}...")
    
    # Build 2D field
    X, Y, n_field, x_arr, y_arr = build_2d_index_field(params)
    dx = x_arr[1] - x_arr[0]
    dy = y_arr[1] - y_arr[0]
    
    # Analyze topology
    topology = analyze_skeleton_topology(n_field, dx, dy)
    
    print(f"  Critical points: {topology['n_maxima']} max, "
          f"{topology['n_minima']} min, {topology['n_saddles']} saddles")
    print(f"  Ridge/valley lengths: {topology['ridge_length']}, {topology['valley_length']} pixels")
    print(f"  Complexity score: {topology['complexity_score']:.1f}")
    
    return topology, (X, Y, n_field)


def compare_good_vs_bad_galaxies(results_file='data/pm_sparc_fits.json'):
    """
    Compare skeleton complexity between good fits and poor fits.
    
    Hypothesis: Good fits have simpler skeleton topology.
    """
    with open(results_file) as f:
        results = json.load(f)
    
    # Categorize galaxies
    good_galaxies = []
    bad_galaxies = []
    
    for name, result in results.items():
        chi2 = result['chi2']
        if chi2 < 5.0:
            good_galaxies.append((name, result))
        elif chi2 > 50.0:
            bad_galaxies.append((name, result))
    
    # Sample a few from each category
    print(f"\nFound {len(good_galaxies)} good fits and {len(bad_galaxies)} poor fits")
    print("\n" + "="*70)
    print("ANALYZING SKELETON COMPLEXITY")
    print("="*70)
    
    good_complexities = []
    bad_complexities = []
    
    # Analyze 5 from each category
    print("\n--- GOOD FITS (χ² < 5) ---")
    for name, result in good_galaxies[:5]:
        disk_params = result['disk']
        medium_params = result['medium']
        params = {**disk_params, **medium_params}
        topology, _ = analyze_galaxy_skeleton(name, params)
        good_complexities.append(topology['complexity_score'])
    
    print("\n--- POOR FITS (χ² > 50) ---")
    for name, result in bad_galaxies[:5]:
        disk_params = result['disk']
        medium_params = result['medium']
        params = {**disk_params, **medium_params}
        topology, _ = analyze_galaxy_skeleton(name, params)
        bad_complexities.append(topology['complexity_score'])
    
    # Statistics
    print("\n" + "="*70)
    print("SKELETON COMPLEXITY COMPARISON")
    print("="*70)
    print(f"Good fits (χ² < 5):")
    print(f"  Mean complexity: {np.mean(good_complexities):.1f}")
    print(f"  Std: {np.std(good_complexities):.1f}")
    print(f"\nPoor fits (χ² > 50):")
    print(f"  Mean complexity: {np.mean(bad_complexities):.1f}")
    print(f"  Std: {np.std(bad_complexities):.1f}")
    
    if np.mean(bad_complexities) > np.mean(good_complexities):
        print(f"\n✓ Hypothesis SUPPORTED: Poor fits show {np.mean(bad_complexities)/np.mean(good_complexities):.1f}x higher complexity")
        print("  → Model insufficient for complex skeleton topologies")
    else:
        print(f"\n✗ Hypothesis NOT supported: Poor fits not systematically more complex")
        print("  → Problem likely elsewhere (data quality, optimizer, etc.)")


def visualize_skeleton_example(galaxy_name, params, output_file=None):
    """Create detailed skeleton visualization for one galaxy."""
    X, Y, n_field, x, y = build_2d_index_field(params, r_max=40.0, n_grid=300)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    
    # Smooth and compute derivatives
    n_smooth = gaussian_filter(n_field, sigma=1.5)
    nx, ny = compute_gradients_2d(n_smooth, dx, dy)
    nxx, nxy, nyy = compute_hessian_2d(n_smooth, dx, dy)
    
    # Eigenvalues
    trace = nxx + nyy
    diff = nxx - nyy
    disc = np.sqrt(0.25*diff*diff + nxy*nxy)
    lam1 = 0.5*trace + disc
    lam2 = 0.5*trace - disc
    
    # Ridge mask
    ridge_mask = (lam1 < -1e-6) & (np.abs(lam1) > np.abs(lam2))
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: n(r) field with skeleton overlay
    ax = axes[0]
    im = ax.contourf(X, Y, n_field, levels=40, cmap='viridis')
    
    # Skeleton overlay
    skel_y, skel_x = np.where(ridge_mask)
    if len(skel_x) > 0:
        ax.plot(x[skel_x], y[skel_y], 'r.', ms=0.5, alpha=0.8, label='Ridge skeleton')
    
    ax.set_title(f"{galaxy_name}: Refractive Index n(r)", fontsize=12, weight='bold')
    ax.set_xlabel("x (kpc)")
    ax.set_ylabel("y (kpc)")
    ax.set_aspect('equal')
    fig.colorbar(im, ax=ax, label="n")
    ax.legend()
    
    # Right: Hessian eigenvalue structure
    ax = axes[1]
    # Color by curvature type
    curvature_type = np.zeros_like(lam1)
    curvature_type[ridge_mask] = -1  # ridges
    curvature_type[(lam2 > 1e-6)] = 1  # valleys/bowls
    curvature_type[(lam1 * lam2 < 0)] = 0  # saddles
    
    im2 = ax.contourf(X, Y, curvature_type, levels=[-1.5, -0.5, 0.5, 1.5], 
                     cmap='RdBu', alpha=0.7)
    ax.contour(X, Y, n_field, levels=20, colors='black', linewidths=0.3, alpha=0.4)
    
    ax.set_title(f"Curvature Topology", fontsize=12, weight='bold')
    ax.set_xlabel("x (kpc)")
    ax.set_ylabel("y (kpc)")
    ax.set_aspect('equal')
    cbar = fig.colorbar(im2, ax=ax, ticks=[-1, 0, 1])
    cbar.ax.set_yticklabels(['Ridge', 'Saddle', 'Valley'])
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"  Saved to {output_file}")
    else:
        plt.show()
    
    plt.close()


def main():
    """Run skeleton analysis on sample galaxies."""
    print("="*70)
    print("SKELETON TOPOLOGY ANALYSIS FOR PM GALAXY FITS")
    print("="*70)
    print("\nHypothesis: Complex skeleton topology correlates with poor fit quality")
    print("Expected: Good fits → simple topology (few critical points)")
    print("          Poor fits → complex topology (many critical points)")
    
    # Run comparison
    compare_good_vs_bad_galaxies()
    
    print("\n" + "="*70)
    print("Analysis complete!")
    print("\nNext: If hypothesis confirmed, need multi-component disk models")
    print("      to capture complex skeleton topologies.")
    print("="*70)


if __name__ == "__main__":
    main()
