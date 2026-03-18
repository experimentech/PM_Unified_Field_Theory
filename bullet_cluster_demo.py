"""Bullet Cluster: GR vs Pushing-Medium flow comparison.

Shows how PM's moving-lens correction creates flow-based lensing patterns
that differ from static GR, particularly for the fast-moving subcluster.

This uses the hybrid approach: vector formalism for the deflection field,
tensor calibration for the moving-lens correction factor.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Activate venv check
venv_path = Path(__file__).parent / "venv"
if not (venv_path / "bin" / "python").exists():
    print("Virtual environment not found. Run: python3 -m venv venv && source venv/bin/activate")
    sys.exit(1)

sys.path.insert(0, str(Path(__file__).parent / "src"))

try:
    from pushing_medium import core as pm
except ImportError as e:
    print(f"Cannot import pushing_medium.core: {e}")
    print("Ensure src/pushing_medium/core.py exists")
    sys.exit(1)

# Physical constants
G = pm.G
c = pm.c
kpc = 3.085677581e19  # metres
M_SUN = 1.98847e30
ARCSEC_PER_RAD = (180.0 / np.pi) * 3600.0

# Calibration from previous work
MU_COEFF = 1.4861356300677034e-27


def plummer_surface_density(M, a_soft, X, Y, x0, y0):
    """Σ(R) = (M a²) / [π (R² + a²)²] for projected Plummer sphere."""
    dx, dy = X - x0, Y - y0
    R2 = dx*dx + dy*dy
    return (M * a_soft**2) / (np.pi * (R2 + a_soft**2)**2)


def plummer_deflection_gr(M, a_soft, X, Y, x0, y0):
    """GR weak-field deflection for Plummer lens."""
    dx, dy = X - x0, Y - y0
    R2 = dx*dx + dy*dy
    pref = 4.0 * G * M / (c**2)
    denom = R2 + a_soft**2
    return pref * dx / denom, pref * dy / denom


def pm_moving_lens_ratio(M, a_soft, v_trans, X, Y, x0, y0, samples=120):
    """Compute PM/GR deflection ratio for moving lens using hybrid approach."""
    if v_trans == 0.0:
        return np.ones_like(X)
    
    R_grid = np.hypot(X - x0, Y - y0)
    b_min = max(0.1 * kpc, float(np.min(R_grid[R_grid > 0])))
    b_max = float(np.max(R_grid))
    
    # Sample impact parameters
    bs = np.linspace(b_min, b_max, samples)
    ratios = np.ones_like(bs)
    z_max = 6.0 * a_soft
    steps = 1200
    
    for i, b in enumerate(bs):
        static_defl = pm.fermat_deflection_static_index(
            M, float(b), mu=MU_COEFF, z_max=z_max, steps=steps
        )
        moving_defl = pm.moving_lens_deflection_numeric(
            M, float(b), mu=MU_COEFF, v_transverse=v_trans, 
            z_max=z_max, steps=steps
        )
        ratios[i] = moving_defl / static_defl if static_defl > 0 else 1.0
    
    # Interpolate to grid
    return np.interp(R_grid.ravel(), bs, ratios, 
                     left=ratios[0], right=ratios[-1]).reshape(R_grid.shape)


def main():
    print("Generating Bullet Cluster comparison...")
    
    # Toy cluster parameters
    M_main = 1.0e15 * M_SUN
    M_sub = 3.0e14 * M_SUN
    a_main = 150.0 * kpc
    a_sub = 80.0 * kpc
    separation = 800.0 * kpc
    v_sub = 3000.0e3  # 3000 km/s
    
    # Grid
    nx = ny = 220
    field_size = 2200.0 * kpc
    xs = np.linspace(-field_size/2, field_size/2, nx)
    ys = np.linspace(-field_size/2, field_size/2, ny)
    XX, YY = np.meshgrid(xs, ys, indexing='xy')
    
    # Cluster centers
    x_main, y_main = -separation/2, 0.0
    x_sub, y_sub = +separation/2, 0.0
    
    # Surface density
    sigma_main = plummer_surface_density(M_main, a_main, XX, YY, x_main, y_main)
    sigma_sub = plummer_surface_density(M_sub, a_sub, XX, YY, x_sub, y_sub)
    sigma_total = sigma_main + sigma_sub
    
    # GR deflection (both static)
    ax_main, ay_main = plummer_deflection_gr(M_main, a_main, XX, YY, x_main, y_main)
    ax_sub, ay_sub = plummer_deflection_gr(M_sub, a_sub, XX, YY, x_sub, y_sub)
    alpha_x_gr = ax_main + ax_sub
    alpha_y_gr = ay_main + ay_sub
    
    # PM deflection (moving subcluster)
    print("Computing PM moving-lens correction...")
    ratio_grid = pm_moving_lens_ratio(M_sub, a_sub, v_sub, XX, YY, x_sub, y_sub)
    alpha_x_pm = ax_main + ax_sub * ratio_grid
    alpha_y_pm = ay_main + ay_sub * ratio_grid
    
    # Magnitudes
    mag_gr = np.hypot(alpha_x_gr, alpha_y_gr)
    mag_pm = np.hypot(alpha_x_pm, alpha_y_pm)
    mag_diff = mag_pm - mag_gr
    mag_frac = np.where(mag_gr > 0, mag_diff / mag_gr, 0.0)
    
    # Convert to arcsec
    mag_gr_arcsec = mag_gr * ARCSEC_PER_RAD
    mag_pm_arcsec = mag_pm * ARCSEC_PER_RAD
    mag_diff_mas = mag_diff * ARCSEC_PER_RAD * 1000.0
    
    # Plot
    extent = [xs[0]/kpc, xs[-1]/kpc, ys[0]/kpc, ys[-1]/kpc]
    fig, axs = plt.subplots(2, 2, figsize=(14, 12))
    
    # Panel 1: Surface density
    ax = axs[0, 0]
    im0 = ax.imshow(sigma_total / (M_SUN / kpc**2), origin='lower', extent=extent,
                    cmap='inferno', norm=LogNorm())
    ax.set_title('Projected Surface Density')
    fig.colorbar(im0, ax=ax, label='Σ [M☉/kpc²]')
    ax.plot([x_main/kpc], [y_main/kpc], 'o', color='cyan', ms=10, label='Main (static)')
    ax.plot([x_sub/kpc], [y_sub/kpc], 'o', color='magenta', ms=10, label='Sub (3000 km/s →)')
    ax.legend()
    
    # Panel 2: GR deflection (points style - discrete samples)
    ax = axs[0, 1]
    im1 = ax.imshow(np.maximum(mag_gr_arcsec, 1e-6), origin='lower', extent=extent,
                    cmap='gray', norm=LogNorm())
    ax.set_title('GR Deflection |α| (discrete point masses)')
    fig.colorbar(im1, ax=ax, label='|α| [arcsec]')
    # Sparse quiver to show "point-like" nature
    step = 15
    ax.quiver(XX[::step, ::step]/kpc, YY[::step, ::step]/kpc,
              alpha_x_gr[::step, ::step], alpha_y_gr[::step, ::step],
              color='white', alpha=0.6, scale=0.0001)
    
    # Panel 3: PM deflection (flow style - streamlines)
    ax = axs[1, 0]
    im2 = ax.imshow(np.maximum(mag_pm_arcsec, 1e-6), origin='lower', extent=extent,
                    cmap='viridis', norm=LogNorm())
    ax.set_title('PM Deflection |α| (flow field)')
    fig.colorbar(im2, ax=ax, label='|α| [arcsec]')
    # Streamlines show continuous flow
    speed_norm = mag_pm / (np.max(mag_pm) + 1e-30)
    lw = 0.6 + 1.8 * speed_norm
    ax.streamplot(XX/kpc, YY/kpc, alpha_x_pm, alpha_y_pm,
                  color='white', linewidth=lw, density=1.3, arrowsize=0.9)
    
    # Panel 4: Difference
    ax = axs[1, 1]
    im3 = ax.imshow(mag_frac * 100, origin='lower', extent=extent,
                    cmap='coolwarm', vmin=-1.0, vmax=1.0)
    ax.set_title('PM effect: (α_PM / α_GR − 1) [%]')
    fig.colorbar(im3, ax=ax, label='Δα / α_GR [%]')
    # Overlay difference vector field
    diff_x, diff_y = alpha_x_pm - alpha_x_gr, alpha_y_pm - alpha_y_gr
    ax.quiver(XX[::step, ::step]/kpc, YY[::step, ::step]/kpc,
              diff_x[::step, ::step], diff_y[::step, ::step],
              color='black', alpha=0.5, scale=0.00002)
    
    for ax in axs.flat:
        ax.set_xlabel('x [kpc]')
        ax.set_ylabel('y [kpc]')
        ax.set_aspect('equal')
    
    fig.suptitle('Bullet Cluster: GR (points) vs PM (flow)', fontsize=15, fontweight='bold')
    fig.tight_layout()
    
    outfile = Path(__file__).parent / "bullet_cluster_pm_vs_gr.png"
    plt.savefig(outfile, dpi=200)
    print(f"✓ Saved: {outfile}")
    print(f"  Max PM correction: {np.max(np.abs(mag_frac))*100:.3f}%")
    print(f"  Max difference: {np.max(mag_diff_mas):.2f} mas")


if __name__ == "__main__":
    main()
