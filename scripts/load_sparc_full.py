#!/usr/bin/env python3
"""
Load full SPARC dataset from .dat rotation curve files.
"""
import sys
from pathlib import Path
import numpy as np

# Venv check
if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    print("❌ Virtual environment not active! Run: source venv/bin/activate")
    sys.exit(1)

sys.path.insert(0, str(Path(__file__).parent / 'legacy' / 'Pushing-Medium' / 'src'))

from galaxy_dynamics.data import RotationCurve

KPC_TO_M = 3.0856775814913673e19

def load_sparc_dat(filepath):
    """Load a single SPARC .dat file."""
    distance_mpc = None
    radii = []
    v_obs = []
    v_err = []
    v_gas = []
    v_disk = []
    v_bul = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Distance'):
                distance_mpc = float(line.split('=')[1].strip().split()[0])
            elif line.startswith('#') or not line:
                continue
            else:
                parts = line.split()
                if len(parts) >= 6:
                    radii.append(float(parts[0]))
                    v_obs.append(float(parts[1]))
                    v_err.append(float(parts[2]))
                    v_gas.append(float(parts[3]))
                    v_disk.append(float(parts[4]))
                    v_bul.append(float(parts[5]))
    
    radii = np.array(radii) * KPC_TO_M
    v_obs = np.array(v_obs) * 1000
    v_err = np.array(v_err) * 1000
    v_gas = np.array(v_gas) * 1000
    v_disk = np.array(v_disk) * 1000
    v_bul = np.array(v_bul) * 1000
    
    # Compute surface densities from velocity contributions
    # v² = G*M/r for a thin disk → Σ ∝ v²*r/G
    # Use placeholder values - the fitting code will work with components
    sigma_star = np.zeros(len(radii))
    sigma_gas = np.zeros(len(radii))
    
    components = {
        'gas': v_gas.tolist(),
        'disk': v_disk.tolist(),
        'bulge': v_bul.tolist(),
    }
    
    return RotationCurve(
        name=filepath.stem.replace('_rotmod', ''),
        radii_m=radii.tolist(),
        v_obs_ms=v_obs.tolist(),
        v_err_ms=v_err.tolist(),
        sigma_star=sigma_star.tolist(),
        sigma_gas=sigma_gas.tolist(),
        meta={'source': 'SPARC'},
        components=components,
        distance_mpc=distance_mpc,
    )

def load_all_sparc(data_dir='data'):
    """Load all SPARC rotation curves from .dat files."""
    data_path = Path(data_dir)
    dat_files = sorted(data_path.glob('*_rotmod.dat'))
    
    galaxies = {}
    for filepath in dat_files:
        try:
            rc = load_sparc_dat(filepath)
            galaxies[rc.name] = rc
        except Exception as e:
            print(f"⚠️  Failed to load {filepath.name}: {e}")
    
    return galaxies

if __name__ == '__main__':
    galaxies = load_all_sparc()
    print(f"Loaded {len(galaxies)} galaxies from SPARC dataset")
    
    # Show sample
    for i, (name, rc) in enumerate(list(galaxies.items())[:5]):
        print(f"\n{name}:")
        print(f"  Distance: {rc.distance_mpc:.2f} Mpc")
        print(f"  Data points: {len(rc.radii_m)}")
        print(f"  Radius range: {rc.radii_m[0]/KPC_TO_M:.2f}-{rc.radii_m[-1]/KPC_TO_M:.2f} kpc")
        print(f"  V_obs range: {rc.v_obs_ms[0]/1000:.1f}-{rc.v_obs_ms[-1]/1000:.1f} km/s")
