#!/usr/bin/env python3
"""
Convert SPARC MassModels file to format compatible with load_sparc_real.

Input: MassModels_Lelli2016c.mrt (3416 lines)
Output: sparc_full.csv (compatible with your loader)

Format:
  Galaxy,Distance_Mpc,R_kpc,V_obs_kms,V_err_kms,V_gas_kms,V_disk_kms,V_bul_kms
"""

import sys
from pathlib import Path

def convert_sparc_mrt_to_csv(input_file, output_file):
    """Convert SPARC .mrt format to CSV."""
    
    print(f"Reading {input_file}...")
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header (find first data line)
    start_idx = 0
    for i, line in enumerate(lines):
        # Data lines start with galaxy name (letters/numbers)
        if line.strip() and line[0].isalpha():
            start_idx = i
            break
    
    print(f"Found data starting at line {start_idx + 1}")
    
    # Parse data
    galaxies = set()
    data_lines = []
    
    for line in lines[start_idx:]:
        if not line.strip() or line.startswith('#'):
            continue
        
        parts = line.split()
        if len(parts) < 9:
            continue
        
        try:
            galaxy = parts[0]
            D = float(parts[1])      # Mpc
            R = float(parts[2])      # kpc
            Vobs = float(parts[3])   # km/s
            eVobs = float(parts[4])  # km/s
            Vgas = float(parts[5])   # km/s
            Vdisk = float(parts[6])  # km/s
            Vbul = float(parts[7])   # km/s
            
            galaxies.add(galaxy)
            data_lines.append(f"{galaxy},{R},{Vobs},{eVobs},{D},{Vgas},{Vdisk},{Vbul}")
        except (ValueError, IndexError):
            continue
    
    print(f"✓ Parsed {len(data_lines)} data points")
    print(f"✓ Found {len(galaxies)} unique galaxies")
    
    # Write CSV
    with open(output_file, 'w') as f:
        f.write("Name,R_kpc,V_obs,Err_V,Dist_Mpc,V_gas,V_disk,V_bulge\n")
        for line in data_lines:
            f.write(line + "\n")
    
    print(f"✓ Written: {output_file}")
    
    # Summary
    print(f"\nGalaxy names (first 20):")
    for name in sorted(galaxies)[:20]:
        print(f"  {name}")
    if len(galaxies) > 20:
        print(f"  ... and {len(galaxies) - 20} more")
    
    return len(galaxies), len(data_lines)

if __name__ == '__main__':
    input_file = sys.argv[1] if len(sys.argv) > 1 else 'data/MassModels_Lelli2016c.mrt'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'data/sparc_full.csv'
    
    n_galaxies, n_points = convert_sparc_mrt_to_csv(input_file, output_file)
    
    print(f"\n{'='*80}")
    print(f"CONVERSION COMPLETE")
    print(f"{'='*80}")
    print(f"  Galaxies: {n_galaxies}")
    print(f"  Data points: {n_points}")
    print(f"  Output: {output_file}")
    print(f"\nReady to run:")
    print(f"  source venv/bin/activate")
    print(f"  python3 batch_fit_sparc.py data/sparc_full.csv")
    print(f"{'='*80}")
