#!/usr/bin/env python3
"""Generate a realistic synthetic galaxy rotation curve for testing."""

import numpy as np
import csv
from pathlib import Path

# Constants
G = 6.67430e-11
c = 299792458.0
KPC_TO_M = 3.0856775814913673e19
KM_TO_M = 1000.0
MU_COEFF = 1.4861356300677034e-27

# True parameters (what we'll try to recover)
M_d_true = 1e41  # kg (~50 billion solar masses)
R_d_true = 10 * KPC_TO_M  # 10 kpc
v_inf_true = 200 * KM_TO_M  # 200 km/s
r_s_true = 15 * KPC_TO_M  # 15 kpc
r_c_true = 2 * KPC_TO_M  # 2 kpc
m_true = 2.0

def pm_velocity(r, M_d, R_d, v_inf, r_s, r_c, m):
    """Generate PM model velocity."""
    if r <= 0:
        return 0.0
    
    # Baryonic
    x = r / R_d
    M_enc = M_d * (1.0 - (1.0 + x) * np.exp(-x))
    a_bar = G * M_enc / (r * r)
    
    # Medium (simplified)
    S = (1.0 + (r_c / r)**m)**(-1.0/m)
    factor = (v_inf**2 / c**2) * np.log(1.0 + r/r_s) * S
    
    # Approximate gradient
    dr = r * 0.001
    r1, r2 = max(r - dr, 1e10), r + dr
    
    Phi1 = -G * M_d / r1 * np.exp(-r1/R_d)
    Phi2 = -G * M_d / r2 * np.exp(-r2/R_d)
    S1 = (1.0 + (r_c / r1)**m)**(-1.0/m)
    S2 = (1.0 + (r_c / r2)**m)**(-1.0/m)
    
    delta_n1 = MU_COEFF * Phi1 * ((v_inf**2 / c**2) * np.log(1.0 + r1/r_s) * S1)
    delta_n2 = MU_COEFF * Phi2 * ((v_inf**2 / c**2) * np.log(1.0 + r2/r_s) * S2)
    
    grad_n = (delta_n2 - delta_n1) / (2.0 * dr)
    a_med = -c**2 * grad_n
    
    a_total = a_bar + a_med
    return np.sqrt(r * a_total) if a_total > 0 else 0.0


# Generate rotation curve
radii_kpc = np.linspace(0.5, 30, 40)  # 40 points from 0.5 to 30 kpc
radii_m = radii_kpc * KPC_TO_M

velocities = []
for r in radii_m:
    v = pm_velocity(r, M_d_true, R_d_true, v_inf_true, r_s_true, r_c_true, m_true)
    # Add 3% noise
    v_noisy = v + np.random.normal(0, 0.03 * v)
    velocities.append(v_noisy)

velocities = np.array(velocities)
errors = velocities * 0.05  # 5% errors

# Write to CSV
output = Path('data/SPARC/TestGalaxy_PM.csv')
with open(output, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Name', 'R_kpc', 'V_obs', 'Err_V'])
    for r_kpc, v_ms, e_ms in zip(radii_kpc, velocities, errors):
        writer.writerow(['TestGalaxy_PM', r_kpc, v_ms/KM_TO_M, e_ms/KM_TO_M])

print(f"Generated: {output}")
print(f"True parameters:")
print(f"  M_d = {M_d_true:.3e} kg")
print(f"  R_d = {R_d_true/KPC_TO_M:.1f} kpc")
print(f"  v_inf = {v_inf_true/KM_TO_M:.0f} km/s")
print(f"  r_s = {r_s_true/KPC_TO_M:.1f} kpc")
print(f"  r_c = {r_c_true/KPC_TO_M:.1f} kpc")
print(f"  m = {m_true:.1f}")
