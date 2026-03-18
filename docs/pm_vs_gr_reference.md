# PM vs GR Reference Values

Computed with `.venv` (PYTHONPATH=src). Inputs: M_sun=1.989e30 kg, b_sun=6.96e8 m, r1=1e11 m, r2=1.2e11 m, D_l=1e9 pc, D_s=2e9 pc, D_ls=D_s-D_l, M_cluster=1e14 M_sun.

| Quantity | PM | GR | Rel. diff |
| --- | --- | --- | --- |
| deflection_angle_point_mass | 8.488869e-06 | 8.488869e-06 | 0.000e+00 |
| shapiro_delay_point_mass (s) | 1.133570e-04 | 1.133570e-04 | 0.000e+00 |
| gravitational_redshift_potential (Δϕ=1e7) | 1.112650e-10 | 1.112650e-10 | 0.000e+00 |
| perihelion_precession (rad/orbit) | 5.020872e-07 | 5.020872e-07 | 0.000e+00 |
| lense_thirring_precession (rad/s) | 6.102516e-15 | 6.102516e-15 | 0.000e+00 |
| einstein_radius_point_mass (rad) | 9.784003e-05 | 9.784003e-05 | 0.000e+00 |
| gw_phase_speed (m/s) | 2.997925e+08 | 2.997925e+08 | 0.000e+00 |
| binary_quadrupole_power (W) | 1.460091e+25 | 1.460091e+25 | 0.000e+00 |
| circular_orbit_energy (J/kg) | -6.637591e+08 | -6.637591e+08 | -0.000e+00 |
| index_deflection_numeric (M_sun,b_sun,z_max=2e10,steps=2500,μ=2G/c^2) | 8.483710e-06 | 8.488869e-06 | -6.078e-04 |

Reproduce:

```bash
. .venv/bin/activate
PYTHONPATH=src python - <<'PY'
from pushing_medium import core as pm
from general_relativity import classical as gr

M_sun = 1.989e30
b_sun = 6.96e8
r1, r2 = 1e11, 1.2e11
D_l, D_s = 1e9*3.086e16, 2e9*3.086e16
D_ls = D_s - D_l
M_cluster = 1e14 * M_sun

cases = []
cases.append(("deflection_angle_point_mass", pm.pm_deflection_angle_point_mass(M_sun, b_sun), gr.deflection_angle_point_mass(M_sun, b_sun)))
cases.append(("shapiro_delay_point_mass", pm.pm_shapiro_delay_point_mass(M_sun, r1, r2, b_sun), gr.shapiro_delay_point_mass(M_sun, r1, r2, b_sun)))
cases.append(("gravitational_redshift_potential", pm.pm_gravitational_redshift_from_potential(1e7), gr.gravitational_redshift_potential(1e7)))
cases.append(("perihelion_precession", pm.pm_perihelion_precession(5.79e10, 0.2056, M_sun), gr.perihelion_precession(5.79e10, 0.2056, M_sun)))
cases.append(("lense_thirring_precession", pm.lense_thirring_precession(7.1e33, 1.2e7), gr.lense_thirring_precession(7.1e33, 1.2e7)))
cases.append(("einstein_radius_point_mass", pm.pm_einstein_radius_point_mass(M_cluster, D_l, D_s, D_ls), gr.einstein_radius_point_mass(M_cluster, D_l, D_s, D_ls)))
cases.append(("gw_phase_speed", pm.pm_gw_phase_speed(), gr.gw_phase_speed()))
cases.append(("binary_quadrupole_power", pm.pm_binary_quadrupole_power(1.4*M_sun, 1.3*M_sun, 1e9), gr.binary_quadrupole_power(1.4*M_sun, 1.3*M_sun, 1e9)))
cases.append(("circular_orbit_energy", pm.pm_circular_orbit_energy(M_sun, 1e11), gr.circular_orbit_energy(M_sun, 1e11)))

for name, pm_val, gr_val in cases:
    rel = (pm_val - gr_val) / gr_val if gr_val != 0 else 0.0
    print(f"{name}: PM={pm_val:.6e}, GR={gr_val:.6e}, rel_diff={rel:.3e}")

mu = 2*pm.G/(pm.c*pm.c)
numeric = pm.index_deflection_numeric(M_sun, b_sun, mu=mu, z_max=2e10, steps=2500)
analytic = pm.pm_deflection_angle_point_mass(M_sun, b_sun)
rel = (numeric - analytic)/analytic
print(f"index_deflection_numeric: numeric={numeric:.6e}, analytic={analytic:.6e}, rel_diff={rel:.3e}")
PY
```
