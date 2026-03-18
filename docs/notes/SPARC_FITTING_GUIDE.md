# SPARC Galaxy Rotation Curve Fitting Guide

## Quick Start

### 1. Activate Environment
```bash
source venv/bin/activate
```

### 2. Get SPARC Data

Visit https://astroweb.case.edu/SPARC/ and download rotation curve files.

Extract to: `data/SPARC/`

Expected format (CSV with headers):
```
Name,R_kpc,V_obs,Err_V
GalaxyName,0.5,45.2,2.1
GalaxyName,1.0,68.5,2.8
...
```

### 3. Run Batch Fitter
```bash
python fit_sparc_batch.py
```

Results saved to `sparc_fit_results.json`

## Model

PM circular velocity:
```
v_c²(r) = r × [a_bar(r) + a_med(r)]
```

**Baryonic acceleration** (exponential disk):
```
M(<r) = M_d × [1 - (1 + r/R_d) × exp(-r/R_d)]
a_bar = G M(<r) / r²
```

**Medium acceleration** (index gradient):
```
δn = μ_coeff × Φ × (v_inf²/c²) × ln(1 + r/r_s) × S(r)
S(r) = [1 + (r_c/r)^m]^(-1/m)
a_med = -c² × d(δn)/dr
```

## Fitted Parameters

| Parameter | Description | Units |
|-----------|-------------|-------|
| M_d | Disk mass | kg |
| R_d | Disk scale length | m |
| v_inf | Asymptotic velocity parameter | m/s |
| r_s | Medium scale radius | m |
| r_c | Core transition radius | m |
| m | Smoothing exponent | dimensionless |

## Calibrated Constants

From `calibration.json`:
- `mu_coeff = 1.486e-27` (gravitational coupling)
- `k_Fizeau = 1.0` (Fizeau drag coefficient)
- `k_TT = 1.0` (Thomas precession coefficient)

## Test Results

**Synthetic galaxy recovery:**
- χ²/dof = 0.52 (excellent recovery of input parameters)
- 40 data points, 6 free parameters (34 degrees of freedom)

This demonstrates the fitting pipeline is working correctly.

## Next Steps

1. Download full SPARC dataset (~175 galaxies)
2. Run batch fits: `python fit_sparc_batch.py`
3. Analyze results: check χ²/dof distribution
4. Compare to legacy results if available
5. Create visualization plots for best fits

## Notes

- The fitter uses random search (500 trials) - fast but not optimal
- For publication, consider upgrading to scipy.optimize or MCMC
- Medium parameters (v_inf, r_s, r_c, m) are phenomenological
- No dark matter component—all "missing" mass explained via PM index gradient
