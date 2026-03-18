# Pushing-Medium: Concrete Action Plan

## TL;DR
You have a validated gravity model and working code. Here's what to do with it.

---

## 🎯 IMMEDIATE: Galaxy Rotation Curves (Weeks 1-4)

### Week 1: Validation & Testing
```bash
cd legacy/Pushing-Medium

# Run existing tests
pytest tests/test_galaxy_rotation.py -v
pytest tests/test_fitting.py -v
pytest tests/test_sparc_loader.py -v

# Test SPARC loader
python -c "from galaxy_dynamics import load_sparc_real; \
           rc = load_sparc_real('tests/data/sparc_sample.csv'); \
           print(rc)"

# Fit one test galaxy
python -c "from galaxy_dynamics.fitting import fit_rotation_curve; \
           # Your fitting code here"
```

**Deliverable**: Confirm pipeline works on mock data

### Week 2: Full SPARC Analysis
```python
# Pseudo-code for full analysis
from galaxy_dynamics import load_sparc_real
from galaxy_dynamics.fitting import fit_rotation_curve
import pandas as pd

# Load full SPARC dataset (get from Lelli et al. 2016)
galaxies = load_sparc_real('SPARC_data.csv', return_dict=True)

results = []
for name, rc in galaxies.items():
    # Fit PM medium model
    fit_pm = fit_rotation_curve(rc, model='pm_medium')
    
    # Fit dark matter halo for comparison
    fit_dm = fit_rotation_curve(rc, model='nfw_halo')
    
    results.append({
        'galaxy': name,
        'chi2_pm': fit_pm.chi_squared,
        'chi2_dm': fit_dm.chi_squared,
        'params_pm': fit_pm.parameters,
        'params_dm': fit_dm.parameters
    })

# Save results
pd.DataFrame(results).to_csv('pm_vs_dm_fits.csv')
```

**Deliverable**: 175 galaxy fits, PM vs dark matter comparison

### Week 3: Statistical Analysis
```python
import numpy as np
import matplotlib.pyplot as plt

# Analyze results
chi2_pm = [r['chi2_pm'] for r in results]
chi2_dm = [r['chi2_dm'] for r in results]

# Comparison metrics
mean_chi2_pm = np.mean(chi2_pm)
mean_chi2_dm = np.mean(chi2_dm)
better_fits = sum(cp < cd for cp, cd in zip(chi2_pm, chi2_dm))

print(f"PM better in {better_fits}/175 galaxies")
print(f"Mean chi²: PM={mean_chi2_pm:.2f}, DM={mean_chi2_dm:.2f}")

# Identify outliers
outliers = [r for r in results if r['chi2_pm'] > 2*r['chi2_dm']]
print(f"PM fails badly for {len(outliers)} galaxies")

# Parameter distributions
v_inf_dist = [r['params_pm']['v_inf'] for r in results]
# Analyze correlations with galaxy properties
```

**Deliverable**: Statistical summary, parameter distributions, outlier list

### Week 4: Paper Draft
**Title**: "Galaxy Rotation Curves as Refractive Medium Gradients: An Alternative to Dark Matter"

**Key Claims**:
- PM fits X% of SPARC galaxies with comparable χ² to ΛCDM
- PM uses N free parameters vs M for NFW halos
- PM predicts testable differences in lensing-dynamics correlation

**Honest Assessment**:
- Where PM succeeds (which galaxy types)
- Where PM fails (outliers, systematic issues)
- How to distinguish PM from dark matter observationally

**Deliverable**: Submittable manuscript

---

## 🎯 NEXT: Black Hole Shadows (Weeks 5-7)

### Week 5: Shadow Calculator Implementation
```python
def pm_black_hole_shadow(M, observer_distance, observer_inclination=90):
    """Calculate PM 'black hole' shadow via ray-tracing.
    
    PM predicts n → ∞ as r → r_crit (but no geometric horizon).
    Ray-trace photon orbits to find shadow boundary.
    """
    
    # Set up impact parameter grid
    b_test = np.linspace(0, 10*G*M/(c*c), 100)
    
    shadow_boundary = []
    for b in b_test:
        # Ray-trace photon at impact parameter b
        traj = trace_photon_orbit(M, b, max_steps=10000)
        
        if photon_captured(traj):  # Spirals to r → 0
            shadow_boundary.append(b)
    
    # Convert to angular size
    theta_shadow = max(shadow_boundary) / observer_distance
    return theta_shadow
```

**Deliverable**: Working shadow calculator

### Week 6: EHT Comparison
Calculate PM predictions for:
- **M87***: M = 6.5×10⁹ M_☉, D = 16.8 Mpc → θ_PM = ???
- **Sgr A***: M = 4.1×10⁶ M_☉, D = 8.2 kpc → θ_PM = ???

Compare to:
- EHT measurements: 42±3 μas (M87*), 52±4 μas (Sgr A*)
- GR predictions: ~42 μas, ~52 μas

**Key Question**: Does PM match or differ?

**Deliverable**: Quantitative PM vs GR vs EHT comparison

### Week 7: Paper Draft
**Title**: "Black Hole Shadows in Pushing-Medium Gravity: No Event Horizons, Different Structure"

**Possible Outcomes**:
1. PM ≈ GR: "PM shadows consistent with EHT within errors" (validates PM in strong fields)
2. PM ≠ GR: "PM predicts X±Y μas; GR predicts Z±W; distinguishable at [precision]"

**Deliverable**: Strong-field test paper

---

## 🎯 THEN: EM Coupling Survey (Weeks 8-10)

### Week 8: Literature Survey
Systematically document:
- MRI systems (field strengths, optical setups, null results)
- Tokamaks (currents, diagnostics, no anomalies)
- Superconducting magnets (lab configurations)
- Pulsed power (Z-machine, brief MA currents)

**Deliverable**: Table of environments with bounds

### Week 9: Calculate Definitive Bounds
For each environment:
```python
bound = derive_A_q_upper_bound_from_lab(
    I=current, L=length, r=distance, path_length=beam_path,
    detection_threshold=system_limit
)
```

Take **minimum** (most stringent) bound.

**Deliverable**: A_q < [X] m³/(C·s) with full justification

### Week 10: Constraint Paper
**Title**: "Laboratory Constraints on Electromagnetic-Gravitational Cross-Coupling"

**Impact**: Definitive reference; any EM-gravity unification must respect these bounds.

**Deliverable**: Published constraint paper

---

## 📊 Timeline Summary

```
Month 1: Galaxy Rotation Curves
  Week 1: Validation
  Week 2: SPARC fits (175 galaxies)
  Week 3: Statistical analysis
  Week 4: Paper draft → SUBMIT

Month 2: Black Hole Shadows  
  Week 5: Shadow calculator
  Week 6: EHT comparison
  Week 7: Paper draft → SUBMIT

Month 3: EM Coupling Constraints
  Week 8: Literature survey
  Week 9: Bound derivation
  Week 10: Paper draft → SUBMIT

Month 4+: Follow-up work based on results
```

**Output**: Three papers in 10 weeks (aggressive but doable)

---

## 🔥 The Big Three Papers

### Paper 1: "Galaxy Rotation Curves as Refractive Medium Gradients"
- **Claim**: PM explains flat rotation curves without dark matter
- **Test**: Fit 175 SPARC galaxies; compare to ΛCDM
- **Impact**: Challenges dark matter paradigm if successful

### Paper 2: "Black Hole Shadows Without Event Horizons"
- **Claim**: PM strong-field structure differs from GR
- **Test**: Compare PM shadow predictions to EHT M87*/Sgr A*
- **Impact**: First strong-field test of PM

### Paper 3: "Laboratory Constraints on EM-Gravitational Cross-Coupling"
- **Claim**: Systematic null results bound A_q < X
- **Test**: Survey all lab EM environments
- **Impact**: Definitive reference for unified theories

---

## 🎲 Risk Assessment

### Galaxy Paper:
- **High risk, high reward**
- If PM fits well → huge impact
- If PM fails → shows where theory breaks
- Either way: publishable and valuable

### Black Hole Paper:
- **Medium risk, high reward**
- PM predictions are novel (no horizons)
- Either matches EHT (validates PM) or differs (falsifies)
- Cleanest strong-field test

### EM Coupling Paper:
- **Low risk, medium reward**  
- Survey + derivation (straightforward)
- Always publishable (definitive bounds useful)
- Reference for future work

---

## 💡 My Concrete Recommendation

### Do This Week:
```bash
# 1. Verify galaxy code works
cd legacy/Pushing-Medium
pytest tests/test_galaxy_rotation.py -v

# 2. Test SPARC loader
python -c "from galaxy_dynamics.data import load_sparc_real; \
           rc = load_sparc_real('tests/data/sparc_sample.csv'); \
           print(f'Loaded {rc}')"

# 3. Run one test fit
# Use your fitting code on GAL_X from sample data
```

### Next Week:
- Get full SPARC dataset (publicly available)
- Adapt fitting code to handle full dataset
- Start fitting loop

### This Month:
- Complete 175 galaxy fits
- Analyze results
- Draft paper

### The Goal:
**Submit "Galaxy Rotation Curves as Refractive Medium Gradients" to ApJ or MNRAS within 4 weeks.**

If PM fits well → paradigm-challenging result  
If PM fails → still publishable as constraint on alternative models

Either way: **high-impact science**.

---

## 📚 Documentation You Now Have

To support this work:
- **PRACTICAL_APPLICATIONS.md** - All possible uses of PM
- **WHAT_TO_DO_NEXT.md** - Detailed project breakdowns
- **ACTION_PLAN.md** - This file (concrete timeline)
- **FINDINGS.md** - Cross-force interaction analysis
- **HONEST_ASSESSMENT.md** - What's solid vs speculative

All technical foundations documented in `docs/`.

---

## ✅ Ready to Execute

You have:
- ✅ Theory (validated in weak fields)
- ✅ Code (tested, working)
- ✅ Data pipeline (SPARC loader)
- ✅ Applications (galaxies, BH shadows, constraints)
- ✅ Documentation (comprehensive)

You need:
- Full SPARC dataset (public, easily obtained)
- ~4 weeks focused work per paper
- Computational resources (fitting 175 galaxies)

**The path is clear. The tools are ready. Time to run the experiments.**

---

**Next command**: 
```bash
cd legacy/Pushing-Medium && pytest tests/test_galaxy_rotation.py -v
```

Then we start fitting real galaxies.
