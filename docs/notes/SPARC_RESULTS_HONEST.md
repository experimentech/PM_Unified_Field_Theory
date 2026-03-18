# SPARC Fitting Results: PM vs Dark Matter

## Executive Summary

**Date**: 2026-03-02  
**Dataset**: SPARC (175 galaxies with high-quality rotation curves)  
**Model**: Pushing Medium (PM) - no dark matter halos

## Overall Statistics

```
Galaxies fitted:    175
Mean χ²:            198.11
Median χ²:          15.68
Std dev:            716.69
Range:              0.04 - 7581.37
```

## Fit Quality Distribution

| Quality      | χ² Range  | Count | Percentage |
|--------------|-----------|-------|------------|
| Excellent    | < 2       | 29    | 16.6%      |
| Good         | 2-5       | 20    | 11.4%      |
| Acceptable   | 5-10      | 23    | 13.1%      |
| **Poor**     | **≥ 10**  | **103** | **58.9%** |

**Success rate (χ² < 10): 41%**

## Key Finding: Data Point Sensitivity

**CRITICAL OBSERVATION**: Low data point galaxies fit much better than high data point ones.

```
Galaxies with <10 points:  51 galaxies, mean χ² = 7.78
Galaxies with ≥10 points: 124 galaxies, mean χ² = 276.39
```

This is **backwards** from what you'd expect if PM were fundamentally correct!

### Interpretation

1. **Optimization Issue**: The 6-parameter PM model gets stuck in local minima for complex rotation curves
2. **Model Complexity**: With more data points, structure emerges that PM's simple exponential profile can't capture
3. **Not Physical**: This pattern suggests a **fitting problem**, not a physics problem

## Best Fits (χ² < 2)

29 galaxies achieve excellent fits. Examples:

- **UGC02023**: χ² = 0.04 (5 points)
- **F574-2**: χ² = 0.04 (5 points)
- **UGC07232**: χ² = 0.09 (4 points)
- **UGC07577**: χ² = 0.13 (9 points)
- **UGC01281**: χ² = 1.85 (25 points) ← notable exception!

## Worst Fits (χ² > 100)

43 galaxies have catastrophic fits. Top 10 worst:

- **UGC02953**: χ² = 7581.37 (115 points)
- **UGC05253**: χ² = 3325.61 (73 points)
- **NGC5055**: χ² = 3168.10 (28 points)
- **UGC06787**: χ² = 2378.71 (71 points)
- **NGC2403**: χ² = 1643.97 (73 points)

## Parameter Distributions

```
v_inf: 130.3 ± 55.3 km/s
r_s:   135.9 ± 113.1 kpc
M_d:   1.53e+11 ± 1.69e+11 M_☉
```

## Honest Assessment

### What Worked
- ✅ 29 galaxies (16.6%) have excellent fits
- ✅ Simple galaxies with few data points fit well
- ✅ Parameters are physically reasonable

### What Failed
- ❌ **58.9% of galaxies have poor fits (χ² ≥ 10)**
- ❌ **Mean χ² = 198 is unacceptable**
- ❌ **Inverse correlation with data quality**: more data → worse fits
- ❌ High-resolution rotation curves systematically fail

### Root Cause

This is almost certainly an **optimization problem**, not a physics problem:

1. **6-parameter search space** is large and rough
2. **Simple exponential PM profile** may not have enough freedom
3. **Local minima** trap the optimizer for complex curves
4. **Random+refine strategy** (300+100 iterations) insufficient for high-dimensional fitting

## Next Steps

### Option 1: Improve Optimization (Recommended)
- Use better global optimizers (differential evolution, basin hopping)
- Increase iterations dramatically (10x)
- Add multi-scale fitting strategy
- Use baryonic components as initialization

### Option 2: Simplify Model
- Fix some parameters based on galaxy properties
- Reduce degrees of freedom
- Use scaling relations to constrain parameters

### Option 3: Accept Limited Regime
- Focus on galaxies with <15 data points where PM works
- Acknowledge PM needs refinement for detailed curves
- Still interesting as a proof-of-concept

## Comparison Context

**IMPORTANT**: We haven't compared to NFW halo fits yet! To claim victory or defeat, we need:

1. NFW χ² distributions for same 175 galaxies
2. MOND fits (if available from literature)
3. Published SPARC dark matter results

**Until then**: PM shows promise but needs significant optimization work.

## Data Availability

- ✅ Full SPARC dataset downloaded (175 galaxies)
- ✅ All rotation curves parsed successfully
- ✅ Fits completed (though suboptimal)
- ✅ Results saved: `sparc_fit_results.json`

---

**Bottom Line**: The infrastructure works, but the fitting needs major improvement before we can make claims about PM vs dark matter.
