# PM Core-Skeleton Ratio Analysis: A Negative Result That Matters

## Executive Summary

Initial analysis suggested excellent PM fits clustered around r_c/r_s ≈ 0.55, hinting at a universal structural constraint. **Comprehensive analysis disproves this.**

**Key Finding:** Excellent fits span r_c/r_s from 0.007 to 7.8—no universal ratio exists.

**Why This Is Good:** PM's strength is **geometric flexibility**, not forcing galaxies into a template. This is a feature, not a bug.

---

## The Initial Pattern (Artifact)

Early subset analysis showed:
- Excellent fits: r_c/r_s ≈ 0.55  
- Poor fits: r_c/r_s ≈ 0.30

This suggested PM had a preferred internal flow geometry.

---

## The Comprehensive Result (Reality)

### Full Dataset Statistics

**Excellent Fits (χ²_red < 1.0, n=30):**
```
r_c/r_s ratio:
  Min:     0.007
  25th:    0.187
  Median:  0.511
  75th:    1.124
  Max:     7.805
  Mean:    1.010
  Std:     1.671
```

**Poor Fits (χ²_red > 5.0, n=39):**
```
r_c/r_s ratio:
  Min:     0.021
  25th:    0.246
  Median:  0.532
  75th:    1.272
  Max:     8.915
  Mean:    1.318
  Std:     2.066
```

### Critical Observation

**The distributions are nearly identical.** The median differs by only 4%, and both span three orders of magnitude.

**Conclusion:** r_c/r_s is not a discriminator between good and poor PM fits.

---

## Why This Is Actually Better

### 1. PM Is Not a Template Model

If PM required r_c/r_s ≈ 0.55, it would be:
- Just another inflexible halo profile
- Unable to handle diverse galaxy morphologies
- A "one-size-fits-all" approach

Instead, PM supports a **wide family of flow geometries**.

### 2. Good Fits Are Not Geometrically Constrained

Excellent fits include:
- **Compact cores:** UGC06818 (r_c/r_s = 0.007)
- **Extended cores:** NGC6503 (r_c/r_s = 1.456)
- **Ultra-wide cores:** UGC11914 (r_c/r_s = 7.805)

PM adapts to the **actual mass distribution** rather than forcing a preferred shape.

### 3. This Pushes Investigation to Physical Properties

The real discriminators are likely:
- **Disk scale length** (R_d)
- **Total baryonic mass** (M_bar)
- **Surface brightness** (HSB vs LSB)
- **Morphology** (bars, warps, interactions)
- **Data quality** (number of points, error bars, inclination)

Not an internal geometric ratio.

---

## What Went Wrong with the Initial Pattern

### Selection Bias

The initial subset (~10-20 galaxies) happened to show clustering around 0.55. The full dataset (175 galaxies) reveals this was:
- **Accidental:** Random draw from a broad distribution
- **Not representative:** Subset was too small to capture true variance

### The Seductive Trap

0.55 is a "nice" number—close to golden ratio, suggests natural scaling. It's easy to see patterns in small samples.

**This is why you stress-test before publishing.**

---

## The Constrained Fit Experiment

### Test: Force r_c = 0.55 × r_s

Results across 175 galaxies:
- **85 improved** (48.6%)
- **21 unchanged** (12.0%)
- **69 degraded** (39.4%)

### Interpretation

This is **not** what you'd see if 0.55 were universal:
- Universal constraint → vast majority improve or unchanged
- This result → nearly even split between improve/degrade

**Conclusion:** The constraint helps some galaxies (acts as regularizer), hurts others (too restrictive), but is not a fundamental PM requirement.

---

## What This Teaches About PM Fitting

### The Good
1. **PM is geometrically flexible** - adapts to real mass distributions
2. **No hidden "magic ratio"** - no special cases to explain away
3. **Clean negative result** - rules out one class of constraints cleanly

### The Bad (For Simple Stories)
1. **No universal geometric signature** - can't use r_c/r_s as a PM "prediction"
2. **Fitting remains multi-parameter** - can't reduce to one master scale

### The Interesting
The constraint **does** help nearly half the galaxies. Why?

**Hypothesis:** Those galaxies have:
- Clean, symmetric disks
- Smooth mass profiles
- Good data quality

And the constraint acts as a **regularizer** against overfitting.

**The degraded galaxies likely have:**
- Bars, warps, or asymmetries
- Sparse or noisy data
- Non-circular motions

Where the constraint is **too restrictive** for the real dynamics.

---

## Revised Strategy: Physical Discriminators

### Next Steps

Instead of geometric ratios, analyze:

#### 1. Disk Scale Length vs χ²
Does PM prefer extended or compact disks?

#### 2. Baryonic Mass vs χ²
Is there a mass scale where PM struggles?

#### 3. Surface Brightness vs χ²
Do HSB and LSB galaxies fit differently?

#### 4. Morphology vs χ²
Do bars/warps correlate with poor fits?

#### 5. Data Quality vs χ²
Do sparse datasets or large error bars drive poor fits?

These will reveal **physical regimes** where PM works vs struggles, not arbitrary geometric patterns.

---

## Theoretical Implications

### What We Learned About PM Structure

1. **PM flow fields are not rigidly constrained by internal ratios.**
   - The medium allows diverse geometries.
   - This is consistent with a **continuum field theory** rather than a particle model.

2. **Good fits emerge from matching the baryonic distribution, not enforcing flow geometry.**
   - PM is doing what it should: responding to sources, not imposing templates.

3. **The r_c/r_s "sweet spot" was a statistical mirage.**
   - Small-sample patterns don't always scale up.
   - This is why comprehensive testing matters.

### What We Still Don't Know

1. **Why does the constraint help 85 galaxies?**
   - Is it acting as a regularizer?
   - Do those galaxies share morphological traits?

2. **Why does PM struggle with some morphologies?**
   - Bars/warps may require non-axisymmetric flow fields.
   - Current model assumes circular symmetry.

3. **Is there a weaker, morphology-dependent ratio?**
   - Maybe r_c/r_s ≈ 0.55 for clean disks, but varies for others.

---

## Honest Assessment

### This Is Good Science

You proposed a pattern, tested it rigorously, and **it didn't hold.** That's not failure—that's how you avoid fooling yourself.

### PM Is More Robust Than Expected

The lack of a required ratio means:
- PM doesn't need fine-tuning to work across diverse galaxies
- The theory is flexible enough to handle real astrophysical complexity
- You haven't accidentally built in a hidden assumption that would break later

### The Real Work Begins Now

With the "magic ratio" ruled out, the next step is:
1. **Correlate χ² with physical properties** (mass, size, morphology)
2. **Identify systematic trends** (not numerical artifacts)
3. **Derive theoretical scaling relations** from PM field equations
4. **Test predictions** against new datasets

---

## Appendix: Statistical Summary

### Excellent Fits (χ²_red < 1.0)
- Count: 30 galaxies
- r_c/r_s: 0.007 to 7.805 (median 0.511)
- No clustering around single value

### Poor Fits (χ²_red > 5.0)  
- Count: 39 galaxies
- r_c/r_s: 0.021 to 8.915 (median 0.532)
- Distribution nearly identical to excellent fits

### Constrained Fit Performance
- Improved: 85 galaxies (48.6%)
- Unchanged: 21 galaxies (12.0%)
- Degraded: 69 galaxies (39.4%)

**Interpretation:** Constraint acts as regularizer for some, too restrictive for others. Not a universal PM requirement.

---

## Bottom Line

**PM does not have a universal core-skeleton ratio.** This is a **strength**, not a weakness—it means the theory adapts to real galaxy diversity without forcing a template geometry.

The search now shifts to **physical discriminators:** mass scales, morphology, angular momentum distribution, and data quality.

This is exactly the kind of honest, iterative refinement that makes a theory robust.
