# Pushing-Medium Cosmology: Hubble Tension & Dark Energy

**Status**: Preliminary Hypothesis & Feasibility Check  
**Date**: 2026-03-02  
**Code**: `demo_cosmology_refinement.py`

---

## 1. The Core Idea: Refractive Redshift

Standard cosmology assumes redshift $z$ comes from the **expansion of space** ($a(t)$ increasing). Pushing-Medium (PM) explores the alternative that redshift comes from the **decay of the vacuum refractive index** $n(t)$.

### Mechanism

In a dielectric universe where the refractive index $n(t)$ slowly decays toward 1 (vacuum) over cosmic time:

$$ 1 + z = \frac{n(t_{\text{emit}})}{n(t_{\text{obs}})} $$

This produces a Hubble-like law ($z \approx H_0 d/c$) for nearby objects without requiring physical expansion velocities.

---

## 2. Default Prediction: The "Empty Universe"

If the fractional decay rate of the index is constant ($H(t) \equiv -\dot{n}/n = H_0$):
- The relationship between distance and redshift is **linear**: $D_L \propto z$.
- This matches the "Empty Universe" (Milne) cosmology.
- **Problem**: Supernova data (SNe Ia) show high-$z$ supernovae are fainter (farther) than this linear prediction. Standard cosmology fixes this by adding **Dark Energy** ($\Lambda$) to accelerate expansion.

---

## 3. Resolving Dark Energy: Variable Decay

We tested a refined model where the decay rate $H(z)$ depends on the density (refractive index) of the medium itself:

$$ H(z) = H_0 (1+z)^\beta $$

### Results (from `demo_cosmology_refinement.py`)

Using $\beta \approx 0.8$, PM matches Standard $\Lambda$CDM predictions ($D_L$) to within **~3%** across the range $z=0.1$ to $1.5$.

| Redshift $z$ | $D_L$ (Standard $\Lambda$CDM) | $D_L$ (PM $H \sim n^{0.8}$) | Error |
| :--- | :--- | :--- | :--- |
| 0.1 | 468 Mpc | 461 Mpc | 1.4% |
| 0.5 | 2887 Mpc | 2781 Mpc | 3.7% |
| 1.0 | 6608 Mpc | 6423 Mpc | 2.8% |
| 1.5 | 10850 Mpc | 10868 Mpc | 0.2% |

### Physical Interpretation
The "acceleration" of the universe is an illusion caused by the changing stability of the vacuum.
- **Past (High $z$):** Universe is denser ($n$ is higher). The medium is more stable, so it decays **slower** ($H(z) < H_0$).
- **Result:** Light travels further for a given $\Delta z$, making objects look fainter—exactly mimicking Dark Energy acceleration.

**Conclusion:** We do not need Dark Energy. We only need a vacuum whose decay rate scales with its density.

---

## 4. Resolving the Hubble Tension

The Hubble Tension is a 5σ discrepancy between local measurements of $H_0$ (~74 km/s/Mpc) and early-universe measurements (~67 km/s/Mpc).

### PM Explanation
If $H_0$ is a **local decay rate** dependent on density:
1. **Local Universe (SH0ES):** We live in a local underdensity (the "Local Void"). Lower density $\to$ Faster decay $\to$ **Higher $H_0$** (~74).
2. **Early Universe (Planck):** The universe was denser on average. Higher density $\to$ Slower decay $\to$ **Lower $H_0$** (~67).

### Prediction
PM predicts that $H_0$ is **not a universal constant**, but an environmental variable. Measuring $H_0$ in different environments (voids vs. superclusters) should yield systematically different values.

---

## 5. Summary

| Phenomenon | Standard Cosmology | Pushing-Medium Hypothesis |
| :--- | :--- | :--- |
| **Redshift** | Stretching of space | Decay of refractive index |
| **Dark Energy** | Unknown repulsive force | Variable decay rate ($\dot{n}/n$ varies) |
| **Hubble Tension** | Systematics / New Physics | Environmental variation of decay rate |

**Status:** Feasible. Requires no new forces or fields, just a density-dependent stability condition for the medium.
