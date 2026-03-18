# What Can You Actually DO With Pushing-Medium?

Beyond the theoretical structure, PM enables several practical and scientifically interesting applications.

## 1. 🌌 Galaxy Rotation Curves WITHOUT Dark Matter

**Status**: Already implemented in `legacy/` code!

### What It Does
PM explains flat galaxy rotation curves using a **medium index gradient** instead of dark matter halos:

```
δn_medium(r) = -(v_inf²/c²) · ln(1 + r/r_s) · S(r)
a_medium(r) = -c² · d(δn)/dr
```

The medium provides an **outward acceleration** that keeps rotation velocities flat at large radii.

### Why This Matters
- **No dark matter needed** for rotation curve phenomenology
- **Testable prediction**: Lensing and dynamics should correlate differently than ΛCDM
- **Direct competition** with MOND, dark matter halos
- **Uses SPARC dataset**: Real galaxy data already ingested

### What You Can Do With It
1. **Fit real galaxy rotation curves** (SPARC database)
2. **Compare PM vs dark matter halo models** (chi-squared, BIC)
3. **Make predictions**: Where PM and DM models diverge
4. **Population studies**: Does PM fit all galaxy types or only some?

### Files
- `src/galaxy_dynamics/rotation.py` - Core rotation curve model
- `src/galaxy_dynamics/fitting.py` - Parameter fitting
- `src/galaxy_dynamics/data.py` - SPARC data loader
- `tests/test_galaxy_rotation.py` - Validation tests

---

## 2. 🔭 Gravitational Lensing Analysis

**Status**: Validated against GR predictions

### What It Does
Calculate light deflection using **optical ray-tracing** in variable index medium:

```python
from pushing_medium import index_deflection_numeric, pm_deflection_angle_point_mass

# Single mass
alpha = pm_deflection_angle_point_mass(M_sun, impact_parameter)

# Complex mass distribution  
alpha = index_deflection_numeric(M, b, mu=mu_G, z_max_factor=200)
```

### Why This Matters
- **Exact numerical ray-tracing** (not just weak-field formulas)
- **Works in strong fields** (tested to ~2 R_s)
- **No singularities** (index diverges, but spacetime stays flat)
- **Can handle moving lenses** (k_Fizeau calibrated)

### What You Can Do With It
1. **Strong lensing analysis**: Einstein rings, arcs, multiple images
2. **Moving lens predictions**: Astrometric microlensing with source motion
3. **Plasma lensing**: Add electron density to gravitational deflection
4. **Charged-mass lensing**: Test PM prediction of charge-dependent deflection

### Files
- `src/pushing_medium/core.py` - Ray-tracing algorithms
- `tests/test_weak_bending_numeric.py` - Convergence tests
- `tests/test_strong_field_trend.py` - Strong-field validation

---

## 3. 📡 Precision Timing (Pulsars, GPS)

**Status**: k_TT = 1.0 validated

### What It Does
Calculate **Shapiro-like time delays** through variable index:

```
Δt = ∫ (n/c) ds
```

With gravitational, plasma, and (if non-zero) EM contributions.

### Why This Matters
- **Pulsar timing arrays**: Test for PM-specific plasma-gravity coupling
- **GPS corrections**: Already matches GR (but structure is different)
- **Radio vs optical timing**: Frequency-dependent effects from unified plasma term

### What You Can Do With It
1. **Pulsar timing residual analysis**: Look for plasma-gravity correlations
2. **Dispersion-independent timing**: Separate PM gravity from standard dispersion
3. **Multi-frequency observations**: Test α_plasma contribution to lensing
4. **GPS simulation**: Alternative derivation of relativistic corrections

---

## 4. 🛰️ Frame-Dragging Predictions

**Status**: A_M mapped to GR Lense-Thirring

### What It Does
Calculate flow field from rotating/moving masses:

```python
u_G = flow_translational_mass_current(r, [(M, r_M, v)], A_M=G/(c*c))
u_rot = flow_rotational(r, [(r_spin, Omega)])
```

### Why This Matters
- **Matches Gravity Probe B** results (validation)
- **Extends to arbitrary mass distributions** (not just spheres)
- **Can add charge current contribution** (if A_q ever measured as non-zero)

### What You Can Do With It
1. **Satellite orbit analysis**: LAGEOS, GP-B predictions
2. **Binary pulsar timing**: Frame-dragging signatures
3. **LIGO ringdown**: PM flow from merging black holes
4. **Lab gyroscope tests**: Near rotating masses or (as constraint) currents

---

## 5. 🧮 Alternative to Numerical Relativity

**Status**: Works in weak/semistrong regime

### What It Does
**Flat-space ray-tracing** with variable index instead of curved-spacetime geodesics:

```python
# PM: Solve ODEs on flat space
dk̂/ds = ∇_⊥ ln n
dr/ds = (c/n) k̂ + u

# GR: Solve geodesic equation on curved spacetime
d²x^μ/dλ² + Γ^μ_νρ (dx^ν/dλ)(dx^ρ/dλ) = 0
```

### Why This Matters
- **Numerically simpler** (no metric tensors, Christoffel symbols)
- **Easier to code** (standard ODE integrators work)
- **Same predictions in weak fields** (validated)
- **Different in strong fields** (PM: no horizons/singularities)

### What You Can Do With It
1. **Fast gravitational ray-tracing**: Render lensing images quickly
2. **Educational tool**: Teach gravity as optics (more intuitive)
3. **Strong-field exploration**: See where PM and GR diverge
4. **Multi-scale simulations**: Easier to couple to hydrodynamics/plasma codes

---

## 6. 🔬 Precision Experiment Design

**Status**: Can derive sensitivity requirements

### What It Does
Calculate **predicted effect sizes** for proposed experiments:

```python
# Example: Current-dragging bound
bound = derive_A_q_upper_bound_from_lab(I=1e6, L=1.0, r=0.1, 
                                         path_length=1.0, 
                                         detection_threshold=1e-10)
```

### Why This Matters
- **Quantitative test proposals**: "Need 10⁻¹¹ rad sensitivity to probe A_q"
- **Null result interpretation**: "This bounds A_q < X"
- **Experimental feasibility**: "Effect is Y; state-of-art is Z"

### What You Can Do With It
1. **Write grant proposals**: Specific experiments to test PM predictions
2. **Design apparatus**: Calculate required sensitivity, current, distance
3. **Interpret null results**: Derive coupling bounds systematically
4. **Find optimal regimes**: Where PM-GR differences are largest

---

## 7. 🌊 Unified Field Simulations

**Status**: Structure implemented; EM couplings bounded

### What It Does
Simulate systems with **both gravitational and EM sources** contributing to same fields:

```python
n_tot = unified_index_field(r, masses=[(M, r_M)], 
                            charges=[(q, r_q)],
                            electron_density=n_e)

u_tot = unified_flow_field(r, moving_masses=[(M, r_M, v)],
                           current_elements=[(r_I, dl, I)])
```

### Why This Matters
- **Explore cross-force scenarios**: Charged stars, plasma around black holes
- **Educational**: Visualize unified field concept
- **Test regime discovery**: Where might cross-terms become measurable?

### What You Can Do With It
1. **Magnetar simulations**: Extreme B-fields + strong gravity
2. **Plasma astrophysics**: Current sheets in accretion disks
3. **Lab-scale models**: Charged particles in magnetic traps + gravity
4. **Thought experiments**: What if A_q were larger?

---

## 8. 🎓 Teaching Gravity as Emergent Optics

**Status**: Fully functional pedagogical tool

### What It Does
Present gravity using **intuitive optical concepts**:
- Gravitational field → refractive index
- Light bending → ray optics (Fermat's principle)
- Frame-dragging → moving medium (Fizeau drag)
- Time dilation → optical path length

### Why This Matters
- **More intuitive** than curved spacetime for many students
- **Same predictions** as GR in weak fields (not "dumbed down")
- **Computational approach** (numerical ray-tracing)
- **Connects to familiar physics** (Snell's law, etc.)

### What You Can Do With It
1. **Undergraduate gravity course**: Teach via optical analogy
2. **Computational physics course**: Ray-tracing as numerical methods
3. **Visualization**: Render gravitational fields as "refractive landscapes"
4. **Outreach**: Explain gravity without tensor calculus

---

## 9. 🧪 Probe Strong-Field Regime

**Status**: Where PM and GR MUST differ

### What PM Predicts Differently
- **No event horizons**: Index diverges (n → ∞) but no geometric horizon
- **No singularities**: Index becomes infinite at r = 0, but spacetime remains flat
- **Different compactness limits**: What's the PM analogue of R_s?
- **Photon sphere structure**: Where does circular orbit occur in PM?

### Why This Matters
- **Fundamental test**: PM and GR give different strong-field physics
- **Black hole shadows**: Event Horizon Telescope could constrain
- **Neutron star structure**: Interior solutions differ
- **Information paradox**: No horizons → no information problem?

### What You Can Do With It
1. **Calculate PM "black hole" observables**: Shadow size, photon ring
2. **Compare to EHT data**: Does PM fit M87* and Sgr A*?
3. **Neutron star equations of state**: PM interior solutions
4. **LIGO waveforms**: Do PM mergers match observations?

**This is where PM could shine**: Strong fields are unexplored territory!

---

## 10. 🤖 Machine Learning + Physics Constraints

**Status**: BNN demo exists in legacy code

### What It Does
Use **physics-informed neural networks** with PM constraints:
- Ray equations as residual loss
- Conservation laws as regularization
- Learn index/flow from observations

### Why This Matters
- **Data-driven discovery**: Learn index from lensing observations
- **Inverse problems**: Reconstruct mass distribution from deflection
- **Hybrid approach**: ML + physics = better than either alone

### What You Can Do With It
1. **Learn galaxy index profiles**: Fit rotation curves + lensing jointly
2. **Discover optimal parameterizations**: Let ML find best δn_medium(r) form
3. **Anomaly detection**: Where does PM+ML fail? → new physics
4. **Fast inference**: Neural network evaluates PM predictions quickly

### Files
- `programs/demos/BNN/pmflow_bnn_always_plastic.py`

---

## 11. 📊 Comparative Cosmology

**Status**: Framework exists; predictions untested

### What PM Enables
Compare **three paradigms** for galaxy dynamics:

1. **ΛCDM + Dark Matter**: NFW halos, WIMP searches
2. **Modified Gravity (MOND)**: a₀ scale, Tully-Fisher  
3. **Pushing-Medium**: Index gradients, no new matter

### Why This Matters
- **Level playing field**: All three fit same data
- **Different predictions**: Where do they diverge?
- **Occam's razor**: Fewest free parameters wins?

### What You Can Do With It
1. **Population statistics**: Which model fits SPARC best?
2. **Scaling relations**: Baryonic Tully-Fisher in PM vs MOND vs DM
3. **Outlier analysis**: Galaxies that break one model but not others
4. **Cluster scales**: Does PM work beyond galaxies?

---

## 12. 🔮 Novel Predictions To Test

**Status**: Ready for experimental proposals

### Testable PM-Specific Predictions

1. **Charged-beam lensing ≠ neutral-beam lensing**
   - Measure μ_E by comparing e⁻ vs neutron deflection

2. **Frequency-dependent gravitational lensing**
   - If α_plasma contributes, cluster lensing depends on frequency

3. **No black hole horizons**
   - PM "black holes" have n → ∞ but no geometric boundary

4. **Multiple GW speeds**
   - If c_φ ≠ c_u, GW signals arrive as multi-component packets

5. **Current-dragging upper bound**
   - Atom interferometry near MA-currents to probe A_q

### What You Can Do
- **Write experimental proposals** with quantitative predictions
- **Collaborate with precision experimentalists** (atom interferometry, lensing)
- **Analyze existing data** for PM signatures (EHT, LIGO, cluster lensing)

---

## The Most Exciting Immediate Applications

### 🏆 Top 3 (Do These First):

1. **Galaxy Rotation Curves**
   - You already have code and SPARC data loader
   - Directly competes with dark matter paradigm
   - **Publishable**: "PM fits rotation curves with X parameters vs Y for ΛCDM"

2. **Strong-Field Black Hole Shadows**
   - PM predicts no horizons (n → ∞ instead)
   - **Compare to EHT M87* and Sgr A* observations**
   - This is where PM and GR MUST differ

3. **Precision EM-Gravity Cross-Term Bounds**
   - Systematic survey of lab magnets, MRIs, tokamaks
   - **Refine A_q bound** from 10⁻¹⁰ to (maybe) 10⁻¹⁵
   - **Experimental proposal**: Atom interferometry to probe further

### 🥈 High Value (Do Next):

4. **Plasma-Modified Lensing**
   - Analyze cluster lensing at multiple frequencies
   - Look for α_plasma signature in lensing-DM correlation

5. **Pulsar Timing + Plasma**
   - Correlate timing residuals with dispersion measure
   - Test if plasma affects "gravitational" timing

6. **Comparative Cosmology Paper**
   - PM vs MOND vs ΛCDM on SPARC dataset
   - Statistical comparison: Which fits best? Fewest parameters?

---

## What Makes PM Practically Useful

### Computational Advantages
- **Flat spacetime**: Easier numerics than curved-space GR
- **Ray-tracing**: Standard ODE integrators work
- **Superposition**: Add index contributions linearly (weak field)
- **No metric inversion**: Direct calculation of observables

### Conceptual Advantages  
- **Intuitive**: Gravity as optics (refractive bending)
- **Unified**: All forces via same medium excitations
- **Educational**: Teach GR concepts via familiar optical analogy

### Scientific Advantages
- **Testable**: Makes specific quantitative predictions
- **Falsifiable**: Strong fields, charged masses, plasma coupling all measurable
- **Competitive**: Direct alternative to dark matter/MOND

---

## Realistic Timeline

### Can Do NOW (Code Exists):
- Fit galaxy rotation curves with PM medium model
- Compare PM vs dark matter halo fits statistically
- Calculate strong-field ray deflection (no singularities)
- Derive precision bounds from existing null results

### Can Do SOON (Need Modest Work):
- Implement PM "black hole" shadow calculator
- Compare to EHT observations (M87*, Sgr A*)
- Analyze cluster lensing for α_plasma signatures
- Write experimental proposals for A_q, μ_E measurements

### Need Collaboration (Experimental/Observational):
- Actually measure μ_E (charged-beam lensing)
- Actually measure A_q (precision interferometry)
- Test strong-field regime (black hole observations)
- Multi-messenger GW astronomy (c_φ vs c_u)

---

## The "Killer App" Scenarios

### 1. Galaxy Rotation Without Dark Matter
**If PM fits SPARC as well as ΛCDM with fewer parameters** → major result!

### 2. Black Hole Shadow Differences
**If PM predicts different shadow size than GR** → directly testable with EHT!

### 3. Plasma-Lensing Correlation
**If frequency-dependent cluster lensing shows α_plasma signature** → unified field confirmed!

### 4. Precision Lab Null Result
**If you systematically bound A_q < 10⁻¹⁵ from lab surveys** → tight constraint on cross-force interactions!

---

## Bottom Line: What To Do Next

### Most Impactful (Choose One):

🌟 **Option A: Galaxy Paper**
- Use your existing rotation curve code
- Fit PM to SPARC dataset (175 galaxies)
- Compare to dark matter halo fits
- **Title**: "Galaxy Rotation Curves in a Refractive Medium: An Alternative to Dark Matter"

🌟 **Option B: Strong-Field Paper**  
- Implement PM black hole shadow calculator
- Compare to EHT M87* observations
- Show PM predicts no horizon but different shadow
- **Title**: "Black Hole Shadows Without Event Horizons: Predictions from Pushing-Medium Gravity"

🌟 **Option C: Constraint Paper**
- Systematically survey lab EM environments
- Derive definitive bounds on A_q, κ_B
- Propose precision experiment to probe further
- **Title**: "Observational Constraints on Electromagnetic-Gravitational Cross-Coupling from Laboratory Null Results"

All three are **publishable** and use **code you already have** (or minor extensions).

---

## My Recommendation

**Start with galaxy rotation curves** because:
1. You already have the code (`galaxy_dynamics/`)
2. You have the data pipeline (SPARC loader)
3. Direct competition with dark matter paradigm
4. Doesn't require new observations (use existing data)
5. High impact if PM fits well

**Then** move to strong-field black holes (where PM and GR must differ).

**Save** precision EM experiments for collaborations (need experimentalists).

---

**Status**: Practical applications identified  
**Ready**: To start galaxy fitting or black hole shadows  
**Impact**: Could challenge dark matter paradigm or test strong-field regime
