# Pushing Medium — Documentation Index

This directory contains the formal documentation for the Pushing Medium (PM) framework.

---

## Directory Structure

```
docs/
├── README.md                       # This file
├── pm-introduction.md              # Accessible prose introduction to PM
├── pm_vs_gr_reference.md           # PM vs GR computed values (10 observables)
├── pushing_medium_paper.md         # Structured summary paper
├── critical_state_computation.md   # Matter/energy transition: φ_crit derivation + results
├── latex/                          # LaTeX source files
├── pdf/                            # Compiled PDFs (one per LaTeX source)
├── figures/                        # PNG figures
├── notes/                          # Working notes and analysis logs
└── LLM_notes_and_conversations/    # Raw LLM session transcripts (reference only)
```

---

## Markdown Documents

| File | Description |
|---|---|
| [pm-introduction.md](pm-introduction.md) | Conceptual introduction: medium, matter, gravity, light; what PM is and is not |
| [pushing_medium_paper.md](pushing_medium_paper.md) | Structured paper: field equations, calibrated constants, GR benchmarks, SPARC results |
| [critical_state_computation.md](critical_state_computation.md) | Matter/energy transition: Landau-Ginzburg derivation, φ_crit = 1, ρ_crit, P_crit, implementation, tests |
| [pm_vs_gr_reference.md](pm_vs_gr_reference.md) | Reference table of PM vs GR values for 10 observables, with reproduction code |

---

## LaTeX Documents (`latex/`)

Compiled PDFs are in `pdf/` with matching filenames.

### Gravity — Core Formulation

| File | Description |
|---|---|
| [gravity-with-effective-metric.tex](latex/gravity-with-effective-metric.tex) | Static PM gravity: index/flow fields, continuous distribution, effective isotropic metric, PPN parameters |
| [gravity-effective-metric-from-Fermats-principle.tex](latex/gravity-effective-metric-from-Fermats-principle.tex) | Effective spacetime metric derived from Fermat's principle in a moving medium |
| [v2-gravity.tex](latex/v2-gravity.tex) | PM Gravity v2: complete formulation + §7 numerical constraints + §8 critical compression |

### Unified Formulations

| File | Description |
|---|---|
| [v1-unified-formulations.tex](latex/v1-unified-formulations.tex) | PM Fluid Theory v1.5: Part I = gravity, Part II = unified EM/plasma/current sources (15 sections) |
| [hamiltonian-formulation-v1.tex](latex/hamiltonian-formulation-v1.tex) | Canonical Hamiltonian; §7 derives Landau-Ginzburg deformation energy and stability → φ_crit = 1 |

### Electromagnetic Sector

| File | Description |
|---|---|
| [em-field-equations-v1.tex](latex/em-field-equations-v1.tex) | PM EM field equations v1: dynamic equations, Maxwell-like emergence, Lagrangian, coupling to gravity |

### Stress-Energy Tensor

| File | Description |
|---|---|
| [stress--energy-tensorv1.1.tex](latex/stress--energy-tensorv1.1.tex) | PM stress-energy tensor v1.1: canonical T^μν, conjugate momenta, energy/momentum density, stress tensor |

### Cross-Force Interactions and Coupling Constants

| File | Description |
|---|---|
| [cross-force-hierarchy-honest.tex](latex/cross-force-hierarchy-honest.tex) | Honest analysis: 4 PM-specific mechanisms, observational hierarchy, A_q and κ_B bounds, κ_B landmine resolved |
| [coupling-constants-reference.tex](latex/coupling-constants-reference.tex) | Complete reference: all coupling constants, dimensional analysis, observational constraints, open questions |
| [pm-gr-cross-force-comparison.tex](latex/pm-gr-cross-force-comparison.tex) | Side-by-side PM vs GR comparison table; quantitative examples |

### Formula Sheet

| File | Description |
|---|---|
| [pm-formula-sheet.tex](latex/pm-formula-sheet.tex) | **Consolidated formula sheet**: all PM equations and constants in one two-column document (fields, Lagrangian, Hamiltonian, stress-energy tensor, EM sector, kinematics, metric, EOS, critical state, calibration constants) |

### Calibration

| File | Description |
|---|---|
| [calibration-constants-derivation.tex](latex/calibration-constants-derivation.tex) | Derivation of μ_G, k_TT, k_Fizeau: physical reasoning, matching procedure, connection to GR |

---

## Key Results

### Critical State (Matter/Energy Transition)
- Critical compression: **φ_crit = 1** (where U''(φ) = 0)
- Critical density: **ρ_crit ≈ 6.25 × 10¹⁷ kg/m³** (≈ 2.7 × nuclear saturation density)
- Critical pressure: **P_crit ≈ 1.78 × 10³⁴ Pa** (PM EOS: c²/2 × ρ_nuc(e−1))
- Theory: `latex/hamiltonian-formulation-v1.tex` §7
- Implementation: `src/pushing_medium/critical_state.py`
- Tests: `tests/test_critical_state.py` (11 tests, all passing)

### Gravitational Calibration
- θ_GR = 1.7516" for solar limb deflection (reproduced by PM)
- k_Fizeau = 1 (exact, by construction)
- μ_G connects PM pressure-gradient force to Newtonian G

### Cross-Force Observational Hierarchy
| Coupling | Value/Bound | Meaning |
|---|---|---|
| A_M | ~10⁻²⁸ m³/(kg·s) | Mass currents couple strongly |
| A_q | <10⁻¹⁰ m³/(C·s) | Charge currents ≳10¹⁷× weaker |
| κ_B | ≈ 0 | Magnetic potential coupling ruled out |

### SPARC Galaxy Fitting (175 galaxies)
- PM median χ² = 15.68, mean = 198
- See `notes/SPARC_RESULTS_HONEST.md` for honest assessment

---

## Working Notes (`notes/`)

| File | Description |
|---|---|
| [notes/SPARC_RESULTS_HONEST.md](notes/SPARC_RESULTS_HONEST.md) | Current honest 175-galaxy SPARC assessment |
| [notes/SPARC_RESULTS_ANALYSIS.md](notes/SPARC_RESULTS_ANALYSIS.md) | Full 175-galaxy analysis with dark-matter comparison |
| [notes/PM_RATIO_CONSTRAINT_FINDINGS.md](notes/PM_RATIO_CONSTRAINT_FINDINGS.md) | Final result: no universal r_c/r_s ratio |
| [notes/PRACTICAL_APPLICATIONS.md](notes/PRACTICAL_APPLICATIONS.md) | PM applications: galaxies, Bullet Cluster, strong-field, gravitational waves |
| [notes/SPARC_FITTING_GUIDE.md](notes/SPARC_FITTING_GUIDE.md) | Quick-start guide for running the SPARC pipeline |
| [notes/HYBRID_SKELETON_APPROACH.md](notes/HYBRID_SKELETON_APPROACH.md) | Hybrid vector/tensor skeleton method for galaxy fitting |
| [notes/SKELETON_TOPOLOGY_FINDINGS.md](notes/SKELETON_TOPOLOGY_FINDINGS.md) | Skeleton topology vs fit quality analysis |
| [notes/COSMOLOGY_FINDINGS.md](notes/COSMOLOGY_FINDINGS.md) | Preliminary hypothesis: refractive redshift (feasibility only, not validated) |

---

**Status**: Gravitational sector validated; EM sector formulated; critical state computed; cross-force constraints derived.
