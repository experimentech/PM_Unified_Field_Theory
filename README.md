# Pushing-Medium: Unified Field Theory

## Overview

This repository documents the extension of the Pushing-Medium (PM) gravitational theory to a unified field framework encompassing gravity, electromagnetism, and plasma effects.

## Core Discovery

**PM predicts cross-force interactions that General Relativity cannot realize**, including:

1. **Unified lensing** from mass, charge, and plasma
2. **Current-induced dragging** from both mass and charge flows  
3. **Dynamic energy transfer** between gravitational and EM sectors

**However**, observational constraints show these interactions are naturally suppressed (charge couplings ≳10¹⁷× weaker than mass couplings), explaining why GR succeeds while PM remains structurally richer.

## Quick Start

### View the Analysis
```bash
python3 em_coupling_resolved.py
```

### Key Results
- Mass currents: **A_M ≈ 7.4×10⁻²⁸ m³/(kg·s)** [validated]
- Charge currents: **A_q < 3×10⁻¹⁰ m³/(C·s)** [constrained by lab null results]
- Magnetic potential: **κ_B ≈ 0** [ruled out]
- Hierarchy: **A_q/A_M < 4×10¹⁷**

## Documentation

### For Theory Understanding
- `FINDINGS.md` - Executive summary (start here)
- `CROSS_FORCE_INTERACTIONS.md` - Detailed interaction analysis
- `docs/README.md` - Complete document index

### For Technical Details
- `docs/cross-force-interactions-beyond-gr.tex` - Main theoretical document
- `docs/coupling-constants-reference.tex` - All coupling constants and constraints
- `docs/pm-gr-cross-force-comparison.tex` - Side-by-side PM vs GR comparison
- `docs/electromagnetic-coupling-constraints.tex` - Derivation of EM bounds

### For Implementation
- `src/pushing_medium/core.py` - Extended unified field implementation
- `src/pushing_medium/__init__.py` - Package interface

## What Makes PM Different From GR

| Aspect | General Relativity | Pushing-Medium |
|--------|-------------------|----------------|
| Gravity | Spacetime curvature (geometric) | Refractive index in medium (optical) |
| EM | Gauge field on curved spacetime | Excitation of same medium |
| Coupling | Only via stress-energy T^μν | Direct via sources S_φ, S_u |
| Cross-terms | None (forces independent) | Unified index and flow fields |
| Structure | Geometric vs gauge fields | Both are medium excitations |

## Key Insight

PM reveals that the **apparent separation of forces** might be an accident of coupling strengths, not fundamental physics. If A_q ~ A_M, we'd have noticed the unified structure long ago. Instead, nature chose A_q ≪ A_M, hiding the unity.

The "no special cases" principle leads to structurally unified equations whose coupling strengths are determined by nature, not theory.

## The Landmine and Resolution

**Problem**: Initial formulation included u_B = κ_B A (magnetic potential coupling), which predicted huge (~10 rad) deflections near lab magnets.

**Solution**: Drop κ_B A term; keep only source coupling S_u = A_M J_M + A_q J_q.

**Lesson**: PM medium couples to **currents** (sources), not **potentials** (derived fields). This avoids double-counting and respects observations.

## Validated Parameters

From `legacy/tests/calibration.json`:
```json
{
  "k_Fizeau": 1.0,
  "k_TT": 1.0,
  "mu_coeff": 1.4861356300677034e-27
}
```

These were empirically derived from solar system tests, binary pulsars, and gravitational lensing observations.

## Status

- **Gravity**: Validated in weak/semistrong regime (passes all GR tests)
- **EM extension**: Structure defined, constraints derived, landmine resolved
- **Plasma**: Structure defined, awaiting precision tests
- **Strong fields**: Unexplored (where PM and GR must diverge)

## Future Tests

1. **Charged-particle lensing**: Constrain μ_E
2. **Plasma cluster lensing**: Frequency-dependent deflection (constrain α_plasma)
3. **Atom interferometry**: Near mega-ampere currents (refine A_q bound)
4. **Magnetar observations**: Extreme B-field environments
5. **Multi-messenger GW**: Test c_φ = c_u assumption

## Philosophy

> "I find it hard to believe that the universe would have special cases."

PM demonstrates this principle by treating all forces as excitations of a single medium with the same mathematical structure. The coupling strengths are empirical properties of the medium, not theoretical necessities.

## Contact / Citation

[Add your information here when ready for publication]

---

**Last Updated**: 2026-03-02  
**Version**: Unified field formulation v1 with EM constraints
