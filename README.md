# Pushing-Medium (PM): Testbench and Theory

## Overview

Pushing-Medium is a scalar-field gravity theory: gravity arises from a compressed
elastic medium with refractive index $n(r) = e^\phi$.  This repository is the
research testbench for PM — quantitative comparison of PM predictions against GR
and observations across 18+ falsification vectors.

**Test suite: 655/655 passing** (as of 2026-04-05)

## What PM is

Pushing-Medium is set in Euclidean space + absolute time.  Gravity is not spacetime
curvature — it is the compression of an elastic medium that pervades space.  The
denser the medium, the slower light travels and the stronger the gravitational
force.  This density variation under gravitational pressure (encoded by the scalar
field $\phi = \ln n$) produces all the observable effects attributed to spacetime
curvature in GR.

| | GR | Pushing-Medium |
|--|----|---------||
| Background | Curved spacetime (Riemann) | Euclidean space + absolute time |
| Gravity mechanism | Spacetime curvature | Medium density variation under gravitational pressure |
| Gravity source | Stress-energy $T^{\mu\nu}$ | Scalar field $\phi = \ln n$ |
| Force law | Geodesic deviation | $\mathbf{a} = +(c^2/2)\nabla\phi$ |
| Light bending | Null geodesics on curved metric | Refraction through denser medium near mass |
| EOS | Arbitrary (TOV) | Self-consistent: $\rho = \rho_\text{nuc}e^\phi$ |

PM reproduces all weak-field GR predictions exactly (PPN $\gamma = \beta = 1$).
Differences appear in the strong-field regime.

## Falsification status

See `docs/LLM_notes_and_conversations/` and `/memories/repo/falsification-status.md`
for the full ledger. Summary:

**Tests PM passes (reproduces GR / observations):**
Light deflection, Shapiro delay, frame dragging, gravitational redshift, Einstein
radius, perihelion precession (Mercury 42.98"/cy), GW propagation speed ($v_{GW}=c$),
GW inspiral waveform at 0PN, SPARC galaxy rotation curves, cosmology fσ8.

**Active falsification challenges:**
- **Mass gap**: GWTC events (GW150914: 30+36 M☉, GW190521: 66+85 M☉) exceed all
  action-derived PM M_max values (9–11.5 M☉). Option-A reaches 30 M☉ but is not
  action-derived.
- **EHT shadow**: M87* and Sgr A* masses are $10^5$–$10^8\times$ PM M_max. Even
  if objects could form, PM shadow ≥ 4.6% larger than GR BH (EHT matches GR to ~5%).
- **Tidal deformability**: GW170817 $\tilde{\Lambda} \leq 800$; PM 1.4\,M$_\odot$ baseline gives $\Lambda \approx 378$ (within
  90\% CI, but 1\textendash2$\sigma$ above the GR-preferred range 70\textendash300). With threshold tuning (f=1.18) $\Lambda \approx 315$.
- **NS radii (NICER)**: PM NS radii 30–50% larger than observed.

## Recent theoretical work

- **α = 2 uniquely selected** — proved β(α) = α/2; three independent derivations
  (PPN, elastic free-energy, linearisation variable). See `scripts/derive_alpha_selection.py`.
- **Option-A derived from action** — $S_\text{full}$ with $\hat{V}(n) = \tfrac{1}{2}n^2 + \ln n$,
  coefficient $A = \kappa\rho_\text{nuc}$ (no free parameter). See `scripts/derive_option_a.py`.
- **Vacuum-stiffened solver** — M_max = 11.5 M☉ (+22% over bare n-field). See
  `src/pushing_medium/stellar_structure.py::solve_pm_star_nfield_stiffened`.

## Quick start

```bash
# Run the test suite
.venv/bin/python -m pytest --tb=short -q

# Run a specific falsification comparison
.venv/bin/python scripts/perihelion_precession_comparison.py
.venv/bin/python scripts/gw_inspiral_comparison.py
.venv/bin/python scripts/eht_shadow_comparison.py

# Compute M-R curves (all solver variants)
.venv/bin/python scripts/pm_mass_radius_curve.py
```

### Key numerical results
- Perihelion (Mercury): **42.98"/cy** (obs: 42.98"/cy) ✓
- GW propagation: **$v_{GW} = c$** exactly ✓
- SPARC galaxies: fitted without dark matter (`results/sparc_fit_results.json`)
- NS M_max (n-field action, vacuum-stiffened): **11.5 M☉**
- NS M_max (Option-A α=+1): **30.0 M☉** (not action-derived)
- GW150914 components: **30.6 + 35.6 M☉** → exceed action-derived M_max

## Documentation

- `docs/pm-extensions-and-open-problems.md` — theoretical gaps, derivations, open problems
- `docs/pm-introduction.md` — PM theory introduction
- `docs/pm_vs_gr_reference.md` — side-by-side PM vs GR comparison
- `WHAT_TO_DO_NEXT.md` — current open problems and suggested next steps

## Source layout

```
src/pushing_medium/
    core.py                — field equations, force law, light bending
    stellar_structure.py   — compact star solvers (5 variants)
    critical_state.py      — RHO_NUC, RHO_CRIT, phase boundary
    cosmology.py           — PM two-phase cosmology
src/galaxy_dynamics/
    fitting.py, rotation.py, data.py, halos.py
src/general_relativity/
    classical.py           — GR reference functions for comparison
tests/                     — 603 tests across 25+ files
scripts/                   — comparison scripts for each falsification vector
```

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

| Domain | Status |
|--------|---------|
| Weak-field gravity | ✅ All GR tests pass |
| GW waveforms (0PN) | ✅ PM = GR exactly |
| Galaxy rotation curves | ✅ SPARC fits done |
| Cosmology fσ8 | ✅ Matches ΛCDM |
| Strong-field (mass gap) | ❌ PM M_max < GWTC masses |
| EHT shadow | ❌ M87*/Sgr A* masses >> PM M_max |
| Tidal deformability | ⚠️ Within 90% CI, mild tension |
| Rotating stars | ⬜ Not yet |
| Cosmological perturbations | ⬜ Not yet |

## Next steps

See `WHAT_TO_DO_NEXT.md` for current priorities.

---

**Last Updated**: 2026-04-05 — 655/655 tests passing
