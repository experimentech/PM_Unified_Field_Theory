# Pushing Medium — Theoretical Extensions and Open Problems

This document collects ideas and predictions that go beyond the verified core of PM.
Items are grouped by their status:

- **Open problem** — a known gap in the current formulation with no agreed resolution
- **Theoretical extension** — an idea that follows naturally from PM's structure but has not been worked out quantitatively
- **Candidate prediction** — something PM should predict differently from GR, but the computation has not yet been done or compared to data

Nothing in this document should be cited as a derived result of the theory.
Cross-references to the formula sheet are given where the relevant verified machinery exists.

---

## 1. Open Problems

### 1.1 The strong-field self-consistency gap

The most significant known problem in PM. The formula sheet's Diagnostic section documents
it in detail; the summary is:

PM's three strong-field ingredients — the Poisson equation, the equation of state
`P = (c²/2)(ρ − ρ_nuc)`, and the deformation potential `U(φ)` — are mathematically
unrelated. In any real medium these must be derivable from one underlying constitutive
relation. They are not, and this mismatch is invisible everywhere the field is weak
(φ ≪ 1) but is exposed in the high-φ interior of compact stars:

- `U(φ)` is absent from the Lagrangian, so the stability cap φ < 1 is imposed by hand
  rather than emerging dynamically.
- The thermodynamic pressure derived from `U(φ)` is inconsistent with the stated EOS
  by a factor of ~4 at lowest order, and is non-monotonic above φ ≈ 0.586.
- Gradient field energy ½|∇φ|² does not appear as a source in the Poisson equation.

**Status:** Three repair strategies (Options a, b, c in the formula sheet) have been
numerically probed. Option (c) — adding `U(φ)` as an additional pressure term — gives
the best observational match. A fully self-consistent derivation from a single Lagrangian
has not been achieved.

**What is needed:** Derive the EOS and the stability condition from a single modified
Lagrangian including a potential `−V(φ)`. The most natural candidate is
`V(φ) = U(φ)/c_φ²` leading to a nonlinear Poisson equation
`∇²φ − U′(φ)/c_φ² = −(8πG/c²)ρ_nuc e^φ`.

---

### 1.2 The preferred-frame argument

The introduction states: "No internal experiment can identify a preferred rest frame
because all instruments and particles are made of the same medium."

This is the right answer but it is asserted, not derived. A compressible medium
ordinarily does have a preferred rest frame — the frame in which the medium is at rest.
The reason PM evades this needs to be demonstrated explicitly. The argument presumably
rests on the fact that all measuring devices, including rulers and clocks, are themselves
configurations of the medium and therefore scale with `n` in exactly the way that cancels
any detectable asymmetry. This needs to be worked out as a formal derivation, not a
verbal claim.

**Status:** Unaddressed. No derivation exists in the codebase or documents.

---

### 1.3 The nature of the nuclear density scale

The reference density `ρ_nuc ≈ 2.3–2.8 × 10¹⁷ kg/m³` appears in the theory as a
phenomenological parameter — the density at which the medium's restoring pressure goes
to zero and matter's surface is defined. Its value is set to match observed neutron-star
properties. But the theory offers no mechanism for why the medium has this particular
reference density rather than any other. This is PM's equivalent of the cosmological
constant problem: the number is known, but its origin is not explained.

**Status:** Open. No proposed mechanism.

---

### 1.4 The cosmological energy-phase density

The cosmological constant Λ corresponds, in the two-phase picture, to a uniform
non-clumping energy-phase density. The compact-star phase transition produces an
energy-phase density at `φ = 1`, which is roughly `ρ_nuc c²` ≈ 2 × 10³⁴ Pa. The
cosmological energy-phase density inferred from Λ is ~10⁻¹⁰ J/m³, some 44 orders of
magnitude smaller. PM does not resolve this — it re-expresses the cosmological constant
problem in medium language without solving it.

**Status:** Open. The two-phase cosmology reproduces ΛCDM expansion history but does not
explain the magnitude of the energy-phase density.

---

### 1.5 Instantaneous vs. retarded flow field

The current formulation uses an instantaneous Poisson equation for the flow field:

```
∇²u = −A_M J_M − A_q J_q
```

This means changes in mass current propagate instantaneously to arbitrarily large
distances — a known problem of instantaneous action in any field theory. Legacy work
explored a wave equation instead:

```
□u_g = κ_J J_TT
```

where `□ = ∂²_t/c² − ∇²` is the d'Alembertian, giving the flow field a finite
propagation speed. This version is causal and structurally closer to the EM vector
potential. The consequences are significant:

- Gravitomagnetic effects would lag behind their source (relevant for binary inspiral
  phasing).
- The flow field would carry energy and be capable of radiating independently of
  the scalar φ field.
- The GW sector would couple directly to `u` dynamics rather than being an
  independent addition.

This was developed in legacy code but appears to have been set aside — possibly
because the instantaneous version gives the correct weak-field limit and is
computationally simpler. Whether the retarded form changes any observational
predictions at current precision is unknown.

**Status:** Explored in legacy documents. Not carried forward into the current
formulation. Should be revisited when the GW sector (inspiral phasing, ringdown) is
developed quantitatively.

---

### 1.6 Well-posedness of the PM initial value problem

Legacy documents contain a formal attempt to cast PM as a symmetric-hyperbolic PDE
system:

```
A⁰ ∂_t U + Aᵏ ∂_k U = R(U)
```

with state vector `U = (n, nv, h_TT, p_TT)` and an explicit energy norm:

```
ε = ½ [ P′(n)/n · n² + |nv|²/n + |p|² + c²|q|² ]
```

A symmetric-hyperbolic system with a positive-definite energy norm admits a
well-posed initial value problem — solutions exist, are unique, and depend
continuously on initial data. If PM's equations have this structure, it means PM
is a properly-defined dynamical theory rather than a collection of static
phenomenological formulas.

This work was apparently not completed or was lost. It is directly relevant to
the mathematical consistency open problem flagged in §4.3.

**Status:** Partially developed in legacy documents. Not verified or completed.
Recovering and finishing this analysis would resolve the well-posedness question
and potentially also the superluminal-mode and hyperbolicity questions simultaneously.

---

## 2. Theoretical Extensions

### 2.1 Matter and energy as phases of the medium (E = mc²)

The two-phase picture (matter phase: localised structured φ-configuration; energy phase:
uniform, non-clumping) provides a natural account of mass-energy equivalence:

- **Rest-mass energy** is the deformation energy `U(φ)` stored in a localised
  φ-configuration. This is already in the Hamiltonian section of the formula sheet.
- **Radioactive decay** corresponds to a metastable φ-configuration relaxing to a
  lower-energy configuration; the difference in `U(φ)` is released as EM or kinetic
  excitations.
- **Nuclear binding energy** is the reduction in total `U(φ)` when two separate
  φ-lumps merge into a single joint configuration.
- **Particle annihilation** corresponds to the erasure of a localised matter-phase
  pattern; the stored energy converts to propagating field excitations (EM, possibly
  gravitational).

This picture is self-consistent at the qualitative level and follows directly from the
medium interpretation. It is not quantitatively developed: PM has no nuclear model, no
quark structure, and no derivation of particle masses from the φ-dynamics. The above
is an interpretive framework, not a set of predictions.

**What is needed to make this quantitative:** A model of what determines the size,
energy, and stability of a minimal φ-lump (the PM analog of a nucleon). This probably
requires solving the nonlinear Poisson equation (see §1.1) in a compact domain.

---

### 2.2 Refractive cosmology and the Hubble tension

Preliminary numerical work (`demo_cosmology_refinement.py`) has shown that a
density-dependent decay rate for the refractive index,
`H(z) = H₀(1+z)^β` with β ≈ 0.8, reproduces ΛCDM luminosity distances to within ~4%
over `z = 0.1–1.5`. This suggests:

- Cosmic redshift is produced by the slow relaxation of `n(t)` toward 1, not by
  spatial expansion.
- The Hubble tension arises because `H₀` is not a universal constant but a local
  decay rate dependent on the medium's density. Observers in an underdensity (the
  Local Void) measure a higher `H₀` than the cosmic average.
- Dark energy is not a new substance; it is the energy-phase density of the medium.

**Status:** Phenomenological agreement with ΛCDM at the ~4% level. Not a derived
prediction from first principles. The decay rate `H(z) ∝ (1+z)^0.8` is fitted, not
derived. The Hubble tension explanation is qualitatively consistent with PM's structure
but has not been worked out as a quantitative prediction. Structure growth (fσ₈) has
been separately matched but with a calibration factor β ≈ 0.8 that needs a first-
principles justification.

---

### 2.3 Particles as topological φ-configurations

The matter-as-medium-configuration idea (§2.1) raises the question of what stabilises
a particle against dispersal: why doesn't a φ-lump just spread out and dilute? In
condensed-matter analogies (solitons, skyrmions, vortex lines), localised configurations
are stabilised by topology — they cannot be continuously deformed into the uniform state.

PM may have an analogous stabilisation mechanism if the φ-field admits topological
defect solutions. This has not been investigated. If it exists, it would provide a
mechanism for: particle number conservation, the discreteness of particle species, and
the minimum mass of a stable PM configuration.

**Status:** Entirely unexplored. No mathematics, no numerics.

---

### 2.4 Gravitational wave echoes and surface modes

GR compact objects (black holes) produce ringdown gravitational waves governed by
quasi-normal modes set by the horizon. PM compact objects have no horizon — they have
a physical surface at `φ → 0` where the matter phase meets ambient medium. This surface
can, in principle, partially reflect gravitational perturbations, producing echoes after
the main merger signal and a surface-mode spectrum rather than horizon QNMs.

Whether PM produces detectable echoes depends on the surface reflectivity, which depends
on the (currently unresolved) strong-field constitutive relation (§1.1). If the surface
is perfectly transmissive, there are no echoes. If it has significant reflectivity, echoes
at a timescale set by the object's light-crossing time might be observable in
high-SNR LIGO/Virgo/KAGRA events.

**Status:** Qualitatively predicted by the no-horizon structure. Not quantitatively
computed. Current GW data does not show echoes at a significant level, which constrains
but does not rule out PM, since the reflectivity is a free parameter until §1.1 is
resolved.

---

## 3. Candidate Predictions (not yet computed)

These are cases where PM's structure implies a specific, falsifiable prediction that
differs from GR, but the computation has not been performed and compared to data.

| Prediction | Observable | Status |
|---|---|---|
| Shadow interior brightness non-zero (no horizon suppression) | EHT M87\*, Sgr A\* | Ray-tracing code exists; prediction not extracted |
| Surface redshift bounded by φ_crit (z_max ≈ e−1 ≈ 1.72) | Neutron-star spectral lines | Formula exists; no comparison to observations |
| lensing shear and rotation curves consistent from same φ-field | Galaxy surveys | Rotation-curve fitting done; lensing shear not computed |
| H₀ varies with local density (voids vs. superclusters) | H₀ tension surveys | Qualitative argument made; no quantitative model |
| No frequency-independent dispersion in GW propagation | LIGO + EM counterparts | Not computed |
| Compact-object mass gap at ~30 M☉ (Option c) | LIGO mass spectrum | Order-of-magnitude match; full calculation needed |
| Lagrange point positions differ from GR/Newtonian prediction | Solar system astrometry, space mission orbit design | Legacy code (`lagrange_1.py`, `lagrange_2.py`) implemented; no quantitative comparison to observed L-point positions done |

---

## 4. Medium Manipulation and "Gravitational Engineering"

This section collects ideas for whether and how the PM medium could be actively
manipulated — whether it is in principle possible to engineer the local compression
state of space.

### 4.1 Electric currents as gravitational flow sources

The PM unified source equation is:

```
∇²u = −A_M J_M − A_q J_q
```

where J_M is a mass current and J_q is a charge current. If A_q ≠ 0, then any
large electric current would generate a flow field u that drags light and matter
— the same mechanism as gravitational frame-dragging, but driven by charge rather
than mass. This is a qualitatively different prediction from GR, which has no
equivalent mechanism at a comparable magnitude.

**Current observational status (from `current_dragging_analysis.py` and
`cross-force-hierarchy-honest.tex`):**
- Laboratory null results from optical systems near large electromagnets constrain
  A_q < 3×10⁻¹⁰ m³/(C·s).
- The dimensional estimate from classical electrodynamics gives
  A_q ~ 2.5×10⁻³⁷ m³/(C·s) — 27 orders of magnitude below the lab bound.
- The magnetic-vector-potential coupling κ_B has been effectively ruled out
  (would produce ~0.3c flow from a 1 MA solenoid, which is not observed).
- The practical working assumption is A_q ≈ 0, but this is an observational
  upper bound, not a theoretical derivation.

**What would make this interesting:** If A_q is non-negligible at some scale,
then a sufficiently large and structured current distribution could create a
localised flow field. The same Fizeau-dragging mechanism that governs light near
rotating stars would apply. This would constitute a laboratory-scale gravitational
effect driven by electromagnetism. The required currents are far beyond current
technology at the dimensional-estimate value of A_q, but the upper bound is much
looser.

**Status:** Quantitatively bounded. Not theoretically derived. `current_dragging_analysis.py`
implements the full comparison between mass-current and charge-current dragging
for magnetars and laboratory solenoids.

---

### 4.2 EM-field modification of the refractive index

The unified index equation is:

```
ln n_tot = μ_G Φ_G + μ_E Φ_E − α_plasma n_e
```

If μ_E ≠ 0, a strong electric field modifies the local refractive index of the
medium. This would mean that a sufficiently powerful electromagnetic field could
locally alter n(r) — compressing or rarefying the medium — in a way that affects
the propagation of both light and matter. In principle this opens the possibility
of:

- **Gravitational lensing by EM fields:** A focused high-intensity EM beam could
  create a local n-gradient that deflects light or slow particles.
- **Medium "valves":** A structured EM field that locally depletes the medium
  (reduces n toward 1) would act as a gravitational optical element — a region
  where light bends less, or where the pressure gradient driving gravitational
  acceleration is reduced.
- **Refractive index engineering:** In the same way that metamaterials engineer
  optical properties using structure, one might in principle engineer the
  gravitational refractive index using structured EM fields — if μ_E is non-zero.

**Current status:** μ_E is dimensionally estimated at ~2.5×10⁻³⁷ m²/C and is
entirely unconstrained by observation. The effect would be negligible for any
presently achievable field strength. This is speculative at a level far beyond
the current-dragging predictions, because it requires not just A_q ≠ 0 but also
a PM-specific coupling not present in standard electrodynamics.

**Status:** No quantitative computation. Entirely unexplored beyond dimensional
estimates.

---

### 4.3 Mathematical consistency: well-posedness and stability

Not yet formally checked:

- **Hyperbolicity:** Are the PM PDEs hyperbolic? If not, they may not admit
  well-posed initial value problems.
- **Superluminal modes:** Does PM permit any characteristic signal speed > c?
  If so, it violates causality.
- **Rotating configurations:** Can PM support stable rotating compact objects?
  Many medium-based theories fail here because rotation introduces centrifugal
  terms that destabilise the constitutive relation.
- **Energy bookkeeping in collapse:** PM predicts no singularity — but what
  happens to the energy as a compact object is pushed toward φ = 1? Does it
  accumulate in a radiation-dominated core stably, or is there a runaway?

**Status:** All four are open. They represent the sort of internal consistency
check that must be completed before PM can be taken as a complete dynamical
theory rather than a static field model.

---

## 5. PM and the Quantum Realm

PM is a classical field theory. However, several structural features of the medium
picture are superficially compatible with quantum mechanical concepts, and the
analogies are worth recording — not as claims, but as pointers to where a future
quantised extension might connect.

### 5.1 The medium and the quantum vacuum

In quantum field theory, the vacuum is not empty — it has a ground-state energy
density from zero-point fluctuations of each field. In PM, free space is not empty
either — it is the medium at its reference state n = 1, ρ = ρ_nuc. These are
structurally analogous: both are "the lowest-energy state of something physical."

The Casimir effect — an attractive force between closely spaced conducting plates,
attributed in QFT to vacuum fluctuations — would, in a PM picture, correspond to
a local modification of the medium's ground state between the plates, creating a
pressure imbalance. Whether PM can derive the Casimir force quantitatively has not
been investigated.

**Status:** Analogy only. No derivation attempted.

---

### 5.2 Wave-particle duality and the medium

In PM, a particle is a localised φ-configuration (§2.1 of this document) and light
is a propagating excitation of the medium. Both are excitations of the same
substrate. This is structurally reminiscent of de Broglie's pilot-wave picture,
where particles are guided by a wave in a physical medium.

The de Broglie wavelength λ = h/p could, in PM language, correspond to the
characteristic spatial scale of a particle's φ-configuration. If the medium has
a natural length scale (set by ρ_nuc and some stiffness parameter), there is a
candidate mechanism for why particles have the wavelengths they do. This is entirely
speculative and undeveloped.

**Status:** Structural analogy. No mathematics.

---

### 5.3 Uncertainty and the medium's fluctuation spectrum

The Heisenberg uncertainty principle ΔxΔp ≥ ℏ/2 could, in a medium picture, be
a statement about the minimum resolvable structure of a φ-configuration: you cannot
localise a particle into a region smaller than the medium's coherence length without
exciting it into a higher-energy state (spreading its momentum). If the medium has
a finite-temperature or finite-resolution granularity, this would naturally produce
uncertainty-like bounds. Whether this can be made quantitative and consistent with
the observed value of ℏ is unknown.

**Status:** Speculative. Mentioned here because the structural fit is suggestive,
not because a derivation exists.

---

### 5.4 Quantisation of PM fields

A Lagrangian density for PM exists (formula sheet §PM Lagrangian). In principle,
canonical quantisation can be applied to any Lagrangian field theory: promote
fields to operators, impose commutation relations. The result would be a quantum
field theory of the PM medium — with creation and annihilation operators for
φ-quanta (gravitational phonons) and u-quanta. Whether the resulting theory is
renormalisable, and whether its particle spectrum matches the known particles of
the Standard Model, is entirely open.

**Status:** The classical Lagrangian exists. Quantisation has not been attempted.

---

## 6. Things PM Does Not (Currently) Address

For completeness — areas where PM has no model at all, and which would need entirely
new development:

- **Quantum mechanics:** PM is a classical field theory. No quantisation scheme exists.
- **The Standard Model:** PM provides no account of electroweak or strong nuclear forces.
- **Particle masses:** Why a proton has the mass it does is not addressed.
- **Spin:** The φ-field is scalar; angular momentum of particles is not accounted for.
- **Charge:** The EM sector (ψ, **w**) describes EM waves but does not derive charge
  quantisation or the electron's charge-to-mass ratio.
