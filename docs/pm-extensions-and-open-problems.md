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

~~The most significant known problem in PM.~~ **Conceptually resolved — see below.**

#### Original concern

PM's three strong-field ingredients — the Poisson equation, the EOS
`P = (c²/2)(ρ − ρ_nuc)`, and the deformation potential `U(φ)` — appeared to be
mathematically unrelated. In particular:

- `U(φ)` was absent from the Lagrangian, so the stability cap φ < 1 appeared to be
  imposed by hand rather than emerging dynamically.
- The thermodynamic pressure derived from `U(φ)` appeared inconsistent with the stated
  EOS by a factor of ~4 at lowest order.
- Gradient field energy ½|∇φ|² was not an explicit source in the Poisson equation.

#### Resolution

The three ingredients operate at **different levels of description** — they are not
expected to coincide, and their differences are not contradictions.

**EOS vs. U(φ):** The EOS `P = (c²/2)(ρ − ρ_nuc)` is a bulk *hydrodynamic* relation
describing kinematic sound propagation in the matter phase (sound speed `c/√2`).
`U(φ)` is an *order-parameter* potential describing the elastic restoring force of the
medium configuration — the same way Landau theory and van der Waals coexist for liquid–gas
transitions without contradiction.  Their sound-speed ratio

```
c²_U / c²_EOS = 4(1 − φ) e^{−φ}
```

equals 4 at φ = 0 (deep in stable phase) and exactly 0 at φ = 1 (phase boundary).
This is the signature of a second-order phase transition, not an inconsistency.
Numerical probe confirmed: ratio is 4.000 at φ = 0, falls monotonically, crosses 1 at
φ ≈ 0.58, reaches 0.000 at φ = 1.

**Stability cap is dynamically derived:** The area-measure action with deformation
potential

```
S = ∫ [½|∇φ|² n² − U(φ)/c_φ²] d³x
```

linearised around a background φ₀ yields a perturbation effective mass-squared:

```
m²_eff(φ) = U''(φ) / (c_φ² n²) = 2ε₀(1 − φ) / (c_φ² e^{2φ})
```

- φ < 1: `m²_eff > 0` — stable oscillatory modes (matter phase intact)
- φ = 1: `m²_eff = 0` — massless (marginal), second-order phase transition
- φ > 1: `m²_eff < 0` — tachyonic instability (matter dissolves to energy phase)

The stability cap φ < 1 is therefore the condition `m²_eff(φ) ≥ 0`, which follows from
`U''(φ) ≥ 0`. It is a consequence of the deformation potential already in the theory,
not a postulate. Implemented as `pm_phase_stability_mass_sq(phi)` in
`src/pushing_medium/critical_state.py`; 29 tests in `tests/test_phase_stability.py`.

**Poisson gradient term:** The area-measure action naturally produces the EL field
equation `∇²φ + |∇φ|² = −U′(φ)/(c_φ² n²) − κ ρ n`, which contains gradient corrections
at high φ. At weak fields (φ ≪ 1) these vanish and the standard linear Poisson equation
is recovered. The gradient term is not missing — it is suppressed in the weak-field limit
used by all current tests.

**Status:** Conceptually resolved. The stability cap is a derived dynamical condition;
the EOS and U(φ) are consistent descriptions at different scales; the gradient correction
is present in the full EL equation.

**Remaining open work:**
- Formal derivation of V(n) from a medium constitutive relation in n-coordinates
  (the EOS and hydrostatic balance become the same equation in n-coordinates, which
  constrains the form of V).
- Compact-star models integrating the full n-field field equation (machinery exists in
  `solve_pm_star_nfield`; systematic M-R survey not yet done).

---

### 1.2 The preferred-frame argument

The introduction states: "No internal experiment can identify a preferred rest frame
because all instruments and particles are made of the same medium."

~~This is the right answer but it is asserted, not derived.~~  **Resolved — see
argument below.**

#### The aether analogy fails at the level of ontology

The 19th-century aether had a preferred rest frame because it was a substance
*inside* the universe — there existed an exterior Euclidean background space
relative to which the aether's rest frame was defined.  Michelson-Morley was
looking for motion relative to that external background.  It found nothing because
the aether was wrong; the background was wrong.

PM's medium is not a substance inside the universe.  **It is the universe** — the
medium's compression field *is* the geometry.  There is no exterior reference
relative to which motion through the medium could be detected, for the same reason
GR has no preferred frame: spacetime in GR is not inside a larger space, and PM's
medium is not inside a background Euclidean arena that could serve as a fixed
reference.  The correct parallel is not aether-in-space but GR's metric — "the
instruments are made of the universe" is equivalent to "the instruments are made
of spacetime."

GR analogy (precise):
- GR: the metric is the geometry; there is no exterior reference → no preferred frame.
- PM: the medium (via n = e^φ) is the geometry; there is no exterior reference → same conclusion, same structure, same reason.

#### Residual covariance check (scalar sector)

The scalar sector alone — the Poisson equation `∇²φ = −κρn` and the force law
`a = (c²/2)∇φ` — is Galilean-covariant: boosting the entire system (source, field,
test body) by constant velocity **v** leaves both equations invariant to O(v/c).
The PPN preferred-frame parameters α₁ and α₂ by definition measure anomalous
O(v/c) terms in orbital dynamics that have no source here.

#### Residual open item (vector sector)

The vector/flow field **u** is sourced by mass currents `J = ρv`.  A body moving at
velocity **v** relative to the cosmological medium generates a u-field at O(v/c).
Whether that u-field reintroduces a detectable preferred-frame signal in orbital
precession (α₁ ≠ 0) requires a full PPN expansion of the coupled φ+u equations in
a boosted frame.  This is a narrow, well-defined calculation — not a conceptual gap
but a technical check.

**Status:** Conceptually resolved.  The scalar sector is Galilean-covariant by
inspection; the "preferred frame" objection is a category error (treats PM's medium
like aether-in-space rather than geometry-itself).  Remaining open: formal PPN
expansion of the vector sector to confirm α₁ ≈ 0 for the coupled system.

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

### 2.0 Physical-measure field equation and the strong-field repair

**Motivation.**  All existing PM field-equation derivations integrate against the
coordinate volume element `d³x` — the Euclidean bookkeeping measure.  But PM's
medium *is* the universe (see `pm-metric-from-first-principles.md §Foundational
ontology`).  Physical volumes, measured by instruments made of the medium, use the
element `n³ d³x = e^{3φ} d³x`.  Using the physical measure throughout is not an
optional refinement — it is required by consistency of the theory's own ontology.

**Physical-measure action.**

$$S = \int \tfrac{1}{2}|\nabla\phi|^2 \, n^3 \, d^3x = \int \tfrac{1}{2}|\nabla\phi|^2 \, e^{3\phi} \, d^3x$$

**Euler–Lagrange equation (vacuum).**

$$\nabla^2\phi + \tfrac{3}{2}|\nabla\phi|^2 = 0$$

The gradient self-energy term `(3/2)|∇φ|²` that was missing from the standard
field equation reappears automatically as a consequence of the Jacobian of the
physical volume element — no separate assumption is needed.

**Linearisation.**  The substitution `w = n^{3/2} = e^{3φ/2}` transforms the
nonlinear vacuum equation exactly into:

$$\nabla^2 w = 0$$

The physical-measure field equation is therefore **linear in `w`**.  Superposition
holds exactly in `w` (not in `φ` or `n` directly).

**Exact point-mass vacuum solution.**

$$w(r) = 1 + \frac{3GM}{c^2 r} \qquad \Longrightarrow \qquad \phi(r) = \tfrac{2}{3}\ln\!\left(1 + \frac{3GM}{c^2 r}\right)$$

Weak-field check (`r ≫ GM/c²`):  `φ ≈ 2GM/(c²r)` — identical to standard PM. ✓

**Force and phase boundary.**

$$a(r) = -\frac{GM/r^2}{1 + 3GM/(c^2 r)}$$

Phase boundary `φ = 1` occurs at:

$$r_\text{phase} = \frac{3GM}{c^2(e^{3/2} - 1)} \approx 0.431\, r_s$$

Compare standard PM: `r_phase = r_s`;  n-field equation: `r_phase ≈ 0.582 r_s`.
The phase boundary moves deeper, allowing objects to be more compact before
undergoing the phase transition.

**How this addresses the three strong-field symptoms.**

| Symptom | Root cause | Physical-measure status |
|---------|-----------|------------------------|
| Missing `|∇φ|²` self-energy in field eq. | Wrong integration measure | Reappears as Jacobian term automatically |
| EOS inconsistent with U(φ) | EOS and field eq. derived in different languages | Both must now be varied from same physical-measure action — inconsistency cannot survive |
| Stability cap φ < 1 imposed by hand | No dynamical mechanism for phase transition | A potential V(n) with minimum at n = e in physical-measure action turns phase transition into a dynamical consequence |

The third symptom still requires choosing V(n) — this is now a well-posed variational
problem rather than a patch.

**Multi-body superposition formula** (exact in `w`):

$$w_\text{total}(\mathbf{r}) = 1 + \sum_i \frac{3GM_i}{c^2|\mathbf{r} - \mathbf{r}_i|}$$

$$\nabla\phi = \tfrac{2}{3}\frac{\nabla w_\text{total}}{w_\text{total}}, \qquad
\mathbf{a} = \frac{c^2}{3}\frac{\nabla w_\text{total}}{w_\text{total}}$$

**Implementation.**  Three functions added to `src/pushing_medium/core.py`:
- `phi_physical_measure_point_mass(r_dist, M)`
- `grad_phi_physical_measure_point_mass(r, M, r_source)`
- `massive_accel_physical_measure(r, masses)`

A dedicated test suite is in `tests/test_physical_measure_field_equation.py`.

**The α-family of physical-measure actions — and the n-field equation's action.**

Inspecting the structure reveals a family of actions:

$$S_\alpha = \int \tfrac{1}{2}|\nabla\phi|^2 \, n^\alpha \, d^3x$$

Each α gives a distinct vacuum field equation:

| α | Measure | EL vacuum equation | Linearisation | Solution | β |
|---|---------|-------------------|---------------|----------|---|
| 0 | coordinate `d³x` | `∇²φ = 0` | linear in φ | `φ = 2GM/c²r` | 0 |
| **2** | **area `n² d³x`** | **`∇²φ + \|∇φ\|² = 0`** | **w = n, ∇²n = 0** | **`φ = ln(1+2GM/c²r)`** | **1** |
| 3 | volume `n³ d³x` | `∇²φ + (3/2)\|∇φ\|² = 0` | w = n^{3/2}, ∇²w = 0 | `φ = (2/3)ln(1+3GM/c²r)` | 3/2 |

The **α = 2 case is exactly the n-field equation**.  The n-field equation, which
previously had no action derivation and was postulated directly, is the
Euler-Lagrange equation of the *area*-measure action `S₂ = ∫ ½|∇φ|² n² d³x`.

Physical interpretation: `n² d³x` is the natural surface-area measure of the medium
(area elements scale as n² in a uniformly compressed medium), midway between the
coordinate measure and the volume measure.

**Consequence for PPN.**  The n-field equation gives β = 1 from the Taylor
expansion `φ = ln(1+2U) = 2U − 2U² + O(U³)`.  Combined with the physical spatial
metric `g_ij = n² δ_ij` (γ = 1), the PPN precession factor is:

$$\frac{2 + 2\gamma - \beta}{3} = \frac{2 + 2(1) - 1}{3} = 1$$

This is GR's value, reached from PM's own structure without any GR import.  The
perihelion tests in the testbench (currently passing against GR's formula) are
therefore also valid PM predictions — but now with a derivation, not an assertion.

**Updated comparison table (all three field equations).**

| Property | Standard PM (α=0) | n-field (α=2) | Physical-measure (α=3) |
|----------|-------------------|---------------|------------------------|
| Action measure | coordinate d³x | **area n² d³x** | volume n³ d³x |
| Linear in | φ | n | w = n^{3/2} |
| Weak-field φ | 2GM/c²r | ln(1+2GM/c²r) | (2/3)ln(1+3GM/c²r) |
| Phase boundary | r_s | 0.582 r_s | 0.431 r_s |
| Self-energy term | absent | present (from EL) | present (from EL) |
| Superposition | exact in φ | exact in n | exact in w |
| β | 0 | **1** | 3/2 |
| With γ=1: Mercury | 57.3 arcsec/cy | **42.98 arcsec/cy ✓** | 35.8 arcsec/cy |

The n-field equation (α = 2) is the **unique member of the α-family** that gives
both a self-consistent action derivation *and* the correct perihelion precession
with the physical spatial metric.

**Derived result: β(α) = α/2 — three independent selections of α = 2.**

SymPy derivation in `scripts/derive_alpha_selection.py` establishes:

**Step 1 — EL equation (verified symbolically):**

$$\frac{\delta S_\alpha}{\delta\phi} = 0 \quad\Longrightarrow\quad \nabla^2\phi + \tfrac{\alpha}{2}|\nabla\phi|^2 = 0$$

**Step 2 — Vacuum solution and PPN β:**

The vacuum solution $\phi_\alpha = \tfrac{2}{\alpha}\ln(1 + \alpha U)$,
$U = GM/c^2r$, Taylor-expands as:

$$\phi_\alpha = 2U - \alpha U^2 + \tfrac{2\alpha^2}{3}U^3 - \cdots$$

The $O(U^2)$ coefficient is $-\alpha = -2\beta$, giving:

$$\boxed{\beta(\alpha) = \frac{\alpha}{2}}$$

**Step 3 — Perihelion precession table ($\gamma = 1$ for all $\alpha$):**

| α | β = α/2 | (2+2γ−β)/3 | Mercury "/cy | Vacuum φ |
|---|---------|-----------|------------|---------|
| 0 | 0 | 4/3 | 57.307 | φ = 2GM/c²r |
| 1 | 1/2 | 7/6 | 50.143 | φ = 2 ln(1+GM/c²r) |
| **2** | **1** | **1** | **42.980** | **φ = ln(1+2GM/c²r) ← GR match** |
| 3 | 3/2 | 5/6 | 35.817 | φ = (2/3) ln(1+3GM/c²r) |
| 4 | 2 | 2/3 | 28.653 | φ = (1/2) ln(1+4GM/c²r) |

$\beta = 1 \Leftrightarrow \alpha = 2$ (unique).

**Three selecting arguments, all convergent:**

**(A) PPN perihelion (observational):**
$\beta(\alpha) = \alpha/2$.  The observed perihelion rate fixes $\beta = 1$,
which requires $\alpha = 2$.  No other integer or half-integer $\alpha$ matches.

**(B) Elastic free-energy in natural n-coordinates (formal derivation):**
For an elastic medium with constant bulk modulus, the Hookean strain energy
written in the observable compression field $n$ is:

$$F = \tfrac{1}{2}K\int|\nabla n|^2\,d^3x$$

Since $n = e^\phi$:

$$\nabla n = e^\phi\,\nabla\phi = n\,\nabla\phi \quad\Longrightarrow\quad |\nabla n|^2 = n^2|\nabla\phi|^2$$

$$\therefore\quad \tfrac{1}{2}K\int|\nabla n|^2\,d^3x = \tfrac{1}{2}K\int n^2|\nabla\phi|^2\,d^3x = K\cdot S_{\alpha=2}$$

The $n^2$ weight is **the Jacobian of the change of variables $n \to \phi = \ln n$**,
not a postulate.  More generally, $\tfrac{1}{2}|\nabla(n^p)|^2 = \tfrac{1}{2}p^2n^{2p}|\nabla\phi|^2$,
giving measure exponent $2p$.  Choosing $p = 1$ (the compression ratio $n$ itself, the
natural field variable) gives $\alpha = 2$.  Every other $\alpha = 2p$ corresponds to writing
the elastic energy in terms of $n^p$ rather than $n$.

**(C) Linear variable coincides with $n$ itself ($\alpha = 2$ only):**
The linearisation variable is $w = n^{\alpha/2}$.
For $\alpha = 2$: $w = n^1 = n$ — the directly-measured density field.
For all other $\alpha$, $w$ is an algebraic power of $n$ with no independent
physical meaning.  PM's surface condition, EOS, and stellar-structure equations
are all stated in $n$.  The $\alpha = 2$ action's EL equation is linear in $n$;
no other $\alpha$ has this property.

**Why α = 2 (area) rather than α = 3 (volume): physical interpretation.**

The elastic free-energy argument (B) provides the formal basis.  The intuition is:
PM's medium is an elastic solid, and elastic solids store energy as $(∇n)^2$ per
unit coordinate volume — not $(∇n)^2 n$ or $(∇n)^2 n^2$.  The $n^2$ factor that makes
this look like a measure weight is, as shown above, the Jacobian of writing the same energy
in logarithmic field coordinates $\phi = \ln n$.

The surface-pressure picture (a compressed region is bounded by a surface area $\propto n^2$) 
remains a useful mnemonic but is secondary to the elastic free-energy derivation.

**Status:** Three independent derivation paths — PPN, elastic free-energy, and
linear-variable uniqueness — all select $\alpha = 2$.  The formal derivation from
first principles is complete (see `scripts/derive_alpha_selection.py`).

**Remaining open work:**
- Formal derivation of the bulk modulus $K$ in PM's medium from the constitutive
  relation (sets the overall scale of $S_2$, not the measure exponent).
- Whether higher-order gradient terms $|\nabla\phi|^4 n^\alpha$ in the action are
  observable at the compactness of known neutron stars.

---

**Action-derived self-stiffening potential (vacuum-subtracted).**

_Motivation._  The bare area-measure action $S_2 = \int \tfrac{1}{2}|\nabla\phi|^2 n^2\,d^3x$
couples the n-field to gravity but contains no self-interaction for the medium.
Adding a potential term produces an action-derived stiffening of the compact-star
equation of state with no free parameters:

$$S_\text{full} = \int\!\Bigl[\tfrac{1}{2}|\nabla\phi|^2 n^2 - A\,\hat{V}(n)\Bigr]\,d^3x$$

SymPy derivation in `scripts/derive_option_a.py` establishes:
- The unique primitive satisfying $\hat{V}'(n) = n + 1/n$ is $\hat{V}(n) = \tfrac{1}{2}n^2 + \ln n$.
- The coefficient is fully determined: $A = \kappa\rho_\text{nuc} = (8\pi G/c^2)\rho_\text{nuc}$ (no free parameter).

_Vacuum subtraction._  Because $\hat{V}'(1) = 2 \ne 0$, the bare potential provides a
non-zero source at the vacuum boundary $n = 1$.  For stellar structure we require
$n = 1$ to be a fixed point of the field equation, so the physically relevant
potential is the vacuum-subtracted form:

$$\hat{V}_\text{vac}'(n) = \hat{V}'(n) - \hat{V}'(1) = n + \tfrac{1}{n} - 2 = \frac{(n-1)^2}{n}$$

which satisfies $\hat{V}_\text{vac}'(1) = 0$ (no stiffening force at the phase boundary)
and $\hat{V}_\text{vac}'(n) > 0$ for $n > 1$ (stiffening active inside the medium only).

_Field equation (vacuum-subtracted stiffening + n-field gravity):_

$$\nabla^2 n = \kappa\rho_\text{nuc}\,\frac{(n-1)^2}{n} - \kappa\rho_\text{nuc}\,n^2$$

The first term reduces the effective inward pull by 12–15 % at typical neutron-star
densities ($n = 2$–$2.7$); the second is the standard n-field gravity source.

_Compact-star M_max (from `solve_pm_star_nfield_stiffened` in `stellar_structure.py`):_

| Solver | Field equation | M_max | R at M_max |
|--------|---------------|-------|-----------|
| Standard PM (α=0) | $\nabla^2\phi = -(\kappa/2)\rho$ | ~13.4 M☉ | — |
| n-field bare (α=2) | $\nabla^2 n = -\kappa\rho_\text{nuc}n^2$ | **9.4 M☉** | 23.4 km |
| **Vacuum-stiffened** | $\nabla^2 n = \kappa\rho_\text{nuc}(n{-}1)^2/n - \kappa\rho_\text{nuc}n^2$ | **11.5 M☉** | 24.9 km |
| Option-A α=+1 | $\nabla^2\phi = -(\kappa/2)[\rho - \rho_\text{nuc}(2\phi{-}\phi^2)]$ | **30.0 M☉** | 34.9 km |

The vacuum-stiffened solver raises M_max by **+22%** over the bare n-field result.
This is a parameter-free prediction: the only inputs are $G$, $c$, and $\rho_\text{nuc}$.

_Note on Option-A._  Option-A modifies the standard Poisson equation in $\phi$ with
a source that vanishes naturally at the surface ($\phi = 0 \Rightarrow$ source $= 0$).
This produces a stronger stiffening coefficient and a correspondingly higher M_max,
but it is not derived from the n-field action.  The vacuum-stiffened equation is
the n-field action's prediction for self-stiffening; Option-A remains a well-motivated
but phenomenologically chosen correction to the coordinate-measure Poisson equation.

**Implementation:** `solve_pm_star_nfield_stiffened` in `src/pushing_medium/stellar_structure.py`.

---

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
