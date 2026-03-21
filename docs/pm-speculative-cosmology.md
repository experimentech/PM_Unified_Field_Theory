# Pushing Medium — Speculative Cosmological Implications

**Status: speculative / exploratory**

This document captures lines of inquiry that emerge naturally from PM's structure
but have not been worked out quantitatively.  They are worth recording because they
are *internally motivated* speculation — each follows from machinery already present
in the theory rather than from a separate hypothesis bolted on.  None should be
cited as a derived result.

Cross-references to verified PM machinery are given where relevant.

---

## 1. The phase transition as a cosmological fixed point

### Background

PM's medium can exist in two phases: a matter phase (structured φ-configurations
below the critical density) and an energy phase (uniform medium at or above φ = 1).
The transition is not a mathematical singularity but a change of state of the medium
itself.  For a compact object this provides a natural collapse floor (see
`pm-extensions-and-open-problems.md §2.2`).

### The speculative step

Scale this up.  Apply the same phase transition physics to a gravitationally
collapsing universe as a whole — a system with enough mass-energy that compression
eventually drives φ → 1 *everywhere* simultaneously.

### What PM's structure implies

1. **Gravity self-extinguishes.**  Gravity in PM requires a density gradient, ∇φ ≠ 0.
   A uniform medium at φ = 1 everywhere has no gradient.  The mechanism that drives
   gravitational attraction ceases to operate — not because it is overcome, but
   because it has no differential to act on.

2. **The flow field u vanishes for the same reason.**  The flow field is sourced by
   mass currents J_M.  Matter in PM is a structured φ-configuration.  In the
   uniform energy phase there are no such structures; matter has dissolved back into
   the medium.  With no matter there are no mass currents, and u → 0 by the same
   logic as ∇φ → 0.

3. **The state is the initial condition of the universe.**  A uniform, critical-density
   medium with no gravity, no structure, and no flow field is precisely the PM
   description of the pre-structure epoch.  There is no singularity — no mathematical
   breakdown — at this point.  It is a phase state, not a boundary.

4. **The uniform state is unstable to fluctuations.**  Small perturbations — thermal,
   quantum, or topological — will break the uniformity.  Some regions dip below φ = 1.
   Gravity re-ignites locally.  Structure begins to form.  The universe *re-emerges*
   from the same medium.

### Consequence

The cosmology is cyclic, not by assumption but as a direct consequence of the phase
transition physics already required for compact objects.  The critical state is both
the endpoint of collapse and the launchpad of re-expansion.  No separate bounce
mechanism, no singularity, no information paradox — the medium simply changes phase
and then differentiates again.

**What is needed to make this rigorous:** A time-dependent field equation governing
how φ evolves globally under self-gravity, and an analysis of the stability of the
φ = 1 uniform state under perturbations.  This is currently inaccessible without
first resolving the strong-field self-consistency gap.

---

## 2. A cosmological Chandrasekhar limit

### Background

The Chandrasekhar limit for white dwarfs is not a measured parameter but a balance
condition: the exact point where electron degeneracy pressure and gravitational
self-compression reach equilibrium.  It is derived from the theory rather than fit
to data.

### The PM analogue

PM's phase transition defines a critical density of the medium, ρ_c, corresponding
to φ = 1.  The question is: what total mass-energy M_universe is required to
self-compress a spherically distributed uniform medium to ρ_c everywhere
simultaneously?

The answer, if it can be derived, would be a fixed point of the theory — not a free
parameter.  Crucially, because matter and energy are different phases of the *same*
field in PM, the total entering this calculation is genuinely the total mass-energy
content.  There are no separate dark matter, baryonic matter, and radiation terms
that must be summed by hand.  The phase transition density does not distinguish
matter from energy — only total medium-density matters.

### The observational hint

In ΛCDM, the fact that the observable universe's Schwarzschild radius 2GM/c² is
close to its actual radius R — i.e. 2GM/c²R ≈ 1 — is treated as a consequence of
the flatness condition and sometimes dismissed as a coincidence.  In PM this
approximate equality would instead be the condition: the universe's total content
is close to, but not over, the PM cosmological limit.  Not a coincidence — a
structural constraint the theory itself enforces.

**What is needed:** A complete variational principle in the strong field.  Without it
this derivation could produce a number but not the confidence that the number *had to*
come out that way.  This is therefore directly downstream of resolving the
strong-field crack (see `pm-extensions-and-open-problems.md §1`).

---

## 3. Hubble tension as a phase transition artefact

### Background

The Hubble tension is the ~5σ discrepancy between H₀ measured from CMB-calibrated
early-universe methods (~67 km/s/Mpc) and from late-universe local distance ladder
measurements (~73 km/s/Mpc).  ΛCDM has no mechanism that would produce a real
difference between these values — both should be measuring the same constant.

### The PM mechanism

PM has a physical transition between the matter phase and the energy phase at a
specific density threshold φ = 1.  The CMB measurements sample the universe when
it was deep in the energy-phase-dominated epoch.  Local distance ladder measurements
sample the late universe, which is matter-phase dominated at galactic scales.

If the medium's effective expansion rate changes character as the energy phase gives
way to matter phase dominance — because the equation of state of each phase is
different — then early and late measurements would be sampling the same medium in
two different states.  The ΛCDM model assumes a single smooth expansion history
parametrised by one constant H₀.  PM's two-phase structure could instead predict
that the measured value is epoch-dependent: not a measurement error, and not a
systematic, but a physical change in the medium's behaviour at the phase transition.

### Why this is particularly interesting

PM does not need to *add* this mechanism.  The two-phase structure is already there
for independent reasons (compact objects, galaxy dynamics, the force law).  The
question is only whether the resulting equation-of-state change at the phase
transition has the right magnitude to account for the observed ~8% offset.

**What is needed:** The full equation of state for the energy phase, including its
contribution to the expansion rate.  The PM cosmology scripts
(`scripts/pm_twophase_cosmology.py`) contain early work on this but the EOS
connection to Hubble tension has not been explicitly computed.

---

## 4. Dark energy and the energy-phase pressure

### Background

The observed accelerating expansion of the universe is explained in ΛCDM by a
cosmological constant Λ, representing a constant vacuum energy density.  It is
the least understood component of ΛCDM — no physical mechanism for Λ is established.

### The PM framing

In PM the energy phase carries energy density that gravitates.  This energy density
is dynamic — it is coupled to φ through the two-phase constitutive relation — and
it exerts pressure.  The question is whether this pressure has the right equation
of state to mimic Λ.

A cosmological constant requires w = P/ρc² = −1.  The energy-phase medium in PM
does not obviously satisfy this, but the equation of state has not been fully
derived from the field equations.  It is possible that the PM energy phase has
w ≈ −1 over the relevant density range, in which case dark energy is not a separate
ingredient but a consequence of the energy phase's thermodynamics.  It is equally
possible that PM needs a separate dark-energy term, in which case it merely renames
the problem.

This is an open calculation, not a closed no.  The basis for hope is that the medium
carries energy with pressure that is already in the theory for other reasons.

**What is needed:** Derive the full equation of state P(ρ) for the energy phase from
PM's two-phase constitutive relation and check whether w ≈ −1 over the range
relevant to late-time cosmological expansion.

---

## 5. The inter-force weak link

### Background

In PM the electromagnetic force and gravitational force are coupled through the
medium — both propagate through the same substrate, and the refractive index n
affects photon propagation.  The coupling is real but very small in the weak-field
regime.

### The speculative question

Does the coupling produce a small but non-zero mixing between the force sectors
that grows in the strong field?  At the phase transition, both gravity and
electromagnetism are operating in the same medium at its critical density.  It is
not clear that the two sectors remain cleanly separable in that regime.

A weak link between forces would be a genuine PM prediction with no GR counterpart,
because GR does not model a physical medium at all.  The size of the coupling, and
whether it is observable at any accessible energy scale, is entirely unknown.

**What is needed:** A careful analysis of the joint PM equations for φ and the EM
field in the strong-field regime, looking for cross-terms that do not vanish at
high φ.

---

## Summary

| Idea | PM basis | Status | Key missing piece |
|------|----------|--------|-------------------|
| Cyclic cosmology via phase transition | Phase transition already in theory for compact objects | Speculative, structurally motivated | Time-dependent global field equation; stability analysis of φ=1 state |
| Cosmological Chandrasekhar limit | Phase transition density + total mass-energy | Speculative, requires strong-field fix | Variational principle; resolution of strong-field crack |
| Hubble tension as phase artefact | Two-phase EOS change at threshold | Speculative, mechanism present | equation of state of energy phase vs expansion rate |
| Dark energy from energy-phase pressure | Energy phase carries energy density with pressure | Open calculation | Derive w = P/ρc² for energy phase |
| Inter-force weak coupling | Both forces in same medium | Highly speculative | Strong-field analysis of joint φ + EM equations |

All five lines of inquiry share the same upstream dependency: **resolving the
strong-field self-consistency gap** (see `pm-extensions-and-open-problems.md §1`).
Until a complete variational principle for the strong field exists, each of these
can be argued for qualitatively but not derived quantitatively.
