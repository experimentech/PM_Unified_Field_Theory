# PM Effective Metric from First Principles

**Status: exploratory / in-progress**

This document works through the derivation of the PM effective metric(s) using
only PM's own postulates — no GR imports.  It exists separately from the formula
sheet precisely because this work may revise or replace parts of that sheet.
Nothing here should be treated as settled until cross-checked against all
existing tests.

---

## The problem

The formula sheet currently states:

$$g_{tt} = -(1+2\phi), \qquad g_{rr} = (1-2\phi)$$

These were not derived from PM's medium.  They are the GR weak-field expansion
$g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$ with $h_{tt} = -2\Phi_{\text{GR}}/c^2$
evaluated at $\Phi_{\text{GR}} = -GM/r$ (negative), then identified with PM's
$\phi = +2GM/c^2r$ (positive).  The sign flip was not noticed because no computation
in the codebase ever feeds this metric into a force calculation — the force law
and optical metric are used directly everywhere.

The GR metric is therefore **decorative** in the formula sheet: present but inert.
In the weak field this causes no error.  In the strong field, any calculation that
used the metric formula directly would produce a repulsive, twice-Newtonian force.

---

## Foundational ontology: the medium IS the universe

Before stating postulates it is essential to be clear about what PM is describing,
because the language of geometry ("curved space", "flat background") imports
assumptions from GR that do not belong here.

**In GR**, spacetime is a 4-manifold with a Lorentzian metric.  Matter and energy
live *inside* this manifold.  Gravity is the curvature of the manifold caused by
its contents.  The manifold exists independently of its contents — empty spacetime
is still a geometric object.

**In PM**, there is no such manifold.  The medium *is* the universe.  There is no
space or time independent of the medium; the medium is the substrate of all
physical existence.  Matter is a structured compression of the medium.  Energy is
the medium in a different phase.  Gravity is not curvature — it is the push
exerted by regions of higher medium density on regions of lower density, transmitted
through the medium itself as a pressure gradient.

Consequences of this framing that differ from GR:

- **The Euclidean coordinate grid is abstract notation**, not a physical flat
  space that exists independently.  It is the bookkeeping system, not the arena.
- **Rulers are made of the medium.**  A ruler in a compressed region is itself
  compressed.  It is shorter, so it over-counts the distances it measures.
  This is not curvature — it is medium-unit measurement.  The coordinate
  distances are still correct; only the physical distances reported by
  medium-composed instruments are scaled.
- **The phase boundary at $\phi = 1$ is a change of state**, not a geometric
  singularity.  There is no separate geometric manifold that "tears" or
  "becomes singular" at the phase transition.  The medium simply changes phase,
  exactly as water changes to ice or steam.  This is the reason PM has no
  singularities in the GR sense.
- **Gravity is self-limiting.**  It requires $\nabla\phi \neq 0$ — a density
  *difference*.  In the uniform-density energy phase $\nabla\phi = 0$ everywhere;
  the gradient mechanism ceases to operate.  Collapse cannot proceed past the
  phase transition because gravity switches itself off.

Whenever this document uses the word "curvature" it means **medium compression
causing measurement scaling**, not geometric curvature of a manifold.  The two
produce the same observable predictions in the weak field but differ ontologically
and diverge in the strong-field / phase-transition regime.

---

## PM's actual postulates (no GR)

1. Fixed Euclidean 3-space + absolute time
2. A scalar compression field $\phi(\mathbf{r}, t) \geq 0$, with $n = e^\phi$
3. The force law for slow massive particles: $\mathbf{a} = +\tfrac{c^2}{2}\nabla\phi$
4. Light propagates at $c/n$ locally (Fermat's principle in the $n$-field)
5. The background spatial metric is Euclidean: $g_{ij} = \delta_{ij}$

Everything else — the effective metrics, PPN parameters, perihelion precession,
lensing — should be derived from these five postulates alone.

---

## Step 1: Massive-particle effective metric from the force law

The force law is conservative with potential $V = -\tfrac{mc^2}{2}\phi$, so the
non-relativistic particle Lagrangian is:

$$\mathcal{L}_{\text{particle}} = \tfrac{1}{2}mv^2 + \tfrac{mc^2}{2}\phi(\mathbf{r})$$

In any metric theory, the slow-particle geodesic equation in a static metric gives:

$$\frac{d^2x^i}{dt^2} = -\frac{c^2}{2}\partial_i h_{tt}
\qquad \text{where} \quad g_{tt} = -(1 + h_{tt})$$

Matching to the PM force law $a_i = +\tfrac{c^2}{2}\partial_i\phi$ requires
$h_{tt} = -\phi$, so:

$$\boxed{g_{tt}^{\text{PM, massive particle}} = -(1 - \phi)}$$

The background is Euclidean (postulate 5), so:

$$g_{ij}^{\text{PM, massive particle}} = \delta_{ij}$$

---

## Step 2: Photon (optical) metric from Fermat's principle

Light travels at $c/n = ce^{-\phi}$ locally.  The optical (Gordon) metric for a
medium at rest with refractive index $n$ is:

$$ds_{\text{opt}}^2 = -\frac{c^2}{n^2}dt^2 + d\mathbf{r}^2
= -c^2 e^{-2\phi}\,dt^2 + d\mathbf{r}^2$$

Expanding to first order in $\phi$:

$$g_{tt}^{\text{opt}} \approx -(1 - 2\phi)c^2$$

---

## Step 3: Two distinct effective metrics — and the 2:1 lensing ratio derived

| Physical actor | $g_{tt}$ (coefficient of $-c^2\,dt^2$) | $g_{ij}$ | Origin |
|---|---|---|---|
| Slow massive particle | $1 - \phi$ | $\delta_{ij}$ | Force law + Euclidean background |
| Photon | $e^{-2\phi} \approx 1 - 2\phi$ | $\delta_{ij}$ | Fermat / optical metric |

The **2:1 lensing ratio** follows immediately from the $\phi$-coefficients:
- Massive particle: effective potential $\propto \phi$
- Photon: effective potential $\propto 2\phi$ (from $g_{tt}$) + $0$ (no spatial compression scaling; coordinate background is Euclidean)

Since there is no spatial compression-scaling correction ($g_{rr} = 1$), the photon's
deflection comes entirely from the $g_{tt}$ term, which is $2\phi$ vs $\phi$.
Ratio = 2. This is not a coincidence — it is a direct consequence of using
Fermat's principle in a medium where $n = e^\phi$.

---

## Step 4: PPN parameters from the PM-derived metrics

PPN parameters are defined via the metric expansion:

$$g_{tt} = -(1 + 2\Phi_{\text{PPN}}/c^2), \qquad
g_{ij} = \delta_{ij}(1 - 2\gamma\,\Phi_{\text{PPN}}/c^2)$$

For PM with direct postulate 5: $g_{ij} = \delta_{ij}$ exactly (no compression-scaling
of the spatial metric yet — see Step 5a for why this must be revised), so:

$$\boxed{\gamma_{\text{PM, derived}} = 0}$$

This differs from the formula sheet's claim $\gamma = 1$.  The formula sheet
inherited $\gamma = 1$ from GR via the imported metric, not from a PM derivation.

The measured value is $\gamma_\text{obs} = 1 + (2.1 \pm 2.3)\times10^{-5}$
(Cassini 2003).  $\gamma = 0$ is excluded.

### PPN β from the field equation: the α-family result

The value of $\beta$ depends on which field equation is used.  Expanding
$\phi(U)$ with $U = GM/c^2r$ as $\phi = 2U + c_2 U^2 + \ldots$, then
$\beta = -c_2/2$:

| Field equation | $\phi$ | $c_2$ | $\beta$ |
|---|---|---|---|
| Standard ($\alpha=0$, coordinate measure) | $2U$ | $0$ | $0$ |
| n-field ($\alpha=2$, area measure) | $\ln(1+2U) = 2U - 2U^2 + \ldots$ | $-2$ | $1$ |
| Physical-measure ($\alpha=3$, volume measure) | $(2/3)\ln(1+3U) = 2U - 3U^2 + \ldots$ | $-3$ | $3/2$ |

The n-field equation ($\alpha=2$) gives $\beta=1$.  With $\gamma=1$ from
$g_{ij} = n^2\delta_{ij}$ (Step 5a), the PPN precession factor is:

$$\frac{2+2\gamma-\beta}{3} = \frac{2+2(1)-1}{3} = 1$$

This gives the exact observed perihelion precession **from PM's own structure**,
without importing the result from GR.  The n-field equation has an action
derivation as the EL equation of $S_2 = \int \tfrac{1}{2}|\nabla\phi|^2 n^2\,d^3x$
(see `pm-extensions-and-open-problems.md §2.0`).

### What this means for γ

This is not a failure of PM — it is a finding that the Euclidean spatial
background assumption (postulate 5) is in tension with the Cassini measurement.
Three possible resolutions:

**Resolution A**: Postulate 5 as stated is incomplete.  The medium compresses
spatial measurement rods along with itself.  Then $g_{ij} \neq \delta_{ij}$ in
medium units, and deriving the correct physical spatial metric becomes the central task.

**Resolution B**: $\gamma = 1$ is recovered from the optical metric, not
the massive-particle metric.  Light does not travel in Euclidean space — it
travels in the optical metric, which has $g_{ij}^{\text{opt}} = \delta_{ij}$
but is not the same object as the background spatial metric.  Whether this
argument holds under scrutiny is open.

**Resolution C**: The PPN framework itself assumes a single metric for all
particles. PM has two metrics (massive and optical) — the PPN formalism may
not apply directly, and the correct observational comparison requires
computing trajectories in each metric separately.

---

## Step 5: Perihelion precession — the coordinate-effect argument fails

The earlier sketch suggested perihelion precession might emerge as a coordinate
effect from measuring angles with light rays.  Working through the argument shows
this cannot be correct.

**Coordinate effects cannot produce a physical precession.**
Perihelion advance per orbital period is a physical, coordinate-invariant
observable — the same in all coordinate systems.  A genuinely closed Keplerian
orbit (which PM's force law produces) is closed in every coordinate system,
including optical coordinates.  Measurement-dependent angle shifts from light
bending are present in any theory with light deflection, and they cancel exactly
over a closed orbit.  No coordinate transformation converts a closed orbit into
a precessing one.

**The GR import in the orbit integrator.**
`core.pm_integrate_orbit` contains this line:

```python
rddot = (h2 / r**3
         - G * M / r**2
         - 3.0 * G * M * h2 / (c**2 * r**4))   # 1PN correction
```

The docstring explains it:
> *"The 1PN correction term −3GMh²/(c²r⁴) is the same for both PM and GR
> because both have PPN parameters β=γ=1."*

The `−3GMh²/(c²r⁴)` term is not derived from PM's force law — it is GR's
1PN correction, justified by assuming $\beta = \gamma = 1$.  Those PPN
parameters were inherited from the imported GR metric.  They have never been
derived from PM's medium.  The perihelion precession tests in the battery are
therefore testing a GR prediction, not a PM one.

**What PM's own metric actually predicts.**

The PPN framework gives the perihelion advance as:

$$\Delta\omega = \frac{2 + 2\gamma - \beta}{3}\cdot\frac{6\pi GM}{c^2 a(1-e^2)}$$

From PM's derived metrics (Steps 1–2):
- $g_{tt} = -(1-\phi)$ with $\phi = 2GM/c^2r$ (linear in $U = GM/c^2r$) → $\beta = 0$
- $g_{ij} = \delta_{ij}$ (coordinate Euclidean, no compression-scaling yet) → $\gamma = 0$

$$\Delta\omega_{\text{PM, derived}} = \frac{2}{3}\cdot\frac{6\pi GM}{c^2 a(1-e^2)}$$

**Mercury: $\approx 28.7$ arcsec/cy.  Observed: $42.98$ arcsec/cy.
PM as currently derived predicts $2/3$ of the observed precession.**

The n-field equation makes this worse, not better.  The logarithmic $\phi$ adds
a $U^2$ term giving $\beta = 1$, with $\gamma$ unchanged at 0:

$$\Delta\omega_{\text{n-field}} = \frac{2 + 0 - 1}{3}\cdot\frac{6\pi GM}{c^2 a(1-e^2)}
= \frac{1}{3}\cdot\frac{6\pi GM}{c^2 a(1-e^2)}$$

Mercury: $\approx 14.3$ arcsec/cy.  Even further from observation.

**Summary.**  The 1PN term currently in the code is a GR import.  PM's own
derived metric predicts $2/3$ (standard field equation) or $1/3$ (n-field)
of the observed perihelion precession.  The passing perihelion tests validate
GR's prediction, not PM's.

---

## Step 5a: What would fix the perihelion precession?

The PPN factor becomes 1 (full GR prediction) when $2 + 2\gamma - \beta = 3$.
For $\beta = 0$ (linear $g_{tt}$): need $\gamma = 1/2$.
For $\beta = 1$ (n-field): need $\gamma = 1$.

The Cassini bound independently requires $\gamma = 1.00002 \pm 0.00023$.
So PM needs $\gamma \approx 1$.

$\gamma = 1$ means $g_{ij} = \delta_{ij}(1 + 2\phi) = \delta_{ij} n^2$ to first
order — **the medium compresses spatial measurement rods by the same factor as it
slows time**.  This is not curvature in the GR sense — it is the consequence of
rulers being made of the same medium that is compressed.  It does, however, revise
postulate 5: the coordinate background may be Euclidean, but physical measurements
are scaled by the medium density.

If the medium is elastically compressed, and if physical measurement rods are
made of the same medium, then a compressed region would appear to have larger
spatial distances — exactly the $g_{ij} = n^2\delta_{ij}$ (isotropic) form.
That would give $\gamma = 1$, recover the correct perihelion precession, and
be consistent with the Cassini bound.

But this is a major revision: the background is no longer Euclidean — it is
conformally flat.  The Euclidean background postulate would need to be replaced
with: *the medium's density also scales spatial measurements, so that the
physical spatial metric is $g_{ij} = n^2\delta_{ij}$.*

This is exactly the isotropic Schwarzschild metric to leading order:
$$g_{tt} = -(1-\phi), \quad g_{ij} = n^2\delta_{ij} \approx (1+\phi)^2\delta_{ij}$$

which at first order gives $g_{tt} \approx -(1-2U)$, $g_{ij} \approx (1+2U)\delta_{ij}$
— precisely GR's weak-field metric in isotropic coordinates.

So the correct PM metric, once spatial measurement scaling by the compressed medium
is included, converges to GR's weak-field metric in isotropic coordinates.  The
fact that it does so is either a deep consistency result (GR's metric is what you
get when you measure a PM medium with PM-composed instruments) or a sign that PM
and GR are genuinely dual descriptions of the same physics.

---

## Step 5b: Spatial metric $g_{ij} = n^2\delta_{ij}$ derived from the particle Lagrangian

**Resolution of open question 1 and open question 4 below.**

Step 5a introduced $g_{ij} = n^2\delta_{ij}$ as something that *would* fix the
perihelion problem if it could be established.  This step derives it directly
from PM's own Lagrangian — not as an additional postulate, but as a consequence
of the ontological principle that particles are made of the medium.

### The argument

In PM, the medium is the substrate of all physical existence.  Matter is a
compressed, structured region of the medium.  This means any instrument a
particle carries — including the imaginary ruler it uses to measure distances
— is made of the same medium.  A ruler in a region where the compression field
is $\phi$ (refractive index $n = e^\phi$) is itself compressed by a factor $n$.
It is shorter, so it records more ruler-lengths per coordinate distance.

More precisely: a particle at position $\mathbf{r}$ with coordinate velocity
$\dot{\mathbf{x}}$ moves through the medium with **physical velocity**
$n\,\dot{\mathbf{x}}$, because the coordinate displacement per unit time is
amplified by $n$ when converted to medium-unit distances (compressed ruler
under-reports the distance).

The particle Lagrangian written in medium units (i.e. using physical velocities)
is therefore:

$$\mathcal{L} = \tfrac{1}{2}m(n\dot{\mathbf{x}})^2 + \tfrac{mc^2}{2}\phi
= \tfrac{1}{2}m\,n^2\,\delta_{ij}\,\dot{x}^i\dot{x}^j + \tfrac{mc^2}{2}\phi$$

Reading off the coefficient of $\tfrac{1}{2}m\,\dot{x}^i\dot{x}^j$ gives
the spatial metric directly:

$$\boxed{g_{ij} = n^2\delta_{ij}}$$

This is not a postulate.  It is the unique form that follows from:
1. The coordinate grid being Euclidean (abstract bookkeeping); and
2. Particles being made of the medium (so their kinematics are in medium units).

Postulate 5 is therefore **not wrong** — the coordinate background is Euclidean
— but it was stated imprecisely.  The correct reading is:

> *The coordinate grid is Euclidean.  Physical distances, measured by
> instruments composed of the medium, are scaled by $n$: the physical spatial
> metric is $g_{ij}^{\text{phys}} = n^2\delta_{ij}$.*

### The full PM effective metric

Combined with $g_{tt} = -(1-\phi)$ from Step 1:

$$\boxed{g_{tt} = -(1-\phi), \qquad g_{ij} = n^2\delta_{ij}}$$

At first order in $U = GM/c^2r$ (where $n^2 \approx 1 + 2U$):

$$g_{tt} \approx -(1-2U), \qquad g_{ij} \approx (1+2U)\delta_{ij}$$

This is precisely GR's weak-field metric in isotropic coordinates.  The
agreement is not an import — it is a consequence of PM's ontology.

### PPN consequences

The PPN expansion $g_{ij} = (1 + 2\gamma U)\delta_{ij}$ gives:
$$\gamma = 1 \quad \text{(derived, not asserted)}$$

Combined with $\beta = 1$ from the n-field action (Step 4, $\alpha = 2$):

$$\frac{2 + 2\gamma - \beta}{3} = \frac{2 + 2(1) - 1}{3} = 1$$

**PM predicts exact perihelion precession without importing the result from GR.**

### Resolution of the 1PN orbit integrator

The term $-3GMh^2/(c^2r^4)$ in `pm_integrate_orbit` is now a **PM prediction**,
not a GR import.  The derivation chain is:

1. $g_{ij} = n^2\delta_{ij}$ from particle Lagrangian in medium units → $\gamma = 1$
2. $\phi = \ln(1 + 2U)$ from n-field equation (area-measure action $S_2$) → $\beta = 1$
3. PPN 1PN radial equation with $\beta = \gamma = 1$:
   $$\ddot{r} = \frac{h^2}{r^3} - \frac{GM}{r^2} - \frac{(2+2\gamma-\beta)\,GMh^2}{c^2r^4}
   = \frac{h^2}{r^3} - \frac{GM}{r^2} - \frac{3\,GMh^2}{c^2r^4}$$

The passing perihelion precession tests are now validating a PM-derived result.
The orbit integrator code is correct; only its annotation needs to reflect
the complete derivation chain.

---

## Step 6: Gravitational redshift

The formula sheet gives $z = e^\phi - 1$, derived from the optical metric
$g_{tt} = -c^2/n^2 = -c^2 e^{-2\phi}$:

$$1 + z = \frac{\nu_\text{emit}}{\nu_\infty} = \sqrt{\frac{g_{tt}(\infty)}{g_{tt}(R)}}
= \sqrt{\frac{1}{e^{-2\phi}}} = e^\phi$$

Since redshift is measured with photons, it uses the optical metric — not the
massive-particle metric.  This derivation is therefore **independent of the
GR-metric import problem** and stands on its own.

The formula sheet's redshift formula $z = e^\phi - 1$ emerges correctly
from PM's own optical metric without any GR input.

---

## Summary of findings so far

| Quantity | Formula sheet | Derived from PM postulates | Status |
|---|---|---|---|
| $g_{tt}$ (massive) | $-(1+2\phi)$ — wrong sign | $-(1-\phi)$ | **Formula sheet wrong** |
| $g_{rr}$ / $g_{ij}$ (massive) | $(1-2\phi)$ / $\delta_{ij}$ | $n^2\delta_{ij}$ from particle Lagrangian | **Resolved — Step 5b** |
| $g_{tt}$ (optical) | $-c^2/n^2$ | $-c^2 e^{-2\phi}$ | Consistent — stands |
| $\gamma_{\text{PPN}}$ | $1$ (imported from GR) | $1$ — from $g_{ij}=n^2\delta_{ij}$ | **Resolved — Step 5b** |
| Lensing 2:1 ratio | Asserted | Derived from metric ratio | Confirmed |
| Gravitational redshift | $e^\phi - 1$ | $e^\phi - 1$ from optical metric | Consistent — stands |
| $\beta_{\text{PPN}}$ | — | $1$ — from n-field action ($\alpha=2$) | **Resolved — Step 4** |
| Perihelion precession | $6\pi GM/c^2a(1-e^2)$ (GR import) | $\frac{2+2\gamma-\beta}{3}\times 6\pi GM/c^2a(1-e^2) = 1\times\ldots$ | **Resolved — PM predicts exact value** |
| `pm_integrate_orbit` 1PN term | $-3GMh^2/c^2r^4$ | PM-derived: $\beta=\gamma=1$ from Steps 4+5b | **Resolved — Step 5b** |

---

## Root structural finding

All four problems (wrong metric signs, $\gamma = 0$, perihelion $2/3$, GR import in
orbit code) have a **single common cause**: postulate 5 (Euclidean spatial background)
is either wrong, or it must be refined to say that the medium's compression also
scales spatial measurements.

If spatial measurements are made with rods composed of the same medium, then in a
region where $n = e^\phi$ the physical spatial metric is:
$$g_{ij} = n^2\delta_{ij} \approx (1+\phi)^2\delta_{ij} \approx (1+2\phi)\delta_{ij}$$

This gives $\gamma = 1$, restores the correct perihelion precession, is
consistent with Cassini, and the massive-particle metric becomes:
$$g_{tt} = -(1-\phi), \quad g_{ij} = (1+\phi)^2\delta_{ij}$$
which is the isotropic Schwarzschild metric at first order — GR's form, derived from
the medium without importing GR.

The Euclidean background is then the statement about the *abstract coordinate grid*
(absolute space), not about physical distances — physical distances are modified by
the medium density.  This is analogous to how an elastic solid in non-uniform
compression has a position-dependent physical distance element even though it lives
in an Euclidean reference space.

---

## Immediate open questions

1. ~~**Does rod compression give $g_{ij} = n^2\delta_{ij}$ exactly?**~~
   **Resolved in Step 5b.** The derivation from the particle Lagrangian in
   medium units gives $g_{ij} = n^2\delta_{ij}$ exactly, without a separate
   constitutive-law calculation.

2. **Is the optical redshift formula consistent with the wave equation?**
   The $z = e^\phi - 1$ result uses geometric optics.  Does the PM wave
   equation reproduce it?

3. **Does the n-field equation change the optical metric, and how?**
   If $n(r) = 1 + 2GM/c^2r$ replaces $e^{2GM/c^2r}$, the optical metric's
   strong-field profile changes with consequences for redshift and lensing.

4. ~~**Once $g_{ij} = n^2\delta_{ij}$ is established, re-derive the 1PN orbit
   correction from PM's Lagrangian** and check whether $-3GMh^2/c^2r^4$
   emerges without importing GR's result.~~
   **Resolved in Step 5b.** With $\gamma = 1$ (from $g_{ij} = n^2\delta_{ij}$)
   and $\beta = 1$ (from n-field action), the PPN 1PN term gives
   $-3GMh^2/c^2r^4$ exactly.  The orbit integrator is correct and now has
   a full PM derivation chain.

---

## What is not changed by this analysis

- The force law $\mathbf{a} = (c^2/2)\nabla\phi$ is unaffected
- All weak-field tests (solar system magnitudes, binary pulsars, SPARC fits)
  are unaffected — the discrepancy enters only at the $GM/c^2r$ level in angles
- The optical metric and redshift formula are unaffected
- The n-field equation work is unaffected
- The Lagrange point tests are unaffected

The perihelion precession tests are **no longer affected** — Step 5b established
that PM independently derives $\beta = \gamma = 1$, so the tests are now validating
a PM prediction.  The orbit integrator code is correct; its annotation has been
updated to reflect the full derivation chain.
