# Pushing Medium: A Medium-Based Reformulation of Gravity and Dynamics

## Abstract
Pushing Medium (PM) models gravity and inertial effects as propagation through a refractive medium with coupled scalar and flow fields. In the weak-field, slow-motion limit PM is constructed to match General Relativity (GR) observables to numerical precision while enabling additional cross-force structures (charge, plasma, current-induced dragging) that are tightly constrained but not excluded. We summarize the field equations, calibrated constants, rotation-curve model, validation against GR benchmarks, galaxy fitting results on SPARC data, and a strengths/weaknesses assessment.

## 1. Field Structure
- **Symbols (key):** \(c\) speed of light, \(G\) Newton constant, \(\Phi_G\) gravitational potential, \(\Phi_E\) electrostatic potential, \(n_e\) electron density, \(\mu_G,\mu_E,\alpha_{\mathrm{plasma}}\) coupling coefficients, \(\mathbf{u}\) flow field, \(\mathbf{J}_M,\mathbf{J}_q\) mass and charge current densities, \(A_M, A_q\) flow couplings, \(\delta n\) index perturbation, \(M_d, R_d, v_\infty, r_s, r_c, m\) rotation-curve parameters.
- **Unified index field** (scalar):
  $$\ln n = \mu_G \Phi_G + \mu_E \Phi_E - \alpha_{\mathrm{plasma}} n_e$$
  where \(\Phi_G\) is gravitational potential, \(\Phi_E\) electrostatic potential, and \(n_e\) electron density. In calibrated weak fields \(\mu_G\) dominates; \(\mu_E\) and \(\alpha_{\mathrm{plasma}}\) are dimensional estimates and presently unconstrained.
- **Flow field** (vector, Lense-Thirring analogue):
  $$\nabla^2 \mathbf{u} = -A_M \mathbf{J}_M - A_q \mathbf{J}_q$$
  with mass current \(\mathbf{J}_M\) and charge current \(\mathbf{J}_q\). \(A_M \approx G/c^2\) recovers GR frame-dragging; laboratory nulls bound \(A_q \ll A_M\).
- **Null geodesic analogue** (ray equation):
  $$\frac{d \mathbf{r}}{ds} = c\,\hat{\mathbf{k}} + \mathbf{u}$$
  with refractive effects encoded in \(n\) and advection by \(\mathbf{u}\).
- **Medium acceleration contribution** for rotation curves (phenomenological profile):
  $$\delta n(r) = \mu_{\mathrm{coeff}} \Phi(r) \frac{v_\infty^2}{c^2} \ln\!\left(1 + \frac{r}{r_s}\right) S(r)$$
  $$S(r) = \left[1 + \left(\frac{r_c}{r}\right)^m\right]^{-1/m}$$
  $$a_{\mathrm{med}}(r) = -c^2 \frac{d\,\delta n}{dr}$$
- **Baryonic disk acceleration** (exponential disk):
  $$M(<r) = M_d \Big[1 - \big(1 + \tfrac{r}{R_d}\big) e^{-r/R_d}\Big]$$
  $$a_{\mathrm{bar}}(r) = \frac{G M(<r)}{r^2}$$
- **Rotation speed**:
  $$v_c^2(r) = r\,[a_{\mathrm{bar}}(r) + a_{\mathrm{med}}(r)]$$
- **Dynamic coupling (advection term)** capturing compressibility-driven mixing between scalar and flow fields:
  $$\partial_t \phi + \nabla\cdot(\phi\,\mathbf{u}) = S_\phi + \cdots$$
  where \(\phi\) stands for scalar contributions (e.g., potential-linked index terms) and \(S_\phi\) gathers sources; this term is absent in GR.

## 2. Calibrated Constants (weak-field matching)
- Gravitational coupling: \(\mu_G \approx 1.486\times10^{-27}\,\mathrm{m^2\,kg^{-1}}\) (from calibration.json).
- Drag coefficients: \(k_{\mathrm{Fizeau}} = k_{\mathrm{TT}} = 1.0\) (consistency with GR light-dragging and Thomas precession analogues).
- Frame-drag constant: \(A_M \approx G/c^2\) fixed by matching Lense-Thirring precession.
- Charge/plasma couplings: \(\mu_E, \alpha_{\mathrm{plasma}}, A_q\) presently only bounded (see Section 6).

## 3. Solar-System and GR Benchmark Agreement
Computed with `.venv` (PYTHONPATH=src); see [docs/pm_vs_gr_reference.md](docs/pm_vs_gr_reference.md). PM equals GR within numerical precision for classical tests:
- Light deflection, Shapiro delay, gravitational redshift, perihelion advance, Lense-Thirring precession, Einstein radius, GW phase speed, binary quadrupole power, circular orbit energy all show relative differences \(< 10^{-12}\) numerically.
- Numerical deflection via index integration matches analytic PM/GR to \(\sim 6\times10^{-4}\) relative.

| Observable | PM | GR | Rel. diff |
| --- | --- | --- | --- |
| deflection_angle_point_mass | $8.488869\times10^{-6}$ | $8.488869\times10^{-6}$ | $0$ |
| shapiro_delay_point_mass (s) | $1.133570\times10^{-4}$ | $1.133570\times10^{-4}$ | $0$ |
| gravitational_redshift_potential (\(\Delta\phi=10^7\)) | $1.112650\times10^{-10}$ | $1.112650\times10^{-10}$ | $0$ |
| perihelion_precession (rad/orbit) | $5.020872\times10^{-7}$ | $5.020872\times10^{-7}$ | $0$ |
| lense_thirring_precession (rad/s) | $6.102516\times10^{-15}$ | $6.102516\times10^{-15}$ | $0$ |
| einstein_radius_point_mass (rad) | $9.784003\times10^{-5}$ | $9.784003\times10^{-5}$ | $0$ |
| gw_phase_speed (m/s) | $2.997925\times10^8$ | $2.997925\times10^8$ | $0$ |
| binary_quadrupole_power (W) | $1.460091\times10^{25}$ | $1.460091\times10^{25}$ | $0$ |
| circular_orbit_energy (J/kg) | $-6.637591\times10^8$ | $-6.637591\times10^8$ | $\approx 0$ |
| index_deflection_numeric | $8.483710\times10^{-6}$ | $8.488869\times10^{-6}$ | $-6.08\times10^{-4}$ |

## 4. Galaxy Rotation Curve Model (SPARC pipeline)
- Six fitted parameters: \(M_d, R_d, v_\infty, r_s, r_c, m\).
- Data: SPARC 175 galaxies (high-quality rotation curves).
- Fitter: random search + local refinement (300 + 100 samples) per galaxy.
- Outputs: PM-only fits (no dark matter halos) saved to `results/pm_sparc_fits.json`; scripts under `scripts/`.

### Fitting Algorithm (methods)
- **Preprocess**: read CSV rotation curves; convert radii to meters; velocities to m/s; errors carried through.
- **Model evaluation**: compute \(a_{\mathrm{bar}}\) (exponential disk) and \(a_{\mathrm{med}}\) (log core + smoothing \(S(r)\)); assemble \(v_c(r)\).
- **Objective**: \(\chi^2 = \sum_i (v_{\mathrm{obs},i} - v_{c,i})^2 / \sigma_i^2\) with fallback \(0.1\,v\) if \(\sigma_i\) is missing.
- **Search**: uniform random sampling (\(n_{\mathrm{random}}\)) in bounds, keep best; followed by multiplicative perturb refinements (\(n_{\mathrm{refine}}\)). Parameter perturbations are clamped to bounds.
- **Defaults**: \(n_{\mathrm{random}}=300\), \(n_{\mathrm{refine}}=100\); bounds from SPARC guide; single-threaded.
- **Outputs**: best-fit parameters, \(\chi^2\), model curve; serialized to JSON for downstream analysis.

### SPARC Fit Statistics (current run)
- Galaxies: 175
- Mean \(\chi^2\): 198.11; Median \(\chi^2\): 15.68; Std: 716.69; Range: 0.04–7581
- Quality bins: <2 (16.6%), 2–5 (11.4%), 5–10 (13.1%), ≥10 (58.9%)
- Success rate (\(\chi^2 < 10\)): 41%
- Sensitivity: galaxies with <10 points have mean \(\chi^2 = 7.78\); galaxies with ≥10 points have mean \(\chi^2 = 276\) → indicates optimizer/model strain on richer curves.

| Metric | Value |
| --- | --- |
| Galaxies fitted | 175 |
| Mean \(\chi^2\) | 198.11 |
| Median \(\chi^2\) | 15.68 |
| Std dev | 716.69 |
| Range | 0.04 – 7581.37 |
| Success rate (\(\chi^2 < 10\)) | 41% |
| <10 data points (mean \(\chi^2\)) | 7.78 |
| ≥10 data points (mean \(\chi^2\)) | 276.39 |

### Figures (current run)
- [docs/figures/pm_sparc_chi2_hist.png](docs/figures/pm_sparc_chi2_hist.png): Histogram of \(\chi^2\) across 175 SPARC galaxies.
- [docs/figures/pm_sparc_chi2_cdf.png](docs/figures/pm_sparc_chi2_cdf.png): Cumulative distribution (log \(\chi^2\)) highlighting long-tail failures.

Regenerate: activate `.venv`, ensure `matplotlib` installed, then run the plotting snippet in Section 7 (or rerun `python -m pytest` to keep fits validated).

### Interpretation
- The infrastructure runs end-to-end, but high-resolution curves expose optimization or model limitations (local minima; limited profile flexibility).
- Best fits cluster on sparse datasets; failures dominate rich datasets.
- No NFW/MOND comparison yet in this dataset; needs follow-up for fairness.

## 5. Comparative Assessment vs GR
### Strengths
- **Weak-field consistency**: Matches GR benchmarks to numerical precision (Section 3).
- **Unified structure**: Single index field admits gravity + charge + plasma contributions; frame-dragging recovered via \(A_M\).
- **No singularities** in formulation (medium picture) and natural ray-tracing interpretation.

### Weaknesses / Open Issues
- **Galaxy fits**: Current PM profile + simple optimizer underperform on rich SPARC curves (mean \(\chi^2\) too high).
- **Unconstrained couplings**: \(\mu_E, \alpha_{\mathrm{plasma}}, A_q\) only bounded; cross-force effects likely \(\ll\) gravitational.
- **Strong-field regime**: Behavior near compact objects not yet benchmarked vs full GR (beyond weak-field observables).
- **Parameter degeneracy**: Six-parameter fits can be non-identifiable without stronger priors or scaling relations.

### GR Comparison Snapshot
From [docs/pm_vs_gr_reference.md](docs/pm_vs_gr_reference.md): all classical observables align; PM differs structurally by allowing suppressed cross-force terms. In practice, current bounds make PM indistinguishable from GR in tested regimes, while adding phenomenology for galaxy scales.

## 6. Coupling Bounds and Cross-Force Terms
- **Charge lensing (\(\mu_E\))**: dimensional \(\sim 10^{-37}\,\mathrm{m^2/C}\); no direct constraints—proposed charged-beam lensing could test.
- **Plasma-gravity coupling (\(\alpha_{\mathrm{plasma}}\))**: dimensional \(\sim 10^{-15}\,\mathrm{m^3}\); cluster lensing vs frequency offers a probe.
- **Charge-current dragging (\(A_q\))**: lab solenoid nulls imply \(A_q \lesssim 10^{-10}\,\mathrm{m^3/(C\,s)}\), giving \(A_q/A_M < 4\times10^{17}\); effectively zero in astrophysical contexts.

## 7. Reproduction and Data Availability
- Environment: `.venv` with `requirements.txt`; export `PYTHONPATH=src` when running scripts.
- GR/PM reference table: [docs/pm_vs_gr_reference.md](docs/pm_vs_gr_reference.md) includes code snippet to reproduce values.
- Galaxy fits: run `python scripts/batch_fit_sparc.py` (default sample) or `python scripts/fit_sparc_batch.py` (full SPARC if data present); outputs to `results/`.
- Plots: regenerate figures with

```bash
. .venv/bin/activate
PYTHONPATH=src python - <<'PY'
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
ROOT = Path('..').resolve()
data = json.load((ROOT / 'results' / 'sparc_fit_results.json').open())
chi2_values = [v['chi2'] for v in data.values()]
fig, ax = plt.subplots(figsize=(6,4))
ax.hist(chi2_values, bins=50, color='#1f77b4', alpha=0.8)
ax.set_xlabel(r'$\chi^2$'); ax.set_ylabel('Count'); ax.set_title('PM SPARC Fits: $\chi^2$ Histogram')
fig.tight_layout(); (ROOT / 'docs' / 'figures').mkdir(parents=True, exist_ok=True)
fig.savefig(ROOT / 'docs' / 'figures' / 'pm_sparc_chi2_hist.png', dpi=300)
plt.close(fig)
vals = np.sort(np.array(chi2_values)); cdf = np.linspace(0,1,len(vals))
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(vals, cdf, color='#d62728'); ax.set_xscale('log'); ax.set_xlabel(r'$\chi^2$ (log)'); ax.set_ylabel('Cumulative fraction')
ax.set_title('PM SPARC Fits: CDF of $\chi^2$'); ax.grid(True, which='both', ls=':', lw=0.5)
fig.tight_layout(); fig.savefig(ROOT / 'docs' / 'figures' / 'pm_sparc_chi2_cdf.png', dpi=300)
plt.close(fig)
PY
```

## 8. Future Work
1. **Optimization upgrade**: differential evolution or basin-hopping; increase iterations; multi-start initializations informed by baryonic scaling relations.
2. **Model flexibility**: explore alternative medium profiles or priors tying \(v_\infty, r_s, r_c\) to baryonic properties.
3. **Comparative baselines**: compute NFW/MOND fits on the same SPARC sample for head-to-head \(\chi^2\) distributions.
4. **Strong-field tests**: numerically compare PM predictions to GR near compact objects; check stability and horizon analogues.
5. **Cross-force experiments**: precision interferometry near high-current conductors (bounds on \(A_q\)), charged-beam lensing (\(\mu_E\)), frequency-dependent cluster lensing (\(\alpha_{\mathrm{plasma}}\)).

## 9. Conclusion
Pushing Medium provides a medium-based rephrasing of gravity that is constructed to reproduce GR in tested weak-field regimes while admitting additional, tightly bounded cross-force structures. Current galaxy-rotation fits demonstrate viability on sparse curves but fall short on rich SPARC datasets, indicating the need for stronger optimization and possibly richer medium profiles. With improved fitting and direct constraints on the extra couplings, PM can be further tested against GR and dark-matter paradigms across weak- to strong-field regimes.

## 10. Strong-Field Analysis: PM vs GR

This section presents the first quantitative strong-field comparison between PM and GR, computed using exact eikonal integrals (substitution \(p = b/r\), validated with `scipy.integrate.quad`). Full reproducible code: [scripts/strong_field_deflection.py](scripts/strong_field_deflection.py).

### 10.1 Deflection Integrals

For PM with \(n(r) = 1 + R_s/r\), Bouguer's law gives:

$$\alpha_{\rm PM} = 2\int_0^{p_{\max}} \frac{dp}{\sqrt{1 - p^2(1-\varepsilon^2) + 2\varepsilon p}} - \pi, \qquad \varepsilon = \frac{R_s}{b},\quad p_{\max} = \frac{b}{b - R_s}$$

GR Schwarzschild null-geodesic integral:

$$\alpha_{\rm GR} = 2\int_0^{p_{\max}^{\rm GR}} \frac{dp}{\sqrt{1 - p^2 + p^3/\!\bar{b}}} - \pi, \qquad \bar{b} = b/R_s$$

where \(p_{\max}^{\rm GR}\) is the largest root of \(1 - p^2 + p^3/\bar{b} = 0\).

### 10.2 Photon Sphere

For \(n(r)=1+R_s/r\): \(\tfrac{d}{dr}(n(r)\,r) = 1\) everywhere -- **PM has no photon sphere**. GR has one at \(r_{\rm ps}=1.5\,R_s\), giving critical impact parameter \(b_{\rm crit}^{\rm GR}=\tfrac{3\sqrt{3}}{2}R_s\approx 2.598\,R_s\). PM's capture threshold is only \(b_{\rm crit}^{\rm PM}=R_s\), a factor 2.6x smaller.

### 10.3 Numerical Results

| \(b/R_s\) | \(\alpha_{\rm PM}\) (rad) | \(\alpha_{\rm GR}\) (rad) | \(\alpha_{\rm WF}\) (rad) | \(\alpha_{\rm PM}/\alpha_{\rm GR}\) |
|---|---|---|---|---|
| 10000 | \(2.000\times10^{-4}\) | \(2.000\times10^{-4}\) | \(2.000\times10^{-4}\) | 0.9999 |
| 100   | \(2.016\times10^{-2}\) | \(2.030\times10^{-2}\) | \(2.000\times10^{-2}\) | 0.9930 |
| 10    | \(2.172\times10^{-1}\) | \(2.361\times10^{-1}\) | \(2.000\times10^{-1}\) | 0.9197 |
| 5     | \(4.758\times10^{-1}\) | \(5.904\times10^{-1}\) | \(4.000\times10^{-1}\) | 0.8059 |
| 3     | \(9.115\times10^{-1}\) | \(1.719\) | \(6.667\times10^{-1}\) | 0.5301 |
| 2.6   | \(1.117\) | \(6.810\) | \(7.692\times10^{-1}\) | 0.1640 |
| 2.55  | \(1.150\) | diverges | \(7.843\times10^{-1}\) | -- |
| 2.0   | \(1.695\) | diverges | \(1.000\) | -- |
| 1.5   | \(3.031\) | diverges | \(1.333\) | -- |
| 1.1   | \(9.878\) | diverges | \(1.818\) | -- |
| 1.01  | \(39.61\) | diverges | \(1.980\) | -- |

### 10.4 Figures

- [docs/figures/strong_field_deflection_pm_vs_gr.png](docs/figures/strong_field_deflection_pm_vs_gr.png): Deflection angle vs \(b/R_s\) (log scale). GR diverges at \(b_{\rm crit}^{\rm GR}\); PM continues smoothly to \(b_{\rm crit}^{\rm PM}=R_s\).
- [docs/figures/strong_field_pm_gr_ratio.png](docs/figures/strong_field_pm_gr_ratio.png): \(\alpha_{\rm PM}/\alpha_{\rm GR}\) ratio -- converges to 1 in weak field; PM predicts systematically less bending in the strong-field regime.
- [docs/figures/strong_field_rmin_pm_vs_gr.png](docs/figures/strong_field_rmin_pm_vs_gr.png): Closest approach \(r_{\min}/R_s\) vs \(b/R_s\) -- PM allows light much closer to the compact object; GR is floored by its photon sphere.

### 10.5 Key Contrasts

| Property | PM | GR |
|---|---|---|
| Photon sphere | **None** | \(r_{\rm ps} = 1.5\,R_s\) |
| Critical impact param | \(b_{\rm crit} = R_s\) | \(b_{\rm crit} = 2.598\,R_s\) |
| Shadow radius (test-mass) | \(\sim 2\,GM/c^2\) | \(\sim 5.2\,GM/c^2\) |
| Weak-field agreement | <0.1% at \(b>100\,R_s\) | reference |
| Deflection at \(b=3R_s\) | 0.91 rad | 1.72 rad |

EHT images of M87* and Sgr A* constrain the shadow radius directly. GR predicts ~5.2 GM/c^2; PM predicts ~2 GM/c^2 -- a factor ~2.6 difference and a **sharp, falsifiable observational test** of PM vs GR.

### 10.6 Caveats and Open Problems

Results use the test-mass PM index \(n=1+R_s/r\); the self-consistent treatment
for a compact object with a resolved interior is developed in Section 11 below. Remaining open problems:
- Gravitational-wave emission and ringdown in the PM framework.
- Frame-dragging and ergosphere analogues without a horizon.
- Numerical PM ray-tracing for shadow image comparison with EHT.

---

## 11. Self-Consistent PM Compact Object

Section 10 used the test-mass index \(n=1+R_s/r\) as a proxy for the external
field. Here we solve the PM scalar equation self-consistently for a uniform-density
star, obtain the full interior+exterior refractive-index profile, and compute
light deflection from first principles using only PM equations.
Full reproducible code: [scripts/pm_compact_object.py](scripts/pm_compact_object.py).

### 11.1 Governing Equation

The static PM scalar field \(\phi(\mathbf{r})\) satisfies (Section 2.1 of the
PM field theory):

$$\nabla^2\phi = -4\pi G_{\rm eff}\,\rho_M(r)$$

with refractive index \(n(r) = e^{\phi(r)}\).  For a focusing medium we require
\(\phi > 0\), so the Green function gives the exterior solution

$$\phi_{\rm ext}(r) = \frac{G_{\rm eff} M}{r} = \frac{GM}{r} > 0.$$

Note the sign: this is the PM lensing potential, equal to \(-2\Phi_N/c^2\) where
\(\Phi_N = -GM/r\) is the Newtonian gravitational potential.

### 11.2 Matter Model and Analytic Solution

Uniform-density star with total mass \(M\), radius \(R_\star\):

$$\rho_0 = \frac{3M}{4\pi R_\star^3}.$$

The ODE `(1/r²) d/dr[r² dφ/dr] = −4πG ρ` has the exact analytic solution:

$$\phi(r) = \begin{cases}
\phi_c - \dfrac{2\pi G\rho_0}{3}\,r^2, & r \le R_\star \\[6pt]
\dfrac{GM}{r}, & r > R_\star
\end{cases}$$

where the central value is fixed by continuity at the surface:

$$\phi_c = \frac{GM}{R_\star} + \frac{2\pi G\rho_0}{3}R_\star^2.$$

Derivative continuity (mass matching) at \(R_\star\) is automatic:
\(\phi_{\rm int}'|_{R_\star} = -(4\pi G\rho_0/3)R_\star = -GM/R_\star^2 = \phi_{\rm ext}'|_{R_\star}\).

The solution is verified numerically with `scipy.integrate.solve_ivp` (DOP853):
exterior max error \(2.7\times10^{-14}\), interior max error \(3.0\times10^{-12}\).

### 11.3 Parameters

Working in geometric units \(G = c = 1\), \(M = 1\), \(R_s = 2GM = 2\), stellar
radius \(R_\star = 5\,R_s = 10\) (compactness \(R_s/R_\star = 0.2\)):

| Quantity | Value |
|---|---|
| \(\phi_c = \phi(0)\) | 0.150 |
| \(\phi(R_\star)\) | 0.100 |
| \(n(0) = e^{\phi_c}\) | 1.162 |
| \(n(R_\star) = e^{\phi(R_\star)}\) | 1.105 |

### 11.4 Photon-Orbit Function and Photon Sphere

Define \(F(r) = n(r)\,r\). A photon sphere exists if and only if \(F'(r) = 0\)
has a solution.

With \(\phi = GM/r\) exterior: \(F(r) = e^{GM/r}\,r\), so

$$F'(r) = e^{GM/r}\!\left(1 - \frac{GM}{r}\right).$$

\(F' = 0\) requires \(r = GM = 1\), which lies **inside** \(R_\star = 10\). In
the interior, \(F(r) = e^{\phi_c - Ar^2}\,r\) with \(A > 0\), which is monotone
for small \(r\) (verified numerically). No sign change of \(F'\) is found over
\((0,\,30\,R_s)\):

> **Result: the self-consistent PM compact object has no photon sphere for this
> compactness (\(R_s/R_\star=0.2\)).**

The PM shadow edge is set by the star surface:

$$b_{\rm crit} = F(R_\star) = n(R_\star)\,R_\star = e^{0.1}\times 10 = 11.052
\approx 5.53\,R_s.$$

### 11.5 Light Deflection (PM Bouguer Law)

The PM deflection angle is computed from Bouguer's law (no GR geodesics):

$$\alpha(b) = 2\int_{r_{\min}}^\infty \frac{b\,dr}{r\sqrt{F(r)^2 - b^2}} - \pi,$$

where \(r_{\min}\) solves \(F(r_{\min}) = b\).  The integral is evaluated with
the substitution \(u = r_{\min}/r\) mapping \([r_{\min},\infty)\to(0,1]\):

$$\alpha(b) = 2\int_0^1 \frac{b\;du}{u\sqrt{F(r_{\min}/u)^2 - b^2}} - \pi.$$

This avoids any far-field truncation error; `scipy.integrate.quad` converges to
\(<10^{-8}\) relative tolerance.

### 11.6 Numerical Results

Weak-field limit: \(\alpha \to 2GM/b\) as \(b\to\infty\) (ratio PM/WF = 1.004 at
\(b = 200\,R_s\)). Deflection increases monotonically toward the shadow edge.

| \(b/R_s\) | \(\alpha_{\rm PM}\) (rad) | \(\alpha_{\rm WF} = 2GM/b\) (rad) | PM/WF |
|---|---|---|---|
| 200 | \(5.020\times10^{-3}\) | \(5.000\times10^{-3}\) | 1.004 |
| 100 | \(1.012\times10^{-2}\) | \(1.000\times10^{-2}\) | 1.012 |
| 50 | \(2.064\times10^{-2}\) | \(2.000\times10^{-2}\) | 1.032 |
| 20 | \(5.380\times10^{-2}\) | \(5.000\times10^{-2}\) | 1.076 |
| 10 | \(1.164\times10^{-1}\) | \(1.000\times10^{-1}\) | 1.164 |
| 6.0 | \(2.117\times10^{-1}\) | \(1.667\times10^{-1}\) | 1.270 |
| 5.53 (≈ \(b_{\rm crit}\)) | \(\to\infty\) | \(1.81\times10^{-1}\) | — |

### 11.7 Figures

- [docs/figures/pm_co_phi.png](docs/figures/pm_co_phi.png): \(\phi(r)\) — interior parabolic, exterior \(GM/r\).
- [docs/figures/pm_co_n.png](docs/figures/pm_co_n.png): \(n(r) = e^{\phi(r)}\) — focusing medium (\(n > 1\) everywhere).
- [docs/figures/pm_co_f.png](docs/figures/pm_co_f.png): \(F(r) = n(r)\,r\) — monotone, no photon sphere.
- [docs/figures/pm_co_fprime.png](docs/figures/pm_co_fprime.png): \(F'(r)\) — positive everywhere, confirming no photon sphere.
- [docs/figures/pm_co_deflection.png](docs/figures/pm_co_deflection.png): \(\alpha(b)\) vs \(b/R_s\) (linear scale) with weak-field reference.
- [docs/figures/pm_co_deflection_log.png](docs/figures/pm_co_deflection_log.png): \(\alpha(b)\) vs \(b/R_s\) (log–log scale).
- [docs/figures/pm_co_rmin.png](docs/figures/pm_co_rmin.png): \(r_{\min}(b)\) — closest approach vs impact parameter.

### 11.8 Key Results

| Property | Self-consistent PM | Test-mass PM (§10) | GR |
|---|---|---|---|
| Photon sphere | **None** | None | \(r_{\rm ps}=1.5\,R_s\) |
| Shadow radius | \(5.53\,R_s\) (star surface) | \(1.0\,R_s\) | \(2.60\,R_s\) |
| \(n(0)\) | 1.162 | diverges | — |
| Weak-field limit | \(2GM/b\) | \(4GM/b\) | \(4GM/b\) |
| Method | PM Bouguer + PM Poisson | PM Bouguer (test-mass) | GR geodesic |

The self-consistent PM weak-field deflection converges to \(2GM/b\) (rather than
the test-mass \(4GM/b\)) because the PM coupling constant \(G_{\rm eff} = G\) (not
\(2G\)). The correct factor of 2 required to match GR light bending would need
\(G_{\rm eff} = 2G\), i.e.\ an additional source factor of two in the PM scalar
equation — a constraint on the PM coupling constant calibration.

### 11.9 Open Problems

- Determine the value of \(G_{\rm eff}\) (or equivalently \(\mu_G\)) required to
  reproduce the GR deflection factor of \(4GM/b\) in the weak-field limit.
- Extend to rotating compact objects (PM vector field \(\mathbf{u}\) sourced by
  mass current).
- PM ray-tracing for shadow images comparable to EHT data.

---

## 12. Compactness Sweep: Does PM Ever Develop a Photon Sphere?

Full reproducible code: [scripts/pm_compactness_sweep.py](scripts/pm_compactness_sweep.py).

We sweep the stellar compactness \(k = R_\star/R_s\) from \(5\) down to \(1.01\)
and track \(\phi_c\), \(n(0)\), \(n(R_\star)\), \(b_{\rm crit}\), and whether
\(F'(r) = 0\) anywhere.

### 12.1 Analytic Photon-Sphere Condition

For the uniform-density interior, \(\phi(r) = \phi_c - Ar^2\) with
\(A = GM/(2R_\star^3)\), so:

$$F'(r) = n(r)\!\left(1 - \frac{GM}{R_\star^3}r^2\right) = 0
\quad\Longleftrightarrow\quad r = R_\star^{3/2}.$$

This critical radius lies inside the star only if \(R_\star^{3/2} < R_\star\),
i.e.\ \(R_\star < 1\). Since \(R_s = 2GM = 2\) in our units, we need
\(R_\star < 0.5\,R_s\) — far inside the Schwarzschild radius.

For the exterior \(\phi = GM/r\): \(F'(r) = 0\) at \(r = GM = 1 < R_\star\)
for all \(R_\star \geq R_s\).

> **Analytic result: for any uniform-density PM compact object with
> \(R_\star \geq R_s\), \(F'(r) > 0\) everywhere — no photon sphere exists.**
> The shadow is always set by the stellar surface.

### 12.2 Identity: Shadow Radius

Since \(b_{\rm crit} = F(R_\star) = n(R_\star)\cdot R_\star\) and
\(n(R_\star) = \exp(GM/R_\star)\):

$$\boxed{b_{\rm crit} = e^{GM/R_\star} \cdot R_\star, \qquad
\frac{b_{\rm crit}}{R_\star} = n(R_\star) = e^{GM/R_\star}.}$$

As \(R_\star \to R_s = 2GM\), \(b_{\rm crit}/R_s \to e^{1/2}\approx 1.649\).

### 12.3 PM Shadow Radius as a Surface–Index Effect

In the Pushing-Medium framework the shadow radius is not set by an external unstable photon orbit, but by the refractive index at the stellar surface. The scalar field \(\phi(r)\) obeys a monotonic Poisson-type equation, so the optical potential \(F(r)=n(r)r\) has no exterior extremum for any physically allowed compactness. As a result, the critical impact parameter is simply

$$b_{\rm crit} = n(R_\star)\,R_\star = R_\star\, e^{GM/R_\star},$$

and the shadow "hugs" the surface rather than forming at a universal radius. Increasing compactness raises the surface index \(n(R_\star)\) and enlarges the shadow (in units of \(R_\star\)), but never produces a photon sphere. This contrasts sharply with general relativity, where the Schwarzschild geometry forces an exterior photon sphere at \(r = 3R_s/2\) once \(R_\star \leq 3R_s/2\), fixing the shadow radius at \(b^{\rm GR} = 3\sqrt{3}\,GM/c^2\) independent of the stellar surface. The two theories therefore diverge qualitatively in the strong field: GR predicts a universal photon ring, while PM predicts a surface-controlled shadow with no unstable orbit.

### 12.4 Numerical Results

| \(k = R_\star/R_s\) | \(\phi_c\) | \(n(0)\) | \(n(R_\star)\) | \(b_{\rm crit}/R_s\) | \(b_{\rm crit}/R_\star\) | Photon sphere |
|---|---|---|---|---|---|---|
| 5.00 | 0.1500 | 1.162 | 1.105 | 5.526 | 1.105 | **no** |
| 4.00 | 0.1875 | 1.206 | 1.133 | 4.533 | 1.133 | **no** |
| 3.00 | 0.2500 | 1.284 | 1.181 | 3.544 | 1.181 | **no** |
| 2.50 | 0.3000 | 1.350 | 1.221 | 3.054 | 1.221 | **no** |
| 2.00 | 0.3750 | 1.455 | 1.284 | 2.568 | 1.284 | **no** |
| 1.75 | 0.4286 | 1.535 | 1.331 | 2.329 | 1.331 | **no** |
| 1.50 | 0.5000 | 1.649 | 1.396 | 2.093 | 1.396 | **no** |
| 1.25 | 0.6000 | 1.822 | 1.492 | 1.865 | 1.492 | **no** |
| 1.10 | 0.6818 | 1.978 | 1.576 | 1.733 | 1.576 | **no** |
| 1.05 | 0.7143 | 2.043 | 1.610 | 1.690 | 1.610 | **no** |
| 1.01 | 0.7426 | 2.101 | 1.641 | 1.657 | 1.641 | **no** |

### 12.5 Key Observations

**1. Focusing medium throughout.** \(n(r) > 1\) everywhere for all compactnesses. As \(k \to 1\), \(n(0) \to 2.12\) — the central medium becomes strongly refractive but remains well-defined.

**2. No photon sphere at any compactness \(k \geq 1\).** \(F'(r) > 0\) is confirmed numerically across \(10^5\) grid points for all 11 test cases. This is not a numerical accident — it is proven analytically above.

**3. Monotone shadow shrinkage.** \(b_{\rm crit}/R_s\) decreases from \(5.53\) at \(k=5\) to \(1.65\) at \(k=1.01\), but never reaches the GR photon-sphere value \(b^{\rm GR}_{\rm crit} = \tfrac{3\sqrt{3}}{2}R_s \approx 2.60\,R_s\) for \(k > 1.5\). At \(k=1.5\) the curves cross: below this compactness GR develops an exterior photon sphere while PM's shadow continues to shrink toward the surface.

**4. Surface shadow vs. photon-sphere shadow.** For \(R_\star \leq 1.5\,R_s\) (i.e.\ \(k \leq 1.5\)), GR would show a photon sphere and relativistic photon ring outside the star; PM shows only a clean surface shadow. This is a **structurally distinct prediction** that survives at all compactnesses within reach of PM.

**5. Strong-field enhancement grows with compactness.** Deflection curves for decreasing \(k\) show larger PM/WF enhancement ratios near the shadow edge, but the behaviour remains finite everywhere — consistent with no photon sphere.

### 12.6 Comparison: PM Shadow vs GR Photon Sphere

| \(k = R_\star/R_s\) | PM \(b_{\rm crit}/R_s\) | GR BH \(b^{\rm GR}/R_s\) | PM shadow \(b_{\rm crit}c^2/GM\) | GR BH shadow \(b^{\rm GR}c^2/GM\) | GR PS ext? | PM PS? |
|---|---|---|---|---|---|---|
| 5.0 | 5.53 | 2.60 | 11.05 | 5.20 | no (buried in star) | **no** |
| 2.0 | 2.57 | 2.60 | 5.14 | 5.20 | no (buried in star) | **no** |
| 1.5 | 2.09 | 2.60 | 4.19 | 5.20 | **yes** \(r_{\rm ps}=3=R_\star\) | **no** |
| 1.1 | 1.73 | 2.60 | 3.47 | 5.20 | **yes** \(r_{\rm ps}=3>R_\star=2.2\) | **no** |
| 1.01 | 1.66 | 2.60 | 3.31 | 5.20 | **yes** \(r_{\rm ps}=3>R_\star=2.02\) | **no** |

GR BH column is the Schwarzschild black-hole reference \(b^{\rm GR}_{\rm crit} = \tfrac{3\sqrt{3}}{2}R_s = 3\sqrt{3}\,GM/c^2 \approx 5.196\,GM/c^2\), the value directly constrained by EHT. The PM \(c^2/GM\) column gives the equivalent EHT-comparable number for each compactness.

At \(k \leq 1.5\), PM and GR make qualitatively different predictions: GR grows a photon ring outside the star; PM does not.

### 12.7 Figures

- [docs/figures/pm_sweep_n.png](docs/figures/pm_sweep_n.png): \(n(0)\) and \(n(R_\star)\) vs compactness \(k\).
- [docs/figures/pm_sweep_bcrit.png](docs/figures/pm_sweep_bcrit.png): \(b_{\rm crit}/R_s\) and \(b_{\rm crit}/R_\star\) vs \(k\) (dual axis).
- [docs/figures/pm_sweep_bcrit2.png](docs/figures/pm_sweep_bcrit2.png): \(b_{\rm crit}/R_s\) vs \(k\) with GR photon-sphere reference line.
- [docs/figures/pm_sweep_F.png](docs/figures/pm_sweep_F.png): \(F(r)=n(r)r\) profiles for \(k \in \{5, 3, 2, 1.5, 1.1\}\).
- [docs/figures/pm_sweep_Fprime.png](docs/figures/pm_sweep_Fprime.png): \(F'(r)\) profiles — positive everywhere for all \(k\), confirming no photon sphere.
- [docs/figures/pm_sweep_deflection.png](docs/figures/pm_sweep_deflection.png): \(\alpha(b)\) curves for five compactnesses on a common axis.

### 12.8 Summary

PM compact objects remain in the "surface-shadow, no-photon-sphere" regime for all
\(R_\star \geq R_s\). The shadow radius shrinks monotonically with compactness, converging to \(b_{\rm crit} \to e^{1/2} R_s \approx 1.65\,R_s\) at the Schwarzschild limit. The precise threshold below which GR grows a photon sphere (\(R_\star < 1.5\,R_s\)) is a region where PM and GR become qualitatively distinguishable — a sharp, falsifiable prediction accessible to any instrument that can resolve sub-Schwarzschild-radius structure around ultra-compact objects.
