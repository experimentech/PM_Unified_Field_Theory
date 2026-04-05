"""
Derrick scaling argument for PM soliton stability
==================================================

Applies the Derrick (1964) theorem to the Pushing Medium field equations to
determine whether stable, static, localised soliton solutions are possible.

Background
----------
PM's vacuum field carries a medium-compression field n(r) = e^φ(r).  A
"nucleon" in PM is expected to be a localised, stable compression (or
rarefaction) of the medium — a soliton held in shape by the balance between
gradient energy (wants to spread out) and potential energy (wants to restore
n → 1, the vacuum).

PM action (α=2 measure, Section 2 of pm-metric-from-first-principles.md):

    S = ∫ [½|∇φ|² n²  −  A·V̂(n)] n^{α-2} d³x   with α = 2

    = ∫ [½|∇φ|² n²  −  A·V̂(n)] d³x

With φ = ln n, |∇φ|² = |∇n|²/n², so the gradient density becomes
½|∇n|²/n² · n² = ½|∇n|².

The vacuum-subtracted potential (fixed point at n=1):

    V̂_vac(n) = ½n² + ln n − ½ − 2(n − 1)

    V̂_vac'(n) = n + 1/n − 2 = (n−1)²/n ≥ 0   (global minimum at n=1)

Static energy functional for a spherically symmetric soliton n(r):

    E[n] = I_grad + I_pot

    I_grad = 4π ∫₀^∞ ½(dn/dr)² r² dr       (always ≥ 0)
    I_pot  = 4π ∫₀^∞ V̂_vac(n(r)) n(r)² r² dr

Derrick scaling
---------------
Under coordinate rescaling  n_λ(r) = n(λr)  the energy transforms as:

    E(λ) = λ^{-1} I_grad  +  λ^{-3} I_pot

For a stationary point (candidate soliton): dE/dλ|_{λ=1} = 0

    ⟹  I_grad + 3 I_pot = 0    [DERRICK CONDITION]

Stability of that stationary point: d²E/dλ²|_{λ=1} > 0

    d²E/dλ² = 2 I_grad + 12 I_pot
             = 2 I_grad + 12(−I_grad/3)   (using Derrick condition)
             = 2 I_grad − 4 I_grad
             = −2 I_grad  <  0            (since I_grad ≥ 0)

CONCLUSION: If a Derrick stationary point exists, it is always a SADDLE in
the scaling direction — the soliton is classically unstable to uniform
dilation/contraction.  PM's classical field equations do not support stable
static solitons.

Sign analysis
-------------
Compression solitons (n > 1 throughout):
    V̂_vac(n) > 0  ⟹  I_pot > 0
    Derrick condition requires I_pot < 0 — VIOLATED.
    No stationary point at all; the soliton simply spreads.

Rarefaction solitons (n < 1 somewhere):
    V̂_vac(n) < 0  ⟹  I_pot can be negative
    Derrick condition CAN be satisfied (stationary point exists)
    But d²E/dλ² = −2 I_grad < 0  → scale-unstable (collapses inward)

Physical meaning
----------------
A stable nucleon analogue in PM requires either:
  (a) A topological charge that prevents the field from unwinding
      (needs non-trivial π₂ or π₃ of the target space in n-field)
  (b) Multiple coupled fields (e.g., n + a charged component)
  (c) A quantum length scale that stabilises the size externally
      (the "quantum foam" route — sets the minimum excitation scale)

All three are outside the classical single-field PM action.  This confirms
that PM's classical sector cannot fix ρ_nuc without one external input.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad


# ---------------------------------------------------------------------------
# Potential
# ---------------------------------------------------------------------------

def v_vac(n: np.ndarray | float) -> np.ndarray | float:
    """
    PM vacuum potential (subtracted so n=1 is fixed point at zero):

        V̂_vac(n) = ½n² + ln(n) − ½ − 2(n − 1)

    Properties:
        V̂_vac(1) = 0
        V̂_vac'(n) = (n−1)²/n ≥ 0  →  global minimum at n = 1
        V̂_vac(n) < 0  for  n ∈ (0, 1)
        V̂_vac(n) > 0  for  n > 1
    """
    n = np.asarray(n, dtype=float)
    return 0.5 * n**2 + np.log(n) - 0.5 - 2.0 * (n - 1.0)


def v_vac_prime(n: np.ndarray | float) -> np.ndarray | float:
    """
    V̂_vac'(n) = n + 1/n − 2 = (n−1)²/n  ≥ 0  for all n > 0.
    """
    n = np.asarray(n, dtype=float)
    return n + 1.0 / n - 2.0


# ---------------------------------------------------------------------------
# Energy integrals — spherically symmetric field n(r) on radial grid
# ---------------------------------------------------------------------------

def energy_integrals(
    n_func,
    r_max: float = 20.0,
    n_pts: int = 4000,
):
    """
    Compute (I_grad, I_pot, E) for a spherically symmetric trial soliton.

    Parameters
    ----------
    n_func : callable
        n(r) — compression field profile.  Must satisfy n(r) → 1 as r → ∞.
    r_max : float
        Outer truncation radius (in units where the soliton width ~ 1).
    n_pts : int
        Number of radial grid points.

    Returns
    -------
    I_grad : float
        Gradient energy  4π ∫ ½(dn/dr)² r² dr  ≥ 0
    I_pot : float
        Potential energy  4π ∫ V̂_vac(n) n² r² dr
    E : float
        Total energy  I_grad + I_pot
    """
    r = np.linspace(1e-6, r_max, n_pts)
    n = np.vectorize(n_func)(r)
    dndr = np.gradient(n, r)

    integrand_grad = 4.0 * np.pi * 0.5 * dndr**2 * r**2
    integrand_pot  = 4.0 * np.pi * v_vac(n) * n**2 * r**2

    dr = r[1] - r[0]
    I_grad = float(np.sum(integrand_grad) * dr)
    I_pot  = float(np.sum(integrand_pot)  * dr)

    return I_grad, I_pot, I_grad + I_pot


def derrick_energy(
    lam: np.ndarray | float,
    I_grad: float,
    I_pot: float,
) -> np.ndarray | float:
    """
    Scaled energy under Derrick rescaling  n_λ(r) = n(λr):

        E(λ) = λ^{-1} I_grad  +  λ^{-3} I_pot
    """
    lam = np.asarray(lam, dtype=float)
    return I_grad / lam + I_pot / lam**3


def derrick_stationary(I_grad: float, I_pot: float) -> float | None:
    """
    Find the scaling λ* that satisfies dE/dλ = 0:

        λ* = sqrt(−3 I_pot / I_grad)   (requires I_pot < 0 and I_grad > 0)

    Returns None if no real stationary point exists.
    """
    if I_grad <= 0:
        return None
    ratio = -3.0 * I_pot / I_grad
    if ratio <= 0.0:
        return None
    return float(np.sqrt(ratio))


def derrick_d2E(lam: float, I_grad: float, I_pot: float) -> float:
    """
    Second derivative d²E/dλ² at λ:

        d²E/dλ² = 2 λ^{-3} I_grad + 12 λ^{-5} I_pot
    """
    return 2.0 * I_grad / lam**3 + 12.0 * I_pot / lam**5


# ---------------------------------------------------------------------------
# Trial soliton families
# ---------------------------------------------------------------------------

def gaussian_soliton(r: float, A: float = 0.3, width: float = 1.0) -> float:
    """n(r) = 1 + A·exp(−r²/width²)  — compression if A > 0, rarefaction if A < 0."""
    return 1.0 + A * np.exp(-r**2 / width**2)


def yukawa_soliton(r: float, A: float = 0.3, width: float = 1.0) -> float:
    """n(r) = 1 + A·exp(−r/width)/r  — Yukawa-type profile (singular at origin, use r>0)."""
    r = max(r, 1e-6)
    return 1.0 + A * np.exp(-r / width) / r


def lorentzian_soliton(r: float, A: float = 0.3, width: float = 1.0) -> float:
    """n(r) = 1 + A·width²/(r² + width²)."""
    return 1.0 + A * width**2 / (r**2 + width**2)


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyse_soliton(
    label: str,
    A: float,
    width: float,
    profile_fn=gaussian_soliton,
) -> dict:
    """Run the full Derrick analysis for one trial soliton."""

    def n_func(r):
        return profile_fn(r, A=A, width=width)

    I_grad, I_pot, E = energy_integrals(n_func)
    lam_star = derrick_stationary(I_grad, I_pot)

    result = {
        "label": label,
        "A": A,
        "width": width,
        "I_grad": I_grad,
        "I_pot": I_pot,
        "E": E,
        "derrick_sum": I_grad + 3.0 * I_pot,   # = 0 at stationary point
        "lam_star": lam_star,
        "d2E_at_lam_star": (
            derrick_d2E(lam_star, I_grad, I_pot)
            if lam_star is not None else None
        ),
        "type": "compression" if A > 0 else "rarefaction",
        "stationary_point_exists": lam_star is not None,
    }
    return result


def print_report(results: list[dict]) -> None:
    print()
    print("=" * 70)
    print("  Derrick Scaling Analysis — PM Soliton Stability")
    print("=" * 70)
    print()
    for r in results:
        print(f"  [{r['label']}]  type={r['type']},  A={r['A']:+.2f},  width={r['width']:.1f}")
        print(f"    I_grad = {r['I_grad']:+.6f}")
        print(f"    I_pot  = {r['I_pot']:+.6f}")
        print(f"    E      = {r['E']:+.6f}")
        print(f"    I_grad + 3·I_pot = {r['derrick_sum']:+.6f}  (0 → stationary point)")
        if r["lam_star"] is not None:
            ls = r["lam_star"]
            d2 = r["d2E_at_lam_star"]
            stable = "STABLE" if d2 > 0 else "UNSTABLE (saddle)"
            print(f"    λ* = {ls:.4f}  →  d²E/dλ² = {d2:+.6f}  → {stable}")
        else:
            print("    No real stationary point  →  soliton spreads to infinity")
        print()

    print("-" * 70)
    print("  CONCLUSION")
    print("-" * 70)
    print()
    n_comp   = sum(1 for r in results if r["type"] == "compression")
    n_rare   = sum(1 for r in results if r["type"] == "rarefaction")
    n_spoint = sum(1 for r in results if r["stationary_point_exists"])
    n_stable = sum(
        1 for r in results
        if r["d2E_at_lam_star"] is not None and r["d2E_at_lam_star"] > 0
    )

    print(f"  Compression solitons: {n_comp} tested.  "
          f"{'None have' if n_comp else '—'} a Derrick stationary point.")
    print(f"  Rarefaction solitons: {n_rare} tested.")
    print(f"  Stationary points found: {n_spoint}")
    print(f"  Classically stable solitons: {n_stable}")
    print()
    print("  Derrick's theorem for PM's α=2 action in D=3:")
    print()
    print("    E(λ) = λ⁻¹ I_grad + λ⁻³ I_pot")
    print()
    print("    Stationary point requires:  I_grad + 3·I_pot = 0")
    print("                                ⟹ I_pot = −I_grad/3  <  0")
    print("                                ⟹ n < 1 somewhere (rarefaction)")
    print()
    print("    Stability check at any stationary point:")
    print("    d²E/dλ²|_{λ*} = 2·I_grad + 12·I_pot")
    print("                  = 2·I_grad − 4·I_grad  =  −2·I_grad  <  0")
    print()
    print("  ⟹  All stationary points are SADDLES in the scaling direction.")
    print("  ⟹  PM's classical field equations admit NO stable static solitons.")
    print()
    print("  Physical interpretation:")
    print("    ρ_nuc cannot be derived from the classical PM action alone.")
    print("    A quantum length scale (nucleon size) must be imported from")
    print("    outside the theory, or a topological/multi-field mechanism added.")
    print()


def main():
    profiles = [
        # --- Compression solitons (A > 0)  →  I_pot > 0, no stationary point ---
        ("Gaussian compression (A=+0.3, w=1.0)", +0.3, 1.0, gaussian_soliton),
        ("Gaussian compression (A=+0.5, w=1.0)", +0.5, 1.0, gaussian_soliton),
        ("Lorentzian compression (A=+0.3, w=1.0)", +0.3, 1.0, lorentzian_soliton),

        # --- Rarefaction solitons (A < 0)  →  I_pot < 0, stationary point exists but unstable ---
        ("Gaussian rarefaction (A=−0.3, w=1.0)", -0.3, 1.0, gaussian_soliton),
        ("Gaussian rarefaction (A=−0.5, w=1.0)", -0.5, 1.0, gaussian_soliton),
        ("Lorentzian rarefaction (A=−0.3, w=1.0)", -0.3, 1.0, lorentzian_soliton),
        # Tune A so I_grad + 3 I_pot ≈ 0 (at Derrick stationary point)
        ("Gaussian rarefaction (A=−0.3, w=3.0)", -0.3, 3.0, gaussian_soliton),
    ]

    results = [analyse_soliton(lbl, A, w, fn) for lbl, A, w, fn in profiles]
    print_report(results)

    # ----- Analytical verification of the key inequality -----
    print("  Analytical check: V̂_vac sign vs n")
    print()
    for n_val, label in [(0.3, "n=0.3"), (0.7, "n=0.7"), (1.0, "n=1.0"),
                          (1.5, "n=1.5"), (2.0, "n=2.0")]:
        v  = float(v_vac(n_val))
        vp = float(v_vac_prime(n_val))
        sign = "< 1 (neg V̂_vac)" if n_val < 1 else ("> 1 (pos V̂_vac)" if n_val > 1 else "= 1 (zero)")
        print(f"    n = {n_val:.1f}  ({sign:22s}): V̂_vac = {v:+.6f},  V̂_vac' = {vp:+.6f}")
    print()


if __name__ == "__main__":
    main()
