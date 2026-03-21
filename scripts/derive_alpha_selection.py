"""Symbolic derivation: selection of α=2 in PM's action measure.

Derives via SymPy:

  1. E-L equation for Sα = ∫ ½|∇φ|² n^α d³x
         →  ∇²φ + (α/2)|∇φ|² = 0  (1-D surrogate; 3D conclusion stated)

  2. Linearisation:  w = n^{α/2}  satisfies ∇²w = 0 when EL = 0
         →  exact, no perturbative approximation

  3. Vacuum solution and PPN β via Taylor expansion of φα(U):
         φα = (2/α) ln(1 + αU),  U = GM/c²r
         →  φα = 2U − αU² + O(U³)
         →  β(α) = α/2

  4. Perihelion factor table for γ=1 (all PM variants share γ=1):
         (2 + 2γ − β)/3 = (4 − α/2)/3
         →  = 1  if and only if  α = 2

  5. Elastic free-energy identity:
         n = exp(φ)  →  ∇n = n ∇φ  →  |∇n|² = n²|∇φ|²
         →  the Hookean strain-energy ½∫|∇n|² d³x is *exactly* S_{α=2}

  6. Uniqueness summary: three independent arguments all select α=2.

Running this script requires only SymPy (included in requirements.txt).

Usage:
    .venv/bin/python scripts/derive_alpha_selection.py
"""

from sympy import (
    symbols, Function, exp, log, series,
    diff, factor, expand, simplify, solve,
    Integer, Rational, Eq,
)

# ── Symbols ───────────────────────────────────────────────────────────────────
x_s     = symbols('x')
U_s     = symbols('U', positive=True)       # Newtonian potential: GM/(c²r)
alp     = symbols('alpha', positive=True)
phi_fn  = Function('phi')

MERCURY_GR = 42.9799    # GR perihelion precession, arcsec/cy

print("=" * 72)
print("Symbolic derivation: action-measure selection in PM (the α-family)")
print("=" * 72)


# ═══════════════════════════════════════════════════════════════════════════════
# 1.  Euler–Lagrange equation for  Sα = ∫ ½(φ')² exp(αφ) dx
# ═══════════════════════════════════════════════════════════════════════════════

print("\n── 1. Euler–Lagrange equation ──────────────────────────────────────────")
print("   Lagrangian density:  L = ½ (φ')² exp(αφ)")

phi_x  = phi_fn(x_s)
phi1   = phi_x.diff(x_s)
phi2   = phi_x.diff(x_s, 2)

L          = Rational(1, 2) * phi1**2 * exp(alp * phi_x)
dL_dphi    = diff(L, phi_x)
dL_dphi1   = diff(L, phi1)
EL_raw     = expand(diff(dL_dphi1, x_s) - dL_dphi)
EL_norm    = simplify(EL_raw / exp(alp * phi_x))

print(f"\n   EL / exp(αφ)  =  {EL_norm}")
print()
print("   ∴  ∇²φ + (α/2)|∇φ|² = 0   [1-D surrogate → 3D spherical result]")


# ═══════════════════════════════════════════════════════════════════════════════
# 2.  Linearisation:  w = n^{α/2}  satisfies ∇²w = 0 exactly when EL = 0
# ═══════════════════════════════════════════════════════════════════════════════

print("\n── 2. Linearisation: ∇²(n^{α/2}) = 0  when  ∇²φ + (α/2)|∇φ|² = 0 ───")

phi_sym  = symbols('phi', real=True)
phi_d1   = symbols("phi_prime",  real=True)   # φ'
phi_d2   = symbols("phi_dprime", real=True)   # φ''

w_expr = exp(alp * phi_sym / 2)                 # w = n^{α/2} = exp(αφ/2)
# d²w/dx² in terms of φ' and φ'' via chain rule:
d2w_chain = (diff(w_expr, phi_sym, 2) * phi_d1**2
             + diff(w_expr, phi_sym)  * phi_d2)
d2w_factored = factor(d2w_chain)

EL_factor = phi_d2 + alp * phi_d1**2 / 2      # EL = 0  ↔  φ'' = -(α/2)(φ')²
d2w_when_EL = d2w_factored.subs(phi_d2, -alp * phi_d1**2 / 2)
d2w_when_EL = simplify(d2w_when_EL)

print(f"\n   d²w/dx² = {simplify(d2w_chain / (alp * w_expr / 2))}"
      f" × (α/2) w")
print(f"           = (α/2) w × [φ'' + (α/2)(φ')²]")
print(f"\n   When EL = 0:")
print(f"   d²w/dx² = {d2w_when_EL}   ✓")
print()
print("   ∴  w = n^{α/2} is the exact linear variable;  ∇²w = 0 in vacuum.")
print("   Point-mass solution:  w(r) = 1 + αGM/(c²r)  for all α.")


# ═══════════════════════════════════════════════════════════════════════════════
# 3.  PPN β from Taylor expansion of vacuum solution
# ═══════════════════════════════════════════════════════════════════════════════

print("\n── 3. PPN β from Taylor expansion of φα(U) ────────────────────────────")
print("   φα(r) = (2/α) ln(1 + αU),   U = GM/c²r  (Newtonian potential)")

phi_vac  = (Integer(2) / alp) * log(1 + alp * U_s)
phi_ser  = series(phi_vac, U_s, 0, 4)

c1 = phi_ser.coeff(U_s, 1)
c2 = phi_ser.coeff(U_s, 2)
c3 = phi_ser.coeff(U_s, 3)

beta_sym = simplify(-c2 / 2)      # φ = 2U - 2β U² + ...  →  β = -c2/2

print(f"\n   Taylor series:  φα = {phi_ser}")
print()
print(f"   O(U¹) coefficient = {c1}             [sets Newton's constant — α-independent]")
print(f"   O(U²) coefficient = {c2}          [= −2β in PPN form φ = 2U − 2β U² + O(U³)]")
print(f"   O(U³) coefficient = {c3}")
print()
print(f"   ∴  β(α)  =  α/2          [verified: {beta_sym}]")


# ═══════════════════════════════════════════════════════════════════════════════
# 4.  Perihelion factor table
# ═══════════════════════════════════════════════════════════════════════════════

print("\n── 4. Perihelion factor table  (γ = 1 for all α: g_ij = n² δ_ij) ─────")
print(f"   Precession per orbit ∝ (2 + 2γ − β)/3  [Will-Nordtvedt formula]")
print()
print(f"   {'α':>4}   {'β':>5}   {'(4−β)/3':>10}   {'Mercury ''/cy':>13}   Notes")
print("   " + "─" * 68)

unique_alpha = solve(Eq(alp / 2, 1), alp)[0]

for a_val in range(5):
    b_val    = Rational(a_val, 2)
    fac      = (4 - b_val) / 3                 # γ=1 → (2+2-β)/3 = (4-β)/3
    mercury  = float(fac) * MERCURY_GR
    sol_form = {
        0: "φ = 2GM/c²r  (linear, no self-energy)",
        1: "φ = 2 ln(1 + GM/c²r)",
        2: "φ = ln(1 + 2GM/c²r)     ← β=1, GR match",
        3: "φ = (2/3) ln(1 + 3GM/c²r)",
        4: "φ = (1/2) ln(1 + 4GM/c²r)",
    }
    marker = "  ◄ UNIQUE GR MATCH" if a_val == 2 else ""
    print(f"   {a_val:>4}   {str(b_val):>5}   {str(fac):>10}   "
          f"{mercury:8.3f} ''/cy{marker}   {sol_form[a_val]}")

print()
print(f"   β = 1  ⟺  α = {unique_alpha}  (unique solution in the integer + half-integer α-family)")
print(f"   GR observed value: 42.980 ± 0.001 arcsec/cy")


# ═══════════════════════════════════════════════════════════════════════════════
# 5.  Elastic free-energy:  ½|∇n|² d³x  =  ½ n²|∇φ|² d³x  =  S_{α=2}
# ═══════════════════════════════════════════════════════════════════════════════

print("\n── 5. Elastic free-energy identity ─────────────────────────────────────")
print("   For an elastic medium with constant bulk modulus K,")
print("   the Hookean strain energy in the compression field n(x) is:")
print()
print("   F_elastic = ½ K ∫ |∇n|² d³x")
print()

phi_v    = symbols('phi', real=True)
nphi_v   = symbols('nabla_phi', real=True)    # |∇φ| placeholder

n_expr      = exp(phi_v)
grad_n_expr = diff(n_expr, phi_v) * nphi_v    # ∇n = (dn/dφ) ∇φ = n ∇φ
elastic_n   = Rational(1, 2) * grad_n_expr**2
elastic_phi = simplify(elastic_n)              # = ½ n² |∇φ|²

print(f"   Chain rule: n = exp(φ)  →  ∇n = exp(φ) ∇φ = n ∇φ")
print(f"   ½ |∇n|²  =  ½ (n ∇φ)²  =  {elastic_phi}")
print()
print("   ∴  ½ K ∫ |∇n|² d³x  =  ½ K ∫ n² |∇φ|² d³x  =  K · S_{α=2}")
print()
print("   The factor n² is the Jacobian of the change of variables n → φ = ln n.")
print("   It is not postulated — it is derived from the chain rule.")
print()

# What α would you get from |∇n^p|² for general power p?
p   = symbols('p', positive=True)
np_expr      = exp(p * phi_v)
grad_np_expr = diff(np_expr, phi_v) * nphi_v
elastic_np   = simplify(Rational(1, 2) * grad_np_expr**2)

print(f"   General power identity: ½|∇(n^p)|² = {elastic_np}")
print(f"   → measure exponent = 2p.  Choosing the compression ratio itself")
print(f"     (p=1, i.e. n directly) gives measure exponent 2 → α = 2.  ✓")


# ═══════════════════════════════════════════════════════════════════════════════
# 6.  Additional distinguishing property: linear variable = n itself (α=2 only)
# ═══════════════════════════════════════════════════════════════════════════════

print("\n── 6. Linearisation variable is n itself only for α=2 ──────────────────")
print()
print(f"   General:  w = n^(α/2).  Values:")
for a_val in range(5):
    b = Rational(a_val, 2)
    if b == 0:
        w_name = "n⁰ = 1  (constant — trivial)"
    elif b == Rational(1, 2):
        w_name = "n^(1/2) = √n"
    elif b == 1:
        w_name = "n¹ = n              ← linear variable IS the density field"
    elif b == Rational(3, 2):
        w_name = "n^(3/2)"
    else:
        w_name = f"n^({b})"
    marker = "  ◄ ONLY α=2" if a_val == 2 else ""
    print(f"   α={a_val}: w = {w_name}{marker}")

print()
print("   For α=2 the linearisation variable w = n is the directly-measurable")
print("   compression field. PM's stellar surface condition is n = 1 (P = 0),")
print("   and n is observable via ρ = ρ_nuc n. No power-law auxiliary needed.")


# ═══════════════════════════════════════════════════════════════════════════════
# 7.  Summary
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("DERIVED RESULT: α=2 is uniquely selected by three independent arguments")
print("=" * 72)
print()
print("  (A) PPN consistency with Mercury perihelion (GR-observed β=1):")
print("      β(α) = α/2  →  β=1  ⟺  α=2.  Every other integer α predicts")
print("      a wrong perihelion rate that is observationally ruled out.")
print()
print("  (B) Elastic free-energy in natural n-coordinates (Hookean medium):")
print("      ½∫|∇n|² d³x  =  ½∫ n²|∇φ|² d³x  =  S_{α=2}.")
print("      The n²-weight is the Jacobian of φ = ln n, not an assumption.")
print("      Any other α corresponds to the elastic energy of n^{α/2}, not n—")
print("      a physically less natural field for PM's compression ratio.")
print()
print("  (C) Linear variable coincides with n itself:")
print("      w_{α=2} = n  (uniquely among the α-family).")
print("      PM's physical observables (density ρ = ρ_nuc n, pressure P ∝ n−1,")
print("      surface condition n=1) are all stated in n, not in an algebraic")
print("      power of n.  The α=2 action is the unique one whose EulerLagrange")
print("      equation is linear directly in the observable field.")
print()
print("  All three constraints select the SAME member of the α-family.  This")
print("  is evidence that α=2 is a structural feature of PM, not a tuning choice.")
print()
print("  Falsification note:")
print("  The n-field (α=2) M_max ≈ 9.5 M☉ < GW150914 component masses (~30 M☉).")
print("  This is a hard prediction: either PM has a distinct energy-phase EOS")
print("  for objects above M_max, or the GW events falsify the α=2 branch.")
print("=" * 72)
