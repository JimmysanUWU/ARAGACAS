"""Exact (finite, rational/algebraic) checks behind 3_HARMONIC_HERSCH.md (Lemma 3.3).  No floating point.

Uses the ATLAS character table of A7 from cover.py (its orthonormality is re-verified there), with cover.py's
names: 14a = 14_(5,2) (the 2-subset module), 14b = 14_(4,3).  Power maps are computed from actual permutations.

  1. 14a is the 2-subset module: chi_14a = pi_2 - pi_1 (fixed 2-subsets minus fixed points), elementwise on A7.
  2. dim (wedge^3 V)^A7 for every real irreducible V, and the decomposition of wedge^2 14a
     (multiplicity-free: 10 + 10b + 15 + 21 + 35).
  3. Induced modules used by the certificate: Ind_{S5}^{A7} 1 = 1 + 6 + 14a (S5 = stabiliser of {5,6}),
     Ind_{A6}^{A7} 1 = 1 + 6, and the multiplicities seen by Q1 = (3^2:4, sgn) and Q2 = (C(16)(23), sgn).
  4. The degree-lattice numbers of section 7: |G| / (exp M(G) lcm(e_i)) for A7 (2,4,7) and L2(13) (2,3,7).
"""
import os
import sys
from itertools import combinations
import sympy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a7 import A7, E, cyc, cycle_type, generated, mul, power, a7_class
from cover import TABLE, CLASSES, SIZES, check_table, character, multiplicities, QUOTIENT_DATA

check_table()
names = list(TABLE)
r7 = cyc(tuple(range(7)))
cls7a = set(a7_class(r7))


def cls(g):
    ct = cycle_type(g)
    if ct == (7,):
        return 7 if g in cls7a else 8
    return CLASSES.index(ct)


reps = {}
for g in A7:
    reps.setdefault(cls(g), g)
assert len(reps) == 9
p2 = [cls(power(reps[i], 2)) for i in range(9)]
p3 = [cls(power(reps[i], 3)) for i in range(9)]


def inner(f, h):
    return sympy.nsimplify(sympy.simplify(sum(SIZES[i] * f[i] * sympy.conjugate(h[i]) for i in range(9)) / 2520))


def decompose(f):
    return {x: inner(f, TABLE[x]) for x in names if inner(f, TABLE[x]) != 0}


# 1. 14a = pi_2 - pi_1, elementwise
pairs = list(combinations(range(7), 2))
for g in A7:
    fix1 = sum(1 for i in range(7) if g[i] == i)
    fix2 = sum(1 for (i, j) in pairs if {g[i], g[j]} == {i, j})
    assert TABLE["14a"][cls(g)] == fix2 - fix1
print("1. chi_14a(g) = #fixed 2-subsets - #fixed points for all 2520 g: 14a = 14_(5,2)")

# 2. wedge^2, wedge^3
real = {"6": TABLE["6"], "10+10b": [a + b for a, b in zip(TABLE["10"], TABLE["10b"])], "14a": TABLE["14a"],
        "14b": TABLE["14b"], "15": TABLE["15"], "21": TABLE["21"], "35": TABLE["35"]}
print("2. real irreducibles:")
for nm, x in real.items():
    w3 = [(x[i] ** 3 - 3 * x[i] * x[p2[i]] + 2 * x[p3[i]]) / 6 for i in range(9)]
    w2 = [(x[i] ** 2 - x[p2[i]]) / 2 for i in range(9)]
    d3 = inner(w3, TABLE["1"])
    print(f"   {nm:7s}: dim (wedge^3)^A7 = {d3};  wedge^2 = {decompose(w2)}")
x = TABLE["14a"]
w2 = [(x[i] ** 2 - x[p2[i]]) / 2 for i in range(9)]
assert decompose(w2) == {"10": 1, "10b": 1, "15": 1, "21": 1, "35": 1}
assert inner(w2, w2) == 5
print("   wedge^2 14_(5,2) = 10 + 10b + 15 + 21 + 35, multiplicity-free (<chi,chi> = 5)")

# 3. induced modules and quotient multiplicities
S5 = character([cyc((0, 1, 2)), cyc((0, 1, 2, 3, 4)), cyc((0, 1), (5, 6))], [1, 1, 1])
A6 = character([cyc((0, 1, 2)), cyc((1, 2, 3, 4, 5))], [1, 1])
assert len(S5) == 120 and len(A6) == 360
assert all(cycle_type(k) != (7,) for k in list(S5) + list(A6))
print("3. Ind_{S5}^{A7} 1 =", {k: v for k, v in multiplicities(S5).items() if v},
      "  Ind_{A6}^{A7} 1 =", {k: v for k, v in multiplicities(A6).items() if v})
assert {k: v for k, v in multiplicities(S5).items() if v} == {"1": 1, "6": 1, "14a": 1}
assert {k: v for k, v in multiplicities(A6).items() if v} == {"1": 1, "6": 1}
for q, (gens, vals) in QUOTIENT_DATA.items():
    print(f"   {q}:", {k: v for k, v in multiplicities(character(gens, vals)).items() if v})

# 4. degree lattice numbers
print("4. |A7|/(exp M(A7) lcm(2,4,7)) =", sympy.Rational(2520, 6 * 28), "; 1296/15 =", sympy.Rational(1296, 15),
      "; |L2(13)|/(exp M lcm(2,3,7)) =", sympy.Rational(1092, 2 * 42), "; 216/13 =", sympy.Rational(216, 13))
print("all exact checks passed")
