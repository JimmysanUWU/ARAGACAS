"""Chapter 6 (and Prop. 1.6): the finite inputs of the side results (output: side_checks_output.txt, ~30 s).

R   ramification transport (§6.1): |C(tau)|, the normalisers N(V0), N(V1) and their amalgam; signature table
D1  the Q2 group C_{A7}((16)(23)) (§6.4)
D4  the involution fixed-point module: isotypic components of the Klein difference and of the lift of (*),
    the Q[G]-modules they generate, and the low-genus quotients carrying the 21- and 35-parts (Prop. 1.6, §6.5)
D5  A7-realisable triangle signatures of genus <= 529 and the pencil-geometry table for gon >= 25 (Ch. 4)
D6  GPT's septic Belyi map C/A6 -> C/A7 with passport [2^2 1^3, 4 2 1, 7] (§6.5)
"""
import numpy as np
from math import comb
from fractions import Fraction as F
from itertools import combinations_with_replacement as cwr
from sympy.combinatorics import Permutation, PermutationGroup
from a7 import *
from chartab import TABLE, idx
from curve_checks import names, a, fix_divisor, PERMS, INV_CLASS, NPTS

G = list(A7)
Gidx = {g: i for i, g in enumerate(G)}

def fixdims(K):
    return {n: round((sum(TABLE[k][idx[x]] for x in K) / len(K)).real) for k, n in enumerate(names)}

def proj(v, n):
    chi = TABLE[names.index(n)]
    w = np.zeros(NPTS, dtype=complex)
    for P, ic in zip(PERMS, INV_CLASS):
        gv = np.empty(NPTS); gv[P] = v; w += chi[ic] * gv
    return (w * chi[0].real / 2520).real

def translates(v):
    M = np.empty((len(PERMS), NPTS))
    for i, P in enumerate(PERMS):
        M[i, P] = v
    return M


c = lambda *cy: cyc(*[tuple(x - 1 for x in t) for t in cy])


def check_transport():
    tau = c((1, 2), (3, 4))
    V0 = [E, tau, c((1, 3), (2, 4)), c((1, 4), (2, 3))]
    V1 = [E, tau, c((1, 2), (5, 6)), c((3, 4), (5, 6))]
    def normalizer(V):
        S = set(V); return [g for g in A7 if {conj(g, v) for v in V} == S]
    N0, N1 = normalizer(V0), normalizer(V1)
    print("|C(tau)| =", len(centralizer(tau)), " |N(V0)| =", len(N0), " |N(V1)| =", len(N1),
          " |<N(V0),N(V1)>| =", len(generated(N0 + N1)))
    z, gam, sig = c((1, 3), (2, 4)), c((5, 6, 7)), c((3, 4), (5, 7))
    print("<z, sigma, gamma sigma, tau> =", len(generated([z, sig, mul(gam, sig), tau])), "= |C(tau)|; tau*x involution for all three:",
          all(cycle_type(mul(tau, x)) == (2, 2, 1, 1, 1) for x in (z, sig, mul(gam, sig))))


def table1():
    # generating triples (x, y, (xy)^-1) of A7 with orders (a, b, c), counted up to conjugation of x
    reps = {}
    for g in A7:
        o = order(g)
        reps.setdefault((o, cycle_type(g)), g)
    for sig in [(2, 4, 7), (3, 3, 5), (2, 5, 7), (3, 3, 6), (3, 4, 4), (2, 6, 7), (3, 3, 7), (2, 7, 7), (3, 4, 5)]:
        a, b, cc = sig
        found = 0
        for (o, ct), x in reps.items():
            if o != a:
                continue
            for y in A7:
                if order(y) == b and order(mul(x, y)) == cc and len(generated([x, y])) == 2520:
                    found += 1
        g = 1 + 1260 * (-2 + sum(1 - 1 / m for m in sig))
        print(f"signature {sig}: genus {g:.0f}, generating pairs with x a fixed class rep: {found}")


def a5_example():
    A5 = [p for p in A7 if p[5] == 5 and p[6] == 6]
    t1, t2, t3, t4 = c((1, 2), (3, 4)), c((1, 3), (2, 4)), c((1, 5), (2, 3)), c((1, 4, 5))
    prod = mul(mul(mul(t1, t2), t3), t4)
    print("A5 example: product", "= 1" if prod == E else "!= 1 (other order?)", " generated order",
          len(generated([t1, t2, t3, t4])))
    for perm in [(t4, t3, t2, t1)]:
        p = E
        for x in perm: p = mul(p, x)
        print("   reversed product = 1:", p == E)



print('R   ramification transport (§6.1)')
check_transport()
a5_example()
table1()

# ---- D1 ------------------------------------------------------------------------------------------
Q2 = generated([(4, 3, 1, 6, 0, 5, 2), (0, 6, 2, 3, 5, 4, 1)])
t2 = cyc((1, 6), (2, 3))
orders = {}
for g in Q2: orders[order(g)] = orders.get(order(g), 0) + 1
centre = [g for g in Q2 if all(mul(g, h) == mul(h, g) for h in Q2)]
print("D1  |Q2| =", len(Q2), " equals C(16)(23):", set(Q2) == set(centralizer(t2)),
      " element orders", dict(sorted(orders.items())), " |centre| =", len(centre))
print("    two elements of order 3 -> normal C3; 9 involutions + 6 of order 4 with a D8 Sylow -> C3 : D8")

# ---- D4 ------------------------------------------------------------------------------------------
tau = a                                                    # (01)(23) in letters 0..6
mu = cyc((0, 2), (1, 3)); tmu = mul(tau, mu)
V1 = [cyc((0, 1), (4, 5)), cyc((2, 3), (4, 5))]
D = fix_divisor(mu) + fix_divisor(tmu) - sum(fix_divisor(x) for x in V1)
goods = [cyc((0, 1), (i, j)) for (i, j) in [(4, 5), (4, 6), (5, 6)]]
goods += [mul(tau, g) for g in goods]
v = sum(fix_divisor(g) for g in goods) - 3 * (fix_divisor(mu) + fix_divisor(tmu))
print("\nD4  |D_delta|^2 by isotype:", {n: round(float(proj(D, n) @ proj(D, n)), 3) for n in names
                                        if np.linalg.norm(proj(D, n)) > 1e-9})
print("    C(tau)-fixed dims = multiplicities in Z[A7/C(tau)]:", fixdims(centralizer(tau)))
for n in ("21", "35"):
    M, Mv = translates(proj(D, n)), translates(proj(v, n))
    r = np.linalg.matrix_rank(M, tol=1e-8); rv = np.linalg.matrix_rank(Mv, tol=1e-8)
    rj = np.linalg.matrix_rank(np.vstack([M, Mv]), tol=1e-8)
    print(f"    {n}: rank Q[G]e(D) = {r}, rank Q[G]e(v) = {rv}, joint {rj}  (one copy each, the same copy)")
L25 = generated([cyc((1, 2, 3, 4, 5)), cyc((1, 6), (2, 5))])      # PSL(2,5) on P^1(F5) = letters 1..6
Q1 = generated([cyc((0, 1, 2)), cyc((3, 4, 5)), cyc((0, 3, 1, 4), (2, 5))])
for K, n, lab in ((L25, "21", "L2(5)"), (Q1, "35", "3^2:4 (Q1 group)")):
    fd = fixdims(K)
    genus = (3 * fd["10"] + 3 * fd["10b"] + 2 * fd["15"] + 2 * fd["21"] + 4 * fd["35"]) // 2
    x = proj(D, n); seen = {}
    for g in G:
        Kg = frozenset(mul(mul(inv(g), k), g) for k in K)
        if Kg in seen: continue
        avg = np.zeros(NPTS)
        for k in Kg:
            gv = np.empty(NPTS); gv[PERMS[Gidx[k]]] = x; avg += gv
        seen[Kg] = np.linalg.norm(avg) > 1e-9
    print(f"    C/{lab}: |K| = {len(K)}, fixed irreducibles {[m for m in names if fd[m]]}, genus {genus};"
          f" {sum(seen.values())} of {len(seen)} conjugates see the {n}-part of D_delta")

# ---- D5 ------------------------------------------------------------------------------------------
print("\nD5  A7-realisable triangle signatures with genus <= 529 (ordered generating triples / |A7|)")
by_order = {}
for g in G: by_order.setdefault(order(g), []).append(g)
reps = {}
for g in G: reps.setdefault(cycle_type(g), g)
csize = {ct: sum(1 for g in G if cycle_type(g) == ct) for ct in reps}
def ncyc(g): return len(cycle_type(g))
for p, q, r in cwr([2, 3, 4, 5, 6, 7], 3):
    chi = 1 - F(1, p) - F(1, q) - F(1, r)
    if chi <= 0 or (2520 * chi).denominator != 1 or (2520 * chi) % 2: continue
    genus = int(2520 * chi) // 2 + 1
    if genus > 529: continue
    total = 0
    for ct, x in reps.items():
        if order(x) != p: continue
        for y in by_order.get(q, []):
            z = inv(mul(x, y))
            if order(z) != r or 21 - ncyc(x) - ncyc(y) - ncyc(z) < 12: continue   # 7-sheeted RH
            if PermutationGroup([Permutation(list(x)), Permutation(list(y))]).order() == 2520:
                total += csize[ct]
    if total:
        print(f"    g = {genus:3d}  {(p, q, r)}  [{total // 2520}]  lambda1 threshold for gon>=25: {48 / (genus - 1):.4f}")

def pi1(d, r=7):
    m1, e1 = divmod(d - 1, r)
    return comb(m1, 2) * r + m1 * (e1 + 1) + (1 if e1 == r - 1 else 0)

def dep_costs(m):
    """minimal singularity cost of 1, 2, 3 distinct dependent thirds: edges {P,Q} with r_P + r_Q = m,
    each centre of multiplicity r costing C(r,2) once (a pencil of (1,1)-forms is fixed by its two base points)"""
    one = min(comb(r, 2) + comb(m - r, 2) for r in range(1, m))
    vee = min(comb(r, 2) + 2 * comb(m - r, 2) for r in range(1, m))
    star = min(comb(r, 2) + 3 * comb(m - r, 2) for r in range(1, m))
    path = min(2 * comb(r, 2) + 2 * comb(m - r, 2) for r in range(1, m))
    tri = 3 * comb(m // 2, 2) if m % 2 == 0 else 10**9
    return one, min(vee, 2 * one), min(star, tri, path, vee + one, 3 * one)

print("\n    gon >= 25 for g >= 336:  m | pair bound | pi1(3m,7) | budget (m-1)^2-336 | cost of 1/2/3 dependent thirds")
for m in range(20, 25):
    pair = max([e * (m // e - 1) ** 2 + (e - 1) * (m - 1) for e in range(2, m // 2 + 1) if m % e == 0], default=None)
    print(f"      {m} | {pair} | {pi1(3 * m)} | {(m - 1) ** 2 - 336} | {dep_costs(m)}")

# ---- D6 ------------------------------------------------------------------------------------------
import sympy as sp
z, T = sp.symbols('z T')
for sg in (1, -1):                       # s = 1 +- 4 sqrt(21)/21, the roots of 21 s^2 - 42 s + 5
    s_ = 1 + sg * 4 * sp.sqrt(21) / 21; t_ = sp.Rational(2, 3) * s_ - sp.Rational(10, 21)
    q = z**7 / 7 - (s_ + 1) * z**6 / 6 + (s_ + t_) * z**5 / 5 - t_ * z**4 / 4
    rs = -128 * (51 * s_ - 65) / 1750329
    ok = [sp.expand(sp.diff(q, z) - z**3 * (z - 1) * (z**2 - s_ * z + t_)) == 0, sp.simplify(q.subs(z, 1)) == 0,
          all(sp.simplify(c) == 0 for c in sp.Poly(sp.expand(q - rs), z).rem(sp.Poly(z**2 - s_ * z + t_, z)).all_coeffs()),
          sp.simplify(sp.discriminant(sp.Poly(sp.expand(1 - q / rs - T), z)) + 7 * (-1 / rs)**6 * T**2 * (T - 1)**4) == 0]
    print(("\nD6  " if sg == 1 else "    ") + f"s = 1 {'+-'[sg < 0]} 4 sqrt21/21:  q' = z^3 (z-1)(z^2-sz+t), q(1) = 0, "
          f"q = r_s at both quadratic critical points, disc(beta - T) = -7 r_s^-6 T^2 (T-1)^4: {all(ok)}")
pr = 17; s0 = 3; inv_ = lambda x: pow(x, -1, pr); t0 = (2 * s0 * inv_(3) - 10 * inv_(21)) % pr
qp = sp.Poly(z**7 * inv_(7) - (s0 + 1) * inv_(6) * z**6 + (s0 + t0) * inv_(5) * z**5 - t0 * inv_(4) * z**4, z, modulus=pr)
f = sp.Poly(1, z, modulus=pr) - qp * inv_(-128 * (51 * s0 - 65) * inv_(1750329) % pr) - sp.Poly(2, z, modulus=pr)
print("    mod 17 (s = 3), beta = 2 factors in degrees", sorted(g.degree() for g, _ in f.factor_list()[1]),
      "-> an element of order 10, not in L2(7) = its own normaliser in S7; so the monodromy is A7")
