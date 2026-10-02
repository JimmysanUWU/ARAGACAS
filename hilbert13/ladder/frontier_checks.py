"""Exact checks behind 6_SIDE_RESULTS.md sections 6.2, 6.6, 6.7 and 1_CURVE.md Prop. 1.8.

D1  structure of the Q2 group C_{A7}((16)(23))
D4  the involution fixed-point module: isotypic components of the Klein difference and of the (*) lift,
    the Q[G]-modules they generate, and the low-genus quotients carrying the 21- and 35-parts
D5  A7-realisable triangle signatures of genus <= 529 and the pencil-geometry table for gon >= 25

About 5 minutes.
"""
import numpy as np
from math import comb
from fractions import Fraction as F
from itertools import combinations_with_replacement as cwr
from sympy.combinatorics import Permutation, PermutationGroup
from a7 import *
from chartab import TABLE, idx
from cusp import names, a, fix_divisor, PERMS, INV_CLASS, NPTS

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
