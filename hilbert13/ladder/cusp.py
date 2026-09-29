"""Isotypic projections of divisors supported on the branch fibres of C -> C/A7.

Points over branch point i = left cosets g<x_i>. A divisor is a vector on these cosets.
If a degree-0 divisor has zero projection onto every isotypic component occurring in H^1(C)
(10, 10b, 15, 21, 35), its class in Jac(C) is torsion (Abel-Jacobi is A7-equivariant and the
isotypic parts of Jac(C) for 1, 6, 14a, 14b vanish).
"""
import numpy as np
from a7 import *
from chartab import TABLE, idx
from triples_data import triples
names = ["1", "6", "10", "10b", "14a", "14b", "15", "21", "35"]
a, b, c = triples[0]

def coset_space(x):
    X = [power(x, k) for k in range(order(x))]
    reps, where = [], {}
    for g in A7:
        if g in where: continue
        cos = [mul(g, h) for h in X]
        for y in cos: where[y] = len(reps)
        reps.append(g)
    return reps, where

SP = [coset_space(x) for x in (a, b, c)]
offs = np.cumsum([0] + [len(r) for r, _ in SP])
NPTS = offs[-1]

def act_perm(g):
    """permutation of all branch points induced by left mult by g"""
    P = np.empty(NPTS, dtype=np.int64)
    for i, (reps, where) in enumerate(SP):
        for j, r in enumerate(reps):
            P[offs[i] + j] = offs[i] + where[mul(g, r)]
    return P

PERMS = [act_perm(g) for g in A7]
INV_CLASS = [idx[inv(g)] for g in A7]

def fix_divisor(g):
    P = act_perm(g)
    v = np.zeros(NPTS)
    v[P == np.arange(NPTS)] = 1
    return v

def isotypic_norms(v):
    out = {}
    for n, chi in zip(names, TABLE):
        d = chi[0].real
        w = np.zeros(NPTS, dtype=complex)
        for g, P, ic in zip(A7, PERMS, INV_CLASS):
            # (g . v)[P[i]] = v[i]
            gv = np.empty(NPTS); gv[P] = v
            w += chi[ic] * gv
        w *= d / 2520
        out[n] = np.linalg.norm(w)
    return out

if __name__ == "__main__":
    tau = a
    mu = cyc((0, 2), (1, 3))          # (13)(24): V4 of type I with tau
    nu = cyc((0, 1), (4, 5))          # (12)(56): V4 of type II with tau
    for lab, m in [("typeI mu=(13)(24)", mu), ("typeII mu=(12)(56)", nu)]:
        tm = mul(tau, m)
        dlt = fix_divisor(m) - fix_divisor(tm)
        print(lab, " deg", dlt.sum(), {k: round(v, 4) for k, v in isotypic_norms(dlt).items()})
    dlt = fix_divisor(tau) - fix_divisor(mu)
    print("Fix(tau)-Fix(mu)", {k: round(v, 4) for k, v in isotypic_norms(dlt).items()})
