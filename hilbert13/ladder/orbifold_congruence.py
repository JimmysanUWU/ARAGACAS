"""Orbifold congruences for symmetric pencils (Chapter 7, Lemma 7.24) [P][X].

A gonal pencil f : C -> P^1 of degree d with class stabiliser K and kernel N descends to Y = C/N (degree k = d/|N|),
equivariantly for Q = K/N acting on P^1 through rho : Q -> PGL_2 = SO(3), and then to
    h : Z = C/K -> P^1/Q   of degree k.
Y -> Z is the normalised pull-back of the Galois cover P^1 -> P^1/Q along h.  If P in Z has h(P) = b and local
degree e, the canonical inertia generator at a point of Y over P is mapped by rho to r_b^e, where r_b is the rotation
by 2 pi / n_b at the corresponding point of P^1.  Hence, for every branch point b of P^1 -> P^1/Q,
    k = sum_{h(P) = b} e_P,   e_P = (residue forced by the inertia class of P) mod n_b,   e_P = 0 mod n_b otherwise.
This script lists the inertia classes of Y -> Z for every stratum (K, N), finds which (b, residue) each class allows
(the rotation axes of rho(q)), and decides by dynamic programming which k are compatible.  Lemma 7.17 (free
divisibility) is the case without inertia.

Points of C with nontrivial stabiliser: over 0, 1, oo of C/A7 they are the cosets g<c> (c of order 2, 4, 7), with
canonical inertia generator g c g^-1; this G-set is independent of the class of the curve.

Usage: python3 orbifold_congruence.py            (prints the inertia data and the strata filter)
"""
import sys
import itertools
from a7 import cyc
from equivariant_spectral import (G, GI, MT, INV, ORD, E0, REPS, closure, normal_subgroups, quotient_orders,
                                  so3_type, describe, genus, strata, ckey, propagate)

C0, C1, CINF = GI[cyc((0, 1), (2, 3))], GI[cyc((0, 1, 2, 3), (4, 5))], GI[cyc((0, 1, 2, 3, 4, 5, 6))]
assert ORD[C0] == 2 and ORD[C1] == 4 and ORD[CINF] == 7


def cosets(c):
    H = closure([c])
    seen, reps = set(), []
    for g in range(len(G)):
        if g in seen:
            continue
        reps.append(g)
        seen.update(MT[g][h] for h in H)
    return reps


# (canonical generator gamma_x) for every point x with nontrivial stabiliser, grouped by branch value of C -> C/A7
POINTS = [MT[MT[g][c]][INV[g]] for c in (C0, C1, CINF) for g in cosets(c)]
assert len(POINTS) == 1260 + 630 + 360


def power(x, a):
    y = E0
    for _ in range(a):
        y = MT[y][x]
    return y


def inertia(K, N):
    """list of (q, multiplicity): canonical inertia generators of Y = C/N -> Z = C/K, one per point of Z, as
    elements of K (representatives of K/N); points of Z are K-orbits of points of C."""
    Ks = set(K)
    out = []
    # K-orbits of points: a point is (coset), represented by its gamma together with a coset label.  Two points with
    # the same gamma can differ (over 0: 12 points per involution), so we orbit the coset set itself.
    for c in (C0, C1, CINF):
        H = closure([c])
        reps = cosets(c)
        label = {}
        for idx, g in enumerate(reps):
            for h in H:
                label[MT[g][h]] = idx
        done = set()
        for idx, g in enumerate(reps):
            if idx in done:
                continue
            orbit = {label[MT[k][g]] for k in K}
            done |= orbit
            gamma = MT[MT[g][c]][INV[g]]
            m = ORD[gamma]
            # K_x = <gamma> cap K, canonical generator gamma^a with a = m / |K_x|
            Kx = [a for a in range(1, m + 1) if m % a == 0 and power(gamma, a) in Ks]
            a = min(Kx)
            q = power(gamma, a)
            if q in N:                                   # trivial inertia in Y -> Z (incl. q = 1)
                continue
            out.append(q)
    return out


def axis_options(K, N, t, qo):
    """for the SO(3)-type t of Q = K/N: the branch points of P^1 -> P^1/Q with their orders, and a function
    opts(q, glob) -> list of (b, residue mod n_b) allowed for a point with inertia q, given a global orientation
    parameter glob (a unit modulo the relevant order).  Also returns the list of global parameters."""
    Kl = sorted(K)
    Nn = set(N)

    def coset_eq(x, y):
        return MT[x][INV[y]] in Nn

    def qpow(x, a):
        return power(x, a)

    def conj_in_Q(x, y):                                 # x ~ y in K/N
        return any(coset_eq(MT[MT[k][x]][INV[k]], y) for k in Kl)

    q_ord = len(K) // len(N)
    if t[0] == "C":
        nn = q_ord
        sigma = next(k for k in Kl if qo[k] == nn)
        def expo(x):                                     # x = sigma^a mod N
            return next(a for a in range(nn) if coset_eq(qpow(sigma, a), x))
        units = [j for j in range(1, nn) if __import__("math").gcd(j, nn) == 1] or [0]   # nn = 1: no inertia
        branch = {"b0": nn, "binf": nn}
        def opts(x, j):
            a = expo(x)
            return [("b0", (a * j) % nn), ("binf", (-a * j) % nn)]
        return branch, opts, units
    if t == "V4":
        invs = [k for k in Kl if qo[k] == 2]
        reps = []
        for x in invs:
            if not any(coset_eq(x, r) for r in reps):
                reps.append(x)
        assert len(reps) == 3
        branch = {f"b{i}": 2 for i in range(3)}
        def opts(x, _):
            i = next(i for i, r in enumerate(reps) if coset_eq(x, r))
            return [(f"b{i}", 1)]
        return branch, opts, [1]
    if t[0] == "D":
        m = int(t[1:])
        rots = [k for k in Kl if qo[k] == m]
        sigma = rots[0]
        R = [qpow(sigma, a) for a in range(m)]
        def rot_expo(x):
            for a in range(m):
                if coset_eq(R[a], x):
                    return a
            return None
        refl = [k for k in Kl if rot_expo(k) is None]
        if m % 2:
            branch = {"pole": m, "b1": 2, "b2": 2}
            def opts(x, j):
                a = rot_expo(x)
                if a is not None:
                    return [("pole", (a * j) % m), ("pole", (-a * j) % m)]
                return [("b1", 1), ("b2", 1)]
        else:
            s0 = refl[0]
            branch = {"pole": m, "bs": 2, "bt": 2}
            def opts(x, j):
                a = rot_expo(x)
                if a is not None:
                    return [("pole", (a * j) % m), ("pole", (-a * j) % m)]
                return [("bs", 1)] if conj_in_Q(x, s0) else [("bt", 1)]
        units = [j for j in range(1, m) if __import__("math").gcd(j, m) == 1]
        return branch, opts, units
    if t == "A4":
        threes = [k for k in Kl if qo[k] == 3]
        s0 = threes[0]
        branch = {"b2": 2, "bv": 3, "bf": 3}
        def opts(x, u):
            o = qo[x]
            if o == 2:
                return [("b2", 1)]
            sgn = 1 if conj_in_Q(x, s0) else -1
            return [("bv", (sgn * u) % 3), ("bf", (-sgn * u) % 3)]
        return branch, opts, [1, 2]
    if t == "S4":
        # double transpositions = squares of elements of order 4
        sq = {MT[k][k] for k in Kl if qo[k] == 4}
        branch = {"b2": 2, "b3": 3, "b4": 4}
        def opts(x, _):
            o = qo[x]
            if o == 2:
                return [("b4", 2)] if any(coset_eq(x, y) for y in sq) else [("b2", 1)]
            if o == 3:
                return [("b3", 1), ("b3", 2)]
            return [("b4", 1), ("b4", 3)]
        return branch, opts, [1]
    if t == "A5":
        fives = [k for k in Kl if qo[k] == 5]
        s0 = fives[0]
        branch = {"b2": 2, "b3": 3, "b5": 5}
        def opts(x, u):
            o = qo[x]
            if o == 2:
                return [("b2", 1)]
            if o == 3:
                return [("b3", 1), ("b3", 2)]
            c = 1 if conj_in_Q(x, s0) else 2           # class of sigma or sigma^2
            return [("b5", (c * u) % 5), ("b5", (-c * u) % 5)]
        return branch, opts, [1, 2]
    raise ValueError(t)


def compatible_k(K, N, t, qo, kmax):
    """set of k <= kmax compatible with the orbifold congruences."""
    pts = inertia(K, N)
    branch, opts, globs = axis_options(K, N, t, qo)
    names = sorted(branch)
    ok = set()
    for glob in globs:
        # group points by option list
        groups = {}
        for x in pts:
            key = tuple(sorted(set(opts(x, glob))))
            groups[key] = groups.get(key, 0) + 1
        # DP over states: per branch point (sum of residues mod n_b, sum of least positive residues capped)
        start = tuple((0, 0) for _ in names)
        states = {start}
        for key, cnt in groups.items():
            new = set()
            for st in states:
                for comp in compositions(cnt, len(key)):
                    s = list(st)
                    for (b, r), c in zip(key, comp):
                        i = names.index(b)
                        nb = branch[b]
                        lpr = r if r else nb
                        s[i] = ((s[i][0] + c * r) % nb, min(s[i][1] + c * lpr, kmax + 1))
                    new.add(tuple(s))
            states = new
        for st in states:
            for k in range(1, kmax + 1):
                if all(k % branch[b] == st[i][0] and st[i][1] <= k for i, b in enumerate(names)):
                    ok.add(k)
    return ok, pts


def compositions(n, parts):
    if parts == 1:
        yield (n,)
        return
    for i in range(n + 1):
        for rest in compositions(n - i, parts - 1):
            yield (i,) + rest


def inertia_summary(K, N, qo, pts):
    hist = {}
    for x in pts:
        hist[qo[x]] = hist.get(qo[x], 0) + 1
    return dict(sorted(hist.items()))


if __name__ == "__main__":
    print("Orbifold congruences (Lemma 7.24) on the strata of Corollary 7.23(2), 25 <= d <= 41.")
    print("inertia of C/N -> C/K: {order in K/N: number of points of C/K}")
    print("  |N| orbits(N)          |K| K/N   g(C/K) inertia                       before              after")
    pairs = {}
    for row in strata(with_groups=True):
        (nn, orbN, gY, lbN, kk, t, d, why, note, K, Nn) = row
        pairs.setdefault((K, Nn), [t, []])[1].append((d, why))
    seen = set()
    for (K, Nn), (t, lst) in pairs.items():
        before = sorted({d for d, w in lst if not w})
        if not before:
            continue
        qo = quotient_orders(K, Nn)
        kmax = 41 // len(Nn)
        ok, pts = compatible_k(K, Nn, t, qo, kmax)
        after = [d for d in before if d // len(Nn) in ok]
        line = (f"  {len(Nn):3d} {str(describe(Nn)[0]):18s} {len(K):4d} {t:4s}  {genus(K):5d}  "
                f"{str(inertia_summary(K, Nn, qo, pts)):29s} {str(before):19s} {after}")
        if line not in seen:
            seen.add(line)
            print(line)
