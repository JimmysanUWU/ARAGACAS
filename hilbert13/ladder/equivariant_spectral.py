"""Equivariant Li-Yau: which symmetric pencils can have low degree (Chapter 7, §7.13) [P][X], numbers [C] from Ch. 2.

Let f : C -> P^1 have degree d, K = stabiliser of the pencil field C(f) in A7, N = kernel of K on the pencil, so that
K/N embeds in PGL_2 and (after conjugation) in SO(3).  Douady-Earle balancing is conformally natural, so the balanced
map x : C -> S^2 is K-equivariant: x o k = R(k) x.  Hence x lies in (L^2(C) (x) R^3)^K, and its component in an
A7-isotype sigma vanishes unless m_sigma(K, N) = <chi_R, sigma|_K> > 0.  Li-Yau then gives 8 pi d >= Lambda Area,
i.e. d >= 67.5 Lambda, where Lambda is the least nonzero eigenvalue whose eigenspace meets an isotype with m > 0.

chi_R(k) depends only on the order of kN in K/N.  For orders 5 and 7 the rotation angle is not determined, but A7
characters are constant on the elements of order 5 and, for every irreducible except 10 and 10b, on the elements of
order 7; and in each cyclic group of rotations of order 5 (resp. 7) the mean of 1 + 2 cos(theta) over the
nonidentity elements is 1/2 (resp. 2/3), also on each of the two A7-classes of 7-cycles.  So the effective values
    ord(kN):  1    2    3    4    5    6    7
    chi_R  :  3   -1    0    1   1/2   2   2/3
give every m_sigma exactly (for 10 + 10b we use alpha + alpha-bar = -1).

The script enumerates all subgroups K of A7 up to conjugacy, all normal N with K/N cyclic, dihedral, A4, S4 or A5,
and prints m_sigma for sigma in 14a, 14b, 21, 6, {10,10b,15,35}, 1, with the free-action divisibility of d.

Certified spectral inputs (Ch. 2-3; classes 0, 1):  spec(C) below 0.55998 is {0, lambda_1}, E_1 = 14a;
isotypes 10, 10b, 15, 35 have eigenvalues >= 0.689027 (Q2), the trivial isotype >= 0.997833 (Q0, nonconstant),
the 6 >= 0.99783 (C/A6).  Classes 12, 14: 10, 10b, 15, 35 >= 0.762988, trivial >= 0.997833, 6 >= 0.99783 (C/A6,
certify_th_output.txt); 14a, 14b, 21 only >= 0.340893 (the 21 has an eigenvalue near 0.38 there, §2.5).

Usage: python3 equivariant_spectral.py
"""
import sys
from fractions import Fraction as Fr
from math import ceil
import itertools
from a7 import A7, E, mul, inv, order, cycle_type
from cover import TABLE, CLASSES

G = list(A7)
GI = {g: i for i, g in enumerate(G)}
n = len(G)
MT = [[GI[mul(x, y)] for y in G] for x in G]
INV = [GI[inv(x)] for x in G]
ORD = [order(x) for x in G]
E0 = GI[E]
FIX = {1: None, 2: 18, 4: 2, 7: 3}                      # fixed points on C by element order; 3, 5, 6 act freely


def col(g):
    ct = cycle_type(g)
    return 7 if ct == (7,) else CLASSES.index(ct)      # 7a/7b merged: only 10, 10b distinguish them


# characters as exact rationals on the merged columns; 10 + 10b combined (alpha + alpha-bar = -1 on 7-cycles)
CHI = {}
for name in ["1", "6", "14a", "14b", "15", "21", "35"]:
    CHI[name] = [Fr(int(TABLE[name][j])) for j in range(8)]
CHI["10+10b"] = [Fr(int(TABLE["10"][j] + TABLE["10b"][j])) for j in range(7)] + [Fr(-1)]
COL = [col(g) for g in G]


def closure(gens):
    S = {E0}
    fr = [E0]
    while fr:
        fr = [y for x in fr for h in gens if (y := MT[x][h]) not in S and not S.add(y)]
    return frozenset(S)


def conj_set(g, S):
    gi = INV[g]
    return frozenset(MT[MT[g][s]][gi] for s in S)


_CK = {}


def ckey(S):
    """canonical representative of the conjugacy class of S (cached)."""
    if S not in _CK:
        if len(S) == n:
            _CK[S] = tuple(range(n))
        else:
            _CK[S] = min(tuple(sorted(conj_set(g, S))) for g in range(n))
    return _CK[S]


# ---- all subgroups up to conjugacy (iterated closure, as in quotient_gonality.py)
GENS = {}


def sub(gens):
    S = closure(gens)
    GENS.setdefault(S, list(gens))
    return S


classes = {}
frontier = []
for i in range(n):
    S = sub([i])
    k = ckey(S)
    if k not in classes:
        classes[k] = S
        frontier.append(S)
while frontier:                                        # extend each class by one element, up to N(S)-conjugacy
    nf = []
    for S in frontier:
        NS = [g for g in range(n) if conj_set(g, S) == S]
        covered = set(S)
        for j in range(n):
            if j in covered:
                continue
            T = sub(GENS[S] + [j])
            k = ckey(T)
            if k not in classes:
                classes[k] = T
                nf.append(T)
            for g in NS:                               # <S, s g j g^-1> = g <S, j> g^-1
                jg = MT[MT[g][j]][INV[g]]
                covered.update(MT[s][jg] for s in S)
    frontier = nf
REPS = sorted(classes.values(), key=lambda S: (len(S), sorted(S)))


def normal_subgroups(K):
    Kl = sorted(K)
    ccl = []                                            # normal closures of single elements
    for x in Kl:
        ccl.append(closure(list({MT[MT[k][x]][INV[k]] for k in Kl})))
    Ns = {frozenset([E0]), K}
    Ns.update(ccl)
    changed = True
    while changed:
        changed = False
        for A, B in itertools.combinations(list(Ns), 2):
            Cc = closure(list(A | B))
            if Cc not in Ns:
                Ns.add(Cc)
                changed = True
    return sorted(Ns, key=len)


def quotient_orders(K, N):
    """order of kN in K/N for every k in K."""
    out = {}
    for k in K:
        e, x = 1, k
        while x not in N:
            x = MT[x][k]
            e += 1
        out[k] = e
    return out


def so3_type(K, N, qo):
    """isomorphism type of K/N if it is a finite subgroup of SO(3), else None."""
    q = len(K) // len(N)
    hist = {}
    for k, e in qo.items():
        hist[e] = hist.get(e, 0) + 1
    hist = {e: c // len(N) for e, c in hist.items()}  # element-order histogram of K/N
    if q in hist:
        return f"C{q}"
    if q % 2 == 0 and q >= 4:
        m = q // 2                                     # dihedral D_m: index-2 cyclic group, all other elements involutions
        if (m == 2 or hist.get(m, 0) > 0) and hist.get(2, 0) == m + (1 if m % 2 == 0 else 0):
            return "V4" if m == 2 else f"D{m}"
    if q == 12 and hist == {1: 1, 2: 3, 3: 8}:
        return "A4"
    if q == 24 and hist == {1: 1, 2: 9, 3: 8, 4: 6}:
        return "S4"
    if q == 60 and hist == {1: 1, 2: 15, 3: 20, 5: 24}:
        return "A5"
    return None


CHIR = {1: Fr(3), 2: Fr(-1), 3: Fr(0), 4: Fr(1), 5: Fr(1, 2), 6: Fr(2), 7: Fr(2, 3)}


def mults(K, qo):
    out = {}
    for name, chi in CHI.items():
        s = sum(CHIR[qo[k]] * chi[COL[k]] for k in K) / len(K)
        assert s.denominator == 1 and s >= 0, (name, s)
        out[name] = int(s)
    return out


def describe(S):
    par = list(range(7))

    def f(x):
        while par[x] != x:
            x = par[x]
        return x
    for h in S:
        for x in range(7):
            par[f(x)] = f(G[h][x])
    orb = sorted(c for c in [sum(1 for x in range(7) if f(x) == r) for r in range(7)] if c)
    hist = {}
    for s in S:
        if ORD[s] > 1:
            hist[ORD[s]] = hist.get(ORD[s], 0) + 1
    return orb, dict(sorted(hist.items()))


def bounds(m, cls01=True):
    """lower bound on d from the certified spectral data."""
    if cls01:
        lam = []
        if m["14a"]:
            lam.append(Fr(340893, 10**6))            # lambda_1 itself (Ch. 2); Ch. 3 gives 25 anyway
        if m["14b"] or m["21"]:
            lam.append(Fr(559982, 10**6))
        if m["10+10b"] or m["15"] or m["35"]:
            lam.append(Fr(689027, 10**6))
        if m["6"] or m["1"]:
            lam.append(Fr(997830, 10**6))
    else:
        lam = []
        if m["14a"] or m["14b"] or m["21"]:
            lam.append(Fr(340893, 10**6))
        if m["10+10b"] or m["15"] or m["35"]:
            lam.append(Fr(762988, 10**6))
        if m["6"] or m["1"]:
            lam.append(Fr(997830, 10**6))
    L = min(lam) if lam else None
    return L, (ceil(Fr(135, 2) * L) if L is not None else None)


if __name__ == "__main__":
    print(f"{len(REPS)} conjugacy classes of subgroups of A7")
    assert sum(1 for _ in REPS) == len(REPS)
    rows = []
    for K in REPS:
        orbK, histK = describe(K)
        for N in normal_subgroups(K):
            qo = quotient_orders(K, N)
            t = so3_type(K, N, qo)
            if t is None:
                continue
            m = mults(K, qo)
            # free-action divisibility: 3 | d if K has elements of order 3 or 6; 9 | d if 9 | |K|; 5 | d if 5 | |K|
            div = 1
            if any(ORD[k] in (3, 6) for k in K):
                div *= 9 if len(K) % 9 == 0 else 3
            if len(K) % 5 == 0:
                div *= 5
            L01, d01 = bounds(m, True)
            L12, d12 = bounds(m, False)
            rows.append((len(K), len(N), t, orbK, histK, m, div, d01, d12))
    print("  |K|  |N|  K/N   orbits(K)        m: 14a 14b 21  6 10s 15 35  1   d | div   d>= (cls 0,1)  d>= (cls 12,14)")
    for (k, nn, t, orbK, histK, m, div, d01, d12) in rows:
        print(f"  {k:4d} {nn:4d}  {t:4s}  {str(orbK):16s}   {m['14a']:3d} {m['14b']:3d} {m['21']:2d} {m['6']:2d} {m['10+10b']:3d}"
              f" {m['15']:2d} {m['35']:2d} {m['1']:2d}   {div:5d}   {d01!s:>10s}   {d12!s:>10s}")


# ------------------------------------------------------------------------------------ strata of Corollary 7.14
FIXN = [{2: 18, 4: 2, 7: 3}.get(o, 0) for o in ORD]


def genus(S):
    return (270 - sum(FIXN[s] for s in S if s != E0)) // (2 * len(S)) + 1


def normaliser(S):
    return frozenset(g for g in range(n) if conj_set(g, S) == S)


def divisors(m):
    return [k for k in range(1, m + 1) if m % k == 0]


def propagate(base=25, maxS=60, maxT=120):
    """[P] lower bounds for gon(C/S) from gon(C) >= base and Lemma 7.12 alone (no numerical inputs):
    a pencil of degree d on C/S either factors through C/<S, j> (j in T - S) or g(C/S) <= m g(C/T) + (m-1)(d-1)."""
    lb = {ckey(S): max(1 if genus(S) == 0 else 2, -(-base // len(S))) for S in REPS}
    tops = {}
    for S in REPS:
        if len(S) > maxS:
            continue
        NS = normaliser(S)
        covered, Ts = set(S), []
        for j in range(n):
            if j in covered:
                continue
            T = closure(GENS[S] + [j])
            if len(T) <= maxT:
                Ts.append((T, {closure(GENS[S] + [x]) for x in T if x not in S}))
            for g in NS:
                jg = MT[MT[g][j]][INV[g]]
                covered.update(MT[s][jg] for s in S)
        tops[ckey(S)] = Ts
    changed = True
    while changed:
        changed = False
        for S in REPS:
            if len(S) > maxS:
                continue
            kS, gS = ckey(S), genus(S)
            for T, K1s in tops[kS]:
                m, gT = len(T) // len(S), genus(T)
                d_cs = max(1, (gS - m * gT - 1) // (m - 1) + 2) if gS > m * gT else 1
                d_f = min(len(K1) // len(S) * lb[ckey(K1)] for K1 in K1s)
                c = min(d_cs, d_f)
                if c > lb[kS]:
                    lb[kS] = c
                    changed = True
    return lb


def strata(dmax=41, dmin=25):
    """surviving (N, K, d) for a gonal pencil of degree dmin <= d <= dmax, classes 0 and 1."""
    LB = propagate()
    allowed = {ckey(S) for S in REPS
               if any(d % len(S) == 0 and d // len(S) >= LB[ckey(S)] for d in range(dmin, dmax + 1))}
    out = []
    for K in REPS:
        for Nn in normal_subgroups(K):
            if ckey(Nn) not in allowed:
                continue
            NG = normaliser(Nn)
            gY = genus(Nn)
            qo = quotient_orders(K, Nn)
            t = so3_type(K, Nn, qo)
            if t is None:
                continue
            m = mults(K, qo)
            _, dspec = bounds(m, True)
            has3 = any(ORD[k] in (3, 6) for k in K)
            for d in range(dmin, dmax + 1):
                if d % len(Nn):
                    continue
                e = d // len(Nn)
                if e < LB[ckey(Nn)]:
                    continue
                why = []
                if has3 and d % 3:
                    why.append("3|d")
                if len(K) % 9 == 0 and d % 9:
                    why.append("9|d")
                if len(K) % 5 == 0 and d % 5:
                    why.append("5|d")
                if dspec is not None and d < dspec:
                    why.append(f"Li-Yau_K>={dspec}")
                moved = len(K) < len(NG)                 # some element of N_G(N) moves the class
                note = ""
                if moved and not why:
                    ks = [k for k in divisors(e) if k < e]   # index of the compositum of two conjugate fields on Y
                    ok, dropped = [], False
                    for k in ks:
                        if k == 1 and gY <= (e - 1) ** 2:
                            ok.append("birational")
                        elif k > 1 and e // k == 2:
                            if len(Nn) * k >= 60:               # Theorem 7.21: elliptic subcovers have degree >= 60
                                ok.append(f"elliptic subcover deg {len(Nn) * k}")
                            else:
                                dropped = True
                        elif k > 1:
                            ok.append(f"subcover deg {len(Nn) * k} genus<={(e // k - 1) ** 2}")
                    if not ok:
                        why.append(f"CS on C/N (g={gY})" + (" + Thm 7.21" if dropped else ""))
                    note = "; ".join(ok)
                out.append((len(Nn), describe(Nn)[0], gY, LB[ckey(Nn)], len(K), t, d, why, note))
    return out


if __name__ == "__main__" and "--strata" in sys.argv:
    print("\nStrata for a gonal pencil of degree 25 <= d <= 41 (classes 0, 1; Corollary 7.23(2)).")
    print("gon(C/N) >= [P]: Lemma 7.12 propagated from gon(C) >= 25 only.  Exclusions: Lemma 7.17 (3|d, 9|d, 5|d),")
    print("Theorem 7.18 (Li-Yau_K), Castelnuovo-Severi for two conjugate pencils on C/N when K < N_G(N), and")
    print("Theorem 7.21 (elliptic subcovers have degree >= 60).")
    print("  |N| orbits(N)        g(C/N) gon>=  |K| K/N   surviving d                     excluded d (reason)")
    rows = {}
    for (nn, orb, gY, lbN, k, t, d, why, note) in strata():
        key = (nn, str(orb), gY, lbN, k, t)
        rows.setdefault(key, []).append((d, why))
    seen = set()
    for key, lst in rows.items():
        surv = sorted({d for d, w in lst if not w})
        exc = sorted({(d, ", ".join(w)) for d, w in lst if w})
        line = (f"  {key[0]:3d} {key[1]:18s} {key[2]:5d} {key[3]:5d}  {key[4]:4d} {key[5]:4s}  {str(surv):30s}  "
                f"{'; '.join(f'{d} ({w})' for d, w in exc)}")
        if line not in seen:
            seen.add(line)
            print(line)
