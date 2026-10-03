"""mu(A7) = 90: EXACT re-verification of the rank steps of 5_ACCESSORY.md section 5.4 (replaces the floating-point eigen_data of
the first floating-point check, in git history at b148338; everything below is integer / rational / polynomial arithmetic).

Setting.  2.A7 < Spin(6) = SU(4) acts on the half-spin V4 with wedge^2 V4 = 6 (standard) and Sym^2 V4 = 10 or 10b.
For g in A7 of order e and a lift s, let mu_a = zeta_{2e}^{k_a} (a = 1..4) be the eigenvalues of s.

 (E1) The exponents k_a are FORCED: {mu_a mu_b : a < b} is the spectrum of g on the 6 (permutation spectrum minus one 1)
      and prod mu_a = 1 (SU(4)).  We solve this exactly.  The solutions are unique up to the lift sign,
      except for 7-cycles, where the two dual half-spin spectra must also be included.
 (E2) In an eigenbasis of s the action Q -> s Q s^T is diagonal, Q_ab -> mu_a mu_b Q_ab.  So the lam-eigenspace of
      Sym^2 V4 is EXACTLY the space of symmetric matrices supported on the pattern P = {(a,b) : mu_a mu_b = lam}, and
      ranks are preserved by the change of basis.  (Sym^2 V4^* is the case lam -> lam^-1; both orientation conventions
      are also lam <-> lam^-1, so testing lam and lam^-1 covers every convention.)
 (E3) For a pattern P (variables x_ab, a <= b in P):
        rank 1 occurs  <=>  P has a loop (a,a)   [v v^T supported in P needs S x S in P];
        rank 3 occurs  <=>  some 3x3 minor m is nonzero on E and, if det is not identically 0 on E, some such m is not in
                            rad(det), i.e. not divisible by every irreducible factor of det (factorisation over Q;
                            radical membership is unchanged by extension to C).
Then the case analysis of that first check is rerun with these exact predicates.
"""
import os, sys, itertools
from fractions import Fraction as F
from math import comb, lcm
from functools import lru_cache
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a7 import A7, mul, inv, order, cycle_type
from signatures_spectrum import curves

G = list(A7)


def pi(n, r):
    k, t = divmod(n - 1, r - 1)
    return (r - 1) * comb(k, 2) + k * t


def representable(x, gens):
    S = [False] * (x + 1); S[0] = True
    for y in range(1, x + 1):
        S[y] = any(y >= g and S[y - g] for g in gens)
    return S[x]


# ------------------------------------------------------------------ (E1) exact lift spectra
@lru_cache(maxsize=None)
def lift_exponents(ct):
    """cycle type -> sorted tuple of exponents k_a mod 2e (one of the two lifts), e = order."""
    e = lcm(*ct)
    L = 2 * e
    # spectrum of g on the permutation module: for each cycle of length l, exponents (L/l) j, j = 0..l-1; remove one 0
    spec = sorted(sum(([(L // l) * j % L for j in range(l)] for l in ct), []))
    spec.remove(0)
    sols = []
    for ks in itertools.combinations_with_replacement(range(L), 4):
        if sum(ks) % L:
            continue
        if sorted((ks[a] + ks[b]) % L for a, b in itertools.combinations(range(4), 2)) == spec:
            sols.append(ks)
    k0 = sols[0]
    shift = lambda ks: tuple(sorted((x + e) % L for x in ks))
    neg = lambda ks: tuple(sorted((-x) % L for x in ks))
    if ct == (7,):
        # two spectra up to sign: S = {1,z,z^2,z^4} (z = zeta_7) and its inverse.  They are the two half-spin modules
        # V4, V4^* (both have wedge^2 = 6); which one a given 7-cycle sees is not decided by wedge^2.  Replacing the
        # spectrum by its inverse is the same as lam -> lam^-1 (pattern(mu^-1, lam) = pattern(mu, lam^-1)), and every
        # other cycle type has a self-dual spectrum (checked below), so testing ONE spectrum against both lam and
        # lam^-1 covers every possibility.
        assert set(sols) == {k0, shift(k0), neg(k0), shift(neg(k0))} and len(set(sols)) == 4, sols
    else:
        assert set(sols) == {k0, shift(k0)}, (ct, sols)
        assert neg(k0) in (k0, shift(k0)), ("spectrum not self-dual", ct)
    return k0, L


# ------------------------------------------------------------------ (E3) exact rank predicates on a pattern
@lru_cache(maxsize=None)
def pattern_ranks(P):
    """P: frozenset of (a,b), a <= b.  Returns (rank1_possible, rank3_possible)."""
    rank1 = any(a == b for a, b in P)
    if not P:
        return False, False
    X = {ab: sp.Symbol(f"x{ab[0]}{ab[1]}") for ab in P}
    M = sp.zeros(4, 4)
    for (a, b), x in X.items():
        M[a, b] = M[b, a] = x
    det = sp.expand(M.det())
    minors = []
    for rows in itertools.combinations(range(4), 3):
        for cols in itertools.combinations(range(4), 3):
            m = sp.expand(M.extract(list(rows), list(cols)).det())
            if m != 0:
                minors.append(m)
    if not minors:
        return rank1, False
    if det == 0:
        return rank1, True
    gens = list(X.values())
    facs = [f for f, _ in sp.factor_list(det, *gens)[1]]
    for m in minors:
        if any(sp.rem(m, f, *gens) != 0 for f in facs):
            return rank1, True              # m does not vanish on some component of {det = 0}: rank 3 occurs there
    return rank1, False                     # every 3-minor vanishes on {det = 0}: no rank-3 element


def pattern(ks, L, lam_exp):
    return frozenset((a, b) for a in range(4) for b in range(a, 4) if (ks[a] + ks[b] - lam_exp) % L == 0)


def eigen_exact(g, k):
    """For both lam = zeta_e^{+k} and zeta_e^{-k}: (no_rank1, no_rank3) must hold for BOTH to be used."""
    ks, L = lift_exponents(tuple(cycle_type(g)))
    e = order(g)
    out = []
    for sgn in (+1, -1):
        lam_exp = (sgn * k * (L // e)) % L
        r1, r3 = pattern_ranks(pattern(ks, L, lam_exp))
        out.append((not r1, not r3))
    return out


# ------------------------------------------------------------------ Sym^2 V4 is irreducible (exact character check)
def sym2_check():
    """Exact character norm by reduction modulo cyclotomic polynomials.

    The squared norm at each cycle type is a rational integer, including
    both dual 7-cycle spectra.  No nsimplify or numerical root recognition.
    """
    tot = 0
    from collections import Counter
    z = sp.Symbol('z')
    cnt = Counter(tuple(cycle_type(g)) for g in G)
    for ct, c in cnt.items():
        ks, L = lift_exponents(ct)
        s = [(ks[a] + ks[b]) % L for a in range(4) for b in range(a, 4)]
        polynomial = sum(z**((p-q)%L) for p in s for q in s)
        val = sp.rem(polynomial,sp.cyclotomic_poly(L,z),z)
        assert val.is_Integer, (ct,val)
        tot += c * val
    return sp.Rational(tot,len(G))


if __name__ == "__main__":
    print("(E1) exact lift exponents (k_a mod 2e), one lift per cycle type:")
    for ct in sorted(set(tuple(cycle_type(g)) for g in G), key=lambda c: (lcm(*c), c)):
        ks, L = lift_exponents(ct)
        print(f"   {ct}: exponents {ks} mod {L}")
    n2 = sym2_check()
    print(f"    <chi_Sym2V4, chi_Sym2V4> = {n2}  (1 = irreducible; it is 10-dimensional, so it is 10 or 10b)")
    assert n2 == 1

    ORD = [2, 3, 4, 5, 6, 7]
    MAXIND = {o: max(7 - len(cycle_type(x)) for x in G if order(x) == o) for o in ORD}
    UNRES = []
    ncases = 0
    for n in (72, 78, 84):
        print(f"\n== n = {n}")
        for h in (0, 1):
            for cnt in range(1, 7):
                for es in itertools.combinations_with_replacement(ORD, cnt):
                    g = 1 + 1260 * (2 * h - 2 + sum(1 - F(1, e) for e in es))
                    N = 2520 // lcm(*es)
                    if not (g.denominator == 1 and 136 <= g <= pi(n, 5) and n % N == 0):
                        continue
                    g = int(g)
                    if sum(MAXIND[e] for e in es) < 12:
                        continue
                    if cnt == 3 and not any(curves(*p) for p in set(itertools.permutations(es))):
                        continue
                    orbits = sorted(set([2520 // e for e in es] + [2520]))
                    forb = [k for k in range(1, 10) if not representable(k * n, orbits)]
                    hmax = max(r for r in range(5, 60) if pi(n, r) >= g) + 1
                    assert min(orbits) > n and hmax <= 12
                    assert ({2, 3, 4, 6} if n == 72 else {2, 3, 4, 7}) <= set(forb)
                    print(f"   {es} h={h} g={g} N={N} h0<={hmax}", "(W = 6 only)" if hmax < 10 else "(W = 6, 6+6, 10, 10b)")
                    if hmax < 10:
                        continue
                    assert h == 0 and cnt == 3 and 4 in forb
                    assert [d for d in range(72, n + 1) if d % N == 0] == [n]
                    for perm in set(itertools.permutations(es)):
                        for (a, b, osz) in curves(*perm):
                            gens = (a, inv(b), mul(b, inv(a)))
                            ee = [order(x) for x in gens]; degs = [2520 // e for e in ee]
                            for n1 in range(ee[0]):
                                for n2_ in range(ee[1]):
                                    rest = n - n1 * degs[0] - n2_ * degs[1]
                                    if rest % degs[2]:
                                        continue
                                    nD = (n1, n2_, rest // degs[2])
                                    ncases += 1
                                    dat = [eigen_exact(gens[i], nD[i]) for i in range(3)]
                                    d2 = F(n, 2)
                                    halph = g > d2 * (d2 - 3) / 6 + 1
                                    # the convention (W* = Sym^2 V4 or its dual; orientation) is GLOBAL: lam_i = a_i^{c n_i}
                                    # with one c = +-1 for all three branches.  Each c must be closed separately.
                                    msg, ok = [], True
                                    for c in (0, 1):
                                        r1 = [ee[i] for i in range(3) if dat[i][c][0]] + (["Halphen"] if halph else [])
                                        r3 = [ee[i] for i in range(3) if dat[i][c][1] and degs[i] > 3 * n]
                                        ok = ok and bool(r1) and bool(r3)
                                        msg.append(f"c={'+-'[c]}: r1 {r1}, r3 {r3}")
                                    strong = (any(all(d[0] for d in dat[i]) for i in range(3)) or halph) and \
                                        any(all(d[1] for d in dat[i]) and degs[i] > 3 * n for i in range(3))
                                    if not ok:
                                        UNRES.append((perm, nD))
                                    print(f"      {[cycle_type(x) for x in gens]} D={nD}: {'; '.join(msg)} -> "
                                          f"{'OK' if ok else 'OPEN'}{' (convention-free)' if strong else ''}")
    print(f"\ncases: {ncases}; unresolved: {UNRES}")
    assert not UNRES
    print("EXACT: every (curve, class, convention) case is closed; with 72 <= mu <= 90: mu(A7) = 90")
