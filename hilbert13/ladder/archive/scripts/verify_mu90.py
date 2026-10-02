"""mu(A7) = 90 (MU90.md): checks for the finite and representation-theoretic inputs.

mu(A7) = least degree of an A7-linearised base-point-free line bundle with h^0 >= 2 on a faithful connected A7-curve.
Known: mu >= 72 (GPT Round 5, Thm 2.2, verified in ROUND5_REVIEW.md) and mu <= 90 (equivariant_rr.py), and the only
lattice-compatible values below 90 are 72, 78, 84 (Round 5, Prop 2.4).  This script checks everything needed to exclude
78, 72 and 84.

 1. Candidates.  If mu = n, the complete series is birational onto Y (minimality); Y has an A7-signature with
    N(Y) | n, g <= pi(n, h^0 - 1), h^0 >= 6.  Enumerate signatures and the h^0 range; the section representation W has no
    trivial summand (an invariant section would vanish on whole orbits, degree n < min orbit).
 2. W contains the standard 6 (any 6-subsystem is base-point free and birational by minimality):
      power sums p_k(x) = sum x_j^k are invariant sections of L^k; if kn is not a sum of orbit sizes, p_k|C = 0.
      n = 72: p2,p3,p4,p6 vanish -> image in {e1=e2=e3=e4=e6=0} = ordered roots of t^7 - u t^2 - v, an irreducible
              curve of degree 144 (monodromy of t^7 - t^2 is S7: 6 simple critical points + a 7-cycle), contradiction.
      n = 84: p2,p3,p4,p7 vanish -> e1..e4 = 0, p7 = 7 e7 -> prod x_j = 0 on C, impossible by transitivity.
 3. W = 10 or 10b (only where h^0 >= 10):  W* = Sym^2 V4 (V4 = half-spin rep of 2.A7 < Spin(6) = SU(4); checked:
      Sym^2 V4 = 10b, wedge^2 V4 = 6).  Points of Y are symmetric 4x4 matrices Q_p.
      (a) det is an A7-invariant quartic; 4n is not a sum of orbit sizes, so det(Q_p) = 0: rank <= 3.
      (b) generic rank 2: wedge^2 Q_p = (v^w)^2 gives a linearised base-point-free 6-system N with 72 <= deg N <= n,
          deg N in N(C)Z, hence deg N = n: excluded by step 2.
      (c) parity lemma: at p over a branch with stabiliser <c>, Q_p is an eigenvector of c with eigenvalue lam = a^{-k}
          (a = rotation, k = multiplicity of the divisor at that fibre; both orientation conventions tested).  If lam is
          not the square of an eigenvalue of the lift of c on V4, every element of that eigenspace has even rank.
          Then rank(Q_p) = 2 on that whole orbit.  Generic rank 1 is impossible (rank would be 1 everywhere); generic
          rank 3 is impossible because the adjugate map gives kappa*O(2) = L^3(-B) with deg B <= 3n < orbit size.
      (d) generic rank 1 where no parity branch exists ((3,5,6), n = 84): phi = v2 o psi, psi: C -> P^3 birational of
          degree n/2 = 42, on no quadric (Sym^2 V4 irreducible), so g <= d(d-3)/6 + 1 = 274 (Halphen; Gruson-Peskine)
          < 379.
"""
import os
import sys
import itertools
from fractions import Fraction as F
from math import comb, lcm
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a7 import A7, mul, inv, order, cycle_type
from equivariant_rr import cls, chars
from signatures_spectrum import curves

G = list(A7)
idx = {g: i for i, g in enumerate(G)}


def pi(n, r):
    k, t = divmod(n - 1, r - 1)
    return (r - 1) * comb(k, 2) + k * t


def representable(x, gens):
    S = [False] * (x + 1); S[0] = True
    for y in range(1, x + 1):
        S[y] = any(y >= g and S[y - g] for g in gens)
    return S[x]


# ------------------------------------------------------------------ the half-spin representation of 2.A7
Q7, _ = np.linalg.qr(np.vstack([np.ones(7), np.eye(7)[:6]]).T)
B6 = Q7[:, 1:7]
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1, -1]).astype(complex)
I2 = np.eye(2)
def kron(*ms):
    out = np.eye(1)
    for m in ms: out = np.kron(out, m)
    return out
gam = [kron(s1, I2, I2), kron(s2, I2, I2), kron(s3, s1, I2), kron(s3, s2, I2), kron(s3, s3, s1), kron(s3, s3, s2)]
chir = (1j) ** 3 * gam[0] @ gam[1] @ gam[2] @ gam[3] @ gam[4] @ gam[5]
w_, V_ = np.linalg.eigh((chir + chir.conj().T) / 2)
Pp = V_[:, w_ > 0]


def lift(g):
    """an element of 2.A7 < SU(4) over g (sign ambiguous), via products of reflections in Clifford(6)."""
    g = list(g); ts = []
    for i in range(7):
        while g[i] != i:
            j = g[i]; ts.append((i, j)); g[i], g[j] = g[j], g[i]
    Sm = np.eye(8, dtype=complex)
    for (i, j) in ts:
        u = B6.T @ (np.eye(7)[i] - np.eye(7)[j]) / np.sqrt(2)
        Sm = Sm @ sum(u[k] * gam[k] for k in range(6))
    return Pp.conj().T @ Sm @ Pp


lifts = {g: lift(g) for g in G}
for name, fn in (("Sym^2", lambda M: (np.trace(M) ** 2 + np.trace(M @ M)) / 2),
                 ("wedge^2", lambda M: (np.trace(M) ** 2 - np.trace(M @ M)) / 2)):
    tr = np.array([fn(lifts[g]) for g in G])
    dec = {nm: (sum(tr[i] * np.conj(chars[nm][cls(g)]) for i, g in enumerate(G)) / 2520) for nm in chars}
    dec = {k: round(v.real, 9) for k, v in dec.items() if abs(v) > 1e-9}
    print(f"0. {name} V4 = {dec}")
    assert dec == ({"10b": 1.0} if name == "Sym^2" else {"6": 1.0})
print("   det(lift) = 1 for all g:", np.allclose([np.linalg.det(lifts[g]) for g in G], 1))


iu = [(i, j) for i in range(4) for j in range(i, 4)]
_rng = np.random.default_rng(0)


def eigen_data(g, k, conv):
    """lam = exp(conv 2 pi i k/e).  Returns (no_rank1, no_rank3): no_rank1 if lam is not the square of an eigenvalue of
    the lift (then every element of the lam-eigenspace of Sym^2 V4 has even rank); no_rank3 if the lam-eigenspace
    contains no rank-3 symmetric matrix with det = 0 (tested on det = 0 along random lines, and on the generic rank)."""
    e = order(g)
    lam = np.exp(conv * 2j * np.pi * k / e)
    ev, P = np.linalg.eig(lifts[g])
    no_rank1 = bool(np.all(np.abs(ev ** 2 - lam) > 1e-6))
    # eigenspace of Q -> s Q s^T in the eigenbasis of s: Q_ab may be nonzero iff ev_a ev_b = lam
    supp = [(a, b) for a in range(4) for b in range(a, 4) if abs(ev[a] * ev[b] - lam) < 1e-6]
    if not supp:
        return no_rank1, True
    def mat(c):
        Q = np.zeros((4, 4), complex)
        for (a, b), x in zip(supp, c): Q[a, b] = Q[b, a] = x
        return Q
    def rk(Q, tol):
        sv = np.linalg.svd(Q, compute_uv=False); return int(np.sum(sv > tol * max(sv[0], 1e-300)))
    d = len(supp)
    generic = max(rk(mat(_rng.normal(size=d) + 1j * _rng.normal(size=d)), 1e-9) for _ in range(6))
    if generic <= 2:
        return no_rank1, True
    if generic == 3:
        return no_rank1, False
    ranks = set()
    for _ in range(60):
        U, W = mat(_rng.normal(size=d) + 1j * _rng.normal(size=d)), mat(_rng.normal(size=d) + 1j * _rng.normal(size=d))
        A = np.linalg.solve(W, U)                    # det(U + tW) = det(W) det(A + t)  ->  t = -eig(A)
        for t in -np.linalg.eigvals(A):
            ranks.add(rk(U + t * W, 1e-6))
    return no_rank1, 3 not in ranks


# ------------------------------------------------------------------ 1-3
ORD = [2, 3, 4, 5, 6, 7]
UNRESOLVED = []
MAXIND = {o: max(7 - len(cycle_type(x)) for x in G if order(x) == o) for o in ORD}
for n in (72, 78, 84):
    print(f"\n== n = {n}:  pi(n,5) = {pi(n, 5)}")
    for h in (0, 1):
        for cnt in range(1, 7):
            for es in itertools.combinations_with_replacement(ORD, cnt):
                g = 1 + 1260 * (2 * h - 2 + sum(1 - F(1, e) for e in es))
                N = 2520 // lcm(*es)
                if not (g.denominator == 1 and 136 <= g <= pi(n, 5) and n % N == 0):
                    continue
                g = int(g)
                if sum(MAXIND[e] for e in es) < 12:          # seven-sheeted Riemann-Hurwitz test
                    continue
                if cnt == 3 and not any(curves(*p) for p in set(itertools.permutations(es))):
                    continue
                orbits = sorted(set([2520 // e for e in es] + [2520]))
                forb = [k for k in range(1, 10) if not representable(k * n, orbits)]
                assert min(orbits) > n
                hmax = max(r for r in range(5, 60) if pi(n, r) >= g) + 1
                assert hmax <= 12          # so W has no 14/15/21/35 and is one of 6, 10, 10b, 6+6
                reps = [W for W in ("6", "10", "10b", "6+6") if {"6": 6, "10": 10, "10b": 10, "6+6": 12}[W] <= hmax]
                line = f"   {es} h={h} g={g} N={N} orbits={orbits} h0<={hmax} W in {reps}; forbidden k<10: {forb}"
                # step 2
                if n == 72:
                    assert {2, 3, 4, 6} <= set(forb); line += "  [6: p2,p3,p4,p6 -> trinomial curve, deg 144]"
                else:
                    assert {2, 3, 4, 7} <= set(forb); line += "  [6: p2,p3,p4,p7 -> e7 = 0]"
                print(line)
                if hmax < 10:
                    continue
                # step 3 (only 3-point signatures reach h0 >= 10 here)
                assert h == 0 and cnt == 3
                assert 4 in forb                                            # (a) det vanishes
                assert [d for d in range(72, n + 1) if d % N == 0] == [n]    # (b) deg N = n
                for perm in set(itertools.permutations(es)):
                    for (a, b, osz) in curves(*perm):
                        gens = (a, inv(b), mul(b, inv(a)))
                        ee = [order(x) for x in gens]; degs = [2520 // e for e in ee]
                        for n1 in range(ee[0]):
                            for n2 in range(ee[1]):
                                rest = n - n1 * degs[0] - n2 * degs[1]
                                if rest % degs[2]:
                                    continue
                                nD = (n1, n2, rest // degs[2])
                                for conv in (+1, -1):
                                    dat = [eigen_data(gens[i], nD[i], conv) for i in range(3)]
                                    big = [degs[i] > 3 * n for i in range(3)]
                                    r1 = [ee[i] for i in range(3) if dat[i][0]]
                                    r3 = [ee[i] for i in range(3) if dat[i][1] and big[i]]
                                    d2 = n // 2
                                    halphen = g > d2 * (d2 - 3) / 6 + 1          # rank 1: psi(C) in P^3, degree n/2, on no quadric
                                    if halphen:
                                        r1 = r1 + ["Halphen"]
                                    status = "OK" if (r1 and r3) else "OPEN(rank 1)" if r3 else "OPEN(rank 3)" if r1 else "OPEN(rank 1,3)"
                                    UNRESOLVED.append((perm, nD, conv)) if status != "OK" else None
                                    print(f"      curve {[cycle_type(x) for x in gens]} D={nD} conv {conv:+d}: rank-1 killed at e={r1}, "
                                          f"rank-3 killed at e={r3}  -> {status}")
print("\nunresolved (curve ordering, class, convention):", sorted(set(UNRESOLVED)))
assert not UNRESOLVED
print("all checks passed: no linearised moving bundle of degree 72, 78 or 84; with 72 <= mu <= 90 (Round 5, equivariant_rr.py): mu(A7) = 90")
