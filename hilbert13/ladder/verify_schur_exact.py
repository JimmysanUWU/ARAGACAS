"""Exact certificate for the degree-60 series and gon(C) <= 42 (Theorem 7.5, Cor. 7.6) [X], independent of twisted_rr.py.

1. The 6 of 3.A7: the ATLAS matrices (3A7G1-Ar6B0, over Z[w], w = E(3)) generate a group of order 7560 whose scalars
   are {1, w, w^2}, and A -> x0, B -> y0 is a well-defined homomorphism onto A7 (checked on every edge of the Cayley graph).
2. Every g in A7 of order prime to 3 has a unique lift of the same order; chi_6 on it is computed exactly in Z[w].
3. Holomorphic Lefschetz, exactly in Q(zeta_84):  tr(g^ | H^0 - H^1) = sum_{p in Fix g} lambda(g^, p) / (1 - a_p^{-1}),
   z central acts by eps(z)(deg - 135); a branch point of type i (e_i = 2, 4, 7) has local datum lambda_i = zeta_{e_i}^{r_i}
   for the order-e_i lift, and deg = 2520 (n + sum r_i/e_i - rho_0) with eps(x1^ x2^ x3^) = e^{2 pi i rho_0}  (Prop. 7.1).
   Validation: O and K_C (eps = 1) give <chi, chi> = 12 = 1 + |H^0(K)|^2 and the right trivial multiplicities.
4. For eps = the central character of the 6 and every local datum of degree 60: the multiplicity of the 6 in chi(L) and
   <chi(L), chi(L)>.  The class with m_6 = 1 has lambda = +1 at all 18 fixed points of an involution tau, and its lift
   acts on the 6 with trace 2, i.e. eigenvalues (+1)^4 (-1)^2.  So |E_-(tau)| is a pencil with the 18 points as base points,
   of moving degree <= 42.
All four (2,4,7) classes.  The orientation is fixed by K_C; the other convention is the complex conjugate.
"""
import itertools
from fractions import Fraction
import numpy as np
import sympy as sp
from a7 import A7, E, mul, inv, order, cyc
from triples_data import triples

# ---------------------------------------------------------------- Z[w] matrices: pairs (P, Q) for P + Q w, w^2 = -1 - w
def zmul(X, Y):
    (A0, A1), (B0, B1) = X, Y
    return (A0 @ B0 - A1 @ B1, A0 @ B1 + A1 @ B0 - A1 @ B1)
I6 = (np.eye(6, dtype=np.int64), np.zeros((6, 6), dtype=np.int64))
def zpow(X, k):
    R = I6
    for _ in range(k): R = zmul(R, X)
    return R
def scalar(k):                                            # w^k I
    return [I6, (np.zeros((6, 6), np.int64), np.eye(6, dtype=np.int64)), (-np.eye(6, dtype=np.int64), -np.eye(6, dtype=np.int64))][k % 3]
def key(X): return X[0].tobytes() + X[1].tobytes()
def zm(rows):                                             # rows of (p, q) pairs
    return (np.array([[e[0] for e in r] for r in rows], np.int64), np.array([[e[1] for e in r] for r in rows], np.int64))
# ATLAS 3A7G1-Ar6B0 (brauer.maths.qmul.ac.uk/Atlas/alt/A7/gap0/3A7G1-Ar6B0.g):
#   A = [[1,0,0,0,0,0],[0,w^2,0,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[-1+w,0,1-w,-w,w,1],[2,0,-1,-1,0,-1]]
#   B = [[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[1,0,0,0,0,0],[-1,1,-w,0,w,1]]
o, one, w, w2 = (0, 0), (1, 0), (0, 1), (-1, -1)
A = zm([[one, o, o, o, o, o], [o, w2, o, o, o, o], [o, o, one, o, o, o], [o, o, o, o, o, one],
        [(-1, 1), o, (1, -1), (0, -1), w, one], [(2, 0), o, (-1, 0), (-1, 0), o, (-1, 0)]])
B = zm([[o, one, o, o, o, o], [o, o, one, o, o, o], [o, o, o, one, o, o], [o, o, o, o, one, o],
        [one, o, o, o, o, o], [(-1, 0), one, (0, -1), o, w, one]])

# ---------------------------------------------------------------- 1. the group and the homomorphism to A7
x0 = cyc((0, 1, 2))
def word(p, q, s):
    r = E
    for ch in s: r = mul(r, {'x': p, 'y': q, 'Y': inv(q)}[ch])
    return r
REL = ['xxx', 'yyyyy', 'xy' * 7, 'xyxYxyxY', 'xYYxyyxYYxyy']
y0 = next(q for q in A7 if order(q) == 5 and all(word(x0, q, r) == E for r in REL))
Binv = zpow(B, 4)
for r in REL:                                             # the ATLAS generators satisfy the relators up to scalars
    M = I6
    for ch in r: M = zmul(M, {'x': A, 'y': B, 'Y': Binv}[ch])
    assert any(key(M) == key(scalar(k)) for k in range(3)), r
elems = {key(I6): (I6, E)}; fr = [key(I6)]
while fr:
    nf = []
    for k_ in fr:
        M, p = elems[k_]
        for G_, g_ in ((A, x0), (B, y0)):
            N, q = zmul(M, G_), mul(p, g_)
            kn = key(N)
            if kn in elems: assert elems[kn][1] == q, "not a homomorphism"
            else: elems[kn] = (N, q); nf.append(kn)
    fr = nf
fib = {}
for M, p in elems.values(): fib.setdefault(p, []).append(M)
assert len(elems) == 7560 and len(fib) == 2520 and all(len(v) == 3 for v in fib.values())
assert sum(1 for v in fib[E] if any(key(v_) == key(scalar(k)) for k in range(3) for v_ in [v])) == 3
print("1. ATLAS 3A7G1-Ar6B0: |<A,B>| = 7560, kernel of the map to A7 = scalars {1,w,w^2}; homomorphism checked on all edges")

def lift(g):                                              # the lift of order order(g) (order prime to 3)
    e = order(g)
    L = [M for M in fib[g] if key(zpow(M, e)) == key(I6)]
    assert len(L) == 1
    return L[0]
def trace(M): return (int(np.trace(M[0])), int(np.trace(M[1])))

# ---------------------------------------------------------------- exact arithmetic in Q(zeta_84)
N = 84; X = sp.symbols('X'); PHI = sp.Poly(sp.cyclotomic_poly(N, X), X, domain='QQ')
def cz(k): return sp.Poly(X ** (k % N), X, domain='QQ').rem(PHI)
def cval(p, q): return (sp.Poly(p, X, domain='QQ') + q * cz(28)).rem(PHI)          # p + q w, w = zeta_84^28
def cconj(P): return sp.Poly(sum(c * X ** ((N - m[0]) % N) for m, c in P.terms()), X, domain='QQ').rem(PHI)
def cinv(P): return sp.Poly(sp.invert(P.as_expr(), PHI.as_expr(), domain='QQ'), X, domain='QQ')
def rational(P):
    P = P.rem(PHI); assert P.degree() <= 0, P
    return Fraction(str(P.as_expr()))

# ---------------------------------------------------------------- Lefschetz data per triple
def lefschetz_data(cls):
    x1, x2, x3 = triples[cls]
    xs, es = (x1, x2, x3), (2, 4, 7)
    assert mul(mul(x1, x2), x3) == E
    pw = [[E] for _ in range(3)]
    for i in range(3):
        for _ in range(es[i] - 1): pw[i].append(mul(pw[i][-1], xs[i]))
    pos = [{pw[i][k]: k for k in range(es[i])} for i in range(3)]
    reps, seen = [], set()
    for g in A7:
        if g in seen or order(g) not in (2, 4, 7): continue
        cl = {mul(mul(inv(h), g), h) for h in A7}; seen |= cl; reps.append((g, len(cl)))
    data = []
    for g, size in reps:
        fixed = {}                                        # (type i, k) -> number of fixed points where g acts as x_i^k
        for h in A7:
            c = mul(mul(inv(h), g), h)
            for i in range(3):
                if c in pos[i]: fixed[(i, pos[i][c])] = fixed.get((i, pos[i][c]), 0) + 1
        fixed = {ik: v // es[ik[0]] for ik, v in fixed.items()}
        data.append((g, size, fixed, trace(lift(g))))
    z0 = zmul(zmul(lift(x1), lift(x2)), lift(x3))
    j = [k for k in range(3) if key(z0) == key(scalar(k))][0]
    return data, j

def lefschetz(data, r, s, deg, eps_order, chi6_exp):
    """returns (<chi, chi>, <chi, chi_6> or None).  r = (r1, r2, r3); s = +-1 orientation (a_i = zeta_{e_i}^s);
    eps_order = 1 (trivial) or 3 (eps = central character of the 6, or its conjugate if chi6_exp = -1)."""
    es = (2, 4, 7); nG = 2520 * eps_order
    tot2 = Fraction(eps_order * (deg - 135) ** 2); tot6 = None
    t6 = sp.Poly(0, X, domain='QQ')
    if eps_order == 3: tot6 = Fraction(eps_order * 6 * (deg - 135))      # sum over the centre of conj(chi6(z)) tr(z)
    for g, size, fixed, (tp, tq) in data:
        tr = sp.Poly(0, X, domain='QQ')
        for (i, k), cnt in fixed.items():
            e = es[i]; lam = cz((N // e) * r[i] * k); den = sp.Poly(1, X, domain='QQ') - cz(-(N // e) * s * k)
            tr = (tr + cnt * lam * cinv(den)).rem(PHI)
        tot2 += eps_order * size * rational((tr * cconj(tr)).rem(PHI))
        if eps_order == 3:
            ch = cval(tp, tq) if chi6_exp == 1 else cconj(cval(tp, tq))
            t6 = (t6 + eps_order * size * cconj(ch) * tr).rem(PHI)
    nrm = tot2 / nG
    m6 = (tot6 + rational(t6)) / nG if eps_order == 3 else None
    return nrm, m6

print("2-4. Lefschetz in Q(zeta_84):")
ok_all = True
for cls in (0, 1, 12, 14):
    data, j = lefschetz_data(cls)
    inv_tr = {order(g): tr for g, _, _, tr in data}
    assert inv_tr[2] == (2, 0), "involution lift does not have trace 2"
    s = 1                       # orientation: K_C (lambda_i = a_i^{-1}) must get degree 270 from Prop. 7.1
    rK = tuple((-s) % e for e in (2, 4, 7))
    degK = 2520 * ((Fraction(rK[0], 2) + Fraction(rK[1], 4) + Fraction(rK[2], 7)) % 1)
    assert degK == 270                                     # (the opposite orientation gives 2250: it is the conjugate convention)
    nO, _ = lefschetz(data, (0, 0, 0), s, 0, 1, 1)
    nK, _ = lefschetz(data, rK, s, 270, 1, 1)
    assert nO == 12 and nK == 12, (nO, nK)
    found = []
    for chi6_exp in (1, -1):
        rho0 = Fraction((j * chi6_exp) % 3, 3)
        for r in itertools.product(range(2), range(4), range(7)):
            if (Fraction(r[0], 2) + Fraction(r[1], 4) + Fraction(r[2], 7) - rho0 - Fraction(60, 2520)) % 1 == 0:
                nrm, m6 = lefschetz(data, r, s, 60, 3, chi6_exp)
                assert nrm.denominator == 1 and m6.denominator == 1        # a virtual character, as it must be
                found.append((chi6_exp, r, nrm, m6))
    good = [f for f in found if f[3] >= 1]
    ok = len(good) >= 1 and all(f[1][0] == 0 and f[1][1] % 2 == 0 for f in good)
    ok_all &= ok
    print(f"   class {cls:2d}: O and K_C give <chi,chi> = 12; degree-60 classes: "
          + "; ".join(f"eps of {'6' if c == 1 else '6bar'}, r = {r}: <chi,chi> = {nr}, mult. of that 6 = {m}" for c, r, nr, m in found)
          + f"  =>  a 6 in H^0 with lambda = +1 at the 18 fixed points of tau: {ok}")
print(f"   the order-2 lift acts on the 6 with trace 2: eigenvalues (+1)^4 (-1)^2")
print(f"==> exact: every (2,4,7) curve has a degree-60 class with 6 in H^0, and gon(C) <= 42: {ok_all}")
