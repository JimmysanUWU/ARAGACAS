"""Schur-twisted equivariant geometry of the (2,4,7) A7-curves (7_TWISTED.md).

An A7-INVARIANT line bundle need not be A7-linearised: its obstruction (Mumford class) lies in
H^2(A7, C^*) = Z/6, and it is then linearised for the Schur cover 6.A7, the centre Z/6 acting by a character eps.
This script computes, for every such twist:
  (1) 2.A7 < SU(4) by Clifford lifts of the standard 6; 3.A7 by coset enumeration of a central lift of the
      presentation A7 = <x,y | x^3, y^5, (xy)^7, (xyxy^-1)^2, (xy^-2xy^2)^2>; 6.A7 = 2.A7 x_{A7} 3.A7 with explicit
      2-cocycles; its 40 classes and character table (Burnside-Dixon, orthogonality checked).
  (2) twisted holomorphic Lefschetz: the virtual 6.A7-character chi(L) = [H^0] - [H^1] for every local datum
      (eigenvalue of the lift (g_i,0,0) of each branch generator on the fibre over its fixed point);
  (3) the exact degree of the twisted class with that local datum:  deg L = 2520 (n + sum_i r_i/e_i - rho_0),
      lambda_i = exp(2 pi i r_i/e_i),  eps(x1 x2 x3) = exp(2 pi i rho_0)   (Proposition 7.2);
  (4) validation: eps = 1 reproduces equivariant_rr.py (B, B+T, K); both orientation conventions;
  (5) the degree-45 obstruction (Theorem 7.4): Molien series of 2.A7 on V4, f14 and f18 vanish on the 210
      involution lines and are coprime, and every degree-45 spin class forces f14 = f18 = 0 on its image;
  (6) the degree-60 class L60 with 6 in H^0 (Theorem 7.5): eigenvalues of an involution at its 18 fixed points,
      the pencil P(E_-) with >= 18 base points -> gon(C) <= 42, gon(C/tau) <= 21; equivariant Plucker check;
  (7) Lemma 7.3 on the linearised degree-90 classes (which power sums of the 6 must vanish);
  (8) a scan of all classes with forced sections (degree <= 270) for pencils with fixed-point base loci;
  (9) the invariants of 3.A7 on its 6 (Molien) and which of them vanish on the degree-60 model;
  (10) the 6 of 3.A7 explicitly (induced from S5 x Z3) and its invariant cubic in eigencoordinates of a 7-element.
Run:  python3 twisted_rr.py > twisted_rr_output.txt     (about 25 seconds)
"""
import os, sys, time, cmath, itertools
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a7 import A7, mul, inv, order, E
from triples_data import triples

T0 = time.time()
G = list(A7); IDX = {g: i for i, g in enumerate(G)}; n = len(G)
assert G[0] == E
P = np.array(G, dtype=np.int64); key = lambda arr: (arr * (7 ** np.arange(7))).sum(-1)
kk = np.array(sorted(int(key(np.array(g))) for g in G)); kmap = {int(key(np.array(g))): i for i, g in enumerate(G)}
vv = np.array([kmap[x] for x in kk])
MUL = np.array([vv[np.searchsorted(kk, key(P[i][P]))] for i in range(n)], dtype=np.int64)   # MUL[g,h] = g o h
INV = np.array([IDX[inv(g)] for g in G])

# ------------------------------------------------------------------ (1a) 2.A7 < SU(4): Clifford lifts
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1]).astype(complex)
I2 = np.eye(2)
kr = lambda *ms: np.kron(np.kron(ms[0], ms[1]), ms[2])
gam = [kr(sx, I2, I2), kr(sy, I2, I2), kr(sz, sx, I2), kr(sz, sy, I2), kr(sz, sz, sx), kr(sz, sz, sy)]
chir = (-1j) ** 3 * gam[0] @ gam[1] @ gam[2] @ gam[3] @ gam[4] @ gam[5]
w_, V_ = np.linalg.eigh(chir); Qh = V_[:, w_ > 0]                       # half-spin space V4
Fb = np.linalg.qr(np.array([np.eye(7)[i] - np.eye(7)[6] for i in range(6)]).T)[0].T
GT = {(i, j): sum(c * g for c, g in zip(Fb @ ((np.eye(7)[i] - np.eye(7)[j]) / np.sqrt(2)), gam))
      for i in range(7) for j in range(i + 1, 7)}
def spin_lift(p):                                                         # product of reflections, transposition-wise
    q, S = list(p), np.eye(8, dtype=complex)
    for i in range(7):
        if q[i] != i:
            k = q[i]; S = S @ GT[min(i, k), max(i, k)]
            t = list(E); t[i], t[k] = k, i; q = [t[x] for x in q]
    return Qh.conj().T @ S @ Qh
U = np.zeros((n, 4, 4), complex); KEY = np.zeros(n, int)
for i, g in enumerate(G):
    u = spin_lift(g); f = u.flatten(); k = int(np.argmax(np.abs(f) > 1e-7)); ph = f[k] / abs(f[k])
    if ph.real < -1e-9 or (abs(ph.real) <= 1e-9 and ph.imag < 0): u = -u
    U[i] = u; KEY[i] = k
ref = U.reshape(n, 16)[np.arange(n), KEY]
C2 = np.zeros((n, n), dtype=np.int64)
for i in range(n):
    prods = np.einsum('ij,hjk->hik', U[i], U).reshape(n, 16); gh = MUL[i]
    v = prods[np.arange(n), KEY[gh]] / ref[gh]
    assert np.allclose(np.abs(v), 1) and np.allclose(v.imag, 0)
    C2[i] = (v.real < 0)
print(f"(1a) 2.A7 < SU(4): spin lifts of all 2520 elements, 2-cocycle built ({time.time() - T0:.0f}s)")

# ------------------------------------------------------------------ (1b) 3.A7 by coset enumeration
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.fp_groups import FpGroup
from sympy.combinatorics import Permutation, PermutationGroup
Fg, x, y, z = free_group("x y z")
WD = {'x': x, 'X': x**-1, 'y': y, 'Y': y**-1}
def fw(word):
    r = Fg.identity
    for ch in word: r = r * WD[ch]
    return r
a_, b_ = (1, 2, 0, 3, 4, 5, 6), (0, 1, 3, 4, 5, 6, 2)                  # x -> (012), y -> (23456)
for word in ('xxx', 'yyyyy', 'xy' * 7, 'xyxY' * 2, 'xYYxyy' * 2):
    r = E
    for ch in word: r = mul(r, {'x': a_, 'X': inv(a_), 'y': b_, 'Y': inv(b_)}[ch])
    assert r == E
pres = FpGroup(Fg, [z**3, x*z*x**-1*z**-1, y*z*y**-1*z**-1, fw('yyyyy'), fw('xy' * 7),
                    fw('xxx'), fw('xyxY' * 2) * z**-1, fw('xYYxyy' * 2)])
Ct = pres.coset_enumeration([y], max_cosets=400000); Ct.compress(); Ct.standardize()
X3, Y3, Z3 = (np.array([row[2 * i] for row in Ct.table]) for i in range(3))
P3 = PermutationGroup([Permutation(list(X3)), Permutation(list(Y3))])
assert len(Ct.table) == 1512 and P3.order() == 7560 and P3.derived_subgroup().order() == 7560
lift = {0: np.arange(len(X3))}; fr = [0]
while fr:
    nf = []
    for gi in fr:
        for ai, Ap in ((IDX[a_], X3), (IDX[b_], Y3)):
            hi = int(MUL[gi, ai])
            if hi not in lift: lift[hi] = Ap[lift[gi]]; nf.append(hi)
    fr = nf
L3 = np.array([lift[i] for i in range(n)]); Zp = [np.arange(len(X3)), Z3, Z3[Z3]]
lookup = -np.ones((n, len(X3)), dtype=np.int64)
for gi in range(n):
    for c in range(3): lookup[gi, Zp[c][L3[gi][0]]] = c
C3 = np.array([lookup[MUL[gi], L3[:, L3[gi][0]]] for gi in range(n)])
assert (C3 >= 0).all()
print(f"(1b) 3.A7: |<x,y>| = 7560, perfect; relator lift (x^3,(xyxy^-1)^2,(xy^-2xy^2)^2) = (1,z,1) ({time.time() - T0:.0f}s)")

# ------------------------------------------------------------------ (1c) 6.A7, classes, characters
N6 = 6 * n
def split(e): return e % n, (e // n) % 2, e // (2 * n)
def gm(e, f):
    g, s, t = split(e); h, u, w = split(f)
    return int(MUL[g, h] + n * ((s + u + C2[g, h]) % 2 + 2 * ((t + w + C3[g, h]) % 3)))
def gi_(e):
    g, s, t = split(e); h = int(INV[g])
    return int(h + n * ((-s - C2[g, h]) % 2 + 2 * ((-t - C3[g, h]) % 3)))
rng = np.random.default_rng(0)
for _ in range(3000):
    p_, q_, r_ = (int(v) for v in rng.integers(N6, size=3)); assert gm(gm(p_, q_), r_) == gm(p_, gm(q_, r_))
gens6 = [IDX[a_], IDX[b_]]; gens6 += [gi_(g) for g in gens6]
cls = -np.ones(N6, int); reps = []
for e0 in range(N6):
    if cls[e0] >= 0: continue
    cc = len(reps); reps.append(e0); cls[e0] = cc; st = [e0]
    while st:
        yy = st.pop()
        for h in gens6:
            zz = gm(gm(h, yy), gi_(h))
            if cls[zz] < 0: cls[zz] = cc; st.append(zz)
kc = len(reps); sizes = np.bincount(cls, minlength=kc)
INV6 = np.array([gi_(e) for e in range(N6)])
A = np.zeros((kc, kc, kc)); xs = np.arange(N6); xin = INV6[xs]
gA, sA, tA = xin % n, (xin // n) % 2, xin // (2 * n)
for l, zr in enumerate(reps):
    h, u, w = split(zr)
    ys = MUL[gA, h] + n * ((sA + u + C2[gA, h]) % 2 + 2 * ((tA + w + C3[gA, h]) % 3))
    np.add.at(A, (cls[xs], cls[ys], l), 1)
evl, Wv = np.linalg.eig(sum(rng.standard_normal() * A[j] for j in range(kc)))
chars = []
for col in Wv.T:
    om_ = col / col[0]; dd = np.sqrt((N6 / np.sum(np.abs(om_) ** 2 / sizes)).real); chars.append(om_ * dd / sizes)
chars = np.array(chars); chars = chars[np.argsort([c[0].real for c in chars])]
assert np.allclose([[np.sum(sizes * c1 * c2.conj()) / N6 for c2 in chars] for c1 in chars], np.eye(kc), atol=1e-6)
deg = [round(c[0].real) for c in chars]
omega = cmath.exp(2j * cmath.pi / 3)
chiZ = lambda eps, s, t: (-1) ** (eps[0] * s) * omega ** (eps[1] * t)
zc = [cls[n * (s + 2 * t)] for t in range(3) for s in range(2)]
CC = [(0 if abs(c[zc[1]] / c[0] - 1) < 1e-6 else 1, [m for m in range(3) if abs(c[zc[2]] / c[0] - omega ** m) < 1e-6][0])
      for c in chars]
print(f"(1c) 6.A7: {kc} classes; irreducible degrees by central character eps = (2-part, 3-part):")
for eps in sorted(set(CC)):
    print(f"      eps {eps}: {[deg[i] for i in range(kc) if CC[i] == eps]}")

# ------------------------------------------------------------------ (2)-(3) twisted Lefschetz and exact degrees
H = np.arange(n); HI = np.array([gi_(h) for h in range(n)])
hg, hs, ht = HI % n, (HI // n) % 2, HI // (2 * n)
def conj_all(yv):
    g, s, t = split(yv)
    a1 = MUL[H, g]; s1 = (s + C2[H, g]) % 2; t1 = (t + C3[H, g]) % 3
    og = MUL[a1, hg]; return og + n * ((s1 + hs + C2[a1, hg]) % 2 + 2 * ((t1 + ht + C3[a1, hg]) % 3))
def theta(gens, lam, eps, sgn):
    TH = np.zeros(N6, complex)
    for i, gi in enumerate(gens):
        e = order(gi); yv = IDX[gi]; pw = [None, yv]
        for j in range(2, e): pw.append(gm(pw[-1], yv))
        for j in range(1, e):
            ys = conj_all(pw[j]); aa = cmath.exp(sgn * 2j * cmath.pi * j / e)
            val = lam[i] ** j / (1 - 1 / aa) / e
            for s in range(2):
                for t in range(3):
                    idx = (ys % n) + n * (((ys // n) % 2 + s) % 2 + 2 * (((ys // (2 * n)) + t) % 3))
                    np.add.at(TH, idx, val * chiZ(eps, s, t))
    return TH
def mults(TH, eps, D):
    TH = TH.copy()
    for s in range(2):
        for t in range(3): TH[n * (s + 2 * t)] = (D - 135) * chiZ(eps, s, t)
    return {i: np.sum(TH * np.conj(chars[i][cls])) / N6 for i in range(kc) if CC[i] == eps}
def classes(t, eps, sgn):
    """every local datum -> (exact degree in [0,2520), local datum, multiplicities)."""
    a, b, c = triples[t]; gens = (a, inv(b), mul(b, a)); es = [order(x) for x in gens]
    cp = []
    for gi in gens:
        yv = IDX[gi]; zz = yv
        for _ in range(order(gi) - 1): zz = gm(zz, yv)
        g_, s_, t_ = split(zz); assert g_ == 0; cp.append((s_, t_))
    _, s0, t0 = split(gm(gm(IDX[gens[0]], IDX[gens[1]]), IDX[gens[2]]))
    rho0 = Fr(eps[0] * s0, 2) + Fr(eps[1] * t0, 3)
    opts = [[m for m in range(6 * e) if abs(cmath.exp(2j * cmath.pi * m / 6) - chiZ(eps, *cp[i])) < 1e-9]
            for i, e in enumerate(es)]
    out = []
    for ld in itertools.product(*opts):
        lam = [cmath.exp(2j * cmath.pi * ld[i] / (6 * es[i])) for i in range(3)]
        sg = 1 if sgn > 0 else -1
        D = int((2520 * (sum(Fr(sg * ld[i], 6 * es[i]) for i in range(3)) - sg * rho0)) % 2520)
        m = mults(theta(gens, lam, eps, sgn), eps, D)
        assert all(abs(v.imag) < 1e-6 and abs(v.real - round(v.real)) < 1e-6 for v in m.values()), (t, eps, ld, D)
        out.append((D, ld, {i: round(v.real) for i, v in m.items()}, gens, lam))
    return sorted(out, key=lambda r: r[0])

# ---- (4) validation: eps = 1 reproduces equivariant_rr.py
print("\n(4) validation with eps = 1 (equivariant_rr.py):")
a, b, c = triples[0]; g0 = (a, inv(b), mul(b, a))
for name, nD, Dg in (("B", (0, -1, 2), 90), ("B+T", (1, -3, 2), 90), ("K", (1, 3, 6), 270)):
    lam = [cmath.exp(2j * cmath.pi * nD[i] / order(g0[i])) for i in range(3)]
    m = mults(theta(g0, lam, (0, 0), +1), (0, 0), Dg)
    print(f"    chi({name}) = {sorted((deg[i], round(v.real)) for i, v in m.items() if abs(v) > 1e-6)}")

print("\n(2)-(3) invariant classes of degree <= 270 (one per local datum; 56 local data per twist):")
TW = {}
for t in (0, 12):
    for sgn in (+1, -1):
        for eps in [(1, 1), (0, 1), (1, 0), (0, 2), (1, 2)]:
            rows = classes(t, eps, sgn); TW[t, eps, sgn] = rows
            low = [(D, {deg[i]: v for i, v in m.items() if v > 0}) for D, ld, m, _, _ in rows if D <= 270]
            print(f"    triple {t:2d} orient {sgn:+d} eps {eps}: degrees = {sorted(set(r[0] % 90 for r in rows))} mod 90; "
                  f"D<=270 (H^0 forced): {[(D, p) for D, p in low]}")
print(f"    ({time.time() - T0:.0f}s)")

# ------------------------------------------------------------------ (5) degree 45: the spinor model does not exist
print("\n(5) degree-45 spin classes: no 2.A7-map to P^3 (Theorem 7.4)")
UU = np.concatenate([U, -U]); Nc = 120
tot = np.zeros(Nc + 1, complex)
for u in UU:
    cf = np.zeros(Nc + 1, complex); cf[0] = 1
    for e_ in np.linalg.eigvals(u):
        for kk_ in range(1, Nc + 1): cf[kk_] += e_ * cf[kk_ - 1]
    tot += cf
dims = np.round((tot / len(UU)).real).astype(int)
print(f"    Molien series of 2.A7 on V4: {{k: dim}} = { {k: int(v) for k, v in enumerate(dims) if v and k <= 40} }")
cvec = rng.standard_normal(4) + 1j * rng.standard_normal(4); rowsU = np.einsum('j,gjk->gk', cvec, U)
finv = lambda k: (lambda v: np.mean((rowsU @ v) ** k))
f8, f12, f14, f18 = finv(8), finv(12), finv(14), finv(18)
inv_lift = [i for i in range(n) if np.allclose(U[i] @ U[i], -np.eye(4))][0]
wl, Vl = np.linalg.eig(U[inv_lift]); vline = Vl[:, np.isclose(wl, 1j)] @ np.array([0.3 + 0.7j, -1.1 + 0.2j])
gen = rng.standard_normal(4)
print(f"    on an involution eigenline |f14|/|f14(gen)| = {abs(f14(vline)) / abs(f14(gen)):.1e}, "
      f"|f18|/|f18(gen)| = {abs(f18(vline)) / abs(f18(gen)):.1e}  (i^14 = i^18 = -1: forced zero)")
pl, ql = rng.standard_normal(4) + 1j * rng.standard_normal(4), rng.standard_normal(4) + 1j * rng.standard_normal(4)
def roots_on_line(f, k):
    ts = np.exp(2j * np.pi * np.arange(k + 1) / (k + 1))
    return np.roots(np.linalg.solve(np.vander(ts, k + 1, increasing=True), [f(pl + s * ql) for s in ts])[::-1])
r14, r18 = roots_on_line(f14, 14), roots_on_line(f18, 18)
print(f"    f14, f18 restricted to a random line: min distance between their roots = "
      f"{min(abs(p1 - p2) for p1 in r14 for p2 in r18):.3f}  (no common factor)")
for (t, eps, sgn), rows in TW.items():
    if eps != (1, 0) or t != 0: continue
    for D, ld, m, gens, lam in rows:
        if D != 45: continue
        es = [order(g) for g in gens]; forced = []
        for k in (14, 18):
            # invariant section of L^k: lambda_i^k = a_i^{o_i};  minimal invariant divisor of degree 45k?
            o = []
            for i in range(3):
                aa = cmath.exp(sgn * 2j * cmath.pi / es[i])
                o.append([oo for oo in range(es[i]) if abs(lam[i] ** k - aa ** oo) < 1e-9][0])
            mindeg = sum(o[i] * 2520 // es[i] for i in range(3))
            forced.append((k, o, mindeg, 45 * k))
        print(f"    orient {sgn:+d} local datum {ld}: (k, orders mod e, min invariant divisor, 45k) = {forced}"
              f" -> f14 = f18 = 0 on the image")
print("    => image in Z(f14) n Z(f18) of degree 252, which contains the 210 involution lines: 45 > 252 - 210 = 42.")

# ------------------------------------------------------------------ (6) degree 60: the P^5 model and gon <= 42
print("\n(6) the degree-60 class with 6 in H^0 (Theorem 7.5)")
for (t, eps, sgn), rows in TW.items():
    if eps[0] != 0: continue
    for D, ld, m, gens, lam in rows:
        if D != 60 or not any(v > 0 for v in m.values()): continue
        i6 = [i for i, v in m.items() if v > 0][0]
        tau = gens[0]; tidx = IDX[tau]
        t2 = [L for L in (tidx + n * (s + 2 * u) for s in range(2) for u in range(3)) if gm(L, L) in (0, n)][0]
        cnt = {}
        for k_, gi in enumerate(gens):
            e = order(gi); yv = IDX[gi]; pw = [None, yv]
            for j in range(2, e): pw.append(gm(pw[-1], yv))
            for hi in range(n):
                for j in range(1, e):
                    if MUL[MUL[hi, pw[j] % n], INV[hi]] != tidx: continue
                    cj = gm(gm(hi, pw[j]), gi_(hi)); _, s_, u_ = split(cj); _, s2, u2 = split(t2)
                    val = lam[k_] ** j * chiZ(eps, (s2 - s_) % 2, (u2 - u_) % 3)
                    kk2 = (e, round(val.real)); cnt[kk2] = cnt.get(kk2, 0) + 1 / e
        mplus = round(((6 + chars[i6][cls[t2]]) / 2).real)
        # Plucker: vanishing sequences at the three branch types
        wsum = 0; seqs = []
        for k_, gi in enumerate(gens):
            e = order(gi); yv = IDX[gi]; pw = [0, yv]
            for j in range(2, 6 * e + 1): pw.append(gm(pw[-1], yv))
            Nn = 6 * e; chi = [chars[i6][cls[pw[j]]] if j else chars[i6][0] for j in range(Nn)]
            mult = [sum(chi[j] * cmath.exp(-2j * cmath.pi * j * r / Nn) for j in range(Nn)) / Nn for r in range(Nn)]
            aa = cmath.exp(sgn * 2j * cmath.pi / e); ords = []
            for r in range(Nn):
                for _ in range(round(mult[r].real)):
                    mu = cmath.exp(2j * cmath.pi * r / Nn)
                    ords.append([oo for oo in range(e) if abs(mu - lam[k_] * aa ** (-oo)) < 1e-6][0])
            seq = []
            for o in sorted(ords):
                while o in seq: o += e
                seq.append(o)
            seqs.append(sorted(seq)); wsum += (2520 // e) * (sum(seq) - 15)
        print(f"    triple {t:2d} orient {sgn:+d} eps {eps}: chi(L60) = "
              f"{sorted((deg[i], v) for i, v in m.items() if v)}; involution lift on the 6: +1 x{mplus}, -1 x{6 - mplus}; "
              f"its 18 fixed points by (branch e, eigenvalue): { {k: round(v) for k, v in cnt.items()} }; "
              f"vanishing sequences {seqs}, Plucker 4410 - {wsum} = {4410 - wsum}")
print("    => all 18 fixed points of tau lie in P(E_+) = P^3; the pencil P(E_-) has them as base points:")
print("       gon(C) <= 60 - 18 = 42, and the pencil is tau-invariant, so gon(C/tau) <= 21.")

# ------------------------------------------------------------------ (7) remark: power sums on the mu = 90 classes
print("\n(7) Lemma 7.3 for the linearised degree-90 classes and the power sums p_k of the 6 (k = 2..7):")
for name, nD in (("B", (0, -1, 2)), ("B+T", (1, -3, 2))):
    res = []
    for k in range(2, 8):
        mind = sum(((k * nD[i]) % e) * 2520 // e for i, e in enumerate((2, 4, 7)))
        res.append(f"p{k}: min {mind} vs {90 * k} -> {'0' if mind > 90 * k else 'allowed'}")
    print(f"    {name}: " + "; ".join(res))

# ------------------------------------------------------------------ (8) scan: pencils with fixed-point base loci
print("\n(8) scan over all classes with forced sections, degree <= 270 (triple 0): best pencil bounds")
for sgn in (+1, -1):
    TW[0, (0, 0), sgn] = classes(0, (0, 0), sgn)
def lift_inv(gidx, eps):
    for s in range(2):
        for u in range(3):
            L = gidx + n * (s + 2 * u); sq = split(gm(L, L))
            if sq[0] == 0 and abs(chiZ(eps, sq[1], sq[2]) - 1) < 1e-9: return L
    for s in range(2):
        for u in range(3):
            L = gidx + n * (s + 2 * u); sq = split(gm(L, L))
            if abs(chiZ(eps, sq[1], sq[2]) + 1) < 1e-9: return L
def fixed_eigs(gens, lam, eps, vidx, vl):
    out = {}
    for k_, gi in enumerate(gens):
        e = order(gi); yv = IDX[gi]; pw = [None, yv]
        for j in range(2, e): pw.append(gm(pw[-1], yv))
        for hi in range(n):
            for j in range(1, e):
                if MUL[MUL[hi, pw[j] % n], INV[hi]] != vidx: continue
                cj = gm(gm(hi, pw[j]), gi_(hi)); _, s_, u_ = split(cj); _, s2, u2 = split(vl)
                val = lam[k_] ** j * chiZ(eps, (s2 - s_) % 2, (u2 - u_) % 3)
                kk3 = complex(round(val.real, 6), round(val.imag, 6)); out[kk3] = out.get(kk3, 0) + 1 / e
    return {k: round(v) for k, v in out.items()}
best = []
for (t, eps, sgn), rows in TW.items():
    if t != 0: continue
    for D, ld, m, gens, lam in rows:
        Wd = [(i, v) for i, v in m.items() if v > 0]
        if D > 270 or not Wd: continue
        tau = gens[0]; tidx = IDX[tau]
        tp = next(g for g in G if order(g) == 2 and g != tau and mul(g, tau) == mul(tau, g))
        lt, ltp = lift_inv(tidx, eps), lift_inv(IDX[tp], eps); lttp = gm(lt, ltp)
        tr = lambda L: sum(mult * chars[i][cls[L]] for i, mult in Wd)
        hW = sum(deg[i] * v for i, v in Wd)
        fe = fixed_eigs(gens, lam, eps, tidx, lt)
        pw = [0, lt]
        for j in range(2, 4): pw.append(gm(pw[-1], lt))
        for mu in {complex(round(v.real, 6), round(v.imag, 6)) for v in (cmath.exp(2j * cmath.pi * r / 4) for r in range(4))}:
            dim = round((sum((tr(pw[j]) if j else hW) * mu ** (-j) for j in range(4)) / 4).real)
            base = sum(cn for lv, cn in fe.items() if abs(lv - mu) > 1e-6)
            if dim >= 2 and base: best.append((D - base - (dim - 2), D, eps, sgn, f"tau-eigenspace {mu:.0f}: dim {dim}, base {base}"))
        if eps[0] == 0 and gm(lt, ltp) == gm(ltp, lt):
            fe2 = fixed_eigs(gens, lam, eps, IDX[tp], ltp); fe3 = fixed_eigs(gens, lam, eps, IDX[mul(tau, tp)], lttp)
            for c1, c2 in itertools.product((1, -1), repeat=2):
                dimc = round(((hW + c1 * tr(lt) + c2 * tr(ltp) + c1 * c2 * tr(lttp)) / 4).real)
                base = sum(cn for lv, cn in fe.items() if abs(lv - c1) > 1e-6) + \
                    sum(cn for lv, cn in fe2.items() if abs(lv - c2) > 1e-6) + \
                    sum(cn for lv, cn in fe3.items() if abs(lv - c1 * c2) > 1e-6)
                if dimc >= 2 and base: best.append((D - base - (dimc - 2), D, eps, sgn, f"Klein char ({c1},{c2}): dim {dimc}, base {base}"))
best.sort(key=lambda r: r[0])
for bb_ in best[:6]:
    print(f"    pencil degree <= {bb_[0]:3d} from the degree-{bb_[1]} class, twist {bb_[2]}, orientation {bb_[3]:+d}: {bb_[4]}")

# ------------------------------------------------------------------ (9) the invariant cubic fourfold
print("\n(9) invariants of 3.A7 on its 6 and their restrictions to the degree-60 model (Lemma 7.3)")
i6 = [i for i in range(kc) if CC[i] == (0, 1) and deg[i] == 6][0]
Nm = 30; molien = np.zeros(Nm + 1, complex)
for c_, rep in enumerate(reps):
    o = 1; zz = rep
    while zz != 0: zz = gm(zz, rep); o += 1
    pw = [0, rep]
    for j in range(2, o): pw.append(gm(pw[-1], rep))
    chi = [chars[i6][cls[pw[j]]] if j else chars[i6][0] for j in range(o)]
    eig = []
    for r in range(o):
        ml = sum(chi[j] * cmath.exp(-2j * cmath.pi * j * r / o) for j in range(o)) / o
        eig += [cmath.exp(2j * cmath.pi * r / o)] * round(ml.real)
    cf = np.zeros(Nm + 1, complex); cf[0] = 1
    for e_ in eig:
        for k_ in range(1, Nm + 1): cf[k_] += e_ * cf[k_ - 1]
    molien += sizes[c_] * cf
md = np.round((molien / N6).real).astype(int)
print(f"    Hilbert series of C[6]^(3.A7): { {k: int(v) for k, v in enumerate(md) if v} }")
for D, ld, m, gens, lam in TW[0, (0, 1), +1]:
    if D == 60 and any(v > 0 for v in m.values()):
        for k in [k for k, v in enumerate(md) if v and 0 < k <= Nm]:
            o = [[oo for oo in range(e) if abs(lam[i] ** k - cmath.exp(2j * cmath.pi * oo / e)) < 1e-9][0]
                 for i, e in enumerate((2, 4, 7))]
            mind = sum(o[i] * 2520 // e for i, e in enumerate((2, 4, 7)))
            print(f"    degree {k:2d} ({md[k]:3d} invariants): least invariant divisor {mind:4d} vs {60 * k:4d} -> "
                  f"{'all vanish on the model' if mind > 60 * k else 'restrictions span <= ' + str(1 + (60 * k - mind) // 2520) + ' dim'}")

# ------------------------------------------------------------------ (10) the cubic fourfold explicitly
print("\n(10) explicit 6 of 3.A7 (induced from S5 x Z3) and its invariant cubic in 7-eigencoordinates")
w3 = cmath.exp(2j * cmath.pi / 3)
def m3(p, q): return (int(MUL[p[0], q[0]]), (p[1] + q[1] + int(C3[p[0], q[0]])) % 3)
def i3(p): h = int(INV[p[0]]); return (h, (-p[1] - int(C3[p[0], h])) % 3)
chi6 = lambda p: chars[i6][cls[p[0] + 2 * n * p[1]]]
def lift_ord(g):
    o = order(g)
    for t in range(3):
        p = (IDX[g], t); zz = p
        for _ in range(o - 1): zz = m3(zz, p)
        if zz == (0, 0): return p
Ssub = {(0, 0)}; fr = [(0, 0)]; gl = [lift_ord((1, 2, 3, 4, 0, 5, 6)), lift_ord((1, 0, 2, 3, 4, 6, 5))]
while fr:
    nf = []
    for p in fr:
        for gg in gl:
            q = m3(p, gg)
            if q not in Ssub: Ssub.add(q); nf.append(q)
    fr = nf
assert len(Ssub) == 120                                          # a complement S5 inside 3.A7
Sg = {p[0]: p for p in Ssub}; coset = {}; creps = []
for p in ((g, t) for g in range(n) for t in range(3)):
    if p in coset: continue
    j = len(creps); creps.append(p)
    for sp in Ssub:
        for u in range(3): coset[m3(p, m3(sp, (0, u)))] = j
def rho21(q):
    M = np.zeros((21, 21), complex)
    for i, r in enumerate(creps):
        qr = m3(q, r); j = coset[qr]; kq = m3(i3(creps[j]), qr)
        M[j, i] = w3 ** ((kq[1] - Sg[kq[0]][1]) % 3)
    return M
ga3, gb3 = (IDX[(1, 2, 0, 3, 4, 5, 6)], 0), lift_ord((0, 1, 3, 4, 5, 6, 2))
mats = {(0, 0): np.eye(21, dtype=complex)}; fr = [(0, 0)]; Ra3, Rb3 = rho21(ga3), rho21(gb3)
while fr:
    nf = []
    for p in fr:
        for gg, R in ((ga3, Ra3), (gb3, Rb3)):
            q = m3(p, gg)
            if q not in mats: mats[q] = mats[p] @ R; nf.append(q)
    fr = nf
Pr = sum(np.conj(chi6(p)) * M for p, M in mats.items()) * 6 / 7560
Uu, sv, _ = np.linalg.svd(Pr); Bb = Uu[:, :6]
assert np.allclose(sv[:6], 1) and np.allclose(sv[6:], 0, atol=1e-9)
R6 = lambda p: Bb.conj().T @ mats[p] @ Bb
assert all(abs(np.trace(R6(p)) - chi6(p)) < 1e-8 for p in list(mats)[:300])
mon3 = list(itertools.combinations_with_replacement(range(6), 3))
Xp = rng.standard_normal((80, 6)) + 1j * rng.standard_normal((80, 6))
mv = lambda X: np.array([[np.prod([x[i] for i in m_]) for m_ in mon3] for x in X])
Mx = np.vstack([mv(Xp @ R6(q).T) - mv(Xp) for q in (ga3, gb3)])
_, svx, vh = np.linalg.svd(Mx); Fc = vh[-1].conj()
print(f"    invariant cubics: null space of dimension {int(np.sum(svx < 1e-8))}")
c7 = lift_ord((1, 2, 3, 4, 5, 6, 0)); h7 = (0, 2, 4, 6, 1, 3, 5)
hl = [(IDX[h7], t) for t in range(3) if m3(m3((IDX[h7], t), (IDX[h7], t)), (IDX[h7], t)) == (0, 0)][0]
ev7, Ev = np.linalg.eig(R6(c7)); kof = [int(round(cmath.phase(e_) / (2 * cmath.pi) * 7)) % 7 for e_ in ev7]
Ev = Ev[:, [kof.index(k) for k in range(1, 7)]]; Vb = np.zeros((6, 6), complex)
for st in (1, 3):
    k = st; v = Ev[:, k - 1] / np.linalg.norm(Ev[:, k - 1])
    for _ in range(3): Vb[:, k - 1] = v; v = R6(hl) @ v; k = (4 * k) % 7
allowed = [m_ for m_ in mon3 if sum(i + 1 for i in m_) % 7 == 0]
Yp = rng.standard_normal((40, 6)) + 1j * rng.standard_normal((40, 6))
Fv = lambda x: np.dot(Fc, [np.prod([x[i] for i in m_]) for m_ in mon3])
cf, *_ = np.linalg.lstsq(np.array([[np.prod([y[i] for i in m_]) for m_ in allowed] for y in Yp]),
                         np.array([Fv(Vb @ y) for y in Yp]), rcond=None)
dct = dict(zip(allowed, cf / cf[allowed.index((0, 1, 3))]))
be, ga_, de = dct[(2, 4, 5)], dct[(0, 0, 4)], dct[(0, 2, 2)]
assert np.isclose(dct[(1, 1, 2)], ga_) and np.isclose(dct[(3, 3, 5)], ga_) and np.isclose(dct[(1, 5, 5)], de)
print("    F = x1x2x4 + b x3x5x6 + g (x1^2x5 + x2^2x3 + x4^2x6) + d (x1x3^2 + x2x6^2 + x4x5^2)")
print(f"    scale-free data: g^3/b = {(ga_**3 / be).real:.12f}  [(23 - 7 sqrt21)/16 = {(23 - 7 * 21**0.5) / 16:.12f}]")
print(f"                     d^3/b^2 = {(de**3 / be**2).real:.12f}  [(23 + 7 sqrt21)/16 = {(23 + 7 * 21**0.5) / 16:.12f}]")
print(f"                     g d/b = {ga_ * de / be:.12f}  [(5/4) exp(-i pi/3) = {1.25 * cmath.exp(-1j * cmath.pi / 3):.12f}]")
print(f"\ndone ({time.time() - T0:.0f}s)")
