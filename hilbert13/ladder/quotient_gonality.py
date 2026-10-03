"""Gonality of the quotients C/K and pencils pulled back from them (Chapter 7, §7.10) [N] + [P].

1. Holomorphic differentials on C as vector-valued weight-2 forms for Delta(2,4,7) in the disk model (order-7 point at 0):
   F(Xz) X'(z) = rho(a) F(z) at the order-2 point, rho(c) a_n = zeta^(n+1) a_n.  One solve per irreducible rho in H^0(K_C).
2. Differentials on Y = C/K: pair F with K-fixed vectors of rho^T.  Canonical points of Y are sampled over the whole curve
   as (w, g) with w in the central disk and g in A7:  x = (u^T rho(g) F(w)).
3. Tests on the canonical curve of Y:
   - Noether: Sym^2 has rank 2g-1 if Y is hyperelliptic, 3g-3 otherwise;
   - Petri: the quadrics through Y have Jacobian rank g-3 at a general point if Y is trigonal, g-2 otherwise;
   - Green-Lazarsfeld nonvanishing (a theorem): Cliff(Y) <= p implies K_{p,2}(Y, K_Y) != 0.  So K_{p,2} = 0 gives gon >= p+3.
4. Castelnuovo-Severi over the subgroup lattice [P]: for S < T of index m, a pencil of degree d on C/S either factors through
   C/<S,j> for some j in T \\ S, or g(C/S) <= m g(C/T) + (m-1)(d-1).  Iterated to a fixed point over all 36 classes of
   subgroups of order <= 120, with the inputs from 3 and the base bound |K| gon(C/K) >= gon(C) >= 25.

Usage: python3 quotient_gonality.py [class]   (class index into triples_data.triples; default 0)
"""
import sys, itertools, math, time
import numpy as np
from a7 import A7, E, mul, inv, order, cyc
from chartab import TABLE, idx as cidx
from triples_data import triples

T0 = time.time()
G = list(A7); GI = {g: i for i, g in enumerate(G)}; n = len(G)
cls = int(sys.argv[1]) if len(sys.argv) > 1 else 0
a, b, c = triples[cls]
zeta = np.exp(2j * np.pi / 7)
d72 = np.arccosh(np.cos(np.pi / 4) / np.sin(np.pi / 7)); P2 = np.tanh(d72 / 2) + 0j     # order-2 point
Tm = lambda p: np.array([[1, p], [np.conj(p), 1]]) / np.sqrt(1 - abs(p) ** 2)
MX = Tm(P2) @ np.diag([1j, -1j]) @ np.linalg.inv(Tm(P2))
Xz = lambda z: (MX[0, 0] * z + MX[0, 1]) / (MX[1, 0] * z + MX[1, 1])
dXz = lambda z: 1.0 / (MX[1, 0] * z + MX[1, 1]) ** 2

def closure(gens):
    S = {E}; fr = [E]
    while fr:
        fr = [z for x in fr for y in gens if (z := mul(x, y)) not in S and not S.add(z)]
    return S

# ---- 1. irreducible modules and forms
def coset_module(H):
    cos, where = [], {}
    for g in G:
        if g in where: continue
        cos.append(g)
        for h in H: where[mul(g, h)] = len(cos) - 1
    def M(g):
        P = np.zeros((len(cos), len(cos)))
        for j, r in enumerate(cos): P[where[mul(g, r)], j] = 1
        return P
    return M
def wedge_module(k):
    subs = list(itertools.combinations(range(7), k)); si = {s: i for i, s in enumerate(subs)}
    def M(g):
        P = np.zeros((len(subs), len(subs)))
        for j, s in enumerate(subs):
            img = [g[x] for x in s]
            inversions = sum(img[u] > img[v] for u in range(k) for v in range(u + 1, k))
            P[si[tuple(sorted(img))], j] = (-1) ** inversions
        return P
    return M
def isotypic_copy(Mfun, chi, deg):
    mats = {g: Mfun(g) for g in G}
    Pr = sum(np.conj(chi[cidx[g]]) * mats[g] for g in G) * deg / 2520
    w, U = np.linalg.eigh((Pr + Pr.conj().T) / 2); B = U[:, w > 0.5]
    assert B.shape[1] == deg
    return {g: B.conj().T @ mats[g] @ B for g in G}
chars = lambda deg: [t for t in TABLE if round(t[0].real) == deg]
L25 = closure([cyc((1, 2, 3, 4, 5)), cyc((1, 6), (2, 5))])
Q1 = closure([cyc((0, 1, 2)), cyc((3, 4, 5)), cyc((0, 3, 1, 4), (2, 5))])
reps = {'15': isotypic_copy(wedge_module(2), chars(15)[0], 15),
        '10': isotypic_copy(wedge_module(3), chars(10)[0], 10), '10b': isotypic_copy(wedge_module(3), chars(10)[1], 10),
        '21': isotypic_copy(coset_module(L25), chars(21)[0], 21), '35': isotypic_copy(coset_module(Q1), chars(35)[0], 35)}

def solve_forms(rho, M=110, r=0.62, seed=1):
    d = rho[a].shape[0]
    ev, V = np.linalg.eig(rho[c]); js = np.round(np.angle(ev) / (2 * np.pi / 7)).astype(int) % 7
    B = {j: np.linalg.qr(V[:, js == j])[0] for j in range(7) if np.any(js == j)}
    cols = [(m, (m + 1) % 7) for m in range(M + 1) if (m + 1) % 7 in B]
    nunk = sum(B[j].shape[1] for _, j in cols)
    rng = np.random.default_rng(seed); npts = int(2.5 * nunk / d) + 40
    zs = P2 + 0.2 * np.sqrt(rng.uniform(0, 1, npts)) * np.exp(2j * np.pi * rng.uniform(0, 1, npts))
    A = np.vstack([np.hstack([(dXz(z) * (Xz(z) / r) ** m * np.eye(d) - (z / r) ** m * rho[a]) @ B[j] for m, j in cols]) for z in zs])
    s, Vh = np.linalg.svd(A, full_matrices=False)[1:]
    sa = s[::-1][:8]; k = int(np.argmax(sa[1:] / sa[:-1])) + 1
    assert sa[k] / sa[k - 1] > 1e5, sa
    def form(row):
        coef, t = [], 0
        for m, j in cols:
            w = B[j].shape[1]; coef.append((m, B[j] @ row[t:t + w])); t += w
        return lambda z: sum(np.outer((np.atleast_1d(z) / r) ** m, v) for m, v in coef)
    return [form(Vh[-1 - i].conj()) for i in range(k)], sa[:k + 1]
forms = {}
print(f"class {cls}.  1. Holomorphic differentials (smallest singular values; a gap > 1e5 separates the forms):")
for name, rho in reps.items():
    forms[name], sv = solve_forms(rho)
    print(f"   rho = {name:3s}: {len(forms[name])} forms   {np.array2string(sv, precision=1)}", flush=True)
assert sum(len(v) * reps[k][a].shape[0] for k, v in forms.items()) == 136

# ---- 2. canonical points of quotients
rngS = np.random.default_rng(11); NS = 600
SW = 0.5 * np.sqrt(rngS.uniform(0, 1, NS)) * np.exp(2j * np.pi * rngS.uniform(0, 1, NS))
SG = [G[t] for t in rngS.integers(0, n, NS)]
FW = {name: [F(SW) for F in Fs] for name, Fs in forms.items()}
def canonical_points(K):
    cols = []
    for name, Fs in forms.items():
        rho = reps[name]; w, U = np.linalg.eig(sum(rho[k].T for k in K) / len(K))
        for us in U[:, np.abs(w - 1) < 1e-8].T:
            for fv in FW[name]:
                cols.append([us @ (rho[SG[p]] @ fv[p]) for p in range(NS)])
    W = np.array(cols).T
    return W / np.linalg.norm(W, axis=1, keepdims=True)
def numrank(s, floor=None):
    s = s / s[0]; k = int(np.argmax(s[:-1] / np.maximum(s[1:], 1e-300))) + 1
    return (k if floor is None or s[k] < floor else len(s)), s[k - 1], (s[k] if k < len(s) else 0)

# ---- 3. Noether, Petri, Koszul
def noether_petri(W):
    g = W.shape[1]; pairs = [(i, j) for i in range(g) for j in range(i, g)]
    S2 = np.array([W[:, i] * W[:, j] for i, j in pairs]).T
    s, Vh = np.linalg.svd(S2)[1:]; rk, lo, hi = numrank(s, 1e-8)
    out = {'g': g, 'rank': rk, 'cut': (lo, hi), 'gon': 2 if rk == 2 * g - 1 else 3}
    if rk == 3 * g - 3 and g >= 5:
        Q = Vh[rk:].conj(); assert np.max(np.abs(S2 @ Q.T)) < 1e-8
        jr = []
        for x in W[:6]:
            J = np.zeros((len(Q), g), complex)
            for k, (i, j) in enumerate(pairs): J[:, i] += Q[:, k] * x[j]; J[:, j] += Q[:, k] * x[i]
            sj = np.linalg.svd(J / np.linalg.norm(J, axis=1, keepdims=True), compute_uv=False)
            jr.append(int(np.sum(sj > 1e-7 * sj[0])))
        out['petri'] = jr
        if min(jr) == g - 2: out['gon'] = 4
    return out
def koszul_Kp2(W, p):
    g = W.shape[1]
    P2_ = np.array([W[:, i] * W[:, j] for i in range(g) for j in range(i, g)]).T
    U2, s2 = np.linalg.svd(P2_, full_matrices=False)[:2]; B2 = U2[:, :int(np.sum(s2 > 1e-9 * s2[0]))]
    P3_ = np.hstack([W[:, [i]] * B2 for i in range(g)])
    U3, s3 = np.linalg.svd(P3_, full_matrices=False)[:2]; B3 = U3[:, :int(np.sum(s3 > 1e-9 * s3[0]))]
    assert B2.shape[1] == 3 * g - 3 and B3.shape[1] == 5 * g - 5
    M1 = [B2.conj().T @ (W[:, [i]] * W) for i in range(g)]
    M2 = [B3.conj().T @ (W[:, [i]] * B2) for i in range(g)]
    def dmat(q, Ms, rin, rout):
        Iq = list(itertools.combinations(range(g), q)); Jq = {J: k for k, J in enumerate(itertools.combinations(range(g), q - 1))}
        D = np.zeros((len(Jq) * rout, len(Iq) * rin), complex)
        for col, I in enumerate(Iq):
            for t, i in enumerate(I):
                r0 = Jq[I[:t] + I[t + 1:]] * rout
                D[r0:r0 + rout, col * rin:(col + 1) * rin] += (-1) ** t * Ms[i]
        return D
    dp = dmat(p, M2, 3 * g - 3, 5 * g - 5); dp1 = dmat(p + 1, M1, g, 3 * g - 3)
    ra = numrank(np.linalg.svd(dp, compute_uv=False)); rb = numrank(np.linalg.svd(dp1, compute_uv=False))
    return dp.shape[1] - ra[0] - rb[0], ra, rb
gen = lambda *gs: closure(list(gs))
TESTS = [('C3^2:2 (order 18)', gen(cyc((0, 1, 2)), cyc((3, 4, 5)), cyc((0, 1), (3, 4))), None),
         ('5:4', gen(cyc((0, 1, 2, 3, 4)), cyc((1, 2, 4, 3), (5, 6))), None),
         ('D12', gen(cyc((0, 1, 2)), cyc((0, 1), (5, 6)), cyc((3, 4), (5, 6))), None),
         ('3:4', gen(cyc((0, 1, 2)), cyc((0, 1), (3, 4, 5, 6))), None),
         ('A4 (on 4 points)', gen(cyc((0, 1, 2)), cyc((0, 1), (2, 3))), None),
         ('A4 (diagonal)', gen(cyc((0, 1, 2), (3, 4, 5)), cyc((0, 3), (1, 4))), None),
         ('C6xC2', gen(cyc((0, 1, 2)), cyc((3, 4), (5, 6)), cyc((3, 5), (4, 6))), None),
         ('D10', gen(cyc((0, 1, 2, 3, 4)), cyc((1, 4), (2, 3))), 2),
         ('C3xC3', gen(cyc((0, 1, 2)), cyc((3, 4, 5))), 2),
         ('D8', gen(cyc((0, 1, 2, 3), (4, 5)), cyc((0, 2), (4, 5))), 3)]
inputs = []
print("\n3. Quotient canonical curves.  Noether: rank Sym^2 = 2g-1 (hyperelliptic) or 3g-3.  Petri: Jacobian rank g-3 (trigonal) or g-2.")
for lab, K, p in TESTS:
    W = canonical_points(K); r = noether_petri(W); g = r['g']; gon = r['gon']
    line = f"   {lab:18s} |K|={len(K):2d} g={g:2d}: rank Sym^2 {r['rank']} (cut {r['cut'][0]:.0e}|{r['cut'][1]:.0e})"
    if 'petri' in r: line += f", Petri ranks {r['petri']}"
    if p:
        for q in range(2, p + 1):
            dim, ra, rb = koszul_Kp2(W, q)
            line += f"; dim K_{q},2 = {dim} (cuts {ra[1]:.0e}|{ra[2]:.0e}, {rb[1]:.0e}|{rb[2]:.0e})"
            if dim == 0: gon = max(gon, q + 3)
    print(line + f"  =>  gon >= {gon}", flush=True)
    inputs.append((K, gon))

# ---- 4. the lattice
MT = np.array([[GI[mul(x, y)] for y in G] for x in G], dtype=np.int32)
INV = np.array([GI[inv(x)] for x in G]); ORD = np.array([order(x) for x in G]); E0 = GI[E]
FIXN = np.array([{2: 18, 4: 2, 7: 3}.get(int(o), 0) for o in ORD])     # fixed points on C of an element of each order
def sub(gens):
    S = {E0}; fr = [E0]
    while fr:
        fr = [y for x in fr for h in gens if (y := int(MT[x, h])) not in S and not S.add(y)]
    S = frozenset(S); GENS.setdefault(S, list(gens)); return S
GENS, _ck = {}, {}
def ckey(S):                                               # canonical representative of the conjugacy class
    if S not in _ck:
        A = np.array(sorted(S))
        _ck[S] = min(map(tuple, np.sort(MT[MT[INV[:, None], A[None, :]], np.arange(n)[:, None]], axis=1)))
    return _ck[S]
def genus(S): return (270 - int(FIXN[list(S)].sum())) // (2 * len(S)) + 1
classes, frontier = {}, []
for i in range(n):
    S = sub([i])
    if ckey(S) not in classes: classes[ckey(S)] = S; frontier.append(S)
while frontier:
    nf = []
    for S in frontier:
        if len(S) > 60: continue
        for j in range(n):
            if j not in S and len(T := sub(GENS[S] + [j])) <= 120 and ckey(T) not in classes:
                classes[ckey(T)] = T; nf.append(T)
    frontier = nf
reps_ = sorted(classes.values(), key=len)
lb = {ckey(S): max(1 if genus(S) == 0 else 2, math.ceil(25 / len(S))) for S in reps_}
for K, gon in inputs:
    k = ckey(sub([GI[x] for x in K])); lb[k] = max(lb[k], gon)
over = {ckey(S): list({T for j in range(n) if j not in S and len(T := sub(GENS[S] + [j])) <= 120}) for S in reps_ if len(S) <= 60}
changed = True
while changed:
    changed = False
    for S in reps_:
        if len(S) > 60: continue
        kS, gS, cand = ckey(S), genus(S), []
        for T in over[kS]:
            m, gT = len(T) // len(S), genus(T)
            d_cs = max(1, (gS - m * gT - 1) // (m - 1) + 2) if gS > m * gT else 1
            d_f = min(len(K1) // len(S) * lb[ckey(K1)] for K1 in over[kS] if K1 <= T)
            cand.append(min(d_cs, d_f))
        if cand and max(cand) > lb[kS]: lb[kS] = max(cand); changed = True
def describe(S):
    par = list(range(7))
    def f(x):
        while par[x] != x: x = par[x]
        return x
    for h in GENS[S]:
        for x in range(7): par[f(x)] = f(G[h][x])
    orb = sorted(np.bincount([f(x) for x in range(7)]).tolist()); orb = [o for o in orb if o]
    hist = {int(o): int(k) for o, k in enumerate(np.bincount(ORD[list(S)], minlength=8)) if k and o > 1}
    return orb, hist
print(f"\n4. Castelnuovo-Severi over the {len(reps_)} classes of subgroups of order <= 120 (rows: |K| <= 60)")
print("   |K|  g(C/K)  gon(C/K) >=  |K|.gon >=        orbits on 7 points   element orders")
for S in reps_:
    if len(S) > 60: continue
    b_ = lb[ckey(S)]; orb, hist = describe(S)
    print(f"   {len(S):3d}  {genus(S):5d}  {b_:8d}     {len(S) * b_:6d} {'' if len(S) * b_ >= 42 else '(<42)':6s}   {str(orb):22s} {hist}")
print(f"\n({time.time() - T0:.0f} s)")
