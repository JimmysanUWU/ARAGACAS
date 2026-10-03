"""Elliptic subcovers of the (2,4,7) A7-curves (Chapter 7, Proposition 7.10).

(1) H^1(C,Z) with its cup product and A7-action, from the dessin of each class: the Delta-complex of 5040
    triangles, a tree-cotree basis, and the Alexander-Whitney cup product. Checks: antisymmetric, unimodular,
    H^1 = 3(10+10b) + 2*15 + 2*21 + 4*35.
(2) For rho = 15, 21 (multiplicity 2, multiplicity space M of rank 2) every v in V_rho(Q) gives a rank-2
    sub-Hodge lattice Lambda_v = (v (x) M) cap H^1(C,Z), i.e. an elliptic subcover C -> E_v of degree
    n(v) = |cup product on a basis of Lambda_v|. With K = A_5 (two points fixed) for 15 and K = L_2(5) for 21,
    (V_rho)^K is a line; its translates span V_rho, and Phi(c) = <sum c_i h_i x_1, sum c_i h_i x_2> is the form
    with n(v) = |Phi(v)| / g(v), where g(v) = [Lambda_v : Z Av + Z Bv] is read off the finite glue H/(L0+L0).
    Output: the K-fixed plane has n = 60 (the Galois quotient C -> C/K), and no plane has n < 60.
"""
import sys, math, itertools
from collections import deque
from fractions import Fraction
import numpy as np, flint
from a7 import A7, mul, inv, cyc, generated
from chartab import TABLE, idx as CIDX
from triples_data import triples

names = ["1", "6", "10", "10b", "14a", "14b", "15", "21", "35"]
G = list(A7); n = len(G); I = {g: i for i, g in enumerate(G)}
K_of = {"15": generated([cyc((0, 1, 2, 3, 4)), cyc((0, 1, 2))]),          # A_5 fixing 5 and 6
        "21": generated([cyc((0, 1, 2, 3, 4)), cyc((0, 5), (1, 4))])}     # L_2(5) on P^1(F_5), fixing 6

def cohomology(a, b):
    """Delta-complex over 0 < 1 < oo: edges E01_g, E1oo_g, E0oo_g; U_g = (E01_g, E1oo_g, E0oo_g), L_g = (E01_g, E1oo_gb, E0oo_ga^-1)."""
    Rb = np.array([I[mul(g, b)] for g in G]); Rai = np.array([I[mul(g, inv(a))] for g in G])
    def cosets(x):
        lab = -np.ones(n, int); k = 0
        for i, g in enumerate(G):
            if lab[i] < 0:
                h = g
                while lab[I[h]] < 0: lab[I[h]] = k; h = mul(h, x)
                k += 1
        return lab, k
    v0, n0 = cosets(a); v1, n1 = cosets(b); vi, ni = cosets(mul(a, b)); nv = n0 + n1 + ni
    tail = np.concatenate([v0, v1 + n0, v0]); head = np.concatenate([v1 + n0, vi + n0 + n1, vi + n0 + n1]); ne = 3 * n
    gi = np.arange(n)
    F = np.vstack([np.stack([gi, n + gi, 2 * n + gi], 1), np.stack([gi, n + Rb, 2 * n + Rai], 1)]); nf = 2 * n
    assert nv - ne + nf == -270
    adj = [[] for _ in range(nv)]
    for e in range(ne): adj[tail[e]].append(e); adj[head[e]].append(e)
    inT = np.zeros(ne, bool); seen = np.zeros(nv, bool); seen[0] = True; order_v = [0]; par_e = -np.ones(nv, int); dq = deque([0])
    while dq:
        u = dq.popleft()
        for e in adj[u]:
            w = head[e] if tail[e] == u else tail[e]
            if not seen[w]: seen[w] = True; inT[e] = True; par_e[w] = e; order_v.append(w); dq.append(w)
    fadj = [[] for _ in range(ne)]
    for f in range(nf):
        for e in F[f]: fadj[e].append(f)
    inTs = np.zeros(ne, bool); fseen = np.zeros(nf, bool); fseen[0] = True; forder = [0]; fpar = -np.ones(nf, int); dq = deque([0])
    while dq:
        f = dq.popleft()
        for e in F[f]:
            if inT[e]: continue
            f2 = [x for x in fadj[e] if x != f][0]
            if not fseen[f2]: fseen[f2] = True; inTs[e] = True; fpar[f2] = e; forder.append(f2); dq.append(f2)
    left = np.where(~inT & ~inTs)[0]; assert len(left) == 272
    sgn = np.array([1, 1, -1])
    Bas = np.zeros((ne, 272), np.int64); Bas[left] = np.eye(272, dtype=np.int64)      # cocycles vanishing on the tree
    for f in reversed(forder[1:]):
        e = fpar[f]; es = F[f]; j = list(es).index(e)
        Bas[e] = -sgn[j] * sum(sgn[i] * Bas[es[i]] for i in range(3) if i != j)
    assert not np.any(Bas[F[:, 0]] + Bas[F[:, 1]] - Bas[F[:, 2]])
    def coords(X):
        f = np.zeros((nv, X.shape[1]), X.dtype)
        for w in order_v[1:]:
            e = par_e[w]; f[w] = f[tail[e]] + X[e] if head[e] == w else f[head[e]] - X[e]
        return X[left] - (f[head[left]] - f[tail[left]])
    Om = Bas[:n].T @ (Bas[n:2 * n] - Bas[n + Rb])                                      # cup product on [C] = sum U - sum L
    def act(h, X):
        p = np.array([I[mul(inv(h), g)] for g in G]); return X[np.concatenate([p, n + p, 2 * n + p])]
    return Bas, coords, Om, act

def fz(M): return flint.fmpz_mat([[int(x) for x in r] for r in M])
def dual_rows(S):
    """basis (columns) of {y in Q^d : S y in Z^m}, S rational of full column rank."""
    S = [[Fraction(x) for x in r] for r in S]; d = len(S[0]); q = 1
    for r in S:
        for x in r: q = q * x.denominator // math.gcd(q, x.denominator)
    rows = [r for r in fz([[x * q for x in r] for r in S]).hnf().tolist() if any(r)]; assert len(rows) == d
    Hi = flint.fmpq_mat(flint.fmpz_mat(rows)).inv()
    return np.array([[Fraction(int(x.p), int(x.q)) * q for x in r] for r in Hi.tolist()], dtype=object)
def kernel_dim_mod(M, p):
    Mp = flint.nmod_mat([[int(x) % p for x in r] for r in M], p); return Mp.ncols() - Mp.rank()
def short_vectors_exist(Gram, bound):
    """is there y != 0 in Z^d with y^T Gram y <= bound?  Exact Fincke-Pohst: Gram = U^T D U with U unit upper
    triangular over Q, so y^T Gram y = sum_i D_i (y_i + sum_{j>i} U_ij y_j)^2, and every branch is pruned by an exact test."""
    d = len(Gram); D = [Fraction(0)] * d; Uu = [[Fraction(int(i == j)) for j in range(d)] for i in range(d)]
    for i in range(d):
        D[i] = Fraction(Gram[i][i]) - sum(D[k] * Uu[k][i] ** 2 for k in range(i))
        assert D[i] > 0, "Gram not positive definite"
        for j in range(i + 1, d):
            Uu[i][j] = (Fraction(Gram[i][j]) - sum(D[k] * Uu[k][i] * Uu[k][j] for k in range(i))) / D[i]
    y = [0] * d
    def rec(i, rem):
        if i < 0:
            return any(y)
        c = -sum(Uu[i][j] * y[j] for j in range(i + 1, d))       # need D_i (y_i - c)^2 <= rem
        s = math.isqrt(math.ceil(rem / D[i])) + 1                 # |y_i - c| <= s
        for yi in range(math.floor(c) - s, math.ceil(c) + s + 1):
            t = D[i] * (yi - c) ** 2
            if t <= rem:
                y[i] = yi
                if rec(i - 1, rem - t): return True
        y[i] = 0; return False
    return rec(d - 1, Fraction(bound))

for cls in (0, 1, 12, 14):
    a, b, c = triples[cls]
    Bas, coords, Om, act = cohomology(a, b)
    chi = {}
    for h in G:
        if CIDX[h] not in chi: chi[CIDX[h]] = int(np.trace(coords(act(h, Bas))))
    sizes = {}
    for h in G: sizes[CIDX[h]] = sizes.get(CIDX[h], 0) + 1
    dec = {nm: round(sum(sizes[k] * chi[k] * np.conj(row[k]) for k in chi).real / n) for nm, row in zip(names, TABLE)}
    print(f"class {cls}: cup product antisymmetric {np.array_equal(Om, -Om.T)}, |det| = {round(abs(np.linalg.det(Om.astype(float))))};"
          f" H^1 = { {k: v for k, v in dec.items() if v} }")
    rng = np.random.default_rng(cls)
    for nm in ("15", "21"):
        d = int(nm); K = K_of[nm]; row = TABLE[names.index(nm)]; ch = {h: round(row[CIDX[h]].real) for h in G}
        X0 = Bas @ rng.integers(-2, 3, (272, 4))
        RK = sum(act(k, X0) for k in K); Cx = coords(sum(ch[inv(h)] * act(h, RK) for h in G))
        xs = [np.array(r, dtype=np.int64) for r in fz(Cx.T).hnf().tolist() if any(r)]; assert len(xs) == 2
        x1, x2 = [x // math.gcd(*[int(t) for t in x]) for x in xs]
        def actc(h, x): return coords(act(h, Bas @ x.reshape(-1, 1)))[:, 0]
        # the K-fixed plane: saturate span(x1, x2) and read its degree
        M2 = np.stack([x1, x2], 1); g2 = 0
        for i, j in itertools.combinations(range(272), 2):
            m_ = int(M2[i, 0] * M2[j, 1] - M2[j, 0] * M2[i, 1])
            if m_: g2 = math.gcd(g2, m_)
        nK = abs(int(x1 @ Om @ x2)) // g2
        hs, cols = [], []
        for h in [G[i] for i in rng.permutation(n)]:
            v = actc(h, x1)
            if np.linalg.matrix_rank(np.array(cols + [v], float).T) > len(cols): hs.append(h); cols.append(v)
            if len(cols) == d: break
        A = np.array(cols, dtype=object).T; B = np.array([actc(h, x2) for h in hs], dtype=object).T
        Phi = A.T.dot(Om.astype(object)).dot(B); assert np.array_equal(Phi, Phi.T)
        W0 = dual_rows(np.vstack([A, B]).tolist())                                 # L0 = L1 cap L2 in c-coordinates
        A0 = np.array([[int(t) for t in r] for r in A.dot(W0)], dtype=object); B0 = np.array([[int(t) for t in r] for r in B.dot(W0)], dtype=object)
        P0 = W0.T.dot(Phi).dot(W0); P0 = -P0 if P0[0, 0] < 0 else P0
        primes = [p for p in (2, 3, 5, 7) if kernel_dim_mod(np.hstack([A0, B0]).tolist(), p)]
        WH = dual_rows(np.hstack([A, B]).tolist()); DH = 1                          # H, to bound the glue exponent
        for t in WH.flatten(): DH = DH * Fraction(t).denominator // math.gcd(DH, Fraction(t).denominator)
        assert all(DH % (q * q) for q in (2, 3, 5, 7)), "glue exponent not squarefree"
        conds = {p: [(0, None)] + [(1, a_ * A0 + b_ * B0) for a_, b_ in [(1, t) for t in range(p)] + [(0, 1)]] + [(2, np.vstack([A0, B0]))] for p in primes}
        below = []
        for combo in itertools.product(*[conds[p] for p in primes]):
            g = 1; rows = [[Fraction(int(i == j)) for j in range(d)] for i in range(d)]
            for p, (dw, Mc) in zip(primes, combo):
                g *= p ** dw
                if Mc is not None: rows += [[Fraction(int(t), p) for t in r] for r in Mc]
            WS = dual_rows(rows); Gm = WS.T.dot(P0).dot(WS)
            den = 1
            for t in Gm.flatten(): den = den * Fraction(t).denominator // math.gcd(den, Fraction(t).denominator)
            Gr = fz([[Fraction(t) * den for t in r] for r in Gm]).lll(rep="gram").tolist()
            if short_vectors_exist([[int(t) for t in r] for r in Gr], 60 * g * den - 1): below.append(combo)
        print(f"   rho = {nm}: glue primes {primes} (exponent | {DH}); Galois plane (V^K, K of order {len(K)}) has degree {nK}; "
              f"planes of degree < 60: {'none' if not below else len(below)}")
