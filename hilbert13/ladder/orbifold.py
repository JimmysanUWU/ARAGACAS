"""Twisted Laplacians on quotients C/K of the (2,4,7) A7-curve, in exact hyperbolic geometry.

Tiling.  C is tiled by 2*2520 copies of the hyperbolic triangle T with angles pi/2 at A (over the
order-2 branch point), pi/4 at B (order 4), pi/7 at C (order 7): tiles U_g, L_g (g in A7), glued
    U_g ~ L_g       along AB,     U_g ~ L_{g b}   along BC,     U_g ~ L_{g a}   along CA
(right multiplications), and A7 acts on the left: k.U_g = U_{kg}.

Chart.  Every tile, U or L, is identified isometrically with ONE reference triangle, placed in the
Klein model (geodesics = straight lines) and parametrised affinely by R = {(s,t): s,t>=0, s+t<=1}:
    x(s,t) = X_A + s (X_B - X_A) + t (X_C - X_A).
Two tiles sharing a side are mirror images across it, so the side gets the same (s,t)-parameters from
both tiles: all gluings are the identity in R-coordinates.  A uniform subdivision of R therefore gives
a conforming mesh of C (and of every quotient), and in R-coordinates the hyperbolic Dirichlet form and
area are EXACTLY
    int grad u . A_R grad u ds dt,     int u^2 w_R ds dt,
    A_R = |det J| J^-1 A_K(x) J^-T,  w_R = |det J| w_K(x),   J = [X_B - X_A, X_C - X_A],
    A_K(x) = (1-|x|^2)^(-1/2) (I - x x^T),   w_K(x) = (1-|x|^2)^(-3/2)      (Klein metric).

Twisted quotients.  For K <= A7 and a character eps: K -> {+1,-1}, the eps-twisted problem is the
Laplacian on {f on C : f(k x) = eps(k) f(x)}, i.e. on the orbifold K\\C with sign-twisted gluings.
Its spectrum is the union of the rho-isotypic spectra of C over the irreps rho with
<rho|_K, eps> != 0 (every such eigenspace contains eps-vectors for K).
"""
import numpy as np
from math import pi, cos, sin, tan, cosh, sinh, acosh, sqrt
from a7 import A7, mul, inv, generated, E


# ---------------------------------------------------------------- reference triangle (Klein model)
def mink(p, q):
    return p[0] * q[0] - p[1] * q[1] - p[2] * q[2]


def reference_triangle(center="mid_bc"):
    """Klein coordinates of A, B, C after moving the chosen centre to the origin."""
    al, be, ga = pi / 2, pi / 4, pi / 7
    c_ab = acosh(cos(ga) / sin(be))      # side AB (opposite C)
    b_ac = acosh(cos(be) / sin(ga))      # side AC (opposite B)
    A = np.array([1.0, 0.0, 0.0])
    B = np.array([cosh(c_ab), sinh(c_ab), 0.0])
    C = np.array([cosh(b_ac), 0.0, sinh(b_ac)])
    if center == "A":
        P = A
    elif center == "mid_bc":
        P = B + C
    elif center == "centroid":
        P = A + B + C
    else:
        wa, wb, wc = center
        P = wa * A + wb * B + wc * C
    P = P / sqrt(mink(P, P))
    g, p = P[0], P[1:]
    Lam = np.eye(3)
    Lam[0, 0] = g
    Lam[0, 1:] = -p
    Lam[1:, 0] = -p
    Lam[1:, 1:] = np.eye(2) + np.outer(p, p) / (1 + g)
    out = []
    for V in (A, B, C):
        W = Lam @ V
        out.append(W[1:] / W[0])
    return out


def metric_R(XA, XB, XC, s, t):
    """(A_R, w_R) at R-coordinates (s,t) (numpy, vectorised over s,t arrays)."""
    J = np.column_stack([XB - XA, XC - XA])
    Ji = np.linalg.inv(J)
    dJ = abs(np.linalg.det(J))
    x = XA[None, :] + np.outer(s, XB - XA) + np.outer(t, XC - XA)
    r2 = (x ** 2).sum(1)
    AK = (np.eye(2)[None] - x[:, :, None] * x[:, None, :]) / np.sqrt(1 - r2)[:, None, None]
    AR = dJ * np.einsum("ij,njk,lk->nil", Ji, AK, Ji)
    wR = dJ * (1 - r2) ** -1.5
    return AR, wR


# ---------------------------------------------------------------- quotient combinatorics
def quotient_tiles(a, b, Kgens, eps_gens=None):
    """Right cosets K g; returns (reps, glue) where glue lists (coset_U, side, coset_L, sign):
    U_{r} is glued to L_{r'} along `side`, and f|L_{r x} = sign * f|L_{r'} transported."""
    K = generated(Kgens) if Kgens else {E}
    # character on K by propagation from generator values
    eps = {E: 1}
    if Kgens:
        vals = eps_gens or [1] * len(Kgens)
        frontier = [E]
        while frontier:
            new = []
            for x in frontier:
                for g, v in zip(Kgens, vals):
                    y = mul(g, x)
                    e = v * eps[x]
                    if y in eps:
                        assert eps[y] == e, "not a character"
                    else:
                        eps[y] = e
                        new.append(y)
            frontier = new
    rep_of, reps = {}, []
    for g in A7:
        if g in rep_of:
            continue
        coset = [mul(k, g) for k in K]
        r = min(coset)
        for h in coset:
            rep_of[h] = r
        reps.append(r)
    reps.sort()
    cid = {r: i for i, r in enumerate(reps)}

    def locate(g):
        r = rep_of[g]
        k = mul(g, inv(r))
        return cid[r], eps[k]

    glue = []
    for i, r in enumerate(reps):
        glue.append((i, "AB", i, 1))
        j, s = locate(mul(r, b))
        glue.append((i, "BC", j, s))
        j, s = locate(mul(r, a))
        glue.append((i, "CA", j, s))
    return reps, glue, len(K)


class SignedUF:
    def __init__(self):
        self.par, self.sgn = {}, {}

    def find(self, x):
        if x not in self.par:
            self.par[x], self.sgn[x] = x, 1
            return x, 1
        s = 1
        path = []
        while self.par[x] != x:
            path.append(x)
            s *= self.sgn[x]
            x = self.par[x]
        # path compression
        acc = s
        for y in path:
            nxt_sign = self.sgn[y]
            self.par[y], self.sgn[y] = x, acc
            acc *= nxt_sign
        return x, s

    def union(self, x, y, s):
        """impose f(x) = s f(y); returns False if this contradicts (forces zero)."""
        rx, sx = self.find(x)
        ry, sy = self.find(y)
        if rx == ry:
            return sx == s * sy
        self.par[rx], self.sgn[rx] = ry, sx * s * sy
        return True


def side_nodes(n, side):
    if side == "AB":
        return [(i, 0) for i in range(n + 1)]
    if side == "CA":
        return [(0, j) for j in range(n + 1)]
    return [(n - j, j) for j in range(n + 1)]


def ref_elements(n):
    els = []
    for i in range(n):
        for j in range(n - i):
            els.append(((i, j), (i + 1, j), (i, j + 1)))
            if i + j + 2 <= n:
                els.append(((i + 1, j), (i + 1, j + 1), (i, j + 1)))
    return els


def build_dofs(n, ntiles_half, glue, kind):
    """kind 'P1': dofs are mesh vertices; 'CR': dofs are mesh edges.
    Returns dof(tile, localkey) -> (index or None, sign)."""
    uf = SignedUF()
    zero = set()
    for (i, side, j, s) in glue:
        tu, tl = 2 * i, 2 * j + 1
        nodes = side_nodes(n, side)
        if kind == "P1":
            keys = [(v,) for v in nodes]
        else:
            keys = [tuple(sorted((nodes[m], nodes[m + 1]))) for m in range(n)]
        for kk in keys:
            if not uf.union((tu, kk), (tl, kk), s):
                zero.add(uf.find((tu, kk))[0])
    return uf, zero


def assemble(n, reps, glue, kind, elem_mats):
    """elem_mats: list over reference elements of (K_e, M_e) 3x3 arrays (same for every tile).
    Returns sparse K, M over the free dofs."""
    import scipy.sparse as sp
    uf, zero = build_dofs(n, len(reps), glue, kind)
    els = ref_elements(n)
    index = {}
    rows, cols, kv, mv = [], [], [], []
    zero_roots = {uf.find(z)[0] for z in zero}
    for tile in range(2 * len(reps)):
        for e, tri in enumerate(els):
            if kind == "P1":
                keys = [(v,) for v in tri]
            else:  # edge opposite vertex m
                keys = [tuple(sorted((tri[(m + 1) % 3], tri[(m + 2) % 3]))) for m in range(3)]
            ids, sg = [], []
            for kk in keys:
                r, s = uf.find((tile, kk))
                if r in zero_roots:
                    ids.append(-1); sg.append(0)
                    continue
                if r not in index:
                    index[r] = len(index)
                ids.append(index[r]); sg.append(s)
            Ke, Me = elem_mats[e]
            for p in range(3):
                if ids[p] < 0:
                    continue
                for q in range(3):
                    if ids[q] < 0:
                        continue
                    rows.append(ids[p]); cols.append(ids[q])
                    kv.append(sg[p] * sg[q] * Ke[p, q]); mv.append(sg[p] * sg[q] * Me[p, q])
    N = len(index)
    K = sp.csr_matrix((kv, (rows, cols)), shape=(N, N))
    M = sp.csr_matrix((mv, (rows, cols)), shape=(N, N))
    return K, M


# ---------------------------------------------------------------- element matrices (numerics)
# degree-4 Dunavant rule on the unit triangle (barycentric weights sum to 1)
_DUN = [(0.223381589678011, 0.108103018168070, 0.445948490915965, 0.445948490915965),
        (0.223381589678011, 0.445948490915965, 0.108103018168070, 0.445948490915965),
        (0.223381589678011, 0.445948490915965, 0.445948490915965, 0.108103018168070),
        (0.109951743655322, 0.816847572980459, 0.091576213509771, 0.091576213509771),
        (0.109951743655322, 0.091576213509771, 0.816847572980459, 0.091576213509771),
        (0.109951743655322, 0.091576213509771, 0.091576213509771, 0.816847572980459)]


def grads(P):
    """gradients of barycentric coordinates of triangle P (3x2); returns (G 2x3, area)."""
    T = np.array([[P[1][0] - P[0][0], P[2][0] - P[0][0]], [P[1][1] - P[0][1], P[2][1] - P[0][1]]])
    area = abs(np.linalg.det(T)) / 2
    Ti = np.linalg.inv(T)            # rows: grad lambda1, grad lambda2
    G = np.zeros((2, 3))
    G[:, 1], G[:, 2] = Ti[0], Ti[1]
    G[:, 0] = -Ti[0] - Ti[1]
    return G, area


def numeric_elem_mats(n, XA, XB, XC, kind):
    out = []
    for tri in ref_elements(n):
        P = [np.array(v, float) / n for v in tri]
        G, area = grads(P)
        lam = np.array([q[1:] for q in _DUN])
        wq = np.array([q[0] for q in _DUN])
        pts = lam @ np.array(P)
        AR, wR = metric_R(XA, XB, XC, pts[:, 0], pts[:, 1])
        Abar = np.einsum("q,qij->ij", wq, AR)
        if kind == "P1":
            Ke = area * G.T @ Abar @ G
            Me = area * np.einsum("q,q,qi,qj->ij", wq, wR, lam, lam)
        else:
            Ke = 4 * area * G.T @ Abar @ G
            psi = 1 - 2 * lam
            Me = area * np.einsum("q,q,qi,qj->ij", wq, wR, psi, psi)
        out.append((Ke, Me))
    return out


def eigs(K, M, k=8, sigma=-0.01):
    import scipy.sparse.linalg as sla
    vals = sla.eigsh(K.tocsc(), k=k, M=M.tocsc(), sigma=sigma, which="LM", return_eigenvectors=False)
    return np.sort(vals)


def clusters(vals, tol=1e-4):
    cl, cur = [], [vals[0]]
    for v in vals[1:]:
        if abs(v - cur[-1]) < tol * max(1, abs(v)):
            cur.append(v)
        else:
            cl.append(cur); cur = [v]
    cl.append(cur)
    return [(round(float(np.mean(x)), 6), len(x)) for x in cl]
