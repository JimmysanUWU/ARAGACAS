"""Certified lower bound for the lowest eigenvalue of an eps-twisted quotient problem on C/K.

Usage:  python3 certify.py <quotient Q1|Q2> <n> <triple index> [sigma_factor]

Chain of inequalities (2_SPECTRAL.md, section 2.4):
 (1) exact hyperbolic problem, written in R-coordinates:  energy A_R(z), mass w_R(z)          [orbifold.py]
 (2) comparison problem with piecewise-constant coefficients A_e = c_e P_e <= A_R(z) and
     w_e >= w_R(z) on each mesh triangle e (rigorous: interval arithmetic, python-flint arb):
         lambda_k(1) >= lambda_k(2)                                             (min-max)
 (3) Crouzeix-Raviart discretisation of (2), eigenvalue lambda_{k,h}:
         lambda_k(2) >= lambda_{k,h} / (1 + C_h^2 lambda_{k,h})                 (Liu / Carstensen-Gedicke)
     C_h^2 = max_e (w_e / lambda_min(A_e)) * kappa_e^2,
     kappa_e^2 = max_{x in E}|P-x|^2/8 + h_e^2/j_{1,1}^2    (Carstensen-Gedicke-Rim 2012, Lemma 2.2)
 (4) lambda_{1,h} > sigma  <=  K - sigma M positive definite, verified by a floating-point sparse
     Cholesky factorisation of (K - sigma M - c I) with the a-priori rounding bound
         |L L^T - B| <= gamma_{k+1} |L||L^T|,  || |L||L^T| ||_2 <= ||L||_1 ||L||_inf,
     k = max number of nonzeros in a row of L; success with c > that bound + assembly errors
     proves K - sigma M > 0 for the exact comparison matrices.
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from flint import arb, arb_mat, ctx
from orbifold import quotient_tiles, ref_elements, build_dofs, grads
from triples_data import triples
from a7 import cyc


def arb_max(xs):
    """a ball containing the maximum of any choice of points from the balls xs (Python's max() is unsound on
    overlapping balls: it keeps the earlier ball whenever '>' is undecided)."""
    m = None
    for x in xs:
        m = x if m is None else m.max(x)
    return m


ctx.prec = 120
U = 2.0 ** -53          # unit roundoff (binary64, round to nearest)

QUOTIENTS = {
    # K = 3^2:4 < A6 (order 36, index 70), eps = its sign character:  irreps 6, 14a, 14b, 15, 21
    "Q1": ([cyc((0, 1, 2)), cyc((3, 4, 5)), cyc((0, 3, 1, 4), (2, 5))], [1, 1, -1]),
    # K = S4 (order 24, index 105), eps = a sign character:  irreps 10, 10b, 15, 35 (twice)
    "Q2": ([(4, 3, 1, 6, 0, 5, 2), (0, 6, 2, 3, 5, 4, 1)], [1, -1]),
}


# ------------------------------------------------------------------ rigorous helpers
def fdown(x):
    """largest double <= every point of the arb ball x"""
    lo = x.lower()
    f = float(lo)
    while arb(f) > lo:
        f = np.nextafter(f, -np.inf)
    return f


def fup(x):
    hi = x.upper()
    f = float(hi)
    while arb(f) < hi:
        f = np.nextafter(f, np.inf)
    return f


def ref_triangle_arb(pqr=(2, 4, 7)):
    """Klein coordinates (arb) of A, B, C: the triangle with angles pi/p, pi/q, pi/r (p = 2), moved by
    the boost that sends the midpoint of BC to the origin (same placement as orbifold.reference_triangle).
    p = 2 uses the original right-angle formulas (bit-identical to certificate.txt); p > 2 uses the general
    hyperbolic law of cosines for angles, with C placed at angle pi/p from the AB axis
    (same placement as signatures_spectrum.ref_triangle)."""
    pi = arb.pi()
    al, be, ga = pi / pqr[0], pi / pqr[1], pi / pqr[2]
    if pqr[0] == 2:
        c_ab = (ga.cos() / be.sin()).acosh()
        b_ac = (be.cos() / ga.sin()).acosh()
        A = [arb(1), arb(0), arb(0)]
        B = [c_ab.cosh(), c_ab.sinh(), arb(0)]
        C = [b_ac.cosh(), arb(0), b_ac.sinh()]
    else:
        c_ab = ((ga.cos() + al.cos() * be.cos()) / (al.sin() * be.sin())).acosh()
        b_ac = ((be.cos() + al.cos() * ga.cos()) / (al.sin() * ga.sin())).acosh()
        A = [arb(1), arb(0), arb(0)]
        B = [c_ab.cosh(), c_ab.sinh(), arb(0)]
        C = [b_ac.cosh(), b_ac.sinh() * al.cos(), b_ac.sinh() * al.sin()]
    P = [B[i] + C[i] for i in range(3)]
    nP = (P[0] ** 2 - P[1] ** 2 - P[2] ** 2).sqrt()
    P = [p / nP for p in P]
    g, p1, p2 = P
    def boost(V):
        # inverse of the boost e0 -> P
        T = g * V[0] - p1 * V[1] - p2 * V[2]
        dot = p1 * V[1] + p2 * V[2]
        X1 = -p1 * V[0] + V[1] + p1 * dot / (1 + g)
        X2 = -p2 * V[0] + V[2] + p2 * dot / (1 + g)
        return (X1 / T, X2 / T)
    return [boost(V) for V in (A, B, C)]


def metric_box(XA, XB, XC, s_lo, s_hi, t_lo, t_hi):
    """arb enclosures of (A_R entries, w_R) over the box [s_lo,s_hi] x [t_lo,t_hi]."""
    s = arb((s_lo + s_hi) / 2, (s_hi - s_lo) / 2)
    t = arb((t_lo + t_hi) / 2, (t_hi - t_lo) / 2)
    J = [[XB[0] - XA[0], XC[0] - XA[0]], [XB[1] - XA[1], XC[1] - XA[1]]]
    dJ = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    dJ = abs(dJ)
    x1 = XA[0] + s * J[0][0] + t * J[0][1]
    x2 = XA[1] + s * J[1][0] + t * J[1][1]
    r2 = x1 ** 2 + x2 ** 2
    q = 1 - r2
    # M = I - x x^T ; A_R = dJ q^(-1/2) J^-1 M J^-T ;  J^-1 = adj(J)/det
    m11, m12, m22 = 1 - x1 ** 2, -x1 * x2, 1 - x2 ** 2
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    Ji = [[J[1][1] / det, -J[0][1] / det], [-J[1][0] / det, J[0][0] / det]]
    M = [[m11, m12], [m12, m22]]
    JM = [[sum(Ji[a][c] * M[c][d] for c in range(2)) for d in range(2)] for a in range(2)]
    S = [[sum(JM[a][d] * Ji[b][d] for d in range(2)) for b in range(2)] for a in range(2)]
    f = dJ / q.sqrt()
    AR = [[f * S[0][0], f * S[0][1]], [f * S[1][0], f * S[1][1]]]
    wR = dJ / (q * q.sqrt())
    return AR, wR


def lam_min_pencil_lower(Q, P):
    """lower bound of lambda_min(P^-1 Q) for every symmetric Q in the arb enclosure (P exact floats, SPD).

    With P = L L^T (Cholesky, in arb), lambda(P^-1 Q) = lambda(W), W = L^-1 Q L^-T symmetric.
    Weyl: lambda_min(W) >= lambda_min(W_mid) - ||W - W_mid||_2, and ||.||_2 <= max row sum of radii.
    (Evaluating the characteristic polynomial on intervals instead would lose sqrt(width).)"""
    P11, P12, P22 = arb(P[0][0]), arb(P[0][1]), arb(P[1][1])
    l11 = P11.sqrt()
    l21 = P12 / l11
    l22 = (P22 - l21 * l21).sqrt()
    # L^-1 = [[1/l11, 0], [-l21/(l11 l22), 1/l22]]
    Li = [[1 / l11, arb(0)], [-l21 / (l11 * l22), 1 / l22]]
    Qs = [[Q[0][0], Q[0][1]], [Q[0][1], Q[1][1]]]
    LQ = [[sum(Li[a][c] * Qs[c][d] for c in range(2)) for d in range(2)] for a in range(2)]
    W = [[sum(LQ[a][d] * Li[b][d] for d in range(2)) for b in range(2)] for a in range(2)]
    w11, w22 = W[0][0], W[1][1]
    w12 = (W[0][1] + W[1][0]) / 2          # W is symmetric; both enclosures contain the true value
    m11, m12, m22 = arb(w11.mid()), arb(w12.mid()), arb(w22.mid())
    r11, r12, r22 = w11.rad(), w12.rad(), w22.rad()
    pert = arb(r11) + arb(r12) + arb(r22)      # >= both row sums of the radius matrix
    tr = m11 + m22
    disc = (m11 - m22) ** 2 + 4 * m12 * m12      # exact up to arb rounding: no interval blow-up
    lam_mid = (tr - disc.sqrt()) / 2
    return (lam_mid - pert).lower()


def element_data(n, XA, XB, XC, sub=3):
    """for each reference element: (P_e float 2x2, c_e, w_e, lamminA_e) with rigorous rounding."""
    import orbifold
    Xf = [np.array([float(x.mid()) for x in X]) for X in (XA, XB, XC)]
    out = []
    for tri in ref_elements(n):
        up = (tri[1][0] == tri[0][0] + 1 and tri[1][1] == tri[0][1])
        i0, j0 = min(v[0] for v in tri), min(v[1] for v in tri)
        cen = np.mean(np.array(tri, float), axis=0) / n
        ARc, _ = orbifold.metric_R(*Xf, np.array([cen[0]]), np.array([cen[1]]))
        P = ARc[0]
        P = [[float(P[0, 0]), float((P[0, 1] + P[1, 0]) / 2)], [float((P[0, 1] + P[1, 0]) / 2), float(P[1, 1])]]
        c_lo, w_hi = None, None
        for p in range(sub):
            for q in range(sub):
                if (up and p + q > sub - 1) or ((not up) and p + q < sub - 1):
                    continue
                s_lo = arb(i0 * sub + p) / (n * sub)
                s_hi = arb(i0 * sub + p + 1) / (n * sub)
                t_lo = arb(j0 * sub + q) / (n * sub)
                t_hi = arb(j0 * sub + q + 1) / (n * sub)
                AR, wR = metric_box(XA, XB, XC, s_lo, s_hi, t_lo, t_hi)
                lo = lam_min_pencil_lower(AR, P)
                hi = wR.upper()
                c_lo = lo if c_lo is None or lo < c_lo else c_lo
                w_hi = hi if w_hi is None or hi > w_hi else w_hi
        c_e = fdown(c_lo)
        w_e = fup(w_hi)
        # lambda_min(P) lower bound
        lp = lam_min_pencil_lower([[arb(P[0][0]), arb(P[0][1])], [arb(P[1][0]), arb(P[1][1])]],
                                  [[1.0, 0.0], [0.0, 1.0]])
        out.append((P, c_e, w_e, lp * arb(c_e)))
    return out


def elem_mats_rigorous(n, edata):
    """CR element matrices from (c_e P_e, w_e): float midpoints and entrywise error bounds."""
    mats = []
    for tri, (P, c_e, w_e, _) in zip(ref_elements(n), edata):
        Pts = [np.array(v, float) for v in tri]          # in units 1/n: G = n * Z
        Z, _ = grads(Pts)                                 # exact small integers here
        Z = np.rint(Z)
        Ke = np.zeros((3, 3)); Kerr = np.zeros((3, 3))
        for a in range(3):
            for b in range(3):
                v = 2 * arb(c_e) * sum(arb(Z[r, a]) * arb(P[r][s]) * arb(Z[s, b]) for r in range(2) for s in range(2))
                Ke[a, b] = float(v.mid())
                Kerr[a, b] = fup(abs(v - arb(Ke[a, b])))
        mv = arb(w_e) / (6 * n * n)
        m = float(mv.mid())
        merr = fup(abs(mv - arb(m)))
        mats.append((Ke, Kerr, m, merr))
    return mats


def assemble_rigorous(n, reps, glue, mats):
    uf, zero = build_dofs(n, len(reps), glue, "CR")
    assert not zero, "a CR dof was identified with its own negative"
    els = ref_elements(n)
    index = {}
    R, Cc, V, Vabs, Verr = [], [], [], [], []
    Md = {}
    for tile in range(2 * len(reps)):
        for e, tri in enumerate(els):
            keys = [tuple(sorted((tri[(m + 1) % 3], tri[(m + 2) % 3]))) for m in range(3)]
            ids, sg = [], []
            for kk in keys:
                r, s = uf.find((tile, kk))
                if r not in index:
                    index[r] = len(index)
                ids.append(index[r]); sg.append(s)
            Ke, Kerr, m, merr = mats[e]
            for p in range(3):
                Md.setdefault(ids[p], []).append((m, merr))
                for q in range(3):
                    R.append(ids[p]); Cc.append(ids[q])
                    V.append(sg[p] * sg[q] * Ke[p, q]); Vabs.append(abs(Ke[p, q])); Verr.append(Kerr[p, q])
    N = len(index)
    K = sp.csr_matrix((V, (R, Cc)), shape=(N, N))
    Kabs = sp.csr_matrix((Vabs, (R, Cc)), shape=(N, N))
    Kerr = sp.csr_matrix((Verr, (R, Cc)), shape=(N, N))
    cnt = sp.csr_matrix((np.ones(len(V)), (R, Cc)), shape=(N, N))
    Mdiag = np.zeros(N); Merr = np.zeros(N); Mabs = np.zeros(N); Mcnt = np.zeros(N)
    for i, lst in Md.items():
        Mdiag[i] = sum(x for x, _ in lst)            # floating sum (a few terms)
        Merr[i] = sum(y for _, y in lst)
        Mabs[i] = sum(abs(x) for x, _ in lst)
        Mcnt[i] = len(lst)
    return K, Kabs, Kerr, cnt, Mdiag, Merr, Mabs, Mcnt


def assemble_rigorous_fast(n, reps, glue, mats):
    """Vectorised version of assemble_rigorous (same matrices up to a permutation of the dofs)."""
    from orbifold import SignedUF, side_nodes
    els = ref_elements(n)
    key = lambda p, q: (p, q) if p < q else (q, p)
    eid = {}
    elem_edges = np.zeros((len(els), 3), dtype=np.int64)
    for e, tri in enumerate(els):
        for m in range(3):
            k = key(tri[(m + 1) % 3], tri[(m + 2) % 3])
            if k not in eid:
                eid[k] = len(eid)
            elem_edges[e, m] = eid[k]
    E = len(eid)
    bnd = {}
    for side in ("AB", "BC", "CA"):
        nodes = side_nodes(n, side)
        for m in range(n):
            bnd[(side, m)] = eid[key(nodes[m], nodes[m + 1])]
    is_b = np.zeros(E, bool)
    is_b[list(bnd.values())] = True
    interior = np.flatnonzero(~is_b)
    T = 2 * len(reps)
    uf = SignedUF()
    for (i, side, j, s) in glue:
        for m in range(n):
            l = bnd[(side, m)]
            assert uf.union((2 * i, l), (2 * j + 1, l), s), "a CR dof was identified with its own negative"
    dof = np.zeros((T, E), dtype=np.int64)
    sgn = np.ones((T, E))
    dof[:, interior] = np.arange(T)[:, None] * len(interior) + np.arange(len(interior))[None, :]
    root_id = {}
    base = T * len(interior)
    bl = np.flatnonzero(is_b)
    for t in range(T):
        for l in bl:
            r, s = uf.find((t, int(l)))
            if r not in root_id:
                root_id[r] = base + len(root_id)
            dof[t, l] = root_id[r]
            sgn[t, l] = s
    N = base + len(root_id)
    Ke = np.array([m[0] for m in mats]); Kerr_e = np.array([m[1] for m in mats])
    me = np.array([m[2] for m in mats]); merr_e = np.array([m[3] for m in mats])
    D = dof[:, elem_edges]                      # (T, nel, 3)
    S = sgn[:, elem_edges]
    rows = np.broadcast_to(D[:, :, :, None], D.shape + (3,)).ravel()
    cols = np.broadcast_to(D[:, :, None, :], D.shape + (3,)).ravel()
    ss = (S[:, :, :, None] * S[:, :, None, :])
    vals = (ss * Ke[None]).ravel()
    vabs = np.broadcast_to(np.abs(Ke)[None], ss.shape).ravel()
    verr = np.broadcast_to(Kerr_e[None], ss.shape).ravel()
    K = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    Kabs = sp.csr_matrix((vabs, (rows, cols)), shape=(N, N))
    Kerr = sp.csr_matrix((verr, (rows, cols)), shape=(N, N))
    cnt = sp.csr_matrix((np.ones(len(vals)), (rows, cols)), shape=(N, N))
    Dm = D.ravel()
    mm = np.broadcast_to(me[None, :, None], D.shape).ravel()
    mr = np.broadcast_to(merr_e[None, :, None], D.shape).ravel()
    Mdiag = np.bincount(Dm, weights=mm, minlength=N)
    Merr = np.bincount(Dm, weights=mr, minlength=N)
    Mabs = np.bincount(Dm, weights=np.abs(mm), minlength=N)
    Mcnt = np.bincount(Dm, minlength=N).astype(float)
    return K, Kabs, Kerr, cnt, Mdiag, Merr, Mabs, Mcnt


def gamma(k):
    return k * U / (1 - k * U)


def certify_trivial(n, tri, sigma=1.0, verbose=True, ab=None, pqr=(2, 4, 7)):
    """Q0: A7-invariant functions (the (2,4,7) orbifold itself, 2 tiles).  Certifies mu_2(Q0) > bound.
    Constants are deflated densely: B + alpha z z^T > 0 with z = M 1, alpha > sigma / (1^T M 1), implies
    K - sigma M > 0 on the M-orthogonal complement of the constants, i.e. lambda_{2,h} > sigma."""
    a, b = ab if ab is not None else triples[tri][:2]
    reps, glue, oK = quotient_tiles(a, b, [a, b], [1, 1])
    XA, XB, XC = ref_triangle_arb(pqr)
    ed = element_data(n, XA, XB, XC)
    j11_lo = arb("3.8317059702075")
    kappa2 = (arb(1) / 8 + arb(2) / (j11_lo * j11_lo)) / (n * n)
    Ch2 = kappa2 * arb_max(arb(e[2]) / e[3] for e in ed)
    mats = elem_mats_rigorous(n, ed)
    K, Kabs, Kerr, cnt, Md, Merr, Mabs, Mcnt = assemble_rigorous(n, reps, glue, mats)
    N = K.shape[0]
    Kd = K.toarray()
    m = float(Md.sum())
    alpha = 2.0 ** int(np.ceil(np.log2(2 * sigma / m)))          # power of two, > sigma/m
    zz = alpha * np.outer(Md, Md)
    Bd = Kd - np.diag(sigma * Md) + zz
    # error bound: assembly errors of K, M; error of z (Merr, summation); rounding of sigma*M, alpha z z^T and sums
    errK = (Kerr + gamma(int(cnt.max())) * Kabs).toarray()
    errM = Merr + gamma(int(Mcnt.max())) * Mabs
    errz = alpha * (np.outer(errM, np.abs(Md) + errM) + np.outer(np.abs(Md), errM)) + 3 * U * np.abs(zz)
    errB = errK + errz + np.diag(sigma * errM + U * sigma * Md) + 2 * U * (np.abs(Kd) + np.abs(zz) + np.abs(Bd))
    eta = float(np.abs(errB).sum(axis=1).max()) * (1 + 1e-6)
    c = 4 * eta
    for attempt in range(8):
        L = np.linalg.cholesky(Bd - c * np.eye(N))      # raises LinAlgError if not (numerically) PD
        L1 = float(np.abs(L).sum(axis=0).max()) * (1 + 2 * gamma(N))
        Linf = float(np.abs(L).sum(axis=1).max()) * (1 + 2 * gamma(N))
        need = gamma(N + 1) * L1 * Linf + U * (float(np.abs(np.diag(Bd)).max()) + c) + eta
        if c > need * 1.01:
            break
        c = 2 * need
    ok = c > need * 1.01
    s = arb(sigma)
    bound = s / (1 + Ch2 * s)
    lam = np.sort(np.linalg.eigvalsh(np.diag(Md ** -0.5) @ Kd @ np.diag(Md ** -0.5)))[:3]
    res = dict(ok=ok, which="Q0", n=n, triple=tri, K_order=oK, dofs=N, lam_h=lam.tolist(), sigma=sigma,
               Ch2=fup(Ch2), eta=eta, chol_shift=c, chol_need=need, bound=fdown(bound) if ok else None)
    if verbose:
        print(res, flush=True)
    return res


def certify(which, n, tri, sigma_factor=0.995, verbose=True, sigma=None, fast=True):
    """sigma=None: sigma = sigma_factor * (numerical lambda_{1,h}); otherwise use the given sigma
    (no eigenvalue solve; the Cholesky either certifies lambda_{1,h} > sigma or fails)."""
    if which == "Q0":
        return certify_trivial(n, tri, verbose=verbose)
    a, b, c = triples[tri]
    Kg, ev = QUOTIENTS[which]
    reps, glue, oK = quotient_tiles(a, b, Kg, ev)
    res = certify_tiles(reps, glue, n, sigma_factor=sigma_factor, sigma=sigma, fast=fast)
    res.update(which=which, triple=tri, K_order=oK)
    if verbose:
        print(res, flush=True)
    return res


def certify_tiles(reps, glue, n, pqr=(2, 4, 7), sigma_factor=0.995, sigma=None, fast=True):
    """Certified lower bound for the lowest eigenvalue of a sign-twisted tile complex without constant
    functions (every coset tile pair U_r, L_r is a copy of the (pi/p, pi/q, pi/r) triangle)."""
    import cvxopt, cvxopt.cholmod as chol
    XA, XB, XC = ref_triangle_arb(pqr)
    ed = element_data(n, XA, XB, XC)
    cmin = min(e[1] for e in ed)
    # C_h^2 = kappa^2 max_e w_e / lambda_min(A_e);  kappa^2 = (1/n)^2/8 + (sqrt2/n)^2/j11^2
    j11_lo = arb("3.8317059702075")          # j_{1,1} = 3.83170597020751231...
    kappa2 = (arb(1) / 8 + arb(2) / (j11_lo * j11_lo)) / (n * n)
    ratio = arb_max(arb(e[2]) / e[3] for e in ed)
    Ch2 = kappa2 * ratio
    mats = elem_mats_rigorous(n, ed)
    asm = assemble_rigorous_fast if fast else assemble_rigorous
    K, Kabs, Kerr, cnt, Md, Merr, Mabs, Mcnt = asm(n, reps, glue, mats)
    N = K.shape[0]
    if sigma is None:
        # numerical lowest eigenvalue of the comparison CR problem (only used to choose sigma)
        Mm = sp.diags(Md)
        lam = sla.eigsh(K.tocsc(), k=3, M=Mm.tocsc(), sigma=-0.01, which="LM", return_eigenvectors=False)
        lam = np.sort(lam)
        sigma = float(lam[0] * sigma_factor)
    else:
        lam = np.array([np.nan])
    # exact-vs-float error of B = K - sigma M  (Gershgorin bound on the spectral norm)
    # K entries: arb conversion errors + summation errors gamma_{cnt} * sum|terms|
    cntmax = int(cnt.max())
    errK = Kerr + gamma(cntmax) * Kabs
    errM = Merr + gamma(int(Mcnt.max())) * Mabs
    B = (K - sp.diags(sigma * Md)).tocsc()
    # rounding of sigma*M_ii and of the diagonal subtraction
    Bd = K.diagonal() - sigma * Md
    errB_diag = sigma * errM + U * sigma * Md + U * np.abs(Bd) + U * np.abs(K.diagonal())
    eta = float((np.asarray(abs(errK).sum(axis=1)).ravel() + errB_diag).max()) * (1 + 1e-6)
    # Cholesky of B - c I
    chol.options["supernodal"] = 2
    c = 4 * eta
    for attempt in range(8):
        Bc = (B - sp.eye(N) * c).tocoo()
        S = cvxopt.spmatrix(Bc.data.tolist(), Bc.row.tolist(), Bc.col.tolist(), (N, N))
        try:
            F = chol.symbolic(S)
            chol.numeric(S, F)
        except ArithmeticError:
            return dict(ok=False, reason="Cholesky failed", n=n, sigma=sigma, lam_h=lam.tolist(), c=c, dofs=N)
        L = chol.getfactor(F)
        Lc = sp.csc_matrix((np.array(L.V).ravel(), np.array(L.I).ravel(), np.array(L.CCS[0]).ravel()), shape=(N, N))
        Lr = Lc.tocsr()
        k = int(np.diff(Lr.indptr).max())
        g = gamma(k + 1)
        # || |L||L^T| ||_2 <= || |L| ||_2^2 <= ||L||_1 ||L||_inf  (floating sums of <= N positive terms,
        # inflated by (1 + 2 gamma_N) to cover their rounding)
        L1 = float(np.asarray(abs(Lc).sum(axis=0)).max()) * (1 + 2 * gamma(N))
        Linf = float(np.asarray(abs(Lr).sum(axis=1)).max()) * (1 + 2 * gamma(N))
        need = g * L1 * Linf + U * (float(np.abs(Bd).max()) + c) + eta
        if c > need * 1.01:
            break
        c = 2 * need
    ok = c > need * 1.01
    # final bound  sigma / (1 + Ch2 sigma), rounded down
    s = arb(sigma)
    bound = s / (1 + Ch2 * s)
    return dict(ok=ok, n=n, dofs=N, lam_h=lam.tolist(), sigma=sigma,
                c_min=cmin, Ch2=fup(Ch2), eta=eta, chol_shift=c, chol_need=need, L_rowmax=k, L1=L1, Linf=Linf,
                bound=fdown(bound) if ok else None)


if __name__ == "__main__":
    which, n, tri = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    sf = float(sys.argv[4]) if len(sys.argv) > 4 else 0.995
    certify(which, n, tri, sf)
