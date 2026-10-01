"""Certified eigenvalue inputs of topological Hersch (FRAMEWORK_TOPOLOGICAL_HERSCH.md section 5, VERIFICATION.md).

For a (2,4,7) triple class t this proves, on top of certificate.txt (lambda_1 >= 0.34089, Q2 >= 0.689, Q0 >= 0.997):

 (A) COUNT on Q1 = (3^2:4, sgn).  The CR comparison pencil (K, M) of certify.py has at most ONE eigenvalue below
     sigma.  Proof: a floating LDL^T factorisation (CHOLMOD simplicial, no pivoting, fill-reducing symmetric
     permutation P) of  B - cI,  B = K - sigma M,  and the backward error bound for symmetric Gaussian elimination
         P(B - cI)P^T + Delta = L D L^T,   |Delta| <= gamma_{k+3} |L||D||L^T|          (k = max nnz per row of L)
     (derivation in VERIFICATION.md).  If c > ||Delta||_2 + assembly/rounding errors, then
     B_exact = L D L^T + (positive definite), so #negative eigenvalues of B_exact <= #negative pivots of D (Sylvester
     + Weyl).  One negative pivot gives lambda_{2,h} >= sigma, and Liu / Carstensen-Gedicke for k = 2 gives
         lambda_2(Q1) >= sigma / (1 + C_h^2 sigma).
     Since Q1 sees 6, 14a, 14b, 15, 21 each exactly once (verify_exact.py), at most one of these isotypes has an
     eigenvalue below that bound, and only once.

 (B) UPPER BOUND on C/S5 (S5 = stabiliser of {5,6}; Ind 1 = 1 + 6 + 14a, multiplicity-free).  Rayleigh-Ritz with
     span{1, v}, v a P1 function, for the UPPER comparison problem (A_e^max >= A_R and w_e^min <= w_R on each cell,
     ball arithmetic): every Rayleigh quotient of the true problem is <= that of the comparison problem, so
         lambda_2(C/S5) <= lambda_max of the 2x2 comparison pencil  =: U.
     The trivial isotype has no nonzero eigenvalue below 0.997 (Q0), so 6 or 14a has an eigenvalue <= U.

 (C) LOWER BOUND on C/A6 (Ind 1 = 1 + 6): the second CR comparison eigenvalue, constants deflated densely as in
     certify.certify_trivial.  It excludes the isotype 6 below the bound.

 Conclusion for the class (with certificate.txt):  E_1 is ONE copy of 14a = 14_(5,2), lambda_1 in [0.34089, U], and
 every other eigenvalue of C is >= min(bound (A), Q2, Q0 second).

Usage: python3 certify_th.py count n t sigma | s5 n t | a6 n t sigma
"""
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from flint import arb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import certify as cf
from certify import (U, gamma, fdown, fup, ref_triangle_arb, element_data, elem_mats_rigorous, assemble_rigorous,
                     assemble_rigorous_fast, metric_box, QUOTIENTS)
from orbifold import quotient_tiles, build_dofs, ref_elements, grads, numeric_elem_mats, reference_triangle
from triples_data import triples
from a7 import cyc

S5_GENS = ([cyc((0, 1, 2)), cyc((0, 1, 2, 3, 4)), cyc((0, 1), (5, 6))], [1, 1, 1])
A6_GENS = ([cyc((0, 1, 2)), cyc((1, 2, 3, 4, 5))], [1, 1])


def ch2_bound(n, ed):
    j11_lo = arb("3.8317059702075")
    kappa2 = (arb(1) / 8 + arb(2) / (j11_lo * j11_lo)) / (n * n)
    return kappa2 * max(arb(e[2]) / e[3] for e in ed)


# ------------------------------------------------------------------------------------------- (A) count
def count_below(reps, glue, n, sigma, pqr=(2, 4, 7)):
    import cvxopt
    import cvxopt.cholmod as chol
    XA, XB, XC = ref_triangle_arb(pqr)
    ed = element_data(n, XA, XB, XC)
    Ch2 = ch2_bound(n, ed)
    mats = elem_mats_rigorous(n, ed)
    K, Kabs, Kerr, cnt, Md, Merr, Mabs, Mcnt = assemble_rigorous_fast(n, reps, glue, mats)
    N = K.shape[0]
    errK = Kerr + gamma(int(cnt.max())) * Kabs
    errM = Merr + gamma(int(Mcnt.max())) * Mabs
    B = (K - sp.diags(sigma * Md)).tocsc()
    Bd = K.diagonal() - sigma * Md
    errB_diag = sigma * errM + U * sigma * Md + U * np.abs(Bd) + U * np.abs(K.diagonal())
    eta = float((np.asarray(abs(errK).sum(axis=1)).ravel() + errB_diag).max()) * (1 + 1e-6)
    chol.options["supernodal"] = 0                 # simplicial LDL^T (up-looking), works for indefinite matrices
    c = 4 * eta
    t0 = time.time()
    for attempt in range(8):
        Bc = (B - sp.eye(N) * c).tocoo()
        S = cvxopt.spmatrix(Bc.data.tolist(), Bc.row.tolist(), Bc.col.tolist(), (N, N))
        F = chol.symbolic(S)
        chol.numeric(S, F)                         # raises ArithmeticError only on an exactly zero pivot
        Lf = chol.getfactor(F)                     # strictly lower part = L (unit diagonal implied), diagonal = D
        Lc = sp.csc_matrix((np.array(Lf.V).ravel(), np.array(Lf.I).ravel(), np.array(Lf.CCS[0]).ravel()), shape=(N, N))
        D = Lc.diagonal().copy()
        Lu = (sp.tril(Lc, k=-1) + sp.eye(N)).tocsc()
        Lr = Lu.tocsr()
        k = int(np.diff(Lr.indptr).max())
        L1 = float(np.asarray(abs(Lu).sum(axis=0)).max()) * (1 + 2 * gamma(N))
        Linf = float(np.asarray(abs(Lr).sum(axis=1)).max()) * (1 + 2 * gamma(N))
        Dmax = float(np.abs(D).max())
        need = gamma(k + 3) * Dmax * L1 * Linf + U * (float(np.abs(Bd).max()) + c) + eta
        if c > need * 1.01:
            break
        c = 2 * need
    ok = c > need * 1.01
    neg = int((D < 0).sum())
    s = arb(sigma)
    bound = s / (1 + Ch2 * s)
    return dict(ok=ok, n=n, dofs=N, sigma=sigma, negative_pivots=neg, Ch2=fup(Ch2), eta=eta, shift=c, need=need,
                L_rowmax=k, L1=L1, Linf=Linf, Dmax=Dmax, min_abs_pivot=float(np.abs(D).min()),
                lambda2_bound=fdown(bound) if (ok and neg <= 1) else None, secs=round(time.time() - t0))


# ------------------------------------------------------------------------------------- (B) upper bound
def lam_max_pencil_upper(Q, P):
    """upper bound of lambda_max(P^-1 Q) for every symmetric Q in the arb enclosure (P floats, SPD)."""
    P11, P12, P22 = arb(P[0][0]), arb(P[0][1]), arb(P[1][1])
    l11 = P11.sqrt(); l21 = P12 / l11; l22 = (P22 - l21 * l21).sqrt()
    Li = [[1 / l11, arb(0)], [-l21 / (l11 * l22), 1 / l22]]
    Qs = [[Q[0][0], Q[0][1]], [Q[0][1], Q[1][1]]]
    LQ = [[sum(Li[a][c] * Qs[c][d] for c in range(2)) for d in range(2)] for a in range(2)]
    W = [[sum(LQ[a][d] * Li[b][d] for d in range(2)) for b in range(2)] for a in range(2)]
    w11, w22, w12 = W[0][0], W[1][1], (W[0][1] + W[1][0]) / 2
    m11, m12, m22 = arb(w11.mid()), arb(w12.mid()), arb(w22.mid())
    pert = arb(w11.rad()) + arb(w12.rad()) + arb(w22.rad())
    lam_mid = (m11 + m22 + ((m11 - m22) ** 2 + 4 * m12 * m12).sqrt()) / 2
    return (lam_mid + pert).upper()


def element_data_upper(n, XA, XB, XC, sub=3):
    """per reference element: (P_e float 2x2, c_e^max, w_e^min) with A_R <= c_e^max P_e and w_R >= w_e^min on e."""
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
        c_hi, w_lo = None, None
        for p in range(sub):
            for q in range(sub):
                if (up and p + q > sub - 1) or ((not up) and p + q < sub - 1):
                    continue
                s_lo = arb(i0 * sub + p) / (n * sub); s_hi = arb(i0 * sub + p + 1) / (n * sub)
                t_lo = arb(j0 * sub + q) / (n * sub); t_hi = arb(j0 * sub + q + 1) / (n * sub)
                AR, wR = metric_box(XA, XB, XC, s_lo, s_hi, t_lo, t_hi)
                hi = lam_max_pencil_upper(AR, P)
                lo = wR.lower()
                c_hi = hi if c_hi is None or hi > c_hi else c_hi
                w_lo = lo if w_lo is None or lo < w_lo else w_lo
        out.append((P, arb(c_hi), arb(w_lo)))
    return out


def s5_upper(n, t, pqr=(2, 4, 7)):
    a, b, _ = triples[t]
    reps, glue, oK = quotient_tiles(a, b, *S5_GENS)
    uf, zero = build_dofs(n, len(reps), glue, "P1")
    assert not zero
    els = ref_elements(n)
    index = {}
    ids = np.zeros((2 * len(reps), len(els), 3), dtype=np.int64)
    for tile in range(2 * len(reps)):
        for e, trg in enumerate(els):
            for p, v in enumerate(trg):
                r, s = uf.find((tile, (v,)))
                assert s == 1
                ids[tile, e, p] = index.setdefault(r, len(index))
    N = len(index)
    # numerical trial function v: second eigenvector of the floating P1 problem
    XAf, XBf, XCf = reference_triangle()
    mats = numeric_elem_mats(n, XAf, XBf, XCf, "P1")
    R, Cc, KV, MV = [], [], [], []
    for tile in range(2 * len(reps)):
        for e in range(len(els)):
            Ke, Me = mats[e]
            for p in range(3):
                for q in range(3):
                    R.append(ids[tile, e, p]); Cc.append(ids[tile, e, q]); KV.append(Ke[p, q]); MV.append(Me[p, q])
    Kf = sp.csc_matrix((KV, (R, Cc)), shape=(N, N)); Mf = sp.csc_matrix((MV, (R, Cc)), shape=(N, N))
    lam, V = sla.eigsh(Kf, k=3, M=Mf, sigma=-1e-3, which="LM")
    o = np.argsort(lam); lam, V = lam[o], V[:, o]
    v = V[:, 1] / np.abs(V[:, 1]).max()
    # rigorous comparison quadratic forms on span{1, v}
    XA, XB, XC = ref_triangle_arb(pqr)
    edu = element_data_upper(n, XA, XB, XC)
    area = arb(1) / (2 * n * n)
    kvv, m11, m1v, mvv = arb(0), arb(0), arb(0), arb(0)
    for e, trg in enumerate(els):
        Pts = [np.array(q, float) for q in trg]
        Z, _ = grads(Pts)
        Z = np.rint(Z)                                   # gradients in units of n (small integers)
        P, chi, wlo = edu[e]
        vl = v[ids[:, e, :]]                             # (tiles, 3) exact floats
        gs = [sum(arb(Z[r, a]) * arb(float(vl[tt, a])) for a in range(3)) * n for r in range(2)
              for tt in range(vl.shape[0])]
        ntl = vl.shape[0]
        ke = arb(0)
        for tt in range(ntl):
            g0, g1 = gs[tt], gs[ntl + tt]
            ke += g0 * g0 * arb(P[0][0]) + 2 * g0 * g1 * arb(P[0][1]) + g1 * g1 * arb(P[1][1])
        kvv += chi * area * ke
        s1 = arb(0); s2 = arb(0)
        for tt in range(ntl):
            a0, a1, a2 = arb(float(vl[tt, 0])), arb(float(vl[tt, 1])), arb(float(vl[tt, 2]))
            ssum = a0 + a1 + a2
            s1 += ssum
            s2 += a0 * a0 + a1 * a1 + a2 * a2 + ssum * ssum
        m11 += wlo * area * ntl
        m1v += wlo * area * s1 / 3
        mvv += wlo * area * s2 / 12
    det = m11 * mvv - m1v * m1v
    assert det > 0
    Ub = kvv * m11 / det
    return dict(n=n, triple=t, K_order=oK, dofs=N, lam_h_numerical=lam.tolist(), U=fup(Ub))


# ---------------------------------------------------------------------------------------- (C) A6 lower
def a6_lower(n, t, sigma, pqr=(2, 4, 7)):
    """second eigenvalue of the CR comparison problem on C/A6 > sigma, constants deflated (as certify_trivial)."""
    a, b, _ = triples[t]
    reps, glue, oK = quotient_tiles(a, b, *A6_GENS)
    XA, XB, XC = ref_triangle_arb(pqr)
    ed = element_data(n, XA, XB, XC)
    Ch2 = ch2_bound(n, ed)
    mats = elem_mats_rigorous(n, ed)
    K, Kabs, Kerr, cnt, Md, Merr, Mabs, Mcnt = assemble_rigorous(n, reps, glue, mats)
    N = K.shape[0]
    Kd = K.toarray()
    m = float(Md.sum())
    alpha = 2.0 ** int(np.ceil(np.log2(2 * sigma / m)))
    zz = alpha * np.outer(Md, Md)
    Bd = Kd - np.diag(sigma * Md) + zz
    errK = (Kerr + gamma(int(cnt.max())) * Kabs).toarray()
    errM = Merr + gamma(int(Mcnt.max())) * Mabs
    errz = alpha * (np.outer(errM, np.abs(Md) + errM) + np.outer(np.abs(Md), errM)) + 3 * U * np.abs(zz)
    errB = errK + errz + np.diag(sigma * errM + U * sigma * Md) + 2 * U * (np.abs(Kd) + np.abs(zz) + np.abs(Bd))
    eta = float(np.abs(errB).sum(axis=1).max()) * (1 + 1e-6)
    c = 4 * eta
    for attempt in range(8):
        L = np.linalg.cholesky(Bd - c * np.eye(N))
        L1 = float(np.abs(L).sum(axis=0).max()) * (1 + 2 * gamma(N))
        Linf = float(np.abs(L).sum(axis=1).max()) * (1 + 2 * gamma(N))
        need = gamma(N + 1) * L1 * Linf + U * (float(np.abs(np.diag(Bd)).max()) + c) + eta
        if c > need * 1.01:
            break
        c = 2 * need
    ok = c > need * 1.01
    s = arb(sigma)
    lam = np.sort(np.linalg.eigvalsh(np.diag(Md ** -0.5) @ Kd @ np.diag(Md ** -0.5)))[:3]
    return dict(ok=ok, n=n, triple=t, K_order=oK, dofs=N, lam_h=lam.tolist(), sigma=sigma, Ch2=fup(Ch2),
                bound=fdown(s / (1 + Ch2 * s)) if ok else None)


if __name__ == "__main__":
    mode, n, t = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    if mode == "count":
        sigma = float(sys.argv[4])
        a, b, _ = triples[t]
        reps, glue, oK = quotient_tiles(a, b, *QUOTIENTS["Q1"])
        print({"which": "Q1 count", "triple": t, **count_below(reps, glue, n, sigma)}, flush=True)
    elif mode == "s5":
        print({"which": "C/S5 upper", **s5_upper(n, t)}, flush=True)
    elif mode == "a6":
        print({"which": "C/A6 lower", **a6_lower(n, t, float(sys.argv[4]))}, flush=True)
