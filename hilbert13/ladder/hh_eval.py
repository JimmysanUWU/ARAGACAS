"""Rigorous (ball-arithmetic) evaluation of the Hejhal trial eigenfunction of a (2,4,7) A7-curve.

The trial function.  Let F0(z) = sum_{m=0}^{M} R_m(|z|) (A_m cos m th + B_m sin m th), z = |z| e^{i th} in the Poincare disk
(order-7 point at 0), A_m, B_m in R^21 (exact binary floats from hejhal_solve.py), with
    R_m(u) = Re[ u^m (1-u^2)^a 2F1(a, m+1/2-it; m+1; u^2) ],  a = 1/2 - it,  lam* = 1/4 + t^2,
the regular solution of (Delta - lam*) (R_m e^{i m th}) = 0 (positive hyperbolic Laplacian, curvature -1).
F0 is an EXACT lam*-eigenfunction on the whole disk.  The certified trial function uses
    Ft(z) = P_V2 (1/7) sum_k P_c^{-k} F0(zeta^k z),
which is exactly C-equivariant (Ft(zeta z) = P_c Ft(z)), V2-valued (V2 = 14_(5,2) inside the 2-subset module R^21),
and still an exact lam*-eigenfunction; its coefficients are enclosed in balls below.
"""
import os, sys, itertools
import numpy as np
from flint import arb, acb, arb_mat, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from triples_data import triples
from a7 import inv, mul

ctx.prec = 106
PAIRS = list(itertools.combinations(range(7), 2)); PIDX = {p: i for i, p in enumerate(PAIRS)}


def perm21_index(p):
    """P_p e_{ij} = e_{p(i)p(j)}: returns sigma with (P_p v)[sigma[r]] = v[r]."""
    return [PIDX[tuple(sorted((p[i], p[j])))] for (i, j) in PAIRS]


def apply_perm(sig, v):
    out = [None] * 21
    for r, s in enumerate(sig):
        out[s] = v[r]
    return out


def PV2_rational():
    """Exact orthogonal projector onto V2 = (span of x_i + x_j)^perp, as a 21x21 list of Fractions."""
    from fractions import Fraction as Fr
    inc = [[1 if i in pr else 0 for i in range(7)] for pr in PAIRS]
    # (Inc^T Inc)^-1 = (1/5)(I - J/12)
    G = [[Fr(1, 5) * ((1 if i == j else 0) - Fr(1, 12)) for j in range(7)] for i in range(7)]
    P = [[Fr(0)] * 21 for _ in range(21)]
    for r in range(21):
        for s in range(21):
            v = sum(inc[r][i] * G[i][j] * inc[s][j] for i in range(7) for j in range(7))
            P[r][s] = (1 if r == s else 0) - v
    return P


def f21(x, a, b, c):
    """2F1(a,b;c;x) (acb); arb's evaluation can return NaN at working precision for large c: retry at higher
    precision (the result is a rigorous enclosure either way)."""
    v = acb(x).hypgeom_2f1(a, b, c)
    p0 = ctx.prec
    while not v.is_finite() and ctx.prec < 1024:
        ctx.prec = ctx.prec * 2
        v = acb(x).hypgeom_2f1(a, b, c)
    ctx.prec = p0
    if not v.is_finite():
        raise ArithmeticError("2F1 evaluation failed")
    return v


class TrialFunction:
    def __init__(self, cls, coef_file):
        d = np.load(coef_file)
        self.lam = arb(float(d["lam"]))                      # lam* : an exact binary float
        A0, B0 = d["A"], d["B"]
        self.M = A0.shape[0] - 1
        a, b, c = triples[cls]
        self.perm = {"a": a, "b": b, "c": c}
        self.t = (self.lam - arb(1) / 4).sqrt()
        self.aa = acb(arb(1) / 2, -self.t)
        # C-symmetrisation and V2 projection of the coefficients (balls)
        PV2 = PV2_rational()
        PV2a = arb_mat(21, 21, [arb(PV2[r][s].numerator) / PV2[r][s].denominator for r in range(21) for s in range(21)])
        cinv = inv(c)
        sig_cinv = perm21_index(cinv)
        two_pi7 = 2 * arb.pi() / 7
        self.A, self.B = [], []
        for m in range(self.M + 1):
            Am = [arb(float(x)) for x in A0[m]]; Bm = [arb(float(x)) for x in B0[m]]
            accA = [arb(0)] * 21; accB = [arb(0)] * 21
            for k in range(7):
                ck, sk = (m * k * two_pi7).cos(), (m * k * two_pi7).sin()
                vA = [Am[r] * ck + Bm[r] * sk for r in range(21)]
                vB = [-Am[r] * sk + Bm[r] * ck for r in range(21)]
                for _ in range(k):                           # apply P_c^{-k}
                    vA = apply_perm(sig_cinv, vA); vB = apply_perm(sig_cinv, vB)
                accA = [x + y for x, y in zip(accA, vA)]; accB = [x + y for x, y in zip(accB, vB)]
            colA = PV2a * arb_mat(21, 1, [x / 7 for x in accA]); colB = PV2a * arb_mat(21, 1, [x / 7 for x in accB])
            self.A.append([colA[r, 0] for r in range(21)]); self.B.append([colB[r, 0] for r in range(21)])
        # coefficient matrix [A_0..A_M | B_0..B_M] (21 x 2(M+1))
        self.Cmat = arb_mat(21, 2 * (self.M + 1), [self.A[m][r] if j < self.M + 1 else self.B[j - self.M - 1][r]
                                                for r in range(21) for j, m in
                                                ((jj, jj if jj < self.M + 1 else jj - self.M - 1) for jj in range(2 * (self.M + 1)))])
        self.absA = [sum(abs(float(x.mid())) + float(x.rad()) for x in self.A[m]) for m in range(self.M + 1)]

    # ---------------------------------------------------------------- radial functions
    def radial(self, u, deriv=False, second=False):
        """For an arb u (point or interval, 0 <= u < 1), with R_m = u^m S_m(u^2):
        R_m(u); if deriv: R_m', Rdiv_m = R_m/u (m >= 1, else 0); if second also R_m'' (from the radial ODE
        R'' = -R'/u + m^2 R/u^2 - 4 lam R/(1-u^2)^2, written without singular terms) and Rdiv_m' = d/du (R_m/u)."""
        x = u * u; one_x = 1 - x
        pre = (acb(one_x).log() * self.aa).exp()           # (1-x)^a
        upow = [arb(1)]
        for m in range(self.M + 1):
            upow.append(upow[-1] * u)
        Rs, dRs, Rdiv, d2Rs, dRdiv = [], [], [], [], []
        four_lam_q = 4 * self.lam / (one_x * one_x)
        for m in range(self.M + 1):
            bm = acb(arb(m) + arb(1) / 2, -self.t); cm = acb(m + 1)
            F1 = f21(x, self.aa, bm, cm)
            S = pre * F1
            Rs.append((S * upow[m]).real)
            if deriv or second:
                F2 = f21(x, self.aa + 1, bm + 1, cm + 1)
                dS = pre * (self.aa * bm / cm * F2 - self.aa * F1 / acb(one_x))      # dS/dx
                val = dS * upow[m + 1] * 2
                if m >= 1:
                    val = val + S * upow[m - 1] * m
                    Rdiv.append((S * upow[m - 1]).real)
                else:
                    Rdiv.append(arb(0))
                dRs.append(val.real)
                if second:
                    # R'' = m(m-1) u^{m-2} S - 2 u^m S' - 4 lam u^m S/(1-x)^2   (S' = dS/dx)
                    v2 = -2 * dS * upow[m] - four_lam_q * S * upow[m]
                    if m >= 2:
                        v2 = v2 + m * (m - 1) * S * upow[m - 2]
                    d2Rs.append(v2.real)
                    # d/du (u^{m-1} S) = (m-1) u^{m-2} S + 2 u^m S'   (m >= 1)
                    if m >= 1:
                        v3 = 2 * dS * upow[m]
                        if m >= 2:
                            v3 = v3 + (m - 1) * S * upow[m - 2]
                        dRdiv.append(v3.real)
                    else:
                        dRdiv.append(arb(0))
        if second:
            return Rs, dRs, Rdiv, d2Rs, dRdiv
        return (Rs, dRs, Rdiv) if deriv else Rs

    def trig(self, th):
        c1, s1 = th.cos(), th.sin()
        cs, ss = [arb(1)], [arb(0)]
        for m in range(1, self.M + 1):
            cs.append((m * th).cos()); ss.append((m * th).sin())
        return cs, ss, c1, s1

    def combine(self, vec):
        col = self.Cmat * arb_mat(2 * (self.M + 1), 1, vec)
        return [col[r, 0] for r in range(21)]

    def value(self, rad, trg):
        Rs = rad if not isinstance(rad, tuple) else rad[0]
        cs, ss, _, _ = trg
        return self.combine([Rs[m] * cs[m] for m in range(self.M + 1)] + [Rs[m] * ss[m] for m in range(self.M + 1)])

    def value_grad(self, rad, trg):
        """F, F_x, F_y (Euclidean disk-coordinate derivatives)."""
        Rs, dRs, Rdiv = rad[:3]
        cs, ss, c1, s1 = trg
        M1 = self.M + 1
        F = self.combine([Rs[m] * cs[m] for m in range(M1)] + [Rs[m] * ss[m] for m in range(M1)])
        Fu = self.combine([dRs[m] * cs[m] for m in range(M1)] + [dRs[m] * ss[m] for m in range(M1)])
        Fthu = self.combine([-m * Rdiv[m] * ss[m] for m in range(M1)] + [m * Rdiv[m] * cs[m] for m in range(M1)])
        Fx = [fu * c1 - ft * s1 for fu, ft in zip(Fu, Fthu)]
        Fy = [fu * s1 + ft * c1 for fu, ft in zip(Fu, Fthu)]
        return F, Fx, Fy

    def polar_derivs(self, rad2, trg):
        """On a polar box (rad2 = radial(U, second=True), trg = trig(Theta)): enclosures of
        dF/du, dF/dth, and of d/du, d/dth of the Euclidean gradient G = (F_x, F_y)."""
        Rs, dRs, Rdiv, d2Rs, dRdiv = rad2
        cs, ss, c1, s1 = trg
        M1 = self.M + 1
        cmb = lambda v, w: self.combine([v[m] * cs[m] for m in range(M1)] + [v[m] * ss[m] for m in range(M1)]) \
            if w is None else self.combine([-w[m] * ss[m] for m in range(M1)] + [w[m] * cs[m] for m in range(M1)])
        mR = [m * Rdiv[m] for m in range(M1)]
        Fu = cmb(dRs, None)
        Fthu = cmb(None, mR)                                       # F_th / u
        Fuu = cmb(d2Rs, None)
        Futh = cmb(None, [m * dRs[m] for m in range(M1)])            # d/dth F_u
        Fthu_th = cmb([-m * m * Rdiv[m] for m in range(M1)], None)   # d/dth (F_th/u)
        Fthu_u = cmb(None, [m * dRdiv[m] for m in range(M1)])         # d/du (F_th/u)
        U = None
        Gx_u = [a * c1 - b * s1 for a, b in zip(Fuu, Fthu_u)]
        Gy_u = [a * s1 + b * c1 for a, b in zip(Fuu, Fthu_u)]
        Gx_t = [a * c1 - b * s1 - c * s1 - d * c1 for a, b, c, d in zip(Futh, Fu, Fthu_th, Fthu)]
        Gy_t = [a * s1 + b * c1 + c * c1 - d * s1 for a, b, c, d in zip(Futh, Fu, Fthu_th, Fthu)]
        return Fu, Fthu, Gx_u, Gy_u, Gx_t, Gy_t

    def at_point(self, z, deriv=False):
        """z: acb point (ball). Returns F (and grads)."""
        u = abs(z); th = z.arg()
        rad = self.radial(u, deriv)
        trg = self.trig(th)
        return self.value_grad(rad, trg) if deriv else self.value(rad, trg)


if __name__ == "__main__":
    import time
    here = os.path.dirname(os.path.abspath(__file__))
    tf = TrialFunction(0, os.path.join(here, "coef_cls0_M90.npz"))
    sys.path.insert(0, "/tmp/claude-0/-home-user-ARAGACAS/32e8ed8f-a61f-5f3c-8c96-fdba0d6af448/scratchpad")
    t0 = time.time()
    z = acb(arb("0.31"), arb("0.12"))
    F, Fx, Fy = tf.at_point(z, deriv=True)
    print("time", time.time() - t0)
    print("F[:3]", F[:3]); print("Fx[:3]", Fx[:3])
    from hejhal2 import radial_and_deriv
    d = np.load(os.path.join(here, "coef_cls0_M90.npz"))
    zz = np.array([0.31 + 0.12j]); u = abs(zz); th = np.angle(zz)
    out = np.zeros((21, 1)); h = 1e-6
    def F0(zz):
        u = np.abs(zz); th = np.angle(zz); out = np.zeros((21, len(zz)))
        for m in range(d["A"].shape[0]):
            R, _ = radial_and_deriv(m, u, float(d["lam"]))
            out += np.outer(d["A"][m], R*np.cos(m*th)) + np.outer(d["B"][m], R*np.sin(m*th))
        return out
    f = F0(zz)[:, 0]; fx = (F0(zz + h)[:, 0] - F0(zz - h)[:, 0]) / (2 * h)
    print("float F0[:3]", f[:3], "fx", fx[:3])
    print("max |Ft - F0|", max(abs(float(F[r].mid()) - f[r]) for r in range(21)), "max rad", max(float(x.rad()) for x in F))
