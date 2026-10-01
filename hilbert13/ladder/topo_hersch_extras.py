"""Topological Hersch (FRAMEWORK_TOPOLOGICAL_HERSCH.md), the finite parts.
  1. [exact, character table] dim (wedge^3 V)^{A7} for the real irreducibles V, and wedge^2 V.
  2. [numerical] the sharpened bracket constant: sup over W = (w1,w2,w3) in E1^3, |W| = 1, of
     F(W) = Q*(w2^w3) + Q*(w3^w1) + Q*(w1^w2), and the decomposable sup q_dec of Q*,
     from the npz written by `topo_hersch.py n t --save`.
  3. [numerical] certification targets for m = 24 on classes 0/1 (Richardson-extrapolated constants).
Usage: python3 topo_hersch_extras.py [th_n8_t0.npz]"""
import os
import sys
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a7 import power
from chartab import TABLE, CL, idx, sizes

# ---------------------------------------------------------------- 1. wedge^3 invariants
names = ["1", "6", "10", "10b", "14a", "14b", "15", "21", "35"]
p2 = [idx[power(c[0], 2)] for c in CL]
p3 = [idx[power(c[0], 3)] for c in CL]


def inner(f):
    return (sizes * f).sum().real / 2520


real = {"6": TABLE[1], "10+10b": TABLE[2] + TABLE[3], "14a (TABLE)": TABLE[4], "14b (TABLE)": TABLE[5],
        "15": TABLE[6], "21": TABLE[7], "35": TABLE[8]}
print("1. real irreducibles of A7 (chartab labels; the lambda_1 isotype 14_(5,2) is TABLE '14b'):")
for nm, x in real.items():
    x2, x3 = x[p2], x[p3]
    w3 = (x ** 3 - 3 * x * x2 + 2 * x3) / 6
    w2 = (x ** 2 - x2) / 2
    dec = " + ".join(f"{int(round(inner(w2 * np.conj(TABLE[i]))))}x{names[i]}"
                     for i in range(9) if round(inner(w2 * np.conj(TABLE[i]))))
    print(f"   {nm:12s} dim {x[0].real:3.0f}: dim (wedge^3)^G = {abs(round(inner(w3)))};  wedge^2 = {dec}")

# ---------------------------------------------------------------- 2. sharpened bracket constant
ratio3, ratio_dec = 0.6143, 0.8191           # values at n = 8, class 0 (recomputed below if an npz is given)
if len(sys.argv) > 1:
    d = np.load(sys.argv[1])
    bw, bV, blocks, res = d["bw"], d["bV"], d["blocks"], d["res"]
    bst = np.zeros(91)
    for (s, e), r in zip(blocks, res):
        bst[s:e] = r[0] * r[2]
    Qs = (bV * bst) @ bV.T
    bmax = bst.max()
    iu = np.triu_indices(14, 1)

    def F_and_grad(x, pairs):
        W = x.reshape(14, 3)
        nrm2 = (W ** 2).sum()
        F, g = 0.0, np.zeros((14, 3))
        for a, b in pairs:
            Om = np.outer(W[:, a], W[:, b]) - np.outer(W[:, b], W[:, a])
            q = Qs @ Om[iu]
            F += Om[iu] @ q
            Gm = np.zeros((14, 14)); Gm[iu] = 2 * q
            Am = Gm - Gm.T
            g[:, a] += Am @ W[:, b]; g[:, b] -= Am @ W[:, a]
        return -F / nrm2 ** 2, -(g / nrm2 ** 2 - 4 * F * W / nrm2 ** 3).reshape(-1)

    rng = np.random.default_rng(7)
    out = {}
    for label, pairs, scale in (("3-plane", [(1, 2), (2, 0), (0, 1)], bmax / 3), ("decomposable", [(0, 1)], bmax / 4)):
        vals = sorted((-minimize(F_and_grad, rng.normal(size=42), args=(pairs,), jac=True, method="BFGS",
                                 options=dict(maxiter=20000, gtol=1e-12)).fun for _ in range(400)), reverse=True)
        out[label] = vals[0] / scale
        print(f"2. {label}: top local maxima / crude = {np.round(np.array(vals[:5]) / scale, 4)}  (400 BFGS starts)")
    ratio3, ratio_dec = out["3-plane"], out["decomposable"]
else:
    print(f"2. (no npz given) using kappa^2 = {ratio3} b*max/3 and q_dec = {ratio_dec} b*max from n = 8, class 0")

# ---------------------------------------------------------------- 3. certification targets
A = 540 * np.pi
bstar, gam = 1.2769e-4, 1.70                 # Richardson from n = 8, 12 (class 0)
Lam = gam ** 2 / A
LAM1_UP = 0.34671                            # P1 value at n = 12 (upper bound up to quadrature)


def margin(lam1, lam_next, c, Lm, m=24):
    eps = 8 * np.pi * m - lam1 * A
    s = eps / (lam_next - lam1)
    lhs = lam1 * (A - s) / 2
    rhs = 3 * np.sqrt(eps) * (A - s) * np.sqrt(c) + np.sqrt(Lm * (A - s)) * np.sqrt(s * (eps + lam1 * s))
    return lhs - rhs, lhs, rhs


def need_lam1(c, L2, t=1.0):
    lo, hi = 0.30, LAM1_UP
    for _ in range(60):
        mid = (lo + hi) / 2
        ok = all(margin(l1, L2, c * t * t, Lam * t * t)[0] > 0 for l1 in np.linspace(mid, LAM1_UP, 9))
        lo, hi = (lo, mid) if ok else (mid, hi)
    return hi


print("3. m = 24, classes 0/1:")
for nm, c in (("crude b*/3", bstar / 3), ("q_dec/3", ratio_dec * bstar / 3), ("sharp kappa^2", ratio3 * bstar / 3)):
    d_, l_, r_ = margin(0.34627, 0.5715, c, Lam)
    print(f"   {nm:14s}: at the true values lhs {l_:.1f}, rhs {r_:.1f};  certified lambda_1 needed: "
          f"{need_lam1(c, 0.56):.5f} (lambda' >= 0.56), {need_lam1(c, 0.55):.5f} (lambda' >= 0.55)")
c = ratio3 * bstar / 3
for L1 in (0.34089, 0.342, 0.344, 0.3455):
    lo, hi = 1.0, 3.0
    for _ in range(50):
        t = (lo + hi) / 2
        ok = all(margin(l1, 0.55, c * t * t, Lam * t * t)[0] > 0 for l1 in np.linspace(L1, LAM1_UP, 9))
        lo, hi = (t, hi) if ok else (lo, t)
    print(f"   certified lambda_1 >= {L1}, lambda' >= 0.55: kappa and sqrt(Lambda) may be inflated by {lo:.3f}x")
