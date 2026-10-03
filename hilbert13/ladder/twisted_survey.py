"""Twisted models and fixed-point pencils on every rigid faithful A7-curve of genus <= 335 (7_TWISTED.md, section 7.8).

For each curve (S7-class of generating triples) and each Schur twist eps (incl. eps = 1), every invariant class of degree
<= 2g-2 whose holomorphic-Lefschetz character forces sections W <= H^0 is tested against every element x with fixed points:
the sections in an eigenspace E_mu of a lift of x vanish at the fixed points whose fibre eigenvalue is not mu, so a pencil
inside E_mu has degree <= D - (those base points) - (dim E_mu - 2).  Prints the best such bound per curve.
Run:  python3 twisted_survey.py > twisted_survey_output.txt   (about 90 s; imports twisted_rr.py)
"""
import sys, io, contextlib, cmath, itertools, time
from fractions import Fraction as Fr
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):
    import twisted_rr as T
from a7 import A7, mul, inv, order, cycle_type
from signatures_spectrum import curves
G, IDX, n, gm, gi_, split, chars, cls, deg, chiZ, MUL, INV, CC = (T.G, T.IDX, T.n, T.gm, T.gi_, T.split, T.chars, T.cls,
                                                                T.deg, T.chiZ, T.MUL, T.INV, T.CC)
kc = len(chars); N6 = 6 * n
def mults(TH, eps, D, gm1):
    TH = TH.copy()
    for s in range(2):
        for t in range(3): TH[n * (s + 2 * t)] = (D - gm1) * chiZ(eps, s, t)
    return {i: np.sum(TH * np.conj(chars[i][cls])) / N6 for i in range(kc) if CC[i] == eps}
def classes_g(gens, eps, sgn, gm1):
    es = [order(x) for x in gens]; cp = []
    for gi in gens:
        yv = IDX[gi]; zz = yv
        for _ in range(order(gi) - 1): zz = gm(zz, yv)
        g_, s_, t_ = split(zz); assert g_ == 0; cp.append((s_, t_))
    _, s0, t0 = split(gm(gm(IDX[gens[0]], IDX[gens[1]]), IDX[gens[2]]))
    rho0 = Fr(eps[0] * s0, 2) + Fr(eps[1] * t0, 3)
    opts = [[m for m in range(6 * e) if abs(cmath.exp(2j * cmath.pi * m / 6) - chiZ(eps, *cp[i])) < 1e-9] for i, e in enumerate(es)]
    out = []
    for ld in itertools.product(*opts):
        lam = [cmath.exp(2j * cmath.pi * ld[i] / (6 * es[i])) for i in range(3)]
        sg = 1 if sgn > 0 else -1
        D = int((2520 * (sum(Fr(sg * ld[i], 6 * es[i]) for i in range(3)) - sg * rho0)) % 2520)
        m = mults(T.theta(gens, lam, eps, sgn), eps, D, gm1)
        assert all(abs(v.imag) < 1e-6 and abs(v.real - round(v.real)) < 1e-6 for v in m.values())
        out.append((D, ld, {i: round(v.real) for i, v in m.items()}, lam))
    return out
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
                key = complex(round(val.real, 6), round(val.imag, 6)); out[key] = out.get(key, 0) + 1 / e
    return {k: round(v) for k, v in out.items()}
def eig_dims(L, Wd, hW):
    """eigenvalue -> multiplicity of the lift L on W = sum of irreps Wd."""
    o = 1; z = L
    while split(z)[0] != 0 or split(z)[1:] != (0, 0):
        z = gm(z, L); o += 1
    pw = [0, L]
    for j in range(2, o): pw.append(gm(pw[-1], L))
    res = {}
    for r in range(o):
        mu = cmath.exp(2j * cmath.pi * r / o)
        d = sum((sum(mult * chars[i][cls[pw[j]]] for i, mult in Wd) if j else hW) * mu ** (-j) for j in range(o)) / o
        if round(d.real): res[complex(round(mu.real, 6), round(mu.imag, 6))] = round(d.real)
    return res
SIGS = [(2, 4, 7), (3, 3, 5), (2, 5, 7), (3, 3, 6), (3, 4, 4), (2, 6, 7), (3, 3, 7), (2, 7, 7), (3, 4, 5), (3, 4, 6), (4, 4, 4)]
t0 = time.time()
for sig in SIGS:
    gg = 1 + Fr(2520, 2) * (1 - sum(Fr(1, e) for e in sig)); gm1 = int(gg) - 1
    for perm in sorted(set(itertools.permutations(sig))):
        cs = curves(*perm)
        if cs: break
    for ci, (a, b, osz) in enumerate(cs):
        gens = (a, inv(b), mul(b, inv(a))); es = [order(x) for x in gens]
        # representatives of elements with fixed points: powers of the generators, one per A7-class
        fixreps = {}
        for gi in gens:
            for j in range(1, order(gi)):
                x = gi
                for _ in range(j - 1): x = mul(x, gi)
                fixreps.setdefault(cycle_type(x), x)
        best = (10 ** 9, None)
        for eps in [(0, 0), (1, 0), (0, 1), (0, 2), (1, 1), (1, 2)]:
            for sgn in (+1,):
                for D, ld, m, lam in classes_g(gens, eps, sgn, gm1):
                    if D > 2 * gm1: continue
                    Wd = [(i, v) for i, v in m.items() if v > 0]
                    if not Wd: continue
                    hW = sum(deg[i] * v for i, v in Wd)
                    if D - 0 >= best[0] + 40: continue
                    for ct, x in fixreps.items():
                        xi = IDX[x]
                        for s in range(2):
                            for u in range(3):
                                Lx = xi + n * (s + 2 * u)
                                fe = fixed_eigs(gens, lam, eps, xi, Lx)
                                for mu, dim in eig_dims(Lx, Wd, hW).items():
                                    base = sum(c for lv, c in fe.items() if abs(lv - mu) > 1e-6)
                                    bd = D - base - max(dim - 2, 0)
                                    if dim >= 2 and bd < best[0]:
                                        best = (bd, f"D={D} eps={eps} W={hW} elt {ct} mu={mu} dim {dim} base {base}")
                                break
                            break
        print(f"{sig} g={gm1 + 1} curve {ci}: best fixed-point pencil degree <= {best[0]}  [{best[1]}]  ({time.time() - t0:.0f}s)", flush=True)
