"""FRAMEWORK_TOPOLOGICAL_HERSCH.md section 6.  Coarse lambda_1 (P1 FEM, Richardson n=4,6) for every faithful A7-curve with triangle signature (p,q,r):
C tiled by 2*2520 hyperbolic triangles with angles pi/p (A), pi/q (B), pi/r (C), glued as in orbifold.py:
U_g~L_g (AB), U_g~L_{gb} (BC), U_g~L_{ga} (CA); vertex cycles have length 2 ord(a), 2 ord(b), 2 ord(b a^-1).
Curves are enumerated up to simultaneous S7-conjugation and (a,b) -> (a^-1,b^-1): the map U_g <-> L_g is an
orientation-reversing isometry X(a,b) -> X(a^-1,b^-1), so the spectra agree.
Usage: python3 signatures_spectrum.py p q r [fast]
  default: elementwise Richardson on n = 4, 6 (checked on 2 4 7: 0.3464 / 0.3598 against 0.34627 / 0.3597);
  fast:    n = 4 only, divided by 1.0105 (the P1 overshoot measured on (2,4,7) and (3,3,5)).
Numerical only.  The fast-mode calibration does NOT transfer to other triangles (it overestimated (4,4,4) by
about 15%); certified values are produced by certify_signatures.py (VERIFICATION.md section 6)."""
import sys, time
from math import pi, cos, sin, cosh, sinh, acosh, sqrt
from itertools import permutations
import numpy as np
import scipy.sparse.linalg as sla
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a7 import A7, mul, inv, generated, order, cycle_type, conj
import orbifold as ob

S7 = list(permutations(range(7)))


def ref_triangle(p, q, r):
    al, be, ga = pi / p, pi / q, pi / r
    c = acosh((cos(ga) + cos(al) * cos(be)) / (sin(al) * sin(be)))     # AB
    b = acosh((cos(be) + cos(al) * cos(ga)) / (sin(al) * sin(ga)))     # AC
    A = np.array([1.0, 0, 0]); B = np.array([cosh(c), sinh(c), 0]); C = np.array([cosh(b), sinh(b) * cos(al), sinh(b) * sin(al)])
    P = B + C; P = P / sqrt(ob.mink(P, P)); g, pp = P[0], P[1:]
    L = np.eye(3); L[0, 0] = g; L[0, 1:] = -pp; L[1:, 0] = -pp; L[1:, 1:] = np.eye(2) + np.outer(pp, pp) / (1 + g)
    return [(L @ V)[1:] / (L @ V)[0] for V in (A, B, C)]


def curves(p, q, r):
    """S7-orbit representatives of generating pairs (a,b), ord a=p, ord b=q, ord(b a^-1)=r, modulo inversion."""
    seen, reps = set(), []
    cand_a = {}
    for x in A7:
        if order(x) == p:
            cand_a.setdefault(cycle_type(x), x)
    for a in cand_a.values():
        for b in A7:
            if order(b) != q or order(mul(b, inv(a))) != r or (a, b) in seen:
                continue
            if len(generated([a, b])) != 2520:
                continue
            orb = set()
            for s in S7:
                si = inv(s)
                for (x, y) in ((a, b), (inv(a), inv(b))):
                    orb.add((mul(mul(s, x), si), mul(mul(s, y), si)))
            seen |= orb
            reps.append((a, b, len(orb)))
    return reps


def lowest(a, b, tri_pts, n, k=24):
    XA, XB, XC = tri_pts
    mats = ob.numeric_elem_mats(n, XA, XB, XC, "P1")
    reps, glue, _ = ob.quotient_tiles(a, b, [], None)
    K, M = ob.assemble(n, reps, glue, "P1", mats)
    vals = sla.eigsh(K.tocsc(), k=k, M=M.tocsc(), sigma=-1e-3, which="LM", return_eigenvectors=False)
    return np.sort(vals)


if __name__ == "__main__":
    p, q, r = map(int, sys.argv[1:4])
    FAST = len(sys.argv) > 4 and sys.argv[4] == "fast"
    g = 1 + 1260 * (1 - 1 / p - 1 / q - 1 / r)
    thr = 48 / (g - 1)
    t0 = time.time()
    cs = curves(p, q, r)
    print(f"({p},{q},{r}) genus {g:.0f}: {len(cs)} curve(s) up to S7 and mirror; threshold lambda_1 > {thr:.4f} "
          f"({time.time()-t0:.0f}s)", flush=True)
    pts = ref_triangle(p, q, r)
    for (a, b, osz) in cs:
        if FAST:                        # n = 4 only; P1 overshoot at n = 4 is about 1.0-1.1% on (2,4,7), (3,3,5)
            l = lowest(a, b, pts, 4) / 1.0105
        else:                           # elementwise Richardson; only meaningful if no crossing between n = 4, 6
            l4 = lowest(a, b, pts, 4); l6 = lowest(a, b, pts, 6)
            l = np.sort((36 * l6 - 16 * l4) / 20)
        cl = ob.clusters(list(l[1:]), tol=3e-3)
        lam1 = cl[0][0]
        print(f"  a={cycle_type(a)} b={cycle_type(b)} orbit {osz}: lambda ({'n=4 / 1.0105' if FAST else 'Richardson 4,6'}) = "
              f"{[(round(v, 4), m) for v, m in cl[:5]]}  lowest {np.round(l[1:4], 4)}  LY gon >= {l[1]*(g-1)/2:.2f}  "
              f"[{'>24 OK' if l[1] > thr else 'LY FAILS'}]  ({time.time()-t0:.0f}s)", flush=True)
