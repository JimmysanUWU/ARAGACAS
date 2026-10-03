"""External validation of the geometry + FEM + certificate pipeline on the Bolza surface.

The Bolza surface (genus 2, y^2 = x^5 - x) is the (2,3,8) triangle curve with group GL(2,3) of order 48,
realised here as permutations of the 8 nonzero vectors of F_3^2.  Its first eigenvalue is known to
high precision:  lambda_1 = 3.8388872588421995  (multiplicity 3; Strohmaier-Uski, Comm. Math. Phys. 317
(2013)).  We check
  (1) conforming P1 eigenvalues on the whole surface converge to it (exact geometry, same code path
      as orbifold.py uses for the A7 curves), and
  (2) the rigorous pipeline (certify.certify_tiles) applied to a sign-twisted quotient that contains the
      lambda_1-isotype produces lower bounds BELOW 3.83888726 that approach it.
Usage: python3 validate_bolza.py
"""
import itertools
import numpy as np
from orbifold import (quotient_tiles, reference_triangle, numeric_elem_mats, assemble, eigs, clusters)
from certify import certify_tiles

LAMBDA1_BOLZA = 3.8388872588421995

VEC = [v for v in itertools.product(range(3), repeat=2) if v != (0, 0)]
IDX = {v: i for i, v in enumerate(VEC)}
N8 = len(VEC)


def mat_perm(M):
    return tuple(IDX[((M[0][0] * x + M[0][1] * y) % 3, (M[1][0] * x + M[1][1] * y) % 3)] for (x, y) in VEC)


E8 = tuple(range(N8))


def mul8(p, q):
    return tuple(p[q[i]] for i in range(N8))


def inv8(p):
    r = [0] * N8
    for i, x in enumerate(p):
        r[x] = i
    return tuple(r)


def order8(p):
    k, q = 1, p
    while q != E8:
        q, k = mul8(q, p), k + 1
    return k


def generated8(gens):
    S, frontier = {E8}, [E8]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = mul8(g, x)
                if y not in S:
                    S.add(y); new.append(y)
        frontier = new
    return S


GL23 = sorted({mat_perm(((a, b), (c, d))) for a, b, c, d in itertools.product(range(3), repeat=4)
               if (a * d - b * c) % 3 != 0})
assert len(GL23) == 48
GROUP = (GL23, mul8, inv8, generated8, E8)


def find_triple():
    for a in GL23:
        if order8(a) != 2:
            continue
        for b in GL23:
            if order8(b) == 3 and order8(mul8(b, a)) == 8 and len(generated8([a, b])) == 48:
                return a, b
    raise ValueError


def main():
    a, b = find_triple()
    print(f"GL(2,3) as permutations of F_3^2 - 0; triple a (order 2), b (order 3), ba (order 8)")
    XA, XB, XC = reference_triangle("mid_bc", pqr=(2, 3, 8))
    reps, glue, _ = quotient_tiles(a, b, None, None, group=GROUP)
    print(f"whole surface: {2 * len(reps)} triangles (pi/2, pi/3, pi/8)")
    for n in (4, 8, 16):
        K, M = assemble(n, reps, glue, "P1", numeric_elem_mats(n, XA, XB, XC, "P1"))
        area = M.sum()
        cl = clusters(eigs(K, M, k=6), 1e-6)
        print(f"  P1 n={n:2d}: area {area:.8f} (4 pi = {4 * np.pi:.8f})   eigenvalues {cl}")
    # a sign-twisted quotient containing the lambda_1 isotype: search cyclic subgroups of even order
    best = None
    for k in GL23:
        o = order8(k)
        if o % 2:
            continue
        r2, g2, _ = quotient_tiles(a, b, [k], [-1], group=GROUP)
        K, M = assemble(8, r2, g2, "P1", numeric_elem_mats(8, XA, XB, XC, "P1"))
        lam = eigs(K, M, k=2)[0]
        if abs(lam - LAMBDA1_BOLZA) < 0.02 and (best is None or o > best[1]):
            best = (k, o, lam, len(r2))
    k, o, lam, nc = best
    print(f"twisted quotient: K = <k>, |K| = {o}, eps(k) = -1, {2 * nc} triangles; lowest P1 eigenvalue (n=8) {lam:.6f}")
    r2, g2, _ = quotient_tiles(a, b, [k], [-1], group=GROUP)
    for n in (8, 16, 32, 64):
        res = certify_tiles(r2, g2, n, pqr=(2, 3, 8), sigma_factor=0.999)
        assert res["ok"]
        ok = res["bound"] < LAMBDA1_BOLZA
        print(f"  certified n={n:2d}: lambda_h {res['lam_h'][0]:.6f}  bound {res['bound']:.6f}  "
              f"< lambda_1(Bolza) = {LAMBDA1_BOLZA:.6f}: {ok}")
        assert ok


if __name__ == "__main__":
    main()
