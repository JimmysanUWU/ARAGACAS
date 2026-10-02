"""Equivariant Plucker test for G-stable linear series on faithful A7-curves with three branch points.

Setting: C -> C/A7 = P^1 branched at three points with monodromy g_1 g_2 g_3 = 1 (orders e_i), L = O(sum n_i D_i)
canonically linearised (D_i = reduced fibre over branch point i), W subset H^0(L) an A7-subrepresentation of dim r+1.

(1) Plucker: sum_p w_p = (r+1)(d + r(g-1)), w_p = sum_i (o_i(p) - i), o_0 < ... < o_r the vanishing orders of W at p.
(2) At p over branch i with stabiliser <c>, c rotating T_pC by a (primitive e-th root), the graded piece of W in
    vanishing order o is L_p (x) (T_p^*)^o, on which c acts by a^{n_i} a^{-o}.  Hence the eigenvalues of c on W
    determine the vanishing orders modulo e: #{i : o_i = m mod e} = mult of eigenvalue a^{n_i - m}.
    The least weight w_min(i) takes the smallest available integers in each residue class; every admissible weight is
    w_min(i) + e t (t >= 0), contributing (2520/e)(w_min + e t) = |orbit| w_min + 2520 t.  Free orbits contribute 2520 w.
(3) Necessary condition:  T := (r+1)(d + r(g-1)) - sum_i (2520/e_i) w_min(i)  is >= 0 and divisible by 2520.

Rotation convention: g_i rotates by exp(2 pi i / e_i) -- the convention validated by equivariant_rr.py; the conjugate
convention is equivalent to swapping 10 <-> 10b, so both W = 10 and W = 10b are always tested.
"""
import os
import sys
import cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a7 import order, power, mul, inv, cycle_type
from equivariant_rr import cls, chars
from signatures_spectrum import curves


def eigen_mults(g, W):
    """multiplicity of exp(2 pi i m / o) as an eigenvalue of g on the representation W = {irrep: mult}."""
    o = order(g)
    vals = [sum(mu * chars[nm][cls(power(g, j))] for nm, mu in W.items()) for j in range(o)]
    out = []
    for m in range(o):
        s = sum(vals[j] * cmath.exp(-2j * cmath.pi * m * j / o) for j in range(o)) / o
        assert abs(s.imag) < 1e-8 and abs(s.real - round(s.real)) < 1e-8
        out.append(round(s.real))
    return out


def w_min(g, n_i, W):
    """least inflection weight at a point with stabiliser <g> (g rotates by exp(2 pi i/e)), fibre exponent n_i."""
    e = order(g)
    mult = eigen_mults(g, W)                  # eigenvalue exp(2 pi i m/e)
    # eigenvalue a^{n_i - o} = exp(2 pi i (n_i - o)/e)  =>  o = n_i - m (mod e)
    orders = []
    for m, mu in enumerate(mult):
        res = (n_i - m) % e
        orders += [res + e * t for t in range(mu)]
    orders.sort()
    return sum(orders) - sum(range(len(orders)))


def test(gens, nD, W, g):
    es = [order(x) for x in gens]
    d = sum(n * (2520 // e) for n, e in zip(nD, es))
    r1 = sum(mu * int(chars[nm][0].real) for nm, mu in W.items())
    total = r1 * (d + (r1 - 1) * (g - 1))
    ws = [w_min(x, n, W) for x, n in zip(gens, nD)]
    T = total - sum((2520 // e) * w for e, w in zip(es, ws))
    return d, r1, total, ws, T


if __name__ == "__main__":
    from triples_data import triples
    # validation 1: canonical series of the (2,4,7) curves (Weierstrass weight g^3 - g)
    for t in (0, 1, 12, 14):
        a, b, c = triples[t]
        gens = (a, inv(b), mul(b, inv(a)))
        for HK in ({"10": 1, "10b": 2, "15": 1, "21": 1, "35": 2}, {"10": 2, "10b": 1, "15": 1, "21": 1, "35": 2}):
            d, r1, tot, ws, T = test(gens, (1, 3, 6), HK, 136)   # K ~ D2 + 3D4 + 6D7 - 2F (F has no fixed points)
            d -= 2 * 2520
            tot = r1 * (d + (r1 - 1) * 135)
            T = tot - sum((2520 // order(x)) * w for x, w in zip(gens, ws))
            print(f"canonical, triple {t}, H0(K)={HK}: total {tot} (g^3-g={136**3-136}), w_min {ws}, T/2520 = {T/2520}")
