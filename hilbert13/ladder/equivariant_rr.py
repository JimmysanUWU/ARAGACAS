"""Equivariant Euler characteristics chi_G(O(D)) = [H^0] - [H^1] on a (2,4,7) A7-curve, for invariant divisors
D = n2 D2 + n4 D4 + n7 D7 (D_e = reduced fibre over the branch point of order e) with the canonical linearisation.

Holomorphic Lefschetz (Atiyah-Bott) on a curve, for g != 1:
    sum_q (-1)^q tr(g | H^q(O(D))) = sum_{p in Fix(g)} a_p^{k_p} / (1 - a_p^{-1}),
a_p = eigenvalue of dg_p on T_pC, k_p = mult_p(D).  (Fibre of O(D) at p: frame w^{-k}, g.w = a^{-1} w, so g acts by a^k.
Checked on P^1 with O and K.)  Points over branch i <-> cosets x<g_i>; g fixes x<g_i> with rotation zeta_{e_i}^j iff
x^{-1} g x = g_i^j; there are |C(g)| [g ~ g_i^j] / e_i such points.  Only the classes of the g_i matter.

Validation: chi(O) = 1 - H^0(K)^*, chi(K) = H^0(K) - 1 must be consistent (integral, H^0(K) >= 0, dim 136).
Application: the two degree-90 linearised classes B = 2D7 - D4 and B + T = D2 - 3D4 + 2D7 (Round 5, section 2).
A positive multiplicity in chi(B) forces h^0(B) > 0, i.e. mu_C(A7) = 90.
"""
import os
import sys
import cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a7 import A7, mul, inv, order, cycle_type, a7_class, cyc, power
from cover import TABLE, CLASSES, SIZES
from triples_data import triples

r7 = set(a7_class(cyc(tuple(range(7)))))


def cls(g):
    ct = cycle_type(g)
    return (7 if g in r7 else 8) if ct == (7,) else CLASSES.index(ct)


reps = {}
for g in A7:
    reps.setdefault(cls(g), g)
cent = [2520 // s for s in SIZES]
chars = {k: [complex(v) for v in row] for k, row in TABLE.items()}


def lefschetz(gens, nD, sign=+1):
    """L[c] for each class c (c = 0 is the identity: returns None there)."""
    L = [None] * 9
    for c in range(1, 9):
        g = reps[c]
        tot = 0
        for gi, k in zip(gens, nD):
            e = order(gi)
            for j in range(1, e):
                if cls(power(gi, j)) != c:
                    continue
                npts = cent[c] / e
                a = cmath.exp(sign * 2j * cmath.pi * j / e)
                tot += npts * a ** k / (1 - 1 / a)
        L[c] = tot
    return L


def decompose(L, dim):
    vals = [dim] + L[1:]
    out = {}
    for name, ch in chars.items():
        m = sum(SIZES[c] * vals[c] * ch[c].conjugate() for c in range(9)) / 2520
        assert abs(m.imag) < 1e-6 and abs(m.real - round(m.real)) < 1e-6, (name, m)
        if round(m.real):
            out[name] = round(m.real)
    return out


if __name__ == "__main__":
    for t in (0, 1, 12, 14):
        a, b, c = triples[t]
        gens = (a, inv(b), mul(b, a))                    # monodromy triple (NOTES 7.3)
        assert [order(x) for x in gens] == [2, 4, 7]
        print(f"triple {t}: 7-class of g_3 = {'7a' if cls(gens[2]) == 7 else '7b'}")
        for sign in (+1, -1):
            cO = decompose(lefschetz(gens, (0, 0, 0), sign), 1 - 136)
            cK = decompose(lefschetz(gens, (1, 3, 6), sign), 270 - 135)   # K = -2F + D2 + 3D4 + 6D7 ~ D2+3D4+6D7-2F
            H0K = {k: v for k, v in cK.items() if k != "1"}
            assert cK.get("1", 0) == -1 and all(v > 0 for v in H0K.values())   # H^1(K) = trivial, H^0(K)^G = 0
            assert sum(v * int(TABLE[k][0]) for k, v in H0K.items()) == 136
            res = {"chi(O)": cO, "H^0(K)": H0K}
            for nm, nD in (("B = 2D7 - D4", (0, -1, 2)), ("B+T = D2 - 3D4 + 2D7", (1, -3, 2)),
                           ("2B = 4D7 - 2D4", (0, -2, 4))):
                deg = 1260 * nD[0] + 630 * nD[1] + 360 * nD[2]
                res[nm] = decompose(lefschetz(gens, nD, sign), deg - 135)
            print(f"  rotation sign {sign:+d}:")
            for k, v in res.items():
                pos = {x: m for x, m in v.items() if m > 0}
                print(f"    {k:22s} {v}" + (f"   POSITIVE PART {pos}" if k.startswith(("B", "2B")) and pos else ""))
