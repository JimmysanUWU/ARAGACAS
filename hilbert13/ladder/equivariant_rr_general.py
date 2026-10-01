"""Equivariant Riemann-Roch test for small linearised moving degrees on faithful A7-curves of any rigid signature
(p,q,r), generalising equivariant_rr.py.  For each curve (signatures_spectrum.curves), each linearised class of
degree n (Pic^G_lin = <D_1,D_2,D_3 | e_1 D_1 = e_2 D_2 = e_3 D_3>, classes indexed by (n_1 mod e_1, n_2 mod e_2)),
compute chi_G(O(D)) by holomorphic Lefschetz with the product-one monodromy triple (a, b^-1, b a^-1).
A positive multiplicity => h^0(D) > 0 => (n < min orbit size) a base-point-free linearised moving bundle of degree n.

Usage: python3 equivariant_rr_general.py n p q r [p q r ...]
"""
import os
import sys
import cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a7 import mul, inv, order, cycle_type, power
from equivariant_rr import cls, cent, chars, SIZES, TABLE
from signatures_spectrum import curves


def lefschetz(gens, nD, sign):
    L = [None] * 9
    for c in range(1, 9):
        tot = 0
        for gi, k in zip(gens, nD):
            e = order(gi)
            for j in range(1, e):
                if cls(power(gi, j)) != c:
                    continue
                a = cmath.exp(sign * 2j * cmath.pi * j / e)
                tot += cent[c] / e * a ** k / (1 - 1 / a)
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
    n = int(sys.argv[1])
    sigs = [tuple(map(int, sys.argv[i:i + 3])) for i in range(2, len(sys.argv), 3)]
    for (p, q, r) in sigs:
        g = 1 + 1260 * (1 - 1 / p - 1 / q - 1 / r)
        g = round(g)
        cs = curves(p, q, r)
        print(f"({p},{q},{r}) genus {g}: {len(cs)} curve(s); degree {n}, min orbit {2520 // max(p, q, r)}", flush=True)
        for (a, b, osz) in cs:
            gens = (a, inv(b), mul(b, inv(a)))
            assert mul(mul(gens[0], gens[1]), gens[2]) == tuple(range(7))
            es = [order(x) for x in gens]
            assert sorted(es) == sorted((p, q, r))
            degs = [2520 // e for e in es]
            # validation on K = sum (e_i - 1) D_i - 2F
            cK = decompose(lefschetz(gens, [e - 1 for e in es], +1), g - 1)
            assert cK.get("1", 0) == -1 and all(v > 0 for k, v in cK.items() if k != "1")
            assert sum(v * int(TABLE[k][0]) for k, v in cK.items() if k != "1") == g
            hits = []
            for n1 in range(es[0]):
                for n2 in range(es[1]):
                    rest = n - n1 * degs[0] - n2 * degs[1]
                    if rest % degs[2]:
                        continue
                    nD = (n1, n2, rest // degs[2])
                    for sign in (+1, -1):
                        ch = decompose(lefschetz(gens, nD, sign), n - g + 1)
                        pos = {k: v for k, v in ch.items() if v > 0}
                        hits.append((nD, sign, ch, pos))
            print(f"  curve a={cycle_type(a)} b={cycle_type(b)} classes {cycle_type(gens[2])}: "
                  f"{len(hits) // 2} linearised class(es) of degree {n}; H^0(K) validated", flush=True)
            for nD, sign, ch, pos in hits:
                print(f"    D = {nD[0]}D{es[0]} + {nD[1]}D{es[1]} + {nD[2]}D{es[2]}  sign {sign:+d}: chi = {ch}"
                      + (f"   ==> h^0 > 0 ({pos})" if pos else ""), flush=True)
