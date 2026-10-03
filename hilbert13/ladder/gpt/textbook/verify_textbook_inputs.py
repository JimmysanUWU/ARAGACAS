"""Exact finite inputs used in the textbook edition.

This checks genus arithmetic, feasible Segre chains, the proper-compositum
maximum, the simultaneous-cluster star, Hilbert deficit polynomials, ordinary
Chevalley--Weil multiplicities, and fixed-divisor projector nonvanishing.
It does not certify FEM assembly, sparse factorization, sampled quotient ranks,
or the classical geometric theorems applying these finite inputs.
Run from any directory with numpy and sympy installed.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def castelnuovo(n, r):
    q, s = divmod(n - 1, r - 1)
    return (r - 1) * comb(q, 2) + q * s


def omega(n):
    out, p = 0, 2
    while n > 1:
        if n % p:
            p += 1
        else:
            n //= p
            out += 1
    return out


def genus_arithmetic():
    assert 2520 * (1 - F(1, 2) - F(1, 4) - F(1, 7)) == 270
    assert 270 == 2 * (2 * 64 - 2) + 18
    assert castelnuovo(60, 5) == 406
    assert castelnuovo(30, 5) == 91
    assert castelnuovo(54, 5) == 325
    assert castelnuovo(36, 5) == 136
    assert castelnuovo(42, 5) == 190
    assert F(135, 2) * F(34089, 100000) > 23
    assert F(355695, 1000000) > F(48, 135)
    print("PASS: Riemann--Hurwitz, Castelnuovo and strict spectral arithmetic")


def finite_pencil_bounds():
    for d in range(2, 25):
        for t in range(3, 2 + omega(d)):
            if t * d >= 2**t - 1:
                assert castelnuovo(t * d, 2**t - 1) <= castelnuovo(3 * d, 7)
    values = [(F(d*d, e) + e*d - 3*d + 1, d, e)
              for d in range(2, 25) for e in range(2, d) if d % e == 0]
    assert max(x[0] for x in values) == 265
    assert [(d, e) for v, d, e in values if v == 265] == [(24, 2), (24, 12)]
    for d in range(18, 25):
        cost = min(comb(r, 2) + comb(d-r, 2) for r in range(d + 1))
        assert 2 * cost > (d - 1)**2 - 266
    stars = [(d, r) for d in range(18, 25) for r in range(d + 1)
             if comb(r, 2) + 33 * comb(d-r, 2) <= (d - 1)**2 - 266]
    assert stars == [(24, 23)]
    for d, expected in [(20, 271), (21, 300), (22, 331), (23, 363), (24, 397)]:
        cost = min(comb(r, 2) + comb(d-r, 2) for r in range(d + 1))
        assert max((d-1)**2 - cost, castelnuovo(3*d, 7)) == expected
    print("PASS: every feasible chain, proper-compositum maximum 265, sole star (24,23)")


def polynomial_bounds():
    # Coefficients in q, with m=6q+r. This certifies all q>=1,
    # rather than extending a finite numerical sample to an infinite claim.
    slopes = [-8, -2, 4, 10, 16, 22]
    constants = [1, 0, 0, 2, 5, 8]
    for r in range(6):
        difference = [18-18, 6*r-6-slopes[r], r*r//2-r+1-constants[r]]
        assert difference == [0, 2, int(r == 2)]
    # Direct deficit sums corroborate the explicit residue polynomials.
    for m in range(8, 80):
        q, r = divmod(m, 6)
        bound = 18*q*q + slopes[r]*q + constants[r]
        deficit = sum(max(0, 3*m - (18*(i//2) + (7 if i % 2 else 1)))
                      for i in range(1, m + 1))
        assert bound == deficit
    import sympy as sp
    s = sp.Symbol("s", integer=True, positive=True)
    for r, expected in [(0, s/2), (1, (s+1)/2), (2, s/2+1)]:
        m = 3*s+r
        lam, eps, mu = (s-1, 8, 2) if r == 0 else (s, 3*r-1, 0)
        pi2 = 9*lam*(lam-1)/2 + lam*(eps+2) + mu
        assert sp.expand(m*m/2-m+1-pi2-expected) == 0
    print("PASS: six Hilbert polynomials and three Petrakiev comparisons for every parameter")


def exact_characters_and_divisors():
    import numpy as np
    import sympy as sp
    from a7 import A7, centralizer, cyc, inv, mul, order, power
    from chartab import EXACT_TABLE, idx
    from triples_data import triples
    from curve_checks import PERMS, fix_divisor, names

    def fixed_dimension(row, g):
        value = sp.simplify(sp.sympify(sum(row[idx[power(g, i)]] for i in range(order(g)))) / sp.Integer(order(g)))
        assert value.is_Integer
        return int(value)

    expected_h1 = {"10": 3, "10b": 3, "15": 2, "21": 2, "35": 4}
    for cls in (0, 1, 12, 14):
        triple = triples[cls]
        result = {}
        for name, row in zip(names, EXACT_TABLE):
            mult = 0 if name == "1" else int(row[idx[tuple(range(7))]]) - sum(fixed_dimension(row, g) for g in triple)
            if mult:
                result[name] = mult
        assert result == expected_h1
    tau = triples[0][0]
    central = centralizer(tau)
    induced = {}
    for name, row in zip(names, EXACT_TABLE):
        value = sp.simplify(sp.sympify(sum(row[idx[g]] for g in central)) / sp.Integer(len(central)))
        assert value.is_Integer
        if value:
            induced[name] = int(value)
    assert induced == {"1": 1, "6": 1, "14a": 2, "14b": 1, "21": 1, "35": 1}
    print("PASS: exact H1 characters in all four classes and the involution permutation module")

    mu = cyc((0, 2), (1, 3))
    small = [cyc((0, 1), (4, 5)), cyc((2, 3), (4, 5))]
    delta = fix_divisor(mu) + fix_divisor(mul(tau, mu)) - sum(fix_divisor(g) for g in small)
    good = [cyc((0, 1), (i, j)) for i, j in [(4, 5), (4, 6), (5, 6)]]
    good += [mul(tau, g) for g in good]
    lifted = sum(fix_divisor(g) for g in good) - 3*(fix_divisor(mu)+fix_divisor(mul(tau, mu)))
    for label, raw in [("Klein difference", delta), ("lifted branch relation", lifted)]:
        vector = raw.astype(np.int64)
        assert np.array_equal(vector, raw) and int(vector.sum()) == 0
        for name in ("21", "35"):
            row = EXACT_TABLE[names.index(name)]
            assert all(isinstance(x, (int, sp.Integer)) for x in row)
            maximum = 2520 * max(abs(int(x)) for x in row) * max(abs(int(x)) for x in vector)
            assert maximum < 2**63  # Proven bound for each accumulated int64 coordinate.
            numerator = np.zeros(len(vector), dtype=np.int64)
            for g, perm in zip(A7, PERMS):
                image = np.empty_like(vector)
                image[perm] = vector
                numerator += int(row[idx[inv(g)]]) * image
            squared = sum(int(x)**2 for x in numerator)  # Python integers, no overflow.
            assert squared > 0
            norm2 = F(squared * int(row[idx[tuple(range(7))]])**2, 2520**2)
            print(f"PASS: {label}, exact {name}-projector squared norm = {norm2}")
    print("The unique 21/35 copies imply common generated modules; no torsion conclusion follows.")


if __name__ == "__main__":
    genus_arithmetic()
    finite_pencil_bounds()
    polynomial_bounds()
    exact_characters_and_divisors()
    print("PASS: textbook finite inputs; geometric and analytical trust bases remain explicit")
