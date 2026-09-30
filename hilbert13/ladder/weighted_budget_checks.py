"""Exact arithmetic supporting GPT_CONTINUATION.md, Sections 1--5.

This verifies numerical inequalities; it does not replace the geometric proofs.
"""
from fractions import Fraction as F
from math import comb


def pi(d, n, alpha):
    q, eps = divmod(d - 1, n - 1 + alpha)
    mu = max(0, (alpha - n + 2 + eps) // 2)
    return (n - 1 + alpha) * comb(q, 2) + q * (eps + alpha) + mu


def A(m):
    return F(m*m, 2) - m + 1


def pair_bound(m):
    vals = [m*m // e + e*m - 3*m + 1 for e in range(2, m // 2 + 1) if m % e == 0]
    return max(vals, default=None)


def dep_costs(m):
    one = min(comb(r, 2) + comb(m-r, 2) for r in range(1, m))
    vee = min(comb(r, 2) + 2*comb(m-r, 2) for r in range(1, m))
    star = min(comb(r, 2) + 3*comb(m-r, 2) for r in range(1, m))
    path = min(2*comb(r, 2) + 2*comb(m-r, 2) for r in range(1, m))
    triangle = 3*comb(m // 2, 2) if m % 2 == 0 else 10**9
    return one, min(vee, 2*one), min(star, triangle, path, vee+one, 3*one)


def h19_genus(m):
    d = 3*m
    total = 0
    i = 1
    while True:
        h = min(d, 18*(i // 2) + (7 if i % 2 else 1))
        if h == d:
            return total
        total += d-h
        i += 1


def main():
    print("Claude genus-336 table: m pair pi1 budget dependent_costs")
    expected = [(181,228,25,(90,120,135)), (148,253,64,(100,133,150)),
                (221,279,105,(110,147,165)), (None,306,148,(121,161,181)),
                (265,335,193,(132,176,198))]
    for m, exp in zip(range(20, 25), expected):
        row = pair_bound(m), pi(3*m,7,1), (m-1)**2-336, dep_costs(m)
        assert row == exp
        print(m, *row)
    print("\nIndependent-triple inequalities (8<=m<=999, plus exact residue formulas)")
    for m in range(8, 1000):
        k, r = divmod(m, 3)
        assert A(m)-pi(3*m,8,2) == F(k+r,2)
        k, r = divmod(m, 6)
        assert int(A(m))-h19_genus(m) == 2*k+(r == 2)
        assert F(7*m*m,16)-F(m,2)+1 <= A(m)
        assert 1+F(3*m*m,8) <= A(m)
    print("all residue formulas and surface comparisons passed")
    print("m=24: A, pi2, H19, type221, type222:",
          A(24),pi(72,8,2),h19_genus(24),F(7*24**2,16)-12+1,1+F(3*24**2,8))
    print("\nGenus-266 star possibilities: m budget minimum_pair_cost allowed_r")
    for m in range(18, 25):
        assert pair_bound(m) is None or pair_bound(m) <= 265
        assert A(m) <= 265
        delta = (m-1)**2-266
        one = dep_costs(m)[0]
        assert delta < 2*one
        allowed = [r for r in range(1,m) if comb(r,2)+33*comb(m-r,2) <= delta]
        assert allowed == ([23] if m == 24 else [])
        print(m,delta,one,allowed)
    assert 2*comb(23,2) > (24-1)**2-266
    assert (25-1)*(25-2)//2 == 276 and 276 % 3 != 1
    assert all(1260 % (3*e) == 0 for e in [2,3,4,5,6,7])
    print("plane degree 25 has genus 276, forbidden by g=1 mod 3")
    print("\nAlgebraic lower bounds in the remaining genera")
    for genus, lower in [(136,18),(169,20),(199,21),(211,22),(241,23)]:
        for m in range(13,lower):
            if genus > (m-1)**2:
                continue
            assert pair_bound(m) is None or genus > pair_bound(m)
            assert genus > A(m)
            delta = (m-1)**2-genus
            assert delta < 2*dep_costs(m)[0]
            edges = 33 if m % 3 == 0 else 40
            allowed = [r for r in range(1,m) if comb(r,2)+edges*comb(m-r,2) <= delta]
            for r in allowed:
                assert r == m-1
                d = 2*m-r
                plane_genus = (d-1)*(d-2)//2
                assert plane_genus != genus or (genus,m,d) == (136,17,18)
                if plane_genus == genus:
                    assert 6*d*d < 2520  # Harui excludes the smooth plane model.
        print(genus,lower)
    print("\nSeven signatures and exact spectral thresholds")
    rows = [((2,4,7),136,4),((3,3,5),169,2),((2,5,7),199,4),
            ((3,3,6),211,2),((3,4,4),211,8),((2,6,7),241,4),((3,3,7),241,4)]
    for sig,g,n in rows:
        assert 1+1260*(1-sum(F(1,e) for e in sig)) == g
        print(sig,g,n,F(48,g-1))
    assert sum(row[2] for row in rows) == 28


if __name__ == "__main__":
    main()
