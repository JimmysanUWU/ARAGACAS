"""Exact check that the twisted quotients Q0, Q1, Q2 see every irreducible representation of A7.

Q_i = (K_i, eps_i):  eps_i-twisted functions on C  <->  Hom_K(eps, .)  so an eigenspace E of the
Laplacian on C (an A7-module) contributes to Q_i iff it contains an irreducible rho with
    m_i(rho) = <rho|_{K_i}, eps_i> = (1/|K_i|) sum_{k in K_i} chi_rho(k) eps_i(k)  >  0.
K_1, K_2 contain no elements of order 7, so only the integer entries of the character table enter and
the computation below is exact (Fractions).  The table is the ATLAS table of A7 (also recomputed
numerically in chartab.py); its orthogonality relations are verified exactly here, and the Frobenius
count  sum_rho m_i(rho) dim(rho) = [A7 : K_i]  is checked as a consistency test.
"""
from fractions import Fraction
import sympy
from a7 import A7, cyc, cycle_type, generated, mul, E

# classes: 1, 2^2, 3, 3^2, 4.2, 5, 3.2^2, 7a, 7b   (sizes below)
CLASSES = [(1, 1, 1, 1, 1, 1, 1), (2, 2, 1, 1, 1), (3, 1, 1, 1, 1), (3, 3, 1), (4, 2, 1), (5, 1, 1),
           (3, 2, 2), "7a", "7b"]
SIZES = [1, 105, 70, 280, 630, 504, 210, 360, 360]
al = (-1 + sympy.sqrt(-7)) / 2
alb = (-1 - sympy.sqrt(-7)) / 2
TABLE = {
    "1":   [1, 1, 1, 1, 1, 1, 1, 1, 1],
    "6":   [6, 2, 3, 0, 0, 1, -1, -1, -1],
    "10":  [10, -2, 1, 1, 0, 0, 1, al, alb],
    "10b": [10, -2, 1, 1, 0, 0, 1, alb, al],
    "14a": [14, 2, 2, -1, 0, -1, 2, 0, 0],
    "14b": [14, 2, -1, 2, 0, -1, -1, 0, 0],
    "15":  [15, -1, 3, 0, -1, 0, -1, 1, 1],
    "21":  [21, 1, -3, 0, -1, 1, 1, 0, 0],
    "35":  [35, -1, -1, -1, 1, 0, -1, 0, 0],
}


def check_table():
    assert sum(SIZES) == 2520
    names = list(TABLE)
    for x in names:
        for y in names:
            s = sum(sz * TABLE[x][i] * sympy.conjugate(TABLE[y][i]) for i, sz in enumerate(SIZES))
            assert sympy.simplify(s - (2520 if x == y else 0)) == 0, (x, y)
    # class sizes against cycle types in A7 (7-cycles: 720 = 360 + 360)
    cnt = {}
    for g in A7:
        ct = cycle_type(g)
        cnt[ct] = cnt.get(ct, 0) + 1
    for ct, sz in zip(CLASSES[:7], SIZES[:7]):
        assert cnt[ct] == sz, ct
    assert cnt[(7,)] == 720


def character(gens, vals):
    K = generated(gens)
    eps = {E: 1}
    frontier = [E]
    while frontier:
        new = []
        for x in frontier:
            for g, v in zip(gens, vals):
                y = mul(g, x)
                e = v * eps[x]
                if y in eps:
                    assert eps[y] == e, "not a homomorphism"
                else:
                    eps[y] = e
                    new.append(y)
        frontier = new
    assert len(eps) == len(K)
    return eps


def multiplicities(eps):
    order = len(eps)
    out = {}
    for name, row in TABLE.items():
        s = Fraction(0)
        for k, e in eps.items():
            ct = cycle_type(k)
            assert ct != (7,), "K contains an element of order 7"
            s += Fraction(int(row[CLASSES.index(ct)]) * e)
        m = s / order
        assert m.denominator == 1 and m >= 0
        out[name] = int(m)
    return out


QUOTIENT_DATA = {
    "Q1": ([cyc((0, 1, 2)), cyc((3, 4, 5)), cyc((0, 3, 1, 4), (2, 5))], [1, 1, -1]),
    "Q2": ([(4, 3, 1, 6, 0, 5, 2), (0, 6, 2, 3, 5, 4, 1)], [1, -1]),
}

if __name__ == "__main__":
    check_table()
    print("character table: exact orthonormality and class sizes verified")
    seen = {"1"}          # Q0 = A7-invariant functions: the trivial representation
    print("Q0: K = A7, eps = 1 -> {1}")
    for q, (gens, vals) in QUOTIENT_DATA.items():
        eps = character(gens, vals)
        m = multiplicities(eps)
        dim = sum(m[x] * TABLE[x][0] for x in m)
        assert dim == 2520 // len(eps), "Frobenius count"
        print(f"{q}: |K| = {len(eps)}, index {2520 // len(eps)}, eps nontrivial: {any(v < 0 for v in eps.values())},"
              f" multiplicities {dict((x, v) for x, v in m.items() if v)}  (Frobenius: sum m dim = {dim})")
        seen |= {x for x, v in m.items() if v}
    assert seen == set(TABLE), set(TABLE) - seen
    print("every irreducible representation of A7 is seen by Q0, Q1 or Q2")
