"""Audit 5.1: an explicit smooth H-invariant (9,9) curve with the same fixed-point data as D.

H = S3 x C2 acts on P1 x P1: S3 diagonally via a3: z -> w z, a2: z -> 1/z, and C2 by the swap.
F = sum c_ij x0^i x1^(9-i) y0^j y1^(9-j) with c_ij in Z, supported on i+j = 0 mod 3 (a3-invariant),
c_ij = c_{9-i,9-j} (a2-invariant), c_ij = c_ji (swap-invariant). All checks are done mod a prime p;
smoothness / squarefreeness / non-vanishing mod p imply the same over Q (open conditions).
"""
import random
from sympy import symbols, Poly, groebner, gcd, diff, GF, factor_list, isprime

x, y = symbols("x y")
p = 1000003
assert isprime(p)
random.seed(20260929)

def build():
    c = {}
    for i in range(10):
        for j in range(10):
            if (i + j) % 3: continue
            key = min((i, j), (j, i), (9 - i, 9 - j), (9 - j, 9 - i))
            if key not in c:
                c[key] = random.randrange(1, 50)
            c[(i, j)] = c[key]
    return {k: v for k, v in c.items() if len(k) == 2}

C = build()
# affine chart x1 = y1 = 1:  f(x, y) = F(x, 1, y, 1)
f = sum(v * x**i * y**j for (i, j), v in C.items())
fP = Poly(f, x, y, modulus=p)

def squarefree_binary(poly_in_x, deg):
    """binary form of degree deg given by its dehomogenization; squarefree incl. the point at infinity"""
    P = Poly(poly_in_x, x, modulus=p)
    if P.degree() < deg - 1:  # double root at infinity
        return False
    g = gcd(P, P.diff(x))
    return g.degree() == 0

report = {}
# 1. smoothness in the chart (other charts follow by symmetry, see notes)
G = groebner([fP.as_expr(), diff(f, x), diff(f, y)], x, y, modulus=p, order="grevlex")
report["singular points in chart"] = "none" if list(G.exprs) == [1] else G.exprs
# 2. transversal to the diagonal: F(x,x) squarefree, degree 18
report["D.Delta squarefree (18 pts)"] = squarefree_binary(f.subs(y, x), 18)
# 3. transversal to the graph of a2 (y = 1/x): x^9 * F(x,1,1,x)
g2 = sum(v * x**i * x**(9 - j) for (i, j), v in C.items())
report["D.Gamma(a2) squarefree (18 pts)"] = squarefree_binary(g2, 18)
# 4. bad involution (a2,a2): fixed pts (+-1, +-1); D contains exactly (1,-1), (-1,1)
val = lambda X, Y: sum(v * X**i * Y**j for (i, j), v in C.items()) % p
report["F(1,1), F(-1,-1), F(1,-1), F(-1,1) mod p"] = (val(1, 1), val(-1, -1), val(1, -1), val(-1, 1))
# 5. C3 fixed points (0,0),(0,inf),(inf,0),(inf,inf) avoided: corner coefficients nonzero
report["corner coeffs c00,c09,c90,c99"] = (C[(0, 0)], C[(0, 9)], C[(9, 0)], C[(9, 9)])
for k, v in report.items():
    print(f"{k}: {v}")
