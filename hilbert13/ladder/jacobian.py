"""Isotypic decomposition of H^1(C) (Chevalley-Weil, convention-free part) and of Jac(D)."""
import numpy as np
from a7 import *
from chartab import TABLE, CL, idx
names = ["1", "6", "10", "10b", "14a", "14b", "15", "21", "35"]
a = cyc((0, 1), (2, 3))
from triples_data import triples
_, b, c = triples[0]
def fixdim(chi, x):
    m = order(x)
    return round(sum(chi[idx[power(x, j)]] for j in range(m)).real / m)
tau = a
print(f"{'rho':>4} {'dim':>4} {'a':>3} {'b':>3} {'c':>3} {'H1':>4} {'mu=H0(K)':>8} {'dim rho^tau':>11} {'contrib to g(D)':>15}")
tot, totD = 0, 0
for n, chi in zip(names, TABLE):
    d = round(chi[0].real)
    fa, fb, fc = fixdim(chi, a), fixdim(chi, b), fixdim(chi, c)
    h1 = 0 if n == "1" else d - fa - fb - fc
    mu = h1 / 2
    ft = fixdim(chi, tau)
    tot += d * mu; totD += ft * mu
    print(f"{n:>4} {d:>4} {fa:>3} {fb:>3} {fc:>3} {h1:>4} {mu:>8} {ft:>11} {ft*mu:>15}")
print("g(C) =", tot, " g(D) =", totD)
