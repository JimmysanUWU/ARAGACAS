"""H-isotypic structure of Jac(D): for each A7-irrep rho, decompose rho^tau as a rep of H=C(tau)/<tau>."""
import numpy as np
from a7 import *
from chartab import TABLE, idx
names = ["1", "6", "10", "10b", "14a", "14b", "15", "21", "35"]
mu_rho = {"10": 1.5, "10b": 1.5, "15": 1, "21": 1, "35": 2}
tau = cyc((0, 1), (2, 3))
Ct = centralizer(tau)
c_bar = cyc((0, 2), (1, 3))            # lifts the central good involution
# H = <c_bar> x S3 with goods = (1,t), bads = (c,t).  Identify via: sign on {5,6,7} part and c-part.
def h_coords(g):
    """return (x in {0,1}, y in S3 as perm of {4,5,6}) with g = c^x * s(y) mod tau."""
    y = tuple(g[i] - 4 for i in (4, 5, 6))
    # s(y) = transposition-type lift: goods are (12)(5x): lift of y odd -> (0 1)*y ; y even -> y
    sy = list(range(7))
    for i in range(3): sy[4 + i] = 4 + y[i]
    sy = tuple(sy)
    if sign(sy) == -1: sy = mul(cyc((0, 1)), sy)
    rest = mul(g, inv(sy))  # in D4 part, even... lies in <tau, c_bar>
    x = 0 if rest in (E, tau) else 1
    assert rest in (E, tau, c_bar, mul(tau, c_bar)), rest
    return x, y
# H irreps: eps^i (x) {triv, sgn, std}
def s3_char(y, which):
    ct = cycle_type(tuple(list(y) + [3, 4, 5, 6])[:3]) if False else None
    fixed = sum(1 for i in range(3) if y[i] == i)
    if which == "triv": return 1
    if which == "sgn": return 1 if fixed in (3, 0) else -1
    if which == "std": return {3: 2, 1: 0, 0: -1}[fixed]
Hirr = [(i, w) for i in (0, 1) for w in ("triv", "sgn", "std")]
lab = {(0, "triv"): "(1,1)", (0, "sgn"): "(1,sgn)", (0, "std"): "(1,std)",
       (1, "triv"): "(e,1)", (1, "sgn"): "(e,sgn)", (1, "std"): "(e,std)"}
for n, chi in zip(names, TABLE):
    # character of rho^tau on h: (chi(h) + chi(h tau))/2 ; inner product over C(tau) (order 24) -> over H
    out = {}
    for (i, w) in Hirr:
        s = 0
        for g in Ct:
            x, y = h_coords(g)
            val = (chi[idx[g]] + chi[idx[mul(g, tau)]]) / 2
            s += val * ((-1) ** (i * x)) * s3_char(y, w)
        m = (s / 24).real  # each H-element counted twice (g, g tau) -> divide by 24 = 2*12
        out[lab[(i, w)]] = round(m)
    print(f"{n:>4}: rho^tau|H =", {k: v for k, v in out.items() if v}, " mu =", mu_rho.get(n, 0))
