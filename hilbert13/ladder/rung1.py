"""Rung 1: fixed points computed from branch orbits for every generating triple class."""
from a7 import *
from triples_data import triples  # noqa
from fractions import Fraction as F

def fix_on_fiber(h, x):
    """# cosets g<x> fixed by h  = # {g : g^-1 h g in <x>} / |<x>|"""
    X = {power(x, k) for k in range(order(x))}
    return sum(1 for g in A7 if conj(inv(g), h) in X) // len(X)

reps = {}
for p in A7:
    key = (cycle_type(p), p in a7_class(cyc(tuple(range(7)))) if order(p) == 7 else None)
    reps.setdefault(key, p)
for (a, b, c) in triples:
    g2 = 1 + F(2520, 2) * (-2 + F(1, 2) + F(3, 4) + F(6, 7))
    print("genus C =", g2)
    for key, h in sorted(reps.items(), key=lambda kv: order(kv[1])):
        if h == E:
            continue
        fx = [fix_on_fiber(h, x) for x in (a, b, c)]
        print(f"  class {key}: fixed points over the three branch points {fx}, total {sum(fx)}")
    break
