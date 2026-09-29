"""Count (2,4,7) generating triples of A7 up to simultaneous conjugation in S7."""
from a7 import *
from itertools import permutations
S7 = list(permutations(range(N)))
inv2 = [p for p in A7 if order(p) == 2]
ord4 = [p for p in A7 if order(p) == 4]
print("involutions", len(inv2), "order4", len(ord4), set(cycle_type(p) for p in ord4))
# fix a representative a of the involution class, enumerate b, set c = (ab)^-1
a = cyc((0, 1), (2, 3))
Ca = centralizer(a)
trips = []
for b in ord4:
    c = inv(mul(a, b))
    if order(c) == 7 and len(generated([a, b])) == 2520:
        trips.append((a, b, c))
print("triples with a fixed:", len(trips), " |C(a)| =", len(Ca))
# total triples = len(trips) * 105 ; orbits under A7-conj = total/2520 ; under S7 = /5040 if free
tot = len(trips) * 105
print("total triples", tot, "A7-orbits", tot / 2520, "S7-orbits", tot / 5040)
# which A7-class of 7-cycle does c lie in
cls = a7_class(cyc(tuple(range(7))))
print("c in class of (0123456):", sum(1 for t in trips if t[2] in cls), "of", len(trips))
