from a7 import *
from itertools import permutations
_a = cyc((0, 1), (2, 3))
triples = []
_seen = set()
for b in [p for p in A7 if order(p) == 4]:
    c = inv(mul(_a, b))
    if order(c) == 7 and len(generated([_a, b])) == 2520:
        triples.append((_a, b, c))
