from a7 import *
from triples_data import triples
a = triples[0][0]
Ca = centralizer(a)
seen, reps = set(), []
for k, (aa, bb, cc) in enumerate(triples):
    if bb in seen: continue
    orb = {conj(g, bb) for g in Ca}
    seen |= orb
    reps.append(k)
print("A7-class representatives (indices into triples):", reps)
