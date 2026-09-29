"""Which triangle signatures (l,m,n) generate A7, and the resulting genus."""
from a7 import *
from fractions import Fraction as F
import random
random.seed(1)
by_order = {}
for p in A7:
    by_order.setdefault(order(p), []).append(p)
orders = sorted(by_order)
res = []
for l in orders:
    for m in orders:
        for n in orders:
            if not (2 <= l <= m <= n): continue
            chi = 1 - F(1, l) - F(1, m) - F(1, n)
            if chi <= 0: continue
            g = 1 + 2520 * chi / 2
            # search: fix x in each class of order l (representatives), y of order m, check xy order n^-1
            found = False
            reps = {a7_class(x) for x in by_order[l]}
            for cls in reps:
                x = next(iter(cls))
                for y in by_order[m]:
                    z = mul(x, y)
                    if order(z) == n and len(generated([x, y])) == 2520:
                        found = True; break
                if found: break
            res.append((g, (l, m, n), found))
for g, s, f in sorted(res):
    if g <= 140: print(s, "genus", g, "generates A7" if f else "-")
