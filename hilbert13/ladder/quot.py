"""Fixed points and quotient genera for subgroups, using Fix counts derived in rung1."""
from a7 import *
from fractions import Fraction as F
from itertools import combinations
GC = 136
def fixC(h):
    if h == E: return None
    o, ct = order(h), cycle_type(h)
    return {2: 18, 4: 2, 7: 3}.get(o, 0)

def genus_quot(K):
    """Riemann-Hurwitz: 2g(C)-2 = |K|(2g_K-2) + sum_{k!=1} Fix(k)."""
    s = sum(fixC(k) for k in K if k != E)
    val = F(2 * GC - 2 - s, len(K))
    return (val + 2) / 2

def fix_on_quotient(w, V):
    """# fixed points on C/V of w (w normalizes V), = (1/|V|) sum_{g in wV} Fix_C(g)."""
    return F(sum(fixC(mul(w, v)) for v in V), len(V))

if __name__ == "__main__":
    tau = cyc((0, 1), (2, 3))
    T = [E, tau]
    Ct = centralizer(tau)
    print("|C(tau)| =", len(Ct), " g(D) =", genus_quot(T))
    # elements of H = C(tau)/<tau>
    seen, H = set(), []
    for h in Ct:
        if h in seen: continue
        seen |= {h, mul(h, tau)}
        H.append(h)
    for h in H:
        if h == E: continue
        oh = order(h); htau = mul(h, tau)
        print(f"  lift {cycle_type(h)}/{cycle_type(htau)}  order on D: {oh if order(mul(h,h)) not in (1,2) or mul(h,h) not in T else (2 if mul(h,h) in T else oh)}"
              f"  Fix_D = {fix_on_quotient(h, T)}")
    # quotient genera of D by subgroups of H (as subgroups of C(tau) containing tau)
    subs = set()
    for x in Ct:
        for y in Ct:
            K = frozenset(generated([x, y, tau]))
            subs.add(K)
    print("subgroups K with tau in K <= C(tau):", len(subs))
    for K in sorted(subs, key=len):
        ords = sorted(order(k) for k in K)
        print(f"   |K|={len(K):3d}  |K/tau|={len(K)//2:2d}  g(C/K)={genus_quot(K)}  orders={ords}")
