"""Chapter 1: the finite facts about the (2,4,7) curves (output: curve_checks_output.txt, ~1 min).

1. triples: counts and the four A7-classes 0, 1, 12, 14 (two S7-classes)
2. fixed points over the branch points, by conjugacy class (Lemma 1.1, Corollary 1.2)
3. D = C/tau: the action of H = C(tau)/tau and the 16 quotient genera (Proposition 1.3)
4. Chevalley-Weil: H^1(C), its tau-invariants (Jac D), and the H-isotypes of Jac(D)
5. the divisor behind (*): isotypic components of Phi_g - 3 Phi_c (Proposition 1.6)
The branch-point divisor tools in part 5 (fix_divisor, isotypic projections) are also used by side_checks.py.
"""
from fractions import Fraction as F
import numpy as np
from a7 import *
from chartab import TABLE, idx
from triples_data import triples

names = ["1", "6", "10", "10b", "14a", "14b", "15", "21", "35"]
a, b, c = triples[0]
FIX = {2: 18, 4: 2, 7: 3}


# ---------------------------------------------------------------- divisors on the branch fibres of C -> C/A7
def coset_space(x):
    X = [power(x, k) for k in range(order(x))]
    reps, where = [], {}
    for g in A7:
        if g in where: continue
        for y in (mul(g, h) for h in X): where[y] = len(reps)
        reps.append(g)
    return reps, where

SP = [coset_space(x) for x in (a, b, c)]
offs = np.cumsum([0] + [len(r) for r, _ in SP])
NPTS = offs[-1]

def act_perm(g):
    P = np.empty(NPTS, dtype=np.int64)
    for i, (reps, where) in enumerate(SP):
        for j, r in enumerate(reps): P[offs[i] + j] = offs[i] + where[mul(g, r)]
    return P

PERMS = [act_perm(g) for g in A7]
INV_CLASS = [idx[inv(g)] for g in A7]

def fix_divisor(g):
    P = act_perm(g); v = np.zeros(NPTS); v[P == np.arange(NPTS)] = 1; return v

def isotypic_norms(v):
    out = {}
    for n, chi in zip(names, TABLE):
        w = np.zeros(NPTS, dtype=complex)
        for P, ic in zip(PERMS, INV_CLASS):
            gv = np.empty(NPTS); gv[P] = v; w += chi[ic] * gv
        out[n] = round(float(np.linalg.norm(w * chi[0].real / 2520)), 4)
    return out


def fixdim(chi, x):
    m = order(x); return round(sum(chi[idx[power(x, j)]] for j in range(m)).real / m)


def main():
    print("1. Triples")
    inv2 = [p for p in A7 if order(p) == 2]; ord4 = [p for p in A7 if order(p) == 4]
    trips = [bb for bb in ord4 if order(inv(mul(a, bb))) == 7 and len(generated([a, bb])) == 2520]
    print(f"   involutions {len(inv2)}, order 4: {len(ord4)}; with a = (01)(23) fixed: {len(trips)} choices of b;"
          f" {len(trips) * 105} triples = {len(trips) * 105 // 2520} A7-classes = {len(trips) * 105 // 5040} S7-classes")
    Ca, seen, reps = centralizer(a), set(), []
    for k, (_, bb, _) in enumerate(triples):
        if bb not in seen: seen |= {conj(g, bb) for g in Ca}; reps.append(k)
    print("   class representatives (indices into triples_data):", reps)

    print("2. Fixed points over the three branch points (Lemma 1.1); genus", 1 + F(2520, 2) * (-2 + F(1, 2) + F(3, 4) + F(6, 7)))
    def fix_on_fiber(h, x):
        X = {power(x, k) for k in range(order(x))}
        return sum(1 for g in A7 if conj(inv(g), h) in X) // len(X)
    cls7 = a7_class(cyc(tuple(range(7)))); creps = {}
    for p in A7: creps.setdefault((cycle_type(p), p in cls7 if order(p) == 7 else None), p)
    for key, h in sorted(creps.items(), key=lambda kv: order(kv[1])):
        if h != E:
            fx = [fix_on_fiber(h, x) for x in (a, b, c)]
            print(f"   {key[0]}{'' if key[1] is None else (' class 1' if key[1] else ' class 2')}: {fx}, total {sum(fx)}")

    print("3. D = C/tau")
    tau = a; T = [E, tau]; Ct = centralizer(tau)
    def genus_quot(K): return (F(270 - sum(FIX.get(order(k), 0) for k in K if k != E), len(K)) + 2) / 2
    print(f"   |C(tau)| = {len(Ct)}, g(D) = {genus_quot(T)}")
    seen = set()
    for h in Ct:
        if h in seen or h == E or h == tau: continue
        seen |= {h, mul(h, tau)}
        fixD = F(sum(FIX.get(order(mul(h, v)), 0) for v in T), 2)
        oD = order(h) if mul(h, h) not in T else (1 if h in T else 2)
        print(f"   lift {cycle_type(h)}: order {oD} on D, {fixD} fixed points")
    subs = {frozenset(generated([x, y, tau])) for x in Ct for y in Ct}
    print(f"   {len(subs)} subgroups K with tau in K <= C(tau); genera g(C/K):",
          sorted((len(K) // 2, int(genus_quot(K))) for K in subs))

    print("4. Chevalley-Weil (multiplicity of rho in H^1(C) = dim rho - sum of fixed dims at the branches)")
    gC = gD = 0
    for n, chi in zip(names, TABLE):
        d = round(chi[0].real); h1 = 0 if n == "1" else d - sum(fixdim(chi, x) for x in (a, b, c)); ft = fixdim(chi, tau)
        gC += d * h1 / 2; gD += ft * h1 / 2
        if h1: print(f"   {n:>4}: H^1 multiplicity {h1}, dim rho^tau {ft}")
    print(f"   g(C) = {gC:.0f}, g(D) = {gD:.0f}")
    cbar = cyc((0, 2), (1, 3))
    def h_coords(g):
        y = tuple(g[i] - 4 for i in (4, 5, 6)); sy = list(range(7))
        for i in range(3): sy[4 + i] = 4 + y[i]
        sy = tuple(sy)
        if sign(sy) == -1: sy = mul(cyc((0, 1)), sy)
        rest = mul(g, inv(sy)); assert rest in (E, tau, cbar, mul(tau, cbar))
        return (0 if rest in (E, tau) else 1), y
    def s3(y, w):
        fixed = sum(1 for i in range(3) if y[i] == i)
        return 1 if w == "1" else ((1 if fixed in (3, 0) else -1) if w == "sgn" else {3: 2, 1: 0, 0: -1}[fixed])
    print("   H-isotypes of rho^tau (H = <e> x S3):")
    for n, chi in zip(names, TABLE):
        out = {}
        for i in (0, 1):
            for w in ("1", "sgn", "std"):
                s = sum((chi[idx[g]] + chi[idx[mul(g, tau)]]) / 2 * (-1) ** (i * h_coords(g)[0]) * s3(h_coords(g)[1], w) for g in Ct)
                if round((s / 24).real): out[("1", "e")[i] + "," + w] = round((s / 24).real)
        print(f"   {n:>4}: {out}")

    print("5. The divisor behind (*)")
    mu = cyc((0, 2), (1, 3)); goods = [cyc((0, 1), ij) for ij in [(4, 5), (4, 6), (5, 6)]]
    goods += [mul(tau, g) for g in goods]
    v = sum(fix_divisor(g) for g in goods) - 3 * (fix_divisor(mu) + fix_divisor(mul(tau, mu)))
    print("   isotypic norms of Phi_g - 3 Phi_c:", {k: x for k, x in isotypic_norms(v).items() if x})


if __name__ == "__main__":
    main()
