"""Finite checks for the review of the 'Ramification transport' synthesis (REVIEW_GPT.md)."""
from a7 import *
c = lambda *cy: cyc(*[tuple(x - 1 for x in t) for t in cy])


def check_transport():
    tau = c((1, 2), (3, 4))
    V0 = [E, tau, c((1, 3), (2, 4)), c((1, 4), (2, 3))]
    V1 = [E, tau, c((1, 2), (5, 6)), c((3, 4), (5, 6))]
    def normalizer(V):
        S = set(V); return [g for g in A7 if {conj(g, v) for v in V} == S]
    N0, N1 = normalizer(V0), normalizer(V1)
    print("|C(tau)| =", len(centralizer(tau)), " |N(V0)| =", len(N0), " |N(V1)| =", len(N1),
          " |<N(V0),N(V1)>| =", len(generated(N0 + N1)))
    z, gam, sig = c((1, 3), (2, 4)), c((5, 6, 7)), c((3, 4), (5, 7))
    print("<z, sigma, gamma sigma, tau> =", len(generated([z, sig, mul(gam, sig), tau])), "= |C(tau)|; tau*x involution for all three:",
          all(cycle_type(mul(tau, x)) == (2, 2, 1, 1, 1) for x in (z, sig, mul(gam, sig))))


def table1():
    # generating triples (x, y, (xy)^-1) of A7 with orders (a, b, c), counted up to conjugation of x
    reps = {}
    for g in A7:
        o = order(g)
        reps.setdefault((o, cycle_type(g)), g)
    for sig in [(2, 4, 7), (3, 3, 5), (2, 5, 7), (3, 3, 6), (3, 4, 4), (2, 6, 7), (3, 3, 7), (2, 7, 7), (3, 4, 5)]:
        a, b, cc = sig
        found = 0
        for (o, ct), x in reps.items():
            if o != a:
                continue
            for y in A7:
                if order(y) == b and order(mul(x, y)) == cc and len(generated([x, y])) == 2520:
                    found += 1
        g = 1 + 1260 * (-2 + sum(1 - 1 / m for m in sig))
        print(f"signature {sig}: genus {g:.0f}, generating pairs with x a fixed class rep: {found}")


def a5_example():
    A5 = [p for p in A7 if p[5] == 5 and p[6] == 6]
    t1, t2, t3, t4 = c((1, 2), (3, 4)), c((1, 3), (2, 4)), c((1, 5), (2, 3)), c((1, 4, 5))
    prod = mul(mul(mul(t1, t2), t3), t4)
    print("A5 example: product", "= 1" if prod == E else "!= 1 (other order?)", " generated order",
          len(generated([t1, t2, t3, t4])))
    for perm in [(t4, t3, t2, t1)]:
        p = E
        for x in perm: p = mul(p, x)
        print("   reversed product = 1:", p == E)


if __name__ == "__main__":
    check_transport()
    a5_example()
    table1()
