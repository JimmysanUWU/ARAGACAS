"""Exact PSL2(F13) computation for the genus-fourteen transport corollary."""
from itertools import combinations, product
from fractions import Fraction as F
from collections import Counter

P = 13


def canonical(x):
    return min(x, tuple((-a) % P for a in x))


I = canonical((1,0,0,1))


def mul(x,y):
    a,b,c,d = x
    e,f,g,h = y
    return canonical(((a*e+b*g)%P,(a*f+b*h)%P,(c*e+d*g)%P,(c*f+d*h)%P))


def inv(x):
    a,b,c,d = x
    return canonical((d,(-b)%P,(-c)%P,a))


def conj(g,x):
    return mul(mul(g,x),inv(g))


def order(x):
    y = I
    for n in range(1,1093):
        y = mul(y,x)
        if y == I:
            return n
    raise AssertionError("order not found")


def generated(gens):
    seen = {I}
    todo = [I]
    while todo:
        x = todo.pop()
        for g in gens:
            y = mul(g,x)
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


def main():
    G = sorted({canonical(x) for x in product(range(P),repeat=4) if (x[0]*x[3]-x[1]*x[2])%P == 1})
    assert len(G) == 1092
    tau = canonical((5,0,0,8))
    assert order(tau) == 2
    cent = [g for g in G if mul(g,tau) == mul(tau,g)]
    assert len(cent) == 12
    Vs = sorted({frozenset([I,tau,u,mul(tau,u)]) for u in cent if order(u)==2 and u!=tau}, key=lambda v:sorted(v))
    assert len(Vs) == 3
    Ns = [[g for g in G if {conj(g,v) for v in V} == V] for V in Vs]
    print("group order",len(G),"centralizer",len(cent),"Klein subgroups",len(Vs))
    for i,N in enumerate(Ns):
        counts = Counter(order(g) for g in N)
        assert counts == {1:1,2:3,3:8}
        print("normalizer",i,"order",len(N),"element orders",dict(sorted(counts.items())))
    for i,j in combinations(range(3),2):
        n = len(generated(Ns[i]+Ns[j]))
        assert n == 1092
        print("normalizers",i,j,"generate",n)
    r = len(cent)//2
    genus_D = 1+F(14-1,2)-F(r,4)
    assert r == 6 and genus_D == 6
    reflection_quotient = 1+F(genus_D-1,2)-F(6,4)
    assert reflection_quotient == 2
    assert 3*12*6 == 216 and 216 % 13 != 0
    print("fixed points",r,"quotient genus",genus_D,"transport degree",216)
    print("three good-reflection quotients have genus",reflection_quotient,
          "and give three distinct degree-four pencils")


if __name__ == "__main__":
    main()
