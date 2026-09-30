"""Exact normal-kernel check for the order-72 pencil stabilizer and Q2 structure.

The maximal-subgroup list and classification of finite PGL2(C) subgroups are
classical inputs; this checks the two finite groups used in the continuation.
"""
from collections import Counter
from a7 import E,cyc,mul,inv,conj,order,generated,centralizer


def normal_subgroups(K):
    principal = {frozenset(generated([conj(g,h) for g in K])) for h in K}
    found = {frozenset([E])}
    todo = list(found)
    while todo:
        H = todo.pop()
        for N in principal:
            J = frozenset(generated(list(H|N)))
            if J not in found:
                found.add(J); todo.append(J)
    return found


def quotient_order(g,N):
    x = E
    for k in range(1,73):
        x = mul(x,g)
        if x in N:
            return k
    raise AssertionError("quotient order not found")


def main():
    H72 = generated([cyc((0,1,2)),cyc((0,1),(2,3)),cyc((4,5,6)),cyc((0,1),(4,5))])
    assert len(H72) == 72
    normals = normal_subgroups(H72)
    prime_to_three = [N for N in normals if len(N)%3]
    assert sorted(len(N) for N in prime_to_three) == [1,4]
    print("order-72 normal subgroup orders:",sorted(len(N) for N in normals))
    for N in sorted(prime_to_three,key=len):
        qo = len(H72)//len(N)
        max_order = max(quotient_order(g,N) for g in H72)
        assert qo in [72,18] and max_order < qo//2
        print("kernel order",len(N),"quotient order",qo,"max element order",max_order,
              "excludes cyclic and dihedral target")
    C3 = generated([cyc((4,5,6))])
    assert C3 in normals
    counts = Counter(quotient_order(g,C3) for g in H72)
    assert counts == {1:3,2:27,3:24,4:18}
    print("normal C3 quotient has S4 orders:",{o:n//3 for o,n in sorted(counts.items())})
    Q2 = generated([(4,3,1,6,0,5,2),(0,6,2,3,5,4,1)])
    tau = cyc((1,6),(2,3))
    assert Q2 == set(centralizer(tau)) and len(Q2)==24
    counts = Counter(order(g) for g in Q2)
    assert counts == {1:1,2:9,3:2,4:6,6:6}
    centre = [g for g in Q2 if all(mul(g,h)==mul(h,g) for h in Q2)]
    assert set(centre) == {E,tau}
    N3 = {E}|{g for g in Q2 if order(g)==3}
    assert all(conj(g,x) in N3 for g in Q2 for x in N3)
    candidates = [generated([x,y]) for x in Q2 if order(x)==4 for y in Q2 if order(y)==2]
    D8 = next(H for H in candidates if len(H)==8 and Counter(order(g) for g in H)=={1:1,2:5,4:2})
    assert len(D8&N3)==1 and {mul(n,d) for n in N3 for d in D8} == Q2
    kernel = [d for d in D8 if all(mul(d,n)==mul(n,d) for n in N3)]
    assert Counter(order(g) for g in kernel)=={1:1,2:3}
    print("Q2 order",len(Q2),"element orders",dict(sorted(counts.items())))
    print("Q2 = C3 semidirect D8; action kernel is V4; centre is <tau>")


if __name__ == "__main__":
    main()
