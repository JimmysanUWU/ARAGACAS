"""Exact generating-triple enumeration through genus 335, without SVD inputs."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from collections import Counter
from sympy.combinatorics import Permutation, PermutationGroup
from a7 import A7,order,cycle_type,inv,mul


EXPECTED = {(2,4,7):(136,4),(3,3,5):(169,2),(2,5,7):(199,4),
            (3,3,6):(211,2),(3,4,4):(211,8),(2,6,7):(241,4),(3,3,7):(241,4),
            (2,7,7):(271,6),(3,4,5):(274,10),(3,4,6):(316,6),(4,4,4):(316,24)}


def main():
    by_order = {}
    representatives = {}
    sizes = Counter(cycle_type(g) for g in A7)
    for g in A7:
        by_order.setdefault(order(g),[]).append(g)
        representatives.setdefault(cycle_type(g),g)
    found = {}
    for p,q,r in combinations_with_replacement([2,3,4,5,6,7],3):
        chi = 1-F(1,p)-F(1,q)-F(1,r)
        genus = 1+1260*chi
        if chi <= 0 or genus.denominator != 1 or genus > 335:
            continue
        total = 0
        for ct,x in representatives.items():
            if order(x) != p:
                continue
            # A first entry of order seven never occurs in this genus range;
            # the only split S7 cycle class therefore needs no separation here.
            assert p != 7
            for y in by_order[q]:
                z = inv(mul(x,y))
                if order(z) != r or 21-len(ct)-len(cycle_type(y))-len(cycle_type(z)) < 12:
                    continue
                if PermutationGroup([Permutation(list(x)),Permutation(list(y))]).order() == 2520:
                    total += sizes[ct]
        if total:
            assert total % 2520 == 0
            found[p,q,r] = (int(genus),total//2520)
    assert found == EXPECTED,found
    for sig,(g,n) in sorted(found.items(),key=lambda row:(row[1][0],row[0])):
        print(sig,"genus",g,"inner classes",n,"strict spectral threshold",F(48,g-1),flush=True)
    max_index = {o:max(7-len(cycle_type(g)) for g in by_order[o]) for o in by_order}
    eliminated = []
    for sig in combinations_with_replacement([2,3,4,5,6,7],4):
        chi = 2-sum(F(1,e) for e in sig)
        genus = 1+1260*chi
        if 1 < genus < 421:
            assert sum(max_index[o] for o in sig) < 12
            eliminated.append(sig)
    print("four-point candidates below 421 fail the degree-seven RH test:",eliminated)
    print("five-point genus >=",1+1260*F(1,2))
    print("positive quotient genus with branching has genus >=",1+1260*F(1,2))
    assert min(g for g,n in found.values()) == 136


if __name__ == "__main__":
    main()
