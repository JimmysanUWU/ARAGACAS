"""Rational Specht-character check of the two fixed-divisor isotypes.

Unlike frontier_checks.py this uses integer Murnaghan--Nakayama characters,
integer conjugation permutations of the 105 involutions, and Fraction norms.
It proves the representation-theoretic assertions, not non-torsion.
"""
from functools import lru_cache
from fractions import Fraction as F
from collections import Counter
import numpy as np
from a7 import A7,E,cyc,mul,inv,conj,cycle_type,order,centralizer,generated


@lru_cache(None)
def partitions(n,limit=7):
    if n == 0:
        return ((),)
    return tuple((a,)+b for a in range(min(n,limit),0,-1) for b in partitions(n-a,a))


def cells(lam):
    return {(i,j) for i,n in enumerate(lam) for j in range(n)}


@lru_cache(None)
def character(lam,ct):
    if not ct:
        return int(not lam)
    k = ct[0]
    diagram = cells(lam)
    ans = 0
    for mu in partitions(sum(lam)-k):
        small = cells(mu)
        if not small <= diagram:
            continue
        strip = diagram-small
        if not strip:
            continue
        if any({(i,j),(i+1,j),(i,j+1),(i+1,j+1)} <= strip for i,j in strip):
            continue
        reached = {next(iter(strip))}
        todo = list(reached)
        while todo:
            i,j = todo.pop()
            for x in [(i-1,j),(i+1,j),(i,j-1),(i,j+1)]:
                if x in strip and x not in reached:
                    reached.add(x); todo.append(x)
        if reached != strip:
            continue
        height = len({i for i,j in strip})
        ans += (-1)**(height-1)*character(mu,ct[1:])
    return ans


REPS = {"1":(7,),"6":(6,1),"14(5,2)":(5,2),"14(4,3)":(4,3),
        "15":(5,1,1),"21":(3,3,1),"35":(4,2,1),"20":(4,1,1,1)}


def main():
    invols = [g for g in A7 if order(g)==2]
    assert len(invols) == 105
    index = {g:i for i,g in enumerate(invols)}
    gidx = {g:i for i,g in enumerate(A7)}
    action = [np.array([index[conj(g,t)] for t in invols]) for g in A7]
    chars = {n:np.array([character(lam,cycle_type(g)) for g in A7],dtype=np.int64) for n,lam in REPS.items()}
    dims = {n:character(lam,(1,)*7) for n,lam in REPS.items()}
    assert dims == {"1":1,"6":6,"14(5,2)":14,"14(4,3)":14,"15":15,"21":21,"35":35,"20":20}
    tau = cyc((0,1),(2,3))
    mu = cyc((0,2),(1,3))
    nu = cyc((0,1),(4,5))
    def vector(items):
        v = np.zeros(105,dtype=np.int64)
        for t,c in items:
            v[index[t]] += c
        return v
    D = vector([(mu,1),(mul(tau,mu),1),(nu,-1),(mul(tau,nu),-1)])
    good = [cyc((0,1),(i,j)) for i,j in [(4,5),(4,6),(5,6)]]
    good += [mul(tau,x) for x in good]
    star = vector([(g,1) for g in good]+[(mu,-3),(mul(tau,mu),-3)])
    def project_numerator(v,name):
        out = np.zeros(105,dtype=np.int64)
        for perm,chi in zip(action,chars[name]):
            out[perm] += int(chi)*v
        return out
    projections = {n:project_numerator(D,n) for n in REPS}
    norms = {n:18*F(dims[n],2520)**2*int(v@v) for n,v in projections.items()}
    expected = {"1":F(0),"6":F(36,5),"14(5,2)":F(144,5),"14(4,3)":F(12),
                "15":F(0),"21":F(18),"35":F(6),"20":F(0)}
    assert norms == expected, norms
    cent = centralizer(tau)
    def fixdims(K):
        result = {}
        for n,chi in chars.items():
            total = sum(int(chi[gidx[g]]) for g in K)
            assert total % len(K) == 0
            result[n] = total // len(K)
        return result
    source_dims = fixdims(cent)
    assert source_dims == {"1":1,"6":1,"14(5,2)":2,"14(4,3)":1,"15":0,"21":1,"35":1,"20":0}
    print("source multiplicities:",source_dims)
    print("Klein squared norms:",{n:str(x) for n,x in norms.items() if x})
    for n in ["21","35"]:
        assert np.any(project_numerator(star,n))
        print("old (*) lift has nonzero",n,"component")
    K_E = generated([cyc((1,2,3,4,5)),cyc((1,6),(2,5))])
    K_S = generated([cyc((0,1,2)),cyc((3,4,5)),cyc((0,3,1,4),(2,5))])
    for K,n,label,expected_seen,expected_count in [(K_E,"21","A5",42,42),(K_S,"35","3^2:4",58,70)]:
        fd = fixdims(K)
        h1 = 3*fd["20"]+2*fd["15"]+2*fd["21"]+4*fd["35"]
        assert h1 == (2 if n == "21" else 4)
        conjugates = {frozenset(conj(g,k) for k in K) for g in A7}
        detected = 0
        v = projections[n]
        for subgroup in conjugates:
            avg = np.zeros(105,dtype=np.int64)
            for k in subgroup:
                avg[action[gidx[k]]] += v
            detected += bool(np.any(avg))
        assert (detected,len(conjugates)) == (expected_seen,expected_count)
        print(label,"order",len(K),"genus",h1//2,"detecting conjugates",detected,"of",len(conjugates))


if __name__ == "__main__":
    main()
