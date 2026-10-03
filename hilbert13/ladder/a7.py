"""Brute-force toolkit for A7 acting on a (2,4,7) curve C.

Permutations of {0..6} are tuples p with p[i] = image of i.
Composition convention: mul(p, q) = p o q (apply q first).
"""
from itertools import permutations
from math import gcd

N = 7
E = tuple(range(N))


def mul(p, q):
    return tuple(p[q[i]] for i in range(N))


def inv(p):
    r = [0] * N
    for i, x in enumerate(p):
        r[x] = i
    return tuple(r)


def sign(p):
    s, seen = 1, [False] * N
    for i in range(N):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def cycle_type(p):
    seen, ct = [False] * N, []
    for i in range(N):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            ct.append(L)
    return tuple(sorted(ct, reverse=True))


def order(p):
    o = 1
    for L in cycle_type(p):
        o = o * L // gcd(o, L)
    return o


def power(p, k):
    r = E
    for _ in range(k % order(p)):
        r = mul(r, p)
    return r


def cyc(*cycles):
    """cyc((0,1),(2,3)) -> permutation."""
    p = list(range(N))
    for c in cycles:
        for i in range(len(c)):
            p[c[i]] = c[(i + 1) % len(c)]
    return tuple(p)


A7 = [p for p in permutations(range(N)) if sign(p) == 1]
assert len(A7) == 2520


def conj(g, h):
    """g h g^-1"""
    return mul(mul(g, h), inv(g))


def centralizer(h, G=A7):
    return [g for g in G if mul(g, h) == mul(h, g)]


# A7-conjugacy classes (7-cycles split into two classes).
def a7_class(h):
    return frozenset(conj(g, h) for g in A7)


def generated(gens):
    S = {E}
    frontier = [E]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = mul(g, x)
                if y not in S:
                    S.add(y)
                    new.append(y)
        frontier = new
    return S
