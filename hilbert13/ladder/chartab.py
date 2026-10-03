"""A7 character table from the exact ATLAS table in cover.py.

EXACT_TABLE contains exact values.  TABLE is its numerical evaluation for
numerical consumers; no eigensolver or rounding is used to identify characters.
The labels of the two conjugate 10s are interchangeable.
"""
import numpy as np
import sympy as sp
from a7 import *
from cover import TABLE as ATLAS_TABLE, CLASSES as ATLAS_CLASSES, check_table

def classes(G):
    left, cls = set(G), []
    while left:
        x = next(iter(left))
        c = {conj(g, x) for g in G}
        cls.append(sorted(c))
        left -= c
    cls.sort(key=lambda c: (order(c[0]), len(c)))
    return cls

CL = classes(A7)
idx = {}
for i, c in enumerate(CL):
    for x in c:
        idx[x] = i
k = len(CL)
sizes = np.array([len(c) for c in CL])
inv_class = [idx[inv(c[0])] for c in CL]

def class_matrix(j):
    """M[i][l] = # (x in C_j, y in C_l... ) : standard a_{j i l} with C_j C_i = sum a_{jil} C_l."""
    M = np.zeros((k, k))
    for i in range(k):
        for x in CL[j]:
            for y in CL[i]:
                M[i, idx[mul(x, y)]] += 1
        # divide: counts over all z in C_l of products; normalize per z
    for l in range(k):
        M[:, l] /= sizes[l]
    return M

check_table()
columns = []
seven_columns = iter((7,8))
for conjugacy_class in CL:
    ct = cycle_type(conjugacy_class[0])
    columns.append(next(seven_columns) if ct == (7,) else ATLAS_CLASSES.index(ct))
EXACT_TABLE = [[row[j] for j in columns] for row in ATLAS_TABLE.values()]


def char_table():
    return [np.array([complex(sp.N(value,17)) for value in row]) for row in EXACT_TABLE]


TABLE = char_table()

def pretty(z):
    if abs(z.imag) < 1e-8:
        return f"{z.real:.3g}"
    return f"{z.real:.3g}{z.imag:+.3g}i"

if __name__ == "__main__":
    print("classes:", [(cycle_type(c[0]), len(c)) for c in CL])
    for chi in TABLE:
        print([pretty(z) for z in chi])
    # orthogonality check
    G = np.array(TABLE)
    gram = (G * sizes) @ G.conj().T / 2520
    print("orthonormal:", np.allclose(gram, np.eye(k)))
