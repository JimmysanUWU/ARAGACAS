"""Character table of A7 by the Burnside/Dixon class-algebra method (numerical, then rounded)."""
import numpy as np
from a7 import *

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

def char_table():
    rng = np.random.default_rng(0)
    Ms = [class_matrix(j) for j in range(k)]
    A = sum(rng.standard_normal() * M for M in Ms)
    w, V = np.linalg.eig(A)
    chars = []
    for t in range(k):
        v = V[:, t]
        v = v / v[0]  # v_i = size_i chi(C_i)/chi(1)
        # chi(1)^2 = |G| / sum_i v_i conj(v_i)/size_i
        s = sum(v[i] * np.conj(v[i]) / sizes[i] for i in range(k)).real
        d = np.sqrt(2520 / s)
        chi = np.array([d * v[i] / sizes[i] for i in range(k)])
        chars.append(chi)
    chars.sort(key=lambda c: (round(c[0].real), round(c[1].real, 3), round(c[-1].imag, 3)))
    return chars

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
