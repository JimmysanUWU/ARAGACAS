"""Numerical Laplace spectrum of the (2,4,7) A7-curve C by P1 finite elements.

C is tiled by 2*2520 hyperbolic triangles with angles pi/2 (at points over a), pi/4 (over b),
pi/7 (over c): the preimages U_g, L_g (g in A7) of the upper/lower half-planes under C -> C/A7.
Gluing: U_g ~ L_g along (0,1); U_g ~ L_{g b} along (1,inf); U_g ~ L_{g a^-1} along (inf,0).
Each triangle is subdivided in the Klein model (straight lines = geodesics); every small triangle
is replaced by the flat triangle with the same hyperbolic side lengths (O(h^2) error).
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from math import pi, cos, sin, acosh, tanh
from a7 import *
from triples_data import triples

a, b, c = triples[int(sys.argv[2]) if len(sys.argv) > 2 else 0]
n = int(sys.argv[1]) if len(sys.argv) > 1 else 6

# reference triangle: A (angle pi/2, over a), B (pi/4, over b), Cv (pi/7, over c)
al, be, ga = pi / 2, pi / 4, pi / 7
cs = acosh((cos(ga) + cos(al) * cos(be)) / (sin(al) * sin(be)))   # AB
bs = acosh((cos(be) + cos(al) * cos(ga)) / (sin(al) * sin(ga)))   # AC
KA, KB, KC = np.array([0.0, 0.0]), np.array([tanh(cs), 0.0]), np.array([0.0, tanh(bs)])

def kpt(i, j):
    return KA + (i / n) * (KB - KA) + (j / n) * (KC - KA)

def hdist(p, q):
    num = 1 - p @ q
    den = np.sqrt((1 - p @ p) * (1 - q @ q))
    return acosh(max(1.0, num / den))

# local element matrices for each small triangle of the reference grid
loc_nodes = [(i, j) for i in range(n + 1) for j in range(n + 1 - i)]
small = []
for i in range(n):
    for j in range(n - i):
        small.append(((i, j), (i + 1, j), (i, j + 1)))
        if i + j + 2 <= n:
            small.append(((i + 1, j), (i + 1, j + 1), (i, j + 1)))
elems = []
for tri in small:
    P = [kpt(*v) for v in tri]
    l = [hdist(P[1], P[2]), hdist(P[0], P[2]), hdist(P[0], P[1])]  # opposite vertices 0,1,2
    s = sum(l) / 2
    A = np.sqrt(max(s * (s - l[0]) * (s - l[1]) * (s - l[2]), 0))
    Kl = np.zeros((3, 3))
    for k in range(3):
        i1, i2 = (k + 1) % 3, (k + 2) % 3
        cot = (l[i1] ** 2 + l[i2] ** 2 - l[k] ** 2) / (4 * A)
        # edge between i1,i2 has weight cot(angle at k)/2
        Kl[i1, i2] -= cot / 2; Kl[i2, i1] -= cot / 2
        Kl[i1, i1] += cot / 2; Kl[i2, i2] += cot / 2
    Ml = A / 12 * (np.ones((3, 3)) + np.eye(3))
    elems.append((tri, Kl, Ml))

# global labelling
Gidx = {g: k for k, g in enumerate(A7)}
def cos_label(g, x):
    return min(Gidx[mul(g, power(x, k))] for k in range(order(x)))
ba = mul(b, a)
binv, ainv = inv(b), inv(a)
labels = {}
def node_id(key):
    if key not in labels:
        labels[key] = len(labels)
    return labels[key]

def global_key(side, g, i, j):
    """side 'U' or 'L'; sheet g; local node (i,j)."""
    if side == "U":
        e01, e1i, ei0 = g, g, g                     # edge labels = the U-sheet they belong to
        vinf = cos_label(g, ba)
    else:
        e01, e1i, ei0 = g, mul(g, binv), mul(g, a)
        vinf = cos_label(mul(g, binv), ba)
    if (i, j) == (0, 0):
        return ("v0", cos_label(g, a))
    if (i, j) == (n, 0):
        return ("v1", cos_label(g, b))
    if (i, j) == (0, n):
        return ("vi", vinf)
    if j == 0:
        return ("e01", Gidx[e01], i)
    if i == 0:
        return ("ei0", Gidx[ei0], j)
    if i + j == n:
        return ("e1i", Gidx[e1i], j)
    return ("int", side, Gidx[g], i, j)

rows, cols, kv, mv = [], [], [], []
for side in ("U", "L"):
    for g in A7:
        ids = {v: node_id(global_key(side, g, *v)) for v in loc_nodes}
        for tri, Kl, Ml in elems:
            I = [ids[v] for v in tri]
            for r in range(3):
                for s_ in range(3):
                    rows.append(I[r]); cols.append(I[s_]); kv.append(Kl[r, s_]); mv.append(Ml[r, s_])
N = len(labels)
K = sp.csr_matrix((kv, (rows, cols)), shape=(N, N))
M = sp.csr_matrix((mv, (rows, cols)), shape=(N, N))
area = M.sum()
print(f"n={n} nodes={N} area={area:.4f} (exact 540*pi={540*pi:.4f})", flush=True)
vals = sla.eigsh(K, k=40, M=M, sigma=-0.05, which="LM", return_eigenvectors=False)
vals = np.sort(vals)
print("eigenvalues:", np.round(vals, 5).tolist())
# group into clusters to read off multiplicities
cl, cur = [], [vals[0]]
for v in vals[1:]:
    if abs(v - cur[-1]) < 1e-4 * max(1, abs(v)): cur.append(v)
    else: cl.append(cur); cur = [v]
cl.append(cur)
print("clusters (value, multiplicity):", [(round(float(np.mean(x)), 5), len(x)) for x in cl])
