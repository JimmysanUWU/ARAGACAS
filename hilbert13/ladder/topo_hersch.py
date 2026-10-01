"""Topological Hersch (FRAMEWORK_TOPOLOGICAL_HERSCH.md, sections 3-4): numerical constants on the full (2,4,7)
curve, P1 elements with n^2 cells per triangle.
  - lambda_1, lambda' (next eigenvalue) and E_1 (14-dimensional);
  - the bracket Gram of {phi_i, phi_j} on wedge^2 E_1 and its dual-norm factors <beta,(Delta-l1)^-1 beta>/|beta|^2
    per isotypic block (Schur: constant on each block), giving b*_sigma;
  - the gradient constant gamma = sqrt(A max_p lambda_max S(p));
  - the worst-case right side of (TH) (unsharpened form) for m = 23..27.
Usage: python3 topo_hersch.py n triple_index [--save]   (memory: about 5 GB at n = 12).
--save writes th_n{n}_t{triple}.npz for topo_hersch_extras.py."""
import sys, time
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orbifold import reference_triangle, quotient_tiles, build_dofs, ref_elements, numeric_elem_mats, metric_R
from triples_data import triples
n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
tri = int(sys.argv[2]) if len(sys.argv) > 2 else 0
t0 = time.time()
XA, XB, XC = reference_triangle(); els = ref_elements(n); mats = numeric_elem_mats(n, XA, XB, XC, "P1")
NE = len(els)
a, b, _ = triples[tri]
greps, glue, _ = quotient_tiles(a, b, [], None)
uf, _ = build_dofs(n, len(greps), glue, "P1")
NT = 2 * len(greps)
gid = {}
ED = np.zeros((NT, NE, 3), dtype=np.int64)
for tile in range(NT):
    for e, trg in enumerate(els):
        for p, v in enumerate(trg):
            r, s = uf.find((tile, (v,)))
            k = gid.get(r)
            if k is None:
                k = gid[r] = len(gid)
            ED[tile, e, p] = k
N = len(gid)
Kel = np.array([Ke for Ke, Me in mats]); Mel = np.array([Me for Ke, Me in mats])
rows = np.repeat(ED[..., :, None], 3, -1).reshape(-1); cols = np.repeat(ED[..., None, :], 3, -2).reshape(-1)
K = sp.csc_matrix((np.broadcast_to(Kel, (NT,) + Kel.shape).reshape(-1), (rows, cols)), shape=(N, N))
M = sp.csc_matrix((np.broadcast_to(Mel, (NT,) + Mel.shape).reshape(-1), (rows, cols)), shape=(N, N))
del rows, cols
print(f"n={n} triple={tri}: {N} dofs, assembled ({time.time()-t0:.0f}s)", flush=True)
ev, V = sla.eigsh(K, k=17, M=M, sigma=-1e-3, which="LM")
o = np.argsort(ev); ev, V = ev[o], V[:, o]
print("lowest eigenvalues:", np.round(ev, 5), flush=True)
lam1 = ev[1:15].mean(); lam_next = ev[15]
E1 = V[:, 1:15]
# element geometry
Gst, AR, wR, ast = [], [], [], []
for trg in els:
    P = [np.array(v, float) / n for v in trg]
    T = np.array([[P[1][0]-P[0][0], P[2][0]-P[0][0]], [P[1][1]-P[0][1], P[2][1]-P[0][1]]])
    Ti = np.linalg.inv(T); Gm = np.zeros((2, 3)); Gm[:, 1], Gm[:, 2] = Ti[0], Ti[1]; Gm[:, 0] = -Ti[0]-Ti[1]
    bc = sum(P) / 3; A_, w_ = metric_R(XA, XB, XC, np.array([bc[0]]), np.array([bc[1]]))
    Gst.append(Gm); AR.append(A_[0]); wR.append(w_[0]); ast.append(abs(np.linalg.det(T)) / 2)
Gst, AR, wR, ast = map(np.array, (Gst, AR, wR, ast))
areah = wR * ast
sgn = np.where(np.arange(NT) % 2 == 0, 1.0, -1.0)
dphi = np.einsum("ejv,tevl->tejl", Gst, E1[ED])                     # (NT, NE, 2, 14) st-gradients
# gradient constant: lambda_max of S = sum_l grad phi_l grad phi_l^T in the metric
Minv = AR / wR[:, None, None]; Lc = np.linalg.cholesky(Minv)
S = np.einsum("teal,tebl->teab", dphi, dphi)
Sm = np.einsum("eia,teab,ejb->teij", Lc, S, Lc)
lmax = np.linalg.eigvalsh(Sm)[..., -1]
trS = np.trace(Sm, axis1=-2, axis2=-1)
print(f"frame energy density: [{trS.min():.5f},{trS.max():.5f}] (14 l1/A = {14*lam1/(540*np.pi):.5f}); "
      f"2 lmax/tr in [{(2*lmax/trS).min():.3f},{(2*lmax/trS).max():.3f}]", flush=True)
del S, Sm
iu = np.triu_indices(14, 1)
def brackets(ts):                                                  # (len ts, NE, 91)
    d = dphi[ts]
    j = d[:, :, 0, :, None] * d[:, :, 1, None, :] - d[:, :, 1, :, None] * d[:, :, 0, None, :]
    return j[..., iu[0], iu[1]] / wR[None, :, None] * sgn[ts][:, None, None]
Bg = np.zeros((91, 91))
for c0 in range(0, NT, 400):
    ts = np.arange(c0, min(NT, c0 + 400)); br = brackets(ts)
    Bg += np.einsum("tep,teq,e->pq", br, br, areah)
bw, bV = np.linalg.eigh(Bg); o = np.argsort(bw)[::-1]; bw, bV = bw[o], bV[:, o]
blocks, start = [], 0
while start < 91:
    end = start
    while end < 91 and abs(bw[end] - bw[start]) < 2e-3 * bw[start]:
        end += 1
    blocks.append((start, end)); start = end
print("bracket Gram blocks:", [(round(float(bw[s]), 8), e - s) for s, e in blocks], flush=True)
shift = lam1 - 1e-4
lu = sla.splu((K - shift * M).tocsc())
Vd = V[:, :15]
res = []
for s_, e_ in blocks:
    rat = []
    for col in (s_, s_ + 1):
        beta = np.zeros((NT, NE))
        for c0 in range(0, NT, 400):
            ts = np.arange(c0, min(NT, c0 + 400)); beta[ts] = brackets(ts) @ bV[:, col]
        f = np.zeros(N); np.add.at(f, ED.reshape(-1, 3), np.repeat((beta * areah[None, :] / 3)[..., None], 3, -1).reshape(-1, 3))
        f -= M @ (Vd @ (Vd.T @ f))
        u = lu.solve(f); u -= Vd @ (Vd.T @ (M @ u))
        rat.append((f @ u) / (beta ** 2 * areah).sum())
    res.append((bw[s_], e_ - s_, float(np.mean(rat))))
    print(f"block dim {e_-s_:2d}: b = {bw[s_]:.4e}  ratio = {np.mean(rat):.4f}  -> mu_eff = {lam1 + 1/np.mean(rat):.4f}"
          f"   b* = {bw[s_]*np.mean(rat):.4e}", flush=True)
if "--save" in sys.argv:
    np.savez(f"th_n{n}_t{tri}.npz", ev=ev, bw=bw, bV=bV, blocks=np.array(blocks), res=np.array([(r[0], r[1], r[2]) for r in res]), lmax=lmax, trS=trS, areah=areah)
A = 540 * np.pi
bstar = max(b_ * r_ for b_, _, r_ in res)
Nst = A * np.sqrt(bstar / 3)
gam = np.sqrt(lmax.max() * A)
print(f"lambda1 = {lam1:.5f}, lambda_next = {lam_next:.5f};  sup||N~||_* <= {Nst:.2f},  sup||grad y||_inf <= {gam:.3f}")
for m in range(23, 28):
    eps = 8 * np.pi * m - lam1 * A
    if eps <= 0:
        print(f"m={m}: excluded by Li-Yau"); continue
    s = eps / (lam_next - lam1); Ez = eps + lam1 * s
    ta, tb, tc = 3 * np.sqrt(eps) * Nst, gam * np.sqrt(s * Ez), Ez / 2
    print(f"m={m}: 4pi m = {4*np.pi*m:6.1f}  worst-case RHS = {ta+tb+tc:6.1f} [{ta:.1f}+{tb:.1f}+{tc:.1f}]"
          f"  -> {'EXCLUDED' if ta+tb+tc < 4*np.pi*m else 'not excluded'}", flush=True)
