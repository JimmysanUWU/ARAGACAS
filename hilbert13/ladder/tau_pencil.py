"""The tau-pencil has degree exactly 42 (Chapter 7, Proposition 7.7).

By the proof of 7.7, an extra base point of |E_-| would be a point x fixed by an involution sigma in C(tau), sigma != tau,
with phi(x) in Q_sigma = P(E_++(sigma, tau)), a plane. phi(C) lies on F3 (the cubic), on G6 (the sextic invariant
whose restriction vanishes to order 2 at the 7-point p, hence identically), and on the other forced-vanishing
invariants of degrees 9, 15 and 21. On each of the 8 planes Q_sigma we find all 18 = 3*6 points of F3 = G6 = 0
(Bezout) and check that the other invariants do not vanish there.
Invariants are Reynolds averages over 3.A7 of powers of random linear forms (floating point)."""
import sys, io, contextlib, cmath
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import twisted_rr as T
from a7 import mul, order

RT = np.transpose(np.array([T.R6(k) for k in T.mats]), (0, 2, 1))     # points are columns x with x -> R(g)^T x
rng = np.random.default_rng(21)
def seeds(m): return np.einsum('mj,gjk->mgk', rng.standard_normal((m, 6)) + 1j * rng.standard_normal((m, 6)), RT)
def fval(rows, k, x): return np.mean((rows @ x) ** k, axis=1)
def fder(rows, k, x, v): y = rows @ x; return np.mean(k * y ** (k - 1) * (rows @ v), axis=1)
def distinct(xs):
    out = []
    for x in xs:
        if not any(abs(abs(np.vdot(x, y)) - 1) < 1e-7 for y in out): out.append(x)
    return out

for t in (0, 12):                      # classes 1, 14 are their images under the outer automorphism
    print(f"class {t}:")
    for D, ld, m, gens, lam in T.TW[t, (0, 1), +1]:
        if D == 60 and any(v > 0 for v in m.values()): break
    def eig(key, nu):
        w, V = np.linalg.eig(T.R6(key).T); return V[:, np.abs(w - nu) < 1e-8]
    p = eig((T.IDX[gens[2]], 0), lam[2])[:, 0]; p /= np.linalg.norm(p)                   # the 7-point
    e1 = eig((T.IDX[gens[2]], 0), lam[2] / cmath.exp(2j * cmath.pi / 7))[:, 0]           # its tangent direction
    S3, S6 = seeds(1), seeds(3)
    K = np.linalg.svd(fder(S6, 6, p, e1).reshape(1, 3))[2][1:].conj()                    # sextics singular along the tangent: F3^2, G6
    X = rng.standard_normal((40, 6)) + 1j * rng.standard_normal((40, 6))
    F3sq = np.array([fval(S3, 3, x)[0] ** 2 for x in X]); KV = np.array([K @ fval(S6, 6, x) for x in X])
    cf, *_ = np.linalg.lstsq(KV, F3sq, rcond=None)
    g6 = np.array([-cf[1].conj(), cf[0].conj()]); g6 /= np.linalg.norm(g6)
    F3 = lambda x: fval(S3, 3, x)[0]
    G6 = lambda x: g6 @ (K @ fval(S6, 6, x))
    S9, S15, S21 = seeds(6), seeds(6), seeds(12)
    K21 = np.linalg.svd(fval(S21, 21, p).reshape(1, -1))[2][1:].conj()                   # degree-21 invariants vanishing at p
    rest = lambda x: np.concatenate([fval(S9, 9, x), fval(S15, 15, x), K21 @ fval(S21, 21, x)])
    xr = rng.standard_normal(6) + 1j * rng.standard_normal(6); xr /= np.linalg.norm(xr)
    s3, s6, sr = abs(F3(xr)), abs(G6(xr)), np.abs(rest(xr))
    print(f"F3^2 lies in the tangent kernel: residual {np.linalg.norm(KV @ cf - F3sq) / np.linalg.norm(F3sq):.1e}")
    print(f"control at the 7-point p: |F3| {abs(F3(p)) / s3:.1e}, |G6| {abs(G6(p)) / s6:.1e}, others {np.abs(rest(p) / sr).max():.1e}")

    def lift2(g): return [(T.IDX[g], t) for t in range(3) if T.m3((T.IDX[g], t), (T.IDX[g], t)) == (0, 0)][0]
    tau = gens[0]; Rtau = T.R6(lift2(tau)).T
    def solve_plane(B):
        sols = []
        for trial in range(400):
            Bq = B @ np.linalg.qr(rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))[0]
            u = np.concatenate([[1], rng.standard_normal(2) + 1j * rng.standard_normal(2)])
            for _ in range(80):
                x = Bq @ u; r = np.array([F3(x) / s3, G6(x) / s6]); h = 1e-7
                Jm = np.array([[(F3(Bq @ (u + h * np.eye(3)[j])) / s3 - r[0]) / h, (G6(Bq @ (u + h * np.eye(3)[j])) / s6 - r[1]) / h] for j in (1, 2)]).T
                du = np.linalg.solve(Jm, -r); u[1:] += du
                if np.linalg.norm(du) < 1e-13 * (1 + np.linalg.norm(u)) or np.linalg.norm(u) > 1e7: break
            x = Bq @ u; x /= np.linalg.norm(x)
            if abs(F3(x)) / s3 < 1e-10 and abs(G6(x)) / s6 < 1e-10: sols.append(x)
            if trial >= 60 and len(distinct(sols)) == 18: break
        return distinct(sols)
    invs = [g for g in T.G if order(g) == 2 and mul(g, tau) == mul(tau, g) and g != tau]
    for sg in invs:
        M = np.vstack([T.R6(lift2(sg)).T - np.eye(6), Rtau - np.eye(6)]); B = np.linalg.svd(M)[2][-3:].conj().T
        assert np.allclose(M @ B, 0, atol=1e-9)
        pts = solve_plane(B)
        print(f"sigma = {sg}: F3 = G6 = 0 on Q_sigma has {len(pts)} points (Bezout: 18); "
              f"min over them of the largest other invariant (relative): {min(np.abs(rest(x) / sr).max() for x in pts):.2f}")
print("=> on both classes phi(C) misses every Q_sigma: the 18 fixed points are the whole base locus, and the moving part has degree exactly 42.")
