"""Equations of the degree-60 model phi(C) in P^5 (Chapter 7, §7.6, class 0) [N].

1. The image of a 4-point: on the line P(E_-1(h)) (h of order 4), the sextic G6 has exactly two double zeros, and they are
   singular points of {G6 = 0}; they are phi(Fix h).
2. Kernels: the invariants of degree 12, 18, 24 vanishing at phi(q), and of degree 21 vanishing at the 7-point p, vanish on phi(C)
   (their restrictions are multiples of s_{2D7}, s_{3D7}, s_{4D7}, s_{2D4}). With F3, G6, the nonics and degree 15 they cut out
   a smooth curve through p (Jacobian rank 4), traced by continuation and spread by the group.
3. Hilbert function of phi(C) in degrees 2..5, the cubic and quartic equations, and a test that they cut out phi(C):
   a random hyperplane meets their zero set in 60 points, and a hyperplane through P(E_+(tau)) in 18 + 42.
4. Klein quadrics (Proposition 7.9): for H = L2(7) the H-invariant quadrics form a line; Q_H vanishes on the 24 seven-points of H,
   to order 5 at p, and over a conjugacy class prod Q_H / S6^5 is constant on phi(C).
"""
import sys, io, contextlib, cmath, itertools
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import twisted_rr as T

RT = np.transpose(np.array([T.R6(k) for k in T.mats]), (0, 2, 1))
rng = np.random.default_rng(21)
def seeds(m): return np.einsum('mj,gjk->mgk', rng.standard_normal((m, 6)) + 1j * rng.standard_normal((m, 6)), RT)
def fval(rows, k, x): return np.mean((rows @ x) ** k, axis=1)
def fder(rows, k, x, v): y = rows @ x; return np.mean(k * y ** (k - 1) * (rows @ v), axis=1)
def kern_at(S, k, x): return np.linalg.svd(fval(S, k, x).reshape(1, -1))[2][1:].conj()
def distinct(xs):
    out = []
    for x in xs:
        if not any(abs(abs(np.vdot(x, y)) - 1) < 1e-7 for y in out): out.append(x)
    return out
for D, ld, m, gens, lam in T.TW[0, (0, 1), +1]:
    if D == 60 and any(v > 0 for v in m.values()): break
def eig(key, nu):
    w, V = np.linalg.eig(T.R6(key).T); return V[:, np.abs(w - nu) < 1e-8]
a7 = cmath.exp(2j * cmath.pi / 7)
p = eig((T.IDX[gens[2]], 0), lam[2])[:, 0]; p /= np.linalg.norm(p)
e1 = eig((T.IDX[gens[2]], 0), lam[2] / a7)[:, 0]
S3, S6 = seeds(1), seeds(3)
K6 = np.linalg.svd(fder(S6, 6, p, e1).reshape(1, 3))[2][1:].conj()
X = rng.standard_normal((40, 6)) + 1j * rng.standard_normal((40, 6))
cf, *_ = np.linalg.lstsq(np.array([K6 @ fval(S6, 6, x) for x in X]), np.array([fval(S3, 3, x)[0] ** 2 for x in X]), rcond=None)
g6 = np.array([-cf[1].conj(), cf[0].conj()]); g6 /= np.linalg.norm(g6)
G6 = lambda x: g6 @ (K6 @ fval(S6, 6, x))

# 1. the 4-point
Q4, _ = np.linalg.qr(eig((T.IDX[gens[1]], 0), lam[1]))
zs = np.exp(2j * np.pi * np.arange(12) / 12)
c6 = np.linalg.lstsq(np.array([[z ** k for k in range(7)] for z in zs]), np.array([G6(Q4 @ np.array([1, z])) for z in zs]), rcond=None)[0][::-1]
dbl = [r for r in np.roots(np.polyder(c6)) if abs(np.polyval(c6, r)) < 1e-8 * np.abs(c6).max()]
s6 = abs(G6(rng.standard_normal(6) / 2.45 + 0j))
def grad(f, x, h=1e-6): return np.array([(f(x + h * np.eye(6)[i]) - f(x - h * np.eye(6)[i])) / (2 * h) for i in range(6)])
x0 = Q4 @ np.array([1, dbl[0]]); x0 /= np.linalg.norm(x0)
print(f"1. G6 on the 4-point line: {len(dbl)} double zeros; |grad G6| there (relative): "
      f"{[f'{np.linalg.norm(grad(G6, Q4 @ np.array([1, r]) / np.linalg.norm(Q4 @ np.array([1, r])))) / s6:.0e}' for r in dbl]}")

# 2. the full system and a branch through p
S9, S15, S21, S12, S18, S24 = seeds(6), seeds(6), seeds(60), seeds(16), seeds(40), seeds(100)
blocks = [(S3, 3, None), (S6, 6, (g6 @ K6).reshape(1, -1)), (S9, 9, None), (S15, 15, None), (S21, 21, kern_at(S21, 21, p)),
          (S12, 12, kern_at(S12, 12, x0)), (S18, 18, kern_at(S18, 18, x0)), (S24, 24, kern_at(S24, 24, x0))]
def f_(x):
    return np.concatenate([(v := np.mean((r @ x) ** k, axis=1)) if K is None else K @ np.mean((r @ x) ** k, axis=1) for r, k, K in blocks])
def J_(x):
    out = []
    for r, k, K in blocks:
        y = r @ x; d = np.einsum('mg,mgj->mj', k * y ** (k - 1), r) / r.shape[1]; out.append(d if K is None else K @ d)
    return np.vstack(out)
xr = rng.standard_normal(6) + 1j * rng.standard_normal(6); xr /= np.linalg.norm(xr)
sc = np.abs(f_(xr)) + 1e-300
ell, ell0 = e1.conj() / np.vdot(e1, e1), p.conj()
def newton_slice(x, c):
    for _ in range(80):
        r = np.concatenate([f_(x) / sc, [ell @ x - c, ell0 @ x - 1]]); Jx = np.vstack([J_(x) / sc[:, None], ell, ell0])
        dx = np.linalg.lstsq(Jx, -r, rcond=None)[0]; x = x + dx
        if np.linalg.norm(dx) < 1e-14 * np.linalg.norm(x): break
    return x, np.linalg.norm(f_(x / np.linalg.norm(x)) / sc)
arc = []; worst = 0
for ts in (np.linspace(0.02, 0.6, 60), 1j * np.linspace(0.02, 0.6, 30)):
    prev = [p / (p.conj() @ p)]
    for t in ts:
        y, res = newton_slice(prev[-1] if len(prev) < 2 else 2 * prev[-1] - prev[-2], t)
        if res > 1e-8: break
        worst = max(worst, res); prev.append(y); arc.append(y / np.linalg.norm(y))
ranks = {int(np.sum((s := np.linalg.svd(J_(y) / sc[:, None], compute_uv=False)) / s[0] > 1e-7)) for y in arc}
print(f"2. the system has a branch through p along its tangent: {len(arc)} points by continuation (residual <= {worst:.0e}), Jacobian ranks {ranks} (4 = curve)")
orb = np.einsum('gij,nj->gni', np.array([T.R6(k).T for k in T.mats]), np.array(arc)).reshape(-1, 6)
orb /= np.linalg.norm(orb, axis=1)[:, None]
sub = orb[rng.choice(len(orb), 3000, replace=False)]

# 3. Hilbert function and equations
def monos(k): return list(itertools.combinations_with_replacement(range(6), k))
def mv(Xs, ms): return np.stack([np.prod(Xs[:, list(mm)], axis=1) for mm in ms], 1)
null = {}
for k, h0 in ((2, "h0(2L) >= 21"), (3, "45 + h0(B+T)"), (4, "105"), (5, "165")):
    A = mv(sub, monos(k)); nrm = np.linalg.norm(A, axis=0); s = np.linalg.svd(A / nrm, compute_uv=False)
    r = int(np.sum(s / s[0] > 1e-9)); _, _, vh = np.linalg.svd(A / nrm); null[k] = vh[r:].conj() / nrm[None, :]
    print(f"3. degree {k}: {len(nrm)} forms, rank {r} on phi(C) [h^0(kL) = {h0}]; {len(nrm) - r} vanish; gap {s[r - 1] / s[0]:.1e} -> {s[r] / s[0] if r < len(s) else 0:.1e}")
M3, M4 = monos(3), monos(4)
def eqs(x): return np.concatenate([mv(x[None], M3)[0] @ null[3].T, mv(x[None], M4)[0] @ null[4].T])
s34 = np.abs(eqs(xr)) + 1e-300
def jac(x, h=1e-7): f0 = eqs(x); return np.array([(eqs(x + h * np.eye(6)[i]) - f0) / h for i in range(6)]).T
def section(B, ntrial):
    mdim = B.shape[1]; sols = []
    for _ in range(ntrial):
        Bq = B @ np.linalg.qr(rng.standard_normal((mdim, mdim)) + 1j * rng.standard_normal((mdim, mdim)))[0]
        u = np.concatenate([[1], rng.standard_normal(mdim - 1) + 1j * rng.standard_normal(mdim - 1)])
        for _ in range(60):
            x = Bq @ u; du = np.linalg.lstsq((jac(x) / s34[:, None] @ Bq)[:, 1:], -eqs(x) / s34, rcond=None)[0]; u[1:] += du
            if np.linalg.norm(du) < 1e-13 * (1 + np.linalg.norm(u)) or np.linalg.norm(u) > 1e7: break
        x = Bq @ u; x /= np.linalg.norm(x)
        if np.linalg.norm(eqs(x) / s34) < 1e-9: sols.append(x)
    return distinct(sols)
pts = section(np.linalg.qr(rng.standard_normal((6, 5)) + 1j * rng.standard_normal((6, 5)))[0], 1500)
print(f"   the cubic and quartics on a random hyperplane: {len(pts)} common zeros (deg phi(C) = 60); on phi(C): "
      f"{sum(np.linalg.norm(f_(x) / sc) < 1e-8 for x in pts)}")
tl = [(T.IDX[gens[0]], t) for t in range(3) if T.m3((T.IDX[gens[0]], t), (T.IDX[gens[0]], t)) == (0, 0)][0]
w, V = np.linalg.eig(T.R6(tl).T); Ep, _ = np.linalg.qr(V[:, np.abs(w - 1) < 1e-8])
ann = np.linalg.svd(Ep.T)[2][4:].conj(); l = rng.standard_normal(2) @ ann
pts = section(np.linalg.svd(l.reshape(1, -1))[2][1:].conj().T, 1500)
fixed = sum(np.linalg.norm(ann @ x) < 1e-7 for x in pts)
print(f"   a hyperplane of the tau-pencil: {len(pts)} points = {fixed} in P(E_+) + {len(pts) - fixed} moving (Proposition 7.7)")

# 4. Klein quadrics
from a7 import A7, mul, inv, order, generated
c7 = gens[2]
Ls = []
for s_ in A7:
    if order(s_) == 2 and order(mul(s_, c7)) == 3:
        Hs = generated([s_, c7])
        if len(Hs) == 168 and not any(set(Hs) == set(K) for K in Ls): Ls.append(Hs)
def lift_same_order(g):
    o = order(g)
    for t in range(3):
        q_ = (T.IDX[g], t); z = q_
        for _ in range(o - 1): z = T.m3(z, q_)
        if z == (0, 0): return q_
M2 = monos(2); Xq = rng.standard_normal((60, 6)) + 1j * rng.standard_normal((60, 6))
def quadrics_of(Hs):
    rows = [mv(Xq @ T.R6(lift_same_order(g)), M2) - mv(Xq, M2) for g in [g for g in Hs if order(g) in (2, 7)][:6]]
    _, sv, vh = np.linalg.svd(np.vstack(rows)); return int(np.sum(sv < 1e-9 * sv[0])), vh[-1].conj()
crow = np.einsum('j,gjk->gk', rng.standard_normal(6) + 1j * rng.standard_normal(6), RT)
S6f = lambda Xs: np.mean((Xs @ crow.T) ** 6, axis=1)
for Hs in Ls:
    k, qH = quadrics_of(Hs)
    OH = distinct([T.R6(lift_same_order(g) if order(g) in (2, 4, 7) else (T.IDX[g], 0)).T @ p for g in Hs])
    OH = [x / np.linalg.norm(x) for x in OH]
    typ = np.median(np.abs(mv(sub, M2) @ qH))
    tv = np.array([abs(ell @ y / (ell0 @ y)) for y in arc[:15]]); qv = np.array([abs(mv((y / np.linalg.norm(y))[None], M2)[0] @ qH) for y in arc[:15]])
    slope = np.polyfit(np.log(tv), np.log(qv), 1)[0]
    cls_ = []
    for g in A7:
        K_ = frozenset(mul(mul(g, h), inv(g)) for h in Hs)
        if K_ not in cls_: cls_.append(K_)
        if len(cls_) == 15: break
    prod = np.prod(np.stack([mv(sub[:400], M2) @ quadrics_of(list(K_))[1] for K_ in cls_]), axis=0) / S6f(sub[:400]) ** 5
    print(f"4. L2(7) containing <c>: {k} invariant quadric; |O_H| = {len(OH)}, max |Q_H| there {max(abs(mv(x[None], M2)[0] @ qH) for x in OH) / typ:.0e}; "
          f"vanishing order at p {slope:.2f}; prod over the class / S6^5 on phi(C): spread {np.std(prod) / abs(np.mean(prod)):.0e}")
