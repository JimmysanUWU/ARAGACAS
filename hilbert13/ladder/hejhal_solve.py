"""Vector-valued Hejhal method for one isotype of a (2,4,7) A7-curve (floating point; produces coef_cls*_M*.npz).

Model.  Poincare disk, order-7 point at 0, order-2 point P2 > 0 on the real axis, order-4 point P4 at angle pi/7.
X = half-turn about P2, Y = quarter-turn about P4, C = rotation by 2pi/7 about 0; X o Y o C = id.  For the triple
(a, b, c) of triples_data (a b c = 1) the curve is H / ker(phi), phi(X) = a, phi(Y) = b, phi(C) = c.
The 14_(5,2)-isotype consists of F : H -> V2 (V2 = 14_(5,2) inside the 2-subset module R^21, P_g = permutation
matrices) with F(delta z) = P_phi(delta) F(z).

Around 0, G = V^* Q14^T F splits by the eigenvalues zeta^j of rho(c); a component with eigenvalue zeta^j is
sum_{m = j mod 7} a_m R_|m|(|z|) e^{i m th}, with R_m the regular radial eigenfunction (hh_eval.py).  Imposing
F(Xz) = P_a F(z) at sample points on both sides of the side through P2 gives a linear system whose smallest singular
value vanishes at an eigenvalue; lambda is located by minimising it (golden section), then the null vector is taken.

Usage: python3 hejhal_solve.py CLASS M [depth]      (default depth 0.65; M = 90 used for the certificate)
"""
import os, sys, itertools
import numpy as np
from scipy.optimize import minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from triples_data import triples

PAIRS = list(itertools.combinations(range(7), 2)); PIDX = {p: i for i, p in enumerate(PAIRS)}
d72 = np.arccosh(np.cos(np.pi / 4) / np.sin(np.pi / 7)); d74 = np.arccosh(np.cos(np.pi / 7) / np.sin(np.pi / 7))
P2 = np.tanh(d72 / 2) + 0j


def Tm(p):
    return np.array([[1, p], [np.conj(p), 1]]) / np.sqrt(1 - abs(p) ** 2)


def act(M, z):
    return (M[0, 0] * z + M[0, 1]) / (M[1, 0] * z + M[1, 1])


MX = Tm(P2) @ np.diag([1j, -1j]) @ np.linalg.inv(Tm(P2))


def perm21(p):
    P = np.zeros((21, 21))
    for (i, j) in PAIRS:
        P[PIDX[tuple(sorted((p[i], p[j])))], PIDX[(i, j)]] = 1
    return P


inc = np.zeros((21, 7))
for (i, j), r in PIDX.items():
    inc[r, i] = inc[r, j] = 1
_q, _ = np.linalg.qr(inc)
_w, _U = np.linalg.eigh(np.eye(21) - _q @ _q.T)
Q14 = _U[:, _w > 0.5]


def radial(m, u, lam, tol=1e-17):
    t = np.sqrt(complex(lam - 0.25)); a = 0.5 - 1j * t; bb = m + 0.5 - 1j * t
    x = u ** 2; s = np.ones_like(u, dtype=complex); c = 1.0 + 0j; xp = np.ones_like(u)
    for n in range(5000):
        c = c * (a + n) * (bb + n) / ((m + 1 + n) * (n + 1)); xp = xp * x; s = s + c * xp
        if n > 5 and np.max(np.abs(c * xp)) < tol * np.max(np.abs(s)):
            break
    return np.real(u ** m * (1 - x) ** a * s)


class Solver:
    def __init__(self, cls, M, depth=0.65, seed=5):
        a, b, c = triples[cls]
        self.Pa, self.Pc = perm21(a), perm21(c)
        rc = Q14.T @ self.Pc @ Q14; ra = Q14.T @ self.Pa @ Q14
        ev, V = np.linalg.eig(rc)
        js = np.round(np.angle(ev) / (2 * np.pi / 7)).astype(int) % 7
        Vu = np.zeros((14, 14), complex)
        for j in set(js):
            idx = np.where(js == j)[0]; q, _ = np.linalg.qr(V[:, idx]); Vu[:, idx] = q
        self.V, self.W, self.M = Vu, Vu.conj().T @ ra @ Vu, M
        self.cols = [(k, m) for k in range(14) for m in range(-M, M + 1) if (m - js[k]) % 7 == 0]
        self.ms = sorted(set(abs(m) for _, m in self.cols))
        rng = np.random.default_rng(seed)
        n = int(3 * len(self.cols) / 14) + 80
        zs = []
        for _ in range(n):
            s = rng.uniform(-0.80, 0.80); d = rng.uniform(-depth, depth)
            Mm = Tm(P2) @ np.array([[np.cosh(s / 2), 1j * np.sinh(s / 2)], [-1j * np.sinh(s / 2), np.cosh(s / 2)]])
            zs.append(act(Mm, -np.tanh(d / 2) + 0j))
        self.zs = np.array(zs)

    def basis(self, z, lam):
        u, th = np.abs(z), np.angle(z)
        RR = {m: radial(m, u, lam) for m in self.ms}
        return np.array([RR[abs(m)] * np.exp(1j * m * th) for (k, m) in self.cols]).T

    def system(self, lam):
        Bz = self.basis(self.zs, lam); Bx = self.basis(act(MX, self.zs), lam)
        A = np.zeros((14 * len(self.zs), len(self.cols)), complex)
        for ci, (k, m) in enumerate(self.cols):
            blk = -np.outer(Bz[:, ci], self.W[:, k]); blk[:, k] += Bx[:, ci]; A[:, ci] = blk.ravel()
        return A

    def smin(self, lam):
        A = self.system(lam); A = A / np.linalg.norm(A, axis=0)
        return np.linalg.svd(A, compute_uv=False)[-1]

    def solve(self, bracket):
        r = minimize_scalar(self.smin, bracket=bracket, tol=1e-13)
        lam = r.x
        A = self.system(lam); nrm = np.linalg.norm(A, axis=0)
        _, sv, Vh = np.linalg.svd(A / nrm); coef = Vh[-1].conj() / nrm
        # the real solution space is one-dimensional: rotate the complex null vector to be real
        zt = np.array([0.3 + 0.1j, -0.2 + 0.25j, 0.1 - 0.4j])
        G = np.zeros((14, len(zt)), complex); B = self.basis(zt, lam)
        for ci, (k, m) in enumerate(self.cols):
            G[k] += coef[ci] * B[:, ci]
        Fz = Q14 @ (self.V @ G)
        coef = coef * np.exp(-0.5j * np.angle(np.sum(Fz ** 2)))
        Am = np.zeros((self.M + 1, 21)); Bm = np.zeros((self.M + 1, 21))
        for ci, (k, m) in enumerate(self.cols):
            v = Q14 @ (self.V[:, k] * coef[ci])
            Am[abs(m)] += np.real(v)
            if m != 0:
                Bm[abs(m)] += -np.sign(m) * np.imag(v)
        # real part (Re and Im are both solutions); normalise the larger one
        return lam, sv[-1], Am, Bm


if __name__ == "__main__":
    cls, M = int(sys.argv[1]), int(sys.argv[2])
    depth = float(sys.argv[3]) if len(sys.argv) > 3 else 0.65
    br = (0.344, 0.3463, 0.349) if cls in (0, 1) else (0.357, 0.3597, 0.362)
    s = Solver(cls, M, depth)
    lam, sm, Am, Bm = s.solve(br)
    print(f"class {cls}  M = {M}  lambda = {lam!r}  smallest singular value {sm:.3e}")
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"coef_cls{cls}_M{M}.npz")
    np.savez(out, lam=lam, A=Am, B=Bm)
    print("saved", out)
