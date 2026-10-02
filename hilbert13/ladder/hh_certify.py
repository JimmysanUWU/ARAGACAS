"""Certificate: gon(C) >= 25 for the (2,4,7) A7-curves of classes 0 and 1, via harmonic Hersch on an exact trial space
built from near-exact eigenfunctions ("vector-valued Hejhal").  Write-up: 3_HARMONIC_HERSCH.md.

Inputs (certified elsewhere, certify_th_output.txt / certificate.txt):
    E_1 = one copy of 14_(5,2),  lambda_1 >= 0.340893,  every eigenvalue on (1 + E_1)^perp is >= b = 0.5599822.

Trial space.  Ft = hh_eval.TrialFunction: an exact lam*-eigenfunction of the hyperbolic disk, V2-valued and exactly
C-equivariant (C = rotation by 2pi/7 about the order-7 point 0).  Centres = Delta-orbit of 0; on each centre g(0) put
F_g(z) = P_phi(g) Ft(g^-1 z) (well defined by C-equivariance).  With h(d) a C^2 cutoff (= 1 for d <= r1, 0 for d >= r2)
and H = sum_c h(d(z,c)),
    Psi(z) = sum_c (h(d(z,c))/H(z)) F_c(z)
is EXACTLY Delta-equivariant and C^2, so its 21 components span a copy E_h of 14_(5,2) in C^2(C) (Schur).  Since
(Delta - lam*) F_c = 0,
    (Delta - lam*) Psi = sum_{c != 0} [Delta, chi_c] (F_c - F_0)       on the sector (where h(d(z,0)) = 1),
so the residual is controlled by the mismatches D_c = F_c - Ft, which are themselves exact lam*-eigenfunctions and are
bounded by circle sampling (Fourier-Legendre coefficients with rigorous aliasing and tail bounds).

Everything below is ball arithmetic (python-flint, 106 bits); floats are used only to choose grids.
Usage:  python3 hh_certify.py CLASS [coef_file]
"""
import os, sys, time, itertools
import numpy as np
from flint import arb, acb, arb_mat, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hh_eval import TrialFunction, perm21_index, apply_perm, PAIRS
from a7 import A7, mul, inv, order
from chartab import CL, TABLE

ctx.prec = 106
HERE = os.path.dirname(os.path.abspath(__file__))
PI = arb.pi()
LAMBDA1_LO = arb("0.340893")         # certificate.txt
B_GAP = arb("0.5599822")             # certify_th_output.txt (count on Q1, n = 96), with Q0 >= 0.997, Q2 >= 0.689
AREA = 540 * PI
TAU = arb(16) / 45                   # 8 pi * 24 / A


def amax(xs):
    m = arb(0)
    for x in xs:
        m = m.max(x)
    return m


def up(x):
    """rigorous upper bound of |x| as an exact arb."""
    return abs(x).upper()


# ------------------------------------------------------------------------------------------------ geometry
class Geometry:
    def __init__(self):
        s7, c7 = (PI / 7).sin(), (PI / 7).cos()
        self.d72 = ((PI / 4).cos() / s7).acosh()
        self.d74 = (c7 / s7).acosh()                       # circumradius of the heptagon star
        self.p2 = (self.d72 / 2).tanh()
        self.q4 = acb((self.d74 / 2).tanh()) * acb(c7, s7)
        self.u_circ = (self.d74 / 2).tanh()
        self.zeta = acb((2 * PI / 7).cos(), (2 * PI / 7).sin())
        self.x0 = (1 + self.p2 ** 2) / (2 * self.p2)

    @staticmethod
    def T(w, xi):
        return (xi + w) / (1 + w.conjugate() * xi)

    @staticmethod
    def Tinv(w, z):
        return (z - w) / (1 - w.conjugate() * z)

    def X(self, z):
        p = acb(self.p2)
        return self.T(p, -self.Tinv(p, z))

    def Y2(self, z):
        return self.T(self.q4, -self.Tinv(self.q4, z))

    def rot(self, z, k):
        return z * self.zeta ** k

    @staticmethod
    def dist(z, w):
        """Interval enclosure of the hyperbolic distance d(z, w) for z, w in balls (upper end +inf if unbounded)."""
        r = abs(z - w) / abs(1 - w.conjugate() * z)
        lo, hi = r.lower().max(arb(0)), r.upper()
        tri = 2 * abs(z).upper().atanh() + 2 * abs(w).upper().atanh()      # d(z,0) + d(0,w); needs |z|,|w| < 1
        assert tri.is_finite()
        dlo = 2 * lo.atanh() if lo < 1 else arb(0)
        dhi = (2 * hi.atanh()).min(tri) if hi < 1 else tri
        return dlo.union(dhi)

    def outside_star(self, zK):
        """True if every point of the ball zK lies strictly beyond some side of the heptagon star (inside the
        side circle |z - x0 e^{2 pi i k/7}| < R0, which excludes 0)."""
        R0 = (1 - self.p2 ** 2) / (2 * self.p2)
        for k in range(7):
            cen = acb(self.x0) * self.zeta ** k
            if abs(zK - cen).upper() < R0:
                return True
        return False

    def rS(self, th):
        c = th.cos()
        return self.x0 * c - (self.x0 ** 2 * c * c - 1).sqrt()


# Moebius matrices for the BFS over centres
def mob_mats(geo):
    def Tm(p):
        s = (1 - abs(p) ** 2).sqrt()
        return [[acb(1) / s, p / s], [p.conjugate() / s, acb(1) / s]]

    def mm(A, B):
        return [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
                [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]

    def inv2(A):
        return [[A[1][1], -A[0][1]], [-A[1][0], A[0][0]]]          # det = 1

    def rot0(th):
        return [[acb(0, th / 2).exp(), acb(0)], [acb(0), acb(0, -th / 2).exp()]]

    p = acb(geo.p2)
    Xm = mm(mm(Tm(p), rot0(PI)), inv2(Tm(p)))
    Y2m = mm(mm(Tm(geo.q4), rot0(PI)), inv2(Tm(geo.q4)))
    Cm = rot0(2 * PI / 7)
    return Xm, Y2m, Cm, mm


def centres_within(geo, D):
    """All points of the Delta-orbit of 0 within hyperbolic distance D of 0 (BFS over the star adjacency graph,
    restricted to centres within D + circumradius: a geodesic from 0 to a centre crosses a chain of stars that pairwise
    share a point, i.e. are side- or corner-adjacent, and whose centres stay within D + circumradius)."""
    Xm, Y2m, Cm, mm = mob_mats(geo)
    nbr = []
    Ck = [[acb(1), acb(0)], [acb(0), acb(1)]]
    for k in range(7):
        nbr.append(mm(Ck, Xm)); nbr.append(mm(Ck, Y2m))
        Ck = mm(Ck, Cm)
    lim = D + geo.d74
    I = [[acb(1), acb(0)], [acb(0), acb(1)]]
    found = {(0.0, 0.0): (acb(0), I)}
    queue = [I]
    while queue:
        g = queue.pop()
        for n in nbr:
            h = mm(g, n)
            c = h[0][1] / h[1][1]
            key = (round(float(c.real.mid()), 6), round(float(c.imag.mid()), 6))
            if key in found:
                continue
            if not (Geometry.dist(c, acb(0)) > lim):           # keep unless certainly beyond the limit
                found[key] = (c, h)
                queue.append(h)
    out = [v[0] for k, v in found.items() if k != (0.0, 0.0) and not (Geometry.dist(v[0], acb(0)) > D)]
    return out


# ------------------------------------------------------------------------------------------------ cutoff
class Cutoff:
    def __init__(self, r1, w):
        self.r1, self.w, self.r2 = arb(r1), arb(w), arb(r1) + arb(w)

    def eval(self, d):
        """bounds (h, |h'|, |h''|) as arbs on the interval d (h(d) = S((r2-d)/w), S quintic smoothstep).
        S is increasing; S' = 30 v^2 and S'' = 60 v (1-2x) with v = x(1-x), so the bounds below are sharp up to
        rounding; S' <= 15/8 and |S''| <= 10/sqrt(3) globally."""
        x = (self.r2 - d) / self.w
        if x.upper() <= 0:
            return arb(0), arb(0), arb(0)
        if x.lower() >= 1:
            return arb(1), arb(0), arb(0)
        lo, hi = x.lower().max(arb(0)), x.upper().min(arb(1))
        Sf = lambda y: y ** 3 * (10 - 15 * y + 6 * y * y)
        h = Sf(lo).union(Sf(hi))
        vlo, vhi = lo * (1 - lo), hi * (1 - hi)
        vmax = arb(1) / 4 if (lo <= arb(1) / 2 and hi >= arb(1) / 2) else vlo.max(vhi)
        S1 = (30 * vmax * vmax).min(arb(15) / 8)
        S2 = (60 * vmax * abs(1 - 2 * lo).max(abs(1 - 2 * hi))).min(10 / arb(3).sqrt() + arb("1e-9"))
        return h, S1.upper() / self.w, S2.upper() / self.w ** 2


# ------------------------------------------------------------------------------------------------ uniform radial bounds
class RadialBounds:
    """For x = u^2 <= xmax and every k >= 0:  g_lo <= R_k(u)/u^k <= G_hi,  |d/du (R_k/u^k)| <= 2u * gp(xmax)."""
    def __init__(self, lam, xmax, need_lo=True):
        t = (lam - arb(1) / 4).sqrt()
        absa = lam.sqrt()                                      # |a| = |1/2 - i t|
        self.absa = absa
        one_x = 1 - xmax
        e1 = one_x ** (-absa) - 1
        self.G_hi = one_x ** (arb(1) / 2 - absa)
        self.g_lo = one_x.sqrt() * ((t * (-one_x.log())).cos() - e1)
        self.gp = 2 * absa * one_x ** (-arb(1) / 2 - absa)
        assert self.g_lo > 0 or not need_lo

    def grad_coef(self, k, ur):
        """sup_{u <= ur} |grad(R_k(u) e^{i m th})|_hyp  (|m| = k)."""
        if k == 0:
            return ur * self.gp
        return k * ur ** (k - 1) * self.G_hi + ur ** (k + 1) * self.gp


# ------------------------------------------------------------------------------------------------ mismatch bounds
def mismatch_sup(tf, geo, kind, cut, log, N=64, u2=arb("0.24"), q=arb("0.6")):
    """Rigorous sup of |D| and |grad D|_hyp over W = {|z| <= u_circ, d(z, c) <= r2}, D = F_c - Ft, c = X(0) or Y^2(0)."""
    if kind == "X":
        cpt = geo.X(acb(0)); sig = perm21_index(tf.perm["a"]); gmap = geo.X
    else:
        cpt = geo.Y2(acb(0)); sig = perm21_index(mul(tf.perm["b"], tf.perm["b"])); gmap = geo.Y2
    uc, ur = q * u2, q * q * u2
    rho_r = 2 * ur.atanh()
    rb = RadialBounds(tf.lam, u2 * u2)
    Gh, gl = rb.G_hi, rb.g_lo
    # crude sup |Ft| on |z| <= U:  |Ft| <= G_hi(U^2) sum_m U^m ||coef_m||
    cn = [((sum(x * x for x in tf.A[m]) + sum(x * x for x in tf.B[m])).sqrt()) for m in range(tf.M + 1)]

    def Fbound(U):
        rbU = RadialBounds(tf.lam, U * U, need_lo=False)
        return rbU.G_hi * sum(cn[m] * U ** m for m in range(tf.M + 1))

    # covering cells (polar), adaptive until circumradius <= rho_r
    cells = []
    nu, nt = 7, 64
    stack = []
    for i in range(nu):
        for j in range(nt):
            u0 = geo.u_circ * i / nu; u1 = geo.u_circ * (i + 1) / nu
            t0 = -PI + 2 * PI * j / nt; t1 = -PI + 2 * PI * (j + 1) / nt
            stack.append((u0, u1, t0, t1))
    while stack:
        u0, u1, t0, t1 = stack.pop()
        U = u0.union(u1); Th = t0.union(t1)
        zK = acb(U * Th.cos(), U * Th.sin())
        if Geometry.dist(zK, cpt).lower() > cut.r2:
            continue
        if geo.outside_star(zK):
            continue                                   # every point of the cell is beyond some side of the star
        um, tm = (u0 + u1) / 2, (t0 + t1) / 2
        w = acb(um * tm.cos(), um * tm.sin())
        w = acb(arb(float(w.real.mid())), arb(float(w.imag.mid())))        # exact float centre
        cr = Geometry.dist(zK, w).upper()
        if cr <= rho_r:
            cells.append(w)
        else:
            if (u1 - u0) * 2 > (t1 - t0) * um:
                stack += [(u0, um, t0, t1), (um, u1, t0, t1)]
            else:
                stack += [(u0, u1, t0, tm), (u0, u1, tm, t1)]
    log(f"  mismatch {kind}: {len(cells)} disks of radius {float(rho_r):.4f}, N = {N}")
    SD, SdD = arb(0), arb(0)
    thk = [2 * PI * k / N for k in range(N)]
    ek = [acb(t.cos(), t.sin()) for t in thk]
    rho2 = 2 * u2.atanh()
    for wi, w in enumerate(cells):
        if wi % 25 == 0:
            log(f"    disk {wi}/{len(cells)}")
        # crude bound M2 on the outer circle
        U1 = ((Geometry.dist(w, acb(0)) + rho2) / 2).tanh()
        U2 = ((Geometry.dist(w, cpt) + rho2) / 2).tanh()
        M2 = Fbound(U1) + Fbound(U2)
        vals = []
        for k in range(N):
            p = Geometry.T(w, uc * ek[k])
            F = tf.at_point(p)
            Fg = tf.at_point(gmap(p))
            Fg = apply_perm(sig, Fg)
            vals.append([Fg[r] - F[r] for r in range(21)])
        Dsup, Dgrad = arb(0), arb(0)
        alias_c = M2 * Gh / gl * 2 / (1 - q ** N)
        for m in range(-(N // 2) + 1, N // 2):
            dh = [sum((vals[k][r] * ek[(-m * k) % N] for k in range(N)), acb(0)) / N for r in range(21)]
            dnorm = (sum(abs(x) ** 2 for x in dh)).sqrt()
            km = abs(m)
            coef = dnorm + alias_c * q ** (N - km)          # >= |d_m| R_km(uc)
            Dsup += coef * ur ** km * Gh / (uc ** km * gl)
            Dgrad += coef * rb.grad_coef(km, ur) / (uc ** km * gl)
        # tail |m| >= N/2
        r_ = q * q
        K = N // 2
        Dsup += 2 * M2 * Gh / gl * r_ ** K / (1 - r_)
        sk = r_ ** K * (K - (K - 1) * r_) / (1 - r_) ** 2           # sum_{k>=K} k r^k
        s0 = r_ ** K / (1 - r_)
        Dgrad += 2 * M2 / gl * (Gh * sk / ur + rb.gp * ur * s0)
        SD, SdD = SD.max(Dsup), SdD.max(Dgrad)
    log(f"  mismatch {kind}: sup|D| <= {float(SD.upper()):.3e},  sup|grad D|_hyp <= {float(SdD.upper()):.3e}")
    return SD, SdD


# ------------------------------------------------------------------------------------------------ isotype projectors
def wedge_projectors():
    """Exact sigma-isotypic projectors on wedge^2 R^21 (basis e<f of 2-subset pairs), sigma in
    {10+10b, 15, 21, 35}: Pi = (dim/2520) sum_g chi(g) wedge^2 P_g, as arb matrices (entries exact rationals)."""
    pairs2 = list(itertools.combinations(range(21), 2)); p2i = {p: i for i, p in enumerate(pairs2)}
    cidx = {}
    for i, c in enumerate(CL):
        for x in c:
            cidx[x] = i
    degs = [int(round(ch[0].real)) for ch in TABLE]
    want = {"10+10b": [i for i, d in enumerate(degs) if d == 10], "15": [degs.index(15)], "21": [degs.index(21)],
            "35": [degs.index(35)]}
    chars = {}
    for key, ids in want.items():
        vals = sum(TABLE[i] for i in ids)
        iv = np.rint(vals.real).astype(int)
        assert np.abs(vals - iv).max() < 1e-8
        chars[key] = iv
    acc = {key: np.zeros((210, 210), dtype=np.int64) for key in want}
    for g in A7:
        sig = perm21_index(g)
        rows, cols, sg = [], [], []
        for j, (e, f) in enumerate(pairs2):
            a, b = sig[e], sig[f]
            if a < b:
                rows.append(p2i[(a, b)]); sg.append(1)
            else:
                rows.append(p2i[(b, a)]); sg.append(-1)
            cols.append(j)
        for key in want:
            ch = chars[key][cidx[g]]
            if ch:
                acc[key][rows, cols] += ch * np.array(sg)
    # sanity: on wedge^2 V2 the four projectors sum to the identity (no 1 or 14 in wedge^2 14_(5,2))
    from hh_eval import PV2_rational
    PV = np.array([[float(x) for x in row] for row in PV2_rational()])
    W2P = np.array([[PV[e, g] * PV[f, h] - PV[e, h] * PV[f, g] for (g, h) in pairs2] for (e, f) in pairs2])
    tot = sum(acc[key] * (int(chars[key][cidx[tuple(range(7))]]) // len(want[key])) / 2520 for key in want)
    assert np.abs(tot @ W2P - W2P).max() < 1e-9
    out = {}
    for key in want:
        dim = int(chars[key][cidx[tuple(range(7))]])           # dimension of the (real) isotype sigma
        d1 = dim // len(want[key])                              # degree of each complex constituent
        Mi = acc[key]
        # exactness check: Pi^2 = Pi with Pi = (d1/2520) Mi, i.e. d1 * Mi @ Mi == 2520 * Mi (integers)
        assert np.array_equal(d1 * (Mi @ Mi), 2520 * Mi)
        out[key] = (dim, arb_mat(210, 210, [arb(int(Mi[r, s]) * d1) / 2520 for r in range(210) for s in range(210)]))
    return out, pairs2


# ------------------------------------------------------------------------------------------------ main
def certify(cls, coef_file, log, nu=48, nt=32, r1="1.37", w="0.25"):
    t_start = time.time()
    tf = TrialFunction(cls, coef_file)
    geo = Geometry()
    cut = Cutoff(r1, w)
    log(f"class {cls}: M = {tf.M}, lam* = {tf.lam.mid()}  (exact binary float), cutoff r1 = {r1}, w = {w}")
    assert geo.d74 < cut.r1
    # centres
    D = cut.r2 + geo.d74
    cents = centres_within(geo, D)
    dists = sorted(float(Geometry.dist(c, acb(0)).mid()) for c in cents)
    log(f"centres within r2 + circumradius = {float(D.mid()):.4f} of 0: {len(cents)} at distances "
        f"{sorted(set(round(x, 4) for x in dists))}")
    c1, cd = geo.X(acb(0)), geo.Y2(acb(0))
    types = []
    for c in cents:
        t = None
        for k in range(7):
            if abs(c - geo.rot(c1, k)) < arb("1e-20"):
                t = ("X", k)
            if abs(c - geo.rot(cd, k)) < arb("1e-20"):
                t = ("Y2", k)
        assert t is not None, "unexpected centre"
        types.append(t)
    # mismatch bounds
    SD = {}
    cache = os.path.join(HERE, f"hh_mismatch_cls{cls}.txt")
    key = f"{os.path.basename(coef_file)} r1={r1} w={w}"
    if os.path.exists(cache) and open(cache).readline().strip() == key:
        lines = open(cache).read().split("\n")[1:]
        for ln in lines:
            if ln.strip():
                kind, a1, a2 = ln.split("|")
                SD[kind] = (arb(a1), arb(a2))
        log(f"  mismatch bounds read from {os.path.basename(cache)} (computed by this script; delete to recompute)")
        for kind in SD:
            log(f"  mismatch {kind}: sup|D| <= {float(SD[kind][0].upper()):.3e},  sup|grad D|_hyp <= {float(SD[kind][1].upper()):.3e}")
    else:
        for kind in ("X", "Y2"):
            SD[kind] = mismatch_sup(tf, geo, kind, cut, log)
        with open(cache, "w") as fh:
            fh.write(key + "\n")
            for kind, (a1, a2) in SD.items():
                fh.write(f"{kind}|{a1.upper().str(30, radius=True)}|{a2.upper().str(30, radius=True)}\n")
    log(f"  [{time.time() - t_start:.0f}s]")
    # box pass over the sector |th| <= pi/7, u <= rS(th)
    proj, pairs2 = wedge_projectors()
    log("  wedge^2 isotype projectors built (exact idempotence checked)")
    us = [geo.u_circ * i / nu for i in range(nu + 1)]
    ts = [-PI / 7 + 2 * PI / 7 * j / nt for j in range(nt + 1)]
    trig_box = [tf.trig(ts[j].union(ts[j + 1])) for j in range(nt)]
    trig_mid = [tf.trig((ts[j] + ts[j + 1]) / 2) for j in range(nt)]
    I_psi_lo, I_res_hi = arb(0), arb(0)
    grad2_max = arb(0)                      # sup conf * |grad Psi|_euc^2   (before dividing by c)
    eps_psi, eps_dpsi = arb(0), arb(0)      # sup |Psi - Ft|, sup |grad(Psi - Ft)|_hyp
    Jc_cols = []                            # centre J values (210 x 2 per box)
    box_info = []                           # (euclid area, Delta_J, |Ft(c)|, |grad Ft(K)|_euc)
    nbox = 0
    for i in range(nu):
        u0, u1 = us[i], us[i + 1]
        U = u0.union(u1); um = (u0 + u1) / 2
        rad_box = tf.radial(U, second=True)
        rad_mid = tf.radial(um, deriv=True)
        conf_hi = (1 - u0 * u0) ** 2 / 4
        for j in range(nt):
            t0, t1 = ts[j], ts[j + 1]
            Th = t0.union(t1)
            rS = geo.rS(Th)
            if not (u0 < rS.upper()):
                continue                                         # box outside the sector
            inside = bool(u1 <= rS.lower())
            nbox += 1
            area_h = (t1 - t0) * 2 * (1 / (1 - u1 * u1) - 1 / (1 - u0 * u0))
            area_e = (t1 - t0) * (u1 * u1 - u0 * u0) / 2
            Fc, Fxc, Fyc = tf.value_grad(rad_mid, trig_mid[j])
            Fu_K, Fthu_K, Gxu_K, Gyu_K, Gxt_K, Gyt_K = tf.polar_derivs(rad_box, trig_box[j])
            nrm = lambda *vs: sum((up(x) ** 2 for v in vs for x in v), arb(0)).sqrt()
            du, dth = u1 - u0, t1 - t0
            dF = (nrm(Fu_K) * du + u1 * nrm(Fthu_K) * dth) / 2          # sup_K |F - F(c)|  (mean value, polar box)
            dG = (nrm(Gxu_K, Gyu_K) * du + nrm(Gxt_K, Gyt_K) * dth) / 2  # sup_K |G - G(c)|, G = Euclidean gradient
            nFc = sum(x * x for x in Fc).sqrt()
            nGc = sum((x * x for x in Fxc + Fyc), arb(0)).sqrt()
            gK = nGc + dG
            # --- partition of unity and residual on K
            zK = acb(U * Th.cos(), U * Th.sin())
            z0 = acb(um * ((t0 + t1) / 2).cos(), um * ((t0 + t1) / 2).sin())
            rK = Geometry.dist(zK, z0).upper()                      # hyperbolic circumradius of K about z0
            act = []
            for c, tp in zip(cents, types):
                d0 = Geometry.dist(z0, c)
                dK = (d0.lower() - rK).max(arb(0)).union(d0.upper() + rK)
                h, h1, h2 = cut.eval(dK)
                if h1 == 0 and h.upper() == 0:
                    continue
                cth = 1 / dK.lower().tanh()
                act.append((h.upper(), h1, h2 + h1 * cth, tp[0]))
            Sh1 = sum((a[1] for a in act), arb(0))
            SLh = sum((a[2] for a in act), arb(0))
            res, e0, e1 = arb(0), arb(0), arb(0)
            for (hu, h1, Lh, kind) in act:
                gchi = h1 + hu * Sh1
                lchi = Lh + 2 * h1 * Sh1 + hu * SLh + 2 * hu * Sh1 * Sh1
                sD, sdD = SD[kind]
                res += lchi * sD + 2 * gchi * sdD
                e0 += hu * sD
                e1 += gchi * sD + hu * sdD
            I_res_hi += area_h * res * res
            eps_psi, eps_dpsi = eps_psi.max(e0), eps_dpsi.max(e1)
            if inside:
                lo = (nFc - dF - e0)
                if lo > 0:
                    I_psi_lo += area_h * lo * lo
            # gradient (Euclidean): |grad Psi|_euc <= |grad Ft|_euc + (2/(1-u^2)) |grad(Psi-Ft)|_hyp
            gpsi = gK + 2 / (1 - u1 * u1) * e1
            grad2_max = grad2_max.max(conf_hi * gpsi * gpsi)
            # --- Jacobian flux J_{ef} = 1/2 (F_e rot grad F_f - F_f rot grad F_e), rot grad = (F_y, -F_x)
            def Jvals(F, Fx, Fy):
                jx, jy = [], []
                for (e, f) in pairs2:
                    jx.append((F[e] * Fy[f] - F[f] * Fy[e]) / 2)
                    jy.append(-(F[e] * Fx[f] - F[f] * Fx[e]) / 2)
                return jx, jy
            jxc, jyc = Jvals(Fc, Fxc, Fyc)
            # |J(z) - J(c)| <= 1/2 (|F(z)-F(c)| sup|G| + |F(c)| |G(z)-G(c)|)    (|J(F,G)| <= |F||G|/2)
            dJ = (dF * gK + nFc * dG) / 2
            # Psi vs Ft:  |J_Psi - J_Ft| <= 1/2 (|Psi - Ft| |grad Psi|_euc + |Ft| |grad(Psi - Ft)|_euc)
            dJ += (e0 * gpsi + (nFc + dF) * 2 / (1 - u1 * u1) * e1) / 2
            Jc_cols.append((jxc, jyc))
            box_info.append((area_e, dJ))
        log(f"  strip {i + 1}/{nu} done [{time.time() - t_start:.0f}s]") if (i + 1) % 8 == 0 else None
    log(f"  boxes: {nbox}")
    # ---------------------------------------------------------------- assemble
    c_lo = 2520 * I_psi_lo / 14                                  # int_C Psi_e Psi_f = c P_V2[e,f]
    eta_s = (I_res_hi / I_psi_lo).sqrt()
    lam_h_hi = tf.lam + eta_s
    lam_h_lo = (tf.lam - eta_s).max(LAMBDA1_LO)
    rho = 2 * eta_s / LAMBDA1_LO.sqrt()
    delta = rho * B_GAP.sqrt() / (B_GAP - lam_h_hi)
    a = (lam_h_lo - rho * (lam_h_hi / (1 - delta * delta)).sqrt()).max(LAMBDA1_LO)
    nu_ = B_GAP - (B_GAP - a) * delta * delta
    Lam_h = grad2_max / c_lo
    Gam = (AREA * Lam_h).sqrt()
    log(f"int_sector |Psi|^2 dA >= {I_psi_lo.lower()},  int_sector |res|^2 dA <= {I_res_hi.upper()}")
    log(f"eta_s <= {float(eta_s.upper()):.3e},  rho <= {float(rho.upper()):.3e},  delta <= {float(delta.upper()):.3e}")
    log(f"lambda_h in [{float(lam_h_lo.lower()):.10f}, {float(lam_h_hi.upper()):.10f}],  lambda_1 >= a = {float(a.lower()):.10f}")
    log(f"eps |Psi - Ft| <= {float(eps_psi.upper()):.2e}, |grad(Psi - Ft)| <= {float(eps_dpsi.upper()):.2e}")
    log(f"Lambda_h <= {float(Lam_h.upper()):.6e},  Gamma_h <= {float(Gam.upper()):.5f}")
    # isotype integrals
    nb = len(Jc_cols)
    Jmat = arb_mat(210, 2 * nb, [Jc_cols[col // 2][col % 2][r] for r in range(210) for col in range(2 * nb)])
    phis = {}
    for key, (dim, P) in proj.items():
        PJ = P * Jmat
        tot = arb(0)
        for bi in range(nb):
            nrm = (sum(PJ[r, 2 * bi] ** 2 + PJ[r, 2 * bi + 1] ** 2 for r in range(210))).sqrt()
            area_e, dJ = box_info[bi]
            tot += area_e * (nrm + dJ) ** 2
        phis[key] = 2520 * tot / (dim * c_lo * c_lo)
        log(f"  isotype {key:>6} (dim {dim}): ||J_omega||^2 / ||omega||^2 <= {float(phis[key].upper()):.4e}")
    phi_max = amax(phis.values())
    B_h = B_GAP / (B_GAP - lam_h_hi) * phi_max
    log(f"B_h = b/(b - lambda_h) * max_sigma phi_sigma <= {float(B_h.upper()):.4e}")
    # ---------------------------------------------------------------- final inequality (Theorem 5.3 form)
    d = TAU - a
    one_m = 1 - a / nu_
    R = (rho + (rho * rho + one_m * d).sqrt()) / one_m
    k = (B_h / 3).sqrt()
    lhs = a * (1 - R * R / nu_) - rho * R
    rhs = 4 * k * AREA.sqrt() * (d + 2 * rho * R).sqrt() + Gam * R * R / nu_.sqrt()
    ok = bool(lhs > rhs)
    log(f"R <= {float(R.upper()):.6f};  left side >= {float(lhs.lower()):.6f};  right side <= {float(rhs.upper()):.6f}")
    log(f"==> class {cls}: degree <= 24 excluded: {ok}   [{time.time() - t_start:.0f}s]")
    return ok


if __name__ == "__main__":
    cls = int(sys.argv[1])
    coef = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, f"coef_cls{cls}_M90.npz")
    nu = int(sys.argv[3]) if len(sys.argv) > 3 else 48
    nt = int(sys.argv[4]) if len(sys.argv) > 4 else 32
    certify(cls, coef, lambda s: print(s, flush=True), nu=nu, nt=nt)
