"""The tau-stratum below 42: no C(tau)-symmetric pencil of degree 30 (Chapter 7, §7.14) [X].

Setting.  tau = (12)(34), K = C(tau) (order 24), N = <tau>, H = K/N = C2 x S3 acting on the pencil as the dihedral
group D6 of P^1.  U = <tau, (13)(24), (567)> is the preimage of the rotation group C6.  The curves
    W = C/(S4 x S3)^+   (3-sets of roots; genus 0),
    T = C/K             (3-set + pairing of the complementary 4-set; genus 3, a cubic cover of W),
    X_rot = C/U         (= T(sqrt disc A), A the cubic of the 3-set; genus 10)
are built over F_p from the septic Belyi model of §6.5 (frobenius_lpolys.model):
  * W: A(z) = z^3 - e1 z^2 + e2 z - e3 divides q(z) - c; eliminating e3 gives a plane model G(e1, e2) = 0; a degree-1
    function w on W (Riemann-Roch) parametrises it, and e1, e2, e3, t become rational functions of w;
  * T: the cubic resolvent of the complementary quartic B(z) = 7 (q(z) - c)/A(z), theta = (b1 + b2)(b3 + b4).
A pencil of the stratum, of degree 2k on C, is a function y on T of degree k with (Lemma 7.24)
    y^*(oo) = F (= 3 P1 + 3 P2 + 3 P3 + 6 E, or with 9 at one P_j),  y^*(0) = A + 2 E_a,  y^*(1) = B + 2 E_b,
A = the nine points with inertia R_a, B = the R_b point, P_j = the z-points.  check(F, A', B') decides over the
algebraic closure of F_p whether such y exists: y and y - lambda range over L(F - A') and L(F - B') + constants;
their norms to the w-line, with the fixed part removed, must be squares of binary forms of the right degree
(Groebner bases, chart by chart).  This is a necessary condition; an empty ideal excludes the configuration.

Modes (p prime, sign = +-1 the two curves; default: all):
  model       genus checks; identification of P1, P2, P3, A, B; odd part of div(disc A) = A + B (so X_rot -> T is
              branched exactly there, g(X_rot) = 10); rationality of the preimages of P_j in X_rot
  control     planted dihedral-norm instance must be found, a wrong one excluded
  fixed       F = 3P + 6P_j (types (3,3,9))
  free        F = 3P + 6Q for every F_p-rational point Q of T (Lemma 7.26 reduces the free type to these)
  degenerate  F = 3P with the ramification points split A' + B' (Q specialising into Fix(iota))
Requires Sage (pip install passagemath-singular passagemath-schemes passagemath-pari passagemath-flint).
Usage: python3 tau_stratum.py p sign [mode ...]
"""
import sys, time, itertools
from sage.all__sagemath_schemes import *
import sympy as sp
from frobenius_lpolys import model


# ------------------------------------------------------------------------------------------------ the models
def build(p, sign):
    rt = int(sp.sqrt_mod(21, p))
    s_, t_, r_, qc = model(p, sign * rt % p)
    K = GF(p)
    R = PolynomialRing(K, ['e1', 'e2', 'e3', 'z'])
    e1, e2, e3, z = R.gens()
    q = sum(K(a) * z**j for j, a in enumerate(qc))
    A = z**3 - e1 * z**2 + e2 * z - e3
    rem = q
    while rem.degree(z) >= 3:
        d = rem.degree(z)
        rem = rem - rem.coefficient({z: d}) * z**(d - 3) * A
    R2, R1, R0 = (rem.coefficient({z: j}) for j in (2, 1, 0))
    a3, b3 = R2.coefficient({e3: 1}), R2.coefficient({e3: 0})
    assert R2 == a3 * e3 + b3                                  # e3 is rational in e1, e2
    G = R2.resultant(R1, e3)
    F1 = FunctionField(K, 'E1'); E1 = F1.gen()
    Y = F1['Y'].gen()
    Gp = sum(F1(G.coefficient({e2: j})(E1, 0, 0, 0)) * Y**j for j in range(G.degree(e2) + 1))
    W = F1.extension(Gp / Gp.leading_coefficient(), 'E2'); E2 = W.gen()
    gW = W.genus()
    pl = [P for P in W.places_finite(1) if P.degree() == 1]
    w = [b for b in pl[0].divisor().basis_function_space() if b.minimal_polynomial('U').degree() > 1][0]
    Fw = FunctionField(K, 'w'); wv = Fw.gen()
    Rw = PolynomialRing(K, ['E', 'U']); Ev, Uv = Rw.gens()
    cs = w.minimal_polynomial('U').list()
    den = lcm([c.denominator() for c in cs])
    Mp = sum(Rw((c * den).numerator()(Ev)) * Uv**j for j, c in enumerate(cs))
    assert Mp.degree(Ev) == 1                                  # w has degree 1 on W
    E1w = -Fw(Mp.coefficient({Ev: 0})(0, wv)) / Fw(Mp.coefficient({Ev: 1})(0, wv))
    Xz = Fw['X'].gen()
    Gw = sum(Fw(G.coefficient({e2: j})(E1w, 0, 0, 0)) * Xz**j for j in range(G.degree(e2) + 1))
    E2w = None
    for cand, _ in Gw.roots():
        good = True
        for P in pl[1:8]:
            k_, fr, to = P.residue_field()
            wP = to(w)
            if cand.denominator()(wP) != 0 and to(E2) != cand.numerator()(wP) / cand.denominator()(wP):
                good = False
        if good:
            E2w = cand
    E3w = -b3(E1w, E2w, 0, 0) / a3(E1w, E2w, 0, 0)
    Cw = R0(E1w, E2w, E3w, 0)                                 # c = r (1 - t)
    Tw = 1 - Cw / K(r_)
    zz = Fw['zz'].gen()
    qw = sum(Fw(K(a)) * zz**j for j, a in enumerate(qc))
    Aw = zz**3 - E1w * zz**2 + E2w * zz - E3w
    Bw, rm = (7 * (qw - Cw)).quo_rem(Aw)
    assert rm == 0 and Bw.degree() == 4 and Bw.leading_coefficient() == 1
    s1, s2, s3, s4 = -Bw[3], Bw[2], -Bw[1], Bw[0]
    Th = Fw['Th'].gen()
    Res = Th**3 - 2 * s2 * Th**2 + (s2**2 + s1 * s3 - 4 * s4) * Th - (s1 * s2 * s3 - s1**2 * s4 - s3**2)
    assert Res.is_irreducible()
    info = dict(gW=gW, deg_t=max(Tw.numerator().degree(), Tw.denominator().degree()))
    return K, qc, Fw, (E1w, E2w, E3w), Tw, Res, info


def resolvent_identities():
    """theta_i = pairing sums for a quartic: e1, e2, e3 of the theta_i in terms of the roots' symmetric functions."""
    import random
    for _ in range(20):
        b = [random.randint(-50, 50) for _ in range(4)]
        th = [(b[0] + b[1]) * (b[2] + b[3]), (b[0] + b[2]) * (b[1] + b[3]), (b[0] + b[3]) * (b[1] + b[2])]
        S1 = sum(b); S2 = sum(b[i] * b[j] for i, j in itertools.combinations(range(4), 2))
        S3 = sum(b[i] * b[j] * b[k] for i, j, k in itertools.combinations(range(4), 3)); S4 = b[0] * b[1] * b[2] * b[3]
        assert sum(th) == 2 * S2
        assert th[0] * th[1] + th[0] * th[2] + th[1] * th[2] == S2**2 + S1 * S3 - 4 * S4
        assert th[0] * th[1] * th[2] == S1 * S2 * S3 - S1**2 * S4 - S3**2


class TModel:
    def __init__(self, p, sign):
        resolvent_identities()
        self.p, self.sign = p, sign
        K, qc, Fw, (E1w, E2w, E3w), Tw, Res, self.info = build(p, sign)
        self.K, self.Fw = K, Fw
        Tff = self.Tff = Fw.extension(Res, 'th')
        th = self.th = Tff.gen()
        self.genus = Tff.genus()
        tT = self.tT = Tff(Tw)
        z = PolynomialRing(K, 'z').gen()
        q = sum(K(a) * z**j for j, a in enumerate(qc))
        quad = q.derivative() // (z**3 * (z - 1))
        quad = quad / quad.leading_coefficient()
        s_, t_ = -quad[1], quad[0]                             # z^2 - s z + t: the double roots of q - r
        e1, e2, e3 = E1w, E2w, E3w
        self.discA = e1**2 * e2**2 - 4 * e2**3 - 4 * e1**3 * e3 + 18 * e1 * e2 * e3 - 27 * e3**2

        def spec(P):                                           # A(z) and theta at the place P
            k_, fr, to = P.residue_field()
            zk = PolynomialRing(k_, 'z').gen()
            v = [to(Tff(f)) for f in (E1w, E2w, E3w)]
            return zk**3 - v[0] * zk**2 + v[1] * zk - v[2], to(th), k_
        D0, D1 = tT.divisor(), (tT - 1).divisor()
        P12, tau_pt, P3, Bpt, Aset = [], None, None, None, []
        for P in D0.support():                                 # over t = 0 (c0 = (12)(34)-type)
            if D0.multiplicity(P) != 1:
                continue
            Aq, thv, k_ = spec(P)
            if Aq.is_squarefree() and Aq.gcd(quad.change_ring(k_)(Aq.parent().gen())) == 1:
                # the 3-set is the three simple roots: stabiliser tau (theta = 4 d1 d2) or z (theta = (d1 + d2)^2)
                if thv == k_(s_**2):
                    P12.append(P)
                elif thv == k_(4 * t_):
                    tau_pt = P
                else:
                    raise ValueError("unexpected theta")
            else:                                              # A -> (z - delta)^2 (z - sigma): inertia R_a
                Aset.append(P)
        for P in D1.support():                                 # over t = 1 (c1 = (1234)(56)-type)
            m = D1.multiplicity(P)
            if m not in (1, 2):
                continue
            Aq, thv, k_ = spec(P)
            if m == 1:                                         # stabiliser C4 inside K: inertia R_b
                Bpt = P
            elif dict(Aq.roots()).get(k_(1), 0) == 2:          # h^2 = z: the third z-point
                P3 = P
            else:                                              # h^2 in R_a, h outside K
                Aset.append(P)
        assert sum(P.degree() for P in P12) == 2 and tau_pt and P3 and Bpt
        assert sum(P.degree() for P in Aset) == 9
        self.P12, self.P3, self.B, self.Aset, self.tau_pt = P12, P3, Bpt, Aset, tau_pt
        self.P = (P12 + [P3]) if len(P12) == 2 else None
        self.Adiv = sum(P.divisor() for P in Aset)
        self.Rp = PolynomialRing(K, 'w')

    def mat(self, f):
        rows = []
        for b in (self.Tff(1), self.th, self.th**2):
            v = (f * b).list()
            rows.append([v[i] for i in range(3)])
        return matrix(self.Fw, rows)

    def check(self, F, Aprime, Bprime):
        """(possible?, info): can some y have y^*(oo) <= F, div(y) + F - A' = 2E_a, div(y - lam) + F - B' = 2E_b?"""
        K, Tff, Rp = self.K, self.Tff, self.Rp
        LF, fromLF, toLF = F.function_space()
        Vb = (F - Aprime).basis_function_space()
        Ub = (F - Bprime).basis_function_space() if Bprime.degree() > 0 else LF.basis()
        if len(Vb) == 0:
            return False, "L(F - A') = 0"
        one = toLF(Tff(1))
        Wsp = LF.subspace([toLF(v) for v in Vb]).intersection(LF.subspace([toLF(u) for u in Ub] + [one]))
        if Wsp.dimension() == 0:
            return False, "V' = 0"
        Mu = matrix(K, [toLF(u) for u in Ub] + [one])
        ys, lams = [], []
        for wv in Wsp.basis():                                 # y = u + lam, u in L(F - B')
            lams.append(Mu.solve_left(wv)[-1])
            ys.append(fromLF(wv))
        nV = len(ys)
        Ms = [self.mat(y) for y in ys]
        den = lcm([e.denominator() for M in Ms for e in M.list()])
        Rx = PolynomialRing(Rp, [f'x{i}' for i in range(nV)])
        xs = Rx.gens()
        Mx = sum(xs[i] * matrix(Rx, 3, 3, [Rx(Rp((e * den).numerator())) for e in Ms[i].list()]) for i in range(nV))
        MxB = Mx - sum(xs[i] * Rx(Rp(den)) * K(lams[i]) for i in range(nV)) * identity_matrix(Rx, 3)

        def variable_part(nrm):
            cs = nrm.coefficients()
            g = cs[0]
            for c in cs[1:]:
                g = gcd(g, c)
            return Rx(sum((c // g) * mono for c, mono in zip(nrm.coefficients(), nrm.monomials())))
        ntA, ntB = variable_part(Mx.det()), variable_part(MxB.det())
        De, Db = F.degree() - Aprime.degree(), F.degree() - Bprime.degree()
        if De % 2 or Db % 2:
            return False, "parity"
        degA = max(c.degree() for c in ntA.coefficients()); degB = max(c.degree() for c in ntB.coefficients())
        if (degA, degB) != (De, Db):                           # a fixed part at w = oo: not decided this way
            return True, f"DEGREE-DROP ({degA}/{De}, {degB}/{Db}): not decided"
        ea, eb = De // 2, Db // 2
        names = [f'x{i}' for i in range(nV)] + [f's{j}' for j in range(ea + 1)] + [f'u{j}' for j in range(eb + 1)]
        Rall = PolynomialRing(K, names, order='degrevlex')
        X = Rall.gens()[:nV]; S = Rall.gens()[nV:nV + ea + 1]; U = Rall.gens()[nV + ea + 1:]

        def eqs(nt_, Sv, D_):                                  # coefficients of nt_ = (sum Sv_j w^j)^2
            out = []
            for j in range(D_ + 1):
                cj = sum(Rall(c[j]) * prod(X[i]**ee for i, ee in enumerate(mono.exponents()[0]))
                         for c, mono in zip(nt_.coefficients(), nt_.monomials()))
                out.append(cj - sum(Sv[a] * Sv[b] for a in range(len(Sv)) for b in range(len(Sv)) if a + b == j))
            return out
        eqA, eqB = eqs(ntA, S, De), eqs(ntB, U, Db)
        info = []
        for chart in range(nV):                                # y in P(V'): charts x_chart = 1, x_i = 0 for i < chart
            fix = [X[i] for i in range(chart)] + [X[chart] - 1]
            RA = PolynomialRing(K, names[:nV + ea + 1], order='degrevlex')
            IA = RA.ideal([RA(str(e_)) for e_ in eqA + fix])
            dA = IA.dimension()
            if dA < 0:
                info.append(f"x{chart}:(a)empty")
                continue
            RXU = PolynomialRing(K, names[:nV] + names[nV + ea + 1:], order='degrevlex')
            if dA == 0:
                J = IA.elimination_ideal(list(RA.gens()[nV:]))
                IB = RXU.ideal([RXU(str(g)) for g in J.gens()] + [RXU(str(e_)) for e_ in eqB + fix])
            else:
                IB = Rall.ideal(eqA + eqB + fix)
            dB = IB.dimension()
            info.append(f"x{chart}:(a)dim{dA},(a+b)dim{dB}")
            if dB >= 0:
                return True, "; ".join(info)
        return False, "; ".join(info)


# ------------------------------------------------------------------------------------------------ the checks
def run(p, sign, modes):
    t0 = time.time()
    M = TModel(p, sign)
    Tff = M.Tff
    zero = Tff(1).divisor()
    print(f"\n=== p = {p}, sign = {sign:+d} ===")
    print(f"  genus W = {M.info['gW']}, deg t(w) = {M.info['deg_t']}, genus T = {M.genus}; "
          f"{len(Tff.places(1))} rational places on T; A = places of degrees {[P.degree() for P in M.Aset]}")
    if M.P is None:
        print("  P1, P2 are conjugate over F_p: use another prime for this curve")
        return
    Pd = [P.divisor() for P in M.P]
    Bd = M.B.divisor()
    if "model" in modes:
        dA = Tff(M.discA).divisor()
        odd = sum((P.divisor() for P in dA.support() if dA.multiplicity(P) % 2), zero)
        print("  odd part of div(disc A) == A + B:", odd == M.Adiv + Bd, " (X_rot -> T branched at these 10 points, "
              "g(X_rot) = 2*3 - 1 + 5 = 10)")
        for j, P in enumerate(M.P):
            f = Tff(M.discA)
            v = f.valuation(P)
            k_, fr, to = P.residue_field()
            lc = to(f / P.local_uniformizer()**v)
            print(f"  P{j + 1}: ord disc A = {v}, leading coefficient a square: {lc.is_square()} "
                  f"(preimages in X_rot rational)")
    if "control" in modes:
        w = Tff(M.Fw.gen())
        set_random_seed(1)
        La = (2 * M.tau_pt.divisor() + Pd[0]).basis_function_space()
        a = sum(M.K.random_element() * f for f in La)
        b = sum(M.K.random_element() * f for f in La[:2])
        for c in range(1, p):
            fib = (w - c).divisor()
            zs = [P for P in fib.support() if fib.multiplicity(P) > 0]
            if len(zs) == 3 and all(P.degree() == 1 for P in zs):
                break
        y = a**2 / (a**2 - b**2 * (w - c))                     # a dihedral norm: y = Tr(g)^2 / 4 Nm(g)
        Dy, D1 = y.divisor(), (y - 1).divisor()
        F = -sum((Dy.multiplicity(P) * P.divisor() for P in Dy.support() if Dy.multiplicity(P) < 0), zero)
        Ap = sum((P.divisor() for P in Dy.support() if Dy.multiplicity(P) > 0 and Dy.multiplicity(P) % 2), zero)
        Bp = sum((P.divisor() for P in D1.support() if D1.multiplicity(P) > 0 and D1.multiplicity(P) % 2), zero)
        ok, info = M.check(F, Ap, Bp)
        print(f"  control: planted y of degree {F.degree()} ({Ap.degree()}, {Bp.degree()} odd points): "
              f"{'found' if ok else 'MISSED'} [{info}]")
        for c2 in range(c + 1, p):
            fib = (w - c2).divisor()
            zs = [P for P in fib.support() if fib.multiplicity(P) > 0]
            if (len(zs) == 3 and all(P.degree() == 1 for P in zs)
                    and not any(P in Bp.support() or P in F.support() for P in zs)):
                break
        ok2, info2 = M.check(F, Ap, sum((P.divisor() for P in zs), zero))
        print(f"  control: same F with a wrong B': {'excluded' if not ok2 else 'NOT excluded'} [{info2}]")
    if "fixed" in modes:
        for mults in [(3, 3, 9), (3, 9, 3), (9, 3, 3)]:
            F = sum((m * D for m, D in zip(mults, Pd)), zero)
            ok, info = M.check(F, M.Adiv, Bd)
            print(f"  fixed type {mults}: {'POSSIBLE' if ok else 'excluded'} [{info}]  ({time.time() - t0:.0f}s)",
                  flush=True)
    if "free" in modes:
        skip = set(M.P) | set(M.Aset) | {M.B}
        Qs = [Q for Q in Tff.places(1) if Q not in skip]
        bad = 0
        for i, Q in enumerate(Qs):
            F = 3 * (Pd[0] + Pd[1] + Pd[2]) + 6 * Q.divisor()
            ok, info = M.check(F, M.Adiv, Bd)
            bad += ok
            if ok:
                print(f"  free point {i}: POSSIBLE [{info}]", flush=True)
        print(f"  free type: {len(Qs)} rational points Q, {'all excluded' if bad == 0 else f'{bad} NOT excluded'}"
              f"  ({time.time() - t0:.0f}s)", flush=True)
    if "degenerate" in modes:
        orbits = M.Aset + [M.B]
        n_ok, n_all = 0, 0
        for r in range(len(orbits) + 1):
            for sub in itertools.combinations(range(len(orbits)), r):
                Ap = [orbits[i] for i in sub]; Bp = [orbits[i] for i in range(len(orbits)) if i not in sub]
                na, nb = sum(P.degree() for P in Ap), sum(P.degree() for P in Bp)
                if na % 2 == 0 or nb % 2 == 0 or na > 9 or nb > 9:
                    continue
                n_all += 1
                ok, info = M.check(3 * (Pd[0] + Pd[1] + Pd[2]), sum((P.divisor() for P in Ap), zero),
                                   sum((P.divisor() for P in Bp), zero))
                n_ok += ok
                if ok:
                    print(f"  degenerate split {[P.degree() for P in Ap]} | {[P.degree() for P in Bp]}: POSSIBLE [{info}]")
        print(f"  degenerate k = 9: {n_all} splits of A + B, {'all excluded' if n_ok == 0 else f'{n_ok} NOT excluded'}"
              f"  ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    p, sign = int(sys.argv[1]), int(sys.argv[2])
    modes = sys.argv[3:] or ["model", "control", "fixed", "degenerate", "free"]
    run(p, sign, modes)
