"""FRAMEWORK_CONFORMAL.md section E (numbers of Proposition E1).  Numerical (non-rigorous) estimate of Lambda^G = sup_h lambda_1(h g) Area(h g) over A7-invariant
conformal factors h on a (2,4,7) A7-curve.  h is piecewise constant on the n^2 cells of the reference
triangle, separately on U- and L-type tiles.  lambda_1(h g) = min over the twisted problems
Q0 (second eigenvalue), Q1, Q2 (first eigenvalues).  P1 elements (upper bounds, fine for estimates)."""
import sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orbifold import (reference_triangle, quotient_tiles, build_dofs, ref_elements, numeric_elem_mats)
from certify import QUOTIENTS
from triples_data import triples

n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
tri = int(sys.argv[2]) if len(sys.argv) > 2 else 0
iters = int(sys.argv[3]) if len(sys.argv) > 3 else 40
XA, XB, XC = reference_triangle()
els = ref_elements(n)
NE = len(els)
mats = numeric_elem_mats(n, XA, XB, XC, "P1")
cell_area = np.array([Me.sum() for (Ke, Me) in mats])          # hyperbolic area of each cell
a, b, _ = triples[tri]


def build(which):
    Kg, eps = ([a, b], [1, 1]) if which == "Q0" else QUOTIENTS[which]
    reps, glue, order = quotient_tiles(a, b, Kg, eps)
    uf, zero = build_dofs(n, len(reps), glue, "P1")
    zero_roots = {uf.find(z)[0] for z in zero}
    index, R, C, KV, MV, CL = {}, [], [], [], [], []
    for tile in range(2 * len(reps)):
        typ = tile % 2
        for e, trg in enumerate(els):
            ids, sg = [], []
            for v in trg:
                r, s = uf.find((tile, (v,)))
                if r in zero_roots:
                    ids.append(-1); sg.append(0); continue
                if r not in index:
                    index[r] = len(index)
                ids.append(index[r]); sg.append(s)
            Ke, Me = mats[e]
            for p in range(3):
                for q in range(3):
                    if ids[p] < 0 or ids[q] < 0:
                        continue
                    R.append(ids[p]); C.append(ids[q])
                    KV.append(sg[p] * sg[q] * Ke[p, q]); MV.append(sg[p] * sg[q] * Me[p, q])
                    CL.append(typ * NE + e)
    N = len(index)
    R, C, CL = np.array(R), np.array(C), np.array(CL)
    K = sp.csc_matrix((KV, (R, C)), shape=(N, N))
    return dict(K=K, R=R, C=C, MV=np.array(MV), CL=CL, N=N, name=which, ntiles=2 * len(reps))


def Mof(P, h):
    return sp.csc_matrix((P["MV"] * h[P["CL"]], (P["R"], P["C"])), shape=(P["N"], P["N"]))


def lowest(P, h, k, skip=0):
    M = Mof(P, h)
    vals, vecs = sla.eigsh(P["K"], k=k + skip, M=M, sigma=-1e-3, which="LM")
    o = np.argsort(vals)
    vals, vecs = vals[o][skip:], vecs[:, o][:, skip:]
    out = []
    for lam, v in zip(vals, vecs.T):
        nrm = v @ (M @ v)
        # d lambda / d h_c = -lambda * v^T M_c v / v^T M v
        contrib = P["MV"] * v[P["R"]] * v[P["C"]]
        g = -lam * np.bincount(P["CL"], weights=contrib, minlength=2 * NE) / nrm
        out.append((lam, g))
    return out


t0 = time.time()
Ps = [build(w) for w in ("Q0", "Q1", "Q2")]
print(f"n={n} triple={tri}: dofs", [P["N"] for P in Ps], f"({time.time()-t0:.0f}s)", flush=True)
area_c = np.concatenate([cell_area, cell_area])      # per class (U cells, L cells), orbifold area
A0 = area_c.sum()                                    # = 2 * 3pi/28


def evaluate(h):
    branches = []
    branches += [("Q0", l, g) for (l, g) in lowest(Ps[0], h, 1, skip=1)]
    branches += [("Q1", l, g) for (l, g) in lowest(Ps[1], h, 3)]
    branches += [("Q2", l, g) for (l, g) in lowest(Ps[2], h, 2)]
    A = (h * area_c).sum()
    return A, branches


# frame function diagnostic at h = 1 from the Q1 eigenvector
h = np.ones(2 * NE)
A, br = evaluate(h)
lam_hyp = min(l for (_, l, _) in br)
print("hyperbolic branches:", [(nm, round(l, 5)) for nm, l, _ in br])
g14 = [g for (nm, l, g) in br if nm == "Q1"][0]
dens = -g14 / (lam_hyp * area_c)              # per-cell density of phi^2 (normalised: sum dens*area = 1)
Fbar = dens * A                               # mean-1 frame function
print(f"frame function: min {Fbar.min():.3f} max {Fbar.max():.3f}  "
      f"mean |Fbar-1| = {(np.abs(Fbar-1)*area_c).sum()/A:.4f}  L2 dev = {np.sqrt(((Fbar-1)**2*area_c).sum()/A):.4f}")
Lam0 = lam_hyp * A / A0
print(f"Lambda(hyp)/Lambda(hyp) = 1.0000   (lambda_1 = {lam_hyp:.5f}, gon bound 67.5*lambda = {67.5*lam_hyp:.3f})")

# max-min ascent in log h (L2-gradient w.r.t. area), steepest direction via min-norm point of active set
def direction(h, A, br):
    vals = np.array([l for (_, l, _) in br])
    lmin = vals.min()
    act = [i for i in range(len(br)) if vals[i] <= lmin * 1.01]
    Gs = np.array([((br[i][2] * h) / br[i][1] + (h * area_c) / A) / area_c for i in act])
    wts = np.ones(len(act)) / len(act)
    for _ in range(300):
        d = wts @ Gs
        j = np.argmin((Gs * area_c) @ d)
        diff = -wts; diff[j] += 1
        dd = diff @ Gs
        denom = (dd * dd * area_c).sum()
        if denom < 1e-20:
            break
        gam = np.clip(-(d * dd * area_c).sum() / denom, 0, 1)
        if gam < 1e-10:
            break
        wts = wts + gam * diff
    d = wts @ Gs
    return d, np.sqrt((d * d * area_c).sum() / A0), [f"{br[i][0]}:{vals[i]:.4f}" for i in act]

eta = np.zeros(2 * NE)
h = np.exp(eta); A, br = evaluate(h)
Lam = min(l for (_, l, _) in br) * A / A0
step = 0.3
for it in range(iters):
    d, dn, names = direction(h, A, br)
    print(f"it {it:2d}: Lambda/Lambda_hyp = {Lam/Lam0:.4f}  lambda_min = {min(l for (_,l,_) in br):.5f}"
          f"  A/A0 = {A/A0:.3f}  h in [{h.min():.3f},{h.max():.3f}]  active {names}  |grad| = {dn:.2e}", flush=True)
    if dn < 1e-7:
        break
    while step > 1e-3:
        eta2 = eta + step * d / dn
        h2 = np.exp(eta2); A2, br2 = evaluate(h2)
        Lam2 = min(l for (_, l, _) in br2) * A2 / A0
        if Lam2 > Lam:
            eta, h, A, br, Lam = eta2, h2, A2, br2, Lam2
            step *= 1.3
            break
        step *= 0.5
    else:
        break
print(f"final Lambda/Lambda_hyp = {Lam/Lam0:.4f}; implied gonality bound on the numerical scale "
      f"{67.5*lam_hyp*Lam/Lam0:.3f}  (needed > 24; certified-baseline equivalent needs ratio > 1.0430)")
np.save(f"eta_n{n}_t{tri}.npy", eta)
