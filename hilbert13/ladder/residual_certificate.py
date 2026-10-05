"""Independent factor-residual certificate for the existing A7 CR comparison problem.

Run all twelve cases, sequentially (about 7 GB peak):
    OPENBLAS_NUM_THREADS=1 python3 residual_certificate.py --out residual_results
Use --which Q0 Q2 for the eight smaller cases, or --triples 12 14 for a subset.
--quick uses n=32 and sigma=.3 for Q1; it does NOT certify the full .34089 bound.

No backward-error theorem for CHOLMOD is assumed. We export its factor and
permutation, compute the residual, and bound every subsequent reduction using
exact binary64 Fractions. JSON records and permutation NPZ files are saved after
each completed case. The geometric/interpolation inputs are those of certify.py;
the coarse Bessel-zero input is proved here by rational Bernstein coefficients.
"""
import argparse
import gc
import hashlib
import json
import math
import platform
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy
import scipy.sparse as sp
import cvxopt
import cvxopt.cholmod as chol
import flint
from flint import arb
import certify as source
from orbifold import quotient_tiles
from triples_data import triples

U = F(1, 2**53)
# Far above the total possible subnormal error for the operations in these cases.
TINY = F(1, 10**280)
SHIFT = 1e-9
SAFE_CH2 = {16: F(1, 400), 32: F(3, 5000), 96: F(3, 50000)}


def q(x):
    x = float(x)
    if not math.isfinite(x):
        raise ArithmeticError("nonfinite floating-point value")
    return F.from_float(x)


def gamma(k):
    ku = int(k) * U
    if ku >= 1:
        raise ArithmeticError("reduction too long for gamma bound")
    return ku / (1 - ku)


def positive_bound(value, operations):
    """Upper bound for a nonnegative exact reduction with this rounded result."""
    return q(value) / (1 - gamma(operations)) + TINY


def matrix_positive_norm(A):
    """Bound both 1- and infinity-norms of the exact nonnegative entries of A."""
    n = max(A.shape)
    if sp.issparse(A):
        A = abs(A)
        row = np.asarray(A.sum(axis=1)).ravel().max(initial=0.0)
        col = np.asarray(A.sum(axis=0)).ravel().max(initial=0.0)
    else:
        A = np.abs(A)
        row = A.sum(axis=1).max(initial=0.0)
        col = A.sum(axis=0).max(initial=0.0)
    return positive_bound(max(row, col), 2 * n + 2)


def bessel_zero_lower_proof():
    """J1(x)/(x/2)>0 for 0<x<=19/5, hence j_(1,1)>19/5.

    With t=x^2/4<=361/100, the series truncated after k=5 is a lower
    bound: the omitted alternating tail starts positive and decreases.
    Every Bernstein coefficient of this degree-five polynomial is positive.
    """
    T = F(361, 100)
    c = [F((-1)**k, math.factorial(k) * math.factorial(k + 1)) for k in range(6)]
    b = [sum(c[k] * T**k * F(math.comb(i, k), math.comb(5, k))
             for k in range(i + 1)) for i in range(6)]
    assert min(b) > 0
    return min(b)


def fraction_record(x):
    return {"numerator": x.numerator, "denominator": x.denominator,
            "display": float(x)}


def element_inputs(n):
    XA, XB, XC = source.ref_triangle_arb()
    ed = source.element_data(n, XA, XB, XC)
    # Bound each Arb ratio BEFORE taking the maximum; ordering overlapping
    # balls is not used to select a rigorous maximum.
    ratios = [q(source.fup(arb(e[2]) / e[3])) for e in ed]
    kappa2 = (F(1, 8) + F(2, 1) / F(19, 5)**2) / n**2
    Ch2 = kappa2 * max(ratios)
    assert Ch2 <= SAFE_CH2[n], (n, Ch2)
    return source.elem_mats_rigorous(n, ed), Ch2


def comparison_matrix(which, n, tri, mats, sigma):
    a, b, _ = triples[tri]
    Kg, signs = ([a, b], [1, 1]) if which == "Q0" else source.QUOTIENTS[which]
    reps, glue, order = quotient_tiles(a, b, Kg, signs)
    asm = source.assemble_rigorous if which == "Q0" else source.assemble_rigorous_fast
    K, Kabs, Kerr, cnt, md, merr, mabs, mcnt = asm(n, reps, glue, mats)
    N = K.shape[0]
    nc = int(cnt.max())
    nm = int(mcnt.max())
    # CSR coalescing of the positive error and magnitude matrices also rounds.
    ek = (matrix_positive_norm(Kerr) + gamma(nc) * matrix_positive_norm(Kabs)) / (1 - gamma(nc))
    em = (q(merr.max()) + gamma(nm) * q(mabs.max())) / (1 - gamma(nm)) + TINY
    sm = sigma * md
    diag = K.diagonal() - sm
    rounding = U / (1 - U) * (q(abs(sm).max()) + q(abs(diag).max())) + TINY
    assembly_inf = ek + q(sigma) * em + rounding
    alpha = None
    if which == "Q0":
        mlo = q(md.sum()) / (1 + gamma(N + 1)) - N * em - TINY
        assert mlo > 0
        alpha = 1.0
        while q(alpha) * mlo <= q(sigma):
            alpha *= 2.0
        zz = alpha * np.outer(md, md)
        B = K.toarray() - np.diag(sm) + zz
        # Md is a rounded approximation to the exact comparison mass vector.
        maxmd = q(abs(md).max())
        sum_md = positive_bound(abs(md).sum(), N + 1)
        rank_error = q(alpha) * em * (N * maxmd + sum_md + N * em)
        outer_round = U / (1 - U) * matrix_positive_norm(zz) + TINY
        final_add = U / (1 - U) * matrix_positive_norm(B) + TINY
        assembly_inf += rank_error + outer_round + final_add
        B = sp.csc_matrix(B)
    else:
        B = (K - sp.diags(sm)).tocsc()
    # This conversion remains valid even if sparse coalescing has produced
    # a nonsymmetric last-bit error: ||E||2 <= ||E||F <= sqrt(N)||E||inf.
    sqrt_N_upper = math.isqrt(N) + (math.isqrt(N)**2 != N)
    assembly_2 = sqrt_N_upper * assembly_inf
    del K, Kabs, Kerr, cnt, md, merr, mabs, mcnt, reps, glue
    gc.collect()
    return B, assembly_inf, assembly_2, order, alpha


def export_factor(target):
    """Recover P before getfactor consumes CHOLMOD's factor object."""
    N = target.shape[0]
    coo = target.tocoo()
    S = cvxopt.spmatrix(coo.data.tolist(), coo.row.tolist(), coo.col.tolist(), (N, N))
    chol.options["supernodal"] = 2
    fac = chol.symbolic(S)
    chol.numeric(S, fac)
    del S, coo
    pvec = cvxopt.matrix(np.arange(N, dtype=float))
    chol.solve(fac, pvec, sys=7)  # P' X = b, so X = P b
    p = np.asarray(pvec).ravel().astype(np.int64)
    pinv_vec = cvxopt.matrix(np.arange(N, dtype=float))
    chol.solve(fac, pinv_vec, sys=8)
    pinv = np.asarray(pinv_vec).ravel().astype(np.int64)
    assert np.array_equal(np.sort(p), np.arange(N))
    assert np.array_equal(p[pinv], np.arange(N))
    exported = chol.getfactor(fac)
    L = sp.csc_matrix((np.asarray(exported.V).ravel(),
                       np.asarray(exported.I).ravel(),
                       np.asarray(exported.CCS[0]).ravel()), shape=(N, N))
    del exported, fac, pvec, pinv_vec
    gc.collect()
    assert np.all(np.isfinite(L.data))
    assert np.all(L.diagonal() > 0)
    assert sp.triu(L, 1).nnz == 0
    return L, p, pinv


def residual_enclosure(L, target, p, block_rows=2048, progress=True):
    N = L.shape[0]
    Lr = L.tocsr()
    k = int(np.diff(Lr.indptr).max())
    l1 = positive_bound(np.asarray(abs(L).sum(axis=0)).max(), 2 * N + 2)
    li = positive_bound(np.asarray(abs(Lr).sum(axis=1)).max(), 2 * N + 2)
    BP = target[p, :][:, p].tocsr()
    LT = Lr.T.tocsc()
    col_sums = np.zeros(N)
    row_max = 0.0
    for start in range(0, N, block_rows):
        end = min(N, start + block_rows)
        R = Lr[start:end] @ LT - BP[start:end]
        row_max = max(row_max, float(np.asarray(abs(R).sum(axis=1)).max(initial=0.0)))
        col_sums += np.asarray(abs(R).sum(axis=0)).ravel()
        if progress and (end == N or (start // block_rows) % 100 == 0):
            print(f"  residual rows {end}/{N}", flush=True)
        del R
    rho = positive_bound(max(row_max, float(col_sums.max())), 2 * N + 2)
    # At most k products enter any sparse dot product. The subtraction of
    # target is enclosed by rho/(1-u); no assertion about factor accuracy.
    delta = rho / (1 - U) + gamma(k + 1) * l1 * li + TINY
    return {"rowmax": k, "L1": l1, "Linf": li, "rho": rho, "delta": delta}


def run_case(which, n, tri, inputs, out, block_rows):
    t0 = time.time()
    sigma = {"Q0": 1.0, "Q1": 0.3409 if n == 96 else 0.3, "Q2": 0.68}[which]
    safe_sigma = {"Q0": F(1), "Q1": F(340899, 10**6) if n == 96 else F(299999, 10**6),
                  "Q2": F(17, 25)}[which]
    assert safe_sigma <= q(sigma)
    mats, actual_ch2_upper = inputs
    print(f"BEGIN triple={tri} {which} n={n}", flush=True)
    B, eta_inf, eta, order, alpha = comparison_matrix(which, n, tri, mats, sigma)
    N = B.shape[0]
    target = (B - SHIFT * sp.eye(N, format="csc")).tocsc()
    shift_error = U / (1 - U) * q(abs(target.diagonal()).max()) + TINY
    print(f"  assembly done: {N} dofs", flush=True)
    L, p, pinv = export_factor(target)
    print(f"  factor exported: {L.nnz} stored entries", flush=True)
    enclosed = residual_enclosure(L, target, p, block_rows)
    needed = enclosed["delta"] + eta + shift_error
    assert q(SHIFT) > needed, (which, tri, float(needed))
    bound = safe_sigma / (1 + SAFE_CH2[n] * safe_sigma)
    perm_file = out / f"triple_{tri}_{which}_permutation.npz"
    np.savez_compressed(perm_file, permutation=p, inverse_permutation=pinv)
    canonical_bytes = p.astype("<i8", copy=False).tobytes()
    record = {"status": "PASS", "which": which, "n": n, "triple": tri,
              "dofs": N, "group_order": order, "factor_stored_entries": L.nnz,
              "rowmax": enclosed.pop("rowmax"), "alpha": alpha,
              "permutation_file": perm_file.name,
              "permutation_sha256_int64_le": hashlib.sha256(canonical_bytes).hexdigest(),
              "seconds": time.time() - t0,
              "versions": {"python": platform.python_version(), "numpy": np.__version__,
                           "scipy": scipy.__version__, "flint": flint.__version__, "cvxopt": cvxopt.__version__}}
    rationals = dict(enclosed, sigma=q(sigma), safe_sigma=safe_sigma, shift=q(SHIFT),
                     assembly_inf=eta_inf, assembly_2=eta, shift_rounding=shift_error,
                     needed=needed, computed_Ch2_upper=actual_ch2_upper,
                     safe_Ch2=SAFE_CH2[n], spectral_bound=bound,
                     bessel_Bernstein_min=bessel_zero_lower_proof())
    record.update({key: fraction_record(value) for key, value in rationals.items()})
    json_file = out / f"triple_{tri}_{which}.json"
    json_file.write_text(json.dumps(record, indent=2) + "\n")
    print(f"PASS triple={tri} {which}: rho={float(enclosed['rho']):.6e} "
          f"delta={float(enclosed['delta']):.6e} needed={float(needed):.6e} "
          f"bound={float(bound):.10f} [{record['seconds']:.1f}s]", flush=True)
    del L, B, target, p, pinv
    gc.collect()
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--which", nargs="+", choices=["Q0", "Q1", "Q2"], default=["Q0", "Q2", "Q1"])
    parser.add_argument("--triples", nargs="+", type=int, default=[0, 1, 12, 14])
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--block-rows", type=int, default=2048)
    parser.add_argument("--out", type=Path, default=Path("residual_results"))
    args = parser.parse_args()
    assert args.block_rows > 0
    args.out.mkdir(parents=True, exist_ok=True)
    bessel_zero_lower_proof()
    for which in args.which:
        n = {"Q0": 16, "Q1": 32 if args.quick else 96, "Q2": 32}[which]
        print(f"Preparing exact coefficient enclosures for {which} n={n}", flush=True)
        inputs = element_inputs(n)
        for tri in args.triples:
            run_case(which, n, tri, inputs, args.out, args.block_rows)
        del inputs
        gc.collect()


if __name__ == "__main__":
    main()
