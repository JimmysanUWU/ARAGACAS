"""Run the full spectral certificate and write certificate.txt.

For every A7-conjugacy class of (2,4,7) generating triples (representatives 0, 1, 12, 14 of
triples_data.py; see curve_checks.py) certify lower bounds for the lowest eigenvalue of
    Q0 (A7-invariant functions, first NONZERO eigenvalue),  Q1,  Q2   (sign-twisted quotients),
which together see every irreducible representation of A7 (cover.py).  Hence
    lambda_1(C) >= min(Q0, Q1, Q2),
and by the Hersch / Yang-Yau / Li-Yau inequality lambda_1 * Area <= 8 pi deg, Area(C) = 540 pi:
    gon(C) >= 67.5 lambda_1(C),     gon(C/<tau>) >= 33.75 lambda_1(C).

Default: Q1 on the mesh n = 96 with sigma = 0.3409 fixed (about 1 min and 7 GB per triple class),
certifying lambda_1 >= 0.34089 and hence gon(C) >= 24.
--quick: Q1 on n = 32 with sigma = 0.999 * (computed discrete eigenvalue): lambda_1 >= 0.33335,
gon(C) >= 23 (about 3 minutes in total, 3 GB).
"""
import sys, time, platform
import numpy as np
from flint import arb
import flint, scipy, cvxopt
import cover
from certify import certify, fdown

REPS = [0, 1, 12, 14]


def main(quick=False, out=None):
    out = out or ("certificate_quick.txt" if quick else "certificate.txt")
    cover.check_table()
    # (quotient, mesh n, sigma_factor, fixed sigma)
    plan = [("Q0", 16, None, None),
            ("Q1", 32, 0.999, None) if quick else ("Q1", 96, None, 0.3409),
            ("Q2", 32, 0.999, None)]
    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)
    log("Spectral certificate for the (2,4,7) A7-curves   (2_SPECTRAL.md)")
    log(f"python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, "
        f"python-flint {flint.__version__}, cvxopt {cvxopt.__version__}")
    log("plan: " + ", ".join(f"{q} n={n}" + (f" sigma={s}" if s else "") for q, n, _, s in plan)
        + "   (each hyperbolic triangle cut into n^2 pieces)")
    log("")
    worst = None
    for tri in REPS:
        best = {}
        for q, n, sf, sig in plan:
            t0 = time.time()
            if q == "Q0":
                r = certify(q, n, tri, verbose=False)
            elif sig is not None:
                r = certify(q, n, tri, verbose=False, sigma=sig)
            else:
                r = certify(q, n, tri, sf, verbose=False)
            assert r["ok"], r
            best[q] = r["bound"]
            lh = r["lam_h"][1 if q == "Q0" else 0]
            lhs = "   (not computed)" if np.isnan(lh) else f"{lh:.6f}"
            log(f"triple {tri:2d}  {q}: n {n:3d}  dofs {r['dofs']:8d}  lambda_h {lhs}  sigma {r['sigma']:.6f}  "
                f"C_h^2 {r['Ch2']:.3e}  Cholesky shift {r['chol_shift']:.2e} (needed {r['chol_need']:.2e})"
                f"  ==> certified >= {r['bound']:.6f}   [{time.time() - t0:.0f}s]")
        lam = min(best.values())
        log(f"triple {tri:2d}  lambda_1(C) >= {lam:.6f}")
        log("")
        worst = lam if worst is None else min(worst, lam)
    L = arb(worst)
    gC = (L * 135 / 2)
    gD = (L * 135 / 4)
    log(f"ALL TRIPLES: lambda_1(C) >= {worst:.6f}")
    log(f"  gon(C)         >= 67.5  * lambda_1 >= {fdown(gC):.4f}  ==>  gon(C) >= {int(np.ceil(fdown(gC)))}")
    log(f"  gon(C/<tau>)   >= 33.75 * lambda_1 >= {fdown(gD):.4f}  ==>  gon(C/<tau>) >= {int(np.ceil(fdown(gD)))}")
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
