"""Run the full spectral certificate and write certificate.txt.

For every A7-conjugacy class of (2,4,7) generating triples (representatives 0, 1, 12, 14 of
triples_data.py; see classreps.py) certify lower bounds for the lowest eigenvalue of
    Q0 (A7-invariant functions, first NONZERO eigenvalue),  Q1,  Q2   (sign-twisted quotients),
which together see every irreducible representation of A7 (cover.py).  Hence
    lambda_1(C) >= min(Q0, Q1, Q2),
and by the Hersch / Yang-Yau / Li-Yau inequality lambda_1 * Area <= 8 pi deg, Area(C) = 540 pi:
    gon(C) >= 67.5 lambda_1(C),     gon(C/<tau>) >= 33.75 lambda_1(C).
Takes a few minutes.  Usage: python3 run_certificate.py [n_Q12=32] [n_Q0=16]
"""
import sys, time, platform
import numpy as np
from flint import arb
import flint, scipy, cvxopt
import cover
from certify import certify, fdown

REPS = [0, 1, 12, 14]


def main(n12=32, n0=16, out="certificate.txt"):
    cover.check_table()
    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)
    log("Spectral certificate for the (2,4,7) A7-curves   (NOTES.md section 7)")
    log(f"python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, "
        f"python-flint {flint.__version__}, cvxopt {cvxopt.__version__}")
    log(f"mesh: n = {n12} for Q1, Q2 and n = {n0} for Q0 (each hyperbolic triangle cut into n^2 pieces)")
    log("")
    worst = None
    for tri in REPS:
        best = {}
        for q, n, sf in (("Q0", n0, None), ("Q1", n12, 0.999), ("Q2", n12, 0.999)):
            t0 = time.time()
            r = certify(q, n, tri, sf, verbose=False) if sf else certify(q, n, tri, verbose=False)
            assert r["ok"], r
            best[q] = r["bound"]
            log(f"triple {tri:2d}  {q}: dofs {r['dofs']:7d}  lambda_h {r['lam_h'][1 if q == 'Q0' else 0]:.6f}  "
                f"sigma {r['sigma']:.6f}  C_h^2 {r['Ch2']:.3e}  Cholesky shift {r['chol_shift']:.2e} "
                f"(needed {r['chol_need']:.2e})  ==> certified >= {r['bound']:.6f}   [{time.time() - t0:.0f}s]")
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
    n12 = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    n0 = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    main(n12, n0)
