"""Certified lambda_1 >= bound for every faithful A7-curve of a rigid signature (p,q,r), and the Li-Yau test
lambda_1 > 48/(g-1)  (=> gon >= 25).  Same chain as certify.py / NOTES section 7, on the (p,q,r) triangle:
  lambda_1(C) = min( mu_2(Q0), mu_1(Q1), mu_1(Q2) )   (cover.py: Q0, Q1, Q2 see every irreducible of A7),
each mu certified by the CR lower bound on the ball-arithmetic comparison problem + verified Cholesky.
Curves: signatures_spectrum.curves (one per S7-conjugation / mirror orbit; both operations are isometries, so
the spectrum, hence the certificate, is the same on the whole orbit).

Usage: python3 certify_signatures.py n p q r [p q r ...]   (n for Q1/Q2; Q0 uses n0 = 16)
"""
import os
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import certify as cf
from certify import certify_tiles, certify_trivial, QUOTIENTS
from orbifold import quotient_tiles
from signatures_spectrum import curves
from a7 import cycle_type

_cache = {}
_element_data = cf.element_data


def element_data_cached(n, XA, XB, XC, sub=3):
    key = (n, sub, str(XA[0]), str(XB[0]), str(XC[0]))
    if key not in _cache:
        _cache[key] = _element_data(n, XA, XB, XC, sub)
    return _cache[key]


cf.element_data = element_data_cached

if __name__ == "__main__":
    n = int(sys.argv[1])
    sigs = [tuple(map(int, sys.argv[i:i + 3])) for i in range(2, len(sys.argv), 3)]
    for (p, q, r) in sigs:
        g = 1 + 1260 * (1 - 1 / p - 1 / q - 1 / r)
        thr = 48 / (g - 1)
        t0 = time.time()
        cs = curves(p, q, r)
        print(f"({p},{q},{r}) genus {g:.0f}, threshold 48/(g-1) = {thr:.5f}, {len(cs)} curve(s) "
              f"[{sum(o for *_, o in cs) // 2520} A7-classes]", flush=True)
        for (a, b, osz) in cs:
            out = {}
            r0 = certify_trivial(16, None, sigma=1.0, verbose=False, ab=(a, b), pqr=(p, q, r))
            out["Q0"] = r0["bound"] if r0["ok"] else None
            for which in ("Q1", "Q2"):
                reps, glue, _ = quotient_tiles(a, b, *QUOTIENTS[which])
                res = certify_tiles(reps, glue, n, pqr=(p, q, r), sigma_factor=0.995)
                out[which] = res["bound"] if res["ok"] else None
                out[which + "_h"] = round(res["lam_h"][0], 5)
            ok = all(out[k] is not None for k in ("Q0", "Q1", "Q2"))
            lam = min(out["Q0"], out["Q1"], out["Q2"]) if ok else None
            verdict = ("CERTIFIED gon >= 25" if ok and lam > thr else "not certified")
            print(f"  a={cycle_type(a)} b={cycle_type(b)} orbit {osz}: Q0 >= {out['Q0']}, Q1 >= {out['Q1']} "
                  f"(h: {out['Q1_h']}), Q2 >= {out['Q2']} (h: {out['Q2_h']})  =>  lambda_1 >= {lam}, "
                  f"Li-Yau gon >= {lam * (g - 1) / 2 if lam else None}  [{verdict}]  ({time.time()-t0:.0f}s)",
                  flush=True)
