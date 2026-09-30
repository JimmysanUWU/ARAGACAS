"""Meaningful exact fixtures for the independent factor-residual checker."""
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp
from residual_certificate import (bessel_zero_lower_proof, export_factor,
                                  residual_enclosure, q, SHIFT)


def exact_residual_inf(L, A):
    L = L.toarray()
    A = A.toarray()
    n = len(L)
    return max(sum(abs(sum(q(L[i, k]) * q(L[j, k]) for k in range(n)) - q(A[i, j]))
                   for j in range(n)) for i in range(n))


def main():
    print("Bessel Bernstein minimum:", bessel_zero_lower_proof())
    # A star forces a nonidentity sparse ordering, so a transposed/reversed
    # permutation convention cannot silently pass the factor comparison.
    n = 16
    A = np.eye(n) * 4
    A[0, 1:] = A[1:, 0] = -0.125
    A = sp.csc_matrix(A)
    L, p, pinv = export_factor(A)
    assert not np.array_equal(p, np.arange(n))
    assert np.array_equal(p[pinv], np.arange(n))
    result = residual_enclosure(L, A, p, block_rows=3, progress=False)
    exact = exact_residual_inf(L, A[p, :][:, p])
    assert exact <= result["delta"]
    assert result["delta"] < q(SHIFT)
    print("permutation:", p.tolist())
    print("exact residual / enclosure:", float(exact), float(result["delta"]))
    bad = L.copy()
    bad[0, 0] += 0.5
    rejected = residual_enclosure(bad, A, p, block_rows=3, progress=False)
    assert rejected["delta"] > q(SHIFT)
    print("corrupted factor rejected:", float(rejected["delta"]))
    # Dense constant deflation produces a different sparsity pattern.
    A0 = sp.csc_matrix([[1.0, -0.5], [-0.5, 1.0]])
    L0, p0, _ = export_factor(A0)
    result0 = residual_enclosure(L0, A0, p0, progress=False)
    assert exact_residual_inf(L0, A0[p0, :][:, p0]) <= result0["delta"] < q(SHIFT)
    print("dense deflation fixture passed")


if __name__ == "__main__":
    main()
