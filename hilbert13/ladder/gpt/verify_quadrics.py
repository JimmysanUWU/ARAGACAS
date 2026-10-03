"""Exact specialization certificates for both exceptional quadratic systems.

For the exceptional six V of 3.A7, Sym^2(V) = V^* + W15.  The
degree-7 multiples of V^*, and the degree-4 multiples of W15, fill
the entire homogeneous piece.  Hence both systems have empty base locus
in characteristic zero.  Run with Python, numpy and sympy installed.
"""
from itertools import combinations_with_replacement
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from a7 import A7, E, cyc, generated, inv, mul, order
from triples_data import triples
from audit_exact import atlas_generators, lifts


def row_basis(matrix, prime):
    rows = np.array(matrix, dtype=np.int64) % prime
    rank = 0
    for column in range(rows.shape[1]):
        candidates = np.flatnonzero(rows[rank:, column])
        if not len(candidates):
            continue
        hit = rank + int(candidates[0])
        rows[[hit, rank]] = rows[[rank, hit]]
        rows[rank, column:] = (rows[rank, column:] * pow(int(rows[rank, column]), -1, prime)) % prime
        if rank + 1 < len(rows):
            rows[rank + 1:, column:] = (
                rows[rank + 1:, column:]
                - rows[rank + 1:, column, None] * rows[rank, column:]
            ) % prime
        rank += 1
        if rank == len(rows):
            break
    return rows[:rank]


def macaulay_rank(quadrics, degree, pairs, prime):
    monomials = list(combinations_with_replacement(range(6), degree))
    index = {monomial: i for i, monomial in enumerate(monomials)}
    factors = list(combinations_with_replacement(range(6), degree - 2))
    relations = np.zeros((len(quadrics) * len(factors), len(monomials)), dtype=np.int64)
    for j, quadric in enumerate(quadrics):
        for k, factor in enumerate(factors):
            for coefficient, pair in zip(quadric, pairs):
                relations[j * len(factors) + k, index[tuple(sorted(pair + factor))]] = coefficient
    return len(row_basis(relations, prime)), len(monomials)


def coset_orbit_sizes(inertia, generators):
    subgroup = generated([inertia])
    representatives, index = [], {}
    for g in A7:
        if g in index:
            continue
        position = len(representatives)
        representatives.append(g)
        for h in subgroup:
            index[mul(g, h)] = position
    actions = [[index[mul(h, g)] for g in representatives] for h in generators]
    unseen = set(range(len(representatives)))
    sizes = []
    while unseen:
        root = unseen.pop()
        queue = [root]
        for i in queue:
            for action in actions:
                j = action[i]
                if j in unseen:
                    unseen.remove(j)
                    queue.append(j)
        sizes.append(len(queue))
    return sorted(sizes)


def main():
    prime = 43
    pairs = list(combinations_with_replacement(range(6), 2))
    (Q, X, Y), _ = atlas_generators()
    data = lifts(Q, X, Y)
    c = cyc(tuple(range(7)))

    def same_order_lift(g):
        for scalar in ((1, 0), (0, 1), (-1, -1)):
            M = tuple(tuple(Q.scalar_multiply(scalar, entry) for entry in row) for row in data[g])
            if Q.power(M, order(g)) == Q.identity:
                return M
        raise AssertionError("No same-order lift")

    # The two Klein subgroups containing this Sylow 7.  Canonical lifts
    # of their order-2/order-7 generators split the central extension;
    # every subgroup Cayley edge is verified, not assumed to split.
    klein = {}
    for s in A7:
        if order(s) != 2 or order(mul(s, c)) != 3:
            continue
        H = frozenset(generated([s, c]))
        if len(H) == 168 and H not in klein:
            klein[H] = s
    assert len(klein) == 2
    for index, (H, s) in enumerate(klein.items()):
        subgroup_lifts = {E: Q.identity}
        queue = [E]
        generators = [(g, same_order_lift(g)) for g in (s, c)]
        for g in queue:
            for h, M in generators:
                gh = mul(g, h)
                actual = Q.multiply(subgroup_lifts[g], M)
                if gh not in subgroup_lifts:
                    subgroup_lifts[gh] = actual
                    queue.append(gh)
                else:
                    assert subgroup_lifts[gh] == actual
        assert set(subgroup_lifts) == set(H)
        norm = (0, 0)
        sym_trace_sum = (0, 0)
        trace_sum = (0, 0)
        for M in subgroup_lifts.values():
            trace = Q.trace(M)
            conjugate = (trace[0] - trace[1], -trace[1])
            term = Q.scalar_multiply(trace, conjugate)
            norm = tuple(a + b for a, b in zip(norm, term))
            square = Q.scalar_multiply(trace, trace)
            power_trace = Q.trace(Q.power(M, 2))
            sym_trace_sum = tuple(a + b + d for a, b, d in zip(sym_trace_sum, square, power_trace))
            trace_sum = tuple(a + b for a, b in zip(trace_sum, trace))
        assert norm == (168, 0) and trace_sum == (0, 0)
        assert sym_trace_sum == (336, 0)
        print(f"Klein subgroup {index}: genuine irreducible six; invariant quadrics = 1", flush=True)
        a, b, _ = triples[0]
        orbit_sizes = [coset_orbit_sizes(inertia, (s, c)) for inertia in (a, b, c)]
        assert orbit_sizes == [[84] * 3 + [168] * 6, [42, 84] + [168] * 3, [24, 168, 168]]
        print(f"  branch-fibre orbit sizes: {orbit_sizes}; unique orbit of size 24", flush=True)
    solutions = [(a, b, c, d) for a in range(6) for b in range(3)
                 for c in range(2) for d in range(1)
                 if 24 * a + 42 * b + 84 * c + 168 * d == 120]
    assert solutions == [(5, 0, 0, 0)]
    print("Klein degree-120 divisor: unique orbit combination = five times the 24-orbit", flush=True)
    for omega in (6, 36):
        assert (omega * omega + omega + 1) % prime == 0
        reduced = {g: np.array([[(a + b * omega) % prime for a, b in row] for row in A], dtype=np.int64)
                   for g, A in data.items()}
        for dual in (False, True):
            matrices = ({g: reduced[inv(g)].T for g in A7} if dual else reduced)
            projector = np.zeros((21, 21), dtype=np.int64)
            for A in matrices.values():
                symmetric_square = np.zeros((21, 21), dtype=np.int64)
                for column, (j, k) in enumerate(pairs):
                    product = np.outer(A[:, j], A[:, k])
                    for row, (a, b) in enumerate(pairs):
                        symmetric_square[row, column] = product[a, b] + (product[b, a] if a != b else 0)
                # conjugate chi(V*) = chi(V); the product with Sym^2(A)
                # is independent of the chosen central lift.
                projector = (projector + int(np.trace(A)) * symmetric_square) % prime
            projector = projector * (6 * pow(2520, -1, prime) % prime) % prime
            assert np.array_equal(projector @ projector % prime, projector)
            six = row_basis(projector.T, prime)
            fifteen = row_basis((np.eye(21, dtype=np.int64) - projector).T, prime)
            assert len(six) == 6 and len(fifteen) == 15
            print(f"p=43, omega={omega}, dual={dual}: projector ranks 6 + 15", flush=True)
            for label, quadrics, degree in (("six", six, 7), ("fifteen", fifteen, 3), ("fifteen", fifteen, 4)):
                rank, dimension = macaulay_rank(quadrics, degree, pairs, prime)
                assert rank == (dimension - 1 if degree == 3 else dimension)
                print(f"  {label} quadratic ideal in degree {degree}: rank {rank}/{dimension}", flush=True)
    print("PASS: both quadratic systems have empty projective base locus; the fifteen ideal has Hilbert function (1,6,6,1)")


if __name__ == "__main__":
    main()
