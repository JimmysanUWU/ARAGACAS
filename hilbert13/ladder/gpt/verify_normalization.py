"""Exact checks for the equivariant normalization-defect argument.

Run from this folder or its parent with Python and sympy installed.
No floating-point arithmetic or numerical eigenvectors are used.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from a7 import A7, E, cyc, generated, inv
from audit_exact import atlas_generators, lifts, character_polynomial


def castelnuovo(degree, dimension):
    quotient, remainder = divmod(degree - 1, dimension - 1)
    return (dimension - 1) * quotient * (quotient - 1) // 2 + quotient * remainder


def multiply(A, B, prime):
    return [[sum(A[i][k] * B[k][j] for k in range(6)) % prime
             for j in range(6)] for i in range(6)]


def column(A, j):
    return [row[j] for row in A]


def line_fixed(A, v, prime):
    image = [sum(a * b for a, b in zip(row, v)) % prime for row in A]
    return all((v[i] * image[j] - v[j] * image[i]) % prime == 0
               for i in range(6) for j in range(i + 1, 6))


def primitive_root(prime):
    return next(x for x in range(2, prime)
                if len({pow(x, k, prime) for k in range(1, prime)}) == prime - 1)


def main():
    bound = castelnuovo(60, 5)
    assert bound == 406
    budget = bound - 136
    assert budget == 270
    orbits = [2520 // e for e in (1, 2, 4, 7)]
    assert min(orbits) == 360 > budget
    assert min(orbits[:-1]) // 2 == 315 > budget
    print(f"Castelnuovo pi(60,5) = {bound}; normalization defect <= {budget}")
    print(f"Orbit sizes (inertia 1,2,4,7): {orbits}")
    print("Any nonimmersed orbit costs >= 360; any collision orbit outside D7 costs >= 315")

    (Q, X, Y), _ = atlas_generators()
    data = lifts(Q, X, Y)  # Includes every exact Cayley-edge check.
    c = cyc(tuple(range(7)))
    subgroup = generated([c])
    matrix = data[c]
    scalar = Q.central_scalar(Q.power(matrix, 7))
    # 7 == 1 (mod 3), so inverse scalar makes the lift have exact order 7.
    inverse_scalar = Q.scalar_multiply(scalar, scalar)
    matrix = tuple(tuple(Q.scalar_multiply(inverse_scalar, entry) for entry in row)
                   for row in matrix)
    assert Q.power(matrix, 7) == Q.identity
    elementary = character_polynomial(Q, matrix, None, None)
    assert elementary == [((-1) ** j, 0) for j in range(1, 7)]
    print("Canonical order-7 lift: characteristic polynomial = Phi_7")

    # Specialize Z[omega,zeta_7] at primes above p=43.  The six spectral
    # projectors remain rank one because p does not divide 7 and zeta_7
    # remains primitive.  If a characteristic-zero projective stabilizer
    # fixes a line, its two-by-two minors also vanish after specialization.
    # Thus a finite-field stabilizer is an UPPER bound, not a heuristic.
    prime = 43
    primitive = primitive_root(prime)
    zeta = pow(primitive, 6, prime)
    omega = pow(primitive, 14, prime)
    assert pow(zeta, 7, prime) == 1 and zeta != 1
    assert (omega * omega + omega + 1) % prime == 0
    identity = [[int(i == j) for j in range(6)] for i in range(6)]

    for root in (omega, pow(omega, 2, prime)):
        reduction = lambda M: [[(a + b * root) % prime for a, b in row] for row in M]
        reduced = {g: reduction(M) for g, M in data.items()}
        for dual in (False, True):
            representation = ({g: [list(row) for row in zip(*reduced[inv(g)])]
                               for g in A7} if dual else reduced)
            S = reduction(matrix)
            if dual:
                inverse_S = identity
                for _ in range(6):
                    inverse_S = multiply(inverse_S, S, prime)
                S = [list(row) for row in zip(*inverse_S)]
            powers = [identity]
            for _ in range(6):
                powers.append(multiply(powers[-1], S, prime))
            sizes = []
            for exponent in range(1, 7):
                eigenvalue = pow(zeta, exponent, prime)
                projector = [[sum(pow(eigenvalue, (-k) % 7, prime) * powers[k][i][j]
                                  for k in range(7)) * pow(7, -1, prime) % prime
                              for j in range(6)] for i in range(6)]
                assert multiply(projector, projector, prime) == projector
                v = next(column(projector, j) for j in range(6)
                         if any(column(projector, j)))
                assert all((v[i] * projector[j][k] - v[j] * projector[i][k]) % prime == 0
                           for i in range(6) for j in range(i + 1, 6) for k in range(6))
                assert all(sum(S[i][j] * v[j] for j in range(6)) % prime == eigenvalue * v[i] % prime
                           for i in range(6))
                stabilizer = {g for g in A7 if line_fixed(representation[g], v, prime)}
                assert stabilizer == subgroup, (root, dual, exponent, len(stabilizer))
                sizes.append(len(stabilizer))
            print(f"p=43, omega={root}, dual={dual}: six eigenline stabilizers {sizes}, each exactly C7")

    print("PASS: exact hypotheses for closed immersion and the degree-42 involution pencil")


if __name__ == "__main__":
    main()
