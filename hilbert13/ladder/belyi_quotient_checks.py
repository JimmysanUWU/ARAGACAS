"""Exact Belyi identities and separation of the two quotient resolvents.

This does not compute an expanded elliptic model or decide either torsion class.
"""
from fractions import Fraction as F
from itertools import combinations
import sympy as S
from flint import nmod_poly


def perfect_matchings(labels):
    if not labels:
        yield ()
        return
    a = labels[0]
    for b in labels[1:]:
        rest = [x for x in labels if x not in [a,b]]
        for matching in perfect_matchings(rest):
            yield ((a,b),)+matching


def mod_fraction(x,p):
    x = F(x)
    return x.numerator * pow(x.denominator,-1,p) % p


def beta_coefficients(p,s):
    t = (2*s*pow(3,-1,p)-10*pow(21,-1,p))%p
    r = (-128*(51*s-65)*pow(1750329,-1,p))%p
    q = [0,0,0,0,-t*pow(4,-1,p),(s+t)*pow(5,-1,p),-(s+1)*pow(6,-1,p),pow(7,-1,p)]
    beta = [(-x*pow(r,-1,p))%p for x in q]
    beta[0] = 1
    return beta


def evaluate(coeffs,z,p):
    result = 0
    for a in reversed(coeffs):
        result = (result*z+a)%p
    return result


def main():
    s,z,X = S.symbols('s z X')
    relation = 21*s*s-42*s+5
    t = 2*s/3-S.Rational(10,21)
    q = z**7/7-(s+1)*z**6/6+(s+t)*z**5/5-t*z**4/4
    r = -S.Rational(128,1750329)*(51*s-65)
    assert S.expand(S.diff(q,z)-z**3*(z-1)*(z*z-s*z+t)) == 0
    assert q.subs(z,0) == q.subs(z,1) == 0
    rem = S.Poly(S.rem(q,z*z-s*z+t,z),z)
    coefficient = S.factor(rem.coeff_monomial(z))
    factors = (7*s*s-7*s+4)*(21*s*s-56*s+40)*relation
    assert S.simplify(coefficient/factors) == -S.Rational(1,129654)
    assert S.rem(rem.coeff_monomial(z),relation,s) == 0
    assert S.rem(rem.coeff_monomial(1)-r,relation,s) == 0
    assert S.rem(s*s-4*t,relation,s) != 0
    a7 = -1/(7*r)
    assert S.simplify(-1/a7*(-1/r)**7+7/r**6) == 0
    lam = 21*s
    integral = X**7-(lam+21)*X**6+36*(lam-6)*X**5-324*(lam-15)*X**4
    assert S.expand(integral.subs(X,18*z)-7*18**7*q) == 0
    assert S.expand(lam*lam-42*lam+105) == 21*relation
    print("critical-value parameter:",coefficient)
    print("critical values: q(0)=q(1)=0; quadratic critical points have q=r")
    print("discriminant: -7(-1/r)^6 T^2(T-1)^4")
    print("integral polynomial scaling and lambda relation passed")
    p = 17
    coeff = beta_coefficients(p,3)
    coeff[0] = (coeff[0]-2)%p
    poly = nmod_poly(coeff,p)
    factors_mod = poly.factor()[1]
    degrees = sorted(int(f.degree()) for f,e in factors_mod for _ in range(e))
    assert degrees == [2,5] and poly.gcd(poly.derivative()).degree() == 0
    print("unramified mod-17 specialization factor degrees:",degrees)
    p,s0,z0 = 457,118,29
    coeff = beta_coefficients(p,s0)
    value = evaluate(coeff,z0,p)
    assert value == 173
    all_roots = [a for a in range(p) if evaluate(coeff,a,p) == value]
    assert all_roots == [29,100,102,274,302,389,390]
    roots = [a for a in all_roots if a != z0]
    matchings = list(perfect_matchings(list(range(6))))
    pentads = [indices for indices in combinations(range(15),5)
               if len({edge for i in indices for edge in matchings[i]}) == 15]
    assert len(matchings) == 15 and len(pentads) == 6
    sigmas = [sum(roots[i]*roots[j] for i,j in M)%p for M in matchings]
    for power in [1,2]:
        assert len({sum(pow(sigmas[i],power,p) for i in Pi)%p for Pi in pentads}) == 1
    elliptic_values = [sum(pow(sigmas[i],3,p) for i in Pi)%p for Pi in pentads]
    assert sorted(elliptic_values) == sorted([128,123,375,88,438,302])
    partitions = [A for A in combinations(range(6),3) if 0 in A]
    def prod_mod(values):
        result = 1
        for a in values:
            result = result*a%p
        return result
    genus_two_values = [(prod_mod(roots[i] for i in A)+prod_mod(roots[i] for i in range(6) if i not in A))%p
                        for A in partitions]
    assert sorted(genus_two_values) == sorted([130,40,334,449,301,390,185,428,395,391])
    assert len(set(elliptic_values)) == 6 and len(set(genus_two_values)) == 10
    print("elliptic pentad values:",elliptic_values)
    print("genus-two partition values:",genus_two_values)
    print("separation witnesses passed; torsion remains open")


if __name__ == "__main__":
    main()
