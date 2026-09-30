"""Symbolic checks of the coordinate-intersection bound and its small-degree ranges."""
import sympy as S
from fractions import Fraction as F
from weighted_budget_checks import pi,h19_genus,A


def main():
    a,b,c,m,x,y,z = S.symbols('a b c m x y z',positive=True)
    M = S.Matrix([[0,c,b],[c,0,a],[b,a,0]])
    numerator = S.Matrix([[-a*a,a*b,a*c],[a*b,-b*b,b*c],[a*c,b*c,-c*c]])
    assert S.simplify(M.det()-2*a*b*c) == 0
    assert M*numerator == S.eye(3)*(2*a*b*c)
    v = S.Matrix([x,y,z])
    quadratic = (v.T*numerator*v)[0]/(2*a*b*c)
    balanced = S.factor(quadratic.subs({x:m,y:m,z:m}))
    expected = m*m*(2*(a*b+a*c+b*c)-a*a-b*b-c*c)/(2*a*b*c)
    assert S.simplify(balanced-expected) == 0
    genus = 1+balanced/2+(a+b+c-6)*m/2
    wb = (m-1)**2-(a+b-1)**2*m*m/(4*a*b)+(a+b-1)*m/2
    assert S.simplify(genus.subs(c,1)-wb) == 0
    print('determinant and inverse identities passed')
    print('general balanced genus bound:',S.factor(genus))
    print('graph weighted budget equals the intersection bound')
    for abc in [(2,1,1),(2,2,1),(2,2,2)]:
        value = S.expand(genus.subs(dict(zip((a,b,c),abc))))
        print(abc,value)
    assert pi(15,7,1) == 9
    assert pi(18,8,0) == 13
    assert pi(21,8,2) == 17
    assert h19_genus(7) == 16
    assert int(1+F(3*6**2,8)) == 14
    assert int(1+F(3*7**2,8)) == 19
    assert A(8)==25
    print('small independent-triple bounds: m=5,6,7,8 -> 9,14,19,25')


if __name__ == '__main__':
    main()
