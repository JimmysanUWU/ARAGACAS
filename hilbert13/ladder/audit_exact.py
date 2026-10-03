"""Independent exact checks added by the October 2026 proof audit.

ATLAS generators (accessed 3 October 2026):
https://brauer.maths.qmul.ac.uk/Atlas/alt/A7/gap0/3A7G1-Ar6B0.g
https://brauer.maths.qmul.ac.uk/Atlas/alt/A7/gap0/2A7G1-Ar4aB0.g

All representation arithmetic is over an explicitly defined quadratic
integer ring.  Cyclotomic calculations use QQ[z]/Phi_84(z), not numerical
root recognition or rounded character values.  Run: python3 audit_exact.py.
"""
from collections import Counter
from fractions import Fraction
from math import comb, isqrt
import itertools
import sympy as sp
from a7 import A7, E, mul, inv, order, cycle_type
from triples_data import triples


class QuadraticMatrices:
    """Matrices over Z[w], w^2=t*w+u; entries are integer pairs (a,b)."""
    def __init__(self, dimension, t, u):
        self.dimension, self.t, self.u = dimension, t, u
        self.zero = (0, 0)
        self.identity = self.matrix([[int(i == j) for j in range(dimension)]
                                     for i in range(dimension)])

    def matrix(self, rows):
        return tuple(tuple((x, 0) if isinstance(x, int) else x for x in row) for row in rows)

    def multiply(self, A, B):
        n, t, u = self.dimension, self.t, self.u
        C = [[(0, 0) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for k in range(n):
                a, b = A[i][k]
                if a == b == 0:
                    continue
                for j in range(n):
                    c, d = B[k][j]
                    if c == d == 0:
                        continue
                    p, q = C[i][j]
                    C[i][j] = (p+a*c+u*b*d, q+a*d+b*c+t*b*d)
        return tuple(tuple(row) for row in C)

    def power(self, A, exponent):
        B = self.identity
        for _ in range(exponent):
            B = self.multiply(B, A)
        return B

    def trace(self, A):
        return (sum(A[i][i][0] for i in range(self.dimension)),
                sum(A[i][i][1] for i in range(self.dimension)))

    def scalar_multiply(self, left, right):
        a,b = left
        c,d = right
        return (a*c+self.u*b*d, a*d+b*c+self.t*b*d)

    def central_scalar(self, A):
        scalar = A[0][0]
        assert all(A[i][j] == (scalar if i == j else self.zero)
                   for i in range(self.dimension) for j in range(self.dimension))
        return scalar


def atlas_generators():
    Q6 = QuadraticMatrices(6, -1, -1)  # w=E(3)
    w, w2 = (0, 1), (-1, -1)
    X6 = Q6.matrix([[1,0,0,0,0,0],[0,w2,0,0,0,0],[0,0,1,0,0,0],
                    [0,0,0,0,0,1],[(-1,1),0,(1,-1),(0,-1),w,1],
                    [2,0,-1,-1,0,-1]])
    Y6 = Q6.matrix([[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],
                    [0,0,0,0,1,0],[1,0,0,0,0,0],[-1,1,(0,-1),0,w,1]])
    Q4 = QuadraticMatrices(4, -1, -2)  # b=(-1+sqrt(-7))/2
    b, B = (0,1), (-1,-1)
    X4 = Q4.matrix([[0,1,0,0],[-1,-1,0,0],[0,0,0,1],[0,0,-1,-1]])
    Y4 = Q4.matrix([[0,1,0,0],[0,0,1,0],[0,B,-1,-1],[1,1,(0,-1),0]])
    return (Q6,X6,Y6), (Q4,X4,Y4)


def lifts(Q, X, Y):
    a = (1,2,0,3,4,5,6)
    b = (0,1,3,4,5,6,2)
    assert Q.power(X,3) == Q.identity
    assert Q.power(Y,5) == Q.identity
    assert Q.power(Q.multiply(X,Y),7) == Q.identity
    data = {E:Q.identity}
    queue = [E]
    for g in queue:
        for h, H in ((a,X),(b,Y)):
            gh = mul(g,h)
            if gh not in data:
                data[gh] = Q.multiply(data[g],H)
                queue.append(gh)
    assert len(data) == 2520
    scalars = set()
    # Every Cayley edge is checked, rather than random associativity tests.
    # Its discrepancy must be one of the stated central scalars.
    for g in A7:
        for h,H in ((a,X),(b,Y)):
            actual = Q.multiply(data[g],H)
            target = data[mul(g,h)]
            allowed = [(1,0),(-1,0)] if Q.dimension == 4 else [(1,0),(0,1),(-1,-1)]
            matches = []
            for scalar in allowed:
                S = Q.matrix([[scalar if i == j else 0 for j in range(Q.dimension)]
                              for i in range(Q.dimension)])
                if actual == Q.multiply(S,target):
                    matches.append(scalar)
            assert len(matches) == 1
            scalars.add(matches[0])
    assert len(scalars) == (2 if Q.dimension == 4 else 3)
    return data


def cyclotomic_field(conductor=84):
    field = sp.QQ.alg_field_from_poly(sp.Poly(sp.cyclotomic_poly(conductor)))
    z = field.dtype([1,0],field.mod.to_list(),field.dom)
    assert z**conductor == field.one
    return field,z


def character_polynomial(Q, matrix, field, root):
    """Elementary symmetric functions from exact Newton identities."""
    powers, running = [], Q.identity
    for j in range(1,Q.dimension+1):
        running = Q.multiply(running,matrix)
        powers.append(Q.trace(running))
    elementary = [(1,0)]
    for j in range(1,Q.dimension+1):
        terms = [Q.scalar_multiply(elementary[j-k],powers[k-1]) for k in range(1,j+1)]
        numerator = tuple(sum((-1)**k*terms[k][i] for k in range(j)) for i in range(2))
        assert all(value % j == 0 for value in numerator)
        elementary.append(tuple(value//j for value in numerator))
    assert elementary[-1] == (1,0)  # ATLAS matrices lie in SL
    return elementary[1:]


def molien(Q, data, field, root, central_order, maximum):
    polynomials = {}
    for matrix in data.values():
        elementary = tuple(character_polynomial(Q,matrix,field,root))
        polynomials[elementary] = polynomials.get(elementary,0)+1
    totals = [(0,0) for _ in range(maximum+1)]
    for elementary,count in polynomials.items():
        complete = [(1,0)]
        for degree in range(1,maximum+1):
            terms = [Q.scalar_multiply(elementary[j-1],complete[degree-j])
                     for j in range(1,min(Q.dimension,degree)+1)]
            complete.append(tuple(sum((-1)**k*term[i] for k,term in enumerate(terms)) for i in range(2)))
        for degree in range(0,maximum+1,central_order):
            totals[degree] = tuple(totals[degree][i]+count*complete[degree][i] for i in range(2))
    dimensions = {}
    for degree,total in enumerate(totals):
        assert total[1] == 0 and total[0] % 2520 == 0
        dimension = total[0]//2520
        assert dimension >= 0
        if dimension:
            dimensions[degree] = dimension
    return dimensions


def scalar_exponent(Q, matrix, central_order):
    scalar = Q.central_scalar(matrix)
    choices = [(1,0),(-1,0)] if central_order == 2 else [(1,0),(0,1),(-1,-1)]
    assert scalar in choices
    return choices.index(scalar)


def twisted_six_rr(Q, data, field, z):
    """Multiplicity of the ATLAS six in degree 60, for all four covers.

    Formula: (d+1-g)*dim(W)/|G| + sum_i 1/e_i sum_{j=1}^{e_i-1}
    lambda_i^j * conjugate(chi_W(s_i^j))/(1-a_i^(-j)).
    Both orientations and both nontrivial central characters are checked.
    """
    results = []
    for cls in (0,1,12,14):
        a,b,_ = triples[cls]
        generators = (a,inv(b),mul(b,a))
        orders = [order(g) for g in generators]
        matrices = [data[g] for g in generators]
        central = [scalar_exponent(Q,Q.power(S,e),3) for S,e in zip(matrices,orders)]
        product = Q.multiply(Q.multiply(matrices[0],matrices[1]),matrices[2])
        central_product = scalar_exponent(Q,product,3)
        traces = []
        for S,e in zip(matrices,orders):
            running = Q.identity
            values = []
            for j in range(1,e):
                running = Q.multiply(running,S)
                values.append(Q.trace(running))
            traces.append(values)
        for orientation in (1,-1):
            successes = []
            for eps in (1,2):
                options = [[(84*k+28*eps*c)//e for k in range(e)]
                           for e,c in zip(orders,central)]
                for datum in itertools.product(*options):
                    degree = (2520*orientation*(sum(Fraction(q,84) for q in datum)
                                                - Fraction(eps*central_product,3))) % 2520
                    assert degree.denominator == 1
                    if degree != 60:
                        continue
                    m = field.convert(Fraction(-75*6,2520))
                    for e,q,powers in zip(orders,datum,traces):
                        for j,(a0,b0) in enumerate(powers,1):
                            # chi_{W_eps} is a+b*w^eps; conjugate it.
                            character_conjugate = field.convert(a0)+field.convert(b0)*z**((-28*eps)%84)
                            m += z**((q*j)%84)*character_conjugate/(1-z**((-orientation*84*j//e)%84))/e
                    coefficients = m.to_list()
                    assert len(coefficients) <= 1
                    multiplicity = int(coefficients[0]) if coefficients else 0
                    assert m == field.convert(multiplicity)
                    if multiplicity > 0:
                        assert multiplicity == 1
                        # Order-two lifts have trace 2 and fibre eigenvalue +1
                        # at both branch types which contain fixed points.
                        for index,power in ((0,1),(1,2)):
                            S = Q.power(matrices[index],power)
                            c = scalar_exponent(Q,Q.power(S,2),3)
                            shift = (-2*c)%3
                            a0,b0 = Q.trace(S)
                            trace = (field.convert(a0)+field.convert(b0)*z**28)*z**(28*shift)
                            assert trace == field.convert(2)
                            assert z**((datum[index]*power+28*eps*shift)%84) == field.one
                        successes.append((eps,datum))
            assert len(successes) == 1, (cls,orientation,successes)
            results.append((cls,orientation,successes[0]))
            print(f'L60 class {cls}, orientation {orientation:+d}: {successes[0]}, six multiplicity 1, all 18 fibres +1')
    return results


def spin45_local(Q, data):
    """Certify the least divisors which force f14=f18=0 in degree 45."""
    for cls in (0,1,12,14):
        a,b,_ = triples[cls]
        generators = (a,inv(b),mul(b,a))
        orders = [order(g) for g in generators]
        matrices = [data[g] for g in generators]
        central = [scalar_exponent(Q,Q.power(S,e),2) for S,e in zip(matrices,orders)]
        product = Q.multiply(Q.multiply(matrices[0],matrices[1]),matrices[2])
        central_product = scalar_exponent(Q,product,2)
        for orientation in (1,-1):
            count = 0
            options = [[(168*k+84*c)//e for k in range(e)] for e,c in zip(orders,central)]
            for datum in itertools.product(*options):
                degree = (2520*orientation*(sum(Fraction(q,168) for q in datum)
                                            - Fraction(central_product,2))) % 2520
                if degree != 45:
                    continue
                count += 1
                least = []
                for power in (14,18):
                    residues = []
                    for q,e in zip(datum,orders):
                        residue = Fraction(orientation*q*power*e,168)
                        assert residue.denominator == 1
                        residues.append(int(residue)%e)
                    minimum = sum(r*2520//e for r,e in zip(residues,orders))
                    assert minimum > 45*power
                    least.append(minimum)
                assert least == [3150,3330]
            assert count == 2
    print('All degree-45 spin data, all classes/orientations: least f14/f18 divisors 3150/3330')


def ordinary_rr(field,z):
    """Ordinary equivariant Riemann--Roch without rounded character sums."""
    from cover import TABLE, CLASSES, check_table
    from a7 import a7_class,cyc,power
    check_table()
    seven_a = a7_class(cyc(tuple(range(7))))
    beta = z**12+z**24+z**48

    def character(name,g):
        ct = cycle_type(g)
        if ct != (7,):
            return field.convert(TABLE[name][CLASSES.index(ct)])
        column = 7 if g in seven_a else 8
        if name == '10':
            return beta if column == 7 else -1-beta
        if name == '10b':
            return -1-beta if column == 7 else beta
        return field.convert(TABLE[name][column])

    for cls in (0,1,12,14):
        a,b,_ = triples[cls]
        generators = (a,inv(b),mul(b,a))
        for label,coefficients,degree in [('B',(0,-1,2),90),('B+T',(1,-3,2),90),('K',(1,3,6),270)]:
            result = {}
            for name,row in TABLE.items():
                m = field.convert(Fraction((degree-135)*row[0],2520))
                for g,n in zip(generators,coefficients):
                    e = order(g)
                    for j in range(1,e):
                        # On Q(sqrt(-7)), conjugation swaps the two 10s.
                        conjugate_name = {'10':'10b','10b':'10'}.get(name,name)
                        m += z**((84*n*j//e)%84)*character(conjugate_name,power(g,j))/(1-z**((-84*j//e)%84))/e
                cs = m.to_list()
                assert len(cs) <= 1
                integer = int(cs[0]) if cs else 0
                assert m == field.convert(integer)
                if integer:
                    result[name] = integer
            if label == 'B':
                assert sorted(result.values()) == [-1,-1] and result.get('35') == -1
            if label == 'B+T':
                assert result.get('6') == -1 and result.get('14a') == result.get('14b') == result.get('21') == -1
                assert (result.get('10',0),result.get('10b',0)) in ((1,0),(0,1))
            if label == 'K':
                assert result.get('1') == -1 and result.get('15') == result.get('21') == 1 and result.get('35') == 2
                assert sorted((result.get('10',0),result.get('10b',0))) == [1,2]
            print(f'Exact ordinary RR, class {cls}, {label}: {result}')


class BiquadraticIntegers:
    """Z[w,b], w^2=-w-1, b^2=-b-2; coordinates (1,w,b,wb)."""
    zero = (0,0,0,0)

    @staticmethod
    def add(x,y):
        return tuple(a+b for a,b in zip(x,y))

    @staticmethod
    def subtract(x,y):
        return tuple(a-b for a,b in zip(x,y))

    @staticmethod
    def scale(n,x):
        return tuple(n*a for a in x)

    @staticmethod
    def conjugate(x):
        a,b,c,d = x
        return (a-b-c+d,-b+d,-c+d,d)

    @staticmethod
    def multiply(x,y):
        def qw(a,b):
            p,q = a
            r,s = b
            return (p*r-q*s,p*s+q*r-q*s)
        u0,u1,v0,v1 = x[:2],x[2:],y[:2],y[2:]
        a,b,c,d = qw(u0,v0),qw(u1,v1),qw(u0,v1),qw(u1,v0)
        return (a[0]-2*b[0],a[1]-2*b[1],c[0]+d[0]-b[0],c[1]+d[1]-b[1])


def central_three_characters(Q6,six,Q4,four):
    """Recover all seven faithful 3.A7 characters from actual tensor powers.

    A virtual character with exact norm 1 and positive dimension is an
    irreducible character.  Subtraction and tensoring known characters
    therefore certify these discoveries, without rounding eigenvectors of
    a class algebra.  Their squared dimensions sum to 2520, proving that
    the central-character sector is complete.
    """
    from cover import TABLE,CLASSES
    B = BiquadraticIntegers
    G = list(A7)
    six_vector = [(*Q6.trace(six[g]),0,0) for g in G]
    ordinary = []
    for name in ('1','6','14a','14b','15','21','35'):
        ordinary.append([(TABLE[name][7 if cycle_type(g)==(7,) else CLASSES.index(cycle_type(g))],0,0,0) for g in G])
    tens = []
    for g in G:
        a = Q4.trace(four[g])
        square = Q4.scalar_multiply(a,a)
        trace_square = Q4.trace(Q4.power(four[g],2))
        numerator = tuple(x+y for x,y in zip(square,trace_square))
        assert all(x%2 == 0 for x in numerator)
        x,y = (v//2 for v in numerator)
        tens.append((x,0,y,0))
    ordinary.extend([tens,[B.conjugate(x) for x in tens]])

    def inner(left,right):
        value = B.zero
        for x,y in zip(left,right):
            value = B.add(value,B.multiply(x,B.conjugate(y)))
        assert value[1:] == (0,0,0) and value[0]%2520 == 0,value
        return value[0]//2520

    lam,sym = [],[]
    for g,x in zip(G,six_vector):
        a = Q6.trace(six[g])
        square = Q6.scalar_multiply(a,a)
        trace_square = Q6.trace(Q6.power(six[g],2))
        for sign,out in ((-1,lam),(1,sym)):
            numerator = tuple(s+sign*t for s,t in zip(square,trace_square))
            assert all(n%2 == 0 for n in numerator)
            out.append(B.conjugate((numerator[0]//2,numerator[1]//2,0,0)))
    # Sym^2(six)^* contains six once; its complement is the second 15.
    sym_complement = [B.subtract(x,y) for x,y in zip(sym,six_vector)]
    known = [six_vector,lam,sym_complement]
    assert all(inner(x,y) == int(i==j) for i,x in enumerate(known) for j,y in enumerate(known))
    products = [[B.multiply(x,y) for x,y in zip(six_vector,ch)] for ch in ordinary]
    while len(known) < 7:
        residuals = []
        for product in products:
            residual = product[:]
            for character in known:
                multiplicity = inner(residual,character)
                residual = [B.subtract(x,B.scale(multiplicity,y)) for x,y in zip(residual,character)]
            if inner(residual,residual):
                residuals.append(residual)
        candidates = residuals+[[B.subtract(x,y) for x,y in zip(a,b)] for a in residuals for b in residuals]
        found = False
        for candidate in candidates:
            if inner(candidate,candidate) != 1:
                continue
            if candidate[0][0] < 0:
                candidate = [B.scale(-1,x) for x in candidate]
            assert candidate[0][1:] == (0,0,0) and candidate[0][0] > 0
            assert all(inner(candidate,k) == 0 for k in known)
            known.append(candidate)
            found = True
            break
        assert found,'extend the list of actual tensor characters'
    assert sorted(k[0][0] for k in known) == [6,15,15,21,21,24,24]
    assert sum(k[0][0]**2 for k in known) == 2520
    print('Exact complete central-3 character sector:',[k[0][0] for k in known])
    return [dict(zip(G,character)) for character in known]


def full_l60_rr(Q,data,characters,field,z):
    """Full virtual character in degree 60 from the exact central-3 sector."""
    B = BiquadraticIntegers
    omega,beta = z**28,z**12+z**24+z**48
    def convert(x):
        a,b,c,d = x
        return field.convert(a)+b*omega+c*beta+d*omega*beta
    for cls in (0,1,12,14):
        a,b,_ = triples[cls]
        generators = (a,inv(b),mul(b,a))
        matrices = [data[g] for g in generators]
        es = [order(g) for g in generators]
        powers = []
        for g,S,e in zip(generators,matrices,es):
            row = []
            h,actual = E,Q.identity
            for j in range(1,e):
                h,actual = mul(h,g),Q.multiply(actual,S)
                phase = None
                for shift,scalar in enumerate(((1,0),(0,1),(-1,-1))):
                    scaled = Q.matrix([[Q.scalar_multiply(scalar,x) for x in line] for line in data[h]])
                    if actual == scaled:
                        phase = shift
                        break
                assert phase is not None
                row.append((h,phase))
            powers.append(row)
        c = [scalar_exponent(Q,Q.power(S,e),3) for S,e in zip(matrices,es)]
        cp = scalar_exponent(Q,Q.multiply(Q.multiply(matrices[0],matrices[1]),matrices[2]),3)
        for orientation in (1,-1):
            found = []
            for eps in (1,2):
                options = [[(84*k+28*eps*t)//e for k in range(e)] for e,t in zip(es,c)]
                for datum in itertools.product(*options):
                    degree = (2520*orientation*(sum(Fraction(q,84) for q in datum)-Fraction(eps*cp,3)))%2520
                    if degree != 60:
                        continue
                    multiplicities = []
                    for ch in characters:
                        m = field.convert(Fraction(-75*ch[E][0],2520))
                        for e,q,row in zip(es,datum,powers):
                            for j,(h,phase) in enumerate(row,1):
                                value = B.conjugate(ch[h]) if eps == 1 else ch[h]
                                conjugate_trace = convert(value)*z**((-28*eps*phase)%84)
                                m += z**((q*j)%84)*conjugate_trace/(1-z**((-orientation*84*j//e)%84))/e
                        coefficients = m.to_list()
                        assert len(coefficients) <= 1
                        integer = int(coefficients[0]) if coefficients else 0
                        assert m == field.convert(integer)
                        multiplicities.append(integer)
                    if multiplicities[0] == 1:
                        result = [(ch[E][0],m) for ch,m in zip(characters,multiplicities) if m]
                        assert sum(d*m for d,m in result) == -75
                        assert [(d,m) for d,m in result if m>0] == [(6,1)]
                        found.append(result)
            assert len(found) == 1
            print(f'Full exact RR L60, class {cls}, orientation {orientation:+d}: {found[0]}')


def main():
    (Q6,X6,Y6),(Q4,X4,Y4) = atlas_generators()
    six = lifts(Q6,X6,Y6)
    four = lifts(Q4,X4,Y4)
    field,z = cyclotomic_field()
    omega = z**28
    beta = z**12+z**24+z**48
    assert omega**2+omega+1 == field.zero
    assert beta**2+beta+2 == field.zero
    for name,Q,data,w in [('3.A7 six',Q6,six,omega),('2.A7 four',Q4,four,beta)]:
        traces = [field.convert(a)+field.convert(b)*w for a,b in map(Q.trace,data.values())]
        # On these quadratic subfields conjugation is w -> -1-w.
        conjugates = [field.convert(a-b)-field.convert(b)*w for a,b in map(Q.trace,data.values())]
        assert sum((c*d for c,d in zip(traces,conjugates)),field.zero) == field.convert(2520)
        print(f'{name}: every Cayley edge central; exact character norm 1')
    spin_molien = molien(Q4,four,field,beta,2,40)
    assert spin_molien == {0:1,8:1,12:1,14:1,16:1,18:1,20:2,22:1,24:3,26:2,28:3,30:3,32:5,34:3,36:6,38:5,40:7}
    print('Exact spin Molien through 40:',spin_molien)
    six_molien = molien(Q6,six,field,omega,3,30)
    assert six_molien == {0:1,3:1,6:3,9:5,12:11,15:18,18:33,21:53,24:86,27:130,30:197}
    print('Exact exceptional-six Molien through 30:',six_molien)
    twisted_six_rr(Q6,six,field,z)
    spin45_local(Q4,four)
    ordinary_rr(field,z)
    characters = central_three_characters(Q6,six,Q4,four)
    full_l60_rr(Q6,six,characters,field,z)


if __name__ == '__main__':
    main()
