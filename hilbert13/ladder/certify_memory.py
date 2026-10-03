"""Memory-efficient certificate for classes 12 and 14.

The n=128 signed Q1 problem has 3,440,640 unknowns.  Replace it with
untwisted quotients and certify their SECOND eigenvalues by LDL inertia:
Ind_S5^G(1)=1+6+14a, Ind_L2(5)^G(1)=1+6+14b+21.
Together with Q0 and Q2 these see every ordinary irreducible.
The largest problem has 1,032,192 unknowns.  The constant function is the
single eigenvalue allowed below sigma; no dense rank-one deflation is used.

Run: python3 certify_memory.py [12|14].  Q0/Q2 bounds are supplied by
run_certificate.py; this file prints and checks the two replacement bounds.
"""
import sys
from fractions import Fraction
from a7 import cyc,cycle_type
from cover import character,TABLE,CLASSES
from orbifold import quotient_tiles
from certify_th import count_below,S5_GENS
from triples_data import triples

L25_GENS = ([cyc((1,2,3,4,5)),cyc((1,6),(2,5))],[1,1])


def verify_induction(generators,expected):
    subgroup = character(*generators)
    assert len(subgroup) in (60,120)
    result = {}
    for name,row in TABLE.items():
        assert all(cycle_type(g)!=(7,) for g in subgroup)
        m = sum(Fraction(row[CLASSES.index(cycle_type(g))]) for g in subgroup)/len(subgroup)
        assert m.denominator == 1 and m >= 0
        if m:
            result[name] = int(m)
    assert result == expected,result


def main(cls):
    assert cls in (12,14)
    verify_induction(S5_GENS,{'1':1,'6':1,'14a':1})
    verify_induction(L25_GENS,{'1':1,'6':1,'14b':1,'21':1})
    a,b,_ = triples[cls]
    bounds = []
    for name,generators,n,sigma in [('S5',S5_GENS,128,0.3557),('L2(5)',L25_GENS,64,0.36)]:
        reps,glue,size = quotient_tiles(a,b,*generators)
        result = count_below(reps,glue,n,sigma)
        print({'class':cls,'quotient':name,'subgroup_order':size,**result},flush=True)
        assert result['ok'] and result['negative_pivots'] == 1,result
        bound = result['lambda2_bound']
        assert bound is not None and Fraction(bound) > Fraction(48,135)
        bounds.append(bound)
    print(f'class {cls}: with Q0 and Q2, lambda_1 >= {min(bounds)} > 48/135; gonality >= 25',flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 12)
