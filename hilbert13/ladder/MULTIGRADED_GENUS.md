# The intersection matrix behind the weighted base-point budget

This continuation isolates the structure behind `GPT_CONTINUATION.md` Sections
2--3. It is a deduction from finite duality, minimal surface resolution,
adjunction, and Hodge index. No claim of bibliographic priority is made.

## 1. Relative canonical divisor on normalization and minimal resolution

Let S be an integral projective Gorenstein surface over C. Let
pi:Y -> S be the composition of normalization and its minimal resolution.
Choose compatible canonical divisors. Then

    K_Y - pi^*K_S <= 0.                            (RC)

The assertion concerns Cartier pullback from S and a divisor on the smooth
surface Y; it makes no Q-Gorenstein assumption on the normalization.

To prove it, finite duality identifies the normalized dualizing module with
the conductor times the local dualizing generator on S. Thus the relative
canonical divisor D=K_Y-pi^*K_S has nonpositive coefficients on every
nonexceptional curve of Y. Write D=A-N with A,N effective and no common
components. Any component of A is exceptional for the resolution.

For every exceptional irreducible E, minimality excludes a smooth rational
(-1) curve. Since E^2 is negative, adjunction yields

    D.E = K_Y.E = 2p_a(E)-2-E^2 >= 0.

Here pi^*K_S.E=0 because K_S is Cartier and E maps to a point. If A were
nonzero, negative definiteness of the exceptional intersection matrix and
nonnegativity of intersections of distinct curves on Y would give

    0 <= D.A = A^2-N.A < 0.

Therefore A=0, proving (RC). The trivial-dualizing lemma in the main
checkpoint is the special case K_S=0. The conductor identity is
Niu--Ulrich, [arXiv:1404.5092](https://arxiv.org/abs/1404.5092), Lemma 3.11.

## 2. A genus bound for any degree vector

Let S be an integral hypersurface of multidegree (a,b,c) in X=(P1)^3,
with a,b,c>0. Let B be an integral curve on S not contained in its singular
locus, and let C be its strict transform on Y as above. Put m_i=F_i.C,
where F_i is the pulled-back coordinate fibre class. These are the degrees
of the coordinate maps on the normalization of B. They need not be equal.

The coordinate intersection matrix is

    M = [[0,c,b],[c,0,a],[b,a,0]],    det(M)=2abc.

It has signature (1,2). Its inverse is

    M^-1 = 1/(2abc) [[-a^2,ab,ac],[ab,-b^2,bc],[ac,bc,-c^2]].

Write v=(m_1,m_2,m_3)^t. The real divisor C0=sum_i (M^-1 v)_i F_i
has the same intersections with all F_i as C. Their span contains the
nef and big class H=F_1+F_2+F_3, with H^2=2(a+b+c)>0. Hence Hodge
index makes C-C0 have nonpositive square, so

    C^2 <= v^t M^-1 v.

Adjunction on X gives K_S=(a-2)F_1+(b-2)F_2+(c-2)F_3 on S.
By (RC), and because C is not a component of its negative relative
canonical divisor,

    C.K_Y <= (a-2)m_1+(b-2)m_2+(c-2)m_3.

The support of that relative divisor maps into the conductor or singular
locus; this is where the hypothesis on B is used. Adjunction on Y and
the nonnegative singularity defect of C now give the general bound

    g(normalization B)
    <= 1 + [v^t M^-1 v
             +(a-2)m_1+(b-2)m_2+(c-2)m_3]/2.      (MG)

This remains valid for nonnormal S and singular normalizations. It does not
assert the same formula for a curve contained in the singular locus.

## 3. Equal coordinate degrees

When m_1=m_2=m_3=m, (MG) becomes

    g <= 1 + m^2/[4abc] * [2(ab+ac+bc)-a^2-b^2-c^2]
             +(a+b+c-6)m/2.

The three relevant cases are:

| surface | quadratic form v^t M^-1 v | genus bound |
|---|---|---|
| (2,1,1) | m^2 | m^2/2-m+1 |
| (2,2,1) | 7m^2/8 | 7m^2/16-m/2+1 |
| (2,2,2) | 3m^2/4 | 3m^2/8+1 |

Thus the previously separate rational-graph and K0 bounds are the same
intersection-theoretic calculation. The graph proof has an additional
advantage: it covers a curve in the singular locus whenever its projection
to P1 x P1 is birational. The general normalization proof requires the
stated exclusion of that locus.

For type (a,b,1), resolving the graph pencil gives sum n_i^2=2ab and
T=sum n_i r_i=(a+b-1)m. Its base-point inequality reads

    g <= (m-1)^2-T^2/(4ab)+T/2.

Expanding this expression gives exactly (MG). The weighted singularity
budget is therefore a concrete blowup realization of the coordinate
intersection matrix and its negative orthogonal complement.

The extremal class for a (2,1,1) surface is also explained. For balanced
v, M^-1 v=(0,m/2,m/2). Thus equality in Hodge index puts the curve
in the numerical class (m/2)(F_2+F_3)=-(m/2)K. For even m=2k,
general smooth curves in |-kK| on a smooth degree-four del Pezzo surface
realize genus m^2/2-m+1 and three pairwise birational degree-m maps.

## 4. The degree-eight borderline and lower-degree bounds

If a balanced curve with pairwise birational projections lies in the
singular locus of an integral (2,2,2) surface, it lies on a nonzero
homogeneous partial derivative of type (2,2,1), up to permutation.
An irreducible surface component containing it has positive coefficients
by pairwise birationality. Its type is therefore (1,1,1), (2,1,1), or
(2,2,1). The first would make the eight product sections dependent.
The graph proof bounds the others, even if they are themselves singular.

For m>=8, their genus bounds are at most A(m)=m^2/2-m+1.
So under the contrary assumption g>A(m) the curve avoids the singular
locus of every relevant (2,2,2) surface. This supplies (MG) at m=8,
where the earlier degree>24 shortcut did not apply. Petrakiev's P8
range and the Hilbert-function argument still apply. Consequently the
sharp independent-triple theorem holds for **m>=8**, with the same bound.

The same method also bounds independent triples in the smaller degrees:

| m | genus upper bound |
|---:|---:|
| 5 | 9 |
| 6 | 14 |
| 7 | 19 |
| 8 | 25 |

For m=6, a ninth section would produce a degree-18 P8 curve, whose
Castelnuovo bound is 13. Above genus 14 this is impossible, giving linear
normality of the P7 model. Its hyperplane section has only 18 points and
therefore cannot have h(2)>=19; an additional (2,2,2) divisor follows.
The graph cases have genus at most 13 and the remaining normalization
bound is floor(1+3*36/8)=14.

For m=7 the P8 pi_2 theorem applies; the three surface bounds have maximum
floor(1+3*49/8)=19, while the no-additional-quadric Hilbert bound is 16.
The same singular-locus exclusion works under g>19.

For m=5 use Eisenbud--Harris in P7 directly: g>pi_1(15,7)=9 forces a
surface of degree six. Its containment in X follows from 15>2*6.
An integral surface of degree six in X containing a pairwise-birational
balanced curve has positive type (1,1,1), hence is a hyperplane section,
contradicting independent products.

These smaller-degree statements concern independent triples. They do not
by themselves exclude degree-six pencils on a genus-fourteen Hurwitz curve.
