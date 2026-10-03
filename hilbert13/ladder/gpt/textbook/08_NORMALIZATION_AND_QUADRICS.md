# 8. Normalization defects and quadratic algebras

## 8.1 Local defects

Let $\nu:C\to Y\subset\mathbb P^r$ be a finite birational map, with $C$ smooth and $Y$ integral, nondegenerate and of degree $n$. Set
$$D_\nu=p_a(Y)-g(C),\qquad B_\nu=\pi(n,r)-g(C).$$
The normalization sequence gives $D_\nu=\sum_y\delta_y\le B_\nu$.

**Lemma 8.1 (first jets; P).** At an image point $y$, let $b_y$ be the number of preimages and $a_y$ the rank of their tangent columns in an affine chart. If $k_y$ branches have zero derivative, then
$$\delta_y\ge2b_y-1-a_y\ge b_y-1+k_y.$$

*Proof.* The completed normalization is $\prod_{i=1}^{b_y}\mathbb C[[t_i]]$. Its quotient by the product of the ideals $(t_i^2)$ has dimension $2b_y$. The image of the completed local ring has one common constant coefficient and an $a_y$-dimensional space of linear coefficients. The normalization quotient surjects onto this first-jet quotient. Finally $a_y\le b_y-k_y$. $\square$

**Theorem 8.2 (orbit bounds; P).** Suppose a finite group acts equivariantly. The nonimmersion set $R$ and the set $S_{\rm col}$ of points in nonsingleton fibres are invariant, finite, and satisfy
$$|R|\le D_\nu\le B_\nu,\qquad |S_{\rm col}|\le2D_\nu\le2B_\nu.$$
Consequently all point orbits larger than $B_\nu$ imply immersion. All point orbits larger than $2B_\nu$ imply a closed embedding. When immersion is already known, separation need only be checked on the union of orbits of size at most $2B_\nu$.

*Proof.* Sum the bound $\delta_y\ge k_y$ for the first inequality. In a collision fibre, $b_y\le2(b_y-1)\le2\delta_y$, proving the second. An injective immersion of a smooth proper curve is a closed embedding. $\square$

For $J=G_y$, with point stabilizers $I_i$ for the $J$-orbits over $y$, one can retain stabilizer information:
$$b_y=\sum_i[J:I_i],\qquad
D_\nu\ge\sum_{[y]}[G:J](2b_y-1-a_y).$$
Higher jets replace the dimension $2b_y$ by $\sum_i m_i$ and the first-jet image rank by the rank of the local functions modulo $(t_i^{m_i})$. This yields a finite local linear-algebra refinement of the same defect inequality.

## 8.2 Application to the degree-60 model

The twisted Lefschetz calculation supplies a degree-60 invariant bundle $L_{60}$ and an exceptional six-dimensional section subspace $V$. Its common base divisor is invariant and has degree at most 60, less than the minimum point-orbit size 360; it is therefore empty.

**Theorem 8.3 (P+X).** The resulting morphism $\varphi:C\to\mathbb P(V^*)$ is a closed embedding of degree 60.

*Proof.* If its degree onto the image were $e\ge2$, the normalized image would have degree $60/e\le30$ and genus at most $\pi(30,5)=91$. Its $A_7$ action is faithful: simplicity makes the kernel either trivial or all of $A_7$, and the latter would force $2520\mid e$. The minimum faithful genus 136 excludes this possibility, so $e=1$.

Arithmetic Castelnuovo gives $p_a(Y)\le\pi(60,5)=406$, hence $D_\varphi\le270$. The only point orbit sizes are $2520,1260,630,360$. Theorem 8.2 excludes nonimmersion and confines any collisions to the single orbit $D_7$.

An order-seven lift on $V$ has characteristic polynomial $\Phi_7$. Exact specialization at $p=43$ proves that each of its six projective eigenlines has stabilizer exactly $C_7$. Since $D_7=G/C_7$, its equivariant map to the orbit of the image eigenline is bijective. This excludes the remaining collisions. $\square$

**Certificate.** Integral ATLAS matrices lie in $\mathbb Z[\omega]$, $\omega^2+\omega+1=0$. At $p=43$, specialize $\omega$ to each of $6,36$ and use a primitive seventh root $\zeta$. The six spectral projectors
$$P_j=\frac17\sum_{k=0}^6\zeta^{-jk}M^k$$
are nonzero and have rank one. Check all 2520 elements against every projector image, both for the representation and its dual. The 24 stabilizer checks return exactly seven elements. Stabilizing a characteristic-zero eigenline is a collection of minor identities preserved by specialization. Thus the finite-field stabilizer is an upper bound; the cyclic subgroup is a lower bound. Nonzero projectors prevent loss of the eigenline on specialization. This is the mathematical implication verified by `verify_normalization.py`.

**Corollary 8.4 (P+X).** The anti-invariant two-dimensional section space of an order-two lift has fixed divisor exactly $R_\tau$, every zero simple. Its moving pencil has degree 42 and descends to degree 21 on $D$.

*Proof.* The lift has section eigenspace dimensions four and two and fixed-point fibre character $+1$. A common anti-invariant zero maps into $\mathbb P(E_+)$, so $\varphi(\tau p)=\varphi(p)$ and injectivity gives $\tau p=p$. Conversely every fixed point is a common zero. At a fixed point the tangent character is $-1$, and in the ambient tangent space the negative directions point towards $E_-$. Immersion supplies an anti-invariant section with a nonzero first derivative. Removing the 18 simple zeros leaves $60-18=42$. Ratios of anti-invariant sections are invariant, hence factor through the degree-two quotient. $\square$

This proves the degrees of these particular pencils. It does not assert that either computes gonality. Immersion proves only the first two vanishing orders $(0,1)$; higher branch orders remain a separate question.

## 8.3 The two quadratic representations

Let $S=\operatorname{Sym}V$. Exact character arithmetic gives
$$S_2=Q_6\oplus Q_{15},\qquad Q_6\simeq V^*,$$
with nonisomorphic irreducible summands, and $S_3^G$ is one-dimensional. Central scalars act on sections, so invariance here means invariance for the relevant Schur cover whenever the central character is nontrivial.

**Proposition 8.5 (P+X).** Both $Q_6$ and $Q_{15}$ have empty projective base locus. Every nonempty invariant closed subset of $\mathbb P(V^*)$ lies on no quadric.

*Proof.* The exact character projector is
$$P_6=\frac6{2520}\sum_{g\in G}\operatorname{tr}(M_g)\operatorname{Sym}^2(M_g),\qquad P_{15}=1-P_6.$$
A change of lift by $\omega$ multiplies the two factors by $\omega$ and $\omega^2$, so the summand is unchanged. In both specializations and both dual conventions, exact Gaussian elimination gives
$$Q_6S_5=S_7\quad(\operatorname{rank}=792),\qquad
Q_{15}S_2=S_4\quad(\operatorname{rank}=126).$$
A full-rank minor modulo a prime is nonzero in characteristic zero. Each ideal therefore contains every form in a positive degree and has empty projective zero set. A nonzero invariant quadratic ideal contains one of the irreducible summands, by multiplicity freeness, and so cannot vanish on a nonempty subset. $\square$

## 8.4 Jacobian and apolar algebras

Let $F\in S_3$ span the invariant cubic, and $F^\vee\in\operatorname{Sym}^3V^*$ span the dual invariant cubic.

**Theorem 8.6 (P+X).** The six derivatives of $F$ span $Q_6$. They are a regular sequence, the cubic hypersurface is smooth, and
$$\operatorname{Hilb}(S/(Q_6),t)=(1+t)^6.$$
The quadrics $Q_{15}$ generate the apolar ideal of $F^\vee$, and
$$\operatorname{Hilb}(S/\operatorname{Ann}(F^\vee),t)=1+6t+6t^2+t^3.$$
The latter quotient is an Artin Gorenstein algebra of socle degree three.

*Proof.* Differentiation gives a nonzero equivariant map $V^*\to S_2$; irreducibility identifies its image with $Q_6$. Empty base locus means the partials have no common projective zero. Their ideal has height six in the Cohen--Macaulay polynomial ring, so the six quadrics are a regular sequence. The complete-intersection series is $(1-t^2)^6/(1-t)^6=(1+t)^6$, and the absence of a gradient zero proves smoothness.

Second differentiation of $F^\vee$ gives a nonzero equivariant map $S_2\to V^*$, necessarily of rank six, so its kernel is $Q_{15}$. In degree three the multiplication $Q_{15}S_1\to S_3$ has rank at least 55 by specialization. It has rank at most 55 because $Q_{15}\otimes V$ has no invariant constituent, whereas $S_3$ has an invariant line. Thus the rank is exactly 55. Degree four is all of $S_4$ by Proposition 8.5. The quotient by the generated quadrics has dimensions $1,6,6,1$ and zero thereafter. The derivative algebra of a nonconical cubic has the same dimensions, and the inclusion $(Q_{15})\subset\operatorname{Ann}(F^\vee)$ is therefore equality. Apolar duality supplies its perfect pairings between complementary degrees and its one-dimensional socle. $\square$

These are ambient ideals. Neither is the ideal of $\varphi(C)$, and neither makes the curve a complete intersection. The general normalization sequence, Castelnuovo bound, Jacobian regular-sequence argument and apolar duality are classical structures.

## 8.5 Klein contact divisors

There are two conjugacy classes of Klein subgroups $H\simeq\mathrm{PSL}_2(7)$, each of size 15. The exceptional six restricts irreducibly, and its invariant quadratic forms form a line, generated by $Q_H$.

**Theorem 8.7 (P+X).** Let $O_H$ be the unique $H$-orbit of 24 seven-points. Then
$$\operatorname{div}(Q_H|_C)=5O_H,\qquad 2L_{60}\sim5O_H.$$

*Proof.* The restriction is nonzero by Proposition 8.5 and has degree 120. Exact subgroup orbit enumeration gives

| Branch fibre | $H$-orbit sizes |
|---|---|
| $D_2$ | $84,84,84,168,168,168,168,168,168$ |
| $D_4$ | $42,84,168,168,168$ |
| $D_7$ | $24,168,168$ |

A point outside the branch fibres has an $H$-orbit of size 168. The only nonnegative solution of $120=24a+42b+84c+168d$ is $a=5$, $b=c=d=0$. This proves the divisor identity. The two Klein subgroups containing a fixed $C_7$ represent the two classes; their exact orbit checks are sufficient by conjugacy. $\square$

For each class the 15 sets $O_H$ partition $D_7$. Thus their quadratic product has divisor $5D_7$. The exact bundle identity $6L_{60}\sim D_7$ identifies this product, up to a nonzero scalar, with the fifth power of the canonical section of $6L_{60}$. Also $5(O_H-O_K)\sim0$. Distinct divisors give nonzero classes once $\operatorname{gon}(C)\ge25$ is used: a principal difference would give a pencil of degree at most 24. Their order is then exactly five. This last nonvanishing depends on the computer-assisted lower bound, unlike the contact multiplicities.

The $H$ representation admits both $V|_H\simeq\operatorname{Sym}^2(3_H)$ and $V|_H\simeq\wedge^2W_H$, with $W_H$ the four-dimensional Weil representation of $\mathrm{SL}_2(7)$. Uniqueness identifies $Q_H$ with the Plücker quadric. Its intersection with the Veronese surface is the Veronese image of Klein's quartic. A Sylow-seven subgroup permutes its three flex eigenlines with exponents proportional to $\{1,2,4\}$; the local characters of $L_{60}$ place the three corresponding fixed points on one of the two conjugate Veronese triangles. Consequently $\varphi(O_H)$ is the Veronese image of the 24 flexes. The equality with a Hessian section of the ambient cubic is an additional numerical identification, and does not follow from these exact divisor identities.
