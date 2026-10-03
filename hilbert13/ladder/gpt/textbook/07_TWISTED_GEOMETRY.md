# 7. Twisted Picard classes and the exceptional model

## 7.1 Central characters

The Schur multiplier of $A_7$ is $\mathbb Z/6$. An invariant line bundle can be linearized for the cover $6.A_7$, with its center acting by a prescribed character. The order of that character is the order of its multiplier class. The corresponding irreducible degree lists are:

| Central-character order | Irreducible degrees in a fixed faithful sector |
|---|---|
| 2 | $4,4,14,14,20,20,36$ |
| 3 | $6,15,15,21,21,24,24$ |
| 6 | $6,6,24,24,36$ |

Their squared dimensions sum to 2520 in each row; their gcds are two, three and six. The order-three sector and the exceptional six have an independent exact certificate in `audit_exact.py`. The complete remaining sector lists are retained representation-theoretic inputs; rounded recognition in `twisted_rr.py` does not independently certify them.

The spin cover is realized on $V_4$ with $\wedge^2V_4$ the ordinary standard six. The triple cover has presentation
$$\begin{aligned}
3.A_7=\langle x,y,z\mid{}&z^3=1,\ [x,z]=[y,z]=1,\ x^3=y^5=(xy)^7=1,\\
&(xyxy^{-1})^2=z,\ (xy^{-2}xy^2)^2=1\rangle.
\end{aligned}$$
It has order 7560 and is perfect; quotienting by $z$ gives $A_7$. The sixfold cover is the fiber product of the double and triple covers. Exact integral ATLAS matrices realize the exceptional six over $\mathbb Z[\omega]$ and the spin four over a quadratic integer ring. Checking every Cayley edge supplies a representation of the intended central extension, beyond checking a few relations.

## 7.2 Local data and degree

Let $x_i$ be branch generators of orders $e_i=2,4,7$, and choose lifts $\hat x_i$. For a central character $\varepsilon$, a local fibre character is $\lambda_i=\exp(2\pi ir_i/e_i)$, with $r_i$ allowed to be rational and
$$\lambda_i^{e_i}=\varepsilon(\hat x_i^{e_i}).$$
Define $\rho_0$ by $\varepsilon(\hat x_1\hat x_2\hat x_3)=\exp(2\pi i\rho_0)$.

**Proposition 7.1 (P+X).** Compatible local data determine twisted invariant classes of degrees
$$\deg L=2520\left(n+\sum_i\frac{r_i}{e_i}-\rho_0\right),\qquad n\in\mathbb Z.$$

*Proof.* Use the central preimage of the triangle group in the universal covering group of $\mathrm{PSL}_2(\mathbb R)$:
$$\widetilde\Delta=\langle c_i\mid c_1^2=c_2^4=c_3^7=c_1c_2c_3=h\rangle.$$
This is a central extension by an infinite cyclic group, not a universal central extension in the perfect-group sense. The triangle group has nontrivial abelianization.

On the upper half-plane a weight $K^w$ has central character $e^{-2\pi iw}$ and local character $e^{-2\pi iw/e_i}$. Tensor with a flat character on the fiber product with the finite Schur preimage, choosing its value on $c_i$ to be $\lambda_i e^{2\pi iw/e_i}$ and its value on the finite center to be $\varepsilon$. The local-power relations give the compatibility conditions above. The product relation gives
$$\frac3{28}w\equiv\sum_i r_i/e_i-\rho_0\pmod1.$$
Its descended degree is $270w$, yielding the formula. Conversely this congruence supplies the compatible character and hence the descended bundle. The product-versus-inverse convention must be fixed consistently; reversing orientation reverses the corresponding residues. $\square$

**Corollary 7.2 (P+X).** Invariant degrees are $15\mathbb Z$, and their residues modulo 90 are fixed by multiplier:

| Multiplier | Trivial | Order two | Order three | Order six |
|---|---|---|---|---|
| Degree modulo 90 | 0 | 45 | $\pm30$ | $\pm15$ |

The enumeration is over all compatible local data, not only effective classes. Theta characteristics have degree 135, providing a consistency check. The order-three twists are not cube roots of $K_C$: the central extension modulo $h^3$ splits, as seen by the corrected lifts $c_1h,c_2h^2,c_3h^2$, so the triangle Euler class is divisible by three.

## 7.3 Twisted Lefschetz and forced sections

For a lift $\hat g$ of a nonidentity element with tangent eigenvalues $a_p$ and fibre eigenvalues $\lambda(\hat g,p)$,
$$\operatorname{tr}(\hat g\mid H^0(L)-H^1(L))
=\sum_{p\in\operatorname{Fix}(g)}\frac{\lambda(\hat g,p)}{1-a_p^{-1}}.$$
At a central element $z$, the trace is $\varepsilon(z)(\deg L-135)$. This determines the virtual section character. A positive irreducible coefficient forces that representation in $H^0$; zero or negative coefficients alone do not determine $H^0$.

**Lemma 7.3 (local-character vanishing; P).** If a homogeneous invariant of degree $k$ restricts nontrivially, its vanishing order $o_i$ at a branch point satisfies
$$\lambda_i^k=a_i^{o_i}.$$
If the sum of its least allowed local orders weighted by $2520/e_i$ exceeds $k\deg L$, the restriction is zero.

*Proof.* Compare the characters of the first nonzero local Taylor coefficient. A nonzero section cannot have more zeros than its degree. $\square$

This refines the orbit-semigroup test by residue information. Congruence fixes orders modulo $e_i$; it does not determine their exact values.

## 7.4 Classes below degree 60

**Theorem 7.4 (P+X+C).** Every invariant bundle of degree less than 60, apart from $\mathcal O_C$, has no section.

*Proof.* Negative degree and nontrivial degree zero have no sections. The invariant lattice leaves only positive degrees 15, 30 and 45. The certified gonality bound 25 and the classical inequality $\operatorname{gon}\le\operatorname{Cliff}+3$ give $\operatorname{Cliff}\ge22$. The bundles are special; Clifford bounds give $h^0\le1$ at 15, $h^0\le5$ at 30 and $h^0\le12$ at 45. The first two bounds are below the minimum faithful dimension in their respective Schur sectors, so no section is possible.

At degree 45 the spin sector forces a copy of $V_4$ or its dual. Its four-dimensional subsystem has no base divisor, since 45 is below the minimum orbit 360. It gives a nondegenerate map to $\mathbb P^3$. The normalized image action is faithful, as its map degree is below $|A_7|$. If nonbirational, the odd map degree divides 45 and is at least three, so the image degree is at most 15, with genus at most $\pi(15,3)=42<136$. Thus the map is birational of degree 45.

Exact Molien arithmetic gives the first spin invariant degrees $8,12,14,16,18$, with one invariant in each and two in degree 20. Lemma 7.3 kills the degree-14 and degree-18 invariants on the curve: their minimum zero degrees are 3150 and 3330, exceeding 630 and 810. Both vanish on the 210 spin eigenspace lines of involutions, since their eigenvalues are $\pm i$ and these polynomial degrees act there by $-1$.

The lines are distinct. In a unitary realization, one eigenspace determines its orthogonal complement and hence the projective involution. A common line would identify the involution. The two invariants are coprime: their gcd is invariant because the cover is perfect. A nonconstant gcd of degree at most 14 has degree 8, 12 or 14; division would give an invariant of degree 6, 2 or 4, all absent. Their complete intersection has degree $14\cdot18=252$, of which the distinct lines use at least 210, leaving at most 42 for a degree-45 curve. Contradiction. $\square$

The spin exclusion is for this degree-45 target model. Spin models on other signatures, or higher-degree kernel maps, are not excluded. The C label is required by the Clifford input, even though the final polynomial calculation is exact.

## 7.5 The degree-60 class

**Theorem 7.5 (P+X).** Every one of the four covers has a degree-60 invariant class of multiplier order three with
$$\chi(L_{60})=V-15-2\cdot21-24.$$
It contains an exceptional six-dimensional section space. Its order-two lift has trace two and eigenspace dimensions four and two; every involution-fixed fibre has character $+1$.

**Exact input.** `audit_exact.py` computes the entire central-three character sector in an algebraic number ring and pairs the twisted Lefschetz character with every constituent, for all four classes and both orientations. For the local datum $(0,2,6)$ the exceptional-six multiplicity is one and the virtual-character norm is seven. It also verifies the lift trace and the fixed-point fibre characters. The positive multiplicity proves existence of the six-dimensional subsystem, not that $h^0(L_{60})=6$.

The base-point-free and birational arguments are given in Theorem 8.3, followed there by the exact embedding proof. Corollary 8.4 supplies the exact degree-42 and quotient degree-21 pencils. Together with Theorem 7.4, this proves $\widetilde\mu(C)=60$; the minimum retains the lower-bound certificate, whereas construction and embedding are independent of it.

The full six-section Plücker degree is
$$6(60+5\cdot135)=4410.$$
Local residues allow minimum sequence $(0,1,2,3,4,6)$ at two- and four-points, of weight one, and $(0,1,2,3,4,5)$ at seven-points. These consume minimum total weight 1890; equality would leave one free orbit of weight 2520. Immersion fixes only the first two orders. The full minimum sequences and the proposed free flex orbit are not yet proved.

## 7.6 The ambient cubic and divisor restrictions

Exact Molien coefficients for the exceptional six begin
$$1+t^3+3t^6+5t^9+11t^{12}+18t^{15}+33t^{18}+53t^{21}+86t^{24}+\cdots.$$
The unique invariant cubic $F_3$ vanishes on the curve by Lemma 7.3. Its smoothness follows independently from the Jacobian-system proof in Chapter 8. This is the exceptional cubic fourfold with symplectic $A_7$ action identified in the work of Laza--Zheng, Yang--Yu--Zhu and Koike. The ambient cubic is prior literature.

In a seven-eigenline frame its form is
$$\begin{aligned}
F_3={}&x_1x_2x_4+\beta x_3x_5x_6
+\gamma(x_1^2x_5+x_2^2x_3+x_4^2x_6)\\
&+\delta(x_1x_3^2+x_2x_6^2+x_4x_5^2),
\end{aligned}$$
with $\gamma^3/\beta$ and $\delta^3/\beta^2$ equal to the two conjugate values $(23\mp7\sqrt{21})/16$. The coefficient matching in `twisted_rr.py` is numerical; the independent rational equation in the cited literature specifies the ambient hypersurface without relying on root recognition. In Koike's coordinates it is
$$\begin{aligned}
2\sum_{i=1}^6x_i^3
&+3\sum_{(i,j)\in\{(1,2),(3,4),(5,6)\}}(x_i^2x_j+x_ix_j^2)\\
&+2(x_2x_3x_5+x_1x_4x_5+x_1x_3x_6+x_2x_4x_6)\\
&+4(x_1x_3x_5+x_2x_4x_5+x_2x_3x_6+x_1x_4x_6)=0.
\end{aligned}$$

Odd-degree invariants vanish on the negative eigenspaces of order-two and order-four lifts. Thus the cubic contains the 105 involution lines and the 315 order-four eigenspace lines. The latter contain the corresponding fixed points.

The exact local-character argument says that nonzero invariant restrictions in degrees $6,12,18,21,24$ have divisors respectively $D_7,2D_7,3D_7,2D_4,4D_7$; it kills every invariant restriction in degrees $3,9,15$. This specifies a possible nonzero restriction and its divisor, not the existence of a particular nonzero polynomial restriction.

The class identity is
$$K_C-3L_{60}=B+T_{\rm tor},\qquad
K_C=D_2+3D_4-8D_7.$$
Consequently $3L_{60}\sim2D_4-3D_7$ and, using $4D_4\sim7D_7$, $6L_{60}\sim D_7$. These are abstract bundle identities. They do not identify a particular ambient Hessian polynomial.

## 7.7 Higher equations and normality

The exact cubic representation is
$$\operatorname{Sym}^3V=1\oplus6\oplus14_a\oplus14_b\oplus21.$$
The invariant cubic is unique in its invariant line. Uniqueness among all cubics through $C$ would additionally require the other four restriction maps to be injective.

Riemann--Roch gives
$$h^0(3L_{60})=45+h^0(B+T_{\rm tor})\ge55,
\quad h^0(4L_{60})=105,\quad h^0(5L_{60})=165.$$
The latter two use degrees 240 and 300 with $K_C$ of degree 270: at degree 240, the residual has degree 30 and is excluded by Theorem 7.4; at degree 300 it has negative degree. Thus the degree-four equality inherits that theorem's C input. The degree-five equality is ordinary nonspecial Riemann--Roch.

Numerical continuation and restriction sampling report ranks $21,55,105,165$ in degrees $2,3,4,5$, a cubic plus 15 additional quartic generators, and 60 points in a random hyperplane section of their common zero set. These ranks and the scheme-generation claim remain N. Rank 21 in degree two now follows independently from Chapter 8's no-quadric theorem. Embedding follows independently from normalization defects. Neither proof upgrades the higher sampled ranks.

Similarly, the sextic $G_6$ is numerically traced from double tangent-line zeros at seven-points; its proposed order-four-point roots, higher equations and Hessian identification remain N. Projective normality requires surjectivity onto the entire $H^0(kL_{60})$ in every relevant degree. It does not follow from an embedding or from the invariant cubic alone.

## 7.8 The degree-90 matrix series

For the ten-dimensional subsystem of $B+T_{\rm tor}$, the associated symmetric spin matrix family has generic rank three or four. Rank one would give the excluded degree-45 spin model. Rank two would give an ordinary-six subsystem in $B$ or $B+T_{\rm tor}$, after the invariant-base-divisor argument. Local characters kill power sums $p_2,p_3,p_5,p_6$ and, for the latter class, $p_7$. The image is then finite or lies in the ordered-root curve of $t^7+ut^3+v$.

For the latter curve, the scale parameter is $s=u^7/v^4$, and the ordered-root cover has degree 5040. Its branch permutations have types a seven-cycle, $(4,3)$ and a transposition; the discriminant is $-v^2(6912u^7+823543v^4)$. Riemann--Hurwitz gives
$$2g-2=5040(-2+6/7+11/12+1/2)=1380,$$
so its genus is 691 and it cannot be dominated by the genus-136 curve. Determining whether the surviving rank-three kernel map is a spin model of degree 135 is open.

## 7.9 Other signatures

For any verified twisted class with section eigenspace $U$ of dimension $r\ge2$, mismatched fixed-point fibre characters force a base divisor of degree $b_U$. Choosing a pencil and imposing $r-2$ additional conditions gives
$$\operatorname{gon}(C)\le\deg L-b_U-(r-2),$$
provided the chosen sections still define a nonconstant map. This is a general sufficient construction, not a lower bound.

The retained `twisted_survey.py` numerical character survey suggests the following bounds. Entries outside the independently exact $(2,4,7)$ construction retain their character-recognition status N unless recertified in exact arithmetic:

| Signature | Genus | Reported best degree | Construction |
|---|---|---|---|
| $(2,4,7)$ | 136 | 42 | Degree-60 triple-cover class; exact |
| $(3,3,7)$ | 241 | 48 | Degree-60 spin class, order-three eigenspace |
| One $(2,7,7)$ curve | 271 | 84 | Degree-90 spin class |
| $(4,4,4)$ | 316 | 95 | Degree-105 order-six class |
| $(2,5,7),(2,6,7),(2,7,7)$ | 199--271 | 104--108 | Degree-120 order-three classes |
| Remaining listed rigid signatures | 169--316 | 192--299 | Higher twisted classes |

The involution/Klein-eigenspace survey through degree 270 found no improvement on 42 by that construction. Exhausting one construction does not exclude unrelated pencils.
