# 9. Jacobians, arithmetic quotients and pencil symmetry

## 9.1 Rational isotypic factors

Chevalley--Weil, checked with exact character sums, gives
$$H^1(C,\mathbb C)=3(10+\overline{10})+2\cdot15+2\cdot21+4\cdot35.$$
The representations $1,6,14_a,14_b$ are absent. Rational group-algebra idempotents yield isogenies
$$\operatorname{Jac}(C)\sim A^{10}\times E_{15}^{15}\times E_{21}^{21}\times S^{35},$$
$$\operatorname{Jac}(D)\sim A^4\times E_{15}^7\times E_{21}^{11}\times S^{17},
\qquad\operatorname{Jac}(T_0)\sim E_{21}\times S.$$
Here $A$ has dimension three and multiplication by $\mathbb Q(\sqrt{-7})$, $E_{15}$ and $E_{21}$ are elliptic, and $S$ has dimension two. The quotients $C/A_5$, with $A_5$ fixing two letters, and $C/L_2(5)$, with $L_2(5)$ transitive on six letters, realize the two elliptic factors. The quotient by $C_3^2:C_4$ realizes the genus-two Jacobian $S$. Isogeny is not an isomorphism of principally polarized varieties, and the notation does not determine endomorphism rings, CM, ranks or conductors.

The permutation representation on involutions has character
$$\mathbb Q[G/C_G(\tau)]=1+6+2\cdot14_a+14_b+21+35.$$
The equivariant fixed-divisor map from this module to degree-zero divisor classes kills the absent Jacobian isotypes after tensoring with $\mathbb Q$. Integral classes in those components are therefore torsion after clearing projector denominators. The assertion concerns degree-zero combinations; a nonzero-degree invariant component is not a Jacobian point.

The Klein difference and the lifted branch relation (2.1) have nonzero projections into the unique copies of 21 and 35 in this permutation module. The exact verifier for this edition computes the integer numerators of the two character projectors, avoiding the floating rank and norm tests of `side_checks.py`. Uniqueness of the copies makes the generated nonzero modules coincide. Passing to Picard classes may still kill them by torsion.

Their quotient images are the named points $P_E\in E_{21}$ and $P_S\in S$. Infinite order of either appropriate detecting image would disprove the branch relation. A 5-adic nonvanishing argument only says that, if the original difference is torsion, its order is divisible by five. It does not prove infinite order. Degree-nine pencils are already excluded independently by ramification transport.

## 9.2 Elliptic subcovers and integral polarization

A primitive rank-two sub-Hodge lattice $\Lambda\subset H^1(C,\mathbb Z)$ corresponds to an optimal elliptic quotient. For an integral basis $e_1,e_2$, its map degree is
$$n=|\langle e_1,e_2\rangle_C|,$$
the restricted cup-product pairing. A nonoptimal map factors further through an elliptic isogeny and has larger degree.

For $\rho=15$ or 21, the multiplicity space has rank two, say $M_\rho$. Every rational plane of the form $v\otimes M_\rho$, $v\in V_\rho(\mathbb Q)$, gives the saturated lattice
$$\Lambda_v=(v\otimes M_\rho)\cap H^1(C,\mathbb Z).$$
The exact dessin construction uses 5040 triangles, a tree--cotree cocycle basis and the Alexander--Whitney cup product. It checks rank 272, antisymmetry and unimodularity. Saturation changes the cup value by a finite glue index:
$$n(v)=|\Phi(v)|/g(v).$$
Here $\Phi$ is the positive rational form from translates of the subgroup-fixed plane, and $g(v)$ records the saturation index. The glue has squarefree exponent with the relevant primes accounted for by exact modular kernels.

**Proposition 9.1 (P+X in the stated plane family).** Every elliptic subcover whose rational Hodge plane has this form in the pure 15 or 21 isotype has degree at least 60. Equality is attained by the two order-60 quotient subgroups above.

**Finite proof.** Enumerate glue congruence types. For each, test whether its exact positive Gram form has a nonzero vector with $|\Phi(v)|<60g(v)$. Rational $LDL^T$ gives positive diagonal pivots, and Fincke--Pohst recursion uses exact partial squared norms for pruning and integer-square-root bounds containing every boundary candidate. Every search is empty. The program `elliptic_subcovers.py` carries this out for all four classes; floating basis proposals do not supply the final lattice inequality. Its exact Hermite-normal-form, saturation and rational enumeration steps provide that implication.

If $\operatorname{End}(E_\rho)=\mathbb Z$, these pure tensor planes exhaust that isotypic family. Extra CM planes, mixed isotypes, or isogenies between the elliptic factors are not excluded by the proposition. Every elliptic map, without a purity assumption, has degree at least 13 once gonality 25 is used, by composing with an elliptic degree-two map to $\mathbb P^1$. The exact pure-plane calculation cannot upgrade this universal lower bound to 60.

## 9.3 A septic Belyi map

Let $21s^2-42s+5=0$ and $t=2s/3-10/21$. Put
$$q_s(z)=\frac{z^7}7-(s+1)\frac{z^6}6+(s+t)\frac{z^5}5-t\frac{z^4}4,
\qquad r_s=-\frac{128(51s-65)}{1750329},$$
$$\beta_s(z)=1-q_s(z)/r_s.$$

**Proposition 9.2 (P+X).** This realizes the degree-seven quotient $C/A_6\to C/A_7$, with passport $[2^21^3,(4,2,1),7]$, over $\mathbb Q(\sqrt{21})$.

*Proof.* Exact differentiation gives $q_s'=z^3(z-1)(z^2-sz+t)$. Reduction modulo the defining polynomial for $s$ gives $q_s(1)=0$ and $q_s=r_s$ at both quadratic critical points. The discriminant of $\beta_s(z)-T$ is
$$-7r_s^{-6}T^2(T-1)^4.$$
Thus the only finite branch values are zero and one, with the stated ramification; infinity is a seven-fold point. All branch permutations are even, so geometric monodromy lies in $A_7$. The remaining proper transitive candidate with this passport is the Klein subgroup. Reduction modulo 17 with $s=3$ gives a specialization with factor degrees two and five, hence an element of order ten, absent from the Klein subgroup and its normalizer. This distinguishes $A_7$. The exact checks occur in `side_checks.py`. $\square$

The equal-critical-value condition factors as
$$ (7s^2-7s+4)(21s^2-56s+40)(21s^2-42s+5). $$
The first factor gives the Klein case. An integral rescaling of the $A_7$ polynomial is
$$P_\lambda=X^7-(\lambda+21)X^6+36(\lambda-6)X^5-324(\lambda-15)X^4,
\quad\lambda^2-42\lambda+105=0.$$
The historical field-of-marked-action statement adjoins $\sqrt{-7}$ as well; this concerns the marking, rather than the displayed unmarked Belyi polynomial.

## 9.4 Quotient resolvents and the arithmetic problem

Over $\beta(X)=\beta(z)$, let $r_1,\ldots,r_6$ be the other roots. For a perfect matching $M$ of these six labels write $\sigma_M=\sum_{ij\in M}r_ir_j$. The historical resolvent for the transitive $A_5$ quotient is
$$R_E(U)=\prod_{\Pi}\left(U-\sum_{M\in\Pi}\sigma_M^3\right),$$
over the six pentads of perfect matchings. The cube avoids the constant lower-power sums. For the $C_3^2:C_4$ quotient it is
$$R_S(V)=\prod_{A\mid A^c}\left(V-\prod_{i\in A}r_i-\prod_{i\in A^c}r_i\right),$$
over the ten unordered splittings into triples. The combinatorial stabilizers identify the intended quotient fields. Turning these definitions into normalized low-genus models, checking the primitive chosen resolvent factor and locating the named divisor images are additional arithmetic work. No Weierstrass model, canonical height or torsion order is supplied here.

The classical arithmeticity of the triangle group and the noncongruence status of this kernel, cited historically, do not imply torsion of these compact elliptic-point divisors. A Manin--Drinfeld statement about congruence cusps cannot be applied merely by analogy.

## 9.5 Class stabilizers

**Proposition 9.3 (P).** Let $f$ compute gonality $d$ and $L=f^*\mathcal O(1)$. Let $K=\{g:g^*L\simeq L\}$ and $N$ be the kernel of its projective action on $H^0(L)$. Then $K/N$ is a finite Möbius group,
$$d\ge|N|\operatorname{gon}(C/N),$$
and the line-bundle orbit has size $[G:K]$.

*Proof.* A gonal bundle has $h^0=2$. Its isomorphisms are unique up to scalar, so $K$ acts on the pencil projectively. The kernel fixes the section ratio, giving descent through $C/N$. Orbit--stabilizer completes the proof. $\square$

For the displayed involution pencil its class stabilizer is $C_G(\tau)$. Indeed $L_{60}$ is invariant and its fixed divisor is $R_\tau$; a different involution with equivalent fixed divisor would give a degree-at-most-18 function, against gonality 25. Distinct involutions have disjoint fixed sets. The exact action on the two anti-invariant sections has scalar kernel precisely $\langle\tau\rangle$. Thus the centralizer quotient acts faithfully on this particular pencil, as a dihedral group of order 12. It need not be a gonal pencil; the class calculation does not assume that.

## 9.6 Quotient propagation

**Lemma 9.4 (P).** For $S_0<T_1\le G$ of index $m$, a degree-$d$ pencil on $C/S_0$ either factors through $C/\langle S_0,j\rangle$ for some $j\in T_1\setminus S_0$, or satisfies
$$g(C/S_0)\le m\,g(C/T_1)+(m-1)(d-1).$$

*Proof.* Apply Castelnuovo--Severi to the pencil and quotient map. A nonbirational common field corresponds, by Galois correspondence for $C\to C/T_1$, to an intermediate subgroup properly containing $S_0$. Normality of $S_0$ in $T_1$ is unnecessary. $\square$

The retained canonical-curve tests have numerical inputs. Max Noether gives multiplication rank $2g-1$ in the hyperelliptic case and $3g-3$ otherwise. Petri's quadric base locus is a surface for a trigonal curve or a plane quintic, and otherwise the curve itself. The latter exception occurs only in genus six and is absent from the listed quotient genera. Green--Lazarsfeld nonvanishing says $\operatorname{Cliff}\le p$ implies $K_{p,2}\ne0$; a verified vanishing therefore gives gonality at least $p+3$ in this range. Sampled analytic differentials and singular values do not verify exact rank or Koszul vanishing.

Conditional on those inputs, subgroup propagation gives:

| Group | Genus | Conditional gonality lower bound | Degree cost |
|---|---|---|---|
| $C_2$ | 64 | 13 | 26 |
| Two $C_3$ classes | 46 | 9 | 27 |
| $C_4$ | 32 | 9 | 36 |
| Two $V_4$ classes | 28 | 7 | 28 |
| $C_5$ | 28 | 9 | 45 |
| $C_6$ | 22 | 8 | 48 |
| Two $S_3$ classes | 19 | 6 or 5 | 36 or 30 |
| $C_7$ | 19 | 4 | 28 |
| $D_8$ | 12 | 6 | 48 |
| $C_3^2$ | 16 | 6 | 54 |
| $D_{10}$ | 10 | 5 | 50 |
| Order-12 types | 7--11 | 4--6 | At least 48 |
| $C_3^2:2$ | 4 | 3 | 54 |
| $C_5:C_4$ | 5 | 4 | 80 |
| $C_7:C_3$ | 7 | 2 | 42 |
| Listed order-24--60 types | 1--4 | 2--3 | At least 48 |

The initial ten numerical inputs are nonhyperellipticity for $C_3^2:2$; nontrigonality for $C_5:C_4,D_{12},C_3:C_4$, two $A_4$ and $C_6\times C_2$; $K_{2,2}=0$ for $D_{10},C_3^2$; and $K_{3,2}=0$ for $D_8$. Group enumeration and propagation are exact, but the resulting stronger table retains N. The universal initial bound $|K|\operatorname{gon}(C/K)\ge25$ retains the spectral C input.

If all these inputs are recertified, a gonal pencil of degree at most 41 has kernel among
$$1,C_2,C_3,C_4,V_4,S_3,C_7,$$
with the two conjugacy types where indicated above. For $N=C_7$, normality gives $K\le C_7:C_3$. For $N=1$, the quotient pencil $C/K\to\mathbb P^1/K$ still has degree $d$, and its fibres over a target branch value have multiplicities divisible by its branching order away from source stabilizers. This is a symmetry classification conditional on the numerical quotient inputs, not a proof of gonality 42.
