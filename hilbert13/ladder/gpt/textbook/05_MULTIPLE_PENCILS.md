# 5. Multiple pencils and multigraded genus

## 5.1 Field generation

**Lemma 5.1 (P).** If rational functions of degrees at most $d$ generate $\mathbb C(C)$, then $g(C)\le(d-1)^2$.

*Proof.* Induct on $d$. A one-generator field is rational. Otherwise choose a maximal proper compositum $E$ of some of the fields, of index $e$ in the full field. Its generating functions have degrees at most $d/e$, so its genus is at most $(d/e-1)^2$ by induction. One further rational field generates the full field with $E$, and Castelnuovo--Severi gives
$$g(C)\le e(d/e-1)^2+(e-1)(d-1)=d^2/e+ed-3d+1.$$
This convex expression is at most $(d-1)^2$ on $2\le e\le d$. Integer floors only strengthen the bound. $\square$

**Corollary 5.2 (P+X).** Every faithful $A_7$-curve has gonality at least 13. If its gonality is at most 24, its gonal pencil's orbit fields generate the full function field.

*Proof.* The orbit compositum is $G$-stable. The action on its curve is faithful: a nontrivial kernel would be all of simple $G$ and have order dividing an index at most the pencil degree. A pencil of degree at most 12 would give a faithful orbit-compositum curve with genus at most 121 by Lemma 5.1, below the minimum 136. If a degree-at-most-24 orbit compositum were proper, its faithful quotient would have a pencil of degree at most 12. $\square$

## 5.2 Orbit sizes

Let $K$ stabilize the pencil field and let $N$ act trivially on that field. Then $N\trianglelefteq K$, $|N|$ divides the pencil degree, and $K/N$ is a finite subgroup of $\mathrm{PGL}_2(\mathbb C)$. Such subgroups are cyclic, dihedral, $A_4,S_4,A_5$.

**Lemma 5.3 (P, with the subgroup classification).** A degree-$m$ pencil on a faithful $A_7$-curve, with $m<60$, has at least 35 conjugate pencil fields, and at least 42 when $3\nmid m$.

*Proof.* Stabilizers of index less than 35 have types $A_7,A_6,\mathrm{PSL}_2(7),S_5$. For a simple stabilizer the kernel is either trivial or the whole group; its embedding in $\mathrm{PGL}_2$ is impossible and the second option would cost at least 168. For $S_5$, the possible kernels contain $A_5$ or are trivial; an embedding is impossible, and a nontrivial kernel costs 60. At index 35 the stabilizer is $(A_4\times C_3):2$. Its elementary abelian subgroup of order nine cannot inject into $\mathrm{PGL}_2$, so the kernel contains a nontrivial subgroup of order three. Thus $3\mid m$ is necessary. The subgroup list has no intervening index. $\square$

The weaker lower bound 15 follows already by excluding the index-seven $A_6$ stabilizer. Its normal kernel has order dividing $m$; simplicity makes it trivial, and then the impossible $A_6$ embedding follows. These are pencil-field stabilizers; a line-bundle stabilizer agrees with them when the complete series has dimension one.

## 5.3 Chains and dependent products

Adjoin conjugate pencil fields one at a time outside the current compositum. Each proper index drops by at least a prime factor, so a generating chain has length $t\le1+\Omega(d)$, where $\Omega$ counts prime factors with multiplicity.

**Lemma 5.4 (P).** If $f_i\notin\mathbb C(f_1,\ldots,f_{i-1})$, the $2^t$ products $\prod_{i\in I}f_i$ are linearly independent.

*Proof.* Write a relation as $Af_t+B=0$. If $A\ne0$, it puts $f_t$ in the previous field; if $A=0$, induction kills the remaining coefficients. $\square$

The Segre product is a birational nondegenerate model of degree $td$ in $\mathbb P^{2^t-1}$. Thus $g\le\pi(td,2^t-1)$. For every feasible chain with $d\le24$, exact finite arithmetic gives $\pi(td,2^t-1)\le\pi(3d,7)$ when $t\ge3$. Feasibility includes $td\ge2^t-1$. This step proves birationality, not an embedding.

**Lemma 5.5 (P).** Suppose a birational pair of degree-$d$ pencils has quadric image $Y$ of bidegree $(d,d)$. If a distinct third pencil has dependent eight products, then
$$\delta(Y)\ge c(d):=\min_{r+s\ge d}\left(\binom r2+\binom s2\right),\qquad r,s\ge0.$$

*Proof.* The relation writes the third map as a ratio of two $(1,1)$ forms. Their common component, if present, is a ruling and cancels to one of the original pencil fields, contrary to distinctness. The complete pencil therefore has a base cluster of length two. Its restriction has degree $2d$ before the fixed part and degree $d$ after it.

For two proper centers the fixed multiplicity is the sum of the curve multiplicities at them, and successive blowup genus losses give the displayed bound. For a tangent base pair there is one proper center and its first infinitely near point in the common tangent direction. The pencil has no further center, since the intersection of two $(1,1)$ members is two. Noether's local intersection formula gives fixed multiplicity at most $r+s$, where $r\ge s$ are the two successive curve multiplicities. The same successive genus losses give $\binom r2+\binom s2$. $\square$

The convex minimum has $r+s=d$ and balanced $r,s$. It equals $k(k-1)$ for $d=2k$ and $k^2$ for $d=2k+1$.

**Corollary 5.6 (P+X).** A faithful $A_7$-curve of gonality $d\le24$ satisfies
$$g\le B^*(d):=\max\{(d-1)^2-c(d),\pi(3d,7)\}.$$

*Proof.* A generating chain of length at least three gives the second term. A length-two chain gives a birational pair; the orbit contains a third. Independent products give the same Segre bound, and dependent products give the first term. $\square$

| $d$ | 20 | 21 | 22 | 23 | 24 |
|---|---|---|---|---|---|
| $(d-1)^2-c(d)$ | 271 | 300 | 331 | 363 | 397 |
| $\pi(3d,7)$ | 261 | 290 | 320 | 352 | 385 |

The original global proof combines the low-genus certificates with $B^*(24)=397$ and the four window signatures. The sharper theorem below removes the window's spectral dependency.

## 5.4 A singular-surface bound

Let $S_0\subset(\mathbb P^1)^3$ be an integral hypersurface of positive multidegree $(a,b,c)$, and $B\subset S_0$ an integral curve not contained in its singular locus. Let $v=(m_1,m_2,m_3)^T$ be its coordinate degrees and
$$M=\begin{pmatrix}0&c&b\\c&0&a\\b&a&0\end{pmatrix}.$$

**Theorem 5.7 (P).** The normalization genus satisfies
$$g(B)\le1+\frac12\left(v^TM^{-1}v+(a-2)m_1+(b-2)m_2+(c-2)m_3\right).$$

*Proof.* Normalize the surface and take its minimal resolution $f:S'\to S_0$. The pullbacks $H_i$ of the coordinate divisors have intersection matrix $M$, and $H=\sum H_i$ is nef and big. For $L=\sum(M^{-1}v)_iH_i$, the strict transform $B'$ has $B'-L$ orthogonal to each $H_i$. Hodge index gives $B'^2\le L^2=v^TM^{-1}v$.

Adjunction on the hypersurface gives $K_{S_0}=(a-2,b-2,c-2)|_{S_0}$. Write $K_{S'}=f^*K_{S_0}+Z$. The nonexceptional coefficients of $Z$ are nonpositive by the conductor formula. For an exceptional irreducible curve $E$, minimality and adjunction give $K_{S'}E=2p_a(E)-2-E^2\ge0$: a rational $(-1)$ curve has been removed. The negative-definite exceptional intersection matrix, with nonnegative off-diagonal entries, and the nonpositive conductor intersections imply that the exceptional coefficients of $Z$ are nonpositive as well. Hence $Z\le0$.

The curve $B'$ lies in neither conductor nor exceptional support, by the singular-locus hypothesis. Therefore $K_{S'}B'\le f^*K_{S_0}B'$. Arithmetic adjunction on $S'$, followed by normalization, gives $2g(B)-2\le B'^2+K_{S'}B'$. Substitute the two bounds. $\square$

For equal coordinate degrees $m$:

| Surface type | Genus bound |
|---|---|
| $(2,1,1)$ | $A(m)=m^2/2-m+1$ |
| $(2,2,1)$ | $7m^2/16-m/2+1$ |
| $(2,2,2)$ | $3m^2/8+1$ |

## 5.5 Graph surfaces and singular-locus exceptions

An integral surface of type $(a,b,1)$ is a graph of a rational map from $\mathbb P^1\times\mathbb P^1$. Resolve its entire base cluster, with pencil multiplicities $n_i$ and projected-curve multiplicities $r_i$. The moving pencil has square zero, so $\sum n_i^2=2ab$. If the first two projections are birational and all three degrees are $m$, then
$$\sum n_ir_i=(a+b-1)m=:T.$$
Successive blowups give
$$g\le(m-1)^2-\frac12\sum(r_i^2-r_i)
\le(m-1)^2-\frac{T^2}{4ab}+\frac T2,$$
using Cauchy--Schwarz and $n_i\ge1$. For $(a,b)=(2,2)$ this recovers the $(2,2,1)$ bound even when the curve lies in the graph's singular locus.

If a curve lies in the singular locus of a $(2,2,2)$ hypersurface, a nonzero first partial derivative gives a hypersurface of type $(1,2,2)$, up to permutation, containing it. Choose an irreducible component containing the curve. A zero-coordinate component is impossible for $m\ge8$ and pairwise birational projections: its two-factor image would have bidegree $(m,m)$, exceeding the available $(2,2)$ degrees. The remaining positive components are graph cases and have the preceding bound. The derivative need not be smooth or irreducible; this component argument is necessary.

## 5.6 The sharp three-pencil theorem

**Theorem 5.8 (P, with Petrakiev and uniform position).** If $C$ has three degree-$m$ pencils, $m\ge8$, every pair is birational, and their eight products are linearly independent, then
$$g(C)\le\lfloor A(m)\rfloor=\left\lfloor\frac{m^2}{2}-m+1\right\rfloor.$$
For every even $m\ge8$ equality is attained.

*Proof.* Suppose $g>A(m)$. The Segre model $B$ is birational, nondegenerate and of degree $3m$ in $\mathbb P^7$. First prove completeness of its eight product sections. A ninth section would give a nondegenerate degree-$3m$ model in $\mathbb P^8$. Petrakiev's Theorem 2.16(a), applied to its arithmetic genus, supplies a surface of degree at most eight. Its numerical threshold is
$$\pi_2(3m,8)=9\binom\lambda2+\lambda(\epsilon+2)+\mu,$$
where $\lambda=\lfloor(3m-1)/9\rfloor$, $\epsilon=3m-1-9\lambda$, $\mu=\max(0,\lfloor(\epsilon-4)/2\rfloor)$. This is at most $A(m)$: for $m=3s,3s+1,3s+2$, the differences are $s/2,(s+1)/2,s/2+1$.

Project the surface to $\mathbb P^7$. Its image is still a surface of degree at most eight: a curve of such degree could not contain $B$, whose degree is at least 24. Each Segre quadric must contain it, because otherwise its intersection degree is at most 16, again too small to contain $B$. Thus it is a divisor on the Segre threefold of degree $2(a+b+c)\le8$. A zero-coordinate type is excluded by birationality of the two-factor projections, since that two-factor image has degree $2m$. Positive types leave $(1,1,1)$, contradicting product independence, or a permutation of $(2,1,1)$, giving $g\le A(m)$ by Theorem 5.7. This proves $h^0$ of the product bundle is exactly eight, hence linear normality.

Let $\Gamma$ be a general hyperplane section, in uniform position, of $3m$ points in $\mathbb P^6$. It has $h_\Gamma(1)=7$. If $h_\Gamma(2)\le18$, linear normality gives $h_B(2)=8+h_\Gamma(2)\le26$. The Segre threefold has 27 sections of multidegree $(2,2,2)$ induced by ambient quadrics; hence $B$ lies on one such hypersurface. Theorem 5.7 gives $3m^2/8+1\le A(m)$ outside its singular locus, and the graph argument gives $7m^2/16-m/2+1\le A(m)$ inside it, for $m\ge8$. Reducible components are treated by the same zero-coordinate exclusion. Therefore $h_\Gamma(2)\ge19$.

Uniform-position subadditivity now gives
$$h_\Gamma(2j)\ge\min(3m,18j+1),\quad j\ge1,
\qquad h_\Gamma(2j+1)\ge\min(3m,18j+7),\quad j\ge0.$$
Castelnuovo's deficit sum bounds arithmetic, and hence geometric, genus. For $m=6q+r$ the summed upper bounds are:

| $r$ | Deficit sum |
|---|---|
| 0 | $18q^2-8q+1$ |
| 1 | $18q^2-2q$ |
| 2 | $18q^2+4q$ |
| 3 | $18q^2+10q+2$ |
| 4 | $18q^2+16q+5$ |
| 5 | $18q^2+22q+8$ |

Each is at most $\lfloor A(m)\rfloor$; the difference is $2q$, except for $r=2$, when it is $2q+1$. This is the final contradiction.

For sharpness, take a smooth $(2,1,1)$ surface, a del Pezzo surface of degree four, and a general smooth curve in $|-kK|$, $k\ge4$. Its three degrees are $m=2k$ and its genus is $2k^2-2k+1=A(2k)$. The two birational surface projections remain birational on the curve. For the degree-two surface projection, choose a curve not invariant under its deck involution; both eigensummands of the anticanonical sections are nonzero, so this is possible. No $(1,1,1)$ section vanishes on the curve, because the residual class has negative intersection with $-K$. All hypotheses hold. $\square$

## 5.7 The algebraic genus-266 threshold

**Theorem 5.9 (P+X).** A faithful $A_7$-curve of genus at least 266 has gonality at least 25.

*Proof.* Suppose its gonality is $d\le24$. By Corollary 5.2 the orbit fields generate the full field. If two distinct fields lay in a proper compositum, enlarge it to a maximal proper compositum of index $e$. Distinctness gives $2\le e<d$, $e\mid d$. Lemma 5.1 and Castelnuovo--Severi give
$$g\le d^2/e+ed-3d+1\le265.$$
The last inequality is exact finite arithmetic for $d\le24$, with maximum at $(d,e)=(24,2),(24,12)$. Thus every pair is birational, and $g\le(d-1)^2$ forces $d\ge18$.

An independent third would give $g\le A(d)\le265$ by Theorem 5.8. Every other pencil is therefore dependent relative to a fixed pair and gives a distinct complete length-two base cluster. Lemma 5.3 supplies at least 33 such clusters. Two disjoint clusters cost more than the available defect $(d-1)^2-266$ for every $18\le d\le24$; at $d=24$ the minimum cost is $2\cdot132=264$, but the available defect is at most 263.

Regard the clusters as edges joining their proper or infinitely near centers. Distinct edges are pairwise intersecting. Four or more such edges form a star: otherwise three form a triangle and there can be no fourth. If its center has multiplicity $r$, its distinct leaves have multiplicity $d-r$, so
$$\binom r2+33\binom{d-r}2\le(d-1)^2-266.$$
Exact integer enumeration leaves only $(d,r)=(24,23)$. An infinitely near center has a fixed proper ancestor as the other vertex of every length-two cluster, contradicting distinctness. The center is proper.

Projecting the quadric from this point gives a birational plane curve of degree $48-23=25$. A singular point of multiplicity at least two would yield, by projection, a pencil of degree at most 23, contradicting the assumed gonality 24. The plane curve is therefore smooth and has genus 276. But every faithful $A_7$-curve has genus congruent to one modulo three. This contradiction proves the theorem. $\square$

The low-genus certificates of Chapters 3--4 cover all genera below 336, so this theorem gives $\gamma(A_7)\ge25$ without the window certificates. The degree-42 example gives $\gamma(A_7)\le42$. The global lower bound retains its low-genus computer-assisted input; only the genus-at-least-266 part is purely algebraic.

## 5.8 A historical Picard-closure mechanism

Under the separate hypotheses $g>257$, gonality 18, all gonal pairs birational and at least seven classes in $W^1_{18}$, the historical manuscript proposes a closure argument. The pencil trick identifies the triple-product kernel with $H^0(A+B-D)$. A kernel of dimension at most one would give a degree-54 model in $\mathbb P^6$ of genus at most 255; hence $A+B-D\in W^1_{18}$ for distinct classes. This suggests a coset $A+H$ in the Picard group. The two-ruling model for $2A$ and $\pi(36,4)=187$ then force $H=H[2]$ once the closure and finiteness arguments are supplied. Summing an $A_7$-stable finite coset yields invariant degree $18\cdot2^a$, incompatible with the relevant free cyclic actions.

This is retained as H: it is unnecessary for Theorem 5.9 and has separate closure, cardinality and completeness obligations. The small historical bounds at genera $136,169,199,211,241$ are superseded by the certified lower bound 25. The ordinary-multiplicity branch of a $(17,17)$ curve with a 16-fold ordinary point gives a smooth plane curve of degree 18 and genus 136; the classical plane-automorphism bound $6\cdot18^2<2520$ excludes that branch only. It does not cover arbitrary infinitely near configurations.
