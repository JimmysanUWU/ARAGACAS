# 10. Conformal variation, other groups and historical reductions

## 10.1 Conformal eigenvalue optimization

Fix a conformal class on a compact surface and a finite conformal action of $G$. Set
$$\Lambda_G(C)=\sup_{h>0,\ h\ G\text{-invariant}}\lambda_1(hg)\int_C h\,dA_g.$$
Dirichlet energy is conformally invariant in dimension two. The balanced degree--energy argument of Theorem 3.1 therefore gives
$$\Lambda_G(C)\le8\pi\operatorname{gon}(C).$$

**Proposition 10.1 (first-eigenspace ceiling; P).** Suppose $E_1(g)^G=0$. For an orthonormal basis $\varphi_1,\ldots,\varphi_q$ of $E_1(g)$, put
$$F=\frac A q\sum_i\varphi_i^2,\qquad A=\operatorname{Area}(g).$$
Then $F$ is invariant, has mean one, and every invariant conformal factor satisfies
$$\lambda_1(hg)\operatorname{Area}(hg)\le\frac{\lambda_1(g)A}{\min_C F}$$
when the denominator is positive.

*Proof.* Each $\varphi_i$ has zero $h\,dA_g$ mean: integration against an invariant measure is an invariant linear functional on $E_1$. Sum their Rayleigh numerators and denominators. The mediant bounds $\lambda_1(hg)$ by $\lambda_1(g)A/\int hF$. Multiplication by $\int h$ proves the assertion. $\square$

The former numerical range $F\approx[0.998,1.001]$ was an uncertified sampling calculation. It gives neither a proved global minimum nor a rigorous percentage ceiling. A nearly constant sum of squares need not have finite spectral support.

**Proposition 10.2 (noncriticality for a hyperbolic triangle action; P).** No invariant positive semidefinite form on $E_1$ has constant squared frame length one. There is an invariant conformal variation strictly increasing $\lambda_1\operatorname{Area}$.

*Proof.* A constant frame length would produce a harmonic map $\Phi:C\to S^N$ with $\Delta\Phi=\lambda_1\Phi$ and constant energy density $|d\Phi|^2=\lambda_1$. Its Hopf differential is invariant. Such a differential descends to a quadratic differential on $\mathbb P^1$ with at most simple poles at three branch values; that line bundle has degree $-4+3=-1$, so the differential is zero. Thus $\Phi$ is a conformal minimal immersion with induced metric $(\lambda_1/2)g$ and curvature $-2/\lambda_1$. Bryant's theorem excluding even local minimal surfaces of constant negative curvature in a finite-dimensional sphere rules this out.

The invariant positive semidefinite forms with normalized integral form a compact convex set. Their squared frame lengths form a compact convex set of smooth invariant functions of mean one, which does not contain the constant one. Separation supplies a smooth invariant function $u$ for which the first variation is positive on every normalized rank-one form in $E_1$, uniformly. The variational description of the lowest eigenvalue branches then makes $\lambda_1((1+tu)g)\operatorname{Area}((1+tu)g)$ strictly increase for sufficiently small positive $t$. This argument also handles a reducible first eigenspace. $\square$

If $E_1$ is absolutely irreducible, Schur's lemma makes the variation especially simple:
$$\left.\frac d{dt}\right|_{0}\bigl[\lambda_1((1+tu)g)\operatorname{Area}((1+tu)g)\bigr]
=\lambda_1\int_C u(1-F)\,dA.$$
Taking $u=1-F$ gives a strictly positive derivative. This establishes strict gain, with no numerical lower bound for its size.

**External existence theorem.** Vinokurov's Corollary 1.4 gives an equivariant harmonic-map maximizer for the first eigenvalue when every orbit has size greater than one. The four covers have orbit sizes $360,630,1260,2520$, so the theorem applies. The maximizing metric may have isolated conical singularities. The triangle Hopf argument makes the harmonic map branched conformal. It does not fix the target representation, the branch support, the optimum or its uniqueness. The precise primary reference is *Conformal optimization of eigenvalues on surfaces with symmetries*, [arXiv:2502.03756](https://arxiv.org/abs/2502.03756), Corollary 1.4; the negative-curvature input above is Bryant's 1985 paper, not the 1980 paper on a different immersion problem.

**Proposition 10.3 (near-eigenspace diagnostic; P).** Suppose $E_1$ is absolutely irreducible and the next eigenvalue is $\lambda_{\rm next}>\lambda_1$. If a degree-$d$ map exists, define
$$\zeta_d=\frac{8\pi d/A-\lambda_1}{\lambda_{\rm next}-\lambda_1}.$$
For the normalized frame density above,
$$\frac1A\int_C|F-1|\,dA\le2\zeta_d+2\sqrt{\zeta_d}.$$

*Proof.* Balance the spherical coordinates $x$ and write $x=z+y$, with $z=P_{E_1}x$. The spectral gap gives $\|y\|^2\le A\zeta_d$. Group averaging $|z|^2$ gives $(\|z\|^2/A)F$ by Schur's lemma. Average the identity $|x|^2=1$ and bound the $|y|^2$ and cross terms in $L^1$. Since $\|z\|^2=A-\|y\|^2$, the three errors are at most $\zeta_d,\zeta_d,2\sqrt{\zeta_d}$. $\square$

The historical quadratic-shape estimate with a numerical gap $0.016$ yielded degree $23.7$ and no integer improvement; it remains N. In an orthonormal tangent frame, if an inverse metric is $LL^T$ and the coordinate gradient Gram matrix is $S$, the physical matrix is $L^TSL$. The retired calculation using $LSL^T$ does not contribute to the present certificates. The harmonic argument in Chapter 4 controls physical gradients directly.

## 10.2 Intrinsic bounds for other finite groups

For a nonabelian simple group $G$, let $g_{\min}(G)$ be its minimum faithful genus. The orbit-field argument gives
$$\min\{|G|,\lceil1+\sqrt{g_{\min}(G)}\rceil\}\le\gamma(G)
\le\left\lfloor\frac{g_{\min}(G)+3}{2}\right\rfloor.$$
For the lower bound, a pencil of degree $d<|G|$ generates with its translates a faithful quotient field; Theorem 5.1 bounds that quotient's genus by $(d-1)^2$. The upper bound is the classical general gonality bound on a minimum-genus curve. These inequalities do not assert a square-root asymptotic.

If $\mu(G)<|G|$, Lemma 6.4 instead gives
$$g_{\min}(G)\le\pi(\mu(G),q(G)-1),$$
where $q(G)$ is the minimum nontrivial genuine representation degree. The linearized degree lattice further restricts the possible values. Subgroup induction is equivalently
$$a(G)=\frac{|G|}{\max_{H\le G}|H|/\mu(H)}.$$

The small-group conclusions retained in the project are
$$\gamma(A_5)=1,\quad a(A_5)=2,\quad\gamma(A_6)=5,\quad a(A_6)=12,
\quad\mu(\mathrm{PSL}_2(7))=a(\mathrm{PSL}_2(7))=4.$$
The icosahedral action supplies the first equality, but has no genuine two-dimensional lift, accounting for the difference between $\gamma$ and $a$. The smooth Valentiner plane sextic has genus ten and gonality five, attaining the minimum-genus lower bound for $A_6$. It is distinct from the singular genus-six Wiman sextic. For $A_6$, linearized degrees are divisible by six and $\pi(6,4)=2$, so $\mu(A_6)\ge12$. The index-six $A_5$ gives the upper bound for $a(A_6)$. A subgroup wall below 12 would require an index below 12: orders 40 and 45 force a normal Sylow five-subgroup and are incompatible with its order-ten normalizer; the order-36 subgroup has no Möbius action and therefore moving degree at least two. The Klein quartic and its canonical plane bundle give the degree-four construction for $\mathrm{PSL}_2(7)$; genus three and $\pi(3,2)=1$ give its lower bound.

**Proposition 10.4 (power-sum curves; P).** For $n\ge5$, the equations
$$\sum_{j=1}^n x_j^r=0\qquad(1\le r\le n-2)$$
define a smooth connected complete-intersection curve in $\mathbb P^{n-1}$. The alternating coordinate action is faithful and genuinely linearizes $\mathcal O(1)$, of degree $(n-2)!$. Consequently
$$\mu(A_n),a(A_n)\le(n-2)!,$$
and, for $n\ge7$, induction from $A_7$ gives
$$a(A_n)\le\min\{(n-2)!,n!/84\}.$$

*Proof.* If a solution had at most $n-2$ distinct nonzero coordinate values $a_1,\ldots,a_k$ with multiplicities $m_j$, the first $k$ equations would be a homogeneous linear system in the nonzero integers $m_j$. Its determinant is $\prod_j a_j\prod_{i<j}(a_j-a_i)\ne0$, a contradiction. The Jacobian rows are proportional to $(1),(x_j),\ldots,(x_j^{n-3})$, so their Vandermonde rank is $n-2$. The curve is smooth; connectedness follows from the complete-intersection Koszul resolution. Bézout gives degree $(n-2)!$. A point with the $n$th roots of unity as coordinates lies on it, and a three-cycle moves that point projectively. Simplicity gives faithfulness. The permutation representation supplies the linearization. Finally $[A_n:A_7]\,60=n!/84$. $\square$

For $n=7$ this has degree 120, $K=\mathcal O(8)$, genus 481 and signature $(3,7,7)$. Its upper bound for $\mu(A_7)$ is superseded by 90. The historical exclusion of attempted minima 78, 102 and 114 is also superseded. The first two have no compatible signatures under Castelnuovo and the degree lattice. At 114 the only candidate is $(3,4,5,7)$, genus 1354, whose section space would be the ordinary six. The invariant product of seven coordinates is nonzero of degree 798, whereas no nonnegative combination of orbit sizes $840,630,504,360,2520$ equals 798. This remains a valid supplementary semigroup obstruction, rather than a necessary step in the current minimum proof.

Cyclic groups and dihedral groups $D_{2m}$ with $m$ odd have projectively faithful genuine two-dimensional representations, so $\mu=1$. For the odd dihedral case, rescale the rotation in a two-dimensional lift so its order is $m$; oddness allows the reflection relation and trivial scalar kernel simultaneously. The odd torus normalizers in $\mathrm{PSL}_2(q)$ therefore give
$$a(\mathrm{PSL}_2(q))\le
\begin{cases}
q(q-1)/2,&q\text{ even or }q\equiv1\pmod4,\\
q(q+1)/2,&q\equiv3\pmod4.
\end{cases}$$
For odd $q$ use the normalizer with odd torus order $(q\pm1)/2$; for even $q$ use order $q+1$. These are upper bounds, with no general exact-value assertion.

## 10.3 The historical cyclic-collision model

This reduction preceded ramification transport. Its consequences are now empty for the four curves, but its distinctions between effectivity, norms and actual collisions remain useful. Use the permutations
$$\tau=(12)(34),\ z=(13)(24),\ \gamma=(567),\ \sigma=(34)(57),\ V_4=\langle\tau,z\rangle,$$
$$J=C_G(\tau),\qquad K=(S_4(1234)\times S_3(567))\cap A_7.$$
The notation $E_{22}$ below avoids confusion with the elliptic factors of Chapter 9.

| Curve | Quotient subgroup | Genus |
|---|---|---|
| $D$ | $\langle\tau\rangle$ | 64 |
| $E_{22}$ | $\langle\tau,\gamma\rangle$ | 22 |
| $Q$ | $V_4$ | 28 |
| $F_{10}$ | $V_4\times\langle\gamma\rangle$ | 10 |
| $B_3$ | $J$ | 3 |
| $Y_4$ | $A_4(1234)\times\langle\gamma\rangle$ | 4 |
| $X_0$ | $K$ | 0 |

After normalization, $D\to E_{22}$ and $Q\to F_{10}$ are cartesian étale degree-three covers, with vertical degree-two maps induced by $z$. Write $\beta:E_{22}\to F_{10}$ for the latter. There is a second cartesian square with $\pi:F_{10}\to Y_4$ étale of degree three, $\phi:B_3\to X_0$ trigonal, and double covers $p:F_{10}\to B_3$, $h:Y_4\to X_0$. The cover $F_{10}\to X_0$ is the $S_3$ closure of $\phi$, and $Q\to Y_4$ is étale with group $C_3^2$.

The involutions in $E_{22},F_{10}$ fix 18 and 10 points. A genus-at-least-five trigonal curve has a unique trigonal pencil, so any involution fixes at most six points; a nonhyperelliptic involution on a hyperelliptic curve fixes at most four. These counts imply gonality at least four, so an effective degree-three class has a unique representative.

Let $\eta\in\operatorname{Jac}(E_{22})[3]$ define $D\to E_{22}$, and choose $\xi\in\operatorname{Jac}(F_{10})[3]$ with $\beta^*\xi=\eta$. For a cyclic étale cover of degree $\ell$, a collision $[U-V]=\eta^j\ne0$ pulls back to distinct equivalent degree-$\ell\deg U$ divisors, producing a pencil after cancellation. Since $\operatorname{gon}(D)\ge10$,
$$W_3(E_{22})\cap(W_3(E_{22})+\eta)=\varnothing.$$
In the hypothetical minimal nine-pencil case, the deck $C_3$ preserves the pencil and acts nontrivially on its target. Its two fixed fibres give ordered effective divisors $U,V$ of degree three with $[U-V]=\eta$. The restriction to at most two pencils follows from the birational-pair and genus-four-core alternatives in Theorem 2.4. This supplies the historical converse in that case; it is not a converse for arbitrary étale covers.

Let $\lambda\in\operatorname{Pic}^3(F_{10})$ be the **actual** double-cover root, with $2\lambda=B_F$ and $\beta^*\lambda=R_z^{E_{22}}$. For $M=\mathcal O(U)$, saturation gives
$$\operatorname{Nm}_\beta M\in\{\lambda,\lambda+\xi,\lambda+2\xi\}.$$
Each effective root has one representative and at most eight degree-three lifts, leaving at most 24 candidates. Effectivity of a norm root alone does not establish a collision.

The $V_4$ cover $E_{22}\to B_3$ has labeled branch divisors
$$\deg\Delta_z=3,\quad\deg\Delta_\sigma=9,\quad\deg\Delta_{z\sigma}=1.$$
Put $A=\phi^*\mathcal O_{X_0}(1)$ and $K_{B_3}\sim A+P_0$. For its simply branched fibres $2R_i+T_i$,
$$R_\phi\sim K_{B_3}+2A,\quad\Delta_\sigma+\Delta_{z\sigma}\sim6A-2K_{B_3},
\quad\operatorname{root}(F_{10}/B_3)=3A-K_{B_3}.$$
Distinguished fibres satisfy
$$\phi^{-1}(x_0)=P_1+P_2+P_\tau,\quad\Delta_z=P_1+P_2+P_3,
\quad\phi^{-1}(b_0)=2P_3+P_t,\quad\Delta_{z\sigma}=P_t.$$
The plane-quartic geometry gives $h^0(\Delta_z)=1$ when $P_3\ne P_0$, and two otherwise.

Let $y_\pm$ lie over $x_0$, $y_*$ be ramified over $b_0$, and $\Theta=y_++y_-+y_*$. The actual roots obey both norms
$$\operatorname{Nm}_p(\lambda+k\xi)=\mathcal O_{B_3}(\Delta_z),\qquad
\operatorname{Nm}_\pi(\lambda+k\xi)=\mathcal O_{Y_4}(\Theta),\qquad k=0,1,2.$$
For the first identity, if $M_2$ is the actual root of $E_{22}/\sigma\to B_3$ and $s_t$ lies over $P_t$, then $\lambda=p^*M_2-s_t$ and $2M_2=\Delta_z+P_t$. Coprime degree-two and degree-three eigensheaf multiplication supplies the second identity. Knowing only their squares would leave a two-torsion ambiguity.

In the rigid case $P_3\ne P_0$, the norms leave four lifted triangles; $\lambda$ is ineffective and a collision pairs the two nonzero twists. This reduces to $[U-\sigma U]=\eta$, two tests up to sign. In the moving case $P_3=P_0$, writing $p^{-1}(P_3)=\{s_+,s_-\}$ leaves **three** compatible divisors
$$2s_++s_t,\qquad s_++s_-+s_t,\qquad2s_-+s_t.$$
The middle can represent $\lambda$, while the outer two can represent its nonzero twists. Discarding those two was the historical error. Untwisted effectivity is equivalent to $M_2\sim P_3+P_t$.

The signed collision is
$$\sum_{i=1}^3\epsilon_i[p_i-\sigma p_i]=\pm\eta,$$
in a Prym of dimension 15. If it held, a function with divisor $3U-3\sigma U$ and the equation $y^3=u$ reconstruct the specified cyclic cover; target actions can be written $\gamma y=\zeta_3y$, $zy=-y$, $\sigma y=1/y$. An $H$-preserved nine-pencil imposes the further relation $3P_\tau\sim P_3+2P_0$, a marked flex condition. The smooth $(9,9)$ audit model in Theorem 2.5 shows why the local $C_2\times S_3$ data alone cannot exclude nine-gonality. Full ramification transport excludes the collision on the actual $A_7$ covers without computing these candidates.

## 10.4 Historical accessory transfer and local obstructions

The old transfer principle remains a useful separate mechanism. Suppose an irreducible $G$-equivariant cover $W\dashrightarrow V$ of degree $d$ dominates a faithful curve $T$, and $G$ contains a Möbius subgroup $M$ with $|M|>d$. Interpolate on a free $M$-orbit in $\mathbb P^1$ to construct equivariant rational curves sweeping the linear representation $V$. Pullback components have degrees summing to $d$; each has a nontrivial stabilizer, since a free orbit of components would already cost at least $|M|$. Some component dominates $T$: constant images fixed by nontrivial subgroups form a finite set and cannot account for a sweeping dominant family.

On such a component, choose generic zero and pole fibres of its map to $\mathbb P^1$ whose images in $T$ are disjoint. Their field norm gives a nonconstant function on $T$ of degree at most $d$. Generic choice matters: a norm can otherwise be constant, as $\operatorname{Nm}_{\mathbb C(x)/\mathbb C(x^2)}((x-1)/(x+1))=1$ shows. This yields $\operatorname{gon}(T)\le d$ in this setting, not a statement about every accessory base.

At primes two and three, the Sylow groups of $A_7$ are $D_8$ and $C_3^2$, each of essential dimension two over $\mathbb C$ by the Karpenko--Merkurjev theorem for $p$-groups. Thus a curve compression after a prime-to-$p$ extension is obstructed at either prime, and the total accessory degree must be divisible by six. The old degree-seventeen argument combined this with the transfer principle and the minimum-genus field bound. All its bounded-degree rungs now follow from $a(A_7)=60$. This local divisibility condition gives no iteration theorem for unrestricted towers.

## 10.5 Specialization and the remaining questions

For a semistable degeneration over a discretely valued field, specialize divisors to the **metric** skeleton, with edge lengths given by node thicknesses. Baker's specialization inequality is $r_C(D)\le r_\Gamma(\operatorname{sp}D)$, hence
$$\operatorname{gon}(\Gamma)\le\operatorname{gon}(C).$$
A regular semistable model provides the appropriate subdivided dual graph. The unweighted dual graph of a contracted stable model need not preserve these lengths or give the needed bound. No stable reduction at seven, skeleton or graph-gonality computation for these four covers has been produced by this repository.

The established curve bounds are $25\le\operatorname{gon}(C)\le42$ and $13\le\operatorname{gon}(D)\le21$; the lower bounds retain the analytical certificate. Exact gonality, normalized elliptic and genus-two quotient models, infinite order of the detecting Jacobian points, certified higher equations, and invariants for actual tower bases remain open. The single-accessory thresholds 60 and 90 do not prove $\mathrm{RD}(A_7)\ge2$ or $\mathrm{RD}(S_7)=3$.
