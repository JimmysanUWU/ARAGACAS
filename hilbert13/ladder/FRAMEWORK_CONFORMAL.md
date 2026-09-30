# Conformal spectral gonality, and answers to the open questions (paper only, no code)

Status tags as in NOTES.md: [P] proved here on paper, [?] conjecture or expectation, [Q] question for GPT.

## A. The framework: beating the Li–Yau barrier by changing the conformal metric

**The barrier.** Section 7 gives $\operatorname{gon}(C)\ge67.5\,\lambda_1$ for the hyperbolic metric. For triple classes 0, 1 this is
$\le23.37$, so hyperbolic Li–Yau can never prove $\operatorname{gon}(C)\ge25$ or $\operatorname{gon}(C/\langle\tau\rangle)\ge13$. That caps the combined
accessory barrier at degree 23.

**The escape.** The Dirichlet energy is conformally invariant, but $\lambda_1$ and the area are not. Any conformal
metric can be used in Hersch's argument.

**Theorem A [P, classical input].** For any positive bounded function $h$ on $C$ and any holomorphic
$\psi:C\to\mathbb P^1$ of degree $d$,
$$8\pi d\ \ge\ \lambda_1(h\,g_{\rm hyp})\int_C h\,dA .$$
*Proof.* Balance $\psi$ against the measure $h\,dA$ (Hersch; the measure has no atoms). Each coordinate $x_i$ then
satisfies $\int|\nabla x_i|^2\ge\lambda_1(hg)\int h x_i^2$. Sum over $i$, use $\sum x_i^2=1$, and use $\sum\int|\nabla x_i|^2=8\pi d$. $\blacksquare$

So $\operatorname{gon}(C)\ge\Lambda^G/8\pi$, where
$$\Lambda^G:=\sup_{h\ G\text{-invariant}}\lambda_1(hg)\,\mathrm{Area}(hg)\ \le\ \Lambda_1^{\rm conf}(C)\le 8\pi\operatorname{gon}(C).$$
Here $\Lambda_1^{\rm conf}$ is the conformal first eigenvalue (Li–Yau conformal volume; El Soufi–Ilias; Petrides).
Restricting to $G$-invariant $h$ keeps the $A_7$-action isometric. The whole certification of section 7 then
runs unchanged: $h$ only multiplies the mass weight $w_R$.

**The frame function.** Let $E_1$ be the first eigenspace. For the $(2,4,7)$ curves it is the absolutely
irreducible $14_a$. Take an orthonormal basis $\varphi_1,\dots,\varphi_{14}$ and put
$$F=\sum_l\varphi_l^2,\qquad \bar F=\tfrac{A}{14}F\quad(\text{mean }1).$$
$F$ does not depend on the basis, and it is $G$-invariant, so it is a function on the $(2,4,7)$ orbifold.

**Proposition B (first variation) [P].** Let $u$ be $G$-invariant and $g_t=(1+tu)g$. Then
$$\frac{d}{dt}\Big|_{0}\lambda_1(g_t)\,\mathrm{Area}(g_t)=\lambda_1\int_C u\,(1-\bar F)\,dA .$$
*Proof.* A $G$-invariant perturbation keeps $E_1$ an irreducible eigenspace, and the next isotype sits at
$\lambda_2\approx0.57$, far away. So the branch is analytic and stays the lowest for small $t$. Hellmann–Feynman
gives $\lambda'=-\lambda_1\int u\varphi^2$ for unit $\varphi\in E_1$. By Schur, $\varphi\mapsto\int u\varphi^2$ is a multiple of the norm on
$E_1$, namely $\frac1{14}\int uF$. Finally $\mathrm{Area}'=\int u$. $\blacksquare$

**Dichotomy [P].** Exactly one of the following holds.
- **(i) $\bar F\equiv1$.** Then $\Phi=\sqrt{A/14}\,(\varphi_1,\dots,\varphi_{14}):C\to S^{13}$ is a $G$-equivariant harmonic map by first
  eigenfunctions. It is automatically conformal: its Hopf differential is a $G$-invariant holomorphic quadratic
  differential, i.e. one on the rigid orbifold $\mathbb P^1(2,4,7)$, where
  $\dim=-3+\lfloor1\rfloor+\lfloor3/2\rfloor+\lfloor12/7\rfloor=0$. So the hyperbolic metric is conformally $\lambda_1$-critical
  (El Soufi–Ilias), and $C$ carries a canonical minimal immersion into $S^{13}$.
- **(ii) Otherwise.** $u=1-\bar F$ strictly increases $\lambda_1\cdot\mathrm{Area}$, so the Li–Yau bound strictly improves.

Either outcome is interesting. When $\dim E_1=3$ (for example Bolza), (i) is impossible: a conformal eigenmap into
$S^2$ with constant energy density would pull back a curvature $+1$ metric. So for Bolza (ii) holds.

**Proposition C (a priori test) [P].** If $C$ has a map of degree $d$, put
$\zeta_d=(8\pi d/A-\lambda_1)/(\lambda_2-\lambda_1)$. Then
$$\frac1A\int_C|\bar F-1|\ \le\ 2\zeta_d+2\sqrt{\zeta_d}.$$
*Proof sketch.*
1. Energy splitting gives $\|x-P_1x\|^2\le\zeta_dA$.
2. Average $|x\circ g|^2=1$ over $G$. By Schur, $\mathrm{avg}_g|P_1x\circ g|^2=\frac{\|P_1x\|^2}{14}F$.
3. The error is $\mathrm{avg}_g(|z|^2-2x\cdot z)\circ g$, whose $L^1$ norm is at most $A(\zeta+2\sqrt\zeta)$. $\blacksquare$

For $d=24$ we get $\zeta\approx0.041$, so the threshold is about $0.49$: a degree-24 map forces $\bar F$ to within
0.49 of 1 in mean absolute deviation. This is weaker than B combined with certification, but it needs no
optimisation.

**A useful identity [P].** For any holomorphic $x:C\to S^2$, harmonicity gives $\Delta x=e\,x$, with
$e=|\nabla x|^2=2\,\mathrm{Jac}$. So $P_1\big((e-\lambda_1)x\big)=0$ exactly. Also the spectral variance of $x$ equals the
spatial variance of the energy density $e$, which vanishes at the $2d+270$ branch points. This is the handle
for sharper versions of C.

**Payoff [?].**
- If $\Lambda^G>8\pi\cdot24$, i.e. a gain of $+2.7\%$ over the hyperbolic metric for classes 0, 1, then
  $\operatorname{gon}(C)\ge25$ and $\operatorname{gon}(C/\langle\tau\rangle)\ge13$. For $D$ the relevant eigenvalue is the $\tau$-invariant one, which
  is still $14_a$ because $14_a^\tau\ne0$.
- A first-order gain needs only $\bar F$ to be non-constant. The gap $\lambda_2-\lambda_1\approx0.23$ leaves room for a
  deformation of size about 10–30%.
- With analogous bounds for all faithful $A_7$-curves of genus $\le(24-1)^2=529$, the accessory barrier would move
  from 23 to 29.

**Classes 12, 14 are different.** There $\lambda_1\approx0.3597$ ($14_a$), but $\lambda_2\approx0.385$ ($21$) is close, so a
conformal deformation soon swaps the two branches. However, hyperbolic Li–Yau alone already has room:
$67.5\cdot0.3597=24.3$. A sharper certificate there (margin 1.2%, as for the $n=96$ run on classes 0, 1) would give
$\operatorname{gon}\ge25$ directly. **Only classes 0, 1 need the conformal trick.**

**Next computation (when code is allowed).**
1. Compute $\bar F$ on the orbifold, which is cheap from the $Q_1$ eigenvector.
2. Iterate $h\leftarrow h(1+t(1-\bar F_h))$.
3. Certify $\lambda_1(hg)\int h$ with the existing pipeline, using $h$ piecewise constant on the mesh.

**Literature position.** Li–Yau (conformal volume), El Soufi–Ilias (criticality), Petrides and Karpukhin–Stern
(conformal maximisers) are known. Not known to us: using the frame-function criterion under symmetry, plus
certified conformal factors, to get **gonality** lower bounds for specific curves.

## B. Answers to the eight questions

**Q2 (transport beyond the saturated case; partial) [P].** Write $x_V=[R_V]\in\mathrm{Pic}^{54}(C)$ and
$\delta=x_{V_0}-x_{V_1}\in J(C)$.

*Lemma (5-adic separation of Klein classes).* $\delta$ has infinite order, or order divisible by 5.

*Proof.* Suppose $n\delta=0$. Then $nx_{V_0}$ is fixed by $N(V_0)$, and also by $N(V_1)$ because it equals $nx_{V_1}$.
These generate $A_7$. A freely acting $C_5$ then forces $5\mid54n$. $\blacksquare$ This is unconditional.

The transport proof is then:
$$\text{9-pencil}\ \Rightarrow\ 24\delta=0,$$
which contradicts the lemma.

A 10-pencil gives instead
$$24\delta=[\mathcal E_{V_1}]-[\mathcal E_{V_0}],\qquad \mathcal E_V=\sum_{h\in C(\tau)}h^*E^{(h)}.$$
Here each $E^{(h)}$ is a $V$-invariant effective divisor of degree 4. On $D$ it is the residual pair $\{q,\bar\mu q\}$:
a node of the symmetric model $(g,g\circ\bar\mu)$ on the diagonal. So the obstruction becomes an incidence question.
Is the specific point $24\delta$ (infinite order or 5-divisible order) in a difference locus of dimension $\le48$
inside the 136-dimensional $J(C)$? That is arithmetic, not a degree count. **Transport is intrinsically a
saturation-plus-one tool.** The spectral bound already gives $\operatorname{gon}(D)\ge12$.

**Q3 (scope of transport) [P, one hypothesis flagged].** The method gives $\operatorname{gon}(C/\langle\tau\rangle)\ge r/2+1$ when:
- (a) $4\nmid r$;
- (b) the saturated good involutions generate $H'\le C(\tau)/\langle\tau\rangle$ with $|H'|\nmid d$ for all $d<r/2$;
- (c) the normalisers of two Klein groups through $\tau$ generate a subgroup containing a freely acting $C_n$ with
  $n\nmid 3|C(\tau)|r$.

It is never more than one step past saturation.

*Example: PSL(2,13) Hurwitz curves, genus 14.* Here $C(\tau)=D_{12}$, $r=6$, and $N(V_4)=A_4$. Two $A_4$'s through $\tau$ lie
in no maximal subgroup ($D_{14}$, $D_{12}$, $A_4$, $13{:}6$), so they generate the group. Order-13 elements act freely,
and $13\nmid216$.

**So the genus-6 involution quotients of the three genus-14 Hurwitz curves have gonality $\ge4$:** not hyperelliptic,
not trigonal. Spectral methods give nothing here; Kim–Sarnak plus Li–Yau gives $<1$.

The same argument gives $\ge8$ for PSL(2,27) (quotient genus 56) and for PSL(2,29) (quotient genus 70). The
maximal-subgroup facts used are standard, but the $A_4$-generation claim should be double-checked. **The niche is
small, non-arithmetic, low-genus cases,** where Castelnuovo–Severi is nearly sharp.

**Q5 (weak points of our certificate; self-audit).** Ranked by residual risk:
1. The trust that CHOLMOD's floating point obeys the textbook backward-error model.
2. Implementation errors in the ball-arithmetic coefficient bounds. The Bolza test and chart invariance mitigate
   this.
3. The twisted-quotient reduction. It is checked exactly, but it is our own code.

Hersch, the Carstensen–Gedicke–Rim constant, and the lower-bound theorem are verified; the last is in Lean.

**Q6 (is $\operatorname{gon}(C)\ge25$ true?)** Hyperbolic Li–Yau cannot decide it. Framework A is a concrete route; it succeeds iff
$\Lambda^G>8\pi\cdot24$. [?] I expect the true gonality to be well above 25, perhaps 40–56, but there is no argument.

**Q7 (the flex condition).** Their $3P_\tau\sim P_3+2P_0$ on the plane quartic $B$ means the following: $P_\tau$ is a flex, and
its flex tangent meets $B$ again at the point $P$ such that the line $PP_3$ is tangent at $P_0$. This follows from
$|K-P|$ being the $g^1_3$'s. Its truth is independent of $\operatorname{gon}\ge10$, which is now proved, and it needs the
equation of $B$. Open, of curiosity value only.

**Q1, Q4, Q8** are GPT's to answer (below). For Q8 (towers): every invariant here (gonality, $\lambda_1$,
$\Lambda^G$, the Picard degree lattice) belongs to *one* curve, and a tower replaces the curve at each step. I have
no candidate for a tower-stable invariant.

## C. Questions for GPT, prioritised

1. **[Q, critical] Lemma 5.5.** Give the exact Petrakiev statements and ranges at $(18,169)$, especially the
   quadric-lifting step. This is the only gate for $\mathrm{ed}_{\mathbb C}(A_7;\le23)>1$.
2. **[Q] Degree 10 and 11.** Given the 5-adic lemma, is there any balancing mechanism that controls
   $[\mathcal E_{V_1}]-[\mathcal E_{V_0}]$ (degree-4 residuals, i.e. nodes on the diagonal)? Or an argument that $24\delta$ lies
   outside the difference locus?
3. **[Q] Frame function.** Do you see an algebraic reason, from triple products of $14_a$ eigenfunctions or from
   Chevalley–Weil, for $\bar F$ to be constant? Case (i) would mean a canonical minimal immersion $C\to S^{13}$, which
   would be remarkable. Case (ii) breaks the Li–Yau barrier.
4. **[Q] Transport applications.** Is "the genus-6 involution quotients of the genus-14 Hurwitz curves are not
   trigonal" known? Can you confirm the $A_4$-generation in PSL(2,13)?
5. **[Q] Beyond 23.** If $\operatorname{gon}\ge25$ is certified in genus 136, which faithful $A_7$-curves of genus $\le529$ would
   still need $\operatorname{gon}\ge25$ for $\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$? Can your pencil geometry handle them?
6. **[Q] Towers.** Can any transport-type class be defined on the versal torsor over every intermediate field of a
   tower, e.g. on a stack, compatibly with the steps $F_i=F_{i-1}E_i'$?

*Section C is superseded by `QUESTIONS_FOR_GPT.md`, which was rewritten after GPT's proof-chain audit (pinned
8b31a2ac). That audit confirmed Lemma 5.5 in Petrakiev's valid ranges, so $\mathrm{ed}_{\mathbb C}(A_7;\le23)>1$ holds
modulo the declared classical inputs. It also independently rechecked the certificate with exact factor residuals.*

## D. New results after the audit (third session) — exact checks in `frontier_checks.py`

**D1. The $Q_2$ group [C].** $C_{A_7}((16)(23))$ has order 24, centre $\langle(16)(23)\rangle$, and element orders
$1{:}1,\ 2{:}9,\ 3{:}2,\ 4{:}6,\ 6{:}6$. Its Sylow 3-subgroup is normal, and its Sylow 2-subgroup is $D_8$. So it is
$C_3\rtimes D_8\cong D_8\times_{C_2}S_3$, where $D_8$ acts on $C_3$ through its quotient by a Klein four-group. It is not $S_4$.

**D2. Pencil orbits [P].** *For $m<60$, a degree-$m$ pencil on a faithful $A_7$-curve has an $A_7$-orbit of size at
least 35. The orbit has size at least 42 if $3\nmid m$.*

*Proof.* Let $K$ be the stabiliser of the pencil $f$ (that is, $f\circ k=M_k\circ f$), and $K_0=\{k:f\circ k=f\}$. Then $f$ factors
through $C/K_0$, so $|K_0|$ divides $m$. Also $K/K_0\subset PGL_2(\mathbb C)$ is cyclic, dihedral, $A_4$, $S_4$ or $A_5$.

The subgroups of $A_7$ of order $>72$ are $A_7$, $A_6$, $L_2(7)$ and $S_5$:
- $A_7$, $A_6$ and $L_2(7)$ are simple and not polyhedral, so $K_0=K$, which is impossible for $m<60$.
- $S_5$ is not polyhedral, and its other quotients have kernel $\supseteq A_5$.

The only subgroup of order 72 is $(A_4\times3){:}2$. Its polyhedral quotients have kernel $C_3$ (quotient $S_4$), $A_4$
(quotient $S_3$) or $A_4\times3$; each needs $3\mid m$. $\blacksquare$

This replaces the audit's "at least eight orbit members". A birational pair therefore leaves at least 33 thirds, not 6.

**D3. The frame function is never constant [P, classical input: Bryant 1985].** Let $(C,g)$ be a closed hyperbolic
surface, and $G$ a finite group of isometries with $C/G$ a triangle orbifold. Then no $G$-invariant positive semidefinite
form $Q$ on the first eigenspace $E_1$ satisfies $Q(\varphi_1,\dots)\equiv1$. In particular $\bar F\not\equiv1$.

*Proof.* Suppose it does. Diagonalise $Q=\sum_i\psi_i^2$ with $\psi_i\in E_1$. This gives a map $\Phi=(\psi_i):C\to S^N$
by first eigenfunctions with $|\Phi|\equiv1$. Then:
- **Constant energy density.** $0=\Delta|\Phi|^2=2\lambda_1-2|\nabla\Phi|^2$, so $|\nabla\Phi|^2\equiv\lambda_1$.
- **Harmonic.** $\Delta\Phi=|\nabla\Phi|^2\Phi$, so $\Phi$ is harmonic into the sphere.
- **Conformal.** Its Hopf differential $\sum(\partial_z\psi_i)^2dz^2$ is holomorphic, and $G$-invariant because $Q$ is. A triangle orbifold carries
  no nonzero invariant holomorphic quadratic differentials ($-3+\sum\lfloor2(1-1/p_i)\rfloor=0$), so $\Phi$ is weakly conformal.
- **Immersion.** Since $|\nabla\Phi|^2>0$ everywhere, $\Phi$ is a conformal minimal immersion.

The induced metric is $\tfrac{\lambda_1}2g$, which has constant curvature $-2/\lambda_1<0$. But R. Bryant, *Minimal surfaces of
constant curvature in $S^n$* (Trans. AMS 290 (1985) 259–271), shows there are no minimal surfaces of constant negative
curvature in $S^n$, even locally. $\blacksquare$

*Consequence.* This resolves the dichotomy of section A. Case (i) never occurs. The hyperbolic metric is not critical
for $\lambda_1\cdot\mathrm{Area}$ among $G$-invariant conformal metrics:
- if $E_1$ is irreducible, the direction $u=1-\bar F$ has derivative $\lambda_1\int(1-\bar F)^2>0$;
- in general, El Soufi–Ilias duality applies to the $G$-averaged form.

Hence $\Lambda^G>540\pi\,\lambda_1$ strictly, for every $(2,4,7)$ $A_7$-curve. This is only first-order. $\operatorname{gon}(C)\ge25$ on classes 0, 1
needs a gain of at least $192/(540\cdot0.34627)-1\approx2.7\%$, which is a computation.

**D4. The involution fixed-point classes are governed by two points [P][C].**

*The fixed-point module.* $A_7$ has one class of involutions, with 105 elements and centraliser $C(\tau)$ of order 24.
$v\mapsto\mathrm{Fix}(v)$ is $A_7$-equivariant. So the $\mathbb Q[A_7]$-span of the fixed divisors is a quotient of
$$\mathbb Q[A_7/C(\tau)]=1+6+14_a+2\cdot14_b+21+35 .$$
Of these, only 21 and 35 occur in $H^1(C)$, since $\mathrm{Jac}(C)\sim A^{10}E_1^{15}E_2^{21}S^{35}$ (section 5 of NOTES).

*Abel–Jacobi.* The map $\mathrm{AJ}\otimes\mathbb Q$ on degree-0 fixed-point divisors is therefore a pair of $G$-maps
$\varphi_{21}:21\to J\otimes\mathbb Q$ and $\varphi_{35}:35\to J\otimes\mathbb Q$. Each is zero or injective, since the target copies are
irreducible. Components in $1,6,14_a,14_b$ are killed by the integral idempotents $2520e_\rho$, so they are torsion.

*The two divisors [C].* Take the Klein difference
$D_\delta=\mathrm{Fix}(v_2)+\mathrm{Fix}(v_3)-\mathrm{Fix}(w_2)-\mathrm{Fix}(w_3)$, with $V_0=\langle(12)(34),(13)(24)\rangle$ and
$V_1=\langle(12)(34),(12)(56)\rangle$. Its squared norms by isotype are $6{:}7.2,\ 14_a{:}12,\ 14_b{:}28.8,\ 21{:}18,\ 35{:}6$.
Both $D_\delta$ and the lift $\Phi_g-3\Phi_c$ of $(\star)$ have nonzero 21- and 35-components, and each generates the unique copy.

*Corollary.*
$$\delta\text{ has infinite order}\iff(\varphi_{21},\varphi_{35})\ne0\iff[\Phi_g-3\Phi_c]\text{ has infinite order}\ \Longrightarrow\ (\star)\ \Longrightarrow\ \operatorname{gon}(D)\ge10 .$$
The 5-adic lemma (section B) still gives $\delta\ne0$ unconditionally. If $\varphi=0$, then $\delta$ is torsion of order divisible by 5.

*Explicit carriers [C].*
- **An elliptic curve.** Let $L_2(5)\cong A_5$ act transitively on 6 letters, inside a point-stabiliser $A_6$.
  $E:=C/L_2(5)$ has genus 1; its fixed irreducibles are $1,6,14_a,21$, so $\mathrm{Jac}(E)\sim E_2$. All 42 conjugates see the
  21-part of $D_\delta$. Hence $\varphi_{21}\ne0$ iff $P_E:=[\pi_*D_\delta]\in E$ has infinite order. For the "if" direction,
  $\pi^*P_E=|K|\,e_K\delta$, and $e_K$ keeps only the 21-part modulo torsion.
- **A genus-2 curve.** $C/(3^2{:}4)$, where $3^2{:}4$ is the $Q_1$ group, has genus 2 and Jacobian $\sim S$. Its fixed irreducibles
  are $1,6,14_a,14_b,35$. 58 of 70 conjugates see the 35-part, giving a point $P_S$ in the same way.

$E$ is a 6-sheeted cover of $C/A_6\cong\mathbb P^1$, whose map to $C/A_7$ is the degree-7 Belyi map with passport
$[2^21^3,\ 4\,2\,1,\ 7]$. So deciding whether $P_E$ is torsion is a finite, explicit computation once a model is known.

*Lens [?].* $\Delta(2,4,7)$ is arithmetic (Takeuchi; invariant trace field $\mathbb Q(\cos\pi/7)$). Congruence quotients of an
arithmetic triangle group are of $PSL_2$/$PGL_2$ type over finite fields, and $A_7$ is not. So $\ker(\Delta\to A_7)$ is
non-congruence, and $\delta$ is a difference of elliptic (CM-type) points. The question resembles Manin–Drinfeld in a
setting where no Hecke argument is available.

**D5. The next accessory frontier: $\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$ [P-sketch, uses the audit's inputs].**

*Reduction.* Follow the audit's page-8 logic.
- The prime-local inputs leave degrees $6,12,18,24$.
- No subgroup index of $A_7$ divides 24.
- The $A_5$ sweep needs $d<60$.

So it suffices that every faithful $A_7$-curve has $\operatorname{gon}\ge25$. Only $g\le(24-1)^2=529$ matters, by the orbit-field argument:
index $\ge2$ would give a faithful curve of genus $\le121$.

*Claim: every faithful $A_7$-curve with $g\ge336$ has $\operatorname{gon}\ge25$.* A degree-$m$ pencil needs $(m-1)^2\ge g$, so
$m\in\{20,\dots,24\}$. The audit's inputs give, for each $m$:

| $m$ | pair bound $\max_e[e(m/e-1)^2+(e-1)(m-1)]$ | $\pi_1(3m,7)$ (Segre gap + Eisenbud–Harris) | budget $(m-1)^2-336$ | min cost of 1 / 2 / 3 dependent thirds |
|---|---|---|---|---|
| 20 | 181 | 228 | 25 | 90 / 120 / 135 |
| 21 | 148 | 253 | 64 | 100 / 133 / 150 |
| 22 | 221 | 279 | 105 | 110 / 147 / 165 |
| 23 | — | 306 | 148 | 121 / 161 / 181 |
| 24 | 265 | 335 | 193 | 132 / 176 / 198 |

*Proof of the claim.*
1. Two distinct orbit pencils are birational, because the pair bound is $<336$.
2. Every independent third is excluded, because $\pi_1(3m,7)<336$.
3. At most two thirds can be dependent: by the table, three dependent thirds cost more than the budget.
4. But D2 supplies at least 33 thirds (the audit's 6 would also suffice). $\blacksquare$

The cost model is the audit's: a dependent third has centres of multiplicities $r+s=m$ costing $\binom r2+\binom s2$. Distinct
thirds have distinct centre pairs, because a pencil of $(1,1)$-forms is fixed by its two base points.

*What remains [C].* All 4-point signatures have $g\ge421$, since the seven-sheeted Riemann–Hurwitz test needs
$\sum\mathrm{ind}\ge12$, so they are covered. The remaining curves are the rigid 3-point signatures with $136\le g\le335$. Each
entry gives [ordered generating triples up to conjugation] and the Yang–Yau threshold $\lambda_1>48/(g-1)$ for
$\operatorname{gon}\ge25$:
- $(2,4,7)$, $g=136$: [4], 0.3556
- $(3,3,5)$, $g=169$: [2], 0.2857
- $(2,5,7)$, $g=199$: [4], 0.2424
- $(3,3,6)$ [2] and $(3,4,4)$ [8], $g=211$: 0.2286
- $(2,6,7)$ [4] and $(3,3,7)$ [4], $g=241$: 0.2000
- $(2,7,7)$, $g=271$: [6], 0.1778
- $(3,4,5)$, $g=274$: [10], 0.1758
- $(3,4,6)$ [6] and $(4,4,4)$ [24], $g=316$: 0.1524

For $g\le289$ this list reproduces the audit's page-7 table exactly.

*The program.* $\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$ reduces to $\operatorname{gon}\ge25$ on these 11 signatures. For $(2,4,7)$:
- classes 12 and 14 need a finer hyperbolic certificate, since $\lambda_1\approx0.3597>0.3556$;
- classes 0 and 1 need the D3 conformal gain.

For the other ten signatures, the options are a spectral certificate on the $(p,q,r)$ tiling, or refined pencil
geometry ($\pi_2$ / Petrakiev), which may reach below 336.
