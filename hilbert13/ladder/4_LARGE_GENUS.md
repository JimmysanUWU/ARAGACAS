# Chapter 4. Large genus: $\gamma(A_7)\ge25$

$\gamma(G)$ is the least gonality of a smooth projective curve with a faithful $G$-action.

**Idea.**
- **Many conjugates.** A pencil of degree $d\le24$ on a faithful $A_7$-curve has at least 15 conjugate pencils, and together they generate the function field.
- **Two pencils** give a model in $\mathbb P^1\times\mathbb P^1$.
- **A third pencil** has two options:
  - it is independent, which gives a birational model in $\mathbb P^7$ with degree $3d$ and bounds the genus by Castelnuovo;
  - it depends on the first two, which costs singularities on the $\mathbb P^1\times\mathbb P^1$ model.
- **Conclusion.** Either way $g\le B^\*(d)\le397$. Every faithful $A_7$-curve of that genus is covered by the certificates of Chapters 2–3.

**Theorem 4.1.** Let $C$ be a faithful $A_7$-curve.
1. $\operatorname{gon}(C)\ge23$, using only the first input below (not the window).
2. If $\operatorname{gon}(C)=d\le24$, then $g(C)\le B^*(d)$, where $B^*(23)=363$, $B^*(24)=397$, and $B^*(d)\le331$ for $d\le22$.
3. **$\gamma(A_7)\ge25$.**

**Inputs.**
- $\operatorname{gon}\ge25$ for all faithful $A_7$-curves of genus $\le335$ (Chapters 2–3) [C].
- The four window signatures, $336\le g\le397$: either the certificates of §4.4 [C], or Theorem 4.9 [P].
- Everything else is algebra [P].

Numbers: `gonality_large_genus.py`.

## 4.1 Conjugate pencils

Let $F=\mathbb C(C)$ and let $f$ be a pencil of degree $d=\operatorname{gon}(C)\le24$. Put $K=\mathbb C(f)$, $\sigma K=\mathbb C(f\circ\sigma^{-1})$ and $H=\mathrm{Stab}(K)$. Claim($d'$) means that no faithful $A_7$-curve has gonality $\le d'$; we argue by strong induction.

**Lemma 4.2 (generation) [P].** Under Claim($\lfloor d/2\rfloor$), the compositum $E$ of the $\sigma K$ equals $F$.

*Proof.* $A_7$ acts on $C_E$, with kernel $1$ or $A_7$. Kernel $A_7$ would force $d\ge2520$. If $[F:E]=k\ge2$, then $C_E$ is faithful with a pencil of degree $d/k\le\lfloor d/2\rfloor$. $\square$

**Lemma 4.3 (orbit size) [P].** $[A_7:H]\ge15$. In particular, there is a third conjugate pencil field.

*Proof.* $H=A_7$ contradicts Lemma 4.2. The only proper subgroup of index $<15$ is $A_6$. Its subgroup fixing $K$ has order dividing $d$, so it is trivial, and then $A_6\hookrightarrow\mathrm{PGL}_2(\mathbb C)$, which is impossible. $\square$

**Lemma 4.4 (multilinear independence) [P].** If $f_i\notin\mathbb C(f_1,\dots,f_{i-1})$ for each $i$, then the $2^t$ products $\prod_{i\in I}f_i$ are linearly independent. (Induct on $t$: a relation $Af_t+B=0$ with $A\ne0$ would put $f_t$ in $\mathbb C(f_1,\dots,f_{t-1})$.)

## 4.2 Genus bounds

**Chains.** Add conjugate pencils outside the current compositum. By Lemma 4.2 this reaches $F$ after $t\le1+\Omega(d)$ steps.
- **$t=2$, a birational pair.** Castelnuovo–Severi gives $g\le(d-1)^2$.
- **$t\ge3$.** The Segre image in $\mathbb P^{2^t-1}$ has degree $td$ and is nondegenerate by Lemma 4.4. So $g\le\pi(td,2^t-1)\le\pi(3d,7)$ in every case that occurs.

**Lemma 4.5 (cost of a dependent third) [P].** Let $(f_1,f_2)$ be a birational pair with image $Y\subset\mathbb P^1\times\mathbb P^1$ of bidegree $(d,d)$, so $\delta(Y)=(d-1)^2-g$. If the 8 products of $(f_1,f_2,f_3)$ are dependent, then
$$\delta(Y)\ge c(d):=\min_{m\ge m',\,m+m'\ge d}\Big[\tbinom m2+\tbinom{m'}2\Big],\qquad c(22,23,24)=110,121,132.$$

*Proof.*
1. **The relation** is $Az_0+Bz_1$, with $A,B$ of bidegree $(1,1)$, so $f_3=(-B:A)|_Y$.
   - Neither vanishes on $Y$, which lies on no $(1,1)$-curve.
   - They are not proportional on $Y$, since $f_3$ is nonconstant.
   - They share no ruling, since $\mathbb C(f_3)\ne\mathbb C(f_i)$.
2. **The fixed part.** $\langle A,B\rangle$ has degree $2d$ on $Y$ and $f_3$ has degree $d$. So the fixed part has degree $d$, over the $A\cdot B=2$ base points.
3. **Two base points.** The fixed part at $b_i$ is $\mathrm{mult}_{b_i}Y$. So $m_1+m_2\ge d$ and $\delta\ge\binom{m_1}2+\binom{m_2}2$.
4. **One base point $b$, where $A$ and $B$ are tangent along $T$.**
   - Every member of $\langle A,B\rangle$ passes through $b$ and through the infinitely near point $b_1$ in direction $T$. It passes through nothing further, since $A\cdot B=2$.
   - By Noether's formula, a general member meets a branch $\gamma$ of $Y$ at $b$ with multiplicity $m_\gamma+m'_\gamma$. Here $m'_\gamma$ is the multiplicity at $b_1$ of the strict transform of $\gamma$, which is 0 if $\gamma$ is transverse to $T$.
   - Summing over branches, the fixed part is at most $m+m'$, where $m=\mathrm{mult}_bY\ge m'=\mathrm{mult}_{b_1}Y$. Moreover $\delta\ge\binom m2+\binom{m'}2$. $\square$

**Corollary 4.6 [P].** $g\le B^*(d):=\max\big((d-1)^2-c(d),\ \pi(3d,7)\big)$. With a birational pair, Lemma 4.3 supplies a third pencil, which is either dependent or gives a birational model in $\mathbb P^7$; otherwise the chain bound applies.

| $d$ | 20 | 21 | 22 | 23 | 24 |
|---|---|---|---|---|---|
| $(d-1)^2-c(d)$ | 271 | 300 | 331 | 363 | 397 |
| $\pi(3d,7)$ | 261 | 290 | 320 | 352 | 385 |

## 4.3 Proof of Theorem 4.1

1. **$d\le22$.** If $g\le335$ this contradicts the base input; if $g\ge336$ it contradicts $B^*\le331$. The induction only calls Claim($\le12$). The same argument gives $g\le B^*(d)$ for $d=23,24$.
2. **$336\le g\le397$.** The faithful curves there have rigid signatures $(3,5,5)$ ($g=337$), $(3,4,7)$ (346), $(3,5,6)$ and $(4,4,5)$ (379). The listing uses three facts:
   - every signature passes the seven-sheeted Riemann–Hurwitz test;
   - a triangle signature needs a generating triple;
   - four branch points give $g\ge421$, and quotient genus $\ge1$ gives $g\ge631$. $\square$

## 4.4 The window signatures [C]

The chain of §2.4 is used (`certify_signatures.py` → `certify_window_output.txt`).

| signature | $g$ | $48/(g-1)$ | curves (classes) | $\lambda_1\ge$ (worst) | Li–Yau $\operatorname{gon}\ge$ |
|---|---|---|---|---|---|
| $(3,5,5)$ | 337 | 0.14286 | 2 (6) | 0.24031 | 40.37 |
| $(3,4,7)$ | 346 | 0.13913 | 4 (16) | 0.23617 | 40.74 |
| $(3,5,6)$ | 379 | 0.12698 | 3 (6) | 0.14941 | 28.24 |
| $(4,4,5)$ | 379 | 0.12698 | 7 (22) | 0.26140 | 49.40 |

Beyond the window, not needed: $(3,5,7)$, $(3,6,6)$, $(3,6,7)$, $(3,7,7)$, $(4,4,6)$ also certify (24 curves). The run stopped at $(4,4,7)$, where the dense $Q_0$ deflation failed numerically. That is not a certified negative.

In Lemma 4.5 the tangential case costs as much as the transversal one.

## 4.5 An algebraic route for $g\ge266$ [P]

This section is due to GPT (DAY 2 §5). The full proof is in `gpt/audit_3oct/MULTIPLE_PENCILS.md`; it is reviewed here, and its finite arithmetic re-checked exactly. It replaces $\pi(3d,7)$ by a sharp bound and gives a second, purely algebraic proof for the window of §4.4.

**Lemma 4.7 (multigraded genus bound).** Let $B$ be an integral curve on an integral surface $S\subset(\mathbb P^1)^3$ of type $(a,b,c)$, with $B\not\subset\mathrm{Sing}\,S$. Let $v=(m_1,m_2,m_3)$ be its coordinate degrees and $M$ the matrix with rows $(0,c,b),(c,0,a),(b,a,0)$. Then
$$g(B)\le1+\tfrac12\Big(v^TM^{-1}v+\textstyle\sum_i(a_i-2)m_i\Big),$$
where $(a_1,a_2,a_3)=(a,b,c)$. For equal degrees $m$ and types $(2,1,1)$, $(2,2,1)$, $(2,2,2)$ this reads $A(m)=m^2/2-m+1$, $7m^2/16-m/2+1$ and $3m^2/8+1$.

*Proof.*
1. On a minimal resolution $S'$ of the normalisation, the pullbacks $H_i$ have intersection matrix $M$, and $\sum H_i$ is nef and big. The Hodge index theorem gives $B'^2\le v^TM^{-1}v$.
2. Write $K_{S'}=f^\*K_S+Z$ with $K_S=(a-2,b-2,c-2)|_S$. The conductor and the minimality of the resolution give $Z\le0$, and $B'$ is not in its support.
3. Adjunction gives $2g-2\le B'^2+K_{S'}\cdot B'\le B'^2+f^\*K_S\cdot B'$.

Types $(a,b,1)$ are graphs of rational maps $\mathbb P^1\times\mathbb P^1\dashrightarrow\mathbb P^1$. Resolving the base cluster gives the same bounds, also when $B\subset\mathrm{Sing}\,S$; a $(2,2,2)$ surface singular along $B$ reduces to this case through a partial derivative. $\square$

**Theorem 4.8 (three pencils, sharp).** Let $C$ carry three pencils of degree $m\ge8$, pairwise birational, whose eight products are linearly independent. Then $g\le\lfloor A(m)\rfloor$. Equality holds for every even $m$: take a smooth curve in $\lvert-\tfrac m2K_S\rvert$ on a smooth $(2,1,1)$ surface $S$, a del Pezzo surface of degree 4.

*Proof.* Suppose $g>A(m)$.
1. **The product model is linearly normal.** It is a nondegenerate curve of degree $3m$ in $\mathbb P^7$. A ninth section would give a model in $\mathbb P^8$. Since $g>\pi_2(3m,8)$, Petrakiev's Theorem 2.16(a) puts that model on a surface of degree $\le8$, which must lie in the Segre threefold. A divisor of type $(1,1,1)$ contradicts independence, and types with a zero entry contradict birationality. That leaves $(2,1,1)$, which gives $g\le A(m)$ by Lemma 4.7.
2. **Quadrics.** Let $\Gamma$ be a general hyperplane section. If $h_\Gamma(2)\le18$, linear normality gives $h_B(2)=8+h_\Gamma(2)\le26<27$. So $B$ lies on a $(2,2,2)$ hypersurface, and Lemma 4.7 gives $g\le A(m)$ for $m\ge8$.
3. **Uniform position.** Otherwise $h_\Gamma(2)\ge19$. Subadditivity gives $h_\Gamma(2j)\ge\min(3m,18j+1)$ and $h_\Gamma(2j+1)\ge\min(3m,18j+7)$. Castelnuovo's sum $\sum_{i\ge1}(3m-h_\Gamma(i))$ is then $\le\lfloor A(m)\rfloor$, checked by six residue polynomials in $m$ mod 6. $\square$

**Theorem 4.9 [P].** A faithful $A_7$-curve of genus $\ge266$ has gonality $\ge25$.

*Proof.* Let $d=\operatorname{gon}(C)\le24$.
1. **Birational pairs.** Faithful $A_7$-curves have gonality $\ge13$: otherwise a pencil's orbit generates the field of a faithful curve of genus $\le(12-1)^2<136$. So every orbit compositum is the whole field. Two conjugate pencils are birational, since otherwise $g\le d^2/e+ed-3d+1\le265$ for a proper divisor $e\mid d$. Castelnuovo–Severi then gives $d\ge18$.
2. **Dependent thirds.** An independent third pencil would give $g\le A(24)=265$ (Theorem 4.8). So every third is dependent, and by Lemma 4.5 it is a $(1,1)$ pencil with a length-two base cluster on the $(d,d)$ model. Distinct pencil fields have distinct clusters, and Lemma 6.1 gives at least 33 of them.
3. **A star.** Two disjoint clusters cost more than $\delta=(d-1)^2-g\le263$. So the clusters pairwise meet, hence form a star with centre multiplicity $r$. The bound $\binom r2+33\binom{d-r}2\le\delta$ leaves only $d=24$, $r=23$.
4. **A plane curve.** Projecting the quadric from the (proper) centre gives a birational plane model of degree $48-23=25$. A singular point would give a pencil of degree $\le23$, so $C$ is a smooth plane curve of genus 276. But every faithful $A_7$-curve has $g\equiv1\pmod3$: in Riemann–Hurwitz, $3\mid2520/e$ for every $e\le7$. Contradiction. $\square$

With Chapters 2–3 for $g\le335$, Theorem 4.9 proves $\gamma(A_7)\ge25$ without §4.4.
