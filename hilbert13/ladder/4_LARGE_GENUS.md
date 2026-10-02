# Chapter 4. Large genus: $\gamma(A_7)\ge25$

$\gamma(G)$ is the least gonality of a smooth projective curve with a faithful $G$-action.

**Theorem 4.1.** Let $C$ be a faithful $A_7$-curve.
1. $\operatorname{gon}(C)\ge23$.
2. If $\operatorname{gon}(C)=d\le24$, then $g(C)\le B^*(d)$, where $B^*(23)=363$, $B^*(24)=397$, and $B^*(d)\le331$ for $d\le22$.
3. **$\gamma(A_7)\ge25$.**

**Inputs.**
- $\operatorname{gon}\ge25$ for all faithful $A_7$-curves of genus $\le335$ (Chapters 2–3) [C].
- The four window signatures, $336\le g\le397$ (§4.4) [C].
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
4. **One base point $b$, with common tangent $T$.**
   - A branch transverse to $T$ contributes $m_\gamma$.
   - A tangent branch contributes $\min(\operatorname{ord}_\gamma A,2m_\gamma)$.
   - By proximity, a branch through $b_1$ follows $A$ or misses $b_2$. So the fixed part is $\le m+m'$, with $m=\mathrm{mult}_bY$ and $m'=\mathrm{mult}_{b_1}Y$. $\square$

**Corollary 4.6 [P].** $g\le B^*(d):=\max\big((d-1)^2-c(d),\ \pi(3d,7)\big)$. With a birational pair, Lemma 4.3 supplies a third pencil, which is either dependent or embeds the curve in $\mathbb P^7$; otherwise the chain bound applies.

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

$c(d)$ agrees with GPT's Round-3 table. Lemma 4.5 adds the tangential case, which costs as much as the transversal one.
