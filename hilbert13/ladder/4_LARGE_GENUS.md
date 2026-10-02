# Chapter 4. Large genus: $\gamma(A_7)\ge25$

$\gamma(G)$ denotes the minimum gonality of a smooth projective curve with a faithful $G$-action.

**Theorem 4.1.** Let $C$ be a smooth projective curve with a faithful $A_7$-action.
1. $\operatorname{gon}(C)\ge23$.
2. If $\operatorname{gon}(C)=d\le24$, then $g(C)\le B^*(d)$, where $B^*(23)=363$, $B^*(24)=397$ and $B^*(d)\le331$ for $d\le22$.
3. **$\gamma(A_7)\ge25$:** every faithful $A_7$-curve has gonality at least 25.

**Inputs.**
- **Base [C].** $\operatorname{gon}\ge25$ for every faithful $A_7$-curve of genus $\le335$ (Chapter 2, §2.5, and Chapter 3).
- **Window [C].** The four rigid signatures with $336\le g\le397$ (§4.4).
- **Algebra [P].** Everything else below is algebra: Castelnuovo–Severi, Castelnuovo's bound, and a singularity count.

Numbers: `gonality_large_genus.py` → `gonality_large_genus_output.txt`.

## 4.1 Conjugate pencils

**Notation.**
- $F=\mathbb C(C)$, and $f:C\to\mathbb P^1$ is a pencil of degree $d=\operatorname{gon}(C)\le24$.
- $K=\mathbb C(f)$, so $[F:K]=d$.
- For $\sigma\in A_7$, the conjugate pencil field is $\sigma K=\mathbb C(f\circ\sigma^{-1})$, and $H=\operatorname{Stab}(K)$.
- Claim($d'$) means: no faithful $A_7$-curve has gonality $\le d'$. We argue by strong induction on $d$.

**Lemma 4.2 (orbit generation) [P].** Assume Claim($\lfloor d/2\rfloor$). Then the compositum $E$ of the fields $\sigma K$ is $F$.

*Proof.*
1. $E$ is $A_7$-stable, so $A_7$ acts on the curve $C_E$ of $E$ with kernel $1$ or $A_7$.
2. If the kernel is $A_7$, then $K\subseteq F^{A_7}$ and $d\ge2520$. So $A_7$ acts faithfully on $C_E$.
3. If $E\ne F$, put $k=[F:E]\ge2$. Then $C_E\to\mathbb P^1$ has degree $d/k\le\lfloor d/2\rfloor$, which contradicts Claim($\lfloor d/2\rfloor$). $\square$

**Lemma 4.3 (orbit size) [P].** $[A_7:H]\ge15$. In particular there is a third conjugate pencil field, distinct from any two given ones.

*Proof.*
- $H=A_7$ would give $E=K\ne F$, which Lemma 4.2 excludes.
- The only proper subgroup of index $<15$ is $A_6$. Suppose $H=A_6$, and let $H_0\trianglelefteq H$ fix $K$ pointwise.
- Then $|H_0|$ divides $d$, so $H_0=1$, since $A_6$ is simple.
- So $A_6$ would embed in $\operatorname{Aut}\mathbb C(t)=\mathrm{PGL}_2(\mathbb C)$, which is impossible. $\square$

**Lemma 4.4 (multilinear independence) [P].** Let $f_1,\dots,f_t\in F$ with $f_i\notin\mathbb C(f_1,\dots,f_{i-1})$. Then the $2^t$ products
$\prod_{i\in I}f_i$ ($I\subseteq\{1,\dots,t\}$) are linearly independent.

*Proof.* Induction on $t$. Write a relation as $Af_t+B=0$, with $A,B$ multilinear in $f_1,\dots,f_{t-1}$.
If $A\ne0$, then $f_t=-B/A$ lies in $\mathbb C(f_1,\dots,f_{t-1})$, a contradiction. So $A=0$, then $B=0$, and induction applies. $\square$

## 4.2 Genus bounds

**Chains.** Start with $f_1=f$ and keep adding a conjugate pencil outside the current compositum. By Lemma 4.2 this reaches $F$ after
$t\le1+\Omega(d)$ pencils, since each step divides $[F:E_i]\mid d$ by at least 2.
- **$t=2$ (a birational pair).** Castelnuovo–Severi with two rational subfields of index $d$ gives $g\le(d-1)^2$.
- **$t\ge3$.** The Segre map of $(f_1,\dots,f_t)$ is birational onto its image in $\mathbb P^{2^t-1}$. The image has degree $td$ and, by Lemma 4.4,
  is nondegenerate. Castelnuovo gives $g\le\pi(td,2^t-1)$, and every case that occurs is $\le\pi(3d,7)$.

A birational pair can still be improved by a third pencil.

**Lemma 4.5 (cost of a dependent third) [P].**
- *Setting.* $(f_1,f_2)$ is a birational pair with image $Y\subset\mathbb P^1\times\mathbb P^1$ of bidegree $(d,d)$, so $\delta(Y)=(d-1)^2-g$.
  $f_3$ is a pencil with $\mathbb C(f_3)\ne\mathbb C(f_1),\mathbb C(f_2)$.
- *Claim.* If the 8 products of $(f_1,f_2,f_3)$ are dependent, then
  $$\delta(Y)\ \ge\ c(d):=\min_{m\ge m',\;m+m'\ge d}\Big[\tbinom m2+\tbinom{m'}2\Big],\qquad c(22,23,24)=110,121,132.$$

*Proof.*
1. **The relation.** A relation is $A(x,y)z_0+B(x,y)z_1$, with $A,B$ of bidegree $(1,1)$, so $f_3=(-B:A)|_Y$.
   - $A,B\not\equiv0$ on $Y$, since $Y$ lies on no $(1,1)$-curve.
   - $A,B$ are not proportional on $Y$, since $f_3$ is nonconstant.
   - They share no component: a common ruling would make $\mathbb C(f_3)$ equal to $\mathbb C(f_1)$ or $\mathbb C(f_2)$.
2. **The fixed part.** On $Y$ the pencil $\langle A,B\rangle$ has degree $2d$, while $f_3$ has degree $d$. So the fixed part has degree $d$,
   supported over the $A\cdot B=2$ base points.
3. **Two base points.** A generic member is transverse to every branch, so the fixed part at $b_i$ is $\operatorname{mult}_{b_i}Y$.
   Hence $m_1+m_2\ge d$, and $\delta\ge\binom{m_1}2+\binom{m_2}2$.
4. **One base point $b$.**
   - Smooth members share a tangent $T$, which is not a ruling direction.
   - A branch $\gamma$ not tangent to $T$ contributes $m_\gamma$.
   - A branch tangent to $T$ contributes $\min(\operatorname{ord}_\gamma A,2m_\gamma)$, where $\operatorname{ord}_\gamma A=m_\gamma+\sum_{k\ge1}m_{b_k}(\gamma)$ over the infinitely near points $b_k$ of $A$.
   - By proximity, a branch through $b_1$ either follows $A$, keeping $m_{b_1}(\gamma)=m_\gamma$, or leaves along the exceptional curve and misses $b_2$.
   - So the fixed part is $\le m+m'$, with $m=\operatorname{mult}_bY$ and $m'=\operatorname{mult}_{b_1}Y\le m$. Again $m+m'\ge d$ and $\delta\ge\binom m2+\binom{m'}2$. $\square$

**Corollary 4.6 [P].** $g\le B^*(d):=\max\big((d-1)^2-c(d),\ \pi(3d,7)\big)$.

*Proof.*
- With a birational pair, Lemma 4.3 supplies a third pencil. Either it is dependent, so $g\le(d-1)^2-c(d)$ by Lemma 4.5, or the triple
  embeds nondegenerately in $\mathbb P^7$ with degree $3d$, so $g\le\pi(3d,7)$.
- Without a birational pair, the chain bound applies. $\square$

| $d$ | 20 | 21 | 22 | 23 | 24 |
|---|---|---|---|---|---|
| $(d-1)^2-c(d)$ | 271 | 300 | 331 | 363 | 397 |
| $\pi(3d,7)$ | 261 | 290 | 320 | 352 | 385 |
| $B^*(d)$ | 271 | 300 | 331 | 363 | 397 |

## 4.3 Proof of Theorem 4.1

1. **Parts 1–2.** Suppose $\operatorname{gon}(C)=d\le22$.
   - If $g\le335$, this contradicts the base input.
   - If $g\ge336$, it contradicts $B^*(d)\le331$.

   The induction in Lemma 4.2 only ever calls Claim($\le12$). The same argument for $d=23,24$ gives $g\le B^*(d)$.
2. **Part 3.** It remains to treat $336\le g\le397$. The faithful $A_7$-curves there (`gonality_large_genus.py`) have the rigid signatures
   - $(3,5,5)$, $g=337$;
   - $(3,4,7)$, $g=346$;
   - $(3,5,6)$ and $(4,4,5)$, $g=379$.

   The listing works as follows:
   - every signature must pass the seven-sheeted Riemann–Hurwitz test;
   - for a triangle signature, an $A_7$-generating triple must exist;
   - four-point signatures start at $g=421$, and quotient genus $\ge1$ at $g=631$.

   §4.4 certifies these four signatures. $\square$

## 4.4 The window signatures [C]

The same certified chain as §2.4 is used (`certify_signatures.py 24 ...` → `certify_window_output.txt`). Every curve passes Li–Yau with
a wide margin.

| signature | $g$ | $48/(g-1)$ | curves (classes) | certified $\lambda_1\ge$ (worst) | Li–Yau gon $\ge$ |
|---|---|---|---|---|---|
| $(3,5,5)$ | 337 | 0.14286 | 2 (6) | 0.24031 | 40.37 |
| $(3,4,7)$ | 346 | 0.13913 | 4 (16) | 0.23617 | 40.74 |
| $(3,5,6)$ | 379 | 0.12698 | 3 (6) | 0.14941 | 28.24 |
| $(4,4,5)$ | 379 | 0.12698 | 7 (22) | 0.26140 | 49.40 |

**Beyond the window** (not needed).
- The same run certified $(3,5,7)$, $(3,6,6)$, $(3,6,7)$, $(3,7,7)$ and $(4,4,6)$, with 24 more curves.
- It stopped at $(4,4,7)$, $g=451$. There the floating dense Cholesky of the $Q_0$ deflation was not positive definite at $n_0=16$.
  That is a failure of the numerics, not a certified negative.

## 4.5 Remarks

- **Cross-check.** $c(d)$ agrees with GPT's Round-3 table (110/121/132). Lemma 4.5 adds the proof in the tangential case, which costs exactly as much as the transversal one.
- **Cited results.**
  - Castelnuovo–Severi (Stichtenoth, Thm 3.11.3);
  - Castelnuovo's bound $\pi(d,r)$;
  - $\delta=\sum\binom{m_p}2$ over infinitely near points, and proximity (Casas-Alvero, *Singularities of Plane Curves*, §3.5);
  - minimum genus 136 and the signature enumeration (exact: `verify_accessory60.py`, `gonality_large_genus.py`).
