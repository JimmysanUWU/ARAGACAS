# Gonality of faithful $A_7$-curves in large genus: an in-repo proof

**Purpose.** This replaces GPT's cited "$g\ge336$ algebra" (Round 3/DAY 2) with a self-contained argument. The tools are
Castelnuovo–Severi, Castelnuovo's bound, and an elementary singularity count. Numbers: `gonality_large_genus.py` →
`gonality_large_genus_output.txt`. Window certificates: `certify_window_output.txt`.

**Theorem.** Let $C$ be a smooth projective curve with a faithful $A_7$-action.
1. $\operatorname{gon}(C)\ge23$, always.
2. If $\operatorname{gon}(C)=d\le24$ then $g(C)\le B^*(d)$, where $B^*(23)=363$ and $B^*(24)=397$.
3. Hence $\operatorname{gon}(C)\ge25$ unless $336\le g\le397$. In that range the faithful $A_7$-curves are those of the rigid signatures
   $(3,5,5)$ ($g=337$), $(3,4,7)$ ($346$), $(3,5,6)$ and $(4,4,5)$ ($379$). These are certified spectrally (§4).
   **So $\gamma(A_7)\ge25$:** every faithful $A_7$-curve has gonality at least 25.

**Base input (certified).** $\operatorname{gon}\ge25$ for every faithful $A_7$-curve of genus $\le335$ (`VERIFICATION.md` §§6–7 and
`HH_CERTIFICATE.md`). These are exactly the eleven rigid signatures. Four-point signatures have $g\ge421$, by the seven-sheeted
test, and $h\ge1$ gives $g\ge631$.

## 1. Set-up

- $F=\mathbb C(C)$, and $f:C\to\mathbb P^1$ a pencil of degree $d=\operatorname{gon}(C)\le24$.
- $K=\mathbb C(f)$, with $[F:K]=d$.
- For $\sigma\in A_7$, the conjugate pencil field is $\sigma K=\mathbb C(f\circ\sigma^{-1})$.
- $H=\operatorname{Stab}(K)$, so the orbit has $[A_7:H]$ elements.

We prove "no faithful $A_7$-curve has gonality $d$" by strong induction on $d$. Write Claim($d'$) for "gonality is never $\le d'$".

**Lemma 1 (orbit generation).** Assume Claim($\lfloor d/2\rfloor$). Then the compositum $E$ of all $\sigma K$ is $F$.

*Proof.*
- $E$ is $A_7$-stable. The kernel of $A_7\to\operatorname{Aut}E$ is $1$ or $A_7$.
- If it is $A_7$, then $K\subseteq F^{A_7}$, so $d\ge2520$. Impossible.
- So $A_7$ acts faithfully on the curve $C_E$.
- If $E\ne F$, put $k=[F:E]\ge2$. Then $k\mid d$, and $C_E\to\mathbb P^1$ (from $K\subseteq E$) has degree $d/k\le\lfloor d/2\rfloor$. This contradicts the claim for $C_E$. $\square$

**Lemma 2 (multilinear independence).** Let $f_1,\dots,f_t\in F$ be nonconstant, with $f_i\notin\mathbb C(f_1,\dots,f_{i-1})$. Then the $2^t$ products
$\prod_{i\in I}f_i$, for $I\subseteq\{1..t\}$, are linearly independent.

*Proof.* Induction on $t$. Write a relation as $A f_t+B=0$, with $A,B$ multilinear in $f_1..f_{t-1}$.
- If $A\ne0$, then $f_t=-B/A\in\mathbb C(f_1..f_{t-1})$, a contradiction.
- So $A=0$, then $B=0$, and induction applies. $\square$

**Lemma 3 (chains).**
- Start with $f_1=f$. Repeatedly add a conjugate pencil not contained in the current compositum.
- By Lemma 1 this ends at $F$ after $t\le1+\Omega(d)$ pencils. Each step divides $[F:E_i]$, a divisor of $d$, by at least 2.
- **If $t=2$** (a *birational pair*): Castelnuovo–Severi with two rational subfields of index $d$ gives $g\le(d-1)^2$.
- **If $t\ge3$:** the Segre map of $(f_1..f_t)$ is birational onto its image in $\mathbb P^{2^t-1}$. By Lemma 2 the image is nondegenerate, and it has degree $td$.
  Castelnuovo gives $g\le\pi(td,2^t-1)$, which is at most $\pi(3d,7)$ in every case occurring.

**Lemma 4 (cost of a dependent third).**
- *Setting.* $(f_1,f_2)$ is a birational pair with image $Y\subset\mathbb P^1\times\mathbb P^1$ of bidegree $(d,d)$, so $\delta(Y)=(d-1)^2-g$.
  $f_3$ is a conjugate pencil with $\mathbb C(f_3)\ne\mathbb C(f_1),\mathbb C(f_2)$.
- *Claim.* If the 8 products of $(f_1,f_2,f_3)$ are dependent, then
  $$\delta(Y)\ \ge\ c(d):=\min_{m\ge m',\,m+m'\ge d}\Big[\tbinom m2+\tbinom{m'}2\Big]\qquad(c(24)=132,\ c(23)=121,\ c(22)=110).$$

*Proof.*
1. **The relation.** A relation is a trilinear $P=A(x,y)z_0+B(x,y)z_1$ with $A,B$ forms of bidegree $(1,1)$, and $f_3=(-B:A)|_Y$.
   - $A,B\not\equiv0$ on $Y$: otherwise $Y$ would lie on a $(1,1)$-curve.
   - $A,B$ are not proportional on $Y$: otherwise $f_3$ is constant.
   - They have no common component: a common ruling would make $f_3$ factor through one projection, so $\mathbb C(f_3)\in\{\mathbb C(f_1),\mathbb C(f_2)\}$.
2. **The fixed part.** The pencil $\langle A,B\rangle$ restricted to $Y$ has degree $2d$, and $f_3$ has degree $d$. So the fixed part has degree $d$,
   supported over the base points ($A\cdot B=2$).
3. **Transversal case** (two distinct base points). A generic member is transverse to every branch, so the fixed part at $b_i$ equals
   $\operatorname{mult}_{b_i}Y$. Hence $m_1+m_2\ge d$, and $\delta\ge\binom{m_1}2+\binom{m_2}2$.
4. **Tangential case** (one base point $b$).
   - *The pencil at $b$.* Smooth members share a tangent $T$, which is not a ruling direction. One member is $D_0=\ell\cup\ell'$, the two rulings through $b$.
   - *Branches not tangent to $T$* contribute $m_\gamma$.
   - *Branches tangent to $T$* contribute $\min(\operatorname{ord}_\gamma A,2m_\gamma)$. Moreover $\operatorname{ord}_\gamma A=m_\gamma+\sum_{k\ge1}m_{b_k}(\gamma)$ along the infinitely near points $b_k$ of $A$.
   - *Proximity.* A branch through $b_1$ that continues along $A$ keeps multiplicity $m_{b_1}(\gamma)=m_\gamma$. One that turns along the exceptional curve does not reach $b_2$.
     So the excess over $m_\gamma$ is at most $m_{b_1}(\gamma)$.
   - *Count.* The fixed part is at most $m+m'$, with $m'=\operatorname{mult}_{b_1}Y\le m$. Hence $m+m'\ge d$ and $\delta\ge\binom m2+\binom{m'}2$. $\square$

**Lemma 5 (orbit size).** $[A_7:H]\ge15$.

*Proof.*
- $H=A_7$ is excluded by Lemma 1, since $E=K\ne F$.
- The only proper subgroup of index $<15$ is $A_6$. Suppose $H=A_6$, and let $H_0$ be the pointwise stabiliser of $K$. Then $H_0\trianglelefteq H$.
- $|H_0|$ divides $d$, because $K\subseteq F^{H_0}$. So $H_0=1$, since $A_6$ is simple.
- Then $A_6$ would embed in $\operatorname{Aut}\mathbb C(t)=\mathrm{PGL}_2(\mathbb C)$, which is impossible.

In particular, a third conjugate pencil field, distinct from $\mathbb C(f_1)$ and $\mathbb C(f_2)$, always exists. $\square$

**Corollary (pair case).** Either $g\le(d-1)^2-c(d)$, or every third conjugate pencil is independent. In the latter case $(f_1,f_2,f_3)$
embeds birationally and nondegenerately in $\mathbb P^7$ with degree $3d$, so $g\le\pi(3d,7)$. Thus $g\le B^*(d)=\max\big((d-1)^2-c(d),\pi(3d,7),\text{chain bounds}\big)$.

## 2. The numbers (`gonality_large_genus_output.txt`)

| $d$ | 20 | 21 | 22 | 23 | 24 |
|---|---|---|---|---|---|
| $(d-1)^2-c(d)$ | 271 | 300 | 331 | 363 | 397 |
| $\pi(3d,7)$ | 261 | 290 | 320 | 352 | 385 |
| $B^*(d)$ | 271 | 300 | 331 | 363 | 397 |

For $d\le22$, $B^*(d)\le331<336$. With the base input, this gives Claim(22) for all genera, and hence (1). For smaller $d$ the
induction in Lemma 1 needs only Claim($\le12$).

## 3. Proof of the Theorem

**Part 1.** If $\operatorname{gon}(C)=d\le22$:
- for $g\le335$ this contradicts the base input;
- for $g\ge336$ it contradicts §1 together with $B^*(d)\le331$.

**Part 2.** Same argument for $d=23,24$, using Claim(12) in Lemma 1.

**Part 3.** The faithful $A_7$-curves with $336\le g\le397$ are listed by `gonality_large_genus.py`:
- the signatures pass the seven-sheeted test;
- for triangle signatures, $A_7$-generating triples exist;
- the one-parameter families start at $g=421$.

## 4. The four window signatures [C]

`certify_signatures.py 24 ...`, output `certify_window_output.txt`. This is the same certified chain as `VERIFICATION.md` §6:
CR lower bounds on $Q_0,Q_1,Q_2$ and a verified Cholesky factorisation. Every curve of the four window signatures passes
Li–Yau with a wide margin:

| signature | $g$ | threshold $48/(g-1)$ | curves (classes) | certified $\lambda_1\ge$ (worst) | Li–Yau $\operatorname{gon}\ge$ (worst) |
|---|---|---|---|---|---|
| $(3,5,5)$ | 337 | 0.14286 | 2 (6) | 0.24031 | 40.37 |
| $(3,4,7)$ | 346 | 0.13913 | 4 (16) | 0.23617 | 40.74 |
| $(3,5,6)$ | 379 | 0.12698 | 3 (6) | 0.14941 | 28.24 |
| $(4,4,5)$ | 379 | 0.12698 | 7 (22) | 0.26140 | 49.40 |

**The same run beyond the window.** It also certified $(3,5,7)$, $(3,6,6)$, $(3,6,7)$, $(3,7,7)$ and $(4,4,6)$ (24 more curves). These are not needed: they lie above 397.
At $(4,4,7)$ ($g=451$, also not needed) the run stopped, because the dense floating Cholesky of the $Q_0$ deflation step was not
positive definite at $n_0=16$. This is a numerical failure of that step, not a certified negative. It does not affect the theorem.

## 5. Remarks

- **Relation to GPT's table.** The dependent-third cost $c(d)$ agrees with GPT's Round-3 numbers (110/121/132). The tangential case, where
  both base points coincide, needs the infinitely-near-point argument in Lemma 4. An earlier version of my own analysis missed it,
  and it costs exactly as much as the transversal case.
- **What is cited.**
  - Castelnuovo–Severi (Stichtenoth Thm 3.11.3);
  - Castelnuovo's genus bound;
  - $\delta=\sum_p\binom{m_p}2$ over infinitely near points, and the proximity equality for branches (Enriques; Casas-Alvero, *Singularities of Plane Curves*, §3.5);
  - the minimum genus 136 and the signature enumeration (exact, `verify_accessory60.py`).
