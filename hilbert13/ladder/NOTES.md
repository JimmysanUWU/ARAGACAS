# A7 gonality ladder — working notes

Status legend: **[P]** proved (paper), **[C]** verified by computer (scripts in this folder),
**[L]** checked in Lean, **[?]** conjecture / heuristic, **[X]** failed approach.

Setup. $C$ smooth connected, faithful $A_7$-action, $C/A_7\cong\mathbb P^1$ branched over three
points with inertia orders $2,4,7$. Equivalently a generating triple $(a,b,c)$ of $A_7$ with
$|a|=2,|b|=4,|c|=7$, $abc=1$. $\tau\in A_7$ an involution, $D=C/\langle\tau\rangle$.

## 0. Group-theoretic facts

- **[C]** $(2,4,7)$ generating triples exist: with $a=(12)(34)$ fixed there are 96 choices of $b$;
  10080 triples in total, forming **4 orbits under $A_7$-conjugation and 2 under $S_7=\mathrm{Aut}(A_7)$**.
  So there are (at least) two non-isomorphic such curves; every argument below uses only
  class-level data, hence applies to all of them (`triples.py`).
- **[P]** Involutions of $A_7$: a single class, cycle type $2^2 1^3$, 105 elements,
  $|C_{A_7}(\tau)|=24$. Elements of order 4: type $(4,2,1)$, 630 elements, $|C(b)|=4$, and
  $b^2$ is an involution. Elements of order 7: two classes of 360, $|C|=7$.
- **[P]** $(2,4,7)$ is not on Singerman's list of non-maximal signatures, so for these curves
  $\mathrm{Aut}(C)=A_7$ (standard for triangle curves with maximal signature). *Only used as
  context; nothing below depends on it.*

## 1. Basic geometry  — **done [P][C]**

Riemann–Hurwitz: $2g(C)-2=2520\,(-2+\tfrac12+\tfrac34+\tfrac67)=2520\cdot\tfrac3{28}=270$, so
$$g(C)=136.$$

Fixed points. Over a branch point with stabilizer $\langle x\rangle$ (fiber $=A_7/\langle x\rangle$), an
element $h$ fixes $\#\{g: g^{-1}hg\in\langle x\rangle\}/|x| = |C(h)|\cdot|h^{A_7}\cap\langle x\rangle|/|x|$ points.
For an involution $h$ (unique class, $|C(h)|=24$):
- over the order-2 point: $\langle a\rangle$ contains 1 involution: $24\cdot1/2=12$;
- over the order-4 point: $\langle b\rangle$ contains exactly one involution $b^2$ (type $2^2$): $24\cdot1/4=6$;
- over the order-7 point: none.

So **every involution has exactly $12+6=18$ fixed points** (confirmed by direct enumeration of
cosets in `rung1.py` for the triple classes). Similarly: order 4 → 2 fixed points (all over the
order-4 point); order 7 → 3; orders 3, 5, 6 → 0.

Riemann–Hurwitz for $C\to D$: $270=2(2g(D)-2)+18$, so
$$g(D)=64.$$

## 2. Geometry inherited by $D$  — **done [P][C]**

Automorphisms of $D$ coming from $A_7$: $H:=N(\langle\tau\rangle)/\langle\tau\rangle=C(\tau)/\langle\tau\rangle$, order 12.
For $\tau=(12)(34)$: $C(\tau)=(D_4\times S_3)\cap A_7$ and
$H\cong\{(x,y)\in V_4\times S_3:\operatorname{sgn}x=\operatorname{sgn}y\}\cong C_2\times S_3$.

**Fixed-point formula on a quotient [P].** If $V\trianglelefteq N$ and $w\in N$, then
$\#\mathrm{Fix}_{C/V}(\bar w)=\frac1{|V|}\sum_{g\in wV}\#\mathrm{Fix}_C(g)$.
(Proof: for a $V$-orbit $O$ fixed by $\bar w$, $\{g\in wV:gp=p\}$ is a coset of $\mathrm{Stab}_V(p)$, so
each fixed orbit contributes $|O|\cdot|\mathrm{Stab}_V(p)|=|V|$.)

**Involution $\mu\neq\tau$ commuting with $\tau$ [P].** $\langle\tau,\mu\rangle\cong V_4$; $\tau\mu$ is an
involution. A point $[p]\in D$ is fixed by $\bar\mu$ iff $\mu p\in\{p,\tau p\}$, i.e.
$p\in\mathrm{Fix}(\mu)\sqcup\mathrm{Fix}(\tau\mu)$. The union is disjoint because a common fixed point would have
stabilizer $\supseteq V_4$, but point stabilizers on a curve are cyclic. For the same reason no
point of $\mathrm{Fix}(\mu)$ is fixed by $\tau$, and $\tau$ preserves $\mathrm{Fix}(\mu)$ (they commute). So the
$18+18=36$ points form 18 free $\tau$-orbits:
$$\#\mathrm{Fix}_D(\bar\mu)=\tfrac{18+18}2=18,\qquad g(D/\bar\mu)=g(C/V_4)=28 .$$
Call these the **good involutions** of $D$. There are 4 in $H$: images of $(13)(24)$ (central in
$H$) and of $(12)(56),(12)(57),(12)(67)$.

**Other elements of $H$ [P][C]** (`quot.py`):
- 3 **bad involutions**: images of $\mu$ of order 4 with $\mu^2=\tau$ (e.g. $(1324)(56)$). Here
  $\mu p=\tau p\iff\mu p=p$, so $\mathrm{Fix}_D=\mathrm{Fix}_C(\mu)$ = 2 points; $g(D/\bar\mu)=32$.
- order 3 and order 6 elements act **freely** on $D$ (no element of order 3 or 6 in $A_7$ has
  fixed points on $C$). $g(D/C_3)=22$.
- Genera of $D/\bar K$ for all 16 subgroups $\bar K\le H$: 64; 28 (×4), 32 (×3); 22; 12 (×3, the
  $V_4$'s); 10, 11, 7 (the three subgroups of order 6); $g(D/H)=g(C/C(\tau))=3$.

Cross-check [C]: $D\to D/H$ has $\sum(|H_p|-1)=4\cdot18+3\cdot2=78$ and
$126=12(2\cdot3-2)+78$. ✓

## 3. First barrier: $\operatorname{gon}(D)\ge 9$  — **done [P]**

**Castelnuovo–Severi (CS).** If $f_i:X\to Y_i$ have degrees $d_i$ and do not factor through a
common map of degree $>1$, then $g(X)\le d_1g(Y_1)+d_2g(Y_2)+(d_1-1)(d_2-1)$.

Let $f:D\to\mathbb P^1$ have degree $d\le 8$ and let $v$ be a good involution, $\pi_v:D\to E_v=D/v$
($g=28$). Since $\deg\pi_v=2$ is prime, a common factor forces $f$ to factor through $\pi_v$. If
not, CS gives $64\le 2\cdot28+(d-1)=55+d$, i.e. $d\ge9$. Hence **$f$ factors through $\pi_v$,
i.e. $f\circ v=f$ exactly** (not merely $f\circ v=M\circ f$).

So $f$ is invariant under the subgroup generated by the four good involutions. That subgroup is
all of $H$: the three images of $(12)(5x)$ generate an $S_3$ (products of two are the 3-cycles
on $\{5,6,7\}$), and the central good involution adds the $C_2$. Hence $f$ factors through
$D\to D/H$, so $12\mid d$, contradicting $d\le 8$. $\blacksquare$

**The distinction $f\circ v=f$ vs. $f\circ v=M\circ f$.** CS yields the first (factorization). The
second (the *pencil* $|f|$ is $v$-stable, with $v$ acting on the target by a Möbius map $M$) is
weaker and is what one gets from uniqueness of a gonal pencil. At $d=9$, $f\circ v=f$ is
impossible for parity reasons ($9$ is odd), so everything must happen at the level of pencils.
