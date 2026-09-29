# A7 gonality ladder — working notes

Status legend: **[P]** proved (paper), **[C]** verified by computer (scripts in this folder),
**[V]** computer-assisted proof (interval arithmetic plus rigorous floating-point error bounds),
**[L]** checked in Lean, **[?]** conjecture / heuristic, **[X]** failed approach.

**Summary (2026-09-29, second session).** Both summits are proved, via a certified spectral gap
(section 7): for every $(2,4,7)$ $A_7$-curve $C$ and every involution $\tau$,
$$\lambda_1(C)\ge0.33335,\qquad \operatorname{gon}(C)\ge23,\qquad \operatorname{gon}(C/\langle\tau\rangle)\ge12 .$$
The algebraic rungs 1–4 are not needed for this; section 5 records why they could not have finished
the job alone (audit theorem 5.1).

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

## 4. The exact threshold $d=9$ — structure theorem **[P]** (modulo Lemma 4.0, standard)

Throughout, $f:D\to\mathbb P^1$ has degree 9 (so $\operatorname{gon}(D)=9$ by rung 3), and for $h\in H$ put
$f_h=f\circ h$. Write $[f]$ for the pencil, i.e. the subfield $k(f)\subset k(D)$.

**Lemma 4.0 (pencils on a $(9,9)$ curve; standard).** On a smooth curve of bidegree $(9,9)$ in
$\mathbb P^1\times\mathbb P^1$ the only $g^1_9$'s are the two rulings. *(Sketch: $K_D=\mathcal O(7,7)|_D$ and
$H^0(\mathcal O(7,7))\cong H^0(K_D)$; a $g^1_9$ is a length-9 scheme failing by one to impose independent
conditions on $|\mathcal O(7,7)|$; residuation along ruling lines puts it on a ruling line. Cf. G. Martens,
*The gonality of curves on a Hirzebruch surface*, Arch. Math. 67 (1996).)*

**4.1 CS equality with every good involution [P].** For a good involution $v$,
$\gcd(\deg\pi_v,\deg f)=\gcd(2,9)=1$, so $(\pi_v,f):D\to E_v\times\mathbb P^1$ is birational onto its image and
CS reads $64\le56+8$: **equality**. Equality in the Hodge-index proof of CS forces the image to be
smooth and $\equiv 9F_1+2F_2$; since $\mathrm{Pic}(E\times\mathbb P^1)=\mathrm{pr}_1^*\mathrm{Pic}(E)\oplus\mathbb Z$,
$$D=\{s_0x^2+s_1xy+s_2y^2=0\}\subset E_v\times\mathbb P^1,\qquad s_i\in H^0(E_v,L_v),\ \deg L_v=9,$$
with $\pi_{v*}\mathcal O_D=\mathcal O\oplus L_v^{-1}$, $L_v^2=\mathcal O(B_v)$ ($B_v$ = the 18 branch points), and
$s_1^2-4s_0s_2$ cutting out $B_v$. So $h^0(E_v,L_v)\ge2$. With $V=\langle s_0,s_1,s_2\rangle$:
- $\dim V=2$ ⟺ $[f\circ v]=[f]$ ⟺ $D=E_v(\sqrt{\ell_1/\ell_2})$, $B_v=\mathrm{div}\,\ell_1+\mathrm{div}\,\ell_2$, $\ell_i\in H^0(L_v)$;
- $\dim V=3$ ⟺ $[f\circ v]\ne[f]$; then $|V|:E_v\to\mathbb P^2$ pulls a conic back to $B_v$ and
  $D=E_v\times_{\mathbb P^2}(\mathbb P^1\times\mathbb P^1)$.

**4.2 Pencil trichotomy [P].** Let $F=k(f_h:h\in H)$, an $H$-stable subfield. For $[f_{h_1}]\ne[f_{h_2}]$ the
field $k(f_{h_1},f_{h_2})$ has index 1 or 3 in $k(D)$ (the index divides 9 and is $<9$). Index 1: CS gives
$64\le(9-1)^2$, equality, so $D$ is a smooth $(9,9)$ curve (the equality case of Farb–Wolfson's Lemma 2.2).
Index 3: $F'=k(f_{h_1},f_{h_2})$ has genus $\le(3-1)^2=4$, and for any $h$ either $f_h\in F'$ or
$k(F',f_h)=k(D)$; the latter is excluded by CS ($64\le3\cdot4+2\cdot8=28$). Hence exactly one of:
- **(a)** $F=k(f)$: the pencil is $H$-invariant, $f\circ h=\rho(h)\circ f$ with $\rho:H\to\mathrm{PGL}_2$;
- **(b)** $F=k(Y)$ has index 3, $g(Y)\le4$, and $H$ acts on $Y$;
- **(c)** $D$ is a smooth $(9,9)$ curve whose rulings are the $[f_h]$ (Lemma 4.0), and $H$ acts on
  $\mathbb P^1\times\mathbb P^1$.

**4.3 Each case exhibits $D\to T$ as a pullback [P].** Let $T=D/H=C/C(\tau)$, of genus 3.
- (a): $|\ker\rho|$ divides 9, so $\ker\rho\in\{1,C_3\}$. $\ker\rho=C_3$ is impossible: then $f=f_3\circ q$ with
  $q:D\to D_3=D/C_3$ ($g=22$) and $\deg f_3=3$; a non-central good involution $w$ acts on $D_3$ with 18
  fixed points and $g(D_3/w)=7$, and CS gives $22\le14+2$ unless $f_3$ factors through $D_3/w$, which is
  impossible since 3 is odd. So $\rho(H)=D_6\subset\mathrm{PGL}_2$ faithfully. As $k(f)\cap k(T)=k(f)^H$,
  Galois theory gives $k(D)=k(f)\,k(T)$: **$D$ is the normalization of the pullback of the $D_6$-cover
  $\mathbb P^1\to\mathbb P^1/D_6$ along some $\varphi:T\to\mathbb P^1$ of degree 9.**
- (b): $H$ acts faithfully on $Y$ (kernel $C_3$ would give $Y=D/C_3$ of genus $22>4$); likewise **$D$ is the
  normalization of $Y\times_{Y/H}T$** with $\varphi:T\to Y/H$ of degree 3.
- (c): all four good involutions swap the rulings (one preserving both is $(\alpha,\beta)$ with $\le4$ fixed
  points, not 18). In suitable coordinates $c̄$ is the swap and $H_0=\langle$bad involutions$\rangle\cong S_3$
  acts diagonally; the non-central good involutions are $(x,y)\mapsto(\alpha y,\alpha x)$ with fixed curve the
  graph $\Gamma_\alpha$, and $c̄$ fixes the diagonal $\Delta$.

**4.4 Branch data of $D\to T$ [P][C].** $D\to T$ has 13 branch points, all with monodromy of order 2:
$t_1,t_2,t_3$ (monodromy $c̄$; images of $\mathrm{Fix}(\mu)\cup\mathrm{Fix}(\tau\mu)$), $g_1,\dots,g_9$ (non-central good
involutions) and $t_b$ (bad involutions). Check: $126=12\cdot4+6\cdot13$.

**4.5 Theorem (necessary condition) [P].** *If $\operatorname{gon}(D)\le9$ then, on the genus-3 curve $T$,*
$$(\star)\qquad g_1+\dots+g_9\ \sim\ 3(t_1+t_2+t_3).$$
*Proof.* Monodromy of a pullback at $t$ is (monodromy at $\varphi(t)$)$^{e_\varphi(t)}$, and $D\to T$ has no
branch points of order 3 or 6.
- (a): $\mathbb P^1\to\mathbb P^1/D_6$ is branched at $y_6$ (rotation $s$, $s^3=c̄$), $y_1$ (the good reflections)
  and $y_2$ (the bad ones). Over $y_6$ every $e$ must be $\equiv0\pmod3$, and the $c̄$-points are those with
  $e\equiv3\pmod6$; the only partition of 9 giving three of them is $3+3+3$, so $\varphi^*(y_6)=3(t_1+t_2+t_3)$.
  All nine good points lie over $y_1$ with odd $e$, so $\varphi^*(y_1)=\sum g_j$.
- (b): the good points force at least 3 good branch points of $Y\to Y/H$; with $g(Y)\le4$, Riemann–Hurwitz
  leaves only $g(Y/H)=0$, $g(Y)=4$, signature $(2,2,2,2,2)$ = one $c̄$-, three good and one bad branch point.
  Then $\varphi^*(y_c)=t_1+t_2+t_3$ and each good fibre is three of the $g_j$: $\sum g_j\sim3\sum t_i$.
- (c): $\Psi=G_1G_2G_3/\Delta^3$ (with $G_i$ the $(1,1)$-forms of the three graphs) has
  $\mathrm{div}\Psi=\sum\Gamma_i-3\Delta$ and $\Psi\circ h=\pm\Psi$, so $\Psi^2$ descends to $\psi$ on $T$. Since $D$ meets
  $\Delta,\Gamma_i$ transversally in exactly the fixed points, $\pi_D^*\mathrm{div}\,\psi=\mathrm{div}(\Psi^2|_D)=\pi_D^*(\sum g_j-3\sum t_i)$. $\blacksquare$

## 5. First summit by algebra alone: reduced to arithmetic **[P][C][?]** — *superseded by section 7*

By 4.5, $\operatorname{gon}(D)\ge10$ follows from $[\sum g_j-3\sum t_i]\ne0$ in $\mathrm{Pic}^0(T)$.

**Jacobians [C]** (`jacobian.py`, Chevalley–Weil; the multiplicity of $\rho$ in $H^1(C)$ is
$\dim\rho-\dim\rho^{a}-\dim\rho^{b}-\dim\rho^{c}$):
$$\mathrm{Jac}(C)\sim A^{10}\times E_1^{15}\times E_2^{21}\times S^{35},\quad
\mathrm{Jac}(D)\sim A^{4}E_1^{7}E_2^{11}S^{17},\quad \mathrm{Jac}(T)\sim E_2\times S,$$
with $E_1,E_2$ elliptic curves, $S$ an abelian surface and $A$ a threefold with $\mathbb Q(\sqrt{-7})$-action
(from $10\oplus\overline{10}$). The irreducibles $6,14_a,14_b$ do **not** occur in $H^1(C)$.

**Consequence [P].** Degree-0 divisors supported on the branch fibres whose isotypic components lie in
$1\oplus6\oplus14_a\oplus14_b$ are automatically torsion in $\mathrm{Jac}(C)$ (Abel–Jacobi is $A_7$-equivariant).
The divisor $\Phi_g-3\Phi_c$ on $C$ lifting $(\star)$ has **nonzero** $21$- and $35$-components
(`star.py`), so $(\star)$ is a genuine relation among branch-point classes in $E_2$ and $S$. Those classes
lie in the Mordell–Weil groups over the field of definition, so $(\star)$ turns on Mordell–Weil ranks
and torsion orders. **No argument using only the $H$-action on $D$ can decide it** (see 5.1).

**5.1 Audit theorem [P][C].** *There is a smooth curve $D'$ of genus 64 with a faithful action of
$H\cong C_2\times S_3$ having exactly the fixed-point data of $D$ (central good involution: 18; three
good involutions: 18 each; three bad involutions: 2 each; elements of order 3 and 6: free), and
$\operatorname{gon}(D')=9$.* Explicitly (`model99.py`): $D'=\{F=0\}\subset\mathbb P^1\times\mathbb P^1$ of bidegree $(9,9)$ with
integer coefficients supported on $i+j\equiv0\pmod3$, symmetric under $x\leftrightarrow y$ and
$(i,j)\mapsto(9-i,9-j)$; $S_3=\langle z\mapsto\omega z,\ z\mapsto1/z\rangle$ acts diagonally and the swap is the central
involution. Verified mod $p=1000003$ (hence over $\mathbb Q$, as these are open conditions): no singular points;
$F(x,x)$ and $F(x,1/x)$ squarefree of degree 18 (the swap and the good involutions have exactly 18
fixed points); $F(1,1),F(-1,-1)\ne0=F(1,-1)=F(-1,1)$ (each bad involution fixes exactly 2 points of $D'$);
the corner coefficients are nonzero (order 3 and 6 elements act freely). A smooth $(9,9)$ curve is
irreducible of genus $64$, with gonality $\le9$ from a ruling and $\ge9$ by the rung-3 argument, which
only uses this fixed-point data.

**Consequence.** Every invariant used in rungs 1–4 (fixed points, genera of all quotients, the
$H$-character of $H^0(K_D)$, which holomorphic Lefschetz determines from the fixed points) is shared
by $D'$. **So any proof of $\operatorname{gon}(C/\langle\tau\rangle)\ge10$ must use information about $C$ beyond
the $H$-equivariant topology of $D$** — for instance the double cover $C\to D$ with its $C(\tau)$-structure,
the $A_7$-structure, or the arithmetic of branch-point classes in 5. Any argument that only combines
Castelnuovo–Severi, fixed-point parity and the $H$-orbit of the gonal pencil proves too much, since it
would also apply to $D'$. (For $D'$ one lands in case (c); analogous $H$-curves realize (a) and (b), since
4.3 shows those configurations are exactly pullbacks of $H$-covers along maps from a genus-3 curve.)

## 6. A spectral route to both summits — **[P] reduction** (input certified in section 7)

**Li–Yau (Invent. Math. 69, 1982).** For a conformal map $X\to S^2$ of degree $d$,
$\lambda_1(X)\,\mathrm{Area}(X)\le 8\pi d$. With the hyperbolic metric, $\mathrm{Area}(C)=2\pi(2g-2)=540\pi$, so
$$\operatorname{gon}(C)\ge 67.5\,\lambda_1(C),\qquad \operatorname{gon}(C/\langle\tau\rangle)\ge 33.75\,\lambda_1(C)$$
(the second by composing a degree-$d$ map on $D$ with $C\to D$). Hence
**$\lambda_1(C)>0.2667\Rightarrow\operatorname{gon}(D)\ge10$** and **$\lambda_1(C)>0.2370\Rightarrow\operatorname{gon}(C)\ge17$**.
Farb–Wolfson's Lemma 2.2 alone gives $\operatorname{gon}(C)\ge13$ ($136>11^2$), rigorously.

**First-session numerics [C, not rigorous, superseded]** (`spectrum.py`: P1 finite elements on the tiling of
$C$ by 5040 hyperbolic $(\pi/2,\pi/4,\pi/7)$ triangles, each small triangle replaced by a flat one with the
same side lengths — this approximation biases $\lambda_1$ low by about 0.3%; see 7.8 for exact-geometry values):

| triple class | $n$ | $\lambda_1$ (mult.) | next |
|---|---|---|---|
| 0 | 3 | 0.34560 (14) | 0.57248 (21) |
| 0 | 5 | 0.34604 (14) | 0.57161 (21) |
| 1 | 3 | 0.34560 (14) | same spectrum as class 0 |
| 12 | 3 | 0.35847 (14) | 0.38548 (21) |
| 14 | 3 | 0.35847 (14) | same spectrum as class 12 |

So there are two isometry classes, with $\lambda_1\approx0.346$ and $\approx0.358$. If these values are
accurate, Li–Yau gives $\operatorname{gon}(C)\ge24$ and $\operatorname{gon}(C/\langle\tau\rangle)\ge12$ for **every** such curve and
every $\tau$, well past both targets. To make this a proof, a *certified lower bound*
$\lambda_1(C)>0.2667$ is needed (e.g. a guaranteed-lower-bound FEM à la Liu–Oishi, adapted to the
hyperbolic metric, on the $A_7$-isotypic twisted problems). There is a 30% margin.

## 7. Both summits: a certified spectral gap — **[P] + [V]**

**Theorem 7.1.** *Let $C$ be any smooth connected complex curve with a faithful $A_7$-action such that
$C/A_7\cong\mathbb P^1$ with three branch points of inertia orders $2,4,7$, carrying its hyperbolic metric
(curvature $-1$). Then*
$$\lambda_1(C)\ \ge\ 0.33335\qquad(\ge0.34625\ \text{for the triple classes 12, 14}),$$
*and consequently, for every involution $\tau\in A_7$,*
$$\operatorname{gon}(C)\ \ge\ 23,\qquad \operatorname{gon}\bigl(C/\langle\tau\rangle\bigr)\ \ge\ 12 .$$
*In particular the first summit ($\ge10$) and the further summit ($\ge17$) both hold.*

For comparison: $\operatorname{gon}(C)\le56$ and $\operatorname{gon}(C/\langle\tau\rangle)\le28$, via $C\to C/P\to\mathbb P^1$ with $P$
of order 8 containing $\tau$ ($g(C/P)=12$, Brill–Noether gonality $\le7$). Numerically (7.8):
$\lambda_1(C)\approx0.34627$ (classes 0, 1) and $\approx0.3597$ (classes 12, 14), in the isotype $14_a$.

The proof is 7.2–7.7. Sections 7.2, 7.4 and 7.6 are paper proofs. Sections 7.5 and 7.7 are finite
computations done rigorously: the metric coefficients are bounded in ball arithmetic (`python-flint`/Arb),
and the floating-point linear algebra carries a-priori error bounds. Code: `orbifold.py`, `certify.py`,
`cover.py`, `run_certificate.py`; output: `certificate.txt`.

### 7.2 Gonality from $\lambda_1$ (Hersch; Yang–Yau 1980; Li–Yau 1982) [P]

*Lemma.* Let $X$ be a compact Riemann surface with a conformal Riemannian metric, and let
$\varphi:X\to\mathbb P^1$ be holomorphic of degree $d\ge1$. Then $\lambda_1(X)\,\mathrm{Area}(X)\le8\pi d$.

*Proof.* Identify $\mathbb P^1$ with the unit sphere $S^2\subset\mathbb R^3$. The push-forward of the area measure of $X$
by $\varphi$ has no atoms. By Hersch's balancing lemma there is therefore a Möbius transformation $\gamma$ with
$\int_X\gamma\circ\varphi\,dA=0\in\mathbb R^3$. Write $\psi=\gamma\circ\varphi=(x_1,x_2,x_3)$. This is again holomorphic of
degree $d$, and each $x_i$ is Lipschitz with mean $0$, so $\lambda_1\int x_i^2\le\int|\nabla x_i|^2$. Summing over $i$
and using $\sum x_i^2=1$ gives $\lambda_1\,\mathrm{Area}(X)\le\sum_i\int|\nabla x_i|^2$. For a conformal map,
$\sum_i|\nabla x_i|^2=2\,\mathrm{Jac}(\psi)$, so the right side equals $2\cdot4\pi d$. $\blacksquare$

The hyperbolic area is $\mathrm{Area}(C)=2\pi(2g-2)=540\pi$, so $\operatorname{gon}(C)\ge67.5\,\lambda_1(C)$. If $f:D\to\mathbb P^1$ has
degree $d$, then $f\circ\pi:C\to\mathbb P^1$ has degree $2d$, so $\operatorname{gon}(C/\langle\tau\rangle)\ge33.75\,\lambda_1(C)$.
With $\lambda_1\ge0.33335$ these give $22.50$ and $11.25$, hence 23 and 12.

*Audit remark.* The audit theorem 5.1 says a proof must use more of $C$ than the $H$-equivariant
topology of $D$. Here the extra input is the conformal structure of $C$, through its hyperbolic metric.
Consistently, the audit curve $D'$, with $\operatorname{gon}(D')=9$, is not $C/\langle\tau\rangle$ for any
$(2,4,7)$-curve $C$.

### 7.3 An exact model of $C$ [P]

**Tiling.** $C$ is tiled by $2\cdot2520$ copies of the hyperbolic triangle $T$ with angles $\pi/2$ at $A$ (over the
order-2 point), $\pi/4$ at $B$ and $\pi/7$ at $C$. Its sides are $AB=\operatorname{arccosh}\frac{\cos(\pi/7)}{\sin(\pi/4)}\approx0.7218$,
$AC=\operatorname{arccosh}\frac{\cos(\pi/4)}{\sin(\pi/7)}\approx1.0704$ and $BC=\operatorname{arccosh}(\cot\frac\pi4\cot\frac\pi7)\approx1.3697$. The tiles
are $U_g,L_g$ for $g\in A_7$, glued by
$U_g\sim L_g$ along $AB$, $U_g\sim L_{gb}$ along $BC$, and $U_g\sim L_{ga}$ along $CA$.
These gluings are right multiplications, and $A_7$ acts on the left.

Around a vertex over $a$, $b$, $c$ there are $4$, $8$, $14$ triangles respectively (for $c$ because $ba$ has
order 7), so the angles sum to $2\pi$. The monodromy triple is $(a,b^{-1},ba)$, whose product is
$a\,b^{-1}\,b\,a=1$. The map $(a,b)\mapsto(a,b^{-1})$ permutes the 96 generating triples with fixed $a$, compatibly
with conjugation by $C(a)$. So the class representatives 0, 1, 12, 14 produce all four curves.

**Chart.** Every tile, $U$ or $L$, is identified isometrically with one placement of $T$ in the Klein
model. The midpoint of $BC$ sits at the origin; this placement is computed in ball arithmetic. The
triangle is parametrised affinely by $R=\{(s,t):s,t\ge0,\ s+t\le1\}$ via
$x=X_A+s(X_B-X_A)+t(X_C-X_A)$. Tiles sharing a side are mirror images across it, so the side gets the same
$(s,t)$ from both tiles, and **every gluing is the identity on the sides of $R$**. In $R$-coordinates the
hyperbolic Dirichlet form and area are exactly
$$\int\nabla u\cdot A_R\nabla u,\quad\int u^2w_R,\qquad A_R=|\det J|\,J^{-1}\frac{I-xx^T}{\sqrt{1-|x|^2}}J^{-T},\quad
w_R=\frac{|\det J|}{(1-|x|^2)^{3/2}},$$
where $J=[X_B-X_A,\ X_C-X_A]$ and the middle factor is the Klein metric. The uniform subdivision of $R$ into
$n^2$ right isosceles triangles is therefore a conforming mesh of $C$ and of every quotient below.

### 7.4 Reduction to three small twisted problems [P][C]

For $K\le A_7$ and a character $\varepsilon:K\to\{\pm1\}$, the $\varepsilon$-twisted problem is the Laplacian on
$\{f:f(kx)=\varepsilon(k)f(x)\}$. It lives on the $2[A_7:K]$ tiles $U_{Kg},L_{Kg}$, with gluings twisted by signs
(`orbifold.quotient_tiles`).

*Lemma.* Let $E$ be an eigenspace of $\Delta_C$ for $\lambda$, and suppose $E$ contains an irreducible $\rho$ with
$\langle\rho|_K,\varepsilon\rangle>0$. Then $E$ contains a nonzero $\varepsilon$-twisted function, so $\lambda$ is an eigenvalue of
the twisted problem. (Restricted to $K$, $E$ contains $\varepsilon$.)

Take the three problems

| | $K$ | $[A_7:K]$ | $\varepsilon$ | irreducibles seen (multiplicity) |
|---|---|---|---|---|
| $Q_0$ | $A_7$ | 1 | trivial | $1$ (use the first **nonzero** eigenvalue) |
| $Q_1$ | $3^2{:}4<A_6$, generated by $(012),(345)$ [$\varepsilon=+1$], $(0314)(25)$ [$\varepsilon=-1$] | 70 | sign | $6,14_a,14_b,15,21$ |
| $Q_2$ | $S_4$, generated by $(04)(1362)$ [$\varepsilon=+1$] and $(16)(45)$ [$\varepsilon=-1$] | 105 | sign | $10,\overline{10},15,35\,(2)$ |

Together they see every irreducible of $A_7$. `cover.py` checks this exactly: $K_1,K_2$ contain no elements of
order 7, so only integer character values enter; it also verifies exact orthonormality of the table and the
Frobenius counts $\sum_\rho m(\rho)\dim\rho=[A_7:K]$. Hence
$$\lambda_1(C)\ \ge\ \min\bigl(\mu_2(Q_0),\ \mu_1(Q_1),\ \mu_1(Q_2)\bigr).$$
Numerical cross-check (P1, $n=2$, all eigenvalues of $C$ below 1.6): every nonzero eigenvalue of $C$ appears
in $Q_1\cup Q_2$ with the predicted multiplicity, and nothing else appears. The reduction replaces the
5040-triangle surface by surfaces of 2, 140 and 210 triangles. $Q_1,Q_2$ have no constant functions, so
their certificates are plain positive-definiteness checks.

### 7.5 Comparison with piecewise-constant coefficients [P][V]

On each mesh triangle $e$, take $P_e=A_R(\text{centroid})$, rounded to doubles. Let $c_e$ be a rigorous lower bound
for $\min_{z\in e}\lambda_{\min}(P_e^{-1}A_R(z))$, and $w_e$ a rigorous upper bound for $\max_e w_R$. Both are computed by ball
arithmetic over $3\times3$ sub-boxes of the cell. For $\lambda_{\min}$ we use the congruence $W=L^{-1}QL^{-T}$ with
$P_e=LL^T$, plus Weyl's inequality for the radius; this loses $O(\text{width})$, not $O(\sqrt{\text{width}})$.
Then $A_e:=c_eP_e\preceq A_R$ and $w_e\ge w_R$ on $e$. By min–max, the eigenvalues of the true problem are $\ge$
those of the comparison problem $(a_c,b_c)$. The loss is first order in $h$: $\min_e c_e=0.95$ at $n=32$.

### 7.6 Guaranteed lower bounds from Crouzeix–Raviart elements [P]

Let $V_h$ be the (twisted) Crouzeix–Raviart space on the mesh, and $\Pi$ the CR interpolation, i.e. edge means
in $R$-coordinates. It is consistent across tile sides because the gluings are the identity there.
1. $a_{c,h}(u-\Pi u,v_h)=0$ for all $v_h\in V_h$. On each $e$, $A_e\nabla v_h$ is constant, so
   $\int_e\nabla(u-\Pi u)\cdot A_e\nabla v_h=\sum_E(A_e\nabla v_h\cdot n_E)\int_E(u-\Pi u)=0$.
2. $\|u-\Pi u\|_{b_c}\le C_h\,|u-\Pi u|_{a_{c,h}}$ with
   $C_h^2=\max_e\kappa_e^2\,w_e/\lambda_{\min}(A_e)$ and $\kappa_e^2=\frac{(1/n)^2}8+\frac{(\sqrt2/n)^2}{j_{1,1}^2}$.
   This uses, for $f$ with zero mean on the edge $E$ opposite $P$,
   $\|f\|_{L^2(e)}\le(\max_{x\in E}|P-x|^2/8+h_e^2/j_{1,1}^2)^{1/2}|f|_{H^1(e)}$ (Carstensen–Gedicke–Rim 2012, Lemma 2.2, as
   used in Carstensen–Gedicke, Math. Comp. 2014, Thm 2.1). Take $E$ to be the hypotenuse.

*Theorem (Liu 2015; Carstensen–Gedicke 2014).* $\lambda_k\ge\lambda_{k,h}/(1+C_h^2\lambda_{k,h})$.

*Proof.* Let $E_k$ be spanned by the first $k$ exact eigenfunctions, and take $u\in E_k$ with $\|u\|_b=1$ and
$t=|u-\Pi u|_{a_h}$. By (1), $|\Pi u|_{a_h}^2=|u|_a^2-t^2\le\lambda_k-t^2$. By (2), $\|\Pi u\|_b\ge1-C_ht$. So
$R(\Pi u)\le(\lambda_k-t^2)/(1-C_ht)^2\le\lambda_k/(1-C_h^2\lambda_k)$; the maximum over $t$ is at $t=C_h\lambda_k$.
$\Pi$ is injective on $E_k$: $\Pi u=0$ would give $1\le C_h\sqrt{\lambda_k}<1$. Min–max then gives
$\lambda_{k,h}\le\lambda_k/(1-C_h^2\lambda_k)$. If instead $C_h^2\lambda_k\ge1$, the claim is trivial. $\blacksquare$

### 7.7 Verifying $\lambda_{1,h}>\sigma$ in floating point [V]

The CR mass matrix is diagonal, since CR basis functions are $L^2$-orthogonal on each triangle. So
$\lambda_{1,h}>\sigma$ iff $B=K_h-\sigma M_h\succ0$.

The exact comparison matrices are built from the double-precision data $(c_e,P_e,w_e)$. Element matrices are
computed in ball arithmetic, and summation errors are bounded by $\gamma_m\sum|\cdot|$. Together these give
$\|B_{\text{exact}}-\hat B\|_2\le\eta$ (a row-sum bound).

Sparse Cholesky (CHOLMOD) of $\hat B-cI$ runs to completion. By the standard backward error bound (Higham,
*Accuracy and Stability*, Thm 10.3), the computed factor satisfies $\tilde L\tilde L^T=Q(\hat B-cI)Q^T+\Delta$ for a
permutation $Q$, with $|\Delta|\le\gamma_{k+1}|\tilde L||\tilde L^T|$, where $k$ is the largest row count of $\tilde L$. Using
$\||\tilde L||\tilde L^T|\|_2\le\|\tilde L\|_1\|\tilde L\|_\infty$,
$$B_{\text{exact}}\succeq\bigl(c-\gamma_{k+1}\|\tilde L\|_1\|\tilde L\|_\infty-u(\max|b_{ii}|+c)-\eta\bigr)I\succ0 .$$
The shift $c$ is chosen twice the needed amount. For $Q_0$ the constants are deflated: $B+\alpha zz^T\succ0$ with
$z=M_h\mathbf 1$ and $\alpha>\sigma/\mathbf 1^TM_h\mathbf 1$ (dense Cholesky, $k=N$). This implies $B\succ0$ on $\mathbf 1^{\perp_M}$, so
$\lambda_{2,h}>\sigma$. Finally, $\lambda\ge\sigma/(1+C_h^2\sigma)$ is evaluated in ball arithmetic and rounded down.

### 7.8 Results

`certificate.txt` (from `python3 run_certificate.py`, about 3 minutes; $n=32$ for $Q_1,Q_2$, $n=16$ for $Q_0$):

| triple class | $Q_0$ | $Q_1$ | $Q_2$ | certified $\lambda_1(C)\ge$ |
|---|---|---|---|---|
| 0, 1 | 0.99783 | **0.33335** | 0.68902 | 0.33335 |
| 12, 14 | 0.99783 | **0.34625** | 0.76298 | 0.34625 |

At $n=32$ the Cholesky shifts are $\sim10^{-10}$, about 2 times the rigorous need, while the true margin
$(\lambda_{1,h}-\sigma)\min M_{ii}$ is $\sim10^{-7}$.

**Non-rigorous numerics with exact geometry** (`orbifold.py`, $Q_1$, class 0):
- conforming P1, which gives upper bounds up to quadrature: $0.346707$ ($n=12$), $0.346378$ ($n=24$);
- plain CR, which approaches from below: $0.345326$ ($n=8$), $0.346029$ ($n=16$);
- both extrapolate to **$\lambda_1(C)\approx0.34627$**.

For class 12, the index-21 quotient $C/S_5$ gives $0.36008$ at $n=12$ (P1). The first eigenfunctions lie in
$14_a$ (the $(5,2)$ representation): they already appear on $C/S_5$, with $S_5$ the stabiliser of a 2-set, and
not on $C/\mathrm{PSL}_2(7)$, which sees only $1+14_b$.

### 7.9 What the spectral theorem does and does not give

- **Rung 5 (first summit)** is settled: $\operatorname{gon}(C/\langle\tau\rangle)\ge12$ for every $C$ and every $\tau$. Theorem 4.5 and $(\star)$
  are no longer needed. Whether $(\star)$ holds on $T$ remains an open arithmetic question of independent
  interest.
- **Rung 6 (further summit)** is settled: $\operatorname{gon}(C)\ge23$. This part is independent of the quotient result. The
  quotient only gives $\operatorname{gon}(C)\ge\operatorname{gon}(D)\ge12$ (norm map, Farb–Wolfson Lemma 2.4), and Farb–Wolfson Lemma 2.2
  gives 13.
- **Limits of the method.** Li–Yau with $\lambda_1\approx0.3463$ cannot give $\operatorname{gon}(D)\ge13$, which would need
  $\lambda_1>0.3556$. $\operatorname{gon}(C)\ge24$ needs a certified $\lambda_1>0.34074$ for classes 0, 1 (true value 0.34627), which
  requires a finer mesh or a second-order comparison. It is already certified for classes 12, 14.
- **Other quotients.** For $K\le A_7$ the same argument gives
  $\operatorname{gon}(C/K)\ge67.5\,\lambda_1^K(C)/|K|\ge22.5/|K|$, where $\lambda_1^K\ge\lambda_1$ is the first eigenvalue of
  $K$-invariant functions.
