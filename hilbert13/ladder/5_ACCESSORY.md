# Chapter 5. Accessory irrationalities: $a(A_7)=60$ and $\mu(A_7)=90$

**Idea.**
- **Compressions are families of divisors.** A curve compression after an accessory is a family of effective divisors on the target curve, parametrised by the base.
- **Over a linear base** the family is a single linearised moving line bundle (Theorem 5.2). So accessory degrees are degrees of linearised series, induced up from subgroups (Theorem 5.3).
- **Two facts decide $A_7$:**
  - *Castelnuovo's bound* forbids low-degree models of high-genus curves. The cheapest route is therefore through the index-15 Klein subgroup: $15\cdot4=60$.
  - *The degree lattice* puts linearised degrees on a $(2,4,7)$ curve in $90\mathbb Z$, and holomorphic Lefschetz shows that degree 90 is attained. So the cheapest series with full monodromy has degree 90.
- **Over a general base** the family is still one linear system, provided the base's Albanese cannot map to $\mathrm{Jac}\,C$. But its Schur obstruction can be cancelled by the base, and exactly by the base's Amitsur subgroup (§5.5).

**Definitions.** Let $K=\mathbb C(V)^G$ for a faithful linear model $V$.
- An *accessory* is a finite extension $F/K$. After it, the monodromy of the torsor may drop to a subgroup $H$.
- A *curve compression* is a faithful $H$-stable curve field $E\subset LF$ with $LF=FE$.
- **$a(G)$** is the least accessory degree that admits a curve compression.
- **$\mu(G)$** is the least degree of a $G$-*linearised*, base-point-free line bundle with $h^0\ge2$ on a faithful connected $G$-curve.

**Theorem 5.1.**
1. $a(A_7)=60$; that is, $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1=\mathrm{ed}_{\mathbb C}(A_7;\le60)$.
2. $\mu(A_7)=90$: compression with connected full $A_7$-monodromy needs exactly accessory degree 90.

**Credits.**
- GPT: part 1, Theorems 5.2–5.3 and Lemmas 5.4–5.5 (DAY 2, `gpt/DAY_2_A7_Mathematical_Reference.pdf`); $\mu\ge72$ (Round 5). All of these are re-derived here.
- GPT (3 October): the open-base gap and the Amitsur theory of §5.5 (`gpt/A7_Amitsur_Correspondence.pdf`), re-derived.
- Ours: $\mu\le90$, the exclusion of 72 and 84, and the first (flawed) version of §5.5.

Computations: `verify_accessory60.py`, `verify_mu90_exact.py`, `equivariant_rr.py`.

## 5.1 From compressions to line bundles

**Theorem 5.2 (correspondence) [P].** An equivariant degree-$d$ cover $W\dashrightarrow V$ with a dominant equivariant map $W\dashrightarrow C$ to a faithful $G$-curve gives $\mu(G)\le d$.

*Proof.*
1. The closure $Z\subset V\times C$ of the joint image is an invariant divisor, so $\mathcal O(Z)$ is canonically linearised.
2. Since $\mathrm{Pic}(V\times C)=\mathrm{Pic}(C)$, all fibres are linearly equivalent: a moving divisor of degree $e\mid d$.
3. Restrict to $\{0\}\times C$ and remove the invariant fixed part. $\square$

**Theorem 5.3 (induction) [P].** $a(G)=\min_{H\le G}[G:H]\,\mu(H)$.

*Proof.*
- **Lower bound.** If $H$ is the monodromy after the accessory, then $[F:K]=[G:H]e$; apply Theorem 5.2 to $H$.
- **Upper bound.** Induce from $H$ and use the incidence variety (DAY 2 Prop. 1.4; Round 5 Prop. 1.1 for disconnected base change). $\square$

**Lemma 5.4 (degree lattice) [P].** A linearised bundle on a faithful $G$-curve has degree in $N_G(C)\mathbb Z$, where $N_G(C)=|G|/\mathrm{lcm}(e_i)$ over the inertia orders $e_i$. (By Hilbert 90 there is an invariant rational section, whose divisor is a sum of orbits.)

**Lemma 5.5 (linearised Castelnuovo) [P].** Let $G$ be simple, and let $L$ be minimal of degree $n<|G|$. Then $|L|$ maps $C$ birationally onto a faithful $Y\subset\mathbb P^r$ with $r\ge q(G)-1$ and $g(Y)\le\pi(n,r)$. Here $q(G)$ is the least nontrivial representation degree; $q(A_7)=6$.

*Proof.*
1. **Faithful.** The kernel of $G$ on the image has order dividing the degree of $C\to Y$, which is $<|G|$. Since $G$ is simple, the kernel is trivial.
2. **Birational.** If $C\to Y$ had degree $k\ge2$, the hyperplane class of $Y$ would be a linearised moving bundle of degree $n/k<n$, contradicting minimality.
3. **Dimension.** $H^0(L)$ is a genuine representation. It contains a nontrivial irreducible, since otherwise $G$ would act trivially on $Y$. So $r+1\ge q(G)$.
4. **Genus.** Castelnuovo bounds the genus of a nondegenerate curve of degree $n$ in $\mathbb P^r$. $\square$

## 5.2 $a(A_7)=60$

**Proposition 5.6 [P][X].** $\mu(A_7)\ge60$.

*Proof.* Inertia orders divide 420, so $6\mid n$, and $n<60$ means $n\le54$. Lemma 5.5 then gives $g\le\pi(54,5)=325$. In genus 136–325 the only lattice-compatible cases are:
- $(2,5,7)$ with $n=36$, where $\pi=136<199$;
- $(3,4,5)$ with $n=42$, where $\pi=190<274$. $\square$

**Subgroups of index $<60$** (ATLAS):

| $H$ | index | $\mu(H)\ge$ | reason | product |
|---|---|---|---|---|
| $A_6$ | 7 | 9 | $\pi(8,4)=5<10=g_{\min}(A_6)$ | 63 |
| $L_2(7)$ | 15 | 4 | $\pi(3,2)=1<3$ | **60** |
| $S_5$ | 21 | 3 | a hyperelliptic involution is central; $S_5\not\subset\mathrm{PGL}_2$ | 63 |
| $(A_4\times C_3){:}2$ | 35 | 2 | $C_3^2\not\subset\mathrm{PGL}_2$ | 70 |
| $A_5$ | 42 | 2 | no 2-dimensional representation | 84 |

**Upper bound.** $L_2(7)$ with the canonical bundle of the Klein quartic gives $15\cdot4=60$. This uses the Farb–Wolfson independence of the linear model; the lower bound needs only §5.1.

**Small groups (Round 5) [P].**
- $a(A_5)=2$: Kronecker–Klein's square root.
- $a(L_2(7))=4$.
- $a(A_6)=12$: $6\mid n$ and $\pi(6,4)=2$; there are no subgroups of order 40 or 45; the one of order 36 is not a Möbius group.

## 5.3 $72\le\mu(A_7)\le90$

**Proposition 5.7 (GPT) [P][X].** $\mu(A_7)\ge72$.

*Proof.*
1. By Proposition 5.6 and the lattice, $n\in\{60,66\}$, with signatures $(2,6,7)$ and $(3,4,7)$.
2. $\pi(60,8)=220$ forces $W=H^0(L)$ to be the standard 6.
3. The $A_6$-fixed vector of $W$ has a zero divisor of degree 60, while $A_6$-orbits have size $\ge120$ (resp. 90). $\square$

**Further exclusions.**
- Degrees 78 and 102 have no lattice signature.
- Degree 114 forces $(3,4,5,7)$ ($g=1354$) and $W=6$, since $\pi(114,6)=1221$. Then $\prod x_j$ is an invariant section of degree 798, which is not a sum of orbit sizes $\{360,504,630,840,2520\}$.
- The power-sum curve $\{\sum x_j^k=0,\ k\le5\}\subset\mathbb P^6$ has genus 481, signature $(3,7,7)$ and $\mathcal O(1)$ of degree 120. So $\mu\le120$.

**Proposition 5.8 [P][X].** On every $(2,4,7)$ curve, $B+T=D_2-3D_4+2D_7$ is a linearised moving bundle of degree 90. So $\mu(A_7)\le90$.

*Proof.*
1. The linearised Picard group is $\langle D_2,D_4,D_7\mid2D_2=4D_4=7D_7\rangle\cong\mathbb Z\oplus\mathbb Z/2$. Its degree-90 classes are $B=2D_7-D_4$ and $B+T$.
2. Holomorphic Lefschetz: for $g\ne1$, $\sum_q(-1)^q\mathrm{tr}(g\mid H^q(\mathcal O(D)))=\sum_{p\in\mathrm{Fix}(g)}a_p^{k_p}/(1-a_p^{-1})$, with $k_p=\mathrm{mult}_pD$. It gives $\chi(B+T)=-6+10-14_a-14_b-21$, so $h^0(B+T)\ge10$.
3. A fixed part would be an invariant divisor of degree $\le90<360$. $\square$

*Checks.*
- $\chi(K)$ gives $H^0(K)=10+2\cdot\overline{10}+15+21+2\cdot35$ (Chevalley–Weil); the other conventions give non-integral multiplicities.
- $\chi(B)=-\overline{10}-35$ and $\chi(2B)=6-\overline{10}+14_a+14_b+21$.
- The other orientation exchanges $10$ and $\overline{10}$.

## 5.4 $\mu(A_7)\ne72,84$

**Candidates [X].** Let $L$ realise $\mu=n$. Then $g\le\pi(n,h^0-1)$ and $N(C)\mid n$. With the seven-sheeted test and generation, this leaves:

| $n$ | signature (genus): $h^0\le$ |
|---|---|
| 72 | $(2,5,7)$ (199): 12; $(3,5,7)$, $(4,5,7)$, $(5,5,7)$: 7, 6, 6 |
| 84 | $(3,4,5)$ (274): 12; $(3,5,6)$ (379): 10; $(4,5,6)$, $(5,5,6)$, $(5,6,6)$, $(5,6,7)$, $(2,2,3,5)$, $(2,2,5,6)$, $(2,3,3,5)$: $\le8$ |

Orbits are larger than $n$, so $W=H^0(L)$ has no trivial summand. Hence $W\in\{6,10,\overline{10},6\oplus6\}$.

**Lemma 5.9 (orbit semigroup) [P].** If $f\in\mathrm{Sym}^kW$ is invariant and $kn$ is not a nonnegative combination of orbit sizes, then $f|_C\equiv0$ (its divisor would be a sum of orbits). Lemma 7.3 sharpens this with local characters.

**Case $W\supseteq6$.** Use coordinates $x_1,\dots,x_7$ with $\sum x_j=0$ and power sums $p_k$.
- **$n=72$.**
  - $p_2,p_3,p_4,p_6$ vanish, since 144, 216, 288, 432 are not sums of orbit sizes.
  - So the image lies in $\{\text{roots of }t^7-ut^2-v\}$, which is irreducible (the monodromy of $t^7-t^2$ is $S_7$) and generically reduced of degree $144>72$.
- **$n=84$.** $p_2,p_3,p_4,p_7$ vanish, and $p_7=7e_7$, so $\prod x_j=0$ on $C$. By transitivity the image lies in every coordinate hyperplane.

**Case $W\in\{10,\overline{10}\}$.** $A_7\subset SO(6)$ lifts to $2.A_7\subset SU(4)$. On the half-spin $V_4$, $\wedge^2V_4=6$ and $\mathrm{Sym}^2V_4=\overline{10}$ (irreducible). So $\varphi_W$ is a family of symmetric $4\times4$ matrices $Q_p$. This occurs only for $(2,5,7)$ at 72 and for $(3,4,5)$, $(3,5,6)$ at 84.

**Lemma 5.10 (parity) [P].** Let $c$ fix $p$, with lift $s$. Then $Q_p$ lies in the $\lambda$-eigenspace of $c$. If $\lambda\ne\mu^2$ for every eigenvalue $\mu$ of $s$, every matrix there has even rank. (In an eigenbasis, $Q_{ab}\ne0$ only if $\mu_a\mu_b=\lambda$, and $\mu\mapsto\lambda/\mu$ is a fixed-point-free involution.)

By Lemma 5.9, $\det Q_p\equiv0$, so the generic rank is $r\le3$.
- **$r=2$.** $\wedge^2Q_p$ maps to $\mathbb P(6)$ by a bundle $N$ with $N^2=L^2(-B')$. So $\deg N=n$ and $W\supseteq6$ for $N$, which is the previous case.
- **$r=3$.** The adjugate gives a kernel map with $\kappa^*\mathcal O(2)=L^3(-B'')$. Some branch eigenspace contains no rank-3 matrix with $\det=0$. So $B''$ contains an orbit of size 1260, 630 or 420, all $>3n$.
- **$r=1$.** Then $Q_p=v_p\otimes v_p$.
  - Lemma 5.10 applies at the $e=2$ branch for $(2,5,7)$, and at $e=4$ for $(3,4,5)$.
  - The two classes of $(3,5,6)$ without a parity branch give $\varphi=v_2\circ\psi$, with $\psi:C\to\mathbb P^3$ of degree 42. Here $\psi$ is birational ($\pi(21,3)=90<136$) and lies on no quadric ($\mathrm{Sym}^2V_4$ is irreducible). Halphen then gives $g\le42\cdot39/6+1=274<379$.

So $\mu\notin\{72,78,84\}$, and $\mu(A_7)=90$. $\square$

**Exactness [X]** (`verify_mu90_exact.py` → `mu90_exact_output.txt`).
- Lift spectra are forced: $\{\mu_a\mu_b\}_{a<b}$ is the spectrum on the 6 and $\prod\mu_a=1$. The solution is unique up to $s\mapsto-s$, and up to $V_4\leftrightarrow V_4^\*$ for 7-cycles.
- Branch eigenspaces are 0/1 pattern spaces of symmetric matrices.
- Rank 1 occurs iff the pattern has a loop. Rank 3 is decided exactly, using the determinant, the $3\times3$ minors and radical membership over $\mathbb Q$.
- All 84 (curve, class) cases close, under both orientations.

## 5.5 Beyond fixed points: Amitsur subgroups [P]

Theorem 5.2 restricts $\mathcal O(Z)$ to a fixed point, and fixed points do not survive a quadratic accessory $t^2=q(v)$ (Round 5 Prop. 4.2). What replaces the fixed point is the *Amitsur subgroup* of the base. An earlier version of this section claimed that $\mathrm{Pic}(B)=0$ on an open base forces linearisation. That is false (Example 5.12); GPT found the gap and supplied the corrected theory (`gpt/A7_Amitsur_Correspondence.pdf`), which is re-derived here.

**Notation.** For a smooth projective $G$-variety $X$ and an invariant class $M\in\mathrm{Pic}(X)^G$, $m_X(M)\in H^2(G,\mathbb C^\*)$ is the obstruction to linearising $M$ (its Mumford class). It is well defined because $\mathcal O(X)^\*=\mathbb C^\*$. The **Amitsur subgroup** is
$$\mathrm{Am}_G(X)=m_X\big(\mathrm{Pic}(X)^G\big)\subseteq H^2(G,\mathbb C^\*).$$
It vanishes if $X$ has a fixed point: on the fibre there, the cocycle becomes a coboundary of scalars.

**Theorem 5.11 (correspondence over a projective base).** Let $X$ be smooth projective, and $C$ a faithful $G$-curve, with $\mathrm{Hom}_G(\mathrm{Alb}\,X,\mathrm{Jac}\,C)=0$. Let $W\dashrightarrow X$ be equivariant and generically finite of degree $d$, and $W\dashrightarrow C$ equivariant and dominant. Then there are invariant classes $M$ on $X$ and $L$ on $C$ with
$$\deg L=e\mid d,\qquad h^0(L)\ge2,\qquad |L|\text{ base-point-free},\qquad m_C(L)=-m_X(M)\in\mathrm{Am}_G(X).$$
In particular $d\ge\tilde\mu(C)\ge\operatorname{gon}(C)$, and $L$ is linearised if $\mathrm{Am}_G(X)=0$.

*Proof.*
1. **The divisor.** Let $Z\subset X\times C$ be the closure of the image of $W$. It is an invariant irreducible divisor, so $\mathcal O(Z)$ is canonically linearised. Its degree $e$ over $X$ divides $d$, since $\mathbb C(X)\subseteq\mathbb C(Z)\subseteq\mathbb C(W)$.
2. **Splitting.** The sequence $0\to\mathrm{Pic}\,X\oplus\mathrm{Pic}\,C\to\mathrm{Pic}(X\times C)\to\mathrm{Hom}(\mathrm{Alb}\,X,\mathrm{Jac}\,C)\to0$ is canonical and equivariant. The image of the invariant class $[Z]$ lies in $\mathrm{Hom}_G=0$, so $\mathcal O(Z)\cong M\boxtimes L$, with $M$ and $L$ invariant classes.
3. **Moving and base-point-free.** The fibres $Z_x\in|L|$ move, because $Z$ dominates $C$. A point $p$ in every $Z_x$ would give $Z=X\times\{p\}$, which does not dominate $C$.
4. **Obstructions.** On the projective $X\times C$, Mumford classes are defined and add: $0=m(\mathcal O(Z))=m_X(M)+m_C(L)$. $\square$

No fixed divisor is removed, so $e\mid d$ holds for $L$ itself.

**Example 5.12 (open bases: units matter).** Theorem 5.11 must be applied on a projective model. On an open base $B$, triviality of $\mathrm{Pic}(B)$ does not force linearisation.
1. **On $\mathbb P^1$.** Let $G=C_2^2$ act on $C=\mathbb P^1$ by $z\mapsto-z$ and $z\mapsto1/z$, and put $W=B=\mathbb P^1\smallsetminus\{0,\infty,\pm1\}\subset C$.
   - This is a cover of degree 1 over a base with $\mathrm{Pic}(B)=0$.
   - Yet $L=\mathcal O(1)$ is not linearised: the lifts to $SL_2$ anticommute. Here $\mu(G)=2$.
2. **Inside the project.** Let $V\subset H^0(L_{60})$ be the $\mathbf 6$ and $I\subset\mathbb P(V)\times C$ the incidence variety. Let $B$ be $\mathbb P(V)$ minus a $G$-orbit of hyperplanes.
   - $B$ is rational with $\mathrm{Pic}(B)=0$.
   - $I|_B\to B$ is an irreducible equivariant cover of degree 60 dominating $C$.
   - $L_{60}$ is not linearised. So $90\mid d$ and $d\ge90$ both fail over such bases.

*The mechanism.* On $B$ the isomorphisms $g^\*M\cong M$ are fixed only up to $\mathcal O(B)^\*$. So the obstruction lives in $H^2(G,\mathcal O(B)^\*)$, and the old argument showed only that $m_C(L)$ dies there. The kernel of $H^2(G,\mathbb C^\*)\to H^2(G,\mathcal O(B)^\*)$ is the image of $H^1(G,\mathcal O(B)^\*/\mathbb C^\*)$. In example 2, take lifts $R_g$ with cocycle $\alpha$ and a linear form $\ell$; the units $u_g=\ell(R_gv)/\ell(v)$ satisfy $u_g(hv)\,u_h(v)=\alpha(g,h)\,u_{gh}(v)$.

**Corollary 5.13 (when $L$ is linearised).** In Theorem 5.11, $L$ is linearised in either of these cases:
- $\mathrm{Am}_G(X)=0$, for instance if $X$ has a fixed point;
- on an invariant open $B\subseteq X$, either $B$ has a fixed point, or every invariant line bundle on $B$ is linearisable and $\mathcal O(B)^\*=\mathbb C^\*$. More generally, the second condition can be weakened to: $H^2(G,\mathbb C^\*)\to H^2(G,\mathcal O(B)^\*)$ is injective.

*Proof.* On $B\times C$ the relation reads $\delta_B(M|_B)+\iota(m_C(L))=0$ in $H^2(G,\mathcal O(B)^\*)$. A fixed point $b$ linearises $L$ directly: restrict $\mathcal O(Z)$ to $\{b\}\times C$. $\square$

**Theorem 5.14 (Amitsur kernel; Hassett–Tschinkel).** For smooth projective $X$ with generically free action,
$$\mathrm{Am}_G(X)=\ker\big(H^2(G,\mathbb C^\*)\to H^2(G,\mathbb C(X)^\*)\big)=\ker\big(H^2(G,\mathbb C^\*)\to\mathrm{Br}(\mathbb C(X)^G)\big).$$
So $\mathrm{Am}_G$ is a stable equivariant birational invariant.

*Proof.*
1. **The kernel.** Use the sequences $1\to\mathbb C^\*\to\mathbb C(X)^\*\to P\to1$ and $0\to P\to\mathrm{Div}\to\mathrm{Pic}\to0$. Hilbert 90 gives $H^1(G,\mathbb C(X)^\*)=0$, and $\mathrm{Div}$ is a permutation module, so $H^1(G,\mathrm{Div})=0$. Hence the kernel is the image of $\mathrm{Pic}(X)^G$.
2. **Stable invariance.** Adding variables with trivial action changes nothing, since $\mathrm{Br}(K)\hookrightarrow\mathrm{Br}(K(t))$. $\square$

**Theorem 5.15 (stable compression formula).** Define two quantities:
- $\mu_A(C)=\min\{\deg L: L\in\mathrm{Pic}(C)^G,\ h^0(L)\ge2,\ m_C(L)\in A\}$ for a subgroup $A\subseteq H^2(G,\mathbb C^\*)$;
- $c^{\rm st}_X(C)$, the least degree of an irreducible equivariant cover dominating $C$ whose base is equivariantly birational to $X\times\mathbb P^N$, with trivial action on $\mathbb P^N$.

Under the hypotheses of Theorem 5.11, with $X$ generically free,
$$c^{\rm st}_X(C)=\mu_{\mathrm{Am}_G(X)}(C).$$

*Proof.*
- **Lower bound.** Apply Theorem 5.11 on a smooth projective model; $\mathrm{Alb}$ and $\mathrm{Am}_G$ are stable invariants.
- **Upper bound.**
  1. Take a minimiser $L$. It is base-point-free, because removing an invariant fixed divisor keeps the multiplier.
  2. Put $V=H^0(L)$, and choose $M$ on $X$ with $m_X(M)=-m_C(L)$. Then $M\otimes V$ is a linearised vector bundle.
  3. By the no-name lemma, $X\times\mathbb P(V)=\mathbb P_X(M\otimes V)$ is equivariantly birational to $X\times\mathbb P^{n-1}$.
  4. Over it, $X\times I_L$ is an irreducible cover of degree $\deg L$ dominating $C$. Here $I_L$ is the incidence variety of $|L|$, a projective bundle over $C$. $\square$

**Corollary 5.16 ($(2,4,7)$ targets).** Let $A=\mathrm{Am}_{A_7}(X)\subseteq\mathbb Z/6$. Corollary 7.2, Theorem 7.4, Theorem 7.5 and Proposition 5.8 give:

| $A$ | lattice of $e=\deg L$ | $c^{\rm st}_X(C)$ |
|---|---|---|
| $0$ | $90\mathbb Z$ | 90 |
| order 2 | $45\mathbb Z$ (45 has no sections) | 90 |
| order 3 | $30\mathbb Z$ (30 has no sections) | 60 |
| $\mathbb Z/6$ | $15\mathbb Z$ | 60 |

1. **Linear base.** It has a fixed point, so $A=0$, and a fixed target needs $90\mid d$, as in §5.1.
2. **The quadratic accessory.**
   - $B=\{t^2=q(v)\}\smallsetminus\{0\}$ has as functions the graded ring of the cone, with $\mathbb C$ in degree 0, so $\mathcal O(B)^\*=\mathbb C^\*$.
   - $\mathrm{Pic}(B)=\mathrm{Pic}(Q^5)/\langle\mathcal O(1)\rangle=0$.
   - By Corollary 5.13, a further compression still needs $d\ge90$.
3. **Any rational base.** $d\ge60$, and this is sharp (Example 5.12.2). Onto any faithful $A_7$-curve, $d\ge\gamma(A_7)\ge25$.

**Theorem 5.17 (Amitsur growth).** Let $Y\dashrightarrow X$ be dominant, equivariant and generically finite of degree $d$, between generically free smooth projective $G$-varieties. Then
$$\mathrm{Am}_G(X)\subseteq\mathrm{Am}_G(Y),\qquad d\cdot\mathrm{Am}_G(Y)\subseteq\mathrm{Am}_G(X).$$

*Proof.* The generic torsor of $Y$ is the base change of that of $X$ along $F=\mathbb C(Y)^G\supseteq K=\mathbb C(X)^G$, with $[F:K]=d$. Apply Theorem 5.14: restriction gives the first inclusion, and $\mathrm{cor}\circ\mathrm{res}=d$ gives the second. $\square$

So the $\ell$-primary parts agree for $\ell\nmid d$. Take a connected-monodromy tower starting from a linear base, and suppose the Albanese hypothesis holds throughout.
- An order-3 class can first appear only at a step of degree divisible by 3, and an order-2 class only at an even step.
- A prefix of degree prime to 3 keeps the threshold 90 for the next compression. This is a lower bound, not 90-divisibility.

**Corollary 5.18 (the Brauer obstruction is cheap).** For the generic $A_7$-torsor over a linear base, a multiplier of order $r\in\{2,3,6\}$ has Brauer index exactly $r$.

*Proof.* Twisting $\mathbb P(V)$ gives a Severi–Brauer variety, so the index divides every projective degree with that multiplier. The gcds of the lists in §7.1 are 2, 3 and 6. The period is $r$, because $\mathrm{Am}=0$. $\square$

So a quadratic and then a cubic extension kill all of $H^2(A_7,\mathbb C^\*)$ while keeping connected full monodromy: $A_7$ has no subgroup of index $\le6$, so an extension of degree $\le6$ is linearly disjoint. Schur–Brauer obstructions alone therefore cannot exclude one-variable towers. This is consistent with $\mu=90$: a cubic step followed by a degree-60 compression has total degree 180.

**Towers: what remains.**
- For a correspondence $Z$ over $X$, the map $x\mapsto[Z_x]$ is an affine map $X\to\mathrm{Pic}^e(C)$. Its linear part $u_Z:\mathrm{Alb}\,X\to\mathrm{Jac}\,C$ is equivariant.
- If $u_Z=0$, Theorem 5.15 applies.
- If not, no bound in terms of gonality or Amitsur subgroups can hold: $X=C$ with $Z=\Delta$ has $e=1$.
- The open problem is which affine maps into $W_e(C)$, lifting to $\mathrm{Sym}^eC$, arise from the actual intermediate bases of a tower.

For a $(2,4,7)$ target, $H^1(C)$ contains none of $1,6,14_a,14_b$. So $u_Z=0$ whenever $H^1(X)$ involves only these.

## 5.6 Trust base and open questions

**Cited.**
- ATLAS subgroup lists.
- Castelnuovo's bound; Halphen–Gruson–Peskine (once).
- Holomorphic Lefschetz (for $\mu\le90$).
- Farb–Wolfson (only for $a\le60$).
- For §5.5: the no-name lemma, the divisorial-correspondence sequence for $\mathrm{Pic}(X\times C)$, and index-period facts for central simple algebras. Theorem 5.14 is in Hassett–Tschinkel.

**Open.**
- **Towers with $u_Z\ne0$** (§5.5). This is the remaining obstacle to $\mathrm{RD}(A_7)>1$ along these lines. The Amitsur part alone is cheap to remove (Corollary 5.18).
- Is $a(A_7)=60=\tilde\mu(C)$ a coincidence? The degree-60 class gives $2L_{60}\sim5\cdot24$ points through the Klein subgroups (Prop. 7.9), while $a$ uses $15\cdot4$ through the same subgroups.
- The arithmetic questions of §6.5.
