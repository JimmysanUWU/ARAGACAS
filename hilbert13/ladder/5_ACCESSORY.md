# Chapter 5. Accessory irrationalities: $a(A_7)=60$ and $\mu(A_7)=90$

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
- Ours: $\mu\le90$, the exclusion of 72 and 84, and §5.5.

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

## 5.5 Beyond fixed points: the twisted correspondence [P]

Theorem 5.2 restricts $\mathcal O(Z)$ to a fixed point. Fixed points do not survive a quadratic accessory $t^2=q(v)$ (Round 5 Prop. 4.2), but the argument does not need them.

**Theorem 5.11.**
- *Setting.*
  - $B$ is a smooth irreducible $G$-variety with a smooth $G$-compactification $X$ such that $\mathrm{Hom}_G(\mathrm{Alb}\,X,\mathrm{Jac}\,C)=0$. This holds, e.g., if $X$ is rational, or if $H^1(X)$ and $H^1(C)$ share no $G$-constituent.
  - $W\to B$ is a $G$-equivariant finite cover of degree $d$, with a dominant $G$-map $W\dashrightarrow C$ to a faithful $G$-curve.
- *Conclusion.* $C$ carries a $G$-invariant class $L$ with $h^0\ge2$ and $\deg L=e\mid d$, whose Mumford class is minus that of an invariant class on $B$. Hence:
  1. $d\ge e\ge\tilde\mu(C)\ge\operatorname{gon}(C)$, where $\tilde\mu(C)$ is the least degree of an invariant class with $h^0\ge2$;
  2. if every invariant line bundle on $B$ is linearisable (e.g. $\mathrm{Pic}(B)=0$, or $B$ has a fixed point), then $L$ is linearised and $d\ge\mu(G)$.

*Proof.*
1. The closure $Z\subset B\times C$ of the image of $W$ is invariant, so $\mathcal O(Z)$ is linearised. Its degree $e$ over $B$ divides $d$.
2. $\mathrm{Pic}(X\times C)=\mathrm{Pic}(X)\oplus\mathrm{Pic}(C)\oplus\mathrm{Hom}(\mathrm{Alb}\,X,\mathrm{Jac}\,C)$, equivariantly. $[\bar Z]$ is invariant, so its correspondence part lies in $\mathrm{Hom}_G=0$. Hence $\mathcal O(Z)\cong p_B^\*M\otimes p_C^\*L$ on $B\times C$.
3. The fibres $Z_b\in|L|$ move, because $W\dashrightarrow C$ is dominant. So $h^0(L)\ge2$.
4. $\sigma^\*L\cong L$, and Mumford classes add, so $m(L)=-m(M)$. Removing the invariant fixed part changes neither. $\square$

**Corollary 5.12.**
1. **Fixed target.** Compression onto a fixed $(2,4,7)$ curve with connected full monodromy needs $90\mid d$ (Lemma 5.4). This answers DAY 2 §9.
2. **After a quadratic accessory.** $B=\{t^2=q(v)\}\smallsetminus\{0\}$ is rational with $\mathrm{Pic}(B)=0$. So a further compression still needs $d\ge90$. What persists is *linearisability*, not the fixed point.
3. **Any rational base.** $d\ge\tilde\mu=60$ with $15\mid e$ onto a $(2,4,7)$ curve (Theorems 7.4–7.5), and $d\ge\gamma(A_7)\ge25$ onto any faithful $A_7$-curve.

So gonality, which plays no role in $a(A_7)$, governs compressions over fixed-point-free rational bases.

**Towers.** The hypothesis is needed stage by stage. For a $(2,4,7)$ target, $H^1(C)$ contains none of $1,6,14_a,14_b$. So the bounds hold over any base whose $H^1$ involves only these. What remains is a base whose $H^1$ shares constituents with $H^1(C)$, where $[Z_b]$ can move along an equivariant map $\mathrm{Alb}\,X\to\mathrm{Jac}\,C$.

## 5.6 Trust base and open questions

**Cited.**
- ATLAS subgroup lists.
- Castelnuovo's bound; Halphen–Gruson–Peskine (once).
- Holomorphic Lefschetz (for $\mu\le90$).
- Farb–Wolfson (only for $a\le60$).

**Open.**
- Towers through a base whose $H^1$ shares a constituent with $H^1(C)$. This is the remaining obstacle to $\mathrm{RD}(A_7)>1$ along these lines.
- Is $a(A_7)=60=\tilde\mu(C)$ a coincidence? The degree-60 class gives $2L_{60}\sim5\cdot24$ points through the Klein subgroups (Prop. 7.9), while $a$ uses $15\cdot4$ through the same subgroups.
- The arithmetic questions of §6.5.
