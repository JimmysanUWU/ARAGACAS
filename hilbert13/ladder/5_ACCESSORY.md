# Chapter 5. Accessory irrationalities: $a(A_7)=60$ and $\mu(A_7)=90$

**Definitions.**
- $K=\mathbb C(V)^{G}$ for a faithful linear model $V$. An *accessory* is a finite extension $F/K$.
- After the accessory, the monodromy of the torsor may drop to a subgroup $H$.
- A *curve compression* is a faithful $H$-stable curve field $E\subset LF$ with $LF=FE$.
- **$a(G)$** is the least accessory degree that allows a curve compression.
- **$\mu(G)$** is the least degree of a $G$-*linearised*, base-point-free line bundle with $h^0\ge2$ on a faithful connected $G$-curve.

**Theorem 5.1.**
1. $a(A_7)=60$; that is, $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$ and $\mathrm{ed}_{\mathbb C}(A_7;\le60)=1$.
2. $\mu(A_7)=90$: compression with connected full $A_7$-monodromy needs exactly accessory degree 90.

**Credits.**
- Part 1, the correspondence theorem and the degree lattice: GPT (DAY 2, §§1–3; `gpt/DAY_2_A7_Mathematical_Reference.pdf`). All are re-derived here.
- $\mu\ge72$ and the power-sum curve: GPT Round 5 (`gpt/A7_Round5_Advances.pdf`).
- $\mu\le90$ and the exclusion of 72 and 84 are ours.

Status tags as in `README.md`; the computations are in `verify_accessory60.py`, `verify_mu90_exact.py` and `equivariant_rr.py`.

## 5.1 From compressions to line bundles

**Theorem 5.2 (correspondence; DAY 2 Thm 1.2) [P].** Take an equivariant degree-$d$ cover $W\dashrightarrow V$ and a dominant equivariant map
$W\dashrightarrow C$ to a faithful $G$-curve. Then $\mu(G)\le d$.

*Proof.*
1. The closure $Z$ of the joint image in $V\times C$ is an invariant Cartier divisor, so $\mathcal O(Z)$ is canonically linearised.
2. $\mathrm{Pic}(V\times C)=\mathrm{Pic}(C)$, so all fibres have one class: an effective divisor of degree $e\mid d$ that moves.
3. Restrict to $\{0\}\times C$, with $0\in V$ fixed, and remove the invariant fixed part. $\square$

**Theorem 5.3 (subgroup induction; DAY 2 Thm 1.5) [P].** $a(G)=\min_{H\le G}[G:H]\,\mu(H)$.

*Proof.*
- *Lower bound.* Let $H$ be the monodromy after the accessory. Then $[F:K]=[G:H]e$, and Theorem 5.2 applies to $H$.
- *Upper bound.* Induce from $H$ and use the incidence variety (DAY 2 Prop. 1.4). The upper bound holds also for disconnected base change (Round 5 Prop. 1.1). $\square$

**Lemma 5.4 (degree lattice; DAY 2 Thm 2.1) [P].** A linearised bundle on a faithful $G$-curve $C$ has degree in $N_G(C)\mathbb Z$, where
$N_G(C)=|G|/\mathrm{lcm}(e_i)$ and the $e_i$ are the inertia orders.

*Proof.* Hilbert 90 gives an invariant rational section. Its divisor is a sum of orbits, of sizes $|G|/e_i$ and $|G|$. $\square$

**Lemma 5.5 (linearised Castelnuovo; DAY 2 Prop. 2.2) [P].** Let $G$ be simple, and let $L$ be minimal of degree $n<|G|$. Then the complete series maps
$C$ birationally onto a faithful $G$-curve $Y\subset\mathbb P^r$, with $r\ge q(G)-1$ and $g(Y)\le\pi(n,r)$.
Here $q(G)$ is the least degree of a nontrivial representation; $q(A_7)=6$.

*Proof.* The kernel on the image has order $\le n<|G|$, so it is trivial. The sections form a nontrivial representation. Minimality makes the map birational. $\square$

## 5.2 $a(A_7)=60$

**Proposition 5.6 [P][X].** $\mu(A_7)\ge60$.

*Proof.*
1. Inertia orders divide 420, so $6\mid N$, and $n<60$ forces $n\le54$.
2. By Lemma 5.5, $g(Y)\le\pi(54,5)=325$.
3. In genus 136–325 the only lattice-compatible cases are
   - $(2,5,7)$ with $n=36$, where $\pi=136<199$;
   - $(3,4,5)$ with $n=42$, where $\pi=190<274$.

   Both are contradictions. $\square$

**Subgroups of index $<60$** (ATLAS):

| $H$ | index | $\mu(H)\ge$ | reason | product |
|---|---|---|---|---|
| $A_6$ | 7 | 9 | $\pi(8,4)=5<10=g_{\min}(A_6)$ | 63 |
| $L_2(7)$ | 15 | 4 | $\pi(3,2)=1<3$ | **60** |
| $S_5$ | 21 | 3 | a hyperelliptic involution is central, and $S_5\not\subset\mathrm{PGL}_2$ | 63 |
| $(A_4\times C_3){:}2$ | 35 | 2 | $C_3^2\not\subset\mathrm{PGL}_2$ | 70 |
| $A_5$ | 42 | 2 | no 2-dimensional representation | 84 |

**Upper bound.** $L_2(7)$ with the canonical bundle of the Klein quartic gives $15\cdot4=60$. This proves Theorem 5.1(1). The upper bound
uses the Farb–Wolfson independence of the linear model; the lower bound needs nothing beyond §5.1.

**Small groups (Round 5 §3) [P].**
- $a(A_5)=2$. This is Kronecker–Klein's square-root accessory for the quintic.
- $a(L_2(7))=4$.
- $a(A_6)=12$:
  - $6\mid n$ and $\pi(6,4)=2$;
  - there are no subgroups of order 40 or 45;
  - the subgroup of order 36 is not a Möbius group.

## 5.3 $72\le\mu(A_7)\le90$

**Proposition 5.7 (GPT Round 5, Thm 2.2) [P][X].** $\mu(A_7)\ge72$.

*Proof.*
- By Prop. 5.6 and the lattice, $n\in\{60,66\}$, with signatures $(2,6,7)$ and $(3,4,7)$.
- $\pi(60,8)=220$ forces the section representation to be exactly the standard 6.
- Its $A_6$-fixed vector has a zero divisor of degree 60, but $A_6$-orbits have size $\ge120$ (resp. 90). $\square$

Every degree in $\{78,102,114\}$ is excluded too.
- 78 and 102 have no lattice signature.
- For 114, the signature is $(3,4,5,7)$ ($g=1354$). $\pi(114,6)=1221$ forces $W=6$, and then $\prod x_j$ is an invariant section of degree 798, which is not a sum of orbit sizes $\{360,504,630,840,2520\}$.

So $\mu(A_7)\in\{72,84,90,96,108,120\}$. The power-sum curve $\{\sum x_j^k=0,\ k\le5\}\subset\mathbb P^6$ (genus 481, signature $(3,7,7)$, $\mathcal O(1)$ of degree 120) shows $\mu\le120$.

**Proposition 5.8 [P][X].** $\mu(A_7)\le90$. On every $(2,4,7)$ curve, $B+T=D_2-3D_4+2D_7$ is a linearised moving bundle of degree 90.

*Proof.*
1. **Linearised Picard group.** It is $\langle D_2,D_4,D_7\mid2D_2=4D_4=7D_7\rangle\cong\mathbb Z\oplus\mathbb Z/2$. Its degree-90 classes are $B=2D_7-D_4$ and $B+T$.
2. **Holomorphic Lefschetz** (Atiyah–Bott). For $g\ne1$,
   $\sum_q(-1)^q\mathrm{tr}(g\mid H^q(\mathcal O(D)))=\sum_{p\in\mathrm{Fix}(g)}a_p^{k_p}/(1-a_p^{-1})$.
   Here $a_p$ is the rotation of $g$ at $p$ and $k_p=\mathrm{mult}_pD$.
3. **The result.** $\chi_{A_7}(B+T)=-6+10-14_a-14_b-21$ (the $10$ becomes $\overline{10}$ under the other orientation), so $h^0(B+T)\ge10$.
4. **No fixed part.** A fixed part would be an invariant divisor of degree $\le90<360=$ the least orbit size. $\square$

*Validation of conventions.*
- $\chi(K)$ gives $H^0(K)=10+2\cdot\overline{10}+15+21+2\cdot35$, matching Chevalley–Weil (§1.4).
- The other conventions give non-integral multiplicities.

For comparison, $\chi(B)=-\overline{10}-35$ decides nothing, and $\chi(2B)=6-\overline{10}+14_a+14_b+21$.

## 5.4 $\mu(A_7)\ne72,84$

**Candidates [X].** Let $L$ realise $\mu=n$ on $C$. Then $g(C)\le\pi(n,h^0-1)$ and $N(C)\mid n$. Combined with the seven-sheeted test and
generation, this leaves:

| $n$ | signature (genus) | $h^0\le$ |
|---|---|---|
| 72 | $(2,5,7)$ (199) | 12 |
| 72 | $(3,5,7)$, $(4,5,7)$, $(5,5,7)$ | 7, 6, 6 |
| 84 | $(3,4,5)$ (274), $(3,5,6)$ (379) | 12, 10 |
| 84 | $(4,5,6)$, $(5,5,6)$, $(5,6,6)$, $(5,6,7)$, $(2,2,3,5)$, $(2,2,5,6)$, $(2,3,3,5)$ | $\le8$ |

Every orbit has size $>n$, so $W=H^0(L)$ has no trivial summand. Hence $W\in\{6,\,10,\,\overline{10},\,6\oplus6\}$.

**Lemma 5.9 (orbit semigroup) [P].** Let $W\subseteq H^0(L)$ be a subrepresentation and $f\in\mathrm{Sym}^kW$ invariant. If $kn$ is not a nonnegative
combination of orbit sizes, then $f$ vanishes on $\varphi_W(C)$.

*Proof.* $f|_C$ is an invariant section of $L^k$, so its divisor is a sum of orbits. $\square$

Lemma 7.3 sharpens this by also using the local characters at branch points.

**Case $W\supseteq6$.** Let $x_1,\dots,x_7$ be the coordinates, with $\sum x_j=0$, and let $p_k$ be the power sums.
- **$n=72$.**
  - $p_2,p_3,p_4,p_6$ vanish, because $144,216,288,432$ are not sums of orbit sizes. So $e_1=e_2=e_3=e_4=e_6=0$.
  - The image lies in $X=\{\text{roots of }t^7-ut^2-v\}$. $X$ is irreducible: $t^7-t^2$ has monodromy $S_7$.
  - $X$ is generically reduced, of degree $1\cdot2\cdot3\cdot4\cdot6=144>72$. So the image cannot be a component.
- **$n=84$.**
  - $p_2,p_3,p_4,p_7$ vanish, and Newton's identities give $p_7=7e_7$. So $\prod x_j=0$ on $C$.
  - By transitivity the image would lie in every coordinate hyperplane.

**Spin model for the 10.** $A_7\subset SO(6)$ lifts to $2.A_7\subset SU(4)$. On the half-spin representation $V_4$, $\wedge^2V_4=6$ and
$\mathrm{Sym}^2V_4=\overline{10}$, which is irreducible (exact character norm 1).

So a curve with $W=10$ is a family of symmetric $4\times4$ matrices $Q_p$. It occurs only for $(2,5,7)$ at 72 and for $(3,4,5)$, $(3,5,6)$ at 84.

**Lemma 5.10 (parity) [P].** Let $c$ fix $p$, with lift $s$ to $2.A_7$. Then $Q_p$ lies in the $\lambda$-eigenspace of $c$. If $\lambda\ne\mu^2$ for every eigenvalue $\mu$ of $s$,
then every matrix in that eigenspace has even rank.

*Proof.* In an eigenbasis of $s$, $Q_{ab}\ne0$ only if $\mu_a\mu_b=\lambda$. Then $\mu\mapsto\lambda/\mu$ is a fixed-point-free involution, so $Q$ is a sum of
blocks $\left(\begin{smallmatrix}0&M\\M^T&0\end{smallmatrix}\right)$. $\square$

**Case $W\in\{10,\overline{10}\}$.** By Lemma 5.9, $\det Q_p\equiv0$, because $4n$ is below the least orbit size. So the generic rank $r$ is 1, 2 or 3.
- **$r=2$.** $\wedge^2Q_p$ defines an equivariant map to $\mathbb P(\wedge^2V_4)=\mathbb P(6)$. Its bundle $N$ has $N^2=L^2(-B')$, so $\deg N=n$ and $W\supseteq6$ for $N$.
  This is the previous case.
- **$r=3$.** The adjugate gives a kernel map with $\kappa^*\mathcal O(2)=L^3(-B'')$, so $\deg B''\le3n$.
  - In every case, some branch eigenspace contains no rank-3 matrix with $\det=0$.
  - So $B''$ contains an orbit of size 1260, 630 or 420. Each exceeds $3n$, a contradiction.
- **$r=1$.** Then $Q_p=v_p\otimes v_p$.
  - For $(2,5,7)$ and $(3,4,5)$, Lemma 5.10 applies at the $e=2$, resp. $e=4$, branch.
  - For the two classes of $(3,5,6)$ without a parity branch, $\varphi=v_2\circ\psi$ with $\psi:C\to\mathbb P^3$ of degree 42.
    - $\psi$ is birational, since $\pi(21,3)=90<136$.
    - $\psi(C)$ lies on no quadric, since $\mathrm{Sym}^2V_4$ is irreducible.
    - Halphen (Gruson–Peskine) gives $g\le42\cdot39/6+1=274<379$.

So $\mu\notin\{72,78,84\}$, and Proposition 5.8 gives $\mu(A_7)=90$. $\square$

**Exactness of the rank steps [X]** (`verify_mu90_exact.py` → `mu90_exact_output.txt`).
- **Lift spectra are forced.** $\{\mu_a\mu_b\}_{a<b}$ is the spectrum of $g$ on the 6, and $\prod\mu_a=1$. The solution is unique up to $s\mapsto-s$, except that 7-cycles also allow $V_4$ versus $V_4^*$.
- **Eigenspaces are pattern spaces.** In an eigenbasis of $s$, $Q\mapsto sQs^T$ is diagonal. So each branch eigenspace is the space of symmetric matrices supported on a 0/1 pattern.
- **Rank tests.**
  - Rank 1 occurs iff the pattern has a loop.
  - Rank 3 is decided by exact polynomial algebra: the determinant, the $3\times3$ minors, and radical membership via factorisation over $\mathbb Q$.
- **Coverage.** All 84 (curve, class) cases close, under both orientation conventions.

## 5.5 Beyond fixed points: the twisted correspondence [P]

Theorem 5.2 restricts $\mathcal O(Z)$ to a fixed point. Round 5 (Prop. 4.2) showed that fixed points do not survive a quadratic
accessory $t^2=q(v)$: its total space keeps full $A_7$ but has no fixed point. The argument does not actually need one.

**Theorem 5.11.**
- *Setting.* $B$ is a smooth irreducible $G$-variety with a smooth $G$-compactification $X$ such that $\mathrm{Hom}_G(\mathrm{Alb}\,X,\mathrm{Jac}\,C)=0$. This holds for instance if $X$ is rational, or if $H^1(X)$ and $H^1(C)$ share no irreducible $G$-constituent. $W\to B$ is a $G$-equivariant finite cover of degree $d$, and $W\dashrightarrow C$ is a dominant $G$-map to a faithful $G$-curve.
- *Conclusion.* $C$ carries a $G$-invariant class $L$ with $h^0(L)\ge2$ and $\deg L=e\mid d$. Its Mumford class is the negative of that of an invariant class on $B$. Hence:
  1. $d\ge e\ge\tilde\mu(C)\ge\operatorname{gon}(C)$, where $\tilde\mu(C)$ is the least degree of an invariant class with $h^0\ge2$;
  2. if every invariant line bundle on $B$ is linearisable (for instance $\mathrm{Pic}(B)=0$, or $B$ has a fixed point), then $L$ is linearised and $d\ge\mu(G)$.

*Proof.*
1. **The divisor.** The closure $Z$ of the image of $W$ in $B\times C$ is an invariant divisor, so $\mathcal O(Z)$ is linearised. Since $W\to Z\to B$, the degree $e$ of $Z\to B$ divides $d$.
2. **Splitting.** $\mathrm{Pic}(X\times C)=\mathrm{Pic}(X)\oplus\mathrm{Pic}(C)\oplus\mathrm{Hom}(\mathrm{Alb}\,X,\mathrm{Jac}\,C)$, $G$-equivariantly. The class of the invariant divisor $\bar Z$ is $G$-invariant, so its correspondence part lies in $\mathrm{Hom}_G=0$. Restricting to $B\times C$, $\mathcal O(Z)\cong p_B^\*M\otimes p_C^\*L$.
3. **The class $L$.** Every fibre $Z_b$ lies in $|L|$, and the fibres move because $W\dashrightarrow C$ is dominant. So $h^0(L)\ge2$.
4. **Invariance.** $\sigma^\*Z=Z$ gives $\sigma^\*L\cong L$. Mumford classes add under $\otimes$, and $\mathcal O(Z)$ has none, so $m(L)=-m(M)$.
5. **Fixed part.** Removing the invariant fixed part of $|L|$ changes neither the Mumford class nor the bound. $\square$

**Corollary 5.12.**
1. **Fixed target.** Compression onto a fixed $(2,4,7)$ curve with connected full $A_7$-monodromy needs $90\mid d$, because linearised degrees there lie in $90\mathbb Z$ (Lemma 5.4). This answers DAY 2 §9.
2. **After a quadratic accessory.** For $B=\{t^2=q(v)\}\smallsetminus\{0\}$, which is rational with $\mathrm{Pic}(B)=0$, a further compression to a curve still needs $d\ge\mu(A_7)=90$. So the obstruction of Round 5 Prop. 4.2 is only apparent: what persists is *linearisability*, not the fixed point.
3. **Any rational base.** For any rational $B$:
   - onto a $(2,4,7)$ curve, $d\ge\tilde\mu=60$ with $15\mid e$ (Theorems 7.4–7.5);
   - onto any faithful $A_7$-curve, $d\ge\gamma(A_7)\ge25$.

   So gonality, which plays no role in $a(A_7)$, governs compressions over fixed-point-free rational bases.

**For towers**, the hypothesis needs only to hold stage by stage. For a $(2,4,7)$ target, $H^1(C)=3(10+\overline{10})+2\cdot15+2\cdot21+4\cdot35$ (Chevalley–Weil), which contains no $1,6,14_a,14_b$. So the bounds hold over any base whose $H^1$ involves only these four representations.
What remains is a base whose $H^1$ shares constituents with $H^1(C)$. There the fibre classes $[Z_b]$ can move along an equivariant map $\mathrm{Alb}\,X\to\mathrm{Jac}\,C$.

## 5.6 Trust base and consequences

**Cited.**
- ATLAS subgroup lists;
- Castelnuovo's bound;
- Halphen–Gruson–Peskine (used once);
- holomorphic Lefschetz (for $\mu\le90$);
- Farb–Wolfson (only for $a(A_7)\le60$).

**Consequences.**
- **The essential-dimension question is settled.** $\mathrm{ed}(A_7;\le59)>1$ needs no gonality input, and the threshold is exactly 60.
- **Gonality returns over fixed-point-free bases** (Corollary 5.12): there the least compression degree onto a faithful $A_7$-curve is bounded below by
  $\gamma(A_7)\in[25,42]$ (Chapters 4 and 7).

**Open.**
- Towers through a base whose $H^1$ shares an irreducible constituent with $H^1(C)$.
- The arithmetic torsion questions (§1.4).
