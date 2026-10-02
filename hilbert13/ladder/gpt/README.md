# GPT's documents: a digest

GPT has supplied five PDFs, all in this folder, listed here in the order they were written. This page states what each one proves, in its sharpest form, and where that result now lives.

**Tags.**
- **[P]/[X]**: re-derived or recomputed in this repository (the chapters carry the proofs).
- **[G]**: GPT's proof only; stated here and not used by any chapter.
- **[S]**: superseded.

| file | date | written against | content | read for |
|---|---|---|---|---|
| `DAY_1C_Verification_Ladder.pdf` | 29 Sep | the Sol/Astra synthesis of 28 Sep | verification note and research charter | the rungs G17–T2 (§3), interfaces (§2.2) |
| `DAY_1B_Proof_Chain_Audit.pdf` | 29 Sep | commit `8b31a2a` | audit of the spectral proof of $\operatorname{gon}\ge24$ | the repair list (§4) |
| `DAY_1A_Research_Checkpoint.pdf` | 30 Sep | PR #2 head `284e1e2` | weighted base points; $\operatorname{gon}\ge25$ for $g\ge266$; the Belyi map | §§2.1, 2.5 |
| `DAY_2_A7_Mathematical_Reference.pdf` | 1 Oct | DAY 1 and PR #2 | consolidated reference; $a(A_7)=60$ | Ch. 5; §§2.1–2.6 |
| `A7_Round5_Advances.pdf` | 1 Oct | commit `1ef3f6f` | answers to Round 5: $\mu\ge72$, small groups, towers, harmonic Hersch | Ch. 3, 5; §2.4 |

DAY 2 consolidates all of DAY 1, so read DAY 2 and Round 5 first. Our questions to GPT are in `../QUESTIONS_FOR_GPT.md`; Rounds 2–7 are in git history.

## 1. Results now proved in the chapters

| result | statement | GPT source | here |
|---|---|---|---|
| ramification transport | $\operatorname{gon}(C/\tau)\ge10$. A nine-pencil saturates the involution determinant. Centraliser averaging gives $T=24[R_V]$ for every Klein $V\ni\tau$, $\langle N(V_0),N(V_1)\rangle=A_7$, and $\deg T=1296\notin15\mathbb Z$ | DAY 1C §3; DAY 2 Thm 4.2 | §6.1 [P][X] |
| the $Q_2$ group | $C_{A_7}((16)(23))\cong C_3\rtimes D_8$, not $S_4$ | DAY 1B p. 3 | §2.3, §6.4 [X] |
| correspondence | an equivariant degree-$d$ cover $W\dashrightarrow V$ dominating a faithful curve gives $\mu(G)\le d$ | DAY 2 Thm 1.2 | Thm 5.2 [P] |
| induction | $a(G)=\min_H[G:H]\,\mu(H)$, including disconnected base change | DAY 2 Thm 1.5; R5 Prop 1.1 | Thm 5.3 [P] |
| degree lattice | linearised degrees lie in $N_G(C)\mathbb Z$, $N_G=\lvert G\rvert/\mathrm{lcm}(e_i)$; invariant classes satisfy $N_G\mid s\deg$, $s=\exp H^2(G,\mathbb C^*)$ | DAY 2 Thm 2.1 | Lemma 5.4; sharp form $15\mathbb Z$ in Cor. 7.2 |
| linearised Castelnuovo | a minimal linearised series embeds in $\mathbb P^{\ge q(G)-1}$ with $g\le\pi(n,q-1)$ | DAY 2 Prop 2.2 | Lemma 5.5 [P] |
| **$a(A_7)=60$** | $\mu(A_7)\ge60$, every subgroup wall below 60, and the Klein quartic of $L_2(7)$ at $15\cdot4$ | DAY 2 Thm 3.1–3.2 | Prop 5.6, Thm 5.1 [P][X] |
| small groups | $a(A_5)=2$, $a(L_2(7))=4$, $a(A_6)=12$ | R5 Prop 3.1 | §5.2 [P] |
| $\mu(A_7)\ge72$ | degree 60 forces $W=6$, and the $A_6$-fixed vector has too few zeros | R5 Thm 2.2 | Prop 5.7 [P][X] |
| the fixed-curve test | on a $(2,4,7)$ curve, $\mu_C\in\{90,180\}$, decided by $h^0(B)$, $h^0(B+T)$ | R5 §2 | Prop 5.8: $h^0(B+T)\ge10$, so 90 |
| no fixed point after a quadratic accessory | $t^2=q(v)$ on the standard 6 | R5 Prop 4.2 | §5.5; Thm 5.11 shows linearisability persists |
| harmonic Hersch | $\lambda_ht+c=4\Theta(z,y,y)+2\Theta(z,z,y)$ from $x\times dx={*dx}$; the trial-space version | R5 Thm 5.1, 5.3 | Thm 3.4 [P], certified in §3.5 |
| bracket bound | $\kappa^2\le B^*/3$, from $\sum_{a<b}\lVert w_a\wedge w_b\rVert^2\le\frac13$ | R5 §5 | Thm 3.4 step 4 |
| non-criticality | an invariant constant-length first-eigenfunction frame would be a minimal immersion of curvature $<0$ in a sphere (Bryant) | DAY 1A p. 6; DAY 2 §7.3 | Prop 6.2 [P] |
| septic Belyi map | explicit model of $C/A_6\to C/A_7$ (§2.5) | DAY 1A p. 5; DAY 2 §8.1 | §6.5; `side_checks.py` D6 [X] |
| torsion together | the Klein difference and the lift of $(\star)$ lie in the same unique 21 and 35 copies; 42/42 and 58/70 conjugates detect them | DAY 1A p. 5; DAY 2 §8.3 | Prop 1.6, §6.5 [X] |

## 2. Results not used in the chapters [G]

### 2.1 Multigraded Castelnuovo theory (DAY 1A pp. 2–3; DAY 2 §5)

- **Weighted base points.** Take a birational pair of degree-$m$ pencils, with image $B\subset\mathbb P^1\times\mathbb P^1$, and a third pencil carried by an integral $(a,b,1)$ surface. Resolving the base cluster gives $\sum n_i^2=2ab$ and $\sum n_ir_i=(a+b-1)m=:T$. Hence
  $$g\le(m-1)^2-\frac{T^2}{4ab}+\frac T2 .$$
- **Intersection-matrix bound (DAY 2 Thm 5.1).** Let $B$ lie on an integral surface of type $(a,b,c)$ in $(\mathbb P^1)^3$, with coordinate degrees $v$, and let $M$ have rows $(0,c,b)$, $(c,0,a)$, $(b,a,0)$. Then
  $$g\le1+\tfrac12\big[v^tM^{-1}v+\textstyle\sum_i(a_i-2)m_i\big].$$
  For equal degrees $m$ the three cases are:
  - $(2,1,1)$: $m^2/2-m+1$;
  - $(2,2,1)$: $7m^2/16-m/2+1$;
  - $(2,2,2)$: $3m^2/8+1$.
- **Sharp independent-triple theorem (Thm 5.2).** Let $m\ge8$, and take three pairwise birational degree-$m$ pencils with independent products. Then $g\le\lfloor m^2/2-m+1\rfloor$, with equality for even $m=2k$ on $\lvert-kK\rvert$ of a $(2,1,1)$ del Pezzo surface. The proof:
  1. linear normality, via Petrakiev's $\pi_2$ applied in $\mathbb P^8$;
  2. $h_\Gamma(2)\ge19$;
  3. uniform position.
- **Thm 5.3: $g\ge266\Rightarrow\operatorname{gon}\ge25$** for faithful $A_7$-curves.
  - Every third pencil is dependent.
  - At least 33 thirds form a star of base clusters, and only $(m,r)=(24,23)$ survives.
  - That case gives a smooth plane curve of degree 25 and genus $276\not\equiv1\pmod3$, a contradiction.
  - *Relation to Ch. 4.* Theorem 4.1 uses the weaker $\pi(3d,7)$ and certifies the window $336\le g\le397$ spectrally (§4.4). Thm 5.3 covers that window algebraically, so it is a second, independent proof.
- **Algebraic gonality bounds.** At $g=136,169,199,211,241$ they are $18,20,21,22,23$, below the spectral bounds of Ch. 2.
- **Picard closure (DAY 2 §5.4).** Suppose $g>257$, the gonality is 18, all pairs are birational and $\lvert W^1_{18}\rvert\ge7$. Then $W^1_{18}$ is a coset of a finite 2-torsion group. An $A_7$-stable set sums to an invariant class of degree $18\cdot2^a$, which a free $C_5$ or $C_7$ excludes.
- **Harui exclusion (DAY 1C §4).** A birational $(17,17)$ model with an ordinary 16-fold point projects to a smooth plane curve of degree 18, so $\lvert\mathrm{Aut}\rvert\le6\cdot18^2<2520$.

### 2.2 Transport beyond saturation (DAY 1C §6; DAY 2 §4.2)

- **Global transport.** Join commuting involutions $\tau,\mu$ when $\tau\mu$ lies in the class. If this graph is connected and $C/\tau$ has a map of degree $r/2$, then $\mathrm{Pic}^{3cr}(C)^G\ne\emptyset$, where $r=\#\mathrm{Fix}$ and $c=\lvert C_G(\tau)\rvert$.
  - For $A_7$ the graph is connected, of degree 8, with 420 edges and 140 Klein triangles.
  - Connectivity is essential. For $A_5$ with signature $(2,2,2,3)$ in genus 6, the graph has five components, and the genus-2 involution quotients have maps of degree 3, although $5\nmid3\cdot4\cdot6$.
- **The open interface.** Above saturation ($m>r$) the determinant leaves a residual:
  $$A_\tau+\mu\cdot A_\tau=R_\mu+R_{\tau\mu}+E_{\tau,\mu},\qquad\deg E_{\tau,\mu}=2(m-r),\qquad T_\tau-T_\mu=[E_{\tau,\mu}]-[E_{\mu,\tau}]\in\mathrm{Jac}(C).$$
  A balancing theorem for these residuals would give $\operatorname{gon}(D)\ge11$ algebraically.
- **The Klein difference.** $\delta=[R_{V_0}-R_{V_1}]$ has infinite order or order divisible by 5.
- **Hurwitz corollary (DAY 2 Cor 4.4).** The genus-6 involution quotients of the three genus-14 $L_2(13)$ Hurwitz curves have gonality exactly 4 (since $13\nmid216$).
- **Retired.** Averaging over all of $A_7$ multiplies away the detecting prime 5; the 24-element centraliser average keeps it.

### 2.3 The torsion/norm model of a nine-pencil (DAY 2 §6) [S]

This was the pre-transport reduction of $\operatorname{gon}(D)=9$ to a finite collision test on the genus-22 curve $C/\langle\tau,(567)\rangle$. It is superseded by §6.1. Two lessons remain:
- the moving case has three compatible divisors, not two;
- a smooth $(9,9)$ curve with the same $C_2\times S_3$ data exists (cf. Theorem 1.7), so the full $A_7$ is indispensable.

### 2.4 Spectral and conformal (DAY 1B; DAY 2 §7; R5 §5)

- **Conformal gain.**
  - $\Lambda_G(C)=\sup_h\lambda_1(hg)\int h\le8\pi\operatorname{gon}(C)$ over invariant conformal factors.
  - The gain is strict: separate the cone of invariant PSD frames. Its size is unmeasured; the ceiling in §6.3 is 0.2%.
  - An equivariant maximiser exists, possibly conical (Vinokurov, Cor. 1.4; all orbits have size $\ge360$).
- **A diagnostic.** If a degree-$d$ map exists and $\zeta_d=(8\pi d/A-\lambda_1)/(\lambda_{\rm next}-\lambda_1)$, then $\frac1A\int\lvert F-1\rvert\le2\zeta_d+2\sqrt{\zeta_d}$.
- **A Bessel bound.** $j_{1,1}>19/5$, proved from the alternating series. The certificate uses a cited decimal; with $19/5$ instead it still gives $\lambda_1\ge0.340893$ (Ch. 2).
- **A sufficient box.** R5 Thm 5.3 asks for $\lambda_h\le0.351$, $\rho\le0.005$, $B_h\le2\cdot10^{-4}$ and $\Gamma_h\le1.9$.
  - Our certified values (§3.5) miss the box in $B_h$ ($4.45\cdot10^{-4}$) and $\Gamma_h$ (2.27).
  - But $\rho\le6\cdot10^{-4}$ is far smaller than 0.005, and (3.1) holds directly.
  - R5 Cor 5.2's box ($B^*\le3\cdot10^{-4}$, $\Gamma\le2$, on the true eigenspace) is not used.
- **Two alternatives Ch. 3 does not need.**
  - The complementary-energy upper bound (R5 (5.4)) for $\langle f,(\Delta-r)^{-1}f\rangle$ via equilibrated fluxes. Ch. 3 uses the Jacobian flux (Lemma 3.5) instead.
  - The frame correction $L^TSL$, not $LSL^T$, which applied to the retired `topo_hersch.py`. `hh_certify.py` bounds $\lvert\nabla\Psi\rvert^2$ directly.

### 2.5 Arithmetic (DAY 1A p. 5; DAY 2 §8)

- **The septic Belyi map** [X].
  - **Definition.** Let $21s^2-42s+5=0$ and $t=\frac23s-\frac{10}{21}$. Put
    $$q_s=\tfrac{z^7}7-(s+1)\tfrac{z^6}6+(s+t)\tfrac{z^5}5-t\tfrac{z^4}4,\qquad r_s=\tfrac{-128(51s-65)}{1750329},\qquad\beta_s=1-q_s/r_s .$$
  - **Ramification.** $q_s'=z^3(z-1)(z^2-sz+t)$, and the passport is $[2^21^3,\,4\,2\,1,\,7]$ over $0,1,\infty$.
  - **Monodromy.** Reduction mod 17 at $s=3$ has an element of order 10, so the monodromy is $A_7$, not $L_2(7)$.
  - **Fields.** The map is defined over $\mathbb Q(\sqrt{21})$, and the marked actions over $\mathbb Q(\sqrt{21},\sqrt{-7})$.
  - **Klein's case.** The equal-critical-value condition factors as $(7s^2-7s+4)(21s^2-56s+40)(21s^2-42s+5)$. The first factor is Klein's $L_2(7)$ case.
  - **Integral form.** $P_\lambda=X^7-(\lambda+21)X^6+36(\lambda-6)X^5-324(\lambda-15)X^4$, with $\lambda^2-42\lambda+105=0$.
- **Resolvents.** Let $r_1,\dots,r_6$ be the other roots of $\beta(X)=\beta(z)$.
  - $C/L_2(5)$ (genus 1). Take $R_E=\prod_\Pi\big(U-\sum_{M\in\Pi}\sigma_M^3\big)$ over the six pentads of perfect matchings, with $\sigma_M=\sum_{ij\in M}r_ir_j$. The cube is needed, because lower powers give constants.
  - $C/3^2{:}4$ (genus 2). Take $R_S=\prod_{A\mid A^c}\big(V-\prod_Ar_i-\prod_{A^c}r_i\big)$ over the ten splittings into triples.
- **Next steps.**
  1. Weierstrass models.
  2. The images of $D_\delta$.
  3. Torsion bounds from good reduction, or a certified canonical height.

### 2.6 Groups and towers (DAY 2 §9; R5 §§3–4)

- **Gonality of $G$-curves.**
  - $\gamma(G)\ge\min(\lvert G\rvert,\lceil1+\sqrt{g_{\min}}\rceil)$ and $\gamma(G)\le\lfloor(g_{\min}+3)/2\rfloor$.
  - $\gamma(A_6)=5$, from the Valentiner smooth plane sextic.
- **Accessory degrees.**
  - $a(G)=\lvert G\rvert/\max_H(\lvert H\rvert/\mu(H))$.
  - $\mu(A_n),\,a(A_n)\le(n-2)!$, from the power-sum complete intersection; and $a(A_n)\le n!/84$ for $n\ge7$, by induction from $A_7$.
  - $a(L_2(q))\le q(q-1)/2$ for $q$ even or $q\equiv1\pmod4$, and $q(q+1)/2$ for $q\equiv3\pmod4$. Use the normaliser of an odd-order torus, which has $\mu=1$.
- **Towers.** Any tower of accessories ending in a curve compression has $\prod d_i\ge60$, and $\ge\mu(A_7)=90$ with connected full monodromy; the step count is unbounded. R5 Prop 4.1 (a base with a fixed point and $\mathrm{Alb}=0$) is generalised by Theorem 5.11.

## 3. Superseded statements [S]

| GPT statement | superseded by |
|---|---|
| $\operatorname{gon}(C)\ge17$ in genus 136, $\ge19$ globally; $\mathrm{ed}_{\mathbb C}(A_7;\le n)>1$ for $n=17,23,29$ | $a(A_7)=60$ (Ch. 5); $\gamma(A_7)\ge25$ (Ch. 4) |
| $\operatorname{gon}(C)\ge24$ and $\operatorname{gon}(C/\tau)\ge12$ | $\ge25$ and $\ge13$ (Ch. 3) |
| $\operatorname{gon}(C)\le56$, via the genus-12 Sylow-2 quotient | $\le42$ (Ch. 7) |
| $72\le\mu(A_7)\le120$, with $\mu\in\{72,84,90,96,108,120\}$ | $\mu=90$ (Ch. 5) |
| "seven signatures remain", "$Q_1$ class 14 remains", "the trial-space bounds are not certified" | Ch. 2–3 certify all of them |
| rungs G17, G18, A18 (DAY 1C §5) | settled by $a(A_7)=60$ and $\gamma\ge25$ |

The rungs **T1** ($\mathrm{RD}(A_7)\ge2$) and **T2** ($\mathrm{RD}(S_7)=3$, algebraic Hilbert 13) remain open. A bounded-accessory obstruction says nothing about unrestricted towers (DAY 1C §7; §5.5 here).

## 4. DAY 1B's repair list

| repair | status |
|---|---|
| 1. call $Q_2$ by its true name | done: $C_3\rtimes D_8$ (§2.3) |
| 2. no unqualified sparse-Cholesky claim | done: Higham's backward error with stored row counts, and the criterion in Lean (§2.4). DAY 1B's exported-factor residual check is an alternative |
| 3. a proved Bessel bound; rational $C_h^2$ | $C_h^2$ is enclosed in ball arithmetic. $j_{1,1}$ is a cited decimal, but $19/5$ would still give 0.340893 (§2.4 above) |
| 4. genuine characters | the ATLAS table, recomputed by Burnside–Dixon (`chartab.py`), with the exact Frobenius count (`cover.py`) |
| 5. the mean-zero subspace for $Q_0$ | min–max over two-dimensional subspaces needs no mean-zero subspace |
| 6. separate the main proof from superseded material | done (chapters, and git history) |
