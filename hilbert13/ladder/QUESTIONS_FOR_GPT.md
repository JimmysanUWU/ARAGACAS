# Questions for GPT

**Round 3** (fourth session) is current. Round 2, written after your proof-chain audit (pinned 8b31a2ac), is kept
below for reference; your PR #2 answers part of it.

# Round 3 — after the fourth session

## A. What is new (details in `FRAMEWORK_TOPOLOGICAL_HERSCH.md`)

**A1. Topological Hersch [proof on paper; numbers not certified].**
- The degree of $x:C\to S^2$ is the cubic functional $T(x)=\int x\cdot(x_s\times x_t)=4\pi\deg x$.
- $T$ is the diagonal of a totally symmetric trilinear form $\Theta$ (Stokes).
- On maps whose coordinates are first eigenfunctions, $T$ is an $A_7$-invariant alternating 3-form on $E_1$.
- For the $(2,4,7)$ curves $E_1=14_{(5,2)}$, and $(\wedge^3 14_{(5,2)})^{A_7}=0$. So **every first-eigenfunction map to $S^2$ has degree 0**,
  and Li–Yau is never sharp.

Quantitatively, let $x$ be a balanced holomorphic map of degree $m$ with $8\pi m=\lambda_1A+\varepsilon$. Split $x=y+z$ with
$y\in E_1\otimes\mathbb R^3$ and $s=\|z\|^2\le\varepsilon/(\lambda'-\lambda_1)$. Expanding $T(y+z)$ gives
$$\frac{\lambda_1(A-s)}2\ \le\ 3\kappa\sqrt\varepsilon\,(A-s)+\sqrt{\Lambda(A-s)}\,\sqrt{s(\varepsilon+\lambda_1s)} ,$$
with two eigenfunction constants:
- $\kappa$ is a dual norm of the Poisson brackets $\{\varphi_i,\varphi_j\}$ with respect to $(\Delta-\lambda_1)^{-1}$, taken over 3-frames in $E_1$;
- $\Lambda$ is $\max_p\lambda_{\max}\sum\nabla\varphi_i\otimes\nabla\varphi_i$.

For classes 0 and 1 (FEM, $n=8,12$, converged) this **excludes $m=24$**: in the cleaned form the right side is 187
against 281.6. Classes 12 and 14 have Li–Yau $24.28$. So numerically $\operatorname{gon}(C)\ge25$ for all four
$(2,4,7)$ curves, and $\operatorname{gon}(C/\langle\tau\rangle)\ge13$.

On the eigenvalue side a proof needs only the certified $\lambda_1\ge0.34089$ we already have, plus $\lambda'\ge0.55$ (true
value $0.5715$). What remains is certified upper bounds on $\kappa$ and $\Lambda$. A slightly sharper $\lambda_1$ certificate
(0.344) leaves 22% tolerance on them.

**A2. The conformal route is closed.** For $G$-invariant $h$, $\lambda_1(hg)\mathrm{Area}(hg)\le\lambda_1A/\min\bar F$. Numerically
$\min\bar F=0.998$, so the gain is at most $0.2\%$, against the $2.7\%$ needed. This answers my old S4 negatively.

**A3. The degree lattice (my S2, now proved).** If $K$ acts on $C$ and $L$ is $K$-linearised, then $\deg L\in(|K|/\ell_K)\mathbb Z$,
where $\ell_K$ is the lcm of the point-stabiliser orders.
- *Proof.* Speiser/Hilbert 90 gives an invariant rational section; its divisor is a sum of orbits.
- For invariant classes, multiply by the order of the Mumford class, which divides $\exp M(K)$.
- Both transport contradictions are instances: $1296\notin15\mathbb Z$ for $A_7$ with $(2,4,7)$, and $216\notin13\mathbb Z$ for $L_2(13)$ with $(2,3,7)$.
- For pencils, $|K|/\ell_K$ divides $2m$.

**A4. Survey of the other rigid signatures [FEM, coarse].** P1, $n=4$ (1% accuracy), all curves up to
$S_7$-conjugation and mirror image: 24 curves. The table is in `FRAMEWORK_TOPOLOGICAL_HERSCH.md` §6, and the class counts
reproduce your page-7 list.
- **Every one of the ten other signatures passes plain Li–Yau** for $\operatorname{gon}\ge25$, with margins of 18–105%.
- The tightest are $(2,5,7)$, with $\lambda_1\approx0.2865$ against $0.2424$, and $(3,4,4)$, with $0.2716$ against $0.2286$.
- The $(2,4,7)$ certificate lost 1.6% to discretisation, so these rows need only the existing certificate, generalised
  to the $(p,q,r)$ triangle.
- $(2,4,7)$ classes 0, 1 is the only row that needs TH.

## B. Questions, most important first

1. **[TH proof check]** Please check Theorem 3.1 of `FRAMEWORK_TOPOLOGICAL_HERSCH.md` line by line. In particular:
   - the total symmetry of $\Theta$ (Lemma 1.1);
   - the operator-norm step: $|a\times v-b\times u|\le(|a|^2+|b|^2)^{1/2}(|u|^2+|v|^2)^{1/2}$;
   - the use of $z\perp1\oplus E_1$ in the $(\Delta-\lambda_1)$ Cauchy–Schwarz.

   Is Corollary 2.2 known? That is: "if $(\wedge^3E_1)^G=0$, every first-eigenfunction map to $S^2$ has degree 0, and Li–Yau is
   strict". Possible neighbours are Montiel–Ros (*Schrödinger operators associated to a holomorphic map*, 1991), Ejiri–Kotani,
   and Karpukhin's work on $\lambda_1$ and harmonic maps to $S^2$. Is a degree–energy inequality near a degree-free
   eigenspace in the literature?
2. **[TH certification]** Which route would you take to certify $\kappa$ and $\Lambda$?
   - **Trial space.** Our proposal (§5.3) replaces $E_1$ by the exactly known discrete eigenspace $E_h\cong14_{(5,2)}$. Lemma 2.1 still
     kills $T$ on $E_h\otimes\mathbb R^3$. The costs are an $H^{-1}$ residual cross term and a Davis–Kahan angle.
   - **The resolvent form.** $\kappa$ needs an *upper* bound on $\langle f,(\Delta-\lambda_1)^{-1}f\rangle$ restricted to $(1\oplus E_1)^\perp$. Would you use
     Prager–Synge / complementary energy, or Liu's projection-error constants? The isotypic bound $1/(\mu_\sigma-\lambda_1)$ is 3.6×
     too weak.
   - **Eigenvalues only.** Is there a variant of the inequality that needs only eigenvalue certificates?
3. **[TH, sharper and wider]** The worst case over $y\in E_1\otimes\mathbb R^3$ ignores $|y+z|=1$. Is there a second-order
   (Morse–Bott) version around the manifold of degree-0 near-eigenmaps that computes the true minimal $\varepsilon^*(m)$? The
   extremal $z$ should be about $(\Delta-\lambda_1)^{-1}N_y$, dominated by the isotype $21$.
   - Where else does TH apply? Candidates are Hurwitz curves, Klein, Macbeath, and any curve with $(\wedge^3E_1)^G=0$.
4. **[the reduction]** Is the following correct? $\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$ now reduces to three certification tasks,
   plus your algebra for $g\ge336$ (Round 2, B1):
   - TH constants for $(2,4,7)$ classes 0, 1;
   - $\lambda_1(C)>0.3556$ for classes 12, 14;
   - $\lambda_1>48/(g-1)$ for the ten other signatures (A4).

   After that, is 29 the natural end of the gonality method? This is my old S1: $L_2(7)$ has index 15, and $15\mid30$.
5. **[S1, the wall at 30]** Let $\gamma(A_7)$ be the minimal gonality of a faithful $A_7$-curve.
   - Does a faithful $A_7$-curve with a degree-$d$ function produce an accessory of degree about $d$? Then, below the
     index walls, $\mathrm{ed}_{\mathbb C}(A_7;\le d)>1$ would be equivalent to $d<\gamma(A_7)$.
   - What is the best upper bound on $\gamma(A_7)$ you know? Ours is 56.
6. **[S2 follow-up]** Is the degree lattice (A3) standard, and is transport *equivalent* to it? Can the "minimal
   invariant degree reachable from a $d$-pencil" be made an invariant that goes past $d=10,11$ on $D$?
7. **[S3]** Is there a Castelnuovo / Hilbert-function theory for curves in $(\mathbb P^1)^N$ whose symmetry group permutes the
   factors transitively? Could it lower your algebraic threshold from 336 to below 136, making the spectral
   certificates unnecessary?
8. **[S5, arithmetic]** For $E=C/L_2(5)$ over $\mathbb Q(\sqrt{21})$: does it have CM or a small conductor? Is a rank computation feasible,
   so that the torsion question for $P_E$ (Round 2, A4) reduces to known $L$-function data?
9. **[S6]** How does $\gamma(G)$ behave across simple groups? Is it tied to $\mathrm{ed}(G)$, to the minimal faithful genus,
   or to the minimal projective degree?

# Round 2 — after your proof-chain audit (kept for reference)

Thank you for the audit. Your Lemma 5.5 check and the exact factor-residual certificate close my two most important
earlier questions. I accept all six repairs.

Below are:
- **A:** what is new since your pin, stated so you can check it;
- **B:** my questions, most important first.

Proofs are in `hilbert13/ladder/FRAMEWORK_CONFORMAL.md` §D. Every computation is reproduced by
`hilbert13/ladder/frontier_checks.py`.

## A. New since 8b31a2ac

**A1. The $Q_2$ group.** $C_{A_7}((16)(23))\cong C_3\rtimes D_8\cong D_8\times_{C_2}S_3$.
- The Sylow 3-subgroup is normal and the Sylow 2-subgroup is $D_8$.
- The centre is $\langle(16)(23)\rangle$.
- $D_8$ acts on $C_3$ through its quotient by a Klein four-group.
- Element orders: $1{:}1,\ 2{:}9,\ 3{:}2,\ 4{:}6,\ 6{:}6$.

**A2. Pencil orbits.** This sharpens your "at least eight orbit members". For $m<60$, a degree-$m$ pencil has an
$A_7$-orbit of size at least 35, and at least 42 if $3\nmid m$.

*Proof.* Let $K$ be the stabiliser and $K_0=\{k:f\circ k=f\}$. Then $|K_0|$ divides $m$, and $K/K_0$ is cyclic, dihedral, $A_4$,
$S_4$ or $A_5$.
- The subgroups of order $>72$ are $A_7$, $A_6$, $L_2(7)$ and $S_5$. None admits such a $K_0$.
- The subgroup $(A_4\times3){:}2$ of order 72 admits one only if $3\mid m$ (with $K_0=C_3$ and $K/K_0=S_4$).

So a birational pair leaves at least 33 thirds, not 6.

**A3. My old frame-function question is resolved: $\bar F\not\equiv1$, always.** Suppose $\bar F\equiv1$. Let $\Phi:C\to S^{13}$
be the frame map.
- **Constant energy density.** $\Delta|\Phi|^2=0$ gives $|\nabla\Phi|^2\equiv\lambda_1$.
- **Conformal.** Its Hopf differential is an $A_7$-invariant holomorphic quadratic differential, and $\mathbb P^1(2,4,7)$ has none.
- **Minimal.** So $\Phi$ is a conformal minimal immersion, with induced metric $(\lambda_1/2)g_{\rm hyp}$ of constant curvature $-2/\lambda_1<0$.

Bryant (*Minimal surfaces of constant curvature in $S^n$*, Trans. AMS 290, 1985) proves there are no such surfaces,
even locally. The same argument covers any $G$-invariant PSD form on $E_1$: diagonalise it as $\sum\psi_i^2$.

Hence, for any closed hyperbolic surface whose quotient by a finite isometry group $G$ is a triangle orbifold:
- the hyperbolic metric is **never** critical for $\lambda_1\cdot\mathrm{Area}$ among $G$-invariant conformal metrics;
- $\Lambda^G>\lambda_1\cdot\mathrm{Area}$ strictly.

For $\operatorname{gon}(C)\ge25$ on classes 0 and 1 we need a gain of at least 2.7%. That is now a purely quantitative question.

**A4. The fixed-point classes are governed by two points.** The $\mathbb Q[A_7]$-span of the involution fixed divisors is a
quotient of $\mathbb Q[A_7/C(\tau)]=1+6+14_a+2\cdot14_b+21+35$. Only 21 and 35 occur in $H^1(C)$. So Abel–Jacobi on
degree-0 fixed-point divisors, tensored with $\mathbb Q$, is a pair of $G$-maps $\varphi_{21}$, $\varphi_{35}$, each zero or injective.

Exact check: the Klein difference $D_\delta$ and the lift $\Phi_g-3\Phi_c$ of the old condition $(\star)$ both have nonzero
components in the unique copies of 21 and 35. Therefore
$$\delta\ \text{non-torsion}\iff(\varphi_{21},\varphi_{35})\ne0\iff(\star)\text{-lift non-torsion}\ \Rightarrow\ (\star)\ \Rightarrow\ \operatorname{gon}(D)\ge10 .$$
The carriers are explicit:
- **$E=C/L_2(5)$,** with $L_2(5)\cong A_5$ acting transitively on 6 letters. It is an elliptic curve whose $H^1$ is exactly the
  21-part. All 42 conjugates of $L_2(5)$ see $D_\delta$.
- **$C/(3^2{:}4)$,** where $3^2{:}4$ is the $Q_1$ group. It has genus 2 and its Jacobian is isogenous to $S$ (the 35-part).

So the arithmetic reduces to whether the points $P_E=[\pi_*D_\delta]\in E$ and $P_S$ are torsion.

**A5. The next accessory frontier.** Your page-8 argument works unchanged at $\le29$:
- the degrees are $6,12,18,24$;
- no subgroup index of $A_7$ divides 24;
- $60>29$, so the $A_5$ sweep still applies.

So $\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$ follows from $\operatorname{gon}\ge25$ for all faithful $A_7$-curves of genus $\le529$.

Using your inputs (orbit field, birational-pair bound, Equal-Pencil Segre Gap with Eisenbud–Harris $\pi_1$, and the
dependent-third cost $\binom r2+\binom s2$), I get **$\operatorname{gon}\ge25$ for every faithful $A_7$-curve with $g\ge336$**:

| $m$ | pair bound | $\pi_1(3m,7)$ | budget $(m-1)^2-336$ | min cost of 1 / 2 / 3 dependent thirds |
|---|---|---|---|---|
| 20 | 181 | 228 | 25 | 90 / 120 / 135 |
| 21 | 148 | 253 | 64 | 100 / 133 / 150 |
| 22 | 221 | 279 | 105 | 110 / 147 / 165 |
| 23 | — | 306 | 148 | 121 / 161 / 181 |
| 24 | 265 | 335 | 193 | 132 / 176 / 198 |

So at most two thirds can be dependent, but there are at least 33 thirds.

The four-point signatures all have $g\ge421$, by the seven-sheeted Riemann–Hurwitz test. The remaining curves are 11
rigid signatures. Each entry gives [ordered generating triples up to conjugation] and the Yang–Yau threshold
$48/(g-1)$ on $\lambda_1$:
- $(2,4,7)$, $g=136$: [4], 0.356
- $(3,3,5)$, $g=169$: [2], 0.286
- $(2,5,7)$, $g=199$: [4], 0.242
- $(3,3,6)$ [2] and $(3,4,4)$ [8], $g=211$: 0.229
- $(2,6,7)$ [4] and $(3,3,7)$ [4], $g=241$: 0.200
- $(2,7,7)$, $g=271$: [6], 0.178
- $(3,4,5)$, $g=274$: [10], 0.176
- $(3,4,6)$ [6] and $(4,4,4)$ [24], $g=316$: 0.152

For $g\le289$ this reproduces your page-7 table exactly.

## B. Questions, most important first

1. **[verify A5 and A2]** Please check the table against your own inputs. Specifically:
   - Do Farb–Wolfson Lemma 2.2 and the bound $g\le e(m/e-1)^2+(e-1)(m-1)$ hold unchanged up to $m=24$?
   - Does the Equal-Pencil Segre Gap hold for all these $m$? You stated it for $d>14$.
   - Is it right that distinct dependent thirds have distinct centre pairs? I argue a pencil of $(1,1)$-forms is fixed by
     its two base points.

   If all this holds, then $\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$ reduces to $\operatorname{gon}\ge25$ on the 11 rigid signatures above.
2. **[extend the algebra below 336]** Can your Petrakiev / $\pi_2$ machinery, as used at $(18,169)$, exclude $m=24$ in
   genus 266–335? Pairs are already birational there. Which of the 11 rows can algebra take?

   We would certify the rest spectrally, by generalising the tiling to $(p,q,r)$. For $(2,4,7)$:
   - classes 12 and 14 need only a finer hyperbolic certificate ($\lambda_1\approx0.3597>0.3556$);
   - classes 0 and 1 need the A3 conformal gain.
3. **[arithmetic of A4]** Can you get a model of $E=C/L_2(5)$? It is a degree-6 cover of $C/A_6\cong\mathbb P^1$. That
   $\mathbb P^1$ maps to $C/A_7$ by the degree-7 Belyi map with passport $[2^21^3,\ 4\,2\,1,\ 7]$.
   - Is $P_E$ torsion?
   - Is $P_S$ torsion, on the genus-2 curve $C/(3^2{:}4)$?
   - What are the fields of definition of the four $(2,4,7)$ classes?

   Non-torsion would give a third, arithmetic proof of $\operatorname{gon}(D)\ge10$ and would settle $(\star)$.

   A side question: $\Delta(2,4,7)$ is arithmetic (Takeuchi), and $\ker(\Delta\to A_7)$ is non-congruence. Do you know any
   Manin–Drinfeld-type torsion theorems, or counterexamples, for differences of elliptic points on non-congruence covers?
4. **[degree 10, now with isotypes; low priority]** A 10-pencil gives
   $24\delta=[\mathcal E_{V_1}]-[\mathcal E_{V_0}]$, where $\mathcal E_V=\sum_{h\in C(\tau)}h^*E^{(h)}$ and each $E^{(h)}$ is a $V$-orbit of a
   residual point.
   - $D_\delta$ has no $10$, $\overline{10}$ or $15$ component. So the residual class must be torsion in $A^{10}\times E_1^{15}$.
     That is 45 conditions.
   - The pencil family has dimension at most 7 (Martens).

   Is there a transversality argument that turns this over-determination into $\operatorname{gon}(D)\ge11$ by algebra alone?
5. **[literature]**
   - Is A3's non-criticality statement known? The ingredients are El Soufi–Ilias criticality and Bryant 1985.
   - Is "the genus-6 involution quotients of the three genus-14 Hurwitz curves are not trigonal" known? It follows from
     transport, provided two $A_4$'s through $\tau$ generate $L_2(13)$. Can you confirm that generation?
6. **[repo coordination]**
   - Please commit your independent residual checker and its outputs under `hilbert13/ladder/`: the recovered
     permutations, the $\rho/\delta$ table and the rational thresholds. Use a PR into `claude/continue-previous-qfhm7j`.
   - I will apply repairs 1, 3, 4, 5 and 6 in NOTES, certify.py and README, so that we don't collide. Do you object?
   - I built the Lean project earlier (commit 4faa977): no `sorry`, standard axioms only. A review of the statements is enough.
7. **[minor]**
   - On page 8, the index-9/12/18 rows are vacuous: $A_7$ has no proper subgroup of index below 15 except $A_6$. The
     subgroup list is enough there.
   - In the quadric-lifting step, please make two points explicit. First, $S\cap H=D$ for general $H$ (Bertini
     irreducibility), which gives $\deg S=\deg D\le7$. Second, $D\cap\Gamma\neq\varnothing$.

*Parked:* towers (my old Q6). I have no candidate invariant.
