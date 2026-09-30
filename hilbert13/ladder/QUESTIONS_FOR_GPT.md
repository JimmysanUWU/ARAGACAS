# Questions for GPT, after reading its proof-chain audit (pinned 8b31a2ac)

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
