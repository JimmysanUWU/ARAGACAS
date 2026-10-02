# Questions for GPT

**Current: Round 7.** Rounds 6, 5, 4, 3 and 2 are kept below for reference.

# Round 7 — your Theorem 5.3 is certified

**New: $\operatorname{gon}\ge25$ on $(2,4,7)$ classes 0, 1, certified** (`HH_CERTIFICATE.md`, `hh_certify.py`, outputs
`hh_certify_output_cls*.txt`). We re-derived your trial-space form (5.3) and certified it with our own trial space. The left side
is $\ge0.3307$ against a right side $\le0.2713$ (class 0); for class 1, $0.3296$ against $0.2795$.

**The method has three new ingredients.**
1. **Vector-valued Hejhal.** The $14_{(5,2)}$-isotype is expanded exactly, as $\sum_mR_m(|z|)(A_m\cos m\theta+B_m\sin m\theta)$
   about the order-7 point, with $R_m$ the regular radial eigenfunction. This gives $\lambda_1$ to 12 digits: $0.346267085404$
   for classes 0, 1 and $0.359671354895$ for classes 12, 14.
2. **Exact equivariance.** A $\Delta$-invariant $C^2$ partition of unity glues $\rho$-translates of the expansion. The
   trial space is then an exact $14_{(5,2)}$. The residual is $\sum_c[\Delta,\chi_c](F_c-F_0)$, and the mismatches are exact
   eigenfunctions, bounded by circle sampling with rigorous aliasing bounds. Result: $\rho\le1.9\cdot10^{-4}$ (class 0).
3. **Jacobian flux** (this replaces your complementary-energy solve). Since $\{u,v\}=\operatorname{div}J$ with
   $J=\tfrac12(u\nabla^\perp v-v\nabla^\perp u)$, we have $\langle f,(\Delta-\lambda_h)^{-1}f\rangle\le\frac b{b-\lambda_h}\|J_f\|^2$. So $B_h\le4.45\cdot10^{-4}$
   from four isotype integrals, with no PDE solve.

**Consequences.**
- $\gamma(A_7)\ge25$, modulo your algebra for $g\ge336$.
- $\operatorname{gon}(C/\langle\tau\rangle)\ge13$ on all $(2,4,7)$ classes.
- A two-sided $\lambda_1\in[0.34610,0.34633]$.

**Questions.**
1. **Audit `HH_CERTIFICATE.md`**, in particular:
   - the Jacobian-flux lemma;
   - the commutator form of the residual, using $h(d(z,0))\equiv1$ on the sector;
   - the uniform radial bounds $g_{\rm lo}\le R_k(u)/u^k\le G_{\rm hi}$ for all $k$, from $|(b_k)_n|\le(k+1)_n$;
   - the mirror identification. The `orbifold.py` tiling curve of a triple is the mirror of $\mathbb H/\ker\varphi$, $\varphi(X,Y,C)=(a,b,c)$, since $\sigma_{BC}\sigma_{AB}=Y^{-1}$.
2. **The large-genus algebra.** To make $\gamma(A_7)\ge25$ unconditional we need a self-contained proof that every faithful
   $A_7$-curve with $g\ge336$ has $\operatorname{gon}\ge25$. Our Round-3 table derives it from four inputs: the birational-pair bound,
   Farb–Wolfson's $g\le e(m/e-1)^2+(e-1)(m-1)$, the Equal-Pencil Segre Gap with Eisenbud–Harris $\pi_1$, and the dependent-third cost.
   Can you write it as one complete proof with exact references? A simpler route for $g>529$ (Castelnuovo–Severi plus subgroup
   structure) would also do.
3. **$\operatorname{gon}\ge26$?** Numerically, (HH) excludes $m=25$ for classes 0, 1 if $\kappa^2=0.614\,B^*/3$, the BFGS maximum; the crude
   $B^*/3$ does not suffice. Is there a certifiable relaxation of $\sup_{\|W\|=1}\sum_{\rm cyc}Q(w_j\wedge w_k)$ on $S^{41}$, using the isotype split
   $Q=\sum_\sigma b_\sigma P_\sigma$? For example an SOS bound, or a bound on the $P_{10+\overline{10}}$ and $P_{15}$ mass of decomposable 2-vectors.
   (Classes 12, 14 stay at 25: their $\lambda'\approx0.386$ is too close.)
4. **Round 6, Q1 (the audit of `MU90.md`) and Q3 (towers)** are still open, if you have not answered them.


# Round 6 — after reviewing your Round-5 answers

**We verified everything in your Round 5** (`ROUND5_REVIEW.md`). Thank you for the bug report; it is fixed. Harmonic Hersch
replaces our inequality.

**New:** $\mu(A_7)\le90$. Holomorphic Lefschetz gives $\chi_{A_7}(B+T)=-6+10-14_a-14_b-21$ on all four $(2,4,7)$ classes. So
$h^0(B+T)\ge10$, and $\mu(A_7)\in\{72,84,90\}$ (`equivariant_rr.py`). On every rigid signature with lattice degree 72 or 84 the same
$\chi$ has no positive part.

1. *(Answered by us since: $\mu(A_7)=90$, `MU90.md` — please audit it.)* **$\mu(A_7)$.** How would you decide degree 72 on the $(2,5,7)$ curves (the unique class $3D_5-4D_7$) and degree 84 on $(3,4,5)$?
   Possible routes: a lower bound on $h^1$, Clifford/Castelnuovo with the 10-dimensional or 6-dimensional section
   representation, or a Brill–Noether-type vanishing for invariant classes. Is $h^0(B)>0$ for the other degree-90 class?
2. **Certification of Thm 5.3.** Any refinement before we implement it? In particular: the S5-quotient trial space and its
   element order, the flux construction across sign-twisted gluings, and how to enclose the hyperbolic metric in
   $\nabla u-p$.
3. **Towers.** Given Prop. 4.2, is there a weaker invariant than a fixed point (for example, a fixed point after a bounded
   further accessory, or the vanishing of an obstruction class in $H^1(G,\mathrm{Pic})$) that does persist?

# Round 5 — after DAY 2

**We verified your accessory theorem.** `ACCESSORY_60.md` and `verify_accessory60.py` on `claude/continue-previous-qfhm7j`.
- **Re-derived on paper:** Thm 1.2, the lower bound of Thm 1.5, Prop. 2.2, Thm 3.1, and every $\mu(H)$ bound.
- **Checked exactly:**
  - the minimum genus is 136; an exhaustive search shows that no triangle signature below 136 is generated by $A_7$;
  - the complete list for genus 136–325, the seven-sheeted exclusions, and the lattice numbers $N$;
  - $\pi(54,5)=325$, $\pi(36,5)=136<199$, $\pi(42,5)=190<274$.
- **Your Thm 2.1** is the same degree lattice as our Lemma 7.1, proved by the same Hilbert 90 argument.

So $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$ holds with no gonality input, and it holds in either convention. We have re-scoped the repo
accordingly.

**Updates for your §7 and §9.** Your reference was built from PR #2 at `284e1e2`. Since then (`VERIFICATION.md`):
- **Certified (Arb + verified factorisations):** $\operatorname{gon}\ge25$ for all 70 $A_7$-classes of the ten rigid signatures other than
  $(2,4,7)$, and for $(2,4,7)$ classes 12, 14 ($\lambda_1\ge0.355696$). Your p. 15 thresholds are therefore met, except for
  $(2,4,7)$ classes 0, 1.
- **$(2,4,7)$ classes 0, 1:** a new inequality, *topological Hersch* (`FRAMEWORK_TOPOLOGICAL_HERSCH.md`), excludes degree 24. Its
  eigenvalue inputs are certified: $E_1$ is one copy of $14_{(5,2)}$, $\lambda_1\in[0.34089,0.36318]$, and $\lambda'\ge0.55998$ via a verified
  $LDL^T$ inertia count. Two eigenfunction constants, $\kappa$ and $\Lambda$, are not yet certified (3.5% tolerance).
- **The conformal route has a quantitative ceiling:** $\Lambda^G\le\lambda_1A/\min\bar F$, which is at most a 0.2% gain numerically
  (`FRAMEWORK_CONFORMAL.md` §E). So Vinokurov's optimiser cannot reach 25 for classes 0, 1.
- **The $n=96$ $Q_1$ certificate** is in `certificate.txt`, with its runner `run_certificate.py`.

**Questions.**
1. **Upper bound.** Please state exactly which result of Farb–Wolfson [1] gives the independence of the faithful linear
   model in Thm 1.5, and check that it covers torsors that become disconnected after the accessory.
2. **$\mu(A_7)$ itself.** It is at least 60. What is the smallest *linearised* moving degree on a faithful $A_7$-curve? Our
   current upper bound is 270, from the canonical bundle of a $(2,4,7)$ curve. Is anything below 270 known or plausible?
   This is the full-monodromy threshold.
3. **Beyond $A_7$.** Does $a(G)=\min_H[G:H]\mu(H)$ together with linearised Castelnuovo compute $a(A_n)$ or $a(L_2(q))$ in
   general? Is there a pattern, such as $a(G)$ attained by a small-genus subgroup curve? For $A_5$ it reproduces
   Kronecker–Klein ($a=2$).
4. **Towers.** Is there any analogue of $\mu$ for iterated accessories, i.e. a fixed point and $\mathrm{Pic}(V\times C)=\mathrm{Pic}(C)$ on
   intermediate varieties? This is the obstacle you name for $\mathrm{RD}(A_7)>1$.
5. **Gonality (now a separate track).** Same as Round 4, Q2: how would you certify $\kappa$ and $\Lambda$? With them, $\operatorname{gon}\ge25$
   holds for every faithful $A_7$-curve, and so $\gamma(A_7)\ge25$.

# Round 4 — distilled

**Where things stand** (`VERIFICATION.md` §1):
- **Certified:** $\operatorname{gon}\ge25$ for all 70 $A_7$-classes of the ten rigid signatures other than $(2,4,7)$, and for $(2,4,7)$
  classes 12 and 14.
- **Classes 0, 1:** $\operatorname{gon}\ge25$ follows from a new inequality, *topological Hersch*
  (`FRAMEWORK_TOPOLOGICAL_HERSCH.md` §3, proof audited). All its eigenvalue inputs are certified:
  - $E_1=14_{(5,2)}$, a single copy;
  - $\lambda_1\in[0.34089,0.36318]$;
  - next eigenvalue $\ge0.55998$.

  Two eigenfunction constants, $\kappa$ and $\Lambda$, are only computed, not certified.

1. **Check the theorem.** The degree of $x:C\to S^2$ is the cubic $T(x)=\int x\cdot(x_s\times x_t)$. Restricted to $E_1\otimes\mathbb R^3$, $T$ is
   an invariant alternating 3-form, so $(\wedge^3E_1)^G=0$ forces it to vanish there. Expanding $T(y+z)$ then gives
   $$\tfrac12\lambda_1(A-s)\le3\kappa\sqrt\varepsilon\,(A-s)+\sqrt{\Lambda(A-s)\,s(\varepsilon+\lambda_1s)},\qquad8\pi m=\lambda_1A+\varepsilon,\quad s\le\varepsilon/(\lambda'-\lambda_1).$$
   - Is the proof (`VERIFICATION.md` §2) correct?
   - Is this inequality, or its corollary "first-eigenfunction maps to $S^2$ have degree 0, so Li–Yau is strict", known?
2. **The last gap: certify $\kappa$ and $\Lambda$ within 3.5%.** The definitions are
   $$\kappa^2=\sup_{W\in E_1^3}\frac{\sum_{\rm cyc}\langle\beta(w_a\wedge w_b),(\Delta-\lambda_1)^{-1}\beta(w_a\wedge w_b)\rangle}{|W|^4},\qquad
   \beta(\omega)=\sum_{i<j}\omega_{ij}\{\varphi_i,\varphi_j\},\qquad \Lambda=\max_p\lambda_{\max}\sum_i\nabla\varphi_i\otimes\nabla\varphi_i .$$
   Computed: $\kappa^2=0.6143\,b^*_{\max}/3$ and $\gamma=\sqrt{\Lambda A}=1.70$. Which route would you take?
   - (a) A trial-space version using the exact P1 eigenspace, with a Davis–Kahan angle, an $H^{-1}$ residual, and an *upper*
     bound for the resolvent form (Prager–Synge? Liu's projection constants?).
   - (b) Eigenfunction enclosures (Plum, Nakao).
   - (c) A variant of the inequality that needs only eigenvalue certificates.

   Separately, the global maximum of the quartic $F(W)/|W|^4$ on $S^{41}$ needs a certificate (SOS/Lasserre?). The cruder
   decomposable bound fails.
3. **The other input.** Please re-verify two things:
   - $\operatorname{gon}\ge25$ for every faithful $A_7$-curve of genus 336–529 (Round 2, A5 table: pair bound, Eisenbud–Harris $\pi_1$,
     dependent-third costs, at least 35 orbit pencils);
   - the transfer "$\operatorname{gon}\ge25$ for $g\le529$ $\Rightarrow\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$".

   With question 2, this would complete the chain.
4. **The wall at 30.** Is 29 a true limit of the gonality method ($[A_7:L_2(7)]=15\mid30$)? Does a degree-$d$ function on a faithful
   $A_7$-curve give an accessory of degree about $d$? What is the best known upper bound for the minimal gonality of a faithful
   $A_7$-curve? Ours is 56.
5. **Degree lattice.** For $K$-linearised $L$, $\deg L\in(|K|/\ell_K)\mathbb Z$, by Speiser. Is this standard, and is transport exactly this
   obstruction?
6. **Arithmetic (optional).** Does $E=C/L_2(5)$ over $\mathbb Q(\sqrt{21})$ have CM or a small conductor? Is $P_E$ torsion? Non-torsion would
   give an arithmetic proof of $(\star)$.

# Round 3 — after the fourth session (superseded by Round 4; numbers in A4 are corrected in `VERIFICATION.md` §6)

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
