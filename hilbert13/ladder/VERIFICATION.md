# Verification report (fifth session)

This report re-checks every claim of `FRAMEWORK_TOPOLOGICAL_HERSCH.md` and of the fourth-session additions to
`FRAMEWORK_CONFORMAL.md`.

**Status tags.**
- **[P✓]** paper proof, re-derived line by line in this session.
- **[X]** exact finite computation (integers, rationals, $\mathbb Q(\sqrt{-7})$), by `verify_exact.py`.
- **[C]** certified. The proof has the same form as `certificate.txt`: ball-arithmetic coefficient bounds (python-flint),
  the Crouzeix–Raviart lower-bound theorem (Liu; Carstensen–Gedicke–Rim), and verified floating-point factorisations
  with a-priori rounding bounds.
- **[N]** numerical only.
- **[cited]** taken from earlier sessions or from GPT, and not re-checked here.

## 1. Bottom line

| statement | status |
|---|---|
| $\operatorname{gon}\ge25$ for every faithful $A_7$-curve of the ten rigid signatures other than $(2,4,7)$ (24 curves up to isometry, all 70 $A_7$-classes) | **[C]** (§6) |
| $\operatorname{gon}\ge25$ for $(2,4,7)$ classes 12, 14 | **[C]** (§7) |
| $(2,4,7)$ classes 0, 1: $\operatorname{gon}\ge24$ | [C] (`certificate.txt`) |
| $(2,4,7)$ classes 0, 1: $\operatorname{gon}\ge25$ | **[C]** (ninth session: harmonic Hersch with a certified trial space, `HH_CERTIFICATE.md`; supersedes the [N] constants of §5) |
| Topological Hersch inequality (Theorem 3.1) and Corollary 2.2 | [P✓] |
| $(\wedge^3V)^{A_7}=0$ exactly for $V=6,14_{(5,2)},14_{(4,3)}$; $\wedge^2 14_{(5,2)}$ is multiplicity-free | [X] |
| Conformal ceiling (Proposition E1) | [P✓]; its numbers are [N] |
| Degree lattice (Lemma 7.1) | [P✓]; the transport degrees 1296 and 216 are [cited] |

**Update (ninth session).** The last gap is closed. `HH_CERTIFICATE.md` certifies GPT's Round-5 harmonic Hersch inequality, in
its trial-space form, for classes 0 and 1. The trial space is built from vector-valued Hejhal expansions. So $\operatorname{gon}\ge25$
now holds **[C]** for every faithful $A_7$-curve of genus $\le335$; for larger genus it rests on GPT's algebra [cited].

**Consequence (updated, sixth session).** The ed question no longer depends on any of this. GPT's DAY-2 theorem
$a(A_7)=60$ (`ACCESSORY_60.md`, verified) gives $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$ with no gonality input. The results here now
serve the separate gonality question. The only gap left for $\operatorname{gon}(C)\ge25$ on every faithful $A_7$-curve is the
pair $\kappa$, $\Lambda$ for $(2,4,7)$ classes 0, 1. GPT's algebra covers $g\ge266$.

## 2. Topological Hersch: proof audit [P✓]

Every step of `FRAMEWORK_TOPOLOGICAL_HERSCH.md` §§1–3 was re-derived. The points that needed care:

1. **Lemma 1.1 (total symmetry).** Write $\Theta=\tfrac12\sum\epsilon_{ijk}\int a_i\,db_j\wedge dc_k$.
   - $(b,c)$: swapping the wedge gives $-1$, and relabelling $\epsilon_{ijk}\to\epsilon_{ikj}$ gives another $-1$.
   - $(a,b)$: Stokes on the exact form $d(a_ib_j\,dc_k)$ gives $-1$, and $\epsilon_{ijk}\to\epsilon_{jik}$ gives another $-1$.
   - All fields are smooth: $x$ is holomorphic, and $y$ consists of eigenfunctions.
2. **Lemma 2.1.** $T(y)=\tfrac12\sum\theta_{ijk}\det(W_i,W_j,W_k)$ (the factor $\tfrac12$ was missing; corrected).
   - $\theta$ is alternating by Stokes.
   - $\theta$ is $G$-invariant because each $g$ is an orientation-preserving diffeomorphism.
   - Hence $\theta\in(\wedge^3E^*)^G=0$.
   - $\int\{\varphi,\psi\}\,dA=\int d\varphi\wedge d\psi=0$, so the brackets are orthogonal to $1\oplus E_1$.
3. **Theorem 3.1.**
   - *Energy.* $\langle\nabla y,\nabla z\rangle=\lambda_1\langle y,z\rangle=0$, $\|y\|^2+\|z\|^2=\int|x|^2=A$, and $E(z)\ge\lambda'\|z\|^2$ on $(1\oplus E_1)^\perp$.
   - *Expansion.* $3\Theta(y,z,z)+T(z)=2\Theta(z,y,z)+\Theta(x,z,z)$ by symmetry.
   - *The $\Theta(x,z,z)$ term.* $\Theta(x,z,z)=\int x\cdot(z_s\times z_t)$, and $|x|=1$ gives $\le E(z)/2$, by conformal invariance of the energy.
   - *The $2\Theta(z,y,z)$ term.* $2\Theta(z,y,z)=\int z\cdot(y_s\times z_t-y_t\times z_s)$. The map $(u,v)\mapsto a\times v-b\times u$ has Gram matrix
     $(|a|^2+|b|^2)I-aa^T-bb^T$. Its top eigenvalue is $|a|^2+|b|^2$, because $aa^T+bb^T$ has rank $\le2$ in $\mathbb R^3$. The
     conformal factors cancel in $\int|z||\nabla y||\nabla z|\,ds\,dt$.
   - *Gradient bound.* $|\nabla y|^2(p)=\sum_k w_k^TJ^TJw_k\le\lambda_{\max}(JJ^T)\|y\|^2$, with $JJ^T=S(p)$.
   - *Dual norm.* The form $q(z)=E(z)-\lambda_1\|z\|^2$ is positive on $(1\oplus E_1)^\perp$. Componentwise Cauchy–Schwarz followed by
     Cauchy–Schwarz over $k$ gives $\sqrt{\sum_kq(z_k)}=\sqrt\varepsilon$.
   - *Components.* The $k$-th component of $N_y$ is $\beta(w_a\wedge w_b)$ with $(a,b,k)$ cyclic, since
     $(W_i\times W_j)_1=(w_2\wedge w_3)_{ij}$.
   - *Combination.* $4\pi m-E(z)/2=\lambda_1(A-s)/2$.
   - *Monotonicity in $s$.* Holds when $6\kappa\sqrt\varepsilon<\lambda_1$ and $s\le A/2$; both are checked numerically for every $\lambda_1$ used.
4. **Remark 3.2.** $Q^*$ is $G$-invariant, because the bracket is equivariant and the resolvent commutes with $G$.
   - Over $\mathbb C$, $\wedge^2 14_{(5,2)}=10+\overline{10}+15+21+35$ [X], multiplicity-free.
   - Over $\mathbb R$, $15$, $21$, $35$ are of real type: a real-valued irreducible of quaternionic type would occur with even
     multiplicity in a real module.
   - $10\oplus\overline{10}$ is one real irreducible of complex type, with $\mathrm{End}_G=\mathbb C=\langle I,J\rangle$ and $J$ skew. So every invariant
     symmetric form is scalar on each block.
   - The SVD reduction uses $F(WO)=F(W)$: $\Omega\mapsto\Omega\,\mathrm{cof}(O)$ with $\mathrm{cof}(O)=\pm O$.
   - The "no decomposable vectors in the 21-block" statement is a dimension count, i.e. a heuristic. The doc now says so.
5. **Corollary 2.2, Proposition E1, Lemma 7.1** are re-derived. Two details:
   - Speiser: a semilinear action of $\mathrm{Gal}(\mathbb C(C)/\mathbb C(C)^K)=K$ on a one-dimensional $\mathbb C(C)$-space has invariant vectors.
   - The Mumford class of $L^{\otimes o}$ is $o$ times that of $L$.

## 3. Exact computations [X] (`verify_exact.py`)

These use the ATLAS table of `cover.py` (orthonormality re-verified symbolically) and power maps computed from
permutations:
- $\chi_{14a}(g)=\#\{\text{fixed 2-subsets}\}-\#\{\text{fixed points}\}$ for all 2520 $g$. So `cover.py`'s $14a$ is $14_{(5,2)}$; chartab's
  label for it is "14b".
- $\dim(\wedge^3V)^{A_7}$: $6{:}0$, $14_{(5,2)}{:}0$, $14_{(4,3)}{:}0$, $15{:}1$, $21{:}1$, $35{:}5$, $10\oplus\overline{10}{:}2$.
- $\wedge^2 14_{(5,2)}=10+\overline{10}+15+21+35$, with $\langle\chi,\chi\rangle=5$.
- $\mathrm{Ind}_{S_5}^{A_7}1=1+6+14_{(5,2)}$ ($S_5$ = stabiliser of $\{5,6\}$) and $\mathrm{Ind}_{A_6}^{A_7}1=1+6$.
- $Q_1$ sees $6,14_{(5,2)},14_{(4,3)},15,21$, each exactly once; $Q_2$ sees $10,\overline{10},15,35^2$.
- $2520/(6\cdot28)=15$ and $1296/15\notin\mathbb Z$; $1092/(2\cdot42)=13$ and $216/13\notin\mathbb Z$.

## 4. Certified eigenvalue inputs for $(2,4,7)$ classes 0, 1 [C] (`certify_th.py`, `certify_th_output.txt`)

| input | method | result (classes 0 and 1) |
|---|---|---|
| $\mu_1(Q_1)=\lambda_1(C)$ | `certificate.txt` ($n=96$) | $\ge0.34089$ |
| $\mu_2(Q_1)$ | **count**: verified $LDL^T$ inertia, $n=96$, $\sigma=0.56$: exactly 1 negative pivot | $\ge0.55998$ |
| $\lambda_2(C/S_5)$ | **upper** comparison problem, P1 Rayleigh–Ritz on $\mathrm{span}\{1,v\}$, $n=24$ | $\le0.36318$ |
| $\mu_1$ of isotype 6 | $C/A_6$, second CR eigenvalue, constants deflated, $n=16$ | $\ge0.99783$ |
| $\mu_1(Q_2)$, $\mu_2(Q_0)$ | `certificate.txt` | $\ge0.68902$, $\ge0.99783$ |

**Deduction.**
- Only one $Q_1$-eigenvalue lies below $0.55998$, and each $Q_1$ isotype occurs once. So exactly one copy of one isotype
  $\rho\in\{6,14_{(5,2)},14_{(4,3)},15,21\}$ has an eigenvalue below $0.55998$, and that eigenvalue is $\lambda_1$.
- $\rho\ne15$, because $Q_2$ sees $15$ and is $\ge0.689$.
- $C/S_5$ forces an eigenvalue $\le0.36318$ in $6$ or $14_{(5,2)}$, so $\rho\in\{6,14_{(5,2)}\}$. The $C/A_6$ bound excludes 6.

So:
- **$E_1$ is exactly one copy of $14_{(5,2)}$**;
- **$\lambda_1\in[0.34089,0.36318]$**;
- **every other eigenvalue is $\ge0.55998$**, i.e. $\lambda'\ge0.55998$.

In particular the hypothesis of Lemma 2.1 holds.

**The count step.** CHOLMOD's simplicial up-looking $LDL^T$ (no pivoting, symmetric fill-reducing permutation $P$) computes
for each row $k$ three things:
- the solve $L\hat y=a_{1:k-1,k}$;
- $\hat l_{kj}=\hat y_j/\hat d_j$;
- $\hat d_k=a_{kk}-\sum_j\hat l_{kj}\hat y_j$.

Standard inner-product and triangular-solve error analysis gives
$$P(B-cI)P^T+\Delta=\hat L\hat D\hat L^T,\qquad|\Delta|\le\gamma_{k+3}\,|\hat L||\hat D||\hat L^T|,$$
where $k$ is the maximum number of nonzeros in a row of $\hat L$. This uses $\hat y_j=\hat l_{kj}\hat d_j(1+\delta)^{-1}$ in both the
off-diagonal and the diagonal relations. Hence
$$\|\Delta\|_2\le\gamma_{k+3}\max|\hat D|\,\|\hat L\|_1\|\hat L\|_\infty .$$
If the shift $c$ exceeds this plus the assembly and rounding errors (`need`), then $B_{\rm exact}=\hat L\hat D\hat L^T+(\text{positive definite})$. Weyl
monotonicity and Sylvester's law then give $\#\{\text{negative eigenvalues of }B_{\rm exact}\}\le\#\{\hat d_i<0\}=1$. So
$\lambda_{2,h}\ge\sigma$, and Liu's bound for $k=2$ gives $\lambda_2(Q_1)\ge\sigma/(1+C_h^2\sigma)$.

| class | need | shift | $\min|\hat d_i|$ |
|---|---|---|---|
| 0 | $6.3\times10^{-9}$ | $1.06\times10^{-8}$ | 0.23 |
| 1 | $1.2\times10^{-8}$ | $2.9\times10^{-8}$ | 0.056 |

There was no pivot growth.

**The upper bound.**
- On each cell, $A_R\preceq c_e^{\max}P_e$ and $w_R\ge w_e^{\min}$ (ball arithmetic). So for every P1 function $u$ the comparison Rayleigh
  quotient bounds the true one from above.
- $\lambda_2\le\max_{\mathrm{span}\{1,v\}}R$ by min–max.
- Since $\nabla1=0$, this maximum is the closed form $k_{vv}m_{11}/(m_{11}m_{vv}-m_{1v}^2)$, evaluated in arb.

## 5. The two uncertified constants [N]

> **Update (Round-5 review, `ROUND5_REVIEW.md`).**
> - GPT's *harmonic Hersch* inequality replaces (TH) here. With it, the rigorous crude bound $\kappa^2\le b^*_{\max}/3$
>   suffices, so the quartic maximum below is no longer needed. The tolerance on $\sqrt{b^*_{\max}}$ and $\Gamma$ becomes 56%.
> - The $\gamma$ values below came from a buggy frame contraction ($LSL^T$ instead of $L^TSL$; now fixed). The
>   corrected values are smaller: $n=8$: $1.635$.
> - A complete sufficient trial-space certificate is GPT's Thm 5.3.

Both are computed for $E_1=14_{(5,2)}$ by P1 elements on the full curve:
- **$\kappa^2=0.6143\,b^*_{\max}/3$**, with $b^*_{\max}\approx1.277\times10^{-4}$ (Richardson on $n=8,12$; the change from $n=8$ to $n=12$ is 0.3%).
  The factor $0.6143$ is the best of 400 BFGS runs, with the top five agreeing to four digits. It is not certified to be the
  global maximum.
- **$\Lambda=\gamma^2/A$**, with $\gamma\approx1.70$ ($n=8{:}\ 1.692$, $n=12{:}\ 1.697$).

With the certified inputs of §4, (TH) excludes $m=24$ for every $\lambda_1\in[0.34089,0.36318]$, provided $\kappa$ and $\sqrt\Lambda$
exceed these values by at most a common factor of **1.035**. With $q_{\rm dec}$ in place of $\kappa$ it fails (0.966): the sharp 3-frame
constant is essential.

To certify, one needs:
1. upper bounds on $\kappa$ and $\Lambda$ that hold for the exact eigenspace, for example by the trial-space route of
   `FRAMEWORK_TOPOLOGICAL_HERSCH.md` §5.3;
2. a certified global maximum for the quartic $F(W)/|W|^4$ on $S^{41}$, for example by an SOS/Lasserre relaxation;
   the cruder $q_{\rm dec}$ is not enough.

A sharper $\lambda_1$ certificate raises the tolerance: $0.342\to1.08$, $0.344\to1.22$.

## 6. The ten other signatures [C] (`certify_signatures.py`, $n=24$; output in `certify_signatures_output.txt`)

Each curve is certified by $\lambda_1\ge\min(\mu_2(Q_0),\mu_1(Q_1),\mu_1(Q_2))$ on its own $(p,q,r)$ tiling. `certify.ref_triangle_arb` is now
general in $p$; for $p=2$ it is unchanged, and it agrees with the float placement to $4\times10^{-16}$. It is then compared with
$48/(g-1)$. All 24 curves (70 $A_7$-classes, matching the audit's counts) pass:

| signature | $g$ | $48/(g-1)$ | certified $\lambda_1\ge$ (worst curve) | Li–Yau $\operatorname{gon}\ge$ (worst) |
|---|---|---|---|---|
| $(3,3,5)$ | 169 | 0.2857 | 0.3530 | 29.65 |
| $(2,5,7)$ | 199 | 0.2424 | 0.2607 | 25.81 |
| $(3,3,6)$ | 211 | 0.2286 | 0.3584 | 37.64 |
| $(3,4,4)$ | 211 | 0.2286 | 0.2481 | 26.05 |
| $(2,6,7)$ | 241 | 0.2000 | 0.2565 | 30.78 |
| $(3,3,7)$ | 241 | 0.2000 | 0.2665 | 31.99 |
| $(2,7,7)$ | 271 | 0.1778 | 0.2126 | 28.71 |
| $(3,4,5)$ | 274 | 0.1758 | 0.2593 | 35.40 |
| $(3,4,6)$ | 316 | 0.1524 | 0.2543 | 40.06 |
| $(4,4,4)$ | 316 | 0.1524 | 0.2599 | 40.94 |

**Correction.** The fourth session's survey (P1, $n=4$, divided by $1.0105$) overestimated several $\lambda_1$. Examples:
- $(4,4,4)$: survey $0.313$, but the comparison CR value at $n=24$ is $0.262$;
- $(2,7,7)$: survey $0.249$, against $0.214$.

The $1.0105$ calibration from $(2,4,7)$ does not transfer to other triangles. The certified table supersedes the survey.
Every certified margin is still at least 7% ($(2,5,7)$: $0.2607$ against $0.2424$).

## 7. $(2,4,7)$ classes 12, 14 [C]

$Q_1$ at $n=128$ with $\sigma=0.3557$ (verified Cholesky, need $1.7\times10^{-10}$ against shift $3.5\times10^{-10}$) gives
$\mu_1(Q_1)\ge0.355696$. `certificate.txt` gives $Q_2\ge0.76298$ and $Q_0\ge0.99783$. So
$$\lambda_1\ge0.355696>48/135=0.355556,\qquad\operatorname{gon}\ge67.5\cdot0.355696=24.0095,$$
i.e. $\operatorname{gon}\ge25$. The margin is small (0.04%) but rigorous. Class 14
gives the same bound (need $1.7\times10^{-10}$ against shift $3.4\times10^{-10}$). Commands and outputs are in `certify_th_output.txt`.

## 8. Trust base (unchanged from `certificate.txt`, plus one item)

- **Classical.** Hersch balancing and the Li–Yau/Yang–Yau inequality (NOTES 7.2); the CR lower-bound theorem (Liu 2015;
  Carstensen–Gedicke–Rim 2012, Lemma 2.2); Courant–Fischer, Weyl and Sylvester.
- **Model.** The exact hyperbolic tiling model (NOTES 7.3, now also for $p\ne2$); IEEE binary64 with round-to-nearest; the
  correctness of python-flint (arb) and CHOLMOD.
- **New.** The $\gamma_{k+3}$ backward-error bound for CHOLMOD's simplicial $LDL^T$ (derived in §4, not taken from a reference).
- **Not verified here.** GPT's algebra for $g\ge336$; the transport degrees 1296 and 216; the transfer
  $\operatorname{gon}\Rightarrow\mathrm{ed}$ (audit, page 8).
