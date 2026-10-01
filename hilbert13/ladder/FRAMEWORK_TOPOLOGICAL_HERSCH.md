# Topological Hersch: the first eigenspace carries no degree

Status tags as in NOTES.md:
- **[P]** proved here;
- **[N]** numerical (P1 finite elements in the exact geometry of `orbifold.py`, not certified);
- **[?]** expectation;
- **[Q]** question.

Code: `topo_hersch.py` (constants and the worst-case bound), `topo_hersch_extras.py` (the $\wedge^3$ table, the
sharpened constant, the certification targets), `signatures_spectrum.py` (the survey in §6). Outputs are in
`topo_hersch_output.txt`.

## 0. Summary

**The barrier.** Hersch/Li–Yau gives $\operatorname{gon}(C)\ge\lambda_1A/8\pi$. For the $(2,4,7)$ curves of classes 0 and 1 this is
$67.5\,\lambda_1\approx23.37$. So no argument that uses only $\lambda_1$ and the area can reach $\operatorname{gon}\ge25$ there.
The conformal-metric escape (`FRAMEWORK_CONFORMAL.md` §A) is numerically dead: see §E of that file, the gain is
$\le0.2\%$.

**The new mechanism [P].**
- The degree of a map $x:C\to S^2$ is a *cubic* functional, $T(x)=\int x\cdot(x_s\times x_t)=4\pi\deg x$.
- Restricted to maps whose three coordinates lie in the first eigenspace $E_1$, $T$ is given by an invariant
  alternating 3-form on $E_1$. For the $(2,4,7)$ curves $E_1\cong14_{(5,2)}$ (the notes' $14_a$), and
  $(\wedge^3 14_{(5,2)})^{A_7}=0$. So **every map with coordinates in $E_1$ has degree 0**.
- A balanced holomorphic map of degree $m$ has energy $8\pi m=\lambda_1A+\varepsilon$. It must leave $E_1$, and
  expanding $T$ shows that it can only gain degree at rate $O(\sqrt\varepsilon)$.
- This forces $\varepsilon$ to be large, i.e. $8\pi m$ to exceed $\lambda_1A$ by a definite amount.

**The result.** For $(2,4,7)$ classes 0 and 1 the inequality excludes $m=24$ with a clear margin (at $n=12$: right
side 233 against 301.6 needed, worst case over all admissible maps). Classes 12 and 14 reach 24.28 by hyperbolic
Li–Yau. Together:
$$\operatorname{gon}(C)\ \ge\ 25\qquad\text{for every }(2,4,7)\ A_7\text{-curve.}$$
This is the first bound past the $23.37$ barrier.

**Rigour** (fifth session; see `VERIFICATION.md`):
- **Classes 12, 14:** certified ($\lambda_1\ge0.355696$).
- **Classes 0, 1:** the theorem is checked on paper. All eigenvalue inputs are certified: $E_1$ is one copy of $14_{(5,2)}$,
  $\lambda_1\in[0.34089,0.36318]$, and $\lambda'\ge0.55998$. Only the two eigenfunction constants $\kappa$ and $\Lambda$ of §3 remain
  numerical. They may exceed their computed values by up to 3.5%.
- **The ten other rigid signatures:** certified to have $\operatorname{gon}\ge25$ (§6).

## 1. The cubic form

Let $X$ be a closed oriented surface with a conformal metric and area $A$. For $a,b,c:X\to\mathbb R^3$ put
$$\Theta(a,b,c)=\tfrac12\sum_{ijk}\epsilon_{ijk}\int_X a_i\,db_j\wedge dc_k,
\qquad T(x)=\Theta(x,x,x)=\int_X x\cdot(x_s\times x_t)\,ds\,dt ,$$
where $(s,t)$ are oriented local coordinates.

**Lemma 1.1 [P].** $\Theta$ is totally symmetric. Hence
$$T(y+z)=T(y)+3\Theta(z,y,y)+3\Theta(y,z,z)+T(z).$$
*Proof.* Symmetry in $(b,c)$: swap the wedge factors and relabel $j\leftrightarrow k$; the two sign changes cancel.
Symmetry in $(a,b)$: Stokes gives $\int a_i\,db_j\wedge dc_k=-\int b_j\,da_i\wedge dc_k$. Relabelling $i\leftrightarrow j$ in
$\epsilon_{ijk}$ gives a second sign. $\blacksquare$

**Lemma 1.2 [P].** For holomorphic $x:X\to\mathbb P^1=S^2$ of degree $m$:
- $T(x)=\int x^*\omega_{S^2}=4\pi m$;
- $E(x)=\int|\nabla x|^2=8\pi m$.

For an arbitrary map, $E(x)-2T(x)=\int|x_t-x\times x_s|^2\ge0$ (Bogomolny).

## 2. Vanishing on the first eigenspace

**Lemma 2.1 [P].** Let a finite group $G$ act on $X$ by orientation-preserving isometries. Let $E\subset C^\infty(X)$ be a
real $G$-subrepresentation with $(\wedge^3E)^G=0$. Then:
- $T\equiv0$ on $E\otimes\mathbb R^3$;
- all Poisson brackets $\{\varphi,\psi\}$ with $\varphi,\psi\in E$ are $L^2$-orthogonal to $E$.

*Proof.* Take an orthonormal basis $\varphi_i$ of $E$ and write $y=\sum_i\varphi_iW_i$ with $W_i\in\mathbb R^3$. Then
$$T(y)=\tfrac12\sum_{ijk}\theta_{ijk}\det(W_i,W_j,W_k),\qquad \theta_{ijk}=\int\varphi_i\,d\varphi_j\wedge d\varphi_k .$$
By Stokes $\theta$ is alternating. It is $G$-invariant because each $g\in G$ is an orientation-preserving
diffeomorphism, so $\int g^*(\varphi\,d\psi\wedge d\chi)=\int\varphi\,d\psi\wedge d\chi$. Hence $\theta\in(\wedge^3E^*)^G=0$. Finally,
$\langle\{\varphi_j,\varphi_k\},\varphi_i\rangle=\theta_{ijk}$. $\blacksquare$

**$A_7$ table [C]** (`topo_hersch_extras.py`, from the character table):

| real irrep $V$ | 6 | $14_{(4,3)}$ | $14_{(5,2)}$ | 15 | 21 | 35 | $10\oplus\overline{10}$ |
|---|---|---|---|---|---|---|---|
| $\dim(\wedge^3V)^{A_7}$ | **0** | **0** | **0** | 1 | 1 | 5 | 2 |

The same computation gives
$$\wedge^2 14_{(5,2)}=10+\overline{10}+15+21+35 ,$$
which is multiplicity-free.

**Corollary 2.2 [P] (Li–Yau is never sharp here).** Suppose $(\wedge^3E_1)^G=0$. Then every map $C\to S^2$ whose
coordinates are first eigenfunctions has degree 0. In particular equality in Li–Yau is impossible. Equality
would force the balanced map to lie in $E_1\otimes\mathbb R^3$. Compare the observation of the last session that the
near-spherical projections $y/|y|$, $y\in E_1\otimes\mathbb R^3$, all have degree 0: this lemma is the reason.

## 3. The inequality

**Setting.**
- $X$ is a closed Riemann surface with a conformal metric of area $A$.
- $E_1$ is the $\lambda_1$-eigenspace, with $T\equiv0$ on $E_1\otimes\mathbb R^3$.
- $\lambda'$ is the next eigenvalue, and $\delta=\lambda'-\lambda_1$.
- $\varphi_1,\dots,\varphi_k$ is an orthonormal basis of $E_1$.

**The two constants.**
- **Bracket constant.** For $\omega\in\wedge^2E_1$ put $\beta(\omega)=\sum_{i<j}\omega_{ij}\{\varphi_i,\varphi_j\}$. Define
  $Q^*(\omega)=\langle\beta(\omega),(\Delta-\lambda_1)^{-1}\beta(\omega)\rangle$, with the resolvent taken on $(1\oplus E_1)^\perp$; this is
  legitimate by Lemma 2.1. Then set
  $$\kappa^2=\sup\Big\{\,Q^*(w_2\wedge w_3)+Q^*(w_3\wedge w_1)+Q^*(w_1\wedge w_2)\ :\ w_a\in E_1,\ \textstyle\sum_a\|w_a\|^2=1\Big\}.$$
- **Gradient constant.** $\Lambda=\max_p\lambda_{\max}S(p)$, where $S(p)=\sum_i\nabla\varphi_i(p)\otimes\nabla\varphi_i(p)$ is a $2\times2$
  matrix in an orthonormal frame.

**Theorem 3.1 [P].** Let $x:X\to\mathbb P^1$ be holomorphic of degree $m$. Put $\varepsilon=8\pi m-\lambda_1A$. Then for some
$s\in[0,\varepsilon/\delta]$,
$$\frac{\lambda_1(A-s)}2\ \le\ 3\kappa\sqrt\varepsilon\,(A-s)+\sqrt{\Lambda(A-s)}\,\sqrt{s(\varepsilon+\lambda_1s)} .\tag{TH}$$
If $6\kappa\sqrt\varepsilon<\lambda_1$ and $\varepsilon/\delta\le A/2$, the right side minus the left side increases with $s$ on
$[0,\varepsilon/\delta]$. Then $m$ is excluded as soon as (TH) fails at $s=\varepsilon/\delta$.

*Proof.*
1. **Split.** Balance $x$ by a Möbius map (Hersch); it stays holomorphic of degree $m$. Write $x=y+z$, where
   $y\in E_1\otimes\mathbb R^3$ is the componentwise $L^2$ projection. Since $x$ has mean zero, $z\perp1\oplus E_1$.
2. **Energy.** With $s=\|z\|^2$:
   - $\|y\|^2=A-s$;
   - $8\pi m=E(x)=\lambda_1\|y\|^2+E(z)$;
   - $E(z)\ge\lambda'\|z\|^2$.

   Hence $E(z)=\varepsilon+\lambda_1s$ and $s\le\varepsilon/\delta$.
3. **Expansion.** By Lemmas 1.1, 1.2 and 2.1,
   $$4\pi m=T(x)=3\Theta(z,y,y)+2\Theta(z,y,z)+\Theta(x,z,z).$$
   This uses $3\Theta(y,z,z)+T(z)=2\Theta(z,y,z)+\Theta(x,z,z)$.
4. **The $\Theta(x,z,z)$ term.** Since $|x|=1$,
   $$|\Theta(x,z,z)|=\Big|\int x\cdot(z_s\times z_t)\Big|\le\int|z_s||z_t|\le E(z)/2 .$$
5. **The $2\Theta(z,y,z)$ term.** $2\Theta(z,y,z)=\int z\cdot(y_s\times z_t-y_t\times z_s)$. Take $a=y_s$, $b=y_t$. The linear
   map $(u,v)\mapsto a\times v-b\times u$ has operator norm exactly $(|a|^2+|b|^2)^{1/2}$: its Gram matrix is
   $(|a|^2+|b|^2)I-aa^T-bb^T$, and $aa^T+bb^T$ is singular in $\mathbb R^3$. Also $|\nabla y(p)|^2\le\lambda_{\max}S(p)\,\|y\|^2$.
   Hence
   $$|2\Theta(z,y,z)|\le\sqrt{\Lambda\|y\|^2}\,\|z\|\,\|\nabla z\| .$$
6. **The $3\Theta(z,y,y)$ term.** $\Theta(z,y,y)=\langle z,N_y\rangle$, where
   $$N_y=\sum_{i<j}\{\varphi_i,\varphi_j\}\,W_i\times W_j .$$
   Its $k$-th component is $\beta(w_a\wedge w_b)$ with $(a,b,k)$ cyclic, where $w_a\in E_1$ is the $a$-th column of $W$.
   Cauchy–Schwarz in the $(\Delta-\lambda_1)$ inner product on $(1\oplus E_1)^\perp$ gives
   $$|\langle z,N_y\rangle|\le\langle z,(\Delta-\lambda_1)z\rangle^{1/2}\,\|N_y\|_*=\sqrt\varepsilon\,\|N_y\|_*,\qquad \|N_y\|_*\le\kappa\|y\|^2 .$$
7. **Combine.** Substitute $4\pi m=(\lambda_1A+\varepsilon)/2$ and $E(z)=\varepsilon+\lambda_1s$. $\blacksquare$

**Remark 3.2 (computing $\kappa$) [P].** $Q^*$ is $G$-invariant on $\wedge^2E_1$. When $\wedge^2E_1$ is multiplicity-free, Schur gives
$Q^*=\sum_\sigma b^*_\sigma P_\sigma$. Then:
- **Crude:** $\kappa^2\le b^*_{\max}/3$, using $e_2(\sigma^2)\le(\sum\sigma^2)^2/3$ for the singular values of $W$.
- **Better:** $\kappa^2\le q_{\rm dec}/3$, where $q_{\rm dec}$ is the maximum of $Q^*$ on *decomposable* unit 2-vectors. Reduce $W$ by its SVD to
  orthogonal columns.

A dimension count predicts that the top isotype contains no decomposable 2-vectors: the cone over
$\mathrm{Gr}(2,14)$ has dimension 25, and the $21$-block has codimension 70 in $\wedge^2E_1$. Then $q_{\rm dec}<b^*_{\max}$. Numerically
$q_{\rm dec}=0.819\,b^*_{\max}$ (§4).

## 4. Numbers for $(2,4,7)$ [N]

Full curve, P1 elements, $n^2$ cells per triangle (`topo_hersch.py`). The bracket blocks are listed as
$(\dim,\ b^*_\sigma=b_\sigma\cdot\text{ratio}_\sigma)$. The "ratio" is the dual-norm factor
$\langle\beta,(\Delta-\lambda_1)^{-1}\beta\rangle/\|\beta\|^2$.

| class, $n$ | $\lambda_1$ | $\lambda'$ | $21$ | $15$ | $35$ | $10+\overline{10}$ | $\gamma=\sqrt{\Lambda A}$ |
|---|---|---|---|---|---|---|---|
| 0, $n=8$ | 0.34725 | 0.57293 | 1.2695e-4 (1.229) | 4.2675e-5 (1.019) | 6.756e-6 (0.491) | 2.559e-6 (0.245) | 1.692 |
| 0, $n=12$ | 0.34671 | 0.57189 | 1.2736e-4 (1.236) | 4.2806e-5 (1.028) | 6.792e-6 (0.496) | 2.568e-6 (0.248) | 1.697 |
| 12, $n=8$ | 0.36058 | 0.38746 | 7.222e-5 (4.772) | 1.030e-5 (0.537) | 2.778e-5 (0.497) | 3.633e-5 (0.605) | 1.770 |

Observations:
- The frame energy density $\sum|\nabla\varphi_i|^2$ lies in $[0.00284,0.00296]$, against a mean of $14\lambda_1/A=0.00286$.
- The frame is nearly conformal: $2\lambda_{\max}/\mathrm{tr}\,S\in[1.012,1.174]$.

**Sharpened constant.** At $n=8$, over 400 BFGS starts, the five best local maxima agree to four digits:
$$\kappa^2=0.6143\cdot b^*_{\max}/3,\qquad q_{\rm dec}=0.8191\cdot b^*_{\max} .$$

**Worst-case exclusion of $m=24$ (class 0).** These use the unsharpened form $\|y\|^2\le A$ and $\kappa^2=b^*_{\max}/3$, with
$4\pi m=301.6$:

| $n$ | right side = $3\sqrt\varepsilon\,\|N\|_*+\gamma\sqrt{sE(z)}+E(z)/2$ | $m=24$ | $m=25$ (314.2) |
|---|---|---|---|
| 8 | $124.3+80.0+17.9=222.2$ | excluded | 479.7, not excluded |
| 12 | $128.5+85.5+19.1=233.1$ | excluded | 489.8, not excluded |

**Extrapolated values.** Richardson on $n=8,12$ gives $b^*_{\max}\approx1.2769\times10^{-4}$ and $\gamma\approx1.70$; the true
$\lambda_1=0.34627$, $\lambda'=0.5715$. In the cleaned form (TH):
- left side $\lambda_1(A-s)/2=281.6$;
- right side $214.4$ (crude $\kappa$), $202.4$ (via $q_{\rm dec}$), $187.0$ (sharp $\kappa$).

**Classes 12, 14.** Here $\lambda'\approx0.386$ is close to $\lambda_1\approx0.3597$, so $s\le\varepsilon/\delta$ is weak, and (TH) adds
little: at $m=25$ the worst-case right side is $902$ against $314$. They do not need it, because Li–Yau already gives
$24.28$. Their bracket spectrum differs from classes 0, 1: the top $b^*$ is in $21$ (ratio $4.77$), because the $21$
eigenvalue is near.

## 5. What a proof needs

**Eigenvalue thresholds for $m=24$.** These are the smallest certified values that make (TH) fail for every
$\lambda_1$ in [certified lower bound, 0.34671], with the extrapolated $\kappa$ and $\Lambda$ (`topo_hersch_extras.py`):

| $\kappa^2$ used | need $\lambda_1\ge$ (with $\lambda'\ge0.56$) | (with $\lambda'\ge0.55$) |
|---|---|---|
| crude $b^*_{\max}/3$ | 0.34246 | 0.34278 |
| $q_{\rm dec}/3$ | 0.34154 | 0.34189 |
| sharp $0.6143\,b^*_{\max}/3$ | **0.34024** | **0.34065** |

The certified $\lambda_1\ge0.34089$ (`certificate.txt`) already clears the sharp row. With the sharp $\kappa$ and $\lambda'\ge0.55$,
$\kappa$ and $\sqrt\Lambda$ can be inflated by a common factor:

| certified $\lambda_1\ge$ | 0.34089 | 0.342 | 0.344 | 0.3455 |
|---|---|---|---|---|
| allowed factor | 1.013 | 1.079 | 1.222 | 1.361 |

So a slightly sharper $\lambda_1$ certificate buys a comfortable 20% tolerance for the eigenfunction constants.

**The items.** *(Fifth session: items 1 and the isotype identification are done, in `VERIFICATION.md` §4. With them,
the tolerance on $\kappa$ and $\sqrt\Lambda$ is 3.5%.)*
1. **$\lambda'\ge0.55$** is the second eigenvalue of $Q_1$ (the $21$). $Q_2\ge0.689$ and $Q_0\ge0.998$ are already certified.
   The CR lower bound (Liu) holds for the $k$-th eigenvalue. The count needs one extra step. Positive
   definiteness of $K-\sigma M+\rho\,(Mv)(Mv)^T$, for an approximate eigenvector $v$, proves that at most one discrete
   eigenvalue is below $\sigma$. This is the same verified-Cholesky machinery.
2. **$\lambda_1\ge0.344$** (optional, for tolerance). The comparison loss of `certify.py` is first order in $h$: 1.55% at
   $n=96$. So $0.344$ would need $n\approx230$ with the present method. A second-order comparison, with piecewise-linear
   coefficient bounds or exact weights with interval quadrature, is the realistic path. The same upgrade gives classes
   12, 14 their $\lambda_1>0.35556$.
3. **$\kappa$ and $\Lambda$: trial-space version [?].** Replace $E_1$ in §3 by the *discrete* eigenspace $E_h$, which is
   exactly known: piecewise linear on a $G$-symmetric mesh, and isomorphic to $14_{(5,2)}$. Then:
   - Lemma 2.1 still gives $T\equiv0$ on $E_h\otimes\mathbb R^3$.
   - $\Lambda_h$ is a finite maximum of piecewise constants, checkable in ball arithmetic.
   - The price is two correction terms: a cross term $\int\nabla y\cdot\nabla z$, bounded by the $H^{-1}$ residual of $E_h$,
     and an eigenspace angle, which enters the bound $E(z)\ge\nu\|z\|^2$ on $(1\oplus E_h)^\perp$. The angle is bounded by
     Davis–Kahan with the certified gap $\delta\ge0.2$.
   - $\kappa_h$ needs an *upper* bound on the resolvent form $\langle f,(\Delta-\lambda_1)^{-1}f\rangle$. The CR/Liu projection-error
     constants give it; the crude isotypic bound $1/(\mu_\sigma-\lambda_1)$ loses a factor of about 3.6 and is not enough.

## 6. The other rigid signatures [C]

**Superseded.** The fourth-session survey (P1 at $n=4$, divided by $1.0105$) overestimated several $\lambda_1$. The calibration
does not transfer between triangles.

**Certified version** (`certify_signatures.py`, `certify_signatures_output.txt`, `VERIFICATION.md` §6).
- Each of the 24 curves is certified by $Q_0,Q_1,Q_2$ on its own $(p,q,r)$ tiling at $n=24$. These cover all 70 $A_7$-classes.
- Every curve satisfies $\lambda_1>48/(g-1)$, hence $\operatorname{gon}\ge25$.
- The tightest rows are $(2,5,7)$ ($\lambda_1\ge0.2607$ against $0.2424$) and $(3,4,4)$ ($0.2481$ against $0.2286$).

Together with GPT's algebraic range $g\ge336$ (`FRAMEWORK_CONFORMAL.md` D5, not re-verified), this reduces
$\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$ to certifying $\kappa$ and $\Lambda$ for $(2,4,7)$ classes 0, 1.

## 7. The invariant degree lattice [P]

**Lemma 7.1.** Let $K$ act faithfully on a curve $C$. Let $\ell_K$ be the lcm of the orders of the point stabilisers.
1. If $L$ is a $K$-linearised line bundle, then $\deg L\in(|K|/\ell_K)\mathbb Z$.
2. If $L$ is only $K$-invariant, then $o\cdot\deg L\in(|K|/\ell_K)\mathbb Z$ for the order $o$ of its Mumford class in
   $H^2(K,\mathbb C^\times)$; $o$ divides $\exp M(K)$.

*Proof.* By Hilbert 90 (Speiser), $H^1(K,\mathbb C(C)^\times)=0$. So the semilinear $K$-action on the one-dimensional
$\mathbb C(C)$-space of rational sections has a nonzero fixed vector $\sigma$. Its divisor is $K$-invariant, hence a
$\mathbb Z$-combination of orbits. An orbit has size $|K|/|K_p|$, which is divisible by $|K|/\ell_K$. For (2), apply (1) to
$L^{\otimes o}$. $\blacksquare$

**Consequences.**
- **Pencils.** For a pencil $f$ with stabiliser $K$, $f^*\mathcal O(2)=f^*\omega_{\mathbb P^1}^{-1}$ is $K$-linearised. So
  $|K|/\ell_K$ divides $2m$.
- **Transport as a lattice violation.** For $G=A_7$ with signature $(2,4,7)$: $|G|/(\exp M\cdot\operatorname{lcm})=2520/(6\cdot28)=15$,
  and the transport class has degree $1296\notin15\mathbb Z$. For $L_2(13)$ with $(2,3,7)$: $1092/(2\cdot42)=13$, and
  $216\notin13\mathbb Z$.

## 8. Negative results recorded for completeness [N]

- **Conformal route.** See `FRAMEWORK_CONFORMAL.md` §E.
- **Quadric gap.** Put $D=\min_{y}\frac1A\int(1-|y|)^2$ over $y\in E_1\otimes\mathbb R^3$ with $\|y\|^2=A$. It is $\approx0.016$, which gives only
  $\operatorname{gon}\ge23.7$ through $\varepsilon\ge\delta A D$. The obstruction is the degree, not the shape.
