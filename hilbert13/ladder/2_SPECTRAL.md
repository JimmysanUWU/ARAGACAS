# Chapter 2. Gonality from a certified spectral gap

**Main results.**
- Every $(2,4,7)$ $A_7$-curve has $\lambda_1(C)\ge0.34089$. Hence $\operatorname{gon}(C)\ge24$ and $\operatorname{gon}(C/\langle\tau\rangle)\ge12$.
- Classes 12 and 14 have $\lambda_1\ge0.355696$. Hence $\operatorname{gon}\ge25$.
- Every faithful $A_7$-curve of the ten other rigid signatures of genus $\le335$ has $\operatorname{gon}\ge25$.
- For classes 0 and 1, the first eigenspace is located exactly; this is the input of Chapter 3.

Status tags: **[P]** paper proof; **[C]** computer-assisted proof (ball arithmetic plus rigorous floating-point error bounds);
**[L]** checked in Lean; **[N]** numerical only.

## 2.1 Gonality from $\lambda_1$

**Lemma 2.1 (Hersch; Yang–Yau 1980; Li–Yau 1982) [P][L].** Let $X$ be a compact Riemann surface with a conformal metric, and
$\varphi:X\to\mathbb P^1$ holomorphic of degree $d$. Then $\lambda_1(X)\operatorname{Area}(X)\le8\pi d$.

*Proof.*
1. Identify $\mathbb P^1$ with $S^2\subset\mathbb R^3$.
2. Hersch's balancing lemma gives a Möbius map $\gamma$ with $\int_X\gamma\circ\varphi=0$.
3. Each coordinate $x_i$ of $\gamma\circ\varphi$ is then admissible, so $\lambda_1\int x_i^2\le\int|\nabla x_i|^2$.
4. Summing over $i$, with $\sum x_i^2=1$ and $\sum|\nabla x_i|^2=2\operatorname{Jac}$ for a conformal map, gives $\lambda_1\operatorname{Area}\le2\cdot4\pi d$. $\square$

With the hyperbolic metric, $\operatorname{Area}(C)=4\pi(g-1)$. So
$$\operatorname{gon}(C)\ge\tfrac12(g-1)\lambda_1(C),\qquad\text{i.e. }\operatorname{gon}(C)\ge25\iff\lambda_1>48/(g-1)\text{ suffices}.$$
For $(2,4,7)$ this reads $\operatorname{gon}(C)\ge67.5\,\lambda_1$ and $\operatorname{gon}(D)\ge33.75\,\lambda_1$; the second holds because a degree-$d$ map on $D$ lifts to degree $2d$ on $C$.

## 2.2 An exact model

**Tiling.** $C$ is tiled by $2|G|$ copies of the hyperbolic triangle with angles $\pi/p,\pi/q,\pi/r$.
- The tiles are $U_g,L_g$ for $g\in G$.
- The gluings are $U_g\sim L_g$, $U_g\sim L_{gb}$ and $U_g\sim L_{ga}$, along the three sides.
- $G$ acts on the left.

**Chart.** Every tile is identified isometrically with one triangle in the Klein model, parametrised affinely by the standard simplex $R$.
- Tiles sharing a side are mirror images across it, so **every gluing is the identity on the sides of $R$**.
- In $R$-coordinates the Dirichlet form and the mass are exactly
  $$\int\nabla u\cdot A_R\nabla u,\qquad\int u^2w_R,$$
  with $A_R=|\det J|J^{-1}\frac{I-xx^T}{\sqrt{1-|x|^2}}J^{-T}$ and $w_R=|\det J|(1-|x|^2)^{-3/2}$.
- A uniform subdivision of $R$ is therefore a conforming mesh of $C$ and of all its quotients (`orbifold.py`).

## 2.3 Reduction to three small twisted problems

For $K\le G$ and a character $\varepsilon:K\to\{\pm1\}$, the $\varepsilon$-twisted problem is the Laplacian on $\{f:f(kx)=\varepsilon(k)f(x)\}$.
It lives on $2[G:K]$ tiles with signed gluings.

**Lemma 2.2 [P].** If an eigenspace of $C$ contains an irreducible $\rho$ with $\langle\rho|_K,\varepsilon\rangle>0$, then its eigenvalue is an eigenvalue of the
$\varepsilon$-twisted problem.

| problem | $K$ | tiles | $\varepsilon$ | irreducibles seen |
|---|---|---|---|---|
| $Q_0$ | $A_7$ | 2 | trivial | $1$ (first nonzero eigenvalue) |
| $Q_1$ | $3^2{:}4$ | 140 | sign | $6,14_a,14_b,15,21$ |
| $Q_2$ | $C_G((16)(23))\cong C_3\rtimes D_8$ | 210 | sign | $10,\overline{10},15,35$ |

Together they see every irreducible (`cover.py`, exact). Hence
$$\lambda_1(C)\ge\min\big(\mu_2(Q_0),\mu_1(Q_1),\mu_1(Q_2)\big).$$

## 2.4 Guaranteed lower bounds

1. **Comparison [C].** On each mesh cell $e$, ball arithmetic gives $A_e=c_eP_e\preceq A_R$ (with $P_e$ the centroid value) and $w_e\ge w_R$. By min–max, the
   eigenvalues of the true problem dominate those of this piecewise-constant problem.
2. **Crouzeix–Raviart bound [P][L]** (Liu 2015; Carstensen–Gedicke 2014). With
   $C_h^2=\max_e\kappa_e^2w_e/\lambda_{\min}(A_e)$ and $\kappa_e^2=\frac{(1/n)^2}8+\frac{(\sqrt2/n)^2}{j_{1,1}^2}$,
   $$\lambda_k\ \ge\ \frac{\lambda_{k,h}}{1+C_h^2\lambda_{k,h}} .$$
   The CR interpolation is consistent across tiles, because the gluings are the identity.
3. **Verified positive definiteness [C][L].** The CR mass matrix is diagonal, so $\lambda_{1,h}>\sigma$ iff $B=K_h-\sigma M_h\succ0$.
   - A floating sparse Cholesky of $\hat B-cI$ runs to completion.
   - Higham's backward error bound gives $|\Delta|\le\gamma_{k+1}|\tilde L||\tilde L^T|$, together with the assembly error $\eta$.
   - So $B_{\rm exact}\succeq\big(c-\gamma_{k+1}\|\tilde L\|_1\|\tilde L\|_\infty-u\max|b_{ii}|-\eta\big)I\succ0$.
   - For $Q_0$, the constants are deflated by a rank-one term.
4. **Counting [C]** (`certify_th.py`). To bound the *second* eigenvalue, use CHOLMOD's simplicial $LDL^T$ without pivoting.
   - It satisfies $P(B-cI)P^T+\Delta=\hat L\hat D\hat L^T$ with $|\Delta|\le\gamma_{k+3}|\hat L||\hat D||\hat L^T|$.
   - Sylvester and Weyl give $\#\{\text{negative eigenvalues of }B\}\le\#\{\hat d_i<0\}$.
5. **Upper bounds [C].** Rayleigh–Ritz on $\operatorname{span}\{1,v\}$ for the reverse comparison problem ($A_R\preceq c_e^{\max}P_e$, $w_R\ge w_e^{\min}$).

## 2.5 Results

**$(2,4,7)$, all classes [C]** (`run_certificate.py` → `certificate.txt`; $Q_1$ at $n=96$, $\sigma=0.3409$):
- $\mu_2(Q_0)\ge0.99783$;
- $\mu_1(Q_1)\ge0.34089$;
- $\mu_1(Q_2)\ge0.68902$ (classes 0, 1) and $0.76298$ (classes 12, 14).

Hence $\lambda_1\ge0.34089$, and so $\operatorname{gon}(C)\ge24$ and $\operatorname{gon}(D)\ge12$. The Lean file `Hilbert13/SpectralCertificate.lean`
checks three pieces: the abstract CR bound, the positive-definiteness criterion with permutation, and the final arithmetic.

**Classes 12, 14 [C]** (`certify_th_output.txt`). $Q_1$ at $n=128$ with $\sigma=0.3557$ gives
$\lambda_1\ge0.355696>48/135$. So $\operatorname{gon}\ge25$.

**Classes 0, 1: the first eigenspace [C]** (`certify_th.py`, `certify_th_output.txt`).
- **Count.** The count on $Q_1$ ($n=96$, $\sigma=0.56$) finds exactly one negative pivot. So $\mu_2(Q_1)\ge0.55998$.
- **Upper bound.** The upper bound on $C/S_5$ gives $\lambda_2(C/S_5)\le0.36318$. Here $\operatorname{Ind}_{S_5}^G1=1+6+14_a$.
- **Excluding the 6.** The CR bound on $C/A_6$, whose functions see only $1+6$, gives $\ge0.99783$.

So **$E_1$ is exactly one copy of $14_a=14_{(5,2)}$**, **$\lambda_1\in[0.34089,0.36318]$**, and **every other eigenvalue is $\ge0.55998$**.
(Only one $Q_1$-isotype can lie below 0.56. $Q_2$ excludes $15$, $C/S_5$ forces $6$ or $14_a$, and $C/A_6$ excludes $6$.)

**The ten other rigid signatures of genus $\le335$ [C]** (`certify_signatures.py 24 ...` → `certify_signatures_output.txt`):
24 curves, 70 classes.

| signature | $g$ | $48/(g-1)$ | certified $\lambda_1\ge$ (worst) | Li–Yau gon $\ge$ |
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

These ten signatures together with $(2,4,7)$ are all faithful $A_7$-curves of genus $\le335$. Four-point signatures have $g\ge421$, by the
seven-sheeted Riemann–Hurwitz test, and quotient genus $h\ge1$ gives $g\ge631$. The four rigid signatures with $336\le g\le397$ are
certified in Chapter 4.

## 2.6 Accurate values and validation

**$(2,4,7)$ [N].** $\lambda_1=0.346267085404$ for classes 0, 1 and $0.359671354895$ for classes 12, 14, with multiplicity 14 (isotype $14_a$).
These are vector-valued Hejhal expansions to 12 digits (`hejhal_solve.py`, Chapter 3). For classes 0, 1 the value is also certified to $2.5\cdot10^{-4}$.

**Bolza surface (external check)** (`validate_bolza.py` → `validation_bolza.txt`). The same pipeline on the $(2,3,8)$ curve of genus 2 certifies
$3.540,\ 3.701,\ 3.772,\ 3.804$ at $n=8,16,32,64$. All are below the known $\lambda_1=3.8388872588$ (Strohmaier–Uski), and they converge to it.

## 2.7 Limits

Li–Yau gives $\operatorname{gon}(C)\ge\lceil67.5\lambda_1\rceil$. For classes 0, 1 that is $24$, because $67.5\cdot0.34627=23.37$.

- **Changing the conformal metric does not help.** The gain is at most $\lambda_1A/\min\bar F$, about $0.2\%$ (§6.3).
- **The way past 24 uses the degree, not the energy.** The degree of a map to $S^2$ is a cubic form that vanishes on the first eigenspace by
  representation theory. Chapter 3 turns this into a proof.
