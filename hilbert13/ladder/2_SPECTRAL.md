# Chapter 2. Gonality from a certified spectral gap

**Idea.**
- **Degree from energy.** A map of degree $d$ to $\mathbb P^1$ gives, after balancing, three functions with total energy $8\pi d$ and total mass equal to the area. So the first eigenvalue bounds the degree from below (Li–Yau).
- **The work** is to certify $\lambda_1$ on a surface of genus 136, in three steps:
  1. an exact hyperbolic model built from 5040 triangles;
  2. three small sign-twisted quotient problems that together see every irreducible representation of $A_7$;
  3. guaranteed lower bounds from Crouzeix–Raviart elements and a Cholesky factorisation with a-priori error bounds.
- **Byproduct.** The same machinery identifies the first eigenspace exactly, and that representation is the input of Chapter 3.

**Results.**
- Every $(2,4,7)$ curve has $\lambda_1\ge0.34089$, so $\operatorname{gon}(C)\ge24$ and $\operatorname{gon}(D)\ge12$.
- Classes 12, 14 have $\lambda_1\ge0.355696$, so $\operatorname{gon}\ge25$.
- The ten other rigid signatures of genus $\le335$ all have $\operatorname{gon}\ge25$.
- For classes 0, 1 the first eigenspace is located exactly; this is the input of Chapter 3.

## 2.1 Gonality from $\lambda_1$

**Lemma 2.1 (Hersch; Yang–Yau; Li–Yau) [P][L].** For a conformal metric on $X$ and holomorphic $\varphi:X\to\mathbb P^1$ of degree $d$: $\lambda_1(X)\,\mathrm{Area}(X)\le8\pi d$.

*Proof.*
1. Hersch's lemma gives a Möbius $\gamma$ with $\int\gamma\circ\varphi=0$ in $S^2\subset\mathbb R^3$.
2. So each coordinate satisfies $\lambda_1\int x_i^2\le\int|\nabla x_i|^2$.
3. Sum over $i$, using $\sum x_i^2=1$ and $\sum|\nabla x_i|^2\,dA=2\varphi^*dA_{S^2}$, whose integral is $8\pi d$. $\square$

With the hyperbolic metric, $\mathrm{Area}=4\pi(g-1)$. So $\operatorname{gon}(C)\ge\frac12(g-1)\lambda_1$, and $\lambda_1>48/(g-1)$ gives $\operatorname{gon}\ge25$. For $(2,4,7)$ this reads $\operatorname{gon}(C)\ge67.5\lambda_1$ and $\operatorname{gon}(D)\ge33.75\lambda_1$, since a map from $D$ lifts to $C$ with twice the degree.

## 2.2 An exact model

$C$ is tiled by $2|G|$ hyperbolic $(\pi/p,\pi/q,\pi/r)$-triangles $U_g,L_g$. The gluings are $U_g\sim L_g,L_{gb},L_{ga}$, and $G$ acts on the left.
- Identify every tile with one Klein-model triangle, parametrised affinely by the simplex $R$. Mirror tiles then glue by the identity on the sides of $R$.
- The Dirichlet form and the mass become exactly $\int\nabla u\cdot A_R\nabla u$ and $\int u^2w_R$, with
  $$A_R=|\det J|\,J^{-1}\frac{I-xx^T}{\sqrt{1-|x|^2}}J^{-T},\qquad w_R=|\det J|(1-|x|^2)^{-3/2}.$$
- So a uniform subdivision of $R$ is a conforming mesh of $C$ and of all its quotients (`orbifold.py`).

## 2.3 Three small twisted problems

For $K\le G$ and $\varepsilon:K\to\{\pm1\}$, the $\varepsilon$-twisted problem is the Laplacian on $\{f(kx)=\varepsilon(k)f(x)\}$. It lives on $2[G:K]$ tiles with signed gluings.

**Lemma 2.2 [P].** If an eigenspace of $C$ contains $\rho$ with $\langle\rho|_K,\varepsilon\rangle>0$, its eigenvalue occurs in the $\varepsilon$-problem. (The $\varepsilon$-isotypic vectors of $\rho|_K$ are $\varepsilon$-twisted eigenfunctions.)

| problem | $K$ | tiles | $\varepsilon$ | irreducibles seen |
|---|---|---|---|---|
| $Q_0$ | $A_7$ | 2 | trivial | $1$ (use $\mu_2$) |
| $Q_1$ | $3^2{:}4$ | 140 | sign | $6,14_a,14_b,15,21$ |
| $Q_2$ | $C_G((16)(23))\cong C_3\rtimes D_8$ | 210 | sign | $10,\overline{10},15,35$ |

They see every irreducible (`cover.py`, exact). Hence $\lambda_1(C)\ge\min(\mu_2(Q_0),\mu_1(Q_1),\mu_1(Q_2))$.

## 2.4 Guaranteed lower bounds [C]

1. **Comparison.** On each cell, ball arithmetic gives $A_e=c_eP_e\preceq A_R$ and $w_e\ge w_R$. By min–max, the true eigenvalues dominate the piecewise-constant ones.
2. **Crouzeix–Raviart [P][L]** (Liu; Carstensen–Gedicke). With $C_h^2=\max_e\kappa_e^2w_e/\lambda_{\min}(A_e)$ and $\kappa_e^2=\frac{(1/n)^2}8+\frac{(\sqrt2/n)^2}{j_{1,1}^2}$,
   $$\lambda_k\ge\frac{\lambda_{k,h}}{1+C_h^2\lambda_{k,h}} .$$
   The CR interpolation is consistent across tiles, because the gluings are the identity.
3. **Positive definiteness [L].** Since $M_h\succ0$ (it is diagonal for CR elements), $\lambda_{1,h}>\sigma$ iff $B=K_h-\sigma M_h\succ0$.
   - A floating Cholesky of $\hat B-cI$ completes.
   - Higham's bound $|\Delta|\le\gamma_{k+1}|\tilde L||\tilde L^T|$, together with the assembly error $\eta$, gives $B\succeq(c-\gamma_{k+1}\|\tilde L\|_1\|\tilde L\|_\infty-u\max|b_{ii}|-\eta)I\succ0$.
   - For $Q_0$, deflate by a rank-one term.
4. **Counting** (`certify_th.py`). Unpivoted $LDL^T$ gives $P(B-cI)P^T+\Delta=\hat L\hat D\hat L^T$ with $|\Delta|\le\gamma_{k+3}|\hat L||\hat D||\hat L^T|$. Sylvester and Weyl then bound the number of eigenvalues of $B$ below $c-\lVert\Delta\rVert_2$ by $\#\{\hat d_i<0\}$.
5. **Upper bounds.** Rayleigh–Ritz on the reverse comparison ($A_R\preceq c_e^{\max}P_e$, $w_R\ge w_e^{\min}$).

## 2.5 Results [C]

**All four $(2,4,7)$ classes** (`run_certificate.py` → `certificate.txt`; $Q_1$ at $n=96$, $\sigma=0.3409$):

| bound | value |
|---|---|
| $\mu_2(Q_0)$ | $\ge0.99783$ |
| $\mu_1(Q_1)$ | $\ge0.34089$ |
| $\mu_1(Q_2)$ | $\ge0.68902$ (classes 0, 1); $\ge0.76298$ (classes 12, 14) |

Lean (`Hilbert13/SpectralCertificate.lean`) checks three pieces: the abstract CR bound, the Cholesky criterion with permutation, and the final arithmetic.

**Classes 12, 14.** $Q_1$ at $n=128$ gives $\lambda_1\ge0.355696>48/135$ (`certify_th_output.txt`).

**Classes 0, 1: the first eigenspace** (`certify_th.py`).
- The count on $Q_1$ ($n=96$, $\sigma=0.56$) finds one negative pivot, so $\mu_2(Q_1)\ge0.55998$.
- On $C/S_5$, Rayleigh–Ritz gives $\lambda_2\le0.36318$, with $\mathrm{Ind}_{S_5}^G1=1+6+14_a$.
- On $C/A_6$ (functions see $1+6$), $\lambda\ge0.99783$, which excludes the 6.

So **$E_1\cong14_a=14_{(5,2)}$ exactly**, $\lambda_1\in[0.34089,0.36318]$, and all other eigenvalues are $\ge0.55998$.

**The ten other rigid signatures** (`certify_signatures.py` → `certify_signatures_output.txt`; 24 curves, 70 classes):

| signature | $g$ | $48/(g-1)$ | $\lambda_1\ge$ (worst) | Li–Yau $\operatorname{gon}\ge$ |
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

With $(2,4,7)$ these are all faithful $A_7$-curves of genus $\le335$: four branch points force $g\ge421$, and quotient genus $\ge1$ forces $g\ge631$.

## 2.6 Accurate values and validation [N]

- **$(2,4,7)$.** $\lambda_1=0.346267085404$ (classes 0, 1) and $0.359671354895$ (12, 14), with multiplicity 14, by vector-valued Hejhal (Chapter 3).
- **Bolza** (`validate_bolza.py`). The pipeline certifies $3.540,3.701,3.772,3.804$ at $n=8,\dots,64$. These are below the known $3.8388872588$ and converge to it.

## 2.7 Limits

For classes 0, 1, Li–Yau gives only $\lceil67.5\cdot0.34627\rceil=24$.
- A better conformal metric gains at most $0.2\%$ (§6.3).
- The way past 24 uses the *degree*, not the energy: it is a cubic form that vanishes on $E_1$ by representation theory (Chapter 3).
