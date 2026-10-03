# 3. Spectral lower bounds

## 3.1 Degree and energy

Use the positive Laplacian $\Delta=-\operatorname{div}\nabla$ and Dirichlet energy $E(u)=\int|\nabla u|^2$. Write $A$ for area and $\lambda_1$ for the first positive eigenvalue.

**Theorem 3.1 (Hersch--Yang--Yau--Li--Yau; P).** Every holomorphic degree-$d$ map from a compact Riemann surface with a conformal metric to $\mathbb P^1$ satisfies
$$\lambda_1A\le8\pi d.$$

*Proof.* Hersch balancing supplies a Möbius transformation after which the three sphere coordinates $x_i$ have mean zero. Their Rayleigh inequalities sum to $\lambda_1\int\sum x_i^2\le\int\sum|\nabla x_i|^2$. The first sum is one, while conformality makes the energy twice the pulled-back spherical area, namely $8\pi d$. $\square$

For curvature $-1$, $A=4\pi(g-1)$, giving $\operatorname{gon}(C)\ge\lambda_1(g-1)/2$. On the genus-136 curves, $A=540\pi$ and the multiplier is $135/2$. A quotient pencil on $D$ lifts with twice the degree, giving the multiplier $135/4$ for that quotient. The strict spectral threshold for excluding degree 24 is $48/(g-1)$; equality would not suffice.

## 3.2 Hyperbolic reference coordinates

The regular cover consists of $2|G|=5040$ triangles with angles $\pi/2,\pi/4,\pi/7$. Label them $U_g,L_g$ and glue $U_g$ to $L_g,L_{gb},L_{ga}$ along the three appropriate sides. The action is left multiplication. Vertex cycles have lengths four, eight and fourteen, so the upstairs metric is smooth.

In the Klein disk, identify each triangle with the same affine reference simplex $R$. For $x=x_A+J(s,t)^T$, the exact energy and mass coefficients are
$$A_R=|\det J|J^{-1}\frac{I-xx^T}{\sqrt{1-|x|^2}}J^{-T},\qquad
w_R=|\det J|(1-|x|^2)^{-3/2}.$$
Adjacent tiles have the same reference coordinates on their common side. Uniform reference subdivision is therefore a genuine mesh of the hyperbolic surface. The method does not replace its geometry by a flat surface.

## 3.3 Representation detection

For $K\le G$ and a sign character $\varepsilon$, restrict to functions satisfying $u(kp)=\varepsilon(k)u(p)$. If an eigenspace contains $\rho$ and $\langle\rho|_K,\varepsilon\rangle>0$, its eigenvalue occurs in this restricted problem.

| Problem | Subgroup | Tiles | Character | Representations detected |
|---|---|---|---|---|
| $Q_0$ | $A_7$ | 2 | trivial | $1$; exclude the constant mode |
| $Q_1$ | $C_3^2\!:\!C_4$ | 140 | sign | $6,14_a,14_b,15,21$ |
| $Q_2$ | $C_G((16)(23))$ | 210 | sign | $10,\overline{10},15,35$ |

The last group has order 24, a normal $C_3$, and Sylow-two subgroup $D_8$; it is $C_3\rtimes D_8$. Its nontrivial center and elements of order six distinguish it from $S_4$. The detected multiplicity of $35$ is two. Exact character sums give dimensions $70$ and $105$, the corresponding subgroup indices, and cover every irreducible. Thus
$$\lambda_1(C)=\min\{\mu_2(Q_0),\mu_1(Q_1),\mu_1(Q_2)\}.$$
Here eigenvalues are listed with multiplicity and $\mu_1(Q_0)=0$.

## 3.4 From a finite matrix to a continuous inequality

On a reference cell, choose constant positive coefficients $A_e\preceq A_R$ and $w_e\ge w_R$, enclosed by ball arithmetic. The comparison problem has smaller energy and larger mass. Let $\Pi u$ be the Crouzeix--Raviart interpolant matching edge means. For the broken comparison energy $a$, constant cell coefficients give
$$a(u-\Pi u,v_h)=0,\qquad a(u,u)=a(\Pi u,\Pi u)+a(u-\Pi u,u-\Pi u).$$
The local interpolation inequality gives
$$\|u-\Pi u\|_b^2\le C_h^2a(u-\Pi u,u-\Pi u),$$
where
$$C_h^2=\max_e\frac{w_e}{\lambda_{\min}(A_e)}\left(\frac1{8n^2}+\frac2{n^2j_{1,1}^2}\right),$$
and $j_{1,1}$ is the first positive zero of $J_1$.

**Lemma 3.2 (P).** If the discrete energy dominates $\sigma\|v_h\|_b^2$, then the comparison continuous energy dominates
$$\frac{\sigma}{1+C_h^2\sigma}\|u\|_b^2.$$

*Proof.* The triangle inequality and the two energy bounds give
$$\|u\|_b\le\sigma^{-1/2}\sqrt{a(\Pi u,\Pi u)}+C_h\sqrt{a(u-\Pi u,u-\Pi u)}
\le\sqrt{\sigma^{-1}+C_h^2}\sqrt{a(u,u)}.$$
Square and rearrange. $\square$

For $Q_0$, use comparison-mean-zero functions and subtract the interpolant's comparison mean. This changes no energy and decreases the error norm. Transfer to the true metric uses weighted variances $\min_c\int(u-c)^2w$: $w_R\le w_e$ orders these variances correctly despite their different mean-zero subspaces. The higher-eigenvalue version follows from the corresponding min--max interpolation theorem; it is an additional classical input to the eigenvalue-count certificate.

The Bessel constant has an exact proof in `certify.py`. The alternating series for $2J_1(x)/x$ has a rational polynomial lower bound on the required interval; positivity follows from exact Bernstein coefficients. It gives $j_{1,1}>3.8317059702075$. The historical weaker bound $j_{1,1}>19/5$ is also sufficient for the conservative coarse certificate.

## 3.5 Positive definiteness and inertia

Let $B=K_h-\sigma M_h$ be the exact discrete matrix and $\widehat B$ its floating approximation, with $\|B-\widehat B\|_2\le\eta$.

**Lemma 3.3 (factor residual; P).** If a permutation $P$ and matrix $L$ satisfy
$$\|LL^T-P(\widehat B-cI)P^T\|_2\le\delta,\qquad c>\delta+\eta,$$
then $B\succeq(c-\delta-\eta)I\succ0$.

*Proof.* $LL^T$ is positive semidefinite; the two perturbations have norm at most $\delta+\eta$. Orthogonal permutation preserves eigenvalues. $\square$

This statement makes no assumption about the algorithm producing $L$. One may certify its residual directly. The current programs instead use a componentwise sparse backward-error estimate of the form
$$|\Delta|\le\gamma_{k+1}|L||L^T|,\qquad
\gamma_j=\frac{ju}{1-ju},\quad u=2^{-53},$$
and $\|\Delta\|_2\le\gamma_{k+1}\|L\|_1\|L\|_\infty$, together with diagonal-shift and assembly errors. Application to the actual sparse update path and the chosen operation count is part of the stated computer-assisted trust base. Successful factorization alone is insufficient.

For inertia, a factorization of a shifted matrix as $L D L^T$ with a bounded perturbation implies, by Sylvester's law and eigenvalue monotonicity, that the exact number of eigenvalues below the threshold is at most the number of negative diagonal pivots. The implemented bound uses $\gamma_{k+3}|L||D||L^T|$. A positive residual margin must include assembly error and the shift actually tested. The trivial problem uses a rank-one deflation in the constant direction.

Conforming Rayleigh--Ritz, with the reverse coefficient comparisons, supplies upper eigenvalue bounds. Lower and upper comparisons play different roles.

## 3.6 Certified inputs and their consequences

The saved certificate supplies the following downward-rounded bounds under this error model:

| Quantity | Bound |
|---|---|
| $\mu_2(Q_0)$ | $0.99783$ |
| $\mu_1(Q_1)$, all four classes | $0.34089$ |
| $\mu_1(Q_2)$, classes 0 and 1 | $0.68902$ |
| $\mu_1(Q_2)$, classes 12 and 14 | $0.76298$ |

Thus $\operatorname{gon}(C)\ge24$ and $\operatorname{gon}(D)\ge12$. For classes 12 and 14, the stronger $Q_1$ certificate gives $\lambda_1\ge0.355695>48/135$ and hence $\operatorname{gon}(C)\ge25$. The displayed value is smaller than the saved raw bound $0.3556959843961548$.

For classes 0 and 1, $Q_1$ has at most one eigenvalue below $0.5599822$. The quotient $C/S_5$ has an upper bound $0.36318$ for its first nonconstant eigenvalue, and
$$\operatorname{Ind}_{S_5}^G1=1+6+14_a.$$
The $C/A_6$ lower bound excludes the six, since its induction character is $1+6$. Coverage by $Q_0,Q_2$ excludes the other constituents. Therefore $E_1$ is exactly one copy of $14_a$, $\lambda_1\in[0.34089,0.36318]$, and the rest of the mean-zero spectrum is at least $0.5599822$. Chapter 4 uses this representation and gap.

The other ten rigid signatures below genus 336 are
$$ (3,3,5),(2,5,7),(3,3,6),(3,4,4),(2,6,7),
(3,3,7),(2,7,7),(3,4,5),(3,4,6),(4,4,4). $$
Exact signature enumeration and the saved sign-quotient certificates exclude degree at most 24 on every corresponding curve. The unneeded window certificates also give:

| Signature | Genus | Worst spectral lower bound | $(g-1)\lambda_1/2$ |
|---|---|---|---|
| $(3,5,5)$ | 337 | $0.24031$ | $>40.37$ |
| $(3,4,7)$ | 346 | $0.23617$ | $>40.73$ |
| $(3,5,6)$ | 379 | $0.14941$ | $>28.23$ |
| $(4,4,5)$ | 379 | $0.26140$ | $>49.40$ |

The later optional failure at $(4,4,7)$ is a failed numerical factorization, not a mathematical negative result. Validation on the Bolza genus-two surface, with known first eigenvalue approximately $3.83888725884$, is a consistency check for geometry and discretization; convergence is not itself a proof for the $A_7$ curves.

## 3.7 Formal scope

`SpectralCertificate.lean` formalizes the abstract interpolation implication, the perturbed-factor positive-definiteness criterion with permutation, and elementary arithmetic transferring a supplied Hersch inequality and lower bound to integer degree bounds 23, 24 and 12. It assumes the analytical and numerical hypotheses. It does not construct the metric, prove Hersch balancing, certify the sparse implementation, or formalize the lower bound 25.
