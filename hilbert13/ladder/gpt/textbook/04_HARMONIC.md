# 4. The cubic degree form and harmonic Hersch

## 4.1 An invariant form carrying degree

For three $\mathbb R^3$-valued fields on an oriented closed surface, define
$$\Theta(u,v,w)=\frac12\sum_{ijk}\epsilon_{ijk}\int u_i\,dv_j\wedge dw_k,
\qquad T(x)=\Theta(x,x,x).$$

**Lemma 4.1 (P).** The form $\Theta$ is symmetric, and for a holomorphic map $x:C\to S^2$ of degree $m$, $T(x)=4\pi m$.

*Proof.* The scalar integral $\int u\,dv\wedge dw$ is alternating in all three entries, by Stokes. Contracting it with the alternating tensor $\epsilon$ makes $\Theta$ symmetric. At equal arguments the integrand is the pulled-back sphere area form. $\square$

**Lemma 4.2 (P+X).** If a $G$-stable function space $E$ satisfies $(\wedge^3E)^G=0$, then $T$ vanishes on $E\otimes\mathbb R^3$.

*Proof.* The scalar integral above is an invariant alternating trilinear form on $E$, hence zero. $\square$

For $E=14_a$, exact character arithmetic gives
$$\dim(\wedge^3E)^G=0,\qquad
\wedge^2E=10\oplus\overline{10}\oplus15\oplus21\oplus35.$$
The power-map character identities used are
$$\chi_{\wedge^2E}(g)=\frac{\chi(g)^2-\chi(g^2)}2,
\quad
\chi_{\wedge^3E}(g)=\frac{\chi(g)^3-3\chi(g)\chi(g^2)+2\chi(g^3)}6.$$
The two-subset permutation character minus the one-subset character identifies the label $14_a$ exactly. The wedge-square has no constants or $14_a$ constituent, so its bracket image lies above the complementary spectral gap.

## 4.2 A trial-space theorem

Let $E_h$ be a nonzero $G$-stable $C^2$ space isomorphic to $14_a$. On this absolutely irreducible space, the energy is a scalar $\lambda_h$ times the mass inner product. Let a balanced holomorphic map be written $x=y+z$ by orthogonal mass projection onto $E_h\otimes\mathbb R^3$. Put
$$t=\|y\|^2,\quad s=\|z\|^2,\quad E_z=E(z),\quad
c=\langle\nabla y,\nabla z\rangle,\quad t+s=A.$$

Assume $a\le\lambda_1,\lambda_h$, the first eigenspace is one copy of $14_a$, and the other mean-zero eigenvalues are at least $b>\lambda_h$. Suppose the trial residual obeys
$$|E(u,v)-\lambda_h\langle u,v\rangle|\le\rho\|u\|\sqrt{E(v)}.$$
For an orthonormal basis $\psi_i$ of $E_h$, set
$$\Lambda_h=\sup_p\lambda_{\max}\sum_i\nabla\psi_i(p)\otimes\nabla\psi_i(p),
\qquad\Gamma_h=\sqrt{A\Lambda_h}.$$
The Poisson bracket map is $\beta(u\wedge v)=\{u,v\}$. On its isotypes define
$$B_h=\sup_{\omega\ne0}
\frac{\langle\beta\omega,(\Delta-\lambda_h)^{-1}\beta\omega\rangle}{\|\omega\|^2}.$$
Only the positive spectral subspaces containing the bracket image enter this definition.

Put
$$\delta=\frac{\rho\sqrt b}{b-\lambda_h},\qquad
\nu=b-(b-a)\delta^2.$$
Require $\delta<1$ and $\nu>\lambda_h$. For a proposed degree bound $m_0$, set $\tau_0=8\pi m_0/A$, $d_0=\tau_0-a>0$, and let $R$ be the positive root of
$$\left(1-\frac a\nu\right)R^2-2\rho R=d_0.$$

**Theorem 4.3 (P).** There is no holomorphic map of degree at most $m_0$ if
$$a\left(1-\frac{R^2}\nu\right)-\rho R
>4\sqrt{B_h/3}\sqrt A\sqrt{d_0+2\rho R}
+\Gamma_h\frac{R^2}{\sqrt\nu}.\tag{4.1}$$

*Proof.* The residual, spectral gap and spectral expansion give $\|P_{E_1}-P_h\|\le\delta$: for high spectral components, $\lambda/(\lambda-\lambda_h)^2\le b/(b-\lambda_h)^2$. Equal finite dimensions convert the one-sided angle estimate to the projection bound. For $z\perp(1\oplus E_h)$, $\|P_{E_1}z\|\le\delta\|z\|$, so $E_z\ge\nu s$.

Balancing gives $E(x)\le A\tau_0$. The residual gives $|c|\le\rho\sqrt{tE_z}$ and $t\ge A-E_z/\nu$. Therefore, with $r=\sqrt{E_z/A}$,
$$a(1-r^2/\nu)-2\rho r+r^2\le\tau_0,$$
and $r\le R$.

Holomorphicity yields $x\times dx=*dx$. Integrating this identity against $dy$ gives
$$2\Theta(x,x,y)=\lambda_ht+c.$$
Since $T(y)=0$, symmetry yields
$$\lambda_ht+c=4\Theta(z,y,y)+2\Theta(z,z,y).\tag{4.2}$$
The three bracket terms satisfy
$$|\Theta(z,y,y)|\le\sqrt{B_h/3}\,t\sqrt{E_z-\lambda_hs}.$$
Indeed resolvent Cauchy--Schwarz applies to the bracket image. Its complementary quadratic form is positive; the $G$-stable orthogonal complement of $E_h$ has lower bound $\nu>\lambda_h$, including the remaining $14_a$ components. For three coefficient vectors $w_i$, the Gram eigenvalues give
$$\sum_{i<j}\|w_i\wedge w_j\|^2\le\frac13\left(\sum_i\|w_i\|^2\right)^2.$$
The gradient frame bound gives $|2\Theta(z,z,y)|\le\sqrt{\Lambda_htsE_z}$.

Finally,
$$E_z-\lambda_hs\le A(\tau_0-\lambda_h)+2\rho\sqrt{tE_z}
\le A(d_0+2\rho R).$$
Insert $t\le A$, $s\le AR^2/\nu$ and $\lambda_ht+c\ge aA(1-R^2/\nu)-\rho AR$ into (4.2), and divide by $A$. The resulting nonstrict inequality contradicts (4.1). $\square$

## 4.3 Flux bounds

**Lemma 4.4 (Jacobian flux; P).** For $J=\tfrac12(u\nabla^\perp v-v\nabla^\perp u)$,
$$\int w\{u,v\}=-\int\nabla w\cdot J,
\qquad\langle\{u,v\},\Delta^{-1}\{u,v\}\rangle\le\|J\|^2.$$
Thus
$$B_h\le\frac b{b-\lambda_h}\max_\sigma\varphi_\sigma,$$
where $\varphi_\sigma=\|J_\omega\|^2/\|\omega\|^2$ on each real wedge-square isotype.

*Proof.* The first identity is integration by parts. It bounds the bracket's dual Dirichlet norm by $\|J\|$. On eigenvalues $\lambda\ge b$, $(\lambda-\lambda_h)^{-1}\le[b/(b-\lambda_h)]\lambda^{-1}$. Schur's lemma reduces the invariant flux quadratic form to four scalars, since the real decomposition is multiplicity-free. $\square$

A second general method is complementary energy. For $q_r(u,v)=E(u,v)-r\langle u,v\rangle$, an approximate solution $u_h$ to $(\Delta-r)u=f$ with residual $\mathcal R$ satisfies
$$\langle f,(\Delta-r)^{-1}f\rangle
=2\langle f,u_h\rangle-q_r(u_h,u_h)+\|\mathcal R\|_{q_r^{-1}}^2.$$
If $|\mathcal R(v)|\le\eta\sqrt{E(v)}$ and the relevant spectrum is at least $b>r$, this is at most
$$2\langle f,u_h\rangle-q_r(u_h,u_h)+\frac{\eta^2}{1-r/b}.$$
This is an upper bound; an uncorrected Rayleigh--Ritz value has the opposite direction. Equilibrated $H(\operatorname{div})$ fluxes can bound residuals: if $-\operatorname{div}p=\lambda_hu$, the energy-dual residual is bounded by $\|\nabla u-p\|$. A divergence defect adds at most its $L^2$ norm divided by $\sqrt a$. These alternatives were proposed in the historical manuscript and are not needed for the present Jacobian-flux certificate.

## 4.4 Exact equivariance from approximate coefficients

In the Poincaré disk centered at a seven-point, regular radial solutions are
$$R_m(u)=\operatorname{Re}\left[u^m(1-u^2)^\alpha
{}_2F_1(\alpha,m+\tfrac12-it;m+1;u^2)\right],
\quad\alpha=\tfrac12-it,\quad\lambda^*=\tfrac14+t^2.$$
A finite Fourier sum with fixed binary coefficients solves the local differential equation exactly. The numerical Hejhal solve locates useful coefficients, with approximate eigenvalues $0.346267085404$ and $0.359671354895$ for the two geometric types; these approximations are N, not spectral certificates.

Average the local sum over the order-seven rotation and project to the rational two-subset constituent. This gives an exactly equivariant local function $\widetilde F$. Translate it to all orbit centers by $F_{g(0)}(z)=\rho(g)\widetilde F(g^{-1}z)$. Blend with a $C^2$ partition of unity
$$\chi_c=\frac{h(d(z,c))}{\sum_{c'}h(d(z,c'))},\qquad
\Psi=\sum_c\chi_cF_c,$$
using a quintic cutoff equal to one for $d\le1.37$ and zero for $d\ge1.62$. The triangle circumradius is less than $1.36005$, so the denominator is positive. The result is exactly equivariant. A positive certified mass proves that its irreducible section space is nonzero.

Since every translated local function solves the same equation, on the base sector
$$ (\Delta-\lambda^*)\Psi=
\sum_{c\ne0}\left[(\Delta\chi_c)D_c-2\nabla\chi_c\cdot\nabla D_c\right],
\qquad D_c=F_c-\widetilde F.$$
This follows by subtracting $\widetilde F$ and using the zero sum of the cutoff derivatives. The unnormalized cutoff being one does not mean its partition weight is one. Mirror reversal changes the order-four rotation to its inverse; this explains the orientation correspondence between this model and the Klein tiling.

## 4.5 Enclosures and the lower bound 25

The certificate treats coefficient bytes as exact data. A proven minimum center separation certifies the finite orbit-center enumeration. Fourier sampling on disks bounds mismatches with explicit aliasing and tail estimates from larger circles and uniform radial envelopes. Whole polar boxes enclose the local functions, first and second derivatives, cutoff derivatives and exact isotypic projectors. Point samples alone would not control a supremum.

If $\eta=\|(\Delta-\lambda^*)\Psi\|/\|\Psi\|$, then $|\lambda_h-\lambda^*|\le\eta$ and $\rho\le2\eta/\sqrt{\lambda_1}$. The mass normalization is $c_h=\int|\Psi|^2/14$; it also scales the gradient supremum and each flux integral. Ball endpoints must be chosen outwards when taking maxima, minima or displaying bounds.

For $m_0=24$, $\tau_0=16/45$. The saved $192\times128$ certificates give:

| Input or inequality | Class 0 | Class 1 |
|---|---|---|
| Lower bound $a$ | $0.346102$ | $0.345770$ |
| Residual upper bound $\rho$ | $1.9\cdot10^{-4}$ | $5.6\cdot10^{-4}$ |
| Gradient upper bound $\Gamma_h$ | $2.273$ | $2.273$ |
| Bracket upper bound $B_h$ | $4.45\cdot10^{-4}$ | $4.45\cdot10^{-4}$ |
| Left side of (4.1), lower bound | $0.33068$ | $0.32959$ |
| Right side of (4.1), upper bound | $0.27130$ | $0.27947$ |

These strict margins prove $\operatorname{gon}(C)\ge25$ for classes 0 and 1 under the spectral and ball-arithmetic trust base. Together with Chapter 3 this applies to all four classes. Lifting quotient pencils gives $\operatorname{gon}(D)\ge13$. The exact embedding gives the upper bounds 42 and 21.

The historical sufficient boxes for the true eigenspace or for a trial space are optional sufficient conditions. Failure to meet an individual box bound does not matter when the full inequality passes. Excluding degree 25 would require a sharper certified bracket optimization; the reported BFGS improvement factor $0.614$ is numerical and does not establish a lower bound 26.
