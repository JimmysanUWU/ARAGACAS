# Chapter 3. Harmonic Hersch: $\operatorname{gon}\ge25$ for the classes 0 and 1

**Theorem 3.1.** A $(2,4,7)$ $A_7$-curve of class 0 or 1 has no holomorphic map of degree $\le24$ to $\mathbb P^1$. Hence
$\operatorname{gon}(C)\ge25$, and $\operatorname{gon}(C/\langle\tau\rangle)\ge13$ for every involution $\tau$.

With Chapter 2, $\operatorname{gon}\ge25$ therefore holds on all four $(2,4,7)$ classes. Li–Yau stops at 24 for classes 0, 1 (§2.7). The extra
input here is the *degree*: it is a cubic form that representation theory kills on the first eigenspace.

## 3.1 The cubic form

For $\mathbb R^3$-valued fields on a surface put
$$\Theta(u,v,w)=\tfrac12\sum_{ijk}\epsilon_{ijk}\int u_i\,dv_j\wedge dw_k,\qquad T(x)=\Theta(x,x,x).$$

**Lemma 3.2 [P].** $\Theta$ is totally symmetric, by Stokes. For a holomorphic $x:X\to S^2$ of degree $m$, $T(x)=4\pi m$.

**Lemma 3.3 (the first eigenspace carries no degree) [P].** Let a finite group $G$ act by orientation-preserving isometries, and let $E$ be a
$G$-stable space of smooth functions with $(\wedge^3E)^G=0$. Then $T\equiv0$ on $E\otimes\mathbb R^3$.

*Proof.* $(u,v,w)\mapsto\int u\,dv\wedge dw$ is an alternating $G$-invariant trilinear form on $E$. $\square$

For $E\cong14_{(5,2)}$ the hypothesis holds: `verify_exact.py` checks $(\wedge^314_{(5,2)})^{A_7}=0$ exactly. Moreover $\wedge^214_{(5,2)}=(10\oplus\overline{10})\oplus15\oplus21\oplus35$ contains
no $1$ and no $14$.

## 3.2 The harmonic Hersch inequality

**Setting.**
- $A=540\pi$, $\tau_0=8\pi\cdot24/A=16/45$.
- $x:C\to S^2$ is holomorphic of degree $\le24$, balanced (Lemma 2.1). So $E(x)\le A\tau_0$ and $\int x=0$.
- $E_h\cong14_{(5,2)}$ is a $G$-stable trial space of $C^2$ functions. By Schur, $E=\lambda_h\|\cdot\|^2$ on $E_h$.
- Split $x=y+z$, with $y$ the projection on $E_h$. Put $t=\|y\|^2$, $s=\|z\|^2$, $E_z=E(z)$ and $c=\langle\nabla y,\nabla z\rangle$.

**Inputs** (Chapter 2, classes 0, 1): $E_1\cong14_{(5,2)}$; $\lambda_1\ge0.340893$; every eigenvalue on $(1\oplus E_1)^\perp$ is $\ge b=0.5599822$.

**Constants of the trial space.**
- **residual** $\rho$: $|\langle\nabla u,\nabla v\rangle-\lambda_h\langle u,v\rangle|\le\rho\|u\|\|\nabla v\|$ for $u\in E_h$;
- **gradient** $\Gamma_h=\sqrt{A\Lambda_h}$, with $\Lambda_h=\sup_p\lambda_{\max}\sum_i\nabla\psi_i\otimes\nabla\psi_i$ over an orthonormal basis $(\psi_i)$ of $E_h$;
- **bracket** $B_h=\sup_\omega\langle\beta(\omega),(\Delta-\lambda_h)^{-1}\beta(\omega)\rangle/\|\omega\|^2$, where $\beta(u\wedge v)=\{u,v\}$ is the Poisson bracket.

**Theorem 3.4 (GPT, Round 5, Thm 5.3; re-derived) [P].** Put $\delta=\rho\sqrt b/(b-\lambda_h)$, $\nu=b-(b-a)\delta^2$ and $d=\tau_0-a$, where $a\le\lambda_1$.
Let $R$ be the positive root of $(1-a/\nu)R^2-2\rho R=d$. If
$$a\Big(1-\frac{R^2}\nu\Big)-\rho R\ >\ 4\sqrt{B_h/3}\,\sqrt A\,\sqrt{d+2\rho R}+\Gamma_h\frac{R^2}{\sqrt\nu},\tag{3.1}$$
then no holomorphic map of degree $\le24$ exists.

*Proof.*
1. **Spectral position.** $\|P_{E_1}-P_h\|\le\delta$, and $E(w)\ge\nu\|w\|^2$ for $w\perp1\oplus E_h$. Also $|\lambda_1-\lambda_h|\le\rho\sqrt{\lambda_h/(1-\delta^2)}$.
2. **Energy.** $\lambda_ht+2c+E_z\le A\tau_0$, $|c|\le\rho\sqrt{tE_z}$ and $s\le E_z/\nu$ give $\sqrt{E_z/A}\le R$.
3. **Harmonicity.** For holomorphic $x$, $x\times dx={*dx}$, so $2\Theta(x,x,y)=\langle\nabla x,\nabla y\rangle=\lambda_ht+c$.
   By Lemmas 3.2–3.3, $T(y)=0$, so $\lambda_ht+c=4\Theta(z,y,y)+2\Theta(z,z,y)$.
4. **Bracket term.** $\Theta(z,y,y)=\langle z,N\rangle$ with $N_i=\{y_j,y_k\}$, which lies in isotypes of spectrum $\ge b$.
   Cauchy–Schwarz in $E-\lambda_h\|\cdot\|^2$, together with $\sum_{\rm cyc}\|w_j\wedge w_k\|^2\le t^2/3$, gives $|\Theta(z,y,y)|\le\sqrt{B_h/3}\,t\sqrt{E_z-\lambda_hs}$.
5. **Gradient term.** $|2\Theta(z,z,y)|\le\sqrt{\Lambda_h\,t\,s\,E_z}$.
6. **Combine.** Insert $E_z-\lambda_hs\le A(d+2\rho R)$, $t\le A$ and $t/A\ge1-R^2/\nu$. $\square$

**Lemma 3.5 (Jacobian flux) [P].** For $C^1$ functions $u,v$ put $J=\tfrac12(u\nabla^\perp v-v\nabla^\perp u)$. Then $\int w\{u,v\}=-\int\nabla w\cdot J$, so
$\langle\{u,v\},\Delta^{-1}\{u,v\}\rangle\le\|J\|^2$. Hence, on isotypes whose spectrum is $\ge b$,
$$B_h\ \le\ \frac b{b-\lambda_h}\,\max_\sigma\varphi_\sigma,\qquad\varphi_\sigma=\|J_\omega\|^2/\|\omega\|^2\ \text{ on the isotype }\sigma\subset\wedge^2E_h .$$

By Schur, $\varphi_\sigma$ is one scalar per isotype, an integral over a fundamental domain. No PDE has to be solved.

## 3.3 An exact trial space from near-exact eigenfunctions

**Vector-valued Hejhal expansion** (`hejhal_solve.py`).
- *Model.* Work in the Poincaré disk, with the order-7 point at $0$. Let $X$, $Y$, $C$ be the rotations by $\pi$, $\pi/2$, $2\pi/7$ about the vertices; they satisfy $XYC=1$.
- *The curve.* It is $\mathbb H/\ker\varphi$ with $\varphi(X,Y,C)=(a,b,c)$. This is the mirror image of the tiling curve of Chapter 2, because $\sigma_{BC}\sigma_{AB}=Y^{-1}$;
  spectra, isotypes and gonality are mirror-invariant.
- *Forms.* A $14_{(5,2)}$-eigenfunction is a map $F:\mathbb H\to V_2$ with $F(\delta z)=P_{\varphi(\delta)}F(z)$. Here $V_2\cong14_{(5,2)}$ sits in the 2-subset permutation module $\mathbb R^{21}$.
- *Expansion.* About $0$,
  $$F_0(z)=\sum_{m=0}^{M}R_m(|z|)(A_m\cos m\theta+B_m\sin m\theta),\qquad R_m(u)=\operatorname{Re}\big[u^m(1-u^2)^a{}_2F_1(a,m+\tfrac12-it;m+1;u^2)\big],$$
  with $a=\frac12-it$ and $\lambda^*=\frac14+t^2$. This is an **exact** eigenfunction on the whole disk.
- *Solving.* Imposing $F(Xz)=P_aF(z)$ at sample points gives a linear system, singular exactly at eigenvalues. This yields
  $\lambda_1=0.346267085404$ (classes 0, 1) and $0.359671354895$ (classes 12, 14) to 12 digits, at $M=90$.

**Exact equivariance by a partition of unity** (`hh_eval.py`, `hh_certify.py`).
1. *Symmetrise.* $\tilde F=P_{V_2}\frac17\sum_kP_c^{-k}F_0(\zeta^kz)$ is exactly $C$-equivariant, $V_2$-valued, and still an exact eigenfunction.
2. *Glue.* At each centre $g(0)$ put $F_{g(0)}(z)=P_{\varphi(g)}\tilde F(g^{-1}z)$. Blend with $\chi_c=h(d(z,c))/\sum_{c'}h(d(z,c'))$, where $h(d)=S((r_2-d)/w)$ with $S$ the quintic smoothstep ($r_1=1.37>R_{\rm circ}=1.36005$, $w=0.25$), so $h\equiv1$ on the fundamental sector.
   Then $\Psi=\sum_c\chi_cF_c$ is exactly $\Delta$-equivariant. Its components span an exact copy $E_h\cong14_{(5,2)}$.
3. *Residual.* Since each $F_c$ is an exact eigenfunction, on the fundamental sector
   $$(\Delta-\lambda^*)\Psi=\sum_{c\ne0}\big[(\Delta\chi_c)D_c-2\nabla\chi_c\cdot\nabla D_c\big],\qquad D_c=F_c-\tilde F .$$
   The mismatches $D_c$ are themselves exact eigenfunctions. By $C$-equivariance they reduce to $D_X$ and $D_{Y^2}$.

## 3.4 Certification (`hh_certify.py`)

All quantities are enclosed in ball arithmetic (python-flint, 106 bits). The floating Hejhal coefficients are treated as exact data.

1. **Mismatch suprema.**
   - Cover the active region by disks of hyperbolic radius $0.1732$.
   - On each disk, sample $D\circ T_w$ at 64 points of the circle $u=0.144$. The DFT gives its Fourier–Legendre coefficients up to aliasing.
   - Aliasing and tail are bounded by a crude sup on a larger circle and the uniform estimates $g_{\rm lo}\le R_k(u)/u^k\le G_{\rm hi}$, valid for all $k$.
2. **Box pass.** Polar boxes on the sector. On each box:
   - mean-value enclosures of $\tilde F$, $\nabla\tilde F$ and $J$ (second derivatives from the radial ODE);
   - cutoff derivatives, from the centre distance $\pm$ the box circumradius;
   - exact rational isotype projectors $\Pi_\sigma=\frac{d_\sigma}{2520}\sum_g\chi_\sigma(g)\wedge^2P_g$, with idempotence checked in integers.
3. **Assembly.**
   - $\eta=\|(\Delta-\lambda^*)\Psi\|/\|\Psi\|$, and $\lambda_h\in[\lambda^*\pm\eta]$.
   - $\rho\le2\eta/\sqrt{\lambda_1}$.
   - $\Lambda_h\le\sup|\nabla\Psi|^2/c$, a trace bound with $c=\int_C|\Psi|^2/14$.
   - $\varphi_\sigma=\frac{2520}{d_\sigma c^2}\int_{\rm sector}|\Pi_\sigma J_\Psi|^2$.

## 3.5 Results

$192\times128$ grid (21778 boxes), about 6.5 minutes per class. Outputs: `hh_certify_output_cls0.txt` and `hh_certify_output_cls1.txt`.

| | class 0 | class 1 |
|---|---|---|
| $\sup\lvert D_X\rvert$, $\sup\lvert D_{Y^2}\rvert$ | $8.7\cdot10^{-10}$, $8.7\cdot10^{-9}$ | $2.7\cdot10^{-9}$, $2.8\cdot10^{-8}$ |
| $\eta$ | $\le5.5\cdot10^{-5}$ | $\le1.7\cdot10^{-4}$ |
| $\rho$ | $\le1.9\cdot10^{-4}$ | $\le5.6\cdot10^{-4}$ |
| $\lambda_1\ge a$ | $0.346102$ | $0.345770$ |
| $\Gamma_h$ | $\le2.273$ | $\le2.273$ |
| $\varphi_\sigma$ ($10{+}\overline{10},15,21,35$), $\times10^{-4}$ | $1.697,1.623,1.076,0.920$ | same |
| $B_h$ | $\le4.45\cdot10^{-4}$ | $\le4.45\cdot10^{-4}$ |
| left / right side of (3.1) | $\ge0.33068$ / $\le0.27129$ | $\ge0.32959$ / $\le0.27947$ |

**So (3.1) holds, which proves Theorem 3.1.** (3.1) would survive inflating $B_h$ and $\Gamma_h$ together by 1.36 (class 0) or 1.29 (class 1).

As a by-product, $\lambda_1\in[0.34610,0.34633]$ for class 0. The certificate gives the lower end; the upper end is $\lambda_h$.

## 3.6 Trust base and reproduction

**Trust base.**
- Classical: Hersch balancing, Stokes, Schur, min–max, and the geodesic-polar expansion of eigenfunctions of $\mathbb H$.
- The ${}_2F_1$ form of the radial solution.
- The certified inputs of Chapter 2.
- python-flint enclosures, including ${}_2F_1$. These are retried at doubled precision when arb returns no finite ball.

**Reproduce:**
```
python3 hejhal_solve.py 0 90; python3 hejhal_solve.py 1 90          # regenerates the committed .npz
python3 hh_certify.py 0 coef_cls0_M90.npz 192 128 > hh_certify_output_cls0.txt
python3 hh_certify.py 1 coef_cls1_M90.npz 192 128 > hh_certify_output_cls1.txt
```
