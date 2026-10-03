# Chapter 3. Harmonic Hersch: $\operatorname{gon}\ge25$ for classes 0 and 1

**Theorem 3.1.** A $(2,4,7)$ curve of class 0 or 1 has no map of degree $\le24$ to $\mathbb P^1$. So $\operatorname{gon}(C)\ge25$ and $\operatorname{gon}(C/\tau)\ge13$ on all four classes (with Chapter 2).

**Idea.**
- **Li–Yau stops at 24 here** (§2.7). It uses only the energy of a balanced map. The degree carries more information: it is the cubic form $T(x)=\Theta(x,x,x)$.
- **The cubic form dies on $E_1$.** The first eigenspace is $14_{(5,2)}$, whose $\wedge^3$ has no invariants. So $T$ vanishes identically on fields built from $E_1$ alone.
- **The squeeze.** The spectral gap forces a map of degree $\le24$ to put almost all of its energy in $E_1$, where it carries no degree. The degree must then come from a small remainder, which is impossible.
- **Making it exact.** Harmonicity of the map turns this into an exact identity (harmonic Hersch). An exact equivariant trial space, certified in ball arithmetic, then closes classes 0 and 1.

## 3.1 The cubic form

For $\mathbb R^3$-valued fields put
$$\Theta(u,v,w)=\tfrac12\sum_{ijk}\epsilon_{ijk}\int u_i\,dv_j\wedge dw_k,\qquad T(x)=\Theta(x,x,x).$$

**Lemma 3.2 [P].** $\Theta$ is totally symmetric (Stokes). For holomorphic $x:X\to S^2$ of degree $m$, $T(x)=4\pi m$, since $T(x)$ is the integral of the pulled-back area form of $S^2$.

**Lemma 3.3 (the first eigenspace carries no degree) [P].** Let $G$ act by orientation-preserving isometries, and let $E$ be a $G$-stable space of functions with $(\wedge^3E)^G=0$. Then $T\equiv0$ on $E\otimes\mathbb R^3$, since $(u,v,w)\mapsto\int u\,dv\wedge dw$ is an alternating invariant trilinear form on $E$.

For $E\cong14_{(5,2)}$ the hypothesis holds exactly (`verify_exact.py`). Also, $\wedge^214_{(5,2)}=(10\oplus\overline{10})\oplus15\oplus21\oplus35$ contains no $1$ and no $14$.

## 3.2 The harmonic Hersch inequality

**Setting.**
- $A=540\pi$ and $\tau_0=8\pi\cdot24/A=16/45$.
- $x$ is holomorphic of degree $\le24$, and balanced, so $E(x)\le A\tau_0$ and $\int x=0$.
- $E_h\cong14_{(5,2)}$ is a $G$-stable $C^2$ trial space, with $E=\lambda_h\|\cdot\|^2$ on it.
- Split $x=y+z$, with $y\in E_h$. Put $t=\|y\|^2$, $s=\|z\|^2$, $E_z=E(z)$ and $c=\langle\nabla y,\nabla z\rangle$.

**Inputs** (Chapter 2): $E_1\cong14_{(5,2)}$; $\lambda_1\ge a$; all other eigenvalues are $\ge b=0.5599822$.

**Constants of $E_h$.**
- **residual** $\rho$: $|\langle\nabla u,\nabla v\rangle-\lambda_h\langle u,v\rangle|\le\rho\|u\|\|\nabla v\|$ for $u\in E_h$;
- **gradient** $\Gamma_h=\sqrt{A\Lambda_h}$, with $\Lambda_h=\sup_p\lambda_{\max}\sum_i\nabla\psi_i\otimes\nabla\psi_i$ over an orthonormal basis $(\psi_i)$;
- **bracket** $B_h=\sup_\omega\langle\beta(\omega),(\Delta-\lambda_h)^{-1}\beta(\omega)\rangle/\|\omega\|^2$, where $\beta(u\wedge v)=\{u,v\}$ is the Poisson bracket.

**Theorem 3.4 (GPT Round 5; re-derived) [P].** Put $\delta=\rho\sqrt b/(b-\lambda_h)$, $\nu=b-(b-a)\delta^2$ and $d=\tau_0-a$. Let $R$ be the positive root of $(1-a/\nu)R^2-2\rho R=d$. If
$$a\Big(1-\frac{R^2}\nu\Big)-\rho R\ >\ 4\sqrt{B_h/3}\,\sqrt A\,\sqrt{d+2\rho R}+\Gamma_h\frac{R^2}{\sqrt\nu},\tag{3.1}$$
then no holomorphic map of degree $\le24$ exists.

*Proof.*
1. **Spectral position.** $\|P_{E_1}-P_h\|\le\delta$; $E(w)\ge\nu\|w\|^2$ on $(1\oplus E_h)^\perp$; and $|\lambda_1-\lambda_h|\le\rho\sqrt{\lambda_h/(1-\delta^2)}$.
2. **Energy.** From $\lambda_ht+2c+E_z\le A\tau_0$, $|c|\le\rho\sqrt{tE_z}$ and $s\le E_z/\nu$, we get $\sqrt{E_z/A}\le R$.
3. **Harmonicity.** For holomorphic $x$, $x\times dx={*dx}$, so $2\Theta(x,x,y)=\lambda_ht+c$. Lemma 3.3 gives $T(y)=0$, hence $\lambda_ht+c=4\Theta(z,y,y)+2\Theta(z,z,y)$.
4. **Bracket term.** $\Theta(z,y,y)=\langle z,N\rangle$ with $N_i=\{y_j,y_k\}$, which lies in isotypes of spectrum $\ge b$. Cauchy–Schwarz and $\sum_{\rm cyc}\|w_j\wedge w_k\|^2\le t^2/3$ give $|\Theta(z,y,y)|\le\sqrt{B_h/3}\,t\sqrt{E_z-\lambda_hs}$.
5. **Gradient term.** $|2\Theta(z,z,y)|\le\sqrt{\Lambda_hts E_z}$.
6. **Combine,** with $E_z-\lambda_hs\le A(d+2\rho R)$, $t\le A$ and $t/A\ge1-R^2/\nu$. $\square$

**Lemma 3.5 (Jacobian flux) [P].** For $J=\tfrac12(u\nabla^\perp v-v\nabla^\perp u)$, $\int w\{u,v\}=-\int\nabla w\cdot J$. So $\langle\{u,v\},\Delta^{-1}\{u,v\}\rangle\le\|J\|^2$, and
$$B_h\le\frac b{b-\lambda_h}\max_\sigma\varphi_\sigma,\qquad\varphi_\sigma=\|J_\omega\|^2/\|\omega\|^2\ \text{on the isotype }\sigma\subset\wedge^2E_h .$$
By Schur, each $\varphi_\sigma$ is one scalar per isotype, an integral over a fundamental domain. No PDE has to be solved.

## 3.3 An exact trial space

**Vector-valued Hejhal** (`hejhal_solve.py`).
- Work in the Poincaré disk with the order-7 point at $0$, and let $X,Y,C$ be the rotations by $\pi,\pi/2,2\pi/7$ ($XYC=1$). Then $C=\mathbb H/\ker\varphi$ with $\varphi(X,Y,C)=(a,b,c)$. This is the mirror image of the tiling curve ($\sigma_{BC}\sigma_{AB}=Y^{-1}$); spectra and gonality are mirror-invariant.
- An eigenfunction is $F:\mathbb H\to V_2\cong14_{(5,2)}\subset\mathbb R^{21}$ (2-subsets), with $F(\delta z)=P_{\varphi(\delta)}F(z)$.
- About 0,
  $$F_0=\sum_{m\le M}R_m(|z|)(A_m\cos m\theta+B_m\sin m\theta),\qquad R_m(u)=\operatorname{Re}\big[u^m(1-u^2)^a{}_2F_1(a,m+\tfrac12-it;m+1;u^2)\big],$$
  with $a=\tfrac12-it$ and $\lambda^*=\tfrac14+t^2$. This is an *exact* eigenfunction on the whole disk.
- Imposing $F(Xz)=P_aF(z)$ at sample points gives $\lambda_1=0.346267085404$ (classes 0, 1) and $0.359671354895$ (12, 14), at $M=90$.

**Exact equivariance by a partition of unity** (`hh_eval.py`).
1. $\tilde F=P_{V_2}\frac17\sum_kP_c^{-k}F_0(\zeta^kz)$ is exactly $C$-equivariant, and still an exact eigenfunction.
2. Put $F_{g(0)}(z)=P_{\varphi(g)}\tilde F(g^{-1}z)$, and blend with $\chi_c=h(d(z,c))/\sum_{c'}h(d(z,c'))$. Here $h(d)=S((r_2-d)/w)$, with $S$ the quintic smoothstep, $r_1=1.37>R_{\rm circ}=1.36005$ and $w=0.25$, so $h\equiv1$ on the fundamental sector. Then $\Psi=\sum\chi_cF_c$ is exactly equivariant, and its components span an exact copy $E_h\cong14_{(5,2)}$.
3. On the sector,
   $$(\Delta-\lambda^*)\Psi=\sum_{c\ne0}\big[(\Delta\chi_c)D_c-2\nabla\chi_c\cdot\nabla D_c\big],\qquad D_c=F_c-\tilde F .$$
   The $D_c$ are exact eigenfunctions, and they reduce to $D_X$, $D_{Y^2}$.

## 3.4 Certification [C] (`hh_certify.py`, 106-bit ball arithmetic)

The Hejhal coefficients are treated as exact data.
0. **Centres.** The orbit points of 0 within $r_2+R_{\rm circ}$ are found by a search over adjacent stars. Distinct orbit points are $\ge2d_{72}$ apart (each star contains the disk of radius $d_{72}$), which certifies the enumeration.
1. **Mismatches.**
   - Cover the active region by disks of radius $0.1732$.
   - On each, the 64-point DFT of $D\circ T_w$ on $u=0.144$ gives Fourier–Legendre coefficients.
   - Aliasing and tail are bounded by a sup on a larger circle and the uniform bounds $g_{\rm lo}\le R_k(u)/u^k\le G_{\rm hi}$.
2. **Boxes.** On each polar box of the sector:
   - mean-value enclosures of $\tilde F$, $\nabla\tilde F$ and $J$ (second derivatives from the radial ODE);
   - cutoff derivatives from the centre distance $\pm$ the box radius;
   - exact rational projectors $\Pi_\sigma=\frac{d_\sigma}{2520}\sum\chi_\sigma(g)\wedge^2P_g$.
3. **Assembly.**
   - $\eta=\|(\Delta-\lambda^*)\Psi\|/\|\Psi\|$, with $\lambda_h\in[\lambda^*\pm\eta]$, and $\rho\le2\eta/\sqrt{\lambda_1}$.
   - $\Lambda_h\le\sup|\nabla\Psi|^2/c$, with $c=\int|\Psi|^2/14$.
   - $\varphi_\sigma=\frac{2520}{d_\sigma c^2}\int_{\rm sector}|\Pi_\sigma J_\Psi|^2$.

## 3.5 Results

$192\times128$ grid (21778 boxes), about 6.5 minutes per class (`hh_certify_output_cls{0,1}.txt`).

| | class 0 | class 1 |
|---|---|---|
| $\sup\lvert D_X\rvert$, $\sup\lvert D_{Y^2}\rvert$ | $8.7\cdot10^{-10}$, $8.7\cdot10^{-9}$ | $2.7\cdot10^{-9}$, $2.8\cdot10^{-8}$ |
| $\eta$; $\rho$ | $5.5\cdot10^{-5}$; $1.9\cdot10^{-4}$ | $1.7\cdot10^{-4}$; $5.6\cdot10^{-4}$ |
| $a\le\lambda_1$ | $0.346102$ | $0.345770$ |
| $\Gamma_h$; $B_h$ | $2.273$; $4.45\cdot10^{-4}$ | the same |
| $\varphi_\sigma$ ($10{+}\overline{10},15,21,35$), $\times10^{-4}$ | $1.697,1.623,1.076,0.920$ | the same |
| (3.1): left $\ge$ / right $\le$ | $0.33068$ / $0.27129$ | $0.32959$ / $0.27947$ |

**So (3.1) holds, proving Theorem 3.1.**
- It survives inflating $B_h$ and $\Gamma_h$ by 1.36 (class 0) or 1.29 (class 1).
- By-product: $\lambda_1\in[0.34610,0.34633]$ for class 0.

**Trust base.**
- Hersch, Stokes, Schur, min–max.
- The geodesic-polar ${}_2F_1$ expansion.
- Chapter 2.
- python-flint enclosures, retried at doubled precision when arb returns no finite ball.

**Reproduce.**
```
python3 hejhal_solve.py 0 90; python3 hejhal_solve.py 1 90          # regenerates the .npz
python3 hh_certify.py 0 coef_cls0_M90.npz 192 128 > hh_certify_output_cls0.txt   # likewise class 1
```
