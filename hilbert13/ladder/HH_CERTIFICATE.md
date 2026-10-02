# $\operatorname{gon}\ge25$ for the $(2,4,7)$ classes 0 and 1: a certified harmonic-Hersch proof

**Theorem.** Let $C$ be a $(2,4,7)$ $A_7$-curve of class 0 or 1 (genus 136). Then $C$ has no holomorphic map of degree
$\le24$ to $\mathbb P^1$. So $\operatorname{gon}(C)\ge25$, and $\operatorname{gon}(C/\langle\tau\rangle)\ge13$ for every involution $\tau$.

**Consequences.**
- With the earlier certificates (`VERIFICATION.md` §§6–7), $\operatorname{gon}(C)\ge25$ now holds for every faithful $A_7$-curve of
  every rigid signature of genus $\le335$. These are the eleven rigid signatures with all their classes.
- Together with GPT's algebra for $g\ge336$ [cited, not re-verified here], this gives $\gamma(A_7)\ge25$: every faithful
  $A_7$-curve has gonality at least 25.
- $\operatorname{gon}(D)\ge13$ for $D=C/\langle\tau\rangle$, on all four $(2,4,7)$ classes and for every involution $\tau$. A degree-12 map
  on $D$ would compose with $C\to D$ to a degree-24 map on $C$. On the user's ladder this replaces "$\operatorname{gon}(D)\ge12$" and
  "$\operatorname{gon}(C)\ge24$".
- The proof also brackets $\lambda_1(C)$ two-sidedly to about $10^{-4}$ (§4). Before, the bracket was $[0.34089,0.36318]$.

**What is new.** The inequality is GPT's harmonic Hersch inequality, Round 5, Theorem 5.3. This document re-derives its
trial-space form (`ROUND5_REVIEW.md`). The certification method is new:

1. **Vector-valued Hejhal expansion** (`hejhal_solve.py`). An eigenfunction in the $14_{(5,2)}$-isotype of $C$ is the same
   thing as a vector-valued Maass form $F:\mathbb H\to V_2$ for the triangle group $\Delta(2,4,7)$, twisted by
   $\rho=14_{(5,2)}\circ\varphi$.
   - The fundamental domain of $\Delta$ is tiny: area $3\pi/14$.
   - One Fourier–Legendre expansion about the order-7 point, $\sum_mR_m(|z|)(A_m\cos m\theta+B_m\sin m\theta)$, is an *exact*
     $\lambda^*$-eigenfunction of the whole disk.
   - Imposing the half-turn relation at sample points pins down $\lambda_1$ to 12 digits:
     $\lambda_1\approx0.346267085404$ (classes 0, 1) and $0.359671354895$ (classes 12, 14).
2. **Exact equivariance by a partition of unity** (`hh_certify.py`). Put a copy of the expansion at every point of the
   $\Delta$-orbit of the order-7 point, transformed by $\rho$, and blend the copies with a $\Delta$-invariant $C^2$ partition
   of unity. This gives a trial function $\Psi$ that is exactly $\Delta$-equivariant. Its components span an exact copy
   $E_h\cong14_{(5,2)}$ in $C^2(C)$.
   - Each copy is an exact eigenfunction. So $(\Delta-\lambda^*)\Psi$ is made only of commutators $[\Delta,\chi_c]$ applied
     to *mismatches* $F_c-F_0$.
   - The mismatches are themselves exact eigenfunctions. They are bounded by circle sampling: Fourier–Legendre
     coefficients with rigorous aliasing and tail bounds.
3. **Jacobian-flux bound for the bracket resolvent** (new; it replaces the complementary-energy solve of Theorem 5.3).
   The Poisson bracket is a divergence: $\int w\{u,v\}\,dA=-\int\nabla w\cdot J$, with $J=\tfrac12(u\nabla^\perp v-v\nabla^\perp u)$.
   Hence $\langle\{u,v\},\Delta^{-1}\{u,v\}\rangle\le\|J\|^2$, and on spectrum $\ge b$
   $$\langle f,(\Delta-\lambda_h)^{-1}f\rangle\le\frac b{b-\lambda_h}\|J_f\|^2 .$$
   So the bracket constant $B_h$ needs only quartic *integrals* of the trial functions, and no PDE solve.
   - Schur reduces $\|J_\omega\|^2$ on $\wedge^2E_h=(10\oplus\overline{10})\oplus15\oplus21\oplus35$ to four scalars.
   - Each scalar is an integral over one fundamental domain of an $A_7$-invariant density, computed with exact rational
     isotype projectors.

Every number below is a ball-arithmetic enclosure (python-flint, 106 bits). Floating point is used only to choose the
expansion coefficients and grids, which then serve as exact data.

## 1. The inequality (Theorem 5.3 form, re-derived)

**Setting and notation.**
- $A=540\pi$ and $\tau=8\pi\cdot24/A=16/45$.
- $x:C\to S^2$ is holomorphic of degree $m\le24$, balanced by a Möbius map. So $E(x)=8\pi m\le A\tau$ and $\int x=0$.
- $E_h\cong14_{(5,2)}$ is the trial space. Schur gives $E=\lambda_h\|\cdot\|^2$ on $E_h$.
- Write $x=y+z$ with $y=P_hx$, and set $t=\|y\|^2$, $s=\|z\|^2$, $E_z=E(z)$, $c=\langle\nabla y,\nabla z\rangle$.

**Inputs (certified elsewhere for classes 0, 1).**
- $E_1$ is one copy of $14_{(5,2)}$.
- $\lambda_1\ge0.340893$ (`certificate.txt`).
- Every eigenvalue on $(1\oplus E_1)^\perp$ is $\ge b=0.5599822$ (`certify_th_output.txt`).

**Steps.**
1. **Residual.** $\rho$ bounds the energy-dual residual: $|\langle\nabla u,\nabla v\rangle-\lambda_h\langle u,v\rangle|\le\rho\|u\|\|\nabla v\|$ for $u\in E_h$.
   - Then $\|P_{E_1}-P_h\|\le\delta:=\rho\sqrt b/(b-\lambda_h)$, and $E(w)\ge\nu\|w\|^2$ for $w\perp1\oplus E_h$, where $\nu=b-(b-a)\delta^2$.
   - Also $|\lambda_1-\lambda_h|\le\rho\sqrt{\lambda_h/(1-\delta^2)}$. This gives the lower bound $a\le\lambda_1$.
2. **Energy.** $\lambda_ht+2c+E_z\le A\tau$ and $|c|\le\rho\sqrt{tE_z}$, with $s\le E_z/\nu$. So $R_z:=\sqrt{E_z/A}\le R$, where
   $(1-a/\nu)R^2-2\rho R=d:=\tau-a$.
3. **Harmonicity.** For holomorphic $x$ we have $x\times dx={*dx}$, so $2\Theta(x,x,y)=\langle\nabla x,\nabla y\rangle=\lambda_ht+c$.
   - Total symmetry of $\Theta$ and $T(y)=0$ give $\lambda_ht+c=4\Theta(z,y,y)+2\Theta(z,z,y)$.
   - $T(y)=0$ holds because $(\wedge^314_{(5,2)})^{A_7}=0$, checked exactly in `verify_exact.py`.
4. **Bracket term.** $\Theta(z,y,y)=\langle z,N\rangle$, where $N_i=\{y_j,y_k\}$ lies in $10,\overline{10},15,21,35$ (spectrum $\ge b$).
   - Cauchy–Schwarz in $E-\lambda_h\|\cdot\|^2$, and $\sum_{\rm cyc}\|w_j\wedge w_k\|^2\le t^2/3$, give $|\Theta(z,y,y)|\le\sqrt{B_h/3}\,t\sqrt{E_z-\lambda_hs}$.
   - Here $B_h=\frac b{b-\lambda_h}\max_\sigma\varphi_\sigma$, where $\varphi_\sigma$ is the value of $\|J_\omega\|^2/\|\omega\|^2$ on the isotype $\sigma$.
5. **Gradient term.** $|2\Theta(z,z,y)|\le\sqrt{\Lambda_ht\,sE_z}$, where $\Lambda_h=\sup_p\lambda_{\max}\sum_i\nabla\psi_i\otimes\nabla\psi_i$ and $\psi_i$ is an orthonormal basis of $E_h$.
6. **Combine.** Using $E_z-\lambda_hs\le A(d+2\rho R)$, $t\le A$, $t/A\ge1-R^2/\nu$ and $\Gamma_h=\sqrt{A\Lambda_h}$:
   $$a\Big(1-\frac{R^2}\nu\Big)-\rho R\ \le\ 4\sqrt{B_h/3}\,\sqrt A\,\sqrt{d+2\rho R}+\Gamma_h\frac{R^2}{\sqrt\nu}.\tag{5.3}$$
   If (5.3) fails, no such $x$ exists.

## 2. The trial space

- **Model.**
  - Poincaré disk, with the order-7 point at $0$, $P_2>0$ on the real axis, and $P_4$ at angle $\pi/7$.
  - $X$ is the half-turn about $P_2$, $Y$ the quarter-turn about $P_4$, and $C$ the rotation by $2\pi/7$. They satisfy $XYC=1$.
  - For the triple $(a,b,c)$ (`triples_data.py`), the curve is $\mathbb H/\ker\varphi$ with $\varphi(X,Y,C)=(a,b,c)$.
  - **This is the mirror image of the tiling curve** of `orbifold.py` for the same triple. We check that
    $\sigma_{BC}\sigma_{AB}=Y^{-1}$, so the tiling has $\varphi_T(Y)=b^{-1}$, and $\varphi_T=\varphi\circ\mathrm{Ad}_{\bar z}$ up to conjugation by $a$.
  - The spectrum, the isotypes, the certified inputs and the gonality are therefore the same for both curves.
- **Expansion.** $\tilde F=P_{V_2}\frac17\sum_kP_c^{-k}F_0(\zeta^kz)$, where $F_0$ has the floating coefficients of `coef_cls*_M90.npz` ($M=90$).
  - $\tilde F$ is an exact $\lambda^*$-eigenfunction, $V_2$-valued, and exactly $C$-equivariant.
  - $V_2\subset\mathbb R^{21}$ is the $14_{(5,2)}$ inside the 2-subset permutation module. $P_{V_2}$ is exact and rational, and the $P_g$ are permutation matrices.
- **Gluing.** Centres are the points $g(0)$; there are 14 within $r_2+R_{\rm circ}$, at distances $2.1408$ and $2.7201$, found by BFS over star adjacency.
  - Put $F_{g(0)}(z)=P_{\varphi(g)}\tilde F(g^{-1}z)$, $h(d)=S((r_2-d)/w)$ with $S$ the quintic smoothstep ($r_1=1.37>R_{\rm circ}=1.36005$, $w=0.25$), and $H=\sum_ch(d(z,c))\ge1$.
  - Then $\Psi=\sum_c(h_c/H)F_c$ is exactly $\Delta$-equivariant and $C^2$, so its components span $E_h\cong14_{(5,2)}$.
- **Residual.** On the fundamental sector, $(\Delta-\lambda^*)\Psi=\sum_{c\ne0}(\Delta\chi_c)D_c-2\nabla\chi_c\cdot\nabla D_c$, where
  $D_c=F_c-\tilde F$. By $C$-equivariance, $D_c$ reduces to $D_X=P_a\tilde F\circ X-\tilde F$ and $D_{Y^2}=P_{b^2}\tilde F\circ Y^2-\tilde F$ on the
  part of the star where the cutoff is active.

## 3. Certification (`hh_certify.py`)

1. **Mismatch bounds.** Cover the active region by disks of hyperbolic radius $0.1732$.
   - On each disk, $D$ is an exact eigenfunction, so $D\circ T_w=\sum d_mR_{|m|}(u)e^{im\theta}$.
   - Sample 64 points on the circle $u_c=0.144$; this gives the DFT coefficients, which equal $d_m$ up to aliasing.
   - A crude bound $M_2$ on the circle $u_2=0.24$ bounds the aliasing and the tail. It uses uniform estimates
     $g_{\rm lo}\le R_k(u)/u^k\le G_{\rm hi}$, valid for all $k$. These come from $|(b_k)_n|\le(k+1)_n$ and $|(a)_n|\le(|a|)_n$.
   - The derivatives are bounded the same way.
2. **Box pass.** Polar boxes on the sector $|\theta|\le\pi/7$, $u\le r_S(\theta)$. On each box:
   - mean-value enclosures of $\tilde F$, $\nabla\tilde F$ and $J$, using second derivatives from the radial ODE;
   - distances to the centres, as the centre value $\pm$ the circumradius, which give $\chi_c$, $|\nabla\chi_c|$, $|\Delta\chi_c|$;
   - the residual, the normalisation (boxes inside the sector only), $\sup|\nabla\Psi|$, and the isotype densities
     $|\Pi_\sigma J|^2$. The projectors are exact rational $\Pi_\sigma=\frac{d_\sigma}{2520}\sum_g\chi_\sigma(g)\wedge^2P_g$; idempotence is checked in integers.
3. **Assembly.**
   - $\eta_s=\|(\Delta-\lambda^*)\Psi\|/\|\Psi\|$, $\lambda_h\in[\lambda^*-\eta_s,\lambda^*+\eta_s]$, $\rho\le2\eta_s/\sqrt{\lambda_1}$.
   - $\Lambda_h\le\sup|\nabla\Psi|^2_{\rm hyp}/c$ (trace bound), where $c=\int_C|\Psi|^2/14$.
   - $\varphi_\sigma=\frac{2520}{d_\sigma c^2}\int_{\rm sector}|\Pi_\sigma J_\Psi|^2$.

## 4. Results

Run on a $192\times128$ polar grid (21778 boxes), about 6.5 minutes per class. Outputs: `hh_certify_output_cls0.txt`, `hh_certify_output_cls1.txt`.

| | class 0 | class 1 |
|---|---|---|
| $\lambda^*$ (exact binary float) | $0.34626708540413$ | $0.34626708539950$ |
| $\sup\lvert D_X\rvert$, $\sup\lvert\nabla D_X\rvert$ | $8.7\cdot10^{-10}$, $3.3\cdot10^{-8}$ | $2.7\cdot10^{-9}$, $1.0\cdot10^{-7}$ |
| $\sup\lvert D_{Y^2}\rvert$, $\sup\lvert\nabla D_{Y^2}\rvert$ | $8.7\cdot10^{-9}$, $3.0\cdot10^{-7}$ | $2.8\cdot10^{-8}$, $9.6\cdot10^{-7}$ |
| $\eta_s=\lVert(\Delta-\lambda^*)\Psi\rVert/\lVert\Psi\rVert$ | $\le5.46\cdot10^{-5}$ | $\le1.65\cdot10^{-4}$ |
| $\rho$, $\delta$ | $\le1.87\cdot10^{-4}$, $\le6.6\cdot10^{-4}$ | $\le5.64\cdot10^{-4}$, $\le1.98\cdot10^{-3}$ |
| $\lambda_h$ | $[0.346212,0.346322]$ | $[0.346102,0.346432]$ |
| $\lambda_1\ge a$ | $0.3461025$ | $0.3457705$ |
| $\Lambda_h$, $\Gamma_h=\sqrt{A\Lambda_h}$ | $\le3.045\cdot10^{-3}$, $\le2.2729$ | $\le3.045\cdot10^{-3}$, $\le2.2729$ |
| $\varphi_\sigma$: $10{+}\overline{10}$, $15$, $21$, $35$ | $1.697,\ 1.623,\ 1.076,\ 0.920$ ($\times10^{-4}$) | same to 4 digits |
| $B_h=\frac b{b-\lambda_h}\max\varphi_\sigma$ | $\le4.448\cdot10^{-4}$ | $\le4.450\cdot10^{-4}$ |
| $R$ | $\le0.157811$ | $\le0.161418$ |
| left side of (5.3) | $\ge0.330681$ | $\ge0.329591$ |
| right side of (5.3) | $\le0.271290$ | $\le0.279466$ |
| **degree $\le24$** | **excluded** | **excluded** |

**Margin.** $B_h$ and $\Gamma_h$ could be inflated together by $1.36$ (class 0) and $1.29$ (class 1) before (5.3) would hold.

**Sanity checks against the earlier numerics.**
- The floating $\Gamma=\sqrt{A\Lambda}$ of the eigenfunction is $1.578$. The certificate uses the trace bound, which is about $\sqrt2$ larger.
- The P1-FEM values were $1.6713$ ($n=4$) and $1.6346$ ($n=8$), after the $L^TSL$ fix; they decrease toward $1.58$.
- The floating $\varphi_\sigma$ are $1.506,1.430,0.947,0.816$ ($\times10^{-4}$). The enclosures exceed them by about 13%, which is the box discretisation.
- The floating $\lVert\beta\rVert^2$ per isotype agrees with the FEM table in `FRAMEWORK_TOPOLOGICAL_HERSCH.md` §4 to 1–2%.
- Classes 0 and 1 give identical frame quantities, as they must if they are mirror images.

## 5. Trust base

- **Classical.** Hersch balancing; Stokes; Schur's lemma; Courant–Fischer.
- **Spectral.** The spectral theorem for the Laplacian on the isotypes; the expansion of an eigenfunction of $\mathbb H$ in
  geodesic polar coordinates.
- **Special function.** The ${}_2F_1$ form of the regular radial solution.
- **Certified inputs.** `certificate.txt` and `certify_th_output.txt`, whose trust base is `VERIFICATION.md` §8.
- **Exact algebra.** $(\wedge^314_{(5,2)})^{A_7}=0$ (`verify_exact.py`). The isotype projectors are checked by exact integer idempotence.
- **Software.** python-flint (arb/acb), including its ${}_2F_1$ enclosures. These are retried at doubled precision when arb
  returns no finite enclosure at 106 bits.

**Reproduce:**
```
python3 hejhal_solve.py 0 90; python3 hejhal_solve.py 1 90      # optional: the .npz files are committed
python3 hh_certify.py 0 coef_cls0_M90.npz 192 128 > hh_certify_output_cls0.txt
python3 hh_certify.py 1 coef_cls1_M90.npz 192 128 > hh_certify_output_cls1.txt
```
