# Open questions for GPT

This file supersedes Rounds 2–7, which are archived in `archive/docs/QUESTIONS_FOR_GPT.md`. Chapter and section numbers refer to this folder
(see `README.md`).

**Since Round 7.**
- $\gamma(A_7)\ge25$ is now proved in-repo (Chapter 4). This answers Round 7 Q2 by our own argument.
- $\mu(A_7)=90$ is now verified in exact arithmetic (Chapter 5, `verify_mu90_exact.py`).
- New Chapter 7: Schur-twisted invariant line bundles give a degree-60 model of $C$ in $\mathbb P^5$ on an explicit $3.A_7$-invariant cubic fourfold,
  and $\operatorname{gon}(C)\le42$. So $25\le\gamma(A_7)\le42$.
- New §5.5 (Theorem 5.11): compression bounds hold over any base with $\mathrm{Hom}_G(\mathrm{Alb},\mathrm{Jac}\,C)=0$, fixed point or not. This answers the fixed-target question
  ($90\mid d$) and the quadratic-accessory obstruction of Round 5 Prop. 4.2.
- Proposition 7.7: the $\tau$-pencil is exactly a base-point-free $g^1_{42}$. Extra base points would be $\tau$-symmetric double points on 8 explicit planes, which $\varphi(C)$ misses.
- Proposition 7.8: $H^1(C,\mathbb Z)$ with cup product is computed from the dessin. Elliptic subcovers of types $15$ and $21$ have degree $\ge60$, with equality for $C/A_5$ and $C/L_2(5)$.

## 1. Audits (most important)

1. **Chapter 3, harmonic Hersch certificate.** Please check:
   - the Jacobian-flux bound, Lemma 3.5;
   - the commutator form of the residual (§3.3, using $h\equiv1$ on the sector);
   - the uniform radial bounds $g_{\rm lo}\le R_k(u)/u^k\le G_{\rm hi}$ for all $k$;
   - the mirror identification $\sigma_{BC}\sigma_{AB}=Y^{-1}$.
2. **Chapter 4, large genus.** Please check:
   - Lemma 4.3 (the orbit size; it needs $A_6\not\subset\mathrm{PGL}_2$);
   - the tangential case of Lemma 4.5 (proximity at the first infinitely near point);
   - that the chain bound $\pi(td,2^t-1)\le\pi(3d,7)$ covers every chain that occurs for $d\le24$.
3. **Chapter 7: Schur-twisted invariant line bundles.** Please check:
   - Proposition 7.1, the degree formula $\deg L=2520(n+\sum r_i/e_i-\rho_0)$ for $\varepsilon$-twisted classes, via $K^w\otimes$(flat) over the universal central extension $\langle c_i\mid c_1^2=c_2^4=c_3^7=c_1c_2c_3\rangle$; in particular the orientation of $c_1c_2c_3=h$ (our conclusions hold for both signs);
   - Theorem 7.4 (degree 45: $f_{14},f_{18}$ and the 210 involution lines; Bézout);
   - Corollary 7.6: the pencil $|E_-(\hat\tau)|\subset\mathbf 6\subseteq H^0(L_{60})$, giving $\operatorname{gon}(C)\le42$.
   - Proposition 7.7, steps 1–4: the parity argument at the fixed points, and the reduction of extra base points to one $C(\tau)$-orbit on the planes $Q_\sigma$.
4. **Chapter 5: §5.4 (the exclusion of 72 and 84) and §5.5.** In particular:
   - the $r=3$ kernel-map step;
   - the use of Halphen for the two $(3,5,6)$ classes;
   - Theorem 5.11, in particular the equivariant splitting of $\mathrm{Pic}(X\times C)$ and the vanishing of the invariant correspondence part.

## 2. Mathematics

1. **$\operatorname{gon}\ge26$.** For classes 0 and 1, (3.1) excludes degree 25 numerically only if $\kappa^2=0.614\,B^*/3$ (the BFGS maximum); the bound $B^*/3$ is too weak.
   Is there a certifiable relaxation of $\sup_{\|W\|=1}\sum_{\rm cyc}Q(w_j\wedge w_k)$ on $S^{41}$, with $Q=\sum_\sigma b_\sigma P_\sigma$? For example:
   - an SOS bound;
   - a bound on the $10{+}\overline{10}$ and $15$ mass of decomposable 2-vectors.

   Classes 12 and 14 stop at 25, because $\lambda'\approx0.386$.
2. **The exact gonality**, now in $[25,42]$. The $\tau$-pencil is sharp at 42 (Prop. 7.7), and elliptic subcovers give $\ge60$ (Prop. 7.8).
   Is there a pencil of degree $<42$ by any other mechanism? What is the true $\gamma(A_7)$?
3. **Towers.** Theorem 5.11 shows that the compression bound survives without fixed points when $\mathrm{Hom}_G(\mathrm{Alb}\,X,\mathrm{Jac}\,C)=0$ for the base. The linearised
   bound $\mu$ persists when $\mathrm{Pic}(B)$ is linearisable (for instance after $t^2=q(v)$); in general the twisted bound $\tilde\mu\ge\operatorname{gon}$ persists.
   What replaces this when $H^1$ of an intermediate total space shares a constituent with $H^1(C)$, so the fibre classes $[Z_b]$ can move in $\mathrm{Pic}(C)$?
   This is now the obstacle to $\mathrm{RD}(A_7)>1$.
4. **Arithmetic (optional).** Are $P_E$ or $P_S$ (§6.7) of infinite order? This needs an explicit model of $E=C/L_2(5)$ through the degree-7 Belyi map
   with passport $[2^21^3,\,4\,2\,1,\,7]$.
