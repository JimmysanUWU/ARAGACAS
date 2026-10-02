# Open questions for GPT

This file supersedes Rounds 2–7, which are in git history (`git show b148338:hilbert13/ladder/archive/docs/QUESTIONS_FOR_GPT.md`). Section numbers refer to this folder.

**New since Round 7.**
- $\gamma(A_7)\ge25$ is proved in-repo (Chapter 4), and $\mu(A_7)=90$ is verified in exact arithmetic (Chapter 5).
- **Theorem 5.11.** Compression bounds hold over any base with $\mathrm{Hom}_G(\mathrm{Alb},\mathrm{Jac}\,C)=0$, with or without a fixed point. This answers the fixed-target question ($90\mid d$) and Round 5 Prop. 4.2.
- **Chapter 7.** Schur-twisted classes give $\varphi:C\hookrightarrow\mathbb P^5$ of degree 60.
  - Its image lies on the Laza–Zheng $A_7$-cubic and is cut out by it and 15 quartics.
  - It carries an exact base-point-free $g^1_{42}$ (Prop. 7.7). So $25\le\gamma(A_7)\le42$.
  - Klein subgroups give quadrics with $\operatorname{div}(Q_H|_C)=5\cdot24$ seven-points (Prop. 7.9).
  - Elliptic subcovers have degree $\ge60$ (Prop. 7.10).

## 1. Audits (most important)

1. **Chapter 3, harmonic Hersch.**
   - Lemma 3.5, the Jacobian flux;
   - the commutator form of the residual (§3.3, with $h\equiv1$ on the sector);
   - the uniform radial bounds $g_{\rm lo}\le R_k(u)/u^k\le G_{\rm hi}$;
   - the mirror identification $\sigma_{BC}\sigma_{AB}=Y^{-1}$.
2. **Chapter 4.**
   - Lemma 4.3, which needs $A_6\not\subset\mathrm{PGL}_2$;
   - the tangential case of Lemma 4.5;
   - that $\pi(td,2^t-1)\le\pi(3d,7)$ covers every chain for $d\le24$.
3. **Chapter 7.**
   - Proposition 7.1: the degree formula via the universal central extension, and the orientation of $c_1c_2c_3=h$ (the conclusions hold for both signs);
   - Theorem 7.4, degree 45: $f_{14}$, $f_{18}$, the 210 lines, Bézout;
   - Proposition 7.7, steps 1–4: parity at the fixed points, and the reduction to one $C(\tau)$-orbit on the planes $Q_\sigma$;
   - Proposition 7.9: the orbit count.
4. **Chapter 5.**
   - §5.4: the $r=3$ kernel map, and Halphen for $(3,5,6)$;
   - Theorem 5.11: the equivariant splitting of $\mathrm{Pic}(X\times C)$.

## 2. Mathematics

1. **The exact gonality in $[25,42]$.** All symmetric constructions stop at 42: fixed-point projections, twisted classes up to degree 270, and elliptic subcovers. Is there a pencil of degree $<42$? What is $\gamma(A_7)$?
2. **A lower bound by degeneration.** Baker's specialisation lemma gives $\operatorname{gon}(C)\ge\operatorname{dgon}(\Gamma)$ for the dual graph of a stable model. At $p=7$ (or 5), $p\,\|\,|A_7|$, the Raynaud–Wewers theory describes the stable reduction of three-point covers. Can you determine $\Gamma$, with its $A_7$-action, and its divisorial gonality?
3. **$\operatorname{gon}\ge26$ for classes 0, 1.** (3.1) excludes 25 only if $\kappa^2=0.614\,B^*/3$ (the BFGS maximum). Is there a certifiable relaxation of $\sup_{\|W\|=1}\sum_{\rm cyc}Q(w_j\wedge w_k)$ on $S^{41}$, for instance SOS? Classes 12, 14 stop at 25 ($\lambda'\approx0.386$).
4. **Towers.** What replaces Theorem 5.11 when $H^1$ of an intermediate total space shares a constituent with $H^1(C)$? This is the obstacle to $\mathrm{RD}(A_7)>1$ along these lines.
5. **Arithmetic (optional).** Are $P_E$ or $P_S$ (§6.5) of infinite order? This needs a model of the degree-7 Belyi map with passport $[2^21^3,\,4\,2\,1,\,7]$.
