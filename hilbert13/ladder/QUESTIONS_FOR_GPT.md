# Open questions for GPT

This file supersedes Rounds 2–7, which are archived in `archive/docs/QUESTIONS_FOR_GPT.md`. Chapter and section numbers refer to this folder
(see `README.md`).

**Since Round 7.**
- $\gamma(A_7)\ge25$ is now proved in-repo (Chapter 4). This answers Round 7 Q2 by our own argument.
- $\mu(A_7)=90$ is now verified in exact arithmetic (Chapter 5, `verify_mu90_exact.py`).

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
3. **Chapter 5, §5.4: the exclusion of 72 and 84.** In particular:
   - the $r=3$ kernel-map step;
   - the use of Halphen for the two $(3,5,6)$ classes.

## 2. Mathematics

1. **$\operatorname{gon}\ge26$.** For classes 0 and 1, (3.1) excludes degree 25 numerically only if $\kappa^2=0.614\,B^*/3$ (the BFGS maximum); the bound $B^*/3$ is too weak.
   Is there a certifiable relaxation of $\sup_{\|W\|=1}\sum_{\rm cyc}Q(w_j\wedge w_k)$ on $S^{41}$, with $Q=\sum_\sigma b_\sigma P_\sigma$? For example:
   - an SOS bound;
   - a bound on the $10{+}\overline{10}$ and $15$ mass of decomposable 2-vectors.

   Classes 12 and 14 stop at 25, because $\lambda'\approx0.386$.
2. **Upper bounds.** We know $\operatorname{gon}(C)\le56$, from the $(2,4,7)$ curve through $C/P$ with $|P|=8$, $g=12$. Is there a faithful $A_7$-curve of gonality below 56?
   What is the true $\gamma(A_7)$?
3. **Fixed target.** For a fixed $(2,4,7)$ target, does compression with connected full $A_7$-monodromy need $90\mid d$?
4. **Towers.** Prop. 4.2 of Round 5 shows that fixed points do not persist through a quadratic accessory. Is there a weaker invariant that does persist? For example:
   - a fixed point after a bounded further accessory;
   - an obstruction class in $H^1(G,\mathrm{Pic})$.

   This is the obstacle to $\mathrm{RD}(A_7)>1$.
5. **Arithmetic (optional).** Are $P_E$ or $P_S$ (§6.7) of infinite order? This needs an explicit model of $E=C/L_2(5)$ through the degree-7 Belyi map
   with passport $[2^21^3,\,4\,2\,1,\,7]$.
