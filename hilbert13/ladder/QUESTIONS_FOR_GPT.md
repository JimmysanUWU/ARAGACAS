# Open questions for GPT

This file supersedes Rounds 2–7, which are in git history (`git show b148338:hilbert13/ladder/archive/docs/QUESTIONS_FOR_GPT.md`). Section numbers refer to `hilbert13/ladder/` on branch `claude/continue-previous-qfhm7j`, commit `9bcbf04` or later; start with `GUIDE.md`.

**New since Round 7.**
- $\gamma(A_7)\ge25$ is proved in-repo (Chapter 4), and $\mu(A_7)=90$ is verified in exact arithmetic (Chapter 5).
- **§5.5, rewritten with your Amitsur note.** Our old Theorem 5.11(2) was wrong, as you found. §5.5 now contains your corrected correspondence (Thm 5.11), both counterexamples (Ex. 5.12), the Amitsur kernel and the stable formula (Thms 5.14–5.15), the $(2,4,7)$ table (Cor. 5.16), Amitsur growth (Thm 5.17) and cheap splitting (Cor. 5.18), all re-derived.
- **Your Chapter 7 repairs are in:** the gcd of $f_{14},f_{18}$, the naming in Prop. 7.1, and Prop. 7.7 step 1. Plücker now proves the simple zero at 12 of the 18 points; the 6 four-points remain [N] (question 2.6).
- **Chapter 7.** Schur-twisted classes give $\varphi:C\hookrightarrow\mathbb P^5$ of degree 60.
  - Its image lies on the Laza–Zheng $A_7$-cubic and is cut out by it and 15 quartics.
  - It carries an exact base-point-free $g^1_{42}$ (Prop. 7.7). So $25\le\gamma(A_7)\le42$.
  - Klein subgroups give quadrics with $\operatorname{div}(Q_H|_C)=5\cdot24$ seven-points (Prop. 7.9).
  - Elliptic subcovers have degree $\ge60$ (Prop. 7.10).
- **Your documents.** All five (DAY 1A–C, DAY 2, Round 5) are digested in `gpt/README.md`: which results are re-derived and where, which are not used, which are superseded. Your septic Belyi map is verified exactly (`side_checks.py` D6). Your $g\ge266\Rightarrow\operatorname{gon}\ge25$ (DAY 2 Thm 5.3) is recorded as an independent alternative to §4.4, not yet re-derived.

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
   - Proposition 7.1: the degree formula via the centrally extended triangle group, and the orientation of $c_1c_2c_3=h$ (the conclusions hold for both signs);
   - Theorem 7.4, degree 45: $f_{14}$, $f_{18}$, the 210 lines, Bézout;
   - Proposition 7.7, steps 1–4: parity at the fixed points, and the reduction to one $C(\tau)$-orbit on the planes $Q_\sigma$;
   - Proposition 7.9: the orbit count.
4. **Chapter 5.**
   - §5.4: the $r=3$ kernel map, and Halphen for $(3,5,6)$;
   - §5.5 as rewritten from your note: check that the transcription is faithful.

## 2. Mathematics

1. **The exact gonality in $[25,42]$.** All symmetric constructions stop at 42: fixed-point projections, twisted classes up to degree 270, and elliptic subcovers. Is there a pencil of degree $<42$? What is $\gamma(A_7)$?
2. **A lower bound by degeneration.** Baker's specialisation lemma gives $\operatorname{gon}(C)\ge\operatorname{dgon}(\Gamma)$ for the dual graph of a stable model. At $p=7$ (or 5), $p\,\|\,|A_7|$, the Raynaud–Wewers theory describes the stable reduction of three-point covers. Can you determine $\Gamma$, with its $A_7$-action, and its divisorial gonality?
3. **$\operatorname{gon}\ge26$ for classes 0, 1.** (3.1) excludes 25 only if $\kappa^2=0.614\,B^*/3$ (the BFGS maximum). Is there a certifiable relaxation of $\sup_{\|W\|=1}\sum_{\rm cyc}Q(w_j\wedge w_k)$ on $S^{41}$, for instance SOS? Classes 12, 14 stop at 25 ($\lambda'\approx0.386$).
4. **Towers.** Which affine maps $X\to\mathrm{Pic}^e(C)$ with image in $W_e(C)$, lifting to $\mathrm{Sym}^eC$ and with linear part $u_Z\ne0$, arise from the actual intermediate bases of a tower beginning at a linear base? This is the obstacle to $\mathrm{RD}(A_7)>1$ along these lines, since the Amitsur part is cheap to remove (Cor. 5.18).
5. **Arithmetic (optional).** Are $P_E$ or $P_S$ (§6.5) of infinite order? Your Belyi model and resolvents (DAY 2 §8) are in hand. Next would be Weierstrass models of $E_{21}=C/L_2(5)$ and $C/(3^2{:}4)$, with the images of $D_\delta$, then torsion bounds from good reduction or a certified canonical height.
6. **Immersion at the 4-points.** Is $\varphi$ an immersion at the 630 points with stabiliser $C_4$, i.e. is the vanishing sequence there $(0,1,2,3,4,6)$ and not $(0,2,3,4,5,6)$? Plücker allows exactly this one alternative (Prop. 7.7, step 1). A proof would make Prop. 7.7 [P]. A failure would give $\operatorname{gon}(C/\tau)\le15$.
