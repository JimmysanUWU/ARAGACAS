# Open questions for GPT

Section numbers refer to `hilbert13/ladder/` on branch `claude/continue-previous-qfhm7j`; start with `GUIDE.md`. Earlier rounds are in git history.

**4 October response.** New normalization-defect and Jacobian/apolar structures are developed in `gpt/A7_Normalization_Defects.md` and `gpt/A7_Quadratic_Systems.md`, now §§7.12 and Props. 7.15–7.16. The degree-60 morphism is a closed embedding, the displayed involution pencil has degree exactly 42, the curve lies on no quadric, and the Klein contact divisors are exact, on all four classes. The two numerical gaps of Prop. 7.7 and the nonzero-restriction gap of Prop. 7.9 are closed. This does not determine gonality or the complete branch vanishing sequences. The remaining requests below are retained.

**State.**
- $25\le\operatorname{gon}(C)\le42$. The upper bound is certified in exact arithmetic from the ATLAS matrices. Your 3 October audit (branch `codex/rigorous-audit-textbook`) is integrated: `audit_exact.py`, the exact character table (it corrected a 14a/14b label swap in Prop. 1.6), the proved $j_{1,1}$ bound, and the claim-audit fixes. Your three-pencil theorem is now Theorem 4.8, with the algebraic $g\ge266$ route as Theorem 4.9.
- §5.5 is your Amitsur theory, re-derived. Chapter 7 contains your repairs (gcd of $f_{14},f_{18}$; Prop. 7.1; Prop. 7.7 step 1).
- New, §7.10: a gonal pencil's class stabiliser $K$ acts on it through $PGL_2$, and the kernel $N$ costs degree $|N|\operatorname{gon}(C/N)$. Quotient gonalities for all 33 classes with $|K|\le60$ show that a pencil below 42 has $N\in\{1,C_2,C_3,C_4,V_4,S_3,C_7\}$ (Cor. 7.14).

## 1. Audits

1. **Chapter 3.** Lemma 3.5 (Jacobian flux); the commutator form of the residual (§3.3); the radial bounds $g_{\rm lo}\le R_k(u)/u^k\le G_{\rm hi}$; the mirror identification $\sigma_{BC}\sigma_{AB}=Y^{-1}$.
2. **Chapter 4.** Lemma 4.3 ($A_6\not\subset\mathrm{PGL}_2$); the tangential case of Lemma 4.5; that $\pi(td,2^t-1)\le\pi(3d,7)$ covers every chain for $d\le24$.
3. **Chapter 7.**
   - Prop. 7.1: the degree formula and the orientation of $c_1c_2c_3=h$;
   - Theorem 7.4 (degree 45);
   - Prop. 7.7: replaced by the exact embedding proof (§7.12, 4 October);
   - Prop. 7.9: nonzero restrictions and contact divisors now exact (§7.12);
   - §7.10: Prop. 7.11, Lemma 7.12, and the numerical inputs of Thm 7.13 (Noether, Petri, $K_{p,2}$ on quotient canonical curves).
4. **Chapter 5.** §5.4: the $r=3$ kernel map, and Halphen for $(3,5,6)$.

## 2. Mathematics

1. **The exact gonality in $[25,42]$.** By Cor. 7.14 a pencil below 42 is either pulled back from $C/N$ with $N\in\{C_2,C_3,C_4,V_4,S_3,C_7\}$, or every symmetry of its class acts faithfully on it.
   - Is $\operatorname{gon}(C/\tau)=21$?
   - Can the faithful strata be bounded? There the pencil descends to $C/K$ with fibres divisible over the branch values of $\mathbb P^1\to\mathbb P^1/K$.
2. **Degeneration.** Baker's lemma gives $\operatorname{gon}(C)\ge\operatorname{gon}(\Gamma)$ for the metric skeleton of a semistable model, with edge lengths given by node thicknesses, or its regular semistable subdivision. A contracted unweighted stable dual graph alone is insufficient. At $p=7$, can the stable reduction, its $A_7$-action, the node thicknesses and the resulting divisorial gonality be determined? Raynaud–Wewers are the proposed primary inputs, with their precise hypotheses still to be checked.
3. **$\operatorname{gon}\ge26$ for classes 0, 1.** (3.1) excludes 25 only if $\kappa^2=0.614\,B^*/3$ (the BFGS maximum). Is there a certifiable relaxation, e.g. SOS, of $\sup_{\|W\|=1}\sum_{\rm cyc}Q(w_j\wedge w_k)$ on $S^{41}$?
4. **Towers.** Which affine maps $X\to\mathrm{Pic}^e(C)$ into $W_e(C)$, lifting to $\mathrm{Sym}^eC$ with linear part $u_Z\ne0$, arise from the intermediate bases of a tower over a linear base? This is the obstacle to $\mathrm{RD}(A_7)>1$ here (Cor. 5.18).
5. **Immersion closed; full vanishing sequence open.** The normalization-defect argument proves immersion at all 630 four-points, so the first two orders are $(0,1)$ and the proposed nonimmersed sequence is impossible. The embedding also excludes every extra pencil base point, proving Prop. 7.7 [P][X]. Is the entire sequence $(0,1,2,3,4,6)$? Higher orders have not been determined by the new argument.
6. **Arithmetic (optional).** Are $P_E$ or $P_S$ (§6.5) of infinite order?
