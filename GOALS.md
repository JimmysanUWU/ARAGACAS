# Goals

Chapter references are to `hilbert13/ladder/` (index: `hilbert13/ladder/README.md`).

## Done

| | result | where |
|---|---|---|
| D1 | the ladder rungs: $g(C)=136$, $g(D)=64$, $\operatorname{gon}(D)\ge9$, the degree-9 structure theorem, the audit curve | Ch. 1 |
| D2 | $\operatorname{gon}(C)\ge25$ on all four $(2,4,7)$ classes, so $\operatorname{gon}(C/\langle\tau\rangle)\ge13$; certified $\lambda_1$ for 24 further curves | Ch. 2–3 |
| D3 | $\gamma(A_7)\ge25$ for every faithful $A_7$-curve ($\ge23$ by algebra alone) | Ch. 4 |
| D4 | $a(A_7)=60$, i.e. $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$, sharp (GPT; verified); $\mu(A_7)=90$ in exact arithmetic; the bounds persist without fixed points (Thm 5.11) | Ch. 5 |
| D5 | Schur-twisted invariant Picard group; a degree-60 model in $\mathbb P^5$ on an explicit cubic fourfold; $\operatorname{gon}(C)\le42$, $\operatorname{gon}(C/\tau)\le21$, sharp for the $\tau$-pencil; elliptic subcovers $\ge60$ | Ch. 7 |

## Open, in priority order

1. **Independent audit** of Chapters 3, 4, 7 and §5.4, by GPT. The questions are written (`QUESTIONS_FOR_GPT.md` §1).
2. **Paper.** Turn the chapters into one self-contained PDF. They are already in textbook order, so this is mostly typesetting.
3. **(Optional) $\operatorname{gon}\ge26$ on classes 0, 1.** It needs a certified sharp $\kappa$. This does not change $\gamma(A_7)$, because classes 12, 14 stop at 25.
4. **The exact gonality**, now in $[25,42]$ (and $[13,21]$ for $C/\tau$). The $\tau$-pencil is exactly 42 and elliptic subcovers give $\ge60$ (§§7.5, 7.9). A better pencil needs a new mechanism.
5. **(Long shot) Towers / RD.** Theorem 5.11 covers bases with $\mathrm{Hom}_G(\mathrm{Alb},\mathrm{Jac}\,C)=0$; open: a base whose $H^1$ shares a constituent with $H^1(C)$.
6. **(Low) Arithmetic.** Is $P_E$ or $P_S$ torsion (§6.7)?

**Budget.** With limited usage, send item 1 to GPT and do item 2. Skip items 3–6 unless a new idea appears.
