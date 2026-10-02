# Goals

Chapter references are to `hilbert13/ladder/` (index: `hilbert13/ladder/README.md`).

## Done

| | result | where |
|---|---|---|
| D1 | $g(C)=136$, $g(D)=64$, $\operatorname{gon}(D)\ge9$ (and $\ge10$ by algebra), degree 9 and the audit curve | Ch. 1, §6.1 |
| D2 | $\operatorname{gon}(C)\ge25$ on all four $(2,4,7)$ classes, so $\operatorname{gon}(C/\tau)\ge13$; certified $\lambda_1$ for 24 further curves | Ch. 2–3 |
| D3 | $\gamma(A_7)\ge25$ for every faithful $A_7$-curve | Ch. 4 |
| D4 | $a(A_7)=60$ (sharp) and $\mu(A_7)=90$; the bounds persist without fixed points (Thm 5.11) | Ch. 5 |
| D5 | Schur-twisted Picard group; the embedding $C\hookrightarrow\mathbb P^5$ of degree 60, cut out by the Laza–Zheng cubic and 15 quartics; an exact $g^1_{42}$, so $\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\tau)\le21$; Klein quadrics; elliptic subcovers $\ge60$ | Ch. 7 |

## Open, in priority order

1. **Independent audit** of Chapters 3, 4, 7 and §§5.4–5.5 by GPT (`QUESTIONS_FOR_GPT.md` §1).
2. **Paper.** Typeset the chapters into one PDF. They are already in textbook order.
3. **The exact gonality**, in $[25,42]$ (and $[13,21]$ for $C/\tau$). Every symmetric mechanism stops at 42. Two candidate new inputs:
   - a lower bound by stable reduction at $p=7$ and graph gonality (§7.10);
   - a non-symmetric pencil on the explicit model.
4. **(Optional) $\operatorname{gon}\ge26$ on classes 0, 1.** This needs a certified sharp $\kappa$ (Ch. 3), and does not change $\gamma(A_7)$.
5. **(Long shot) Towers / RD**: bases whose $H^1$ shares a constituent with $H^1(C)$ (§5.5).
6. **(Low) Arithmetic**: are $P_E$, $P_S$ torsion (§6.5)?

**Budget.** With limited usage, send item 1 to GPT and do item 2. Skip items 3–6 unless a new idea appears.
