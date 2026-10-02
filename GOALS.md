# Goals and targets (compact list, 2026-10-02)

Branch `claude/continue-previous-qfhm7j`, draft PR #1. Details live in `HANDOFF.md` §5 and `hilbert13/ladder/*.md`.

## Done (frozen; do not reopen)

- **D1. The ladder rungs.**
  - $g(C)=136$, $g(D)=64$.
  - $\operatorname{gon}(D)\ge9$, with a structure theorem at degree 9.
  - First summit $\operatorname{gon}(D)\ge10$ and further summit $\operatorname{gon}(C)\ge17$; both are superseded by D2.
- **D2. $\operatorname{gon}(C)\ge25$ on all four $(2,4,7)$ classes**, certified (`HH_CERTIFICATE.md`, `VERIFICATION.md` §7). Hence $\operatorname{gon}(C/\langle\tau\rangle)\ge13$.
- **D3. $\operatorname{gon}\ge25$ for every faithful $A_7$-curve of genus $\le335$**, certified (`VERIFICATION.md` §6).
- **D4. $a(A_7)=60$**, i.e. $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$ (GPT, verified; `ACCESSORY_60.md`). The published bound is $\le6$.
- **D5. $\mu(A_7)=90$**: exact full-monodromy compression degree (ours; `MU90.md`).
- **D6. Certified $\lambda_1$ data.** All 24 curves of the other rigid signatures; $(2,4,7)$ to 12 digits numerically and certified to $2.5\cdot10^{-4}$.

## Open targets, in priority order

1. **T1. Make $\gamma(A_7)\ge25$ unconditional.** Write and check an in-repo proof that $\operatorname{gon}\ge25$ for every faithful $A_7$-curve with
   $g\ge336$. It replaces GPT's cited algebra: the pair bound, Farb–Wolfson's genus bound, the Segre gap and the dependent-third cost.
   Pure algebra, no compute. *Highest value per token.*
2. **T2. Paper-style write-up.** Consolidate D2–D5 into one self-contained document (`PAPER.md`): statements, proofs and certificate
   descriptions, with the literature comparison of `LITERATURE_CHECK.md`. *Needed before anything is shared.*
3. **T3. Independent audit** of `HH_CERTIFICATE.md` and `MU90.md`, by GPT (Round 7 Q1, Round 6 Q1). *External; costs us nothing.*
4. **T4. Housekeeping.** The repairs the GPT audit asked for (e.g. the $S_4$ misnaming of $C_3\rtimes D_8$ in NOTES; NOTES/certify/README
   repairs 1, 3, 4, 5, 6). Small.
5. **T5. (Optional) $\operatorname{gon}\ge26$ on classes 0, 1.** Needs a certified sharp $\kappa$. It does not change $\gamma(A_7)$, because classes 12, 14 stay at 25.
6. **T6. (Long shot) Towers / RD.** No candidate invariant yet. Park it unless a new idea appears.
7. **T7. (Low) Arithmetic.** Torsion of $P_E\in C/L_2(5)$ and $P_S$; fields of definition of the classes.

**Budget note.** With limited usage, do T1 then T2, and send T3 to GPT. Skip T5–T7.
