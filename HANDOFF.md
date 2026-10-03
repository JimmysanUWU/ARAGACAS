# HANDOFF

Status and open problems. The mathematics: [`GUIDE.md`](GUIDE.md).

## Status

| chapter | established |
|---|---|
| 1, §6.1 | $g(C)=136$, $g(D)=64$; $\operatorname{gon}(D)\ge9$, and $\ge10$ by algebra |
| 2–4 | $\operatorname{gon}(C)\ge25$ on all four classes, so $\operatorname{gon}(C/\tau)\ge13$; $\gamma(A_7)\ge25$; algebraically, $g\ge266\Rightarrow\operatorname{gon}\ge25$ (the three-pencil theorem, §4.5) |
| 5 | $a(A_7)=60$ (sharp), $\mu(A_7)=90$; for a fixed target, the stable threshold over a smooth projective generically free base with the global Albanese vanishing hypothesis is set by its Amitsur subgroup |
| 7 | the twisted Picard group. The degree-60 six-section map is a closed embedding, and the displayed pencil has degree exactly 42 (quotient degree 21), all [P][X]. No quadrics and Klein contact divisors are exact. The quadratic systems are Jacobian/apolar ideals (§7.12). The proposed quartic ideal and Hessian identification remain [N]. Pure 15/21 elliptic subcovers $\ge60$. A pencil below 42 has kernel $N\in\{1,C_2,C_3,C_4,V_4,S_3,C_7\}$, conditional on the quotient rank inputs |

## Open, in priority order

1. **Audit** of Ch. 3, 4, 7 and §§5.4–5.5 (`QUESTIONS_FOR_GPT.md` §1).
2. **Paper:** typeset the chapters into one PDF.
3. **The exact gonality** in $[25,42]$.
   - Remove $C_7$ and $S_3$ from Cor. 7.14: $K_{3,2}=0$ on $C/C_7$ and $K_{4,2}=0$ on $C/S_3$ (genus 19; artinian reduction, matrices of size $\sim10^4$).
   - The core cases: $\operatorname{gon}(C/\tau)=21$? Pencils with $N=1$?
   - Stable reduction at $p=7$ with graph gonality (§7.11).
   - The entire branch vanishing sequences and free flex orbit; immersion itself and the exact base divisor are now proved (§7.12).
4. $\operatorname{gon}\ge26$ on classes 0, 1 via a certified sharp $\kappa$ (Ch. 3).
5. Towers and RD: correspondences with $u_Z\ne0$ (§5.5).
6. Are $P_E$, $P_S$ torsion (§6.5)?

**4 October checkpoint.** `gpt/A7_Normalization_Defects.md` and `gpt/A7_Quadratic_Systems.md` develop the normalization-defect and Jacobian/apolar structures. Their exact scripts and saved outputs remove the old embedding, no-quadric and Proposition 7.7/7.9 numerical dependencies. The remaining audits and open problems above are not claimed closed.

## Rules

- Work on `claude/continue-previous-qfhm7j` (draft PR #1). GPT's branches `codex/a7-weighted-budget-residuals` (PR #2), `codex/a7-quadric-envelopes-and-residuals` and `codex/rigorous-audit-textbook` (the 3 October audit; reviewed and integrated) are not to be merged, pushed to or commented on without the user's go-ahead.
- GPT documents go in `hilbert13/ladder/gpt/`, digested in `gpt/README.md`; every claim is re-derived or recomputed before use.
- No model identifiers in files or commits. Usage is limited: prefer small decisive steps.
- Superseded notes and scripts are in git history only (up to `b148338`).

## Setup

```sh
pip install numpy scipy sympy python-flint cvxopt pypdf
cd hilbert13/ladder && python3 cover.py          # smoke test, 1 s
# Lean 4 + Mathlib
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o elan-init.sh
bash elan-init.sh -y --default-toolchain leanprover/lean4:stable && export PATH=$HOME/.elan/bin:$PATH
cd hilbert13 && lake exe cache get && lake build   # Lean/Mathlib v4.34.1
```

Notes: FEM jobs use up to 7 GB (one at a time). `curl` reaches arXiv and the ATLAS; use `pypdf` for PDFs. `hh_eval.f21` retries `acb_hypgeom_2f1` at doubled precision when a ball is non-finite.
