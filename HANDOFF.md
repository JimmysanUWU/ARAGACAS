# HANDOFF

Status and open problems. The mathematics: [`GUIDE.md`](GUIDE.md).

The [textbook edition](hilbert13/ladder/gpt/textbook/README.md) is being pushed by mathematical sections. It supplies complete statements and proofs across the source collection and keeps exact, analytical, numerical and historical dependencies explicit. Its final source index and PDF assembly are in progress.

## Status

| chapter | established |
|---|---|
| 1, §6.1 | $g(C)=136$, $g(D)=64$; $\operatorname{gon}(D)\ge9$, and $\ge10$ by algebra |
| 2–4 | $\operatorname{gon}(C)\ge25$ on all four classes, so $\operatorname{gon}(C/\tau)\ge13$; $\gamma(A_7)\ge25$; algebraically, $g\ge266\Rightarrow\operatorname{gon}\ge25$ (the three-pencil theorem, §4.5) |
| 5 | $a(A_7)=60$ (sharp), $\mu(A_7)=90$; for a fixed target, the stable threshold over a smooth projective generically free base with the global Albanese vanishing hypothesis is set by its Amitsur subgroup |
| 7 | the twisted Picard group. The degree-60 six-section map is a closed embedding, and the displayed pencil has degree exactly 42 (quotient degree 21), all [P][X]. No quadrics and Klein contact divisors are exact. The quadratic systems are Jacobian/apolar ideals (§7.12). The proposed quartic ideal and Hessian identification remain [N]. **Every** elliptic subcover has degree $\ge60$ [P][X] (§7.13, Frobenius on $A,S,E_{15},E_{21}$). The kernel $C_7$ is impossible [P][X][C], so a pencil below 42 has kernel $N\in\{1,C_2,C_3,C_4,V_4,S_3\}$, conditional on the quotient rank inputs; the unconditional list is in Cor. 7.23(2) |

## Open, in priority order

1. **Audit.**
   - Done on 4 October (`hilbert13/ladder/AUDIT_4OCT.md`): every proof re-derived and every exact script and certificate rerun, with no errors found.
   - An independent human review of Chapters 3, 4, 7 and §§5.4–5.5 is still the publication bar.
2. **Paper:** typeset the chapters into one PDF.
3. **The exact gonality** in $[25,42]$.
   - $C_7$ is removed (Thm 7.22). Remove $S_3$: only $d=30,36$ remain (Cor. 7.23). This needs $\operatorname{gon}(C/S_3)\ge7$, e.g. $K_{4,2}=0$ on the genus-19 quotients, or a Frobenius argument as for $C_7$.
   - Make Theorem 7.13 exact where it matters: Cor. 7.23(2) shows which quotient bounds ($C_5$, $C_6$, $D_8$, $D_{10}$, $D_{12}$, $5{:}4$) the [N] list relies on.
   - The core cases: $\operatorname{gon}(C/\tau)=21$? Pencils with $N=1$?
   - Stable reduction at $p=7$ with graph gonality (§7.11).
   - The entire branch vanishing sequences and free flex orbit; immersion itself and the exact base divisor are now proved (§7.12).
4. $\operatorname{gon}\ge26$ on classes 0, 1 via a certified sharp $\kappa$ (Ch. 3). On classes 12, 14 the harmonic Hersch route is closed: a $21$-eigenvalue near $0.38$ sits just above $\lambda_1$ (§2.5).
5. Towers and RD: correspondences with $u_Z\ne0$ (§5.5).
6. Are $P_E$, $P_S$ torsion (§6.5)?

**4 October, second session.** This session ran an independent audit and the following new mathematics (`hilbert13/ladder/AUDIT_4OCT.md`; §2.5, §7.13):
- free-symmetry divisibility (Lemma 7.17);
- equivariant Li–Yau (Thm 7.18);
- Frobenius on the Jacobian factors at $p=37,43,67$ (Prop. 7.20);
- every elliptic subcover has degree $\ge60$ (Thm 7.21);
- no kernel $C_7$ (Thm 7.22).

**4 October checkpoint.** `gpt/A7_Normalization_Defects.md` and `gpt/A7_Quadratic_Systems.md` develop the normalization-defect and Jacobian/apolar structures. Their exact scripts and saved outputs remove the old embedding, no-quadric and Proposition 7.7/7.9 numerical dependencies. The remaining audits and open problems above are not claimed closed.

## Rules

- Work on `claude/continue-previous-qfhm7j` (draft PR #1). The 4 October second session worked on `claude/practical-bardeen-b3zki8`, in a draft PR into that branch. GPT's branches `codex/a7-weighted-budget-residuals` (PR #2), `codex/a7-quadric-envelopes-and-residuals` and `codex/rigorous-audit-textbook` (the 3 October audit; reviewed and integrated) are not to be merged, pushed to or commented on without the user's go-ahead.
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
