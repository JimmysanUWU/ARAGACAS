# HANDOFF: Hilbert 13 and the $A_7$ gonality ladder

This file holds the live status and the working rules. To learn the mathematics, read [`GUIDE.md`](GUIDE.md) first: reading order, notation, dependencies and verification.

## Where things are

| | |
|---|---|
| repository | `JimmysanUWU/ARAGACAS`, branch `claude/continue-previous-qfhm7j`, draft PR #1 |
| GPT's branches | `codex/a7-weighted-budget-residuals` (PR #2) and `codex/a7-quadric-envelopes-and-residuals`: do not merge, push or comment without the user's go-ahead |
| mathematics | `hilbert13/ladder/1_CURVE.md` … `7_TWISTED.md` (textbook order) |
| Lean | `hilbert13/Hilbert13/Superposition.lean` (`hilbert13/README.md`) and `SpectralCertificate.lean` |
| history | git only: superseded notes and scripts up to `b148338`, session logs up to `d5d99c4` |
| web stub | root `index.html`, `src/`, `css/`, `dist/`, `tsconfig.json`: the original 2022 ARAGACAS stub on `master`, unrelated. The user asked for its removal on this branch; the deletion was blocked by the permission classifier and awaits the user's explicit go-ahead |

## Status

**Done.**
- **Ch. 1, §6.1.** $g(C)=136$, $g(D)=64$, $\operatorname{gon}(D)\ge9$ ($\ge10$ by algebra); degree 9 and the audit curve.
- **Ch. 2–4.** $\operatorname{gon}(C)\ge25$ on all four $(2,4,7)$ classes, so $\operatorname{gon}(C/\tau)\ge13$; also $\gamma(A_7)\ge25$.
- **Ch. 5.** $a(A_7)=60$ (sharp) and $\mu(A_7)=90$; the bounds persist without fixed points (Theorem 5.11).
- **Ch. 7.**
  - Schur-twisted Picard group.
  - $C\hookrightarrow\mathbb P^5$ of degree 60, cut out by the Laza–Zheng cubic and 15 quartics.
  - An exact $g^1_{42}$, so $\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\tau)\le21$.
  - Klein quadrics; elliptic subcovers have degree $\ge60$.

**Open, in priority order.**
1. **Audit** of Ch. 3, 4, 7 and §§5.4–5.5 by GPT (`QUESTIONS_FOR_GPT.md` §1).
2. **Paper.** Typeset the chapters into one PDF.
3. **The exact gonality**, in $[25,42]$. Every symmetric mechanism stops at 42. New inputs would be:
   - stable reduction at $p=7$ with graph gonality (§7.10);
   - a non-symmetric pencil on the explicit model.
4. **(Optional)** $\operatorname{gon}\ge26$ on classes 0, 1, via a certified sharp $\kappa$ (Ch. 3).
5. **(Long shot)** Towers and RD: bases whose $H^1$ shares a constituent with $H^1(C)$ (§5.5).
6. **(Low)** Are $P_E$, $P_S$ torsion (§6.5)?

**Budget.** The user's usage is limited. Send item 1 to GPT and do item 2. Do items 3–6 only on a new idea, in small decisive steps.

## Workflow

- **GPT** runs in parallel. Its PDFs go in `hilbert13/ladder/gpt/`, and each new one is digested into `gpt/README.md`. Every claim is re-derived [P] or recomputed [X] before use. Questions go in `QUESTIONS_FOR_GPT.md`, paste-ready.
- **Commits.**
  - Use the trailers from the session instructions.
  - Put no model identifiers in commits, PRs or files.
  - Push with `git push -u origin claude/continue-previous-qfhm7j`, retrying on network errors.

## Setup

```sh
pip install numpy scipy sympy python-flint cvxopt pypdf
cd hilbert13/ladder && python3 cover.py          # smoke test, 1 s
# Lean 4 + Mathlib (needs release.lean-lang.org, github.com and the Mathlib cache)
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o elan-init.sh
bash elan-init.sh -y --default-toolchain leanprover/lean4:stable && export PATH=$HOME/.elan/bin:$PATH
cd hilbert13 && lake exe cache get && lake build   # Lean/Mathlib v4.34.1
```

## Infrastructure notes

- **Papers.** WebFetch is blocked for arxiv.org; `curl` works. Use `pypdf`, because `pdfminer` breaks on this image.
- **Memory.** FEM jobs use up to 7 GB. Run one per shell call.
- **cvxopt.** `cholmod.options['supernodal']=2` forces supernodal $LL^T$. `getfactor` returns the factor of a CHOLMOD-chosen permutation; the certificates use only its row counts and norms.
- **python-flint.** `acb_hypgeom_2f1` can return a non-finite ball at large $m$; `hh_eval.f21` retries at doubled precision.
- **Killing jobs.** `pkill -f pattern` also kills the calling shell if the pattern is in its own command line. Use `pgrep -f "[p]attern" | xargs -r kill` in a separate call.
