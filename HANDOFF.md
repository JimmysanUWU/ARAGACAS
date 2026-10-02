# HANDOFF: Hilbert 13 and the $A_7$ gonality ladder

**Read this first, then:**
- `hilbert13/ladder/README.md`: the results, chapter map, scripts and reproduction;
- `GOALS.md`: what is done and what is open.

## Where the work lives

| | |
|---|---|
| repository | `JimmysanUWU/ARAGACAS` |
| branch | `claude/continue-previous-qfhm7j` (draft PR https://github.com/JimmysanUWU/ARAGACAS/pull/1) |
| other open PR | GPT's PR #2. Do not merge it or comment on it without the user's go-ahead. |
| proofs | `hilbert13/ladder/1_CURVE.md` … `6_SIDE_RESULTS.md` |
| Lean | `hilbert13/Hilbert13/Superposition.lean` (single superpositions; see `hilbert13/README.md`) and `SpectralCertificate.lean` |
| history | `hilbert13/ladder/archive/` (superseded notes and scripts) and git history (session logs up to commit `d5d99c4`) |
| web | the root `index.html`, `src/`, `css/`, `dist/` are an unrelated web stub |

## Setup

```sh
pip install numpy scipy sympy python-flint cvxopt pypdf
cd hilbert13/ladder && python3 cover.py          # smoke test, 1 s

# Lean 4 + Mathlib. The hosts release.lean-lang.org, github.com and the Mathlib cache must be allowed.
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o elan-init.sh
bash elan-init.sh -y --default-toolchain leanprover/lean4:stable
export PATH=$HOME/.elan/bin:$PATH
cd hilbert13 && lake exe cache get && lake build   # pinned to Lean/Mathlib v4.34.1
```

## Workflow

- **Collaboration.** The user runs GPT in parallel.
  - GPT's documents arrive as PDFs in `hilbert13/ladder/gpt/`.
  - Each of its claims is re-derived on paper [P] or recomputed [X] before it is used.
  - Questions for GPT go in `hilbert13/ladder/QUESTIONS_FOR_GPT.md`, paste-ready.
- **Status tags** are used throughout; see the ladder README.
- **Commits.**
  - End every commit message with the trailers given in the session's instructions.
  - Never put model identifiers in commits, PRs or files.
  - Push with `git push -u origin claude/continue-previous-qfhm7j`, retrying on network errors.
- **Budget.** The user's usage is limited, so prefer small, decisive steps.

## Infrastructure notes

- **Fetching papers.** `WebFetch` is blocked for arxiv.org and some publishers, but `curl` works. `pdfminer` breaks on this image (cryptography/cffi), so use `pypdf`.
- **Memory.** Large FEM jobs use up to 7 GB. Run one per shell call, or they run out of memory.
- **cvxopt.**
  - `cholmod.options['supernodal']=2` forces supernodal $LL^T$.
  - `getfactor` returns the factor of a CHOLMOD-chosen permutation; the certificates use only its row counts and norms.
- **python-flint.** `acb_hypgeom_2f1` can return a non-finite ball at large $m$. `hh_eval.f21` retries at doubled precision.
