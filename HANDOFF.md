# HANDOFF — Hilbert 13 / A7 gonality ladder

Updated at the end of the second session (2026-09-29). **Read this first**, then
`hilbert13/ladder/NOTES.md` (all proofs; section 7 is new) and `hilbert13/ladder/LIT.md`.

## 0. Where the work lives

Repository `JimmysanUWU/ARAGACAS`, branch `claude/continue-previous-qfhm7j`, draft PR
https://github.com/JimmysanUWU/ARAGACAS/pull/1. The first session's 6 commits (originally on
`claude/hopeful-allen-my45u6`, transferred via `aragacas-ladder.bundle`) are the start of this branch.

## 1. Environment setup

```sh
# Python: group theory, character table, FEM, certified spectral gap
pip install numpy scipy sympy python-flint cvxopt pypdf
cd hilbert13/ladder
python3 cover.py              # exact irrep cover check (1 s)
python3 run_certificate.py    # the full spectral certificate (about 3 min, ~3 GB RAM) -> certificate.txt

# Lean 4 + Mathlib (needs release.lean-lang.org, github.com and the Mathlib cache hosts allowed)
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o elan-init.sh
bash elan-init.sh -y --default-toolchain leanprover/lean4:stable
export PATH=$HOME/.elan/bin:$PATH
cd hilbert13 && lake exe cache get && lake build     # the project pins Lean/Mathlib v4.34.1
```

## 2. Files

| path | what |
|---|---|
| `hilbert13/Hilbert13/Superposition.lean` | first mini-project (Lean, no `sorry`): single superpositions $g(\varphi(x)+\psi(y))$ |
| `hilbert13/README.md` | paper proofs for that mini-project |
| `hilbert13/ladder/NOTES.md` | **main notes**: rungs 1–6 and section 7 (the certified spectral gap) |
| `hilbert13/ladder/LIT.md` | literature map |
| `a7.py`, `triples*.py`, `classreps.py`, `chartab.py` | $A_7$ toolkit, the 96 triples with $a=(01)(23)$ fixed, class reps 0, 1, 12, 14, character table |
| `rung1.py`, `quot.py`, `mingenus.py`, `jacobian.py`, `hdecomp.py`, `cusp.py`, `star.py`, `model99.py` | rungs 1–5 (see NOTES) |
| **`orbifold.py`** | exact hyperbolic model of $C$ and of the twisted quotients $K\backslash C$: Klein-model chart, identity gluings, P1/CR finite elements (numerics) |
| **`certify.py`** | rigorous lower bound for one quotient: ball-arithmetic coefficients, CR lower-bound theorem, verified Cholesky |
| **`cover.py`** | exact check that $Q_0,Q_1,Q_2$ see all irreducibles of $A_7$ |
| **`run_certificate.py`, `certificate.txt`** | the full certificate and its output |
| `spectrum.py` | first-session FEM with a flat approximation (superseded by `orbifold.py`) |

## 3. The task (the user's ladder)

$C$: smooth connected complex curve with faithful $A_7$-action, $C/A_7\cong\mathbb P^1$ branched at three points
with inertia orders 2, 4, 7; $\tau$ an involution, $D=C/\langle\tau\rangle$.
1. $g(C)$, the fixed points of an involution, $g(D)$.
2. Involutions $\mu\ne\tau$ commuting with $\tau$ on $D$.
3. $\operatorname{gon}(D)\ge9$; $f\circ v=f$ vs $f\circ v=M\circ f$.
4. The degree-9 case in detail.
5. **First summit:** $\operatorname{gon}(D)\ge10$ for all $C$ and $\tau$.
6. **Further summit:** $\operatorname{gon}(C)\ge17$, with a fresh audit.

Standing instruction: "Aim for QED of the ladder — or your own invention of a mathematical structure /
novel framework." Lean formalization where feasible.

## 4. Results

- Rungs 1–4 [P][C]: as in the first session (NOTES §1–4). $g(C)=136$, $18$ fixed points, $g(D)=64$,
  $H\cong C_2\times S_3$ on $D$, $\operatorname{gon}(D)\ge9$, structure theorem at degree 9.
- §5 [P][C]: by algebra alone the first summit reduces to the arithmetic condition $(\star)$ on $T=C/C(\tau)$.
  The audit curve $D'$ shows that the $H$-topology of $D$ cannot decide it.
- **§7 [P]+[V] — both summits proved.** For every $(2,4,7)$ $A_7$-curve and every involution $\tau$:
  $$\lambda_1(C)\ge0.33335,\qquad\operatorname{gon}(C)\ge23,\qquad\operatorname{gon}(C/\langle\tau\rangle)\ge12 .$$
  The route has four steps:
  1. Hersch/Yang–Yau: $\lambda_1\mathrm{Area}\le8\pi\deg$, so $\operatorname{gon}(C)\ge67.5\lambda_1$ and $\operatorname{gon}(D)\ge33.75\lambda_1$.
  2. $\lambda_1(C)\ge\min(\mu_2(Q_0),\mu_1(Q_1),\mu_1(Q_2))$, where $Q_1=(3^2{:}4,\mathrm{sgn})$ and $Q_2=(S_4,\mathrm{sgn})$ are sign-twisted
     quotients (140 and 210 triangles) that, together with the $A_7$-invariant problem $Q_0$, see every
     irreducible representation.
  3. Each quotient gets a guaranteed Crouzeix–Raviart lower bound (Liu / Carstensen–Gedicke). The
     comparison uses piecewise-constant coefficients bounded in ball arithmetic, in an exact Klein-model
     chart whose gluings are the identity.
  4. A verified sparse Cholesky with a-priori rounding bounds.

  Upper bounds for context: $\operatorname{gon}(C)\le56$, $\operatorname{gon}(D)\le28$.
- Numerics (exact geometry): $\lambda_1(C)\approx0.34627$ (classes 0, 1) and $\approx0.3597$ (classes 12, 14), isotype
  $14_a$.

## 5. Suggested next steps

1. **Independent re-verification** of the certificate: rerun with another sparse Cholesky
   (scikit-sparse / CHOLMOD LDLᵀ, or an interval Cholesky), or on another machine. Cross-check the
   Carstensen–Gedicke–Rim constant against the original 2012 paper. The notes quote it from
   CG 2014, which was read.
2. **$\operatorname{gon}(C)\ge24$**: needs a certified $\lambda_1>0.34074$ for classes 0, 1 (true value 0.34627). The loss is
   the $O(h)$ piecewise-constant comparison (at $n=32$ the certified value is 0.3334). Options are
   $n\approx128$ (3.4M dofs for $Q_1$), or keeping the exact weight $w_R$ in the mass form with interval
   quadrature.
3. **Lean**: formalize the finite pieces: the arithmetic of Riemann–Hurwitz and fixed points, the final
   step "$\lambda_1\ge0.33335\Rightarrow$ gonality bounds", and possibly the abstract CR lower-bound lemma in
   finite-dimensional form.
4. $(\star)$ on $T$ is now of independent interest only.
5. Verify the citation for Lemma 4.0 (Martens 1996). No longer needed for the summits.

## 6. Infrastructure notes

- `WebFetch` is blocked for arxiv.org and some publishers; `curl` works. `pdfminer` breaks on this image
  (cryptography/cffi), so use `pypdf`.
- In `cvxopt`, `cholmod.options['supernodal']=2` forces supernodal $LL^T$. `getfactor` returns the factor of a
  CHOLMOD-chosen permutation; the certificate needs only its row counts and norms.
- Do not run several large FEM jobs in one shell call (OOM).
