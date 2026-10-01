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
python3 run_certificate.py    # full spectral certificate (about 6 min, 7 GB RAM) -> certificate.txt
python3 run_certificate.py --quick   # coarser (3 min, 3 GB): gon(C) >= 23 -> certificate_quick.txt

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
| **`run_certificate.py`, `certificate.txt`, `certificate_quick.txt`** | the full certificate and its outputs |
| **`hilbert13/Hilbert13/SpectralCertificate.lean`** | Lean, no `sorry`: abstract CR lower bound, perturbed-Cholesky criterion (with permutation), final arithmetic |
| **`validate_bolza.py`, `validation_bolza.txt`** | external check of the whole pipeline on the Bolza surface (known $\lambda_1=3.83888726$) |
| **`REVIEW_GPT.md`, `review_checks.py`** | rigorous review of the Sol/Astra *Ramification transport* synthesis and *Verification ladder* (main theorem verified correct) |
| `spectrum.py` | first-session FEM with a flat approximation (superseded by `orbifold.py`) |
| **`hilbert13/ladder/FRAMEWORK_CONFORMAL.md`** | conformal spectral gonality; §D: results after GPT's audit (D1–D5) |
| **`hilbert13/ladder/QUESTIONS_FOR_GPT.md`** | current, paste-ready questions for GPT |
| **`frontier_checks.py`** | exact checks for §D: $Q_2$ structure, fixed-point Abel–Jacobi carriers, signatures up to genus 529 |
| **`hilbert13/ladder/FRAMEWORK_TOPOLOGICAL_HERSCH.md`** | fourth session: topological Hersch (degree obstruction on $E_1$), certification targets, signature survey, degree lattice |
| **`topo_hersch.py`, `topo_hersch_extras.py`, `topo_hersch_output.txt`** | TH constants on the full curve (FEM), $\wedge^3$ table, sharpened $\kappa$, thresholds; outputs |
| **`signatures_spectrum.py`** | coarse $\lambda_1$ for every faithful $A_7$-curve of a given triangle signature (general $(p,q,r)$ tiling) |
| **`hilbert13/ladder/VERIFICATION.md`** | fifth session: claim-by-claim verification report (proof audit, exact checks, certificates, what is still numerical) |
| **`verify_exact.py`** | exact character computations behind TH ($\wedge^3$, $\wedge^2$, induced modules, quotient multiplicities) |
| **`certify_th.py`, `certify_th_output.txt`** | certified TH eigenvalue inputs: verified-$LDL^T$ count on $Q_1$, $C/S_5$ upper bound, $C/A_6$ lower bound |
| **`certify_signatures.py`, `certify_signatures_output.txt`** | certified $\lambda_1$ for all 24 curves of the ten other rigid signatures (general $(p,q,r)$ in `certify.ref_triangle_arb`) |

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
  $$\lambda_1(C)\ge0.34089,\qquad\operatorname{gon}(C)\ge24,\qquad\operatorname{gon}(C/\langle\tau\rangle)\ge12 .$$
  The route has four steps:
  1. Hersch/Yang–Yau: $\lambda_1\mathrm{Area}\le8\pi\deg$, so $\operatorname{gon}(C)\ge67.5\lambda_1$ and $\operatorname{gon}(D)\ge33.75\lambda_1$.
  2. $\lambda_1(C)\ge\min(\mu_2(Q_0),\mu_1(Q_1),\mu_1(Q_2))$, where $Q_1=(3^2{:}4,\mathrm{sgn})$ and $Q_2=(S_4,\mathrm{sgn})$ are sign-twisted
     quotients (140 and 210 triangles) that, together with the $A_7$-invariant problem $Q_0$, see every
     irreducible representation.
  3. Each quotient gets a guaranteed Crouzeix–Raviart lower bound (Liu / Carstensen–Gedicke). The
     comparison uses piecewise-constant coefficients bounded in ball arithmetic, in an exact Klein-model
     chart whose gluings are the identity.
  4. A verified sparse Cholesky with a-priori rounding bounds.

  The Lean file checks the logic of steps 3–4 and the final arithmetic. For classes 0 and 1, 24 and 12
  are the limits of the Hersch/Yang–Yau route.

  Upper bounds for context: $\operatorname{gon}(C)\le56$, $\operatorname{gon}(D)\le28$.
- Numerics (exact geometry): $\lambda_1(C)\approx0.34627$ (classes 0, 1) and $\approx0.3597$ (classes 12, 14), isotype
  $14_a$.

## 5. Suggested next steps

**Fifth session (2026-10-01): verification** (`VERIFICATION.md`).
- **Certified:** $\operatorname{gon}\ge25$ for all 70 $A_7$-classes of the ten rigid signatures other than $(2,4,7)$, and for $(2,4,7)$
  classes 12, 14 ($Q_1$ at $n=128$: $\lambda_1\ge0.355696$).
- **Classes 0, 1:** the TH theorem is audited [P]. Its eigenvalue inputs are certified: $E_1$ is one copy of $14_{(5,2)}$,
  $\lambda_1\in[0.34089,0.36318]$, and $\lambda'\ge0.55998$, via a new verified $LDL^T$ inertia count.
- **Only $\kappa$ and $\Lambda$ remain numerical**, with 3.5% tolerance. So $\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$ reduces to those two
  constants, plus GPT's $g\ge336$ algebra [cited].
- **Correction.** The fourth-session survey overestimated several $\lambda_1$; the certified table supersedes it.
- **Next.** Certify $\kappa$ and $\Lambda$ (TH doc §5.3, `QUESTIONS_FOR_GPT.md` Round 4, Q2). Optionally sharpen $\lambda_1$ for classes
  0, 1 (0.344 gives 22% tolerance).

**Fourth session (2026-09-30): topological Hersch** (`FRAMEWORK_TOPOLOGICAL_HERSCH.md`).
- **New mechanism [P].** $\deg x$ is the cubic $T(x)=\int x\cdot(x_s\times x_t)$. On $E_1\otimes\mathbb R^3$ it is an invariant alternating 3-form,
  and $(\wedge^3 14_{(5,2)})^{A_7}=0$, so first-eigenfunction maps to $S^2$ have degree 0.
- **The inequality [P].** Expanding $T(y+z)$ gives an inequality (TH) that forces $8\pi m-\lambda_1A$ to be large.
- **Result [N].** For classes 0, 1, (TH) excludes $m=24$ (right side 187 against 281.6 needed). So numerically
  $\operatorname{gon}(C)\ge25$ for all four $(2,4,7)$ curves, and $\operatorname{gon}(D)\ge13$.
- **Conformal route closed [P+N].** The gain is at most $\lambda_1A/\min\bar F$, i.e. $\le0.2\%$ (`FRAMEWORK_CONFORMAL.md` §E).
- **Degree lattice lemma [P]** (TH doc §7): transport is a lattice violation.
- **Survey [N]** (TH doc §6): every other rigid signature passes plain Li–Yau with a wide margin.
- **Next.** Certify $\lambda'\ge0.55$ (second eigenvalue of $Q_1$, with a rank-one deflation for the count). Certify upper
  bounds for $\kappa$ and $\Lambda$ (trial-space route, TH doc §5). Run the $\lambda_1$ certificate on the ten other
  signatures, which needs `certify.py` generalised to the $(p,q,r)$ triangle of `signatures_spectrum.py`.

**Third session (2026-09-30), after GPT's proof-chain audit (pinned 8b31a2ac).** The audit confirms 24/12, verifies
Lemma 5.5, and so gives $\mathrm{ed}_{\mathbb C}(A_7;\le23)>1$ modulo classical inputs. It also asks for six repairs (its page 9);
these are not yet applied. New results are in `FRAMEWORK_CONFORMAL.md` §D:
- the frame function is never constant (via Bryant 1985), so the conformal ascent strictly improves Li–Yau;
- $\delta$ and $(\star)$ are controlled by two points, on the elliptic curve $C/L_2(5)$ and the genus-2 curve $C/(3^2{:}4)$;
- pencil orbits have size at least 35;
- $\operatorname{gon}\ge25$ for all faithful $A_7$-curves with $g\ge336$, so $\mathrm{ed}_{\mathbb C}(A_7;\le29)>1$ reduces to 11 rigid signatures.

0a. **Conformal spectral gonality** (`FRAMEWORK_CONFORMAL.md`, paper-level). Replace the hyperbolic metric by a
   $G$-invariant conformal factor $h$ and certify $\lambda_1(hg)\int h>8\pi\cdot24$. That gives $\operatorname{gon}(C)\ge25$ and $\operatorname{gon}(C/\langle\tau\rangle)\ge13$.
   Test first whether the frame function $\bar F$ of $E_1=14_a$ is non-constant; the ascent direction is $1-\bar F$.
   For classes 12, 14 a finer hyperbolic certificate may already suffice ($67.5\cdot0.3597=24.3$).
0. **Close $\mathrm{ed}_{\mathbb C}(A_7;\le23)>1$.** This needs $\operatorname{gon}\ge19$ for the $(3,3,5)$ curves (genus 169). Either verify
   Lemma 5.5 of the synthesis (Petrakiev ranges), or run the spectral certificate on the $(3,3,5)$ curves; they need
   $\lambda_1>0.2143$. Running it on all eight Table 1 signatures would bypass their §5–6 entirely. First generalize
   `orbifold.reference_triangle` to angle $\pi/p$ at $A$ with $p\ne2$.

1. **Independent re-verification** of the certificate: rerun on another machine, or with a factorization
   that exposes its permutation (scikit-sparse). That would allow an a-posteriori residual check
   $\|L L^T-QBQ^T\|$ computed by our own code. `cvxopt` reorders internally, so the certificate currently
   relies on the a-priori backward error bound. The Carstensen–Gedicke–Rim constant was checked
   against the 2012 paper (Lemma 2.2, with proof).
2. Done this session: $\operatorname{gon}(C)\ge24$ ($n=96$ on $Q_1$). Next integers via this route: $\operatorname{gon}(C)\ge25$ is
   possible only for classes 12, 14. It needs a certified $\lambda_1>0.35556$, which would take $n\approx128$ or
   a second-order comparison (exact weight $w_R$ with interval quadrature).
3. **Lean** (done: `SpectralCertificate.lean`). Possible extensions:
   - the $k$-th eigenvalue version, which needs min–max, not in Mathlib;
   - Hersch's balancing lemma;
   - finite $A_7$ checks by `decide`/`native_decide`, such as the cover multiplicities.
4. $(\star)$ on $T$ is now of independent interest only.
5. Verify the citation for Lemma 4.0 (Martens 1996). No longer needed for the summits.

## 6. Infrastructure notes

- `WebFetch` is blocked for arxiv.org and some publishers; `curl` works. `pdfminer` breaks on this image
  (cryptography/cffi), so use `pypdf`.
- In `cvxopt`, `cholmod.options['supernodal']=2` forces supernodal $LL^T$. `getfactor` returns the factor of a
  CHOLMOD-chosen permutation; the certificate needs only its row counts and norms.
- Do not run several large FEM jobs in one shell call (OOM).
