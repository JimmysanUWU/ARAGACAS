# The $A_7$ gonality ladder

Curves with an $A_7$-action, their gonality, and the accessory-irrationality problem for the general septic.

## Main results

$C$ denotes a $(2,4,7)$ $A_7$-curve (genus 136), and $\tau$ an involution.

| statement | status | where |
|---|---|---|
| $g(C)=136$, the minimum; $g(C/\langle\tau\rangle)=64$; $\operatorname{gon}(C/\langle\tau\rangle)\ge9$ | [P][X] | Ch. 1 |
| $\lambda_1(C)\ge0.34089$, so $\operatorname{gon}(C)\ge24$; classes 12, 14: $\operatorname{gon}\ge25$ | [C][L] | Ch. 2 |
| classes 0, 1: $\operatorname{gon}(C)\ge25$ (harmonic Hersch) | [P][C] | Ch. 3 |
| **every $(2,4,7)$ curve: $\operatorname{gon}(C)\ge25$ and $\operatorname{gon}(C/\langle\tau\rangle)\ge13$** | [P][C] | Ch. 2–3 |
| **$\gamma(A_7)\ge25$**: every faithful $A_7$-curve has gonality $\ge25$ (always $\ge23$ by algebra alone) | [P][C] | Ch. 4 |
| **$a(A_7)=60$**: $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$, sharp (GPT; verified) | [P][X] | Ch. 5 |
| **$\mu(A_7)=90$**: connected full-monodromy compression needs exactly degree 90 | [P][X] | Ch. 5 |
| the compression bounds hold without fixed points, over bases with $\mathrm{Hom}_G(\mathrm{Alb},\mathrm{Jac}\,C)=0$; over such bases gonality bounds the degree from below | [P] | §5.5 |
| $\operatorname{gon}(C/\langle\tau\rangle)\ge10$ by pure algebra (Sol/Astra) | [P][X] | §6.1 |
| Schur-twisted invariant classes have degrees $15\mathbb Z$, none with sections below degree 60; a birational $3.A_7$-model $C\to\mathbb P^5$ of degree 60 on the Laza–Zheng $A_7$-cubic fourfold | [P][X] | Ch. 7 |
| **$\operatorname{gon}(C)\le42$, $\operatorname{gon}(C/\langle\tau\rangle)\le21$** (previously 56, 28); hence $25\le\gamma(A_7)\le42$ | [P][X] | Ch. 7 |
| the $\tau$-pencil is exactly a base-point-free $g^1_{42}$; elliptic subcovers (types $15$, $21$) have degree $\ge60$, from $H^1(C,\mathbb Z)$ computed exactly | [P][X] | §7.5, §7.9 |

**Tags.**
- **[P]** paper proof.
- **[X]** finite computation: exact arithmetic, or floating-point character sums whose integrality is checked.
- **[C]** computer-assisted proof: ball arithmetic, or floating point with rigorous a-priori error bounds.
- **[L]** checked in Lean (`../Hilbert13/SpectralCertificate.lean`).
- **[N]** numerical only.

## Chapters

| file | contents |
|---|---|
| [`1_CURVE.md`](1_CURVE.md) | group data, fixed points, $D=C/\langle\tau\rangle$, $\operatorname{gon}(D)\ge9$, degree-9 structure, the audit curve |
| [`2_SPECTRAL.md`](2_SPECTRAL.md) | Li–Yau, exact hyperbolic model, twisted quotients, CR lower bounds, certified results for 11 signatures |
| [`3_HARMONIC_HERSCH.md`](3_HARMONIC_HERSCH.md) | the cubic-form lemma, harmonic Hersch, vector-valued Hejhal, certificate for classes 0, 1 |
| [`4_LARGE_GENUS.md`](4_LARGE_GENUS.md) | conjugate pencils, Castelnuovo bounds, dependent-third cost, the window signatures |
| [`5_ACCESSORY.md`](5_ACCESSORY.md) | correspondence, degree lattice, $a(A_7)=60$, $\mu(A_7)=90$ |
| [`6_SIDE_RESULTS.md`](6_SIDE_RESULTS.md) | algebraic $\operatorname{gon}(D)\ge10$, pencil orbits, conformal ceiling, negative results, open arithmetic |
| [`7_TWISTED.md`](7_TWISTED.md) | Schur-twisted invariant line bundles, twisted Lefschetz, the $\mathbb P^5$ model, $\operatorname{gon}\le42$ and its sharpness for the pencil, elliptic subcovers |
| [`LITERATURE.md`](LITERATURE.md) | comparison with the literature, and references |
| [`QUESTIONS_FOR_GPT.md`](QUESTIONS_FOR_GPT.md) | open questions and audit requests |
| `gpt/` | GPT's source documents (DAY 2 reference, Round-5 answers) |
| [`archive/`](archive/README.md) | superseded notes, reviews and scripts |

## Scripts, claims and outputs

| script | claim | output |
|---|---|---|
| `a7.py`, `chartab.py`, `classreps.py`, `triples.py`, `triples_data.py` | $A_7$ toolkit; the 4 classes of $(2,4,7)$ triples | — |
| `rung1.py`, `quot.py`, `mingenus.py` | fixed points, quotient genera, minimum genus (§1.1–1.2) | stdout |
| `jacobian.py`, `hdecomp.py`, `cusp.py`, `star.py` | Jacobian decomposition, $(\star)$ isotypes (§1.4) | stdout |
| `model99.py` | the audit curve $D'$ (Thm 1.9) | stdout |
| `orbifold.py` | exact Klein-chart model and twisted FEM (§2.2) | — |
| `cover.py` | $Q_0,Q_1,Q_2$ see every irreducible (§2.3) | stdout |
| `certify.py`, `run_certificate.py` | $\lambda_1\ge0.34089$ (§2.5) | `certificate.txt`, `certificate_quick.txt` |
| `certify_th.py` | $E_1=14_a$, $\lambda'\ge0.55998$; classes 12, 14: $\lambda_1\ge0.355696$ | `certify_th_output.txt` |
| `signatures_spectrum.py`, `certify_signatures.py` | the ten other signatures (§2.5) and the window (§4.4) | `certify_signatures_output.txt`, `certify_window_output.txt` |
| `validate_bolza.py` | external check on the Bolza surface (§2.6) | `validation_bolza.txt` |
| `verify_exact.py` | $(\wedge^314_a)^{A_7}=0$ and other exact character facts (Lemma 3.3) | stdout |
| `hejhal_solve.py`, `hh_eval.py`, `hh_certify.py` | harmonic Hersch certificate (Ch. 3) | `coef_cls*_M90.npz`, `hh_certify_output_cls{0,1}.txt` |
| `gonality_large_genus.py` | $B^*(d)$ table and window listing (Ch. 4) | `gonality_large_genus_output.txt` |
| `verify_accessory60.py` | $a(A_7)=60$ inputs (§5.2) | stdout |
| `equivariant_rr.py` | holomorphic Lefschetz: $h^0(B+T)\ge10$ (Prop. 5.8) | stdout |
| `verify_mu90_exact.py` | exclusion of 72 and 84, exact (§5.4) | `mu90_exact_output.txt` |
| `twisted_rr.py` | $2.A_7$, $3.A_7$, $6.A_7$; twisted Lefschetz; Theorems 7.4–7.5, Corollary 7.6, the cubic fourfold (Ch. 7) | `twisted_rr_output.txt` |
| `twisted_survey.py` | the same on all 26 rigid curves of genus $\le335$ (§7.8) | `twisted_survey_output.txt` |
| `tau_pencil.py` | the $\tau$-pencil has no base points beyond the 18 (Prop. 7.7) | `tau_pencil_output.txt` |
| `elliptic_subcovers.py` | $H^1(C,\mathbb Z)$ with cup product from the dessin; elliptic subcovers have degree $\ge60$ (Prop. 7.8) | `elliptic_subcovers_output.txt` |
| `review_checks.py`, `frontier_checks.py` | finite inputs of §6.1, §6.6, §6.7 | stdout |

## Reproduction

```sh
pip install numpy scipy sympy python-flint cvxopt pypdf
cd hilbert13/ladder
python3 cover.py                       # 1 s
python3 run_certificate.py             # ~6 min, 7 GB RAM -> certificate.txt
python3 certify_th.py count 96 0 0.56 # likewise t = 1; s5 24 t; a6 16 t 1.0 (commands in certify_th_output.txt)
python3 certify_signatures.py 24 3 3 5 2 5 7 3 3 6 3 4 4 2 6 7 3 3 7 2 7 7 3 4 5 3 4 6 4 4 4
python3 hh_certify.py 0 coef_cls0_M90.npz 192 128   # ~6.5 min per class
python3 hh_certify.py 1 coef_cls1_M90.npz 192 128
python3 gonality_large_genus.py
python3 verify_mu90_exact.py
python3 twisted_rr.py > twisted_rr_output.txt           # ~25 s
python3 twisted_survey.py > twisted_survey_output.txt   # ~90 s
python3 tau_pencil.py > tau_pencil_output.txt           # ~1 min
python3 elliptic_subcovers.py > elliptic_subcovers_output.txt   # ~80 s
```

Run large FEM jobs one at a time, because several in one shell can run out of memory.

## Trust base

- **Classical theorems.**
  - Castelnuovo–Severi, Castelnuovo's bound, Halphen–Gruson–Peskine;
  - Hersch balancing, Li–Yau, min–max, Schur;
  - holomorphic Lefschetz, Chevalley–Weil, Coppens–Martens ($\operatorname{gon}\le\mathrm{Cliff}+3$);
  - Bryant (1985), only in §6.3.
- **Finite-element bounds.** The Crouzeix–Raviart lower bound (Liu; Carstensen–Gedicke) and Higham's Cholesky backward error. Both are checked in Lean in abstract form.
- **Software.** python-flint (arb) ball arithmetic, and CHOLMOD for the factorisations whose errors are bounded a priori.
- **External inputs.** GPT's DAY 2 and Round 5 theorems, all re-derived (Ch. 5); the ATLAS subgroup lists.
