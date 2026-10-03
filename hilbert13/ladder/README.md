# The $A_7$ gonality ladder

Curves with an $A_7$-action, their gonality, and the accessory-irrationality problem for the general septic.

The [textbook edition](gpt/textbook/README.md) distills the complete source collection into a unified mathematical treatment. Its status labels distinguish exact inputs, analytical certificates, numerical proposals and historical reductions.

## Main results

$C$ is a $(2,4,7)$ $A_7$-curve (genus 136, the minimum), $\tau$ an involution, and $D=C/\langle\tau\rangle$.

| statement | status | where |
|---|---|---|
| **$25\le\operatorname{gon}(C)\le42$ and $13\le\operatorname{gon}(D)\le21$** | | |
| — lower: $\lambda_1\ge0.34089$ with Li–Yau; classes 12, 14 directly, classes 0, 1 by harmonic Hersch | [P][C]; specified logical cores [L] | Ch. 2–3 |
| — upper: the sections of a degree-60 Schur-twisted class give a base-point-free pencil of degree exactly 42, by the embedding and exact normalization-defect argument | [P][X] | Ch. 7, §7.12 |
| **$25\le\gamma(A_7)\le42$**: the least gonality of a faithful $A_7$-curve. Genus $\le335$ is spectral; genus $\ge266$ is algebraic (the three-pencil theorem) | [P][C] | Ch. 4, 7 |
| **$a(A_7)=60$**: $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$, sharp (GPT; verified) | [P][X] | Ch. 5 |
| **$\mu(A_7)=90$**: compression with connected full monodromy needs exactly degree 90 | [P][X] | Ch. 5 |
| over a smooth projective generically free base with $\mathrm{Hom}_G(\mathrm{Alb},\mathrm{Jac}\,C)=0$, the least stable compression degree is $\mu_{\mathrm{Am}(X)}(C)$: 90 or 60 for $(2,4,7)$ targets, by the Amitsur subgroup | formula [P]; numerical threshold [P][X][C] | §5.5 |
| invariant degrees are $15\mathbb Z$; below 60 only $\mathcal O_C$ has a section; Mumford class in $\mathbb Z/6$ | lattice [P][X]; section exclusion [P][X][C] | §§7.2–7.4 |
| $\varphi:C\hookrightarrow\mathbb P^5$ is a degree-60 closed embedding on the Laza–Zheng cubic, with no quadrics. Being cut out by the cubic and 15 quartics remains [N] | embedding/no quadrics [P][X]; ideal [N] | §§7.6, 7.12 |
| Klein subgroup quadrics have exact divisor $5\times24$ seven-points; $2L_{60}\sim5O_H$, $6L_{60}\sim D_7$. Their Plücker/Veronese interpretation is in Prop. 7.9 | [P][X] | Prop. 7.9, §7.12 |
| the two exceptional quadratic systems are the Jacobian complete intersection of the cubic and the dual cubic's apolar ideal, with quotient Hilbert series $(1+t)^6$ and function $(1,6,6,1)$ | [P][X] | §7.12 |
| a gonal pencil's class stabiliser acts on it through $PGL_2$, with kernel $N$ costing $|N|\operatorname{gon}(C/N)$; all quotients $C/K$ with $|K|\le60$ are priced, so a pencil below 42 has $N\in\{1,C_2,C_3,C_4,V_4,S_3,C_7\}$ | [P][N] | §7.10 |
| elliptic subcovers with pure rational tensor Hodge planes in types 15, 21 have degree $\ge60$; $H^1(C,\mathbb Z)$ exactly from the dessin. Mixed or extra CM planes are not covered | [P][X] | §7.9 |
| $\operatorname{gon}(D)\ge9$ by Castelnuovo–Severi; $\ge10$ by ramification transport | [P][X] | §1.3, §6.1 |

**Tags.**
- **[P]** paper proof.
- **[X]** finite computation in integers, rationals, finite fields or specified algebraic number rings. Rounded floating-point character recognition retains [N] unless replaced by an independent exact calculation.
- **[C]** computer-assisted proof: ball arithmetic, or floating point with a-priori error bounds.
- **[L]** checked in Lean (`../Hilbert13/SpectralCertificate.lean`).
- **[N]** numerical only.
- **[G]** GPT's proof only, not re-derived here (`gpt/README.md` §2).

## Chapters

| file | contents |
|---|---|
| [`1_CURVE.md`](1_CURVE.md) | group data, fixed points, $D$, $\operatorname{gon}(D)\ge9$, degree 9 and the audit curve |
| [`2_SPECTRAL.md`](2_SPECTRAL.md) | Li–Yau, exact hyperbolic model, twisted quotients, certified $\lambda_1$ for 11 signatures |
| [`3_HARMONIC_HERSCH.md`](3_HARMONIC_HERSCH.md) | the cubic-form lemma, harmonic Hersch, vector-valued Hejhal, certificate for classes 0, 1 |
| [`4_LARGE_GENUS.md`](4_LARGE_GENUS.md) | conjugate pencils, Castelnuovo bounds, cost of a dependent third, the window |
| [`5_ACCESSORY.md`](5_ACCESSORY.md) | correspondence, degree lattice, $a=60$, $\mu=90$, the twisted correspondence |
| [`6_SIDE_RESULTS.md`](6_SIDE_RESULTS.md) | algebraic $\operatorname{gon}(D)\ge10$, pencil orbits, conformal ceiling, negative tests, open arithmetic |
| [`7_TWISTED.md`](7_TWISTED.md) | Schur-twisted classes, the $\mathbb P^5$ model and its equations, $\operatorname{gon}\le42$, Klein quadrics, elliptic subcovers, pencils by symmetry type |
| [`LITERATURE.md`](LITERATURE.md) | the literature and references |
| [`QUESTIONS_FOR_GPT.md`](QUESTIONS_FOR_GPT.md) | audit requests and open questions |
| [`gpt/README.md`](gpt/README.md) | digest of GPT's five PDFs (DAY 1A–C, DAY 2, Round 5): what each proves, where it is used, what is superseded |

## Scripts

Every claim has a script and a saved output, `NAME_output.txt` unless noted.

| script | claim |
|---|---|
| `a7.py`, `chartab.py`, `triples_data.py` | $A_7$ toolkit: permutations, the exact ATLAS character table, the four classes of triples |
| `curve_checks.py` | Ch. 1: triples, fixed points, quotients of $D$, Jacobian, $(\star)$ isotypes |
| `model99.py` | the audit curve (Thm 1.7) |
| `orbifold.py`, `cover.py` | exact Klein-chart model; $Q_0,Q_1,Q_2$ see every irreducible (§§2.2–2.3) |
| `certify.py`, `run_certificate.py` | $\lambda_1\ge0.34089$ (§2.5) → `certificate.txt` |
| `certify_th.py` | $E_1=14_a$, $\lambda'\ge0.55998$; classes 12, 14: $\lambda_1\ge0.355695$ |
| `signatures_spectrum.py`, `certify_signatures.py` | other signatures (§2.5) and the window (§4.4) → `certify_{signatures,window}_output.txt` |
| `validate_bolza.py` | the Bolza surface (§2.6) → `validation_bolza.txt` |
| `verify_exact.py` | $(\wedge^314_a)^{A_7}=0$ (Lemma 3.3) |
| `hejhal_solve.py`, `hh_eval.py`, `hh_certify.py` | harmonic Hersch (Ch. 3) → `coef_cls*_M90.npz`, `hh_certify_output_cls{0,1}.txt` |
| `gonality_large_genus.py` | $B^*(d)$ and the window (Ch. 4) |
| `verify_accessory60.py`, `equivariant_rr.py`, `verify_mu90_exact.py` | $a=60$; $h^0(B+T)\ge10$; the exclusion of 72, 84 (Ch. 5) → `mu90_exact_output.txt` for the last |
| `side_checks.py` | §6.1, §6.4, §6.5 (with GPT's Belyi map), Prop. 1.6 and the Ch. 4 table |
| `twisted_rr.py`, `twisted_survey.py` | $2.A_7,3.A_7,6.A_7$, twisted Lefschetz, §§7.1–7.6; the 26 rigid curves (§7.8) |
| `audit_exact.py` | exact (GPT's audit): the ATLAS $\mathbf 6$ of $3.A_7$ and $V_4$ of $2.A_7$; Molien series; the central-3 character sector; $\chi(L_{60})$ with $\mathbf 6\subseteq H^0$ and eigenvalue $+1$ at the 18 fixed points (Thm 7.5, Cor. 7.6); degree-45 data (Thm 7.4); Prop. 5.8 |
| `tau_pencil.py` | no base points beyond the 18 (Prop. 7.7) |
| `p5_curve.py` | equations of $\varphi(C)$, embedding, the fibre $18+42$, Klein quadrics (Props. 7.8–7.9) |
| `elliptic_subcovers.py` | $H^1(C,\mathbb Z)$ with cup product; subcovers $\ge60$ (Prop. 7.10) |
| `quotient_gonality.py` | differentials on $C/K$; Noether, Petri, $K_{p,2}$; Castelnuovo–Severi over the subgroup lattice (Thm 7.13) |
| `gpt/verify_normalization.py` | exact order-7 eigenline stabilizers and normalization-defect arithmetic: embedding and exact degree-42 pencil (Prop. 7.15, 7.7) |
| `gpt/verify_quadrics.py` | exact quadratic-system multiplication ranks and Klein restriction characters: no quadrics, Jacobian/apolar ideals, Klein divisors (Props. 7.16, 7.9) |

## Reproduction

```sh
pip install numpy scipy sympy python-flint cvxopt pypdf
cd hilbert13/ladder
python3 curve_checks.py > curve_checks_output.txt  # 15 s
python3 side_checks.py > side_checks_output.txt    # 30 s
python3 run_certificate.py                         # ~6 min, 7 GB RAM
python3 certify_th.py count 96 0 0.56              # see certify_th_output.txt for the other runs
python3 certify_signatures.py 24 3 3 5 2 5 7 3 3 6 3 4 4 2 6 7 3 3 7 2 7 7 3 4 5 3 4 6 4 4 4 > certify_signatures_output.txt   # <= 12 min
python3 certify_signatures.py 24 3 4 7 3 5 5 3 5 6 3 5 7 3 6 6 3 6 7 3 7 7 4 4 5 4 4 6 4 4 7 > certify_window_output.txt       # <= 25 min; stops at (4,4,7) (§4.4)
python3 hh_certify.py 0 coef_cls0_M90.npz 192 128  # ~6.5 min per class; likewise class 1
for f in cover verify_exact gonality_large_genus verify_accessory60 equivariant_rr model99; do python3 $f.py > ${f}_output.txt; done
python3 verify_mu90_exact.py > mu90_exact_output.txt
python3 twisted_rr.py > twisted_rr_output.txt               # 25 s
python3 twisted_survey.py > twisted_survey_output.txt       # 90 s
python3 audit_exact.py > audit_exact_output.txt             # 8 s
python3 tau_pencil.py > tau_pencil_output.txt               # 1 min
python3 p5_curve.py > p5_curve_output.txt                   # 7 min
python3 elliptic_subcovers.py > elliptic_subcovers_output.txt   # 80 s
python3 quotient_gonality.py 0 > quotient_gonality_output.txt  # 9 min; class 12 gives the same table
python3 gpt/verify_normalization.py > gpt/normalization_check.txt
python3 gpt/verify_quadrics.py > gpt/quadrics_check.txt
```

Times are historical estimates for an idle 4-core machine; run the FEM jobs one at a time (7 GB each). Saved outputs are retained source records. The textbook verification manifest records which computations were rerun for this edition.

## Trust base

- **Classical theorems.**
  - Castelnuovo–Severi, Castelnuovo's bound, Halphen–Gruson–Peskine;
  - Hersch balancing, Li–Yau, min–max, Schur;
  - holomorphic Lefschetz, Chevalley–Weil, Coppens–Martens;
  - Bryant (§6.3 only).
- **Finite-element bounds.** Crouzeix–Raviart (Liu; Carstensen–Gedicke) and Higham's Cholesky backward error, both checked in Lean in abstract form.
- **Software.** python-flint ball arithmetic; CHOLMOD with a-priori error bounds.
- **External inputs.** GPT's theorems, re-derived wherever used (`gpt/README.md` §1); the ATLAS character table, subgroup lists and matrices `3A7G1-Ar6B0`, `2A7G1-Ar4aB0` (checked in `cover.py`, `audit_exact.py`).
