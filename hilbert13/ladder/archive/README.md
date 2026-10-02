# Archive

This folder holds superseded material, kept for provenance. Every result here that is still valid has been restated in the
chapters (`../1_CURVE.md` to `../6_SIDE_RESULTS.md`). Cross-references inside these files use the old file names:

| old name | now |
|---|---|
| `NOTES.md` | `1_CURVE.md`, `2_SPECTRAL.md` |
| `HH_CERTIFICATE.md` | `3_HARMONIC_HERSCH.md` |
| `GONALITY_LARGE_GENUS.md` | `4_LARGE_GENUS.md` |
| `MU90.md`, `ACCESSORY_60.md` | `5_ACCESSORY.md` |
| `LIT.md`, `LITERATURE_CHECK.md` | `LITERATURE.md` |

## `docs/`

| file | what it was | superseded by |
|---|---|---|
| `VERIFICATION.md` | fifth-session audit of topological Hersch and the certificates | Ch. 2 (§2.5) and Ch. 3 |
| `FRAMEWORK_TOPOLOGICAL_HERSCH.md` | the first degree-based inequality, its numerics, and the degree lattice | Ch. 3 (Lemma 3.3), §6.4, Lemma 5.4 |
| `FRAMEWORK_CONFORMAL.md` | the conformal-metric route, and the pencil-orbit and $Q_2$ results | §6.2, §6.3, §6.6, §6.7 |
| `ACCESSORY_60.md` | verification of GPT's $a(A_7)=60$ | §5.1–5.2 |
| `ROUND5_REVIEW.md` | claim-by-claim review of GPT's Round-5 answers; $\mu\le90$ | §5.3, Ch. 3 |
| `REVIEW_GPT.md` | review of the Sol/Astra *Ramification transport* synthesis | §6.1 |
| `LITERATURE_CHECK.md` | 2 October state-of-the-art check | `LITERATURE.md` |
| `QUESTIONS_FOR_GPT.md` | question rounds 2–7, as sent | `QUESTIONS_FOR_GPT.md` (open questions only) |

## `scripts/`

These scripts import the toolkit in `../..` (`a7.py`, `orbifold.py`, ...). Run them from `hilbert13/ladder` as
`PYTHONPATH=. python3 archive/scripts/NAME.py`.

| file | what it was | superseded by |
|---|---|---|
| `verify_mu90.py`, `mu90_output.txt` | first $\mu=90$ check; the rank steps used floating point | `verify_mu90_exact.py` |
| `topo_hersch.py`, `topo_hersch_extras.py`, `topo_hersch_output.txt` | constants of topological Hersch (P1 FEM, not certified) | `hh_certify.py` |
| `conformal.py` | frame function $\bar F$ and the conformal ascent (§6.3 numerics) | — (negative result) |
| `equivariant_plucker.py` | equivariant Plücker test (§6.5) | — (negative result) |
| `equivariant_rr_general.py` | holomorphic Lefschetz on other signatures (no positive part at 72 or 84) | §5.4 |
| `spectrum.py` | first-session FEM with a flat approximation | `orbifold.py` |
