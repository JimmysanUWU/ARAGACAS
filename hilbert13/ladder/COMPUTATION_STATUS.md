# Recovery computation status

The interrupted environment was reset. Its numerical export files did not
survive. The continuation proof and independent checker were reconstructed
from the session checkpoint and the unchanged source commit
`dbcd7ea788d0b1d59b92bedb75e72af5b8af6969`.

A fresh sequential run of all twelve cases is in progress. Eight Q0/Q2 cases
and Q1 class 0 have completed. The large Q1 classes 1,12,14 are not yet included
as completed in this initial checkpoint. The eventual JSON records and
permutations under `residual_results/` supersede the rounded pre-reset table.
The original `certificate.txt` is Claude's certificate, not a substitute for
the independent residual records.

Completed exact checks:

- `weighted_budget_checks.py`: the genus-336 table, residue formulas for the
  sharp independent-triple bound, genus-266 star budget, seven target thresholds.
- `signature_checks.py`: all eleven signatures through genus 335, their exact
  generating-triple counts, the four-point genus bound.
- `belyi_quotient_checks.py`: critical-value identities, discriminant identity,
  integral polynomial scaling, mod-17 monodromy witness, mod-457 quotient separation.
- `fixed_divisor_exact.py`: integer Specht characters, Fraction projection norms,
  both quotient carriers and their detecting conjugates.
- `psl13_transport_check.py`: three Klein subgroups, their A4 normalizers, and
  generation by every pair.
- `residual_checker_tests.py`: a nonidentity sparse permutation, exact rational
  residual comparison, a deliberately corrupted factor that is rejected, dense
  constant deflation, and the Bessel-zero bound by rational Bernstein coefficients.

Reproduce the numerical run from this directory with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 residual_certificate.py --out residual_results
```

Run cases sequentially: the full Q1 factor is large. A `--quick` result uses
sigma=.3 and cannot be presented as the .34089 certificate. Completion is
reported only after the case's JSON and permutation files have been written.
