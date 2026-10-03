# Hilbert's thirteenth problem and curves with an alternating action

This edition reorganizes the mathematics of the entire repository by definitions, hypotheses, statements and proofs. The research chapters, programs, saved outputs and historical manuscripts remain the source records. Numbered results in this edition use their own numbering; the source correspondence will be recorded in the verification appendix.

All algebraic curves are smooth, projective and connected over $\mathbb C$ unless an image, normalization or different ground field is specified. A group action is faithful unless stated otherwise. Arithmetic genus and geometric genus are distinguished throughout.

The edition is being committed in mathematical sections. The completed edition will include the analytic arguments, accessory theory, twisted geometry, arithmetic models, historical alternatives and a file-by-file source map, together with a typeset PDF.

## Reading the statements

| Label | Meaning |
|---|---|
| P | A mathematical argument is supplied, with its hypotheses and named classical inputs. |
| X | A finite calculation in integers, rationals, finite fields or specified algebraic number rings supplies an input. |
| C | An analytical calculation uses explicit enclosures and a stated floating-point or interval-arithmetic error model. |
| N | Uncertified numerical evidence, including rounded character recognition, sampled ranks and numerical root counts. |
| H | A historical argument or proposal that has not been incorporated into the established chain. |

An assertion retains every status required by its dependencies. A formal theorem proves its displayed implication; it does not formalize an application whose hypotheses are supplied outside the prover. The term *exact* never means rounding a floating-point result to the nearest integer.

## Sections available in this checkpoint

- [1. Function classes and curve invariants](01_FOUNDATIONS.md)
- [2. The alternating curves and ramification transport](02_CURVES.md)
- [3. Spectral lower bounds](03_SPECTRAL.md)
- [4. The cubic degree form and harmonic Hersch](04_HARMONIC.md)
- [5. Multiple pencils and multigraded genus](05_MULTIPLE_PENCILS.md)
- [6. Linearized moving degree and accessories](06_ACCESSORIES.md)
- [7. Twisted Picard classes and the exceptional model](07_TWISTED_GEOMETRY.md)
- [8. Normalization defects and quadratic algebras](08_NORMALIZATION_AND_QUADRICS.md)
- [9. Jacobians, arithmetic quotients and pencil symmetry](09_JACOBIANS_AND_SYMMETRY.md)
- [10. Conformal variation, other groups and historical reductions](10_OTHER_FRAMEWORKS.md)

`verify_textbook_inputs.py` checks the exact arithmetic and character inputs introduced by this edition. The source index, complete verification record and typeset PDF are the remaining assembly work.

The general structures developed here are the orbit bounds for normalization defects and the Jacobian/apolar decomposition of the exceptional quadratic representation. Their classical ingredients and their specific finite certificates are stated separately.
