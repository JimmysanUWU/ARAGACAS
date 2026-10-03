# Independent mathematical audit

Audit date: 3 October 2026. Initial baseline: c60db100bc2dbd2f95cba80c75f47bee56de6bef. The audit incorporates research head de8482869fff4f5c8dcf10333b6e0975cde5e805 on claude/continue-previous-qfhm7j, including its additional exact Schur certificate.

This report covers every numbered statement in the seven current chapters, the computational claims supporting them, the Lean project's actual scope, and the retained families of claims in the six historical documents. Historical proposals are not silently promoted to theorems.

This is a recovery checkpoint of the audit findings. The full corrected source tree, additional exact programs, raw logs, manifest, and 61-page manuscript were prepared in the workspace. The execution service then reported environment_offline before their upload. They are **not** included in this checkpoint, and the existing upstream programs must not be confused with the repaired local snapshot. The conclusions below distinguish analytical correctness, reproduced calculations, and outstanding certification obligations.

## Proof standard

The algebraic results a(A7)=60 and mu(A7)=90 survive the audit, under the generic-torsor convention below. The exact upper bound gon(C)<=42 survives independently of the numerical projective equations. The lower bound 25 has a valid analytical argument and was reproduced using repaired certificates, subject to the stated finite-element and sparse floating-point error model. It is not a fully formalized theorem.

Several headline formulations were too strong. The exact result is a base-point-free morphism birational onto a degree-60 image. Closed immersion, the proposed ideal, and equality 42 for the displayed pencil retain numerical inputs. The elliptic bound concerns specified pure rational Hodge planes. Quotient gonality and the seven possible pencil kernels retain uncertified rank inputs. The old open-base linearisation claim is false and has been replaced by Amitsur theory.

| Label | Meaning |
|---|---|
| P | Mathematical proof with explicit hypotheses |
| X | Exact integer, rational, polynomial or algebraic-number calculation |
| C | Computer-assisted certificate with an explicit analytical and arithmetic trust base |
| N | Floating-point evidence without a complete error certificate |
| H | Historical argument outside the established chain |
| S | Superseded or false in its original formulation |

Dependencies retain their weakest required status. Rounding a character inner product to an integer is N. A successful numerical run does not certify rank, root counts, immersion, or scheme structure.

**Accessory convention.** After monodromy drops, the target is generically free on each component and compresses the actual base-changed torsor with its subgroup H. A faithful permutation of disconnected rational components is insufficient. Essential-dimension comparisons use the generic torsor of a faithful linear representation and this component convention.

## Chapter 1

| Statement | Verdict | Reason / limitation |
|---|---|---|
| Four inner Nielsen classes, two outer types, genus 136 | P+X | Permutation enumeration, generating triples and Riemann–Hurwitz reproduce. Minimum genus uses an exhaustive signature and seven-sheeted test, not new classification theory. |
| Lemma 1.1 | P | Count cosets fixed by an element. |
| Corollary 1.2 | P+X | Involution fixed points include 12 from order-2 inertia and 6 from order-4 inertia. |
| Proposition 1.3 | P+X | Centralizer order 24, quotient H order 12; good/bad involution and quotient counts reproduce. |
| Theorem 1.4 | P | Castelnuovo–Severi forces pencils of degree <=8 to factor through every good involution. |
| Theorem 1.5 | P with classical input | Index-3 and index-1 arguments check. Uniqueness of the gonal rulings on a smooth balanced quadric curve uses Martens. |
| Proposition 1.6 | P+X for representation/isogeny data | Does not determine torsion or canonical height of the named arithmetic points. |
| Theorem 1.7 | P+X | Finite-field smoothness establishes a nonempty characteristic-zero smooth locus. The same local centralizer data allow gonality 9; those data alone cannot exclude it. |

## Chapter 2

| Statement | Verdict | Reason / limitation |
|---|---|---|
| Lemma 2.1 | P | Hersch balancing and Rayleigh quotients. Lean records an abstract implication, not the construction of a balanced holomorphic map. |
| Klein metric and signed tiling | P+X | Common reference mesh and metric tensor check. |
| Lemma 2.2 / representation coverage | P+X | Exact ordinary ATLAS characters and subgroup sums reproduce. Q2 is C3 semidirect D8, not S4. |
| CR lower-bound transfer | P subject to interpolation estimate | Piecewise constant anisotropic comparison has the needed energy orthogonality. The Bessel lower bound received an exact alternating-series/Bernstein proof. |
| Sparse Cholesky / LDL | C | Reproduced under the asserted sparse backward-error model. An implementation-level justification of CHOLMOD's update path and the chosen gamma_k estimate remains in the trust base. Lean assumes the matrix perturbation estimate. |
| lambda1>=0.34089, four classes | C | Reproduced at n=96 using repaired interval operations and actual tested shifts. |
| Classes 12,14: gon>=25 | C | Original n=128 Q1 solves exceeded available memory. Untwisted quotient replacements, transported by an exact odd conjugator, succeeded. Raw bound 0.355695984... must be rounded down, e.g. to 0.35569. |
| E1=14a, next eigenvalue>=0.55998, classes 0,1 | C+X | Q1 inertia, S5 upper and A6 lower bounds isolate one copy of 14a. |
| Ten other signatures and four window signatures | C+X | Required cases reproduced. Optional unused failures do not enter the chain. |
| Hejhal eigenvalues / Bolza | N | Consistency tests, separate from the certified enclosure argument. |

## Chapter 3: harmonic Hersch

| Statement | Verdict | Reason / limitation |
|---|---|---|
| Lemma 3.2 | P | Alternation in the scalar integral and tensor give symmetry of Theta. |
| Lemma 3.3 | P+X | Exact character arithmetic gives (wedge^3 14a)^G=0; the required wedge-square isotypes lie above b. |
| Theorem 3.4 | P | Requires the stated mean-zero projection, a<=lambda1 and lambda_h, b>lambda_h, delta<1, nu>lambda_h, residual and bracket bounds. |
| Lemma 3.5 | P | Consistent orientation makes div J the Poisson bracket. Its dual Dirichlet norm is <=||J||. Spectral comparison contributes b/(b-lambda_h). |
| Residual commutator | P | Translated local solutions satisfy the same exact equation. Sum of cutoff derivatives is zero; subtracting the sector solution gives (Delta chi)D-2 grad chi.grad D. h0=1 does not imply chi0=1. |
| Uniform radial bounds | P+C | Hypergeometric coefficient estimates give uniform envelopes in angular frequency. Positivity and cosine range must hold on every evaluation disk. |
| Mirror identification | P | Reversing the mirrors inverts their rotation product: sigma_BC sigma_AB=Y^-1. Both orientations were checked in exact local-character calculations. |
| Trial space | C+X | Coefficients are exact binary data, not certified eigenfunctions. Averaging/blending gives exact equivariance; positive certified mass proves the 14a copy is nonzero. |
| Theorem 3.1 | P+C | Final 192x128 runs pass: class 0 has lhs 0.330681 > rhs 0.271290; class 1 has 0.329591 > 0.279466. The coarse class-1 run fails. |

Rounded orbit-centre keys were replaced by proven minimum separation of order-7 centres on the bounded disk. Cache provenance includes coefficient bytes, relevant source hashes, precision and radii. The repaired CLI exits with failure when the decisive inequality fails.

## Chapter 4

| Statement | Verdict | Reason / limitation |
|---|---|---|
| Lemma 4.2 | P | Simplicity and induction on degree establish compositum generation. Image action is faithful whenever degree<2520. |
| Lemma 4.3 | P; exposition repaired | Kernel in A6 is normal and has order dividing d; simplicity kills it. Classification of finite subgroups of PGL2 excludes A6. |
| Lemma 4.4 | P | Induction in the newly adjoining function establishes product independence. |
| Chain Castelnuovo bound | P+X | Only feasible td>=2^t-1 count. All t<=1+Omega(d), d<=24, satisfy the comparison. A birational image need not be an embedding. |
| Lemma 4.5 | P | Complete (1,1) pencils have a length-two cluster, proper or tangential. Generic members are smooth at its centres; proximity gives the same convex cost. |
| Corollary 4.6 / Theorem 4.1 | P+C+X | Bounds 331,363,397 and signature enumeration reproduce. Original route uses four window certificates. |
| Stronger genus>=266 route | P+X with classical inputs | Re-derived using canonical comparison, singular graph/222 cases, Petrakiev in P8, uniform position, simultaneous clusters and exact star arithmetic. Its low-genus complement still needs certificates. |

## Chapter 5

| Statement | Verdict | Reason / limitation |
|---|---|---|
| Theorem 5.2 | P | Joint-image closure is Cartier. Picard class splits over the linear base. Restriction of the canonical linearisation to the origin works even when that fibre lies in the divisor. |
| Theorem 5.3 | P | Lower bound uses an H-component. Universal upper bound twists a genuinely linearised bundle and selects a closed point, avoiding a representation-choice shortcut. |
| Lemmas 5.4–5.5 | P | Hilbert 90 gives the degree lattice. Normalizing a nonbirational image would contradict minimal moving linearised degree. |
| Proposition 5.6 | P+X with ATLAS inputs | Signature/subgroup/Castelnuovo calculations reproduce; no numerical projective model is used. |
| Proposition 5.7 | P+X | Lattice excludes 66,78,102 without candidates; A6-fixed section excludes degree-60 candidates. |
| Proposition 5.8 | P+X | Ordinary Lefschetz characters checked in all four classes and both orientations. Positive multiplicity 10 guarantees h0>=10; equality is not asserted. |
| Lemmas 5.9–5.10 | P | Local fibre characters and eigenvalue parity, not element order alone. |
| Exclusion of 72 and 84 | P+X | All 84 exact pattern cases close. Sym2 spin character independently checked by cyclotomic reduction. |
| Rank-2 and rank-3 constructions | P | Compound of rank-two symmetric matrix is rank-one symmetric on wedge2; square root lies in ordinary 6. Adjugate rank-one lies in spin square. Invariant base divisors exceed the whole degree. |
| Halphen (3,5,6), degree 42 | P with classical input | Faithful minimum genus and pi(21,3)=90 force birationality. Irreducible Sym2(V4) excludes quadrics; no-quadric genus bound 274<379. |
| Theorem 5.11 | P | Canonical quotient Pic(XxC)->Hom(Alb X,Jac C) is equivariant. External-product factors are uniquely determined and invariant. Global Hom_G=0 is essential. |
| Example 5.12 | P | Units on an open base can trivialize a constant Schur cocycle despite Pic(B)=0. |
| Theorems 5.14–5.15 | P | Amitsur invariance/no-name are existing results. Stable compression formula requires Hom_G=0. |
| Theorem 5.17 / Corollary 5.18 | P with exact/external degrees | Generic indices 2,3,6 follow from degree gcds, Am=0 and restriction/corestriction. Quadratic and cubic steps remove obstructions while retaining full monodromy. |
| Tower assertion u_Z=0 implies full stable formula | S; repaired locally | Vanishing for one correspondence gives its lower bound. Minimizing over all correspondences requires the global Hom_G hypothesis. |

## Chapter 6

| Statement | Verdict | Reason / limitation |
|---|---|---|
| Algebraic gon(D)>=10 | P+X | Saturated determinant, undivided centralizer average, normalizer generation and free C5 give contradiction. |
| Lemma 6.1 | P with subgroup input | Kernel degree and finite PGL2 quotients classify large pencil stabilizers, giving orbit>=35. |
| Proposition 6.2 | P with Bryant | Constant frame would give a minimal immersion with constant negative curvature; Bryant's primary author summary excludes it. Not in main gonality-25 chain. |
| Proposition 6.3 | P after adding E1^G=0 | Schur orthogonality for the new invariant weight needs the mean-zero hypothesis. |
| Numerical ceiling / route closed | N | No certified global frame minimum. Categorical closure claim withdrawn. Near constancy does not establish a finite spectral band. |
| Retired topological Hersch / shape optimization | Historical P conditional on bounds / N | Superseded diagnostics. |
| Plücker tests | X for checked residues; N for numerical character extraction | Passing an obstruction proves neither existence nor uniqueness. |
| PE,PS | Open | 5-adic nonzero does not prove infinite order. |
| Septic Belyi map | P+X | Exact derivative, branch values, discriminant, reduction and monodromy checks reproduce. |
| Arithmeticity / noncongruence aside | External; unused | Needs specific references before inclusion in a paper. |

## Chapter 7

| Statement | Verdict | Reason / limitation |
|---|---|---|
| Proposition 7.1 | P+X | Extension through the universal cover of PSL2(R) is not the triangle group's universal central extension. Product and inverse correspond to opposite orientations. |
| Corollary 7.2 | P+X | Local datum enumeration, centre order 6 and ordinary degree lattice 90. Reversing orientation exchanges signs of order-3/order-6 residues. |
| Theorem 7.4 | P+X+C | Clifford input depends on lower gonality. Degree-45 invariant polynomials, local orders, 210 lines and coprimality have independent exact checks. Exception: the trivial bundle has a section. Only degree-45 spinor models are excluded. |
| Theorem 7.5 | P+X | Exact ATLAS exceptional 6 and complete central-3 sector give chi(L60)=6-15-2*21-24 in all four classes and both orientations. Six-dimensional subrepresentation gives exact base-point-free birational morphism; h0=6 is not claimed. |
| Corollary 7.6 | P+X | Involution trace 2 and all 18 fibre characters are exact. Moving pencil has degree <=42 and descends to quotient degree <=21. |
| Proposition 7.7 | P conditional on N | Parity/Plücker force simple zeros at twelve 2-points. Immersion at six 4-points and off-fixed-point centralizer-orbit exclusion remain numerical. Immersion alone is insufficient. |
| Proposition 7.8 | N | Sampled ranks and continued roots do not certify embedding, ideal, or scheme structure. |
| Invariant cubic | P+X | Unique invariant cubic vanishes by local characters; ambient exceptional cubic is prior literature. Uniqueness among all cubics containing C still uses restriction rank. |
| Proposition 7.9 | P conditional on nonzero restriction N | Exactly one Klein orbit 24; other orbit sizes 42,84,168. Degree 120 forces 5*24 if quadric restriction is nonzero. |
| Klein differences / Hessian | Conditional / N | Linear equivalences require nonzero restrictions. Relation 6L60~D7 is exact; Hessian identification is numerical. |
| Proposition 7.10 | X in its stated scope | Exact cup form and rational-LDL enumeration establish bound 60 for pure 15/21 rational Hodge planes. Mixed isotypes and additional CM planes are not excluded. |
| Proposition 7.11 | P | Minimal gonality forces h0(M)=2, since otherwise M(-p) gives a cheaper pencil. Stabilizer quotient embeds in PGL2. |
| Lemma 7.12 | P | Galois correspondence works even for S nonnormal in T. |
| Theorem 7.13 / Corollary 7.14 | P conditional on N | Group propagation is exact; Noether/Petri/Koszul ranks are numerical. Petri must include the plane-quintic exception. |
| Constructed pencil stabilizer | Conditional on Proposition 7.7 | Full base divisor then identifies class/field stabilizers. Constructed pencil is not proved gonal. |
| Stable reduction at 5 or 7 | Open | Baker uses metric skeleton with node thicknesses, or appropriate regular semistable subdivision; arbitrary unweighted stable graph is insufficient. |

## Computational repairs and reproduction

Local work added exact literal ATLAS checks in quadratic integer rings, Molien dimensions, full central-3 Lefschetz decomposition, both orientation conventions, exact rational-LDL elliptic enumeration and finite pencil-bound checks. It also checked integer overflow before every int64 multiplication in the new upstream Schur program.

Local certificate repairs used Arb enclosure maxima/minima rather than comparisons of overlapping balls. Factor loops record the shift actually used, not a hypothetical untested retry. The LDL routine can shrink an overly conservative valid shift and refactor. Hersch keys, cache provenance and failure exit status were repaired.

The evidence index prepared locally contains **50 run records**, versions and final-snapshot source/data hashes. These hashes were not captured separately for every earlier intermediate run. Failed coarse Hersch and memory-exhausted jobs were retained. The final Hersch inequalities above, exact 84-case mu90 result, exact ATLAS sector calculation, and all four elliptic searches succeeded. Numerical projective/quotient runs retained N status.

Large certificate runs should be sequential. Class-12 transport uses exact conjugator (0,1,2,3,4,6,5) to class-14 geometry. Subgroup induction characters are 1+6+14a for S5 and 1+6+14b+21 for L2(5). With Q0/Q2 these cover all ordinary representations. Replacement second-eigenvalue bounds were 0.3556959843961548 and 0.3599833853937989, above 48/135.

**Lean.** Source inspection finds no sorry/admit in the two files. Formal scope is three elementary superposition statements plus abstract implications in the spectral certificate. It does not include CHOLMOD implementation, FEM assembly, analytic Hersch estimates, or a full proof of gamma(A7)>=25. The toolchain became unavailable after an earlier environment reset; a fresh build is **not independently reproduced**.

## Historical material and open frontier

Six historical PDFs were extracted and read. DAY1A–C are checkpoints/verification charters. DAY2's accessory correspondence, subgroup induction and transport have current proofs. Round5's harmonic identity and mu>=72 are integrated. The Amitsur note corrects the earlier open-base error.

The stronger multigraded genus theorem and genus>=266 algebraic route were re-derived. Historical Picard closure, residual transport, Hurwitz-group and conformal-attainment proposals have separate hypotheses and should not be promoted by association. Torsion/norm collision models are superseded by the degree-9 exclusion. Balancing laws above saturation remain proposals.

Exact gonality in [25,42], bad-prime stable reduction, a certified lower bound 26, towers with overlapping H1 constituents and infinite order of PE/PS remain open.

Publication-critical work is delimited: independent review of the algebraic chain, implementation-level sparse error justification or outward residual certificates, and exact rank/root/scheme certificates for stronger projective claims. This checkpoint records the assessment; it does not supply the unavailable full repaired source package.
