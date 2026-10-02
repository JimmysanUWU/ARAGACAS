# Literature check: where the project stands (2026-10-02)

This compares the repository's results with the published literature as of 1 October 2026. Searches covered
arXiv/EMS/AMS/Springer for:
- resolvent degree (Hilbert 13);
- essential dimension with accessory irrationalities;
- gonality of curves with large automorphism group;
- certified Laplace spectra.

The earlier novelty check is `LIT.md` (29 September).

## 1. Hilbert's 13th problem itself (resolvent degree): not moved

- It is still unknown whether $\mathrm{RD}(n)>1$ for any $n$, so Hilbert's conjectures $\mathrm{RD}(6)=2$, $\mathrm{RD}(7)=3$, $\mathrm{RD}(8)=4$ are wide open.
- Recent work:
  - Reichstein, *Hilbert's 13th problem for algebraic groups*, Enseign. Math. 2025: RD is $\le5$ for connected groups, and RD$>1$ is open for every group.
  - Edens–Reichstein, Doc. Math. 2025: in characteristic $p$ all three conjectures can fail.
  - Upper bounds on $\mathrm{RD}(n)$: Sutherland; Heberle–Sutherland, NYJM 2023; Gómez-Gonzáles–Sutherland–Wolfson, J. Algebra 2024.
- **Our work does not touch RD.** Towers of accessories remain open (`QUESTIONS_FOR_GPT.md`, Round 6 Q3).

## 2. The accessory-irrationality version (Klein/Hilbert; Farb–Wolfson 2025): large, sharp advance

**Published state of the art.** Farb–Wolfson, *Essential dimension relative to branched covers of degree at most n*,
arXiv:2510.22786 (26 October 2025), Corollary 1.7:
- $\mathrm{ed}_{\mathbb C}(A_7;\le6)>1$, i.e. there is no reduction of the septic to one parameter after an accessory sextic;
- $\mathrm{ed}_{\mathbb C}(G;\le2)>1$ for every simple $G\ne A_5$;
- for $\mathrm{PSL}_2(\mathbb F_p)$, $n$ up to about $p-1$.

The method is Castelnuovo–Severi with the genus bound $g\le(n-1)^2$. No follow-ups were found.

| quantity | published | this project | status |
|---|---|---|---|
| $\mathrm{ed}_{\mathbb C}(A_7;\le n)>1$ | $n\le6$ (FW 2025) | $n\le59$, sharp: $a(A_7)=60$ | theorem by GPT (DAY 2), verified here (`ACCESSORY_60.md`) |
| same, keeping connected full $A_7$-monodromy | — | exactly $\mu(A_7)=90$ | **ours** (`MU90.md`, `verify_mu90.py`) |
| $a(A_6)$, $a(L_2(7))$, $a(A_5)$ | $\ge5$ and $\ge3$ from FW's formula; $a(A_5)=2$ by Klein | $12$, $4$, $2$ | GPT Round 5, verified here |

For $A_7$ this is a tenfold improvement on a bound published three weeks before our work, with an *exact* value. The new
mechanism is linearised Castelnuovo, Hilbert-90 degree lattices and subgroup induction (DAY 2). For $\mu$ we added
orbit-semigroup base loci and the spin geometry $2.A_7\subset SU(4)$.

**Caveats.**
- None of this is refereed or public.
- The value 60 uses the convention in which the monodromy may drop after the accessory. The bound $\le59$ holds in either convention.
- The upper bound uses Farb–Wolfson's model independence.

## 3. Gonality of $A_7$-curves: new

- **Before this project**, the best available lower bound for the genus-136 $(2,4,7)$ curve was $\operatorname{gon}\ge13$, from Farb–Wolfson's lemma with Conder's genus bound (`LIT.md`).
- **Now:**
  - $\operatorname{gon}\ge25$ is certified on all four classes, so $\operatorname{gon}(C/\langle\tau\rangle)\ge13$.
  - $\operatorname{gon}\ge25$ holds for every faithful $A_7$-curve of genus $\le335$.
  - $\gamma(A_7)\ge25$, modulo GPT's large-genus algebra (not yet re-verified).
- **Methods.** Spectral gonality bounds of Li–Yau/Yang–Yau type are classical: Abramovich 1996 for modular curves, and Cornelissen–Kato–Kool for Drinfeld curves.
  - **Going beyond Li–Yau** through representation theory appears new. The degree is a cubic form on $E_1\otimes\mathbb R^3$, killed by $(\wedge^3E_1)^G=0$; we used it as topological Hersch, and GPT as harmonic Hersch. No precedent was found.
  - Confidence is moderate: the search was by keyword, not exhaustive.

## 4. Certified spectra: new application of a known paradigm

- **Known.**
  - Hejhal's method (1990s) and Strömberg's congruence and vector-valued extensions.
  - Quasimode certification: Booker–Strömbergsson–Venkatesh 2006; the database of rigorous Maass forms of Seymour-Howell et al., arXiv:2502.01442, which covers $\Gamma_0(N)$ only.
  - Bootstrap *upper* bounds on $\lambda_1$ (Kravchuk–Mazáč–Pal 2024: Bolza and Klein quartic).
  - Rigorous $\lambda_1$ exists only in low genus.
- **Ours.**
  - $\lambda_1=0.346267085404$ (12 digits, numerical) and certified two-sided to $2.5\cdot10^{-4}$, for a **genus-136
    non-congruence** cover of the arithmetic group $\Delta(2,4,7)$, in a finite-group-twisted, cocompact, vector-valued setting.
  - Earlier: a CR/Cholesky certificate for 24 more curves of other signatures.
- **Technically new pieces.** These are modest, and the paradigm is BSV's:
  - exact equivariance by a partition of unity of expansions;
  - the Jacobian-flux bound for bracket resolvents.

## 5. Correction found by this check

`LIT.md` said "(2,4,7) is not arithmetic". **That is wrong.** $(2,4,7)$ is in Takeuchi's list of compact arithmetic triangle
groups (Nugent–Voight, arXiv:1510.04637, §6.1.1). The entry is now fixed. Nothing in the proofs depended on it.

## 6. One-line assessment

**Accessory version of Hilbert 13 for $A_7$.** It is settled sharply: $60$, and $90$ with full monodromy. This improves the
published bound ($6$) by an order of magnitude.

**Gonality of $A_7$-curves.** A new certified lower bound of 25, against the previous 13.

**Hilbert 13 proper (RD).** Untouched; it is open worldwide.

All of this is publishable once it is written up and independently audited.

## Sources

- Farb–Wolfson, arXiv:2510.22786 — https://arxiv.org/abs/2510.22786
- Farb–Wolfson, *Resolvent degree, Hilbert's 13th problem and geometry* — https://arxiv.org/abs/1803.04063
- Reichstein, *Hilbert's 13th problem for algebraic groups* — https://arxiv.org/abs/2204.13202
- Edens–Reichstein, *Hilbert's 13th problem in prime characteristic* — https://arxiv.org/abs/2406.15954
- Heberle–Sutherland, *Upper bounds on resolvent degree via Sylvester's obliteration algorithm* — https://arxiv.org/abs/2110.08670
- Sutherland, *Upper bounds on resolvent degree and its growth rate* — https://arxiv.org/abs/2107.08139
- Cornelissen–Kato–Kool, *A combinatorial Li–Yau inequality* — https://arxiv.org/abs/1211.2681
- Abramovich, *A linear lower bound on the gonality of modular curves* — https://arxiv.org/abs/alg-geom/9609012
- Seymour-Howell et al., *A database of rigorous Maass forms* — https://arxiv.org/abs/2502.01442
- Kravchuk–Mazáč–Pal, *Automorphic spectra and the conformal bootstrap* — https://arxiv.org/abs/2111.12716
- Nugent–Voight, *On the arithmetic dimension of triangle groups* — https://arxiv.org/abs/1510.04637
