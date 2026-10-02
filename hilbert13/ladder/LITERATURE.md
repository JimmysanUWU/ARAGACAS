# Literature

**Searches.** 29 September 2026 (novelty) and 2 October 2026 (state of the art).
- Sources searched: web, arXiv abstracts, and EMS/AMS/Springer.
- Not searched: MathSciNet or zbMATH.

**Tags.**
- **[read]**: the full text was read.
- **[abstract]**: the abstract or a summary only.
- **[unverified]**: cited from memory; check before relying on it.

## 1. Where the project stands

| question | published | this project |
|---|---|---|
| Hilbert 13 proper: is $\mathrm{RD}(n)>1$ for some $n$? | open for every $n$ | not touched |
| $\mathrm{ed}_{\mathbb C}(A_7;\le n)>1$ | $n\le6$ (Farb–Wolfson 2025) | $n\le59$, sharp: $a(A_7)=60$ (GPT, DAY 2; verified here, Ch. 5) |
| same, with connected full $A_7$-monodromy | — | exactly $\mu(A_7)=90$ (ours, Ch. 5) |
| $a(A_6)$, $a(L_2(7))$, $a(A_5)$ | $\ge5$, $\ge3$ from FW; $a(A_5)=2$ (Klein) | $12$, $4$, $2$ (GPT Round 5; verified) |
| gonality of the genus-136 $(2,4,7)$ curves | $[13,56]$ (FW Lemma 2.2 with Conder's genus; a quotient map) | $[25,42]$ (Ch. 2–3, 7) |
| $\gamma(A_7)$, the least gonality of a faithful $A_7$-curve | $\ge13$ | $[25,42]$ (Ch. 4, 7) |
| certified $\lambda_1$ of a closed hyperbolic surface | genus $\le7$ | genus 136 and 24 further $A_7$-curves (Ch. 2) |

**Assessment.**
- **Accessory version of Hilbert 13 for $A_7$.** It is settled sharply: 60, and 90 with full monodromy.
  This improves the published bound (6) by an order of magnitude.
- **Gonality.** The genus-136 curves have gonality in $[25,42]$, against $[13,56]$ before.
- **Methods.** Going past Li–Yau by representation theory (Lemma 3.3: the degree is a cubic form killed by $(\wedge^3E_1)^G=0$) appears to have no
  precedent. The certification paradigm is that of Booker–Strömbergsson–Venkatesh. Our additions are modest:
  - exact equivariance by a partition of unity;
  - the Jacobian-flux bound.
- **Schur-twisted models** (Ch. 7). Holomorphic Lefschetz for line bundles linearised only for $6.A_7$, with an exact degree formula from the
  universal central extension of the triangle group. We found no prior use of Schur-multiplier twists to build low-degree models of symmetric curves.
- **Caveats.** None of this is refereed. The value $a=60$ uses the convention in which the monodromy may drop after the accessory; the bound
  $\le59$ holds under either convention.

## 2. Hilbert 13 and essential dimension

- **B. Farb, J. Wolfson, *Essential dimension relative to branched covers of degree at most n*, arXiv:2510.22786 (Oct 2025) [read].** The backbone.
  - **Main Theorem 1.6.** $\mathrm{ed}_k(G;\le n)>1$ if three conditions hold:
    - no proper subgroup of $G$ has index $\le n$;
    - some $M\hookrightarrow\mathrm{PSL}_2(k)$ has $|M|>n$;
    - $G$ acts on no curve of genus $\le(n-1)^2$.
  - **Cor. 1.7.** $\mathrm{ed}_{\mathbb C}(A_7;\le6)>1$, and $\mathrm{ed}_{\mathbb C}(G;\le2)>1$ for every simple $G\ne A_5$.
  - **Tools.**
    - Lemma 2.2: if functions of degree $n$ generate $k(C)$, then $g\le(n-1)^2$.
    - Lemma 2.4: the gonality of a quotient is at most the gonality of the cover.
    - Lemma 3.1: $\mathrm{ed}=1$ produces a faithful curve with a function of degree $\le n$.
  - **Footnote 4 after (1.2).** Cites FKW Ex. 4.6 and Lemma 4.9; this was verified.
- **B. Farb, J. Wolfson, *Resolvent degree, Hilbert's 13th problem and geometry*, Enseign. Math. 65 (2019), arXiv:1803.04063 [abstract].**
  It gives $\mathrm{RD}(A_7)\le3$.
- **Farb–Kisin–Wolfson** [abstract].
  - Papers: Compositio 2021, Math. Ann. 2023, Duke 2024.
  - Their methods give essential dimension at $p$.
  - FW Remark 1.9 says these methods cannot give the $A_7$ result.
- **Z. Reichstein, *Hilbert's 13th problem for algebraic groups*, Enseign. Math. 2025, arXiv:2204.13202 [abstract].** RD is $\le5$ for connected groups.
- **Edens–Reichstein, Doc. Math. 2025, arXiv:2406.15954 [abstract].** In characteristic $p$ the conjectures can fail.
- **Upper bounds on $\mathrm{RD}(n)$ [abstract].**
  - Sutherland, arXiv:2107.08139;
  - Heberle–Sutherland, NYJM 2023, arXiv:2110.08670;
  - Gómez-Gonzáles–Sutherland–Wolfson, J. Algebra 2024.

## 3. The curve and its group

- **Conder: the strong symmetric genus of $A_7$ is 136 [abstract]** (e.g. arXiv:1310.3871).
  - $A_7$ has 4 regular maps of genus 136, which matches `triples.py`.
  - The value is rechecked by `verify_accessory60.py`.
- **D. Singerman, *Finitely maximal Fuchsian groups*, J. London Math. Soc. 1972 [unverified].** $(2,4,7)$ is not on the list of non-maximal signatures. So $\mathrm{Aut}(C)=A_7$.
- **K. Takeuchi, arithmetic triangle groups, 1977.** $(2,4,7)$ **is arithmetic**, with invariant trace field $\mathbb Q(\cos\frac{2\pi}7)$.
  - Source: Nugent–Voight, arXiv:1510.04637, §6.1.1.
  - The kernel of $\Delta(2,4,7)\to A_7$ is non-congruence.
  - So no Selberg-type bound applies, and the certified $\lambda_1>1/4$ is a fact about this cover.
  - An earlier note here wrongly said "not arithmetic". No proof depended on it.

## 4. Algebraic gonality tools

- **H. Stichtenoth, *Algebraic Function Fields and Codes*.** Thm 3.11.3 (Castelnuovo–Severi); Cor. 3.11.4.
- **R. Accola, *On the Castelnuovo–Severi inequality for Riemann surfaces*, Kodai Math. J. 29 (2006) [abstract].** Refinements and the equality case.
- **G. Martens, *The gonality of curves on a Hirzebruch surface*, Arch. Math. 67 (1996) [unverified].**
  Used in Thm 1.6(c): the only $g^1_9$'s on a smooth $(9,9)$ curve are the rulings.
- **N. Ishii, *Coverings over d-gonal curves*, Tsukuba J. Math. 16 (1992) [unread].**
- **E. Casas-Alvero, *Singularities of Plane Curves*, §3.5.** Infinitely near points and proximity (Lemma 4.5).
- **L. Gruson, C. Peskine (Halphen's bound) for space curves not on a quadric.** Used in §5.4.
- **M. Coppens, G. Martens, *Secant spaces and Clifford's theorem*, Compositio Math. 78 (1991) [unverified wording].** $\operatorname{gon}\le\mathrm{Cliff}+3$ (Theorem 7.4).
- **I. Schur (1911); the exceptional triple covers $3.A_6$, $3.A_7$** (ATLAS). $H^2(A_7,\mathbb C^\*)=\mathbb Z/6$; $3.A_7\subset SL_6(\mathbb C)$ (Ch. 7).
- **R. Laza, Z. Zheng, *Automorphisms and periods of cubic fourfolds*, Math. Z. (2022), arXiv:1905.11547 [abstract].** The 34 symplectic groups,
  six maximal. $A_7$ is realised by exactly two smooth cubic fourfolds (Ch. 7, §7.6).
- **Yang, Yu, Zhu (2023) [unverified]; K. Koike, *Cubic fourfolds with symplectic automorphisms*, arXiv:2409.08448 [read: Example 2.1].**
  The $3.A_7$-invariant cubic fourfold, with an equation over $\mathbb Q$. Our degree-60 model of the genus-136 curve lies on it.
- **D. Mumford, *On the equations defining abelian varieties I*, Invent. Math. 1 (1966).** Theta groups and the Mumford class of an invariant line bundle (Ch. 7).
- **M. F. Atiyah, R. Bott, *A Lefschetz fixed point formula for elliptic complexes II*, Ann. Math. 88 (1968).** The holomorphic Lefschetz formula (Ch. 5, 7).

## 5. Spectral gonality and certified eigenvalues

**Li–Yau.**
- **J. Hersch, C. R. Acad. Sci. Paris 270 (1970) [unverified wording].** The balancing lemma.
- **P. C. Yang, S.-T. Yau, Ann. SNS Pisa 7 (1980); P. Li, S.-T. Yau, Invent. Math. 69 (1982) [unverified wording].**
  These give $\lambda_1\mathrm{Area}\le8\pi d$ (Lemma 2.1).
- **Earlier applications.**
  - D. Abramovich, IMRN 1996 (modular curves);
  - Ellenberg–Hall–Kowalski, Duke 2012;
  - Cornelissen–Kato–Kool, arXiv:1211.2681, and Amini–Kool (graph versions).

  All of these rely on automorphic eigenvalue bounds (Selberg, Kim–Sarnak).
- **R. Bryant, *Minimal surfaces of constant curvature in $S^n$*, Trans. AMS 290 (1985).** Used in Prop. 6.2.

**Guaranteed lower bounds.**
- **C. Carstensen, J. Gedicke, *Guaranteed lower bounds for eigenvalues*, Math. Comp. 83 (2014) [read: Thm 2.1 and its proof].**
  - $\kappa^2=1/8+j_{1,1}^{-2}$.
  - The single-triangle estimate (also Carstensen–Gedicke–Rim, J. Comput. Math. 2012).
- **R. Laugesen, B. Siudeja, J. Differential Equations 249 (2010).** The Neumann Poincaré constant $h_T/j_{1,1}$ on triangles.
- **X. Liu, Appl. Math. Comput. 267 (2015); X. Liu, S. Oishi, SIAM J. Numer. Anal. 51 (2013) [unverified wording].**
  The abstract bound $\lambda_k\ge\lambda_{k,h}/(1+C_h^2\lambda_{k,h})$.
- **N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., Thm 10.3.** The componentwise Cholesky backward error.
- **S. M. Rump, BIT 46 (2006) [unverified wording].** Shift-and-Cholesky verification of positive definiteness.

**Eigenvalues of specific surfaces.**
- Strohmaier–Uski, CMP 317 (2013): Bolza $\lambda_1=3.8388872588$, used in `validate_bolza.py`.
- Jenni 1984 (Bolza).
- Fortier Bourque–Petri, arXiv:2111.14699 (Klein quartic).
- Lee, arXiv:2311.02632 (Fricke–Macbeath, genus 7).
- Kravchuk–Mazáč–Pal, arXiv:2111.12716 (bootstrap upper bounds).

**Maass forms and quasimodes.**
- Hejhal's method;
- Strömberg's vector-valued extension;
- Booker–Strömbergsson–Venkatesh 2006;
- Seymour-Howell et al., arXiv:2502.01442 (only $\Gamma_0(N)$).

**Context.**
- Magee–Naud–Puder and Hide–Magee: spectral gaps of random covers.
- Monk–Naud, ICM 2026 survey, arXiv:2601.13988: no deterministic certified computations.

## 6. Novelty search: findings

- **Prior work on these curves.** None was found on the gonality or the spectrum of the $(2,4,7)$ $A_7$-curves. Only their existence and census are known.
- **Prior certified $\lambda_1$.** It exists only in genus $\le7$, by symmetry reduction, particular solutions or the trace formula.
  - None of these uses guaranteed-lower-bound finite elements on a closed hyperbolic surface.
  - None uses sign-twisted quotients that cover all irreducibles.
- **Follow-ups.** None to Farb–Wolfson 2025 was found.
- **Checked, not relevant.**
  - OpenAI, *Ten advances in mathematics and theoretical computer science* (Aug 2026): none concerns Hilbert 13.
  - *A Differential-Algebraic Solution to Hilbert's 13th Problem* (Cambridge Engage, Jan 2026): not refereed, not evaluated.

**Re-fetching.** `curl -sL -o fw.pdf https://arxiv.org/pdf/2510.22786v1` works through the proxy, while WebFetch to arxiv.org was blocked.
Extract the text with `pdfminer.six`.
