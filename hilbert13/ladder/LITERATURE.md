# Literature

Searches were made on 29 September 2026 (novelty) and 2 October 2026 (state of the art), over the web, arXiv abstracts and EMS/AMS/Springer. MathSciNet and zbMATH were not searched.

Tags:
- **[read]**: the full text was read.
- **[abstract]**: the abstract or a summary only.
- **[unverified]**: cited from memory; check before relying on it.

## 1. Where the project stands

| question | published | this project |
|---|---|---|
| Hilbert 13 proper: $\mathrm{RD}(n)>1$ for some $n$? | open for every $n$ | not touched |
| $\mathrm{ed}_{\mathbb C}(A_7;\le n)>1$ | $n\le6$ (Farb–Wolfson 2025) | $n\le59$, sharp: $a(A_7)=60$ (GPT, verified; Ch. 5) |
| the same with connected full monodromy | — | exactly $\mu(A_7)=90$ (Ch. 5) |
| $a(A_6)$, $a(L_2(7))$, $a(A_5)$ | $\ge5$, $\ge3$ (FW); $a(A_5)=2$ (Klein) | 12, 4, 2 (GPT Round 5; verified) |
| gonality of the genus-136 $(2,4,7)$ curves | $[13,56]$ (FW Lemma 2.2 with Conder's genus; a quotient) | $[25,42]$ (Ch. 2–3, 7) |
| $\gamma(A_7)$, the least gonality of a faithful $A_7$-curve | $\ge13$ | $[25,42]$ (Ch. 4, 7) |
| certified $\lambda_1$ of a closed hyperbolic surface | genus $\le7$ | genus 136 and 24 further $A_7$-curves (Ch. 2) |

**Assessment.**
- **Accessory Hilbert 13 for $A_7$.** It is settled sharply: 60, and 90 with full monodromy. The published bound was 6.
- **Spectral methods.** Going past Li–Yau by representation theory appears new (Lemma 3.3: the degree is a cubic form killed by $(\wedge^3E_1)^G=0$). The certification follows Booker–Strömbergsson–Venkatesh; our additions are exact equivariance by a partition of unity and the Jacobian-flux bound.
- **Schur-twisted models** (Ch. 7). We found no prior use of Schur-multiplier twists, with an exact degree formula from the centrally extended triangle group, to build low-degree models of symmetric curves. The resulting link is also new: the genus-136 curves are degree-60 curves cut out of the Laza–Zheng $A_7$-cubic fourfold by quartics.
- **Elliptic subcovers** (§7.9). Optimal maps $C\to E$ correspond to primitive rank-2 sub-Hodge lattices, with degree equal to the restricted cup product; this is classical (Kani; Birkenhake–Lange). New is only the exact $H^1(C,\mathbb Z)$ from the dessin and the lattice search.
- **Caveats.** Nothing here is refereed. The value $a=60$ uses the convention in which the monodromy may drop after the accessory; the bound $\le59$ holds in either convention.

## 2. Hilbert 13 and essential dimension

- **B. Farb, J. Wolfson, *Essential dimension relative to branched covers of degree at most n*, arXiv:2510.22786 (2025) [read].** The backbone.
  - **Theorem 1.6.** $\mathrm{ed}_k(G;\le n)>1$ if three conditions hold:
    - no proper subgroup of $G$ has index $\le n$;
    - some $M\hookrightarrow\mathrm{PSL}_2(k)$ has $|M|>n$;
    - $G$ acts on no curve of genus $\le(n-1)^2$.
  - **Cor. 1.7.** $\mathrm{ed}_{\mathbb C}(A_7;\le6)>1$.
  - **Tools.**
    - Lemma 2.2: if functions of degree $n$ generate $k(C)$, then $g\le(n-1)^2$.
    - Lemma 2.4: the gonality of a quotient is at most that of the cover.
    - Lemma 3.1: $\mathrm{ed}=1$ gives a faithful curve with a function of degree $\le n$.
  - Footnote 4 cites FKW Ex. 4.6 and Lemma 4.9; checked.
- **B. Farb, J. Wolfson, Enseign. Math. 65 (2019), arXiv:1803.04063 [abstract].** $\mathrm{RD}(A_7)\le3$.
- **Farb–Kisin–Wolfson**, Compositio 2021, Math. Ann. 2023, Duke 2024 [abstract]. Essential dimension at $p$. FW Remark 1.9 says these methods cannot give the $A_7$ result.
- **Z. Reichstein, Enseign. Math. 2025, arXiv:2204.13202 [abstract].** $\mathrm{RD}\le5$ for connected groups.
- **Edens–Reichstein, Doc. Math. 2025, arXiv:2406.15954 [abstract].** The conjectures can fail in characteristic $p$.
- **Upper bounds on $\mathrm{RD}(n)$ [abstract].**
  - Sutherland, arXiv:2107.08139;
  - Heberle–Sutherland, NYJM 2023;
  - Gómez-Gonzáles–Sutherland–Wolfson, J. Algebra 2024.

## 3. The curve and its group

- **Conder: the strong symmetric genus of $A_7$ is 136 [abstract]** (arXiv:1310.3871). There are 4 regular maps, matching `curve_checks.py`; the genus is rechecked by `verify_accessory60.py`.
- **D. Singerman, *Finitely maximal Fuchsian groups*, J. London Math. Soc. 1972 [unverified].** $(2,4,7)$ is finitely maximal, so $\mathrm{Aut}(C)=A_7$.
- **K. Takeuchi (1977), via Nugent–Voight, arXiv:1510.04637, §6.1.1.**
  - $\Delta(2,4,7)$ is arithmetic, with trace field $\mathbb Q(\cos\frac{2\pi}7)$.
  - $\ker(\Delta\to A_7)$ is non-congruence, so no Selberg-type bound applies.
- **R. Laza, Z. Zheng, Math. Z. 2022, arXiv:1905.11547 [abstract].** Symplectic automorphism groups of cubic fourfolds: 34 groups, six maximal. $A_7$ acts on exactly two smooth cubic fourfolds.
- **Yang–Yu–Zhu (2023) [unverified]; K. Koike, arXiv:2409.08448 [read: Ex. 2.1].** The $3.A_7$-invariant cubic fourfold over $\mathbb Q$ (§7.6).

## 4. Algebraic tools

- **H. Stichtenoth, *Algebraic Function Fields and Codes*, Thm 3.11.3.** Castelnuovo–Severi.
- **R. Accola, Kodai Math. J. 29 (2006) [abstract].** The equality case of Castelnuovo–Severi.
- **G. Martens, Arch. Math. 67 (1996) [unverified].** The only $g^1_9$'s on a smooth $(9,9)$ curve are the rulings (Thm 1.5(c)).
- **E. Casas-Alvero, *Singularities of Plane Curves*, §3.5.** Proximity (Lemma 4.5).
- **L. Gruson, C. Peskine.** Halphen's bound (§5.4).
- **M. Coppens, G. Martens, Compositio 78 (1991) [unverified wording].** $\operatorname{gon}\le\mathrm{Cliff}+3$ (Theorem 7.4).
- **B. Hassett, Y. Tschinkel, *Torsors and stable equivariant birational geometry*, Nagoya Math. J. 250 (2023) [unverified; cited by GPT].** Amitsur subgroups, the unit sequence on Picard-trivial opens, the no-name lemma (Theorem 5.14).
- **J. Milnor, *On the 3-dimensional Brieskorn manifolds* (1975) [unverified].** The centrally extended triangle groups $\langle c_i\mid c_1^p=c_2^q=c_3^r=c_1c_2c_3\rangle$ (Prop. 7.1).
- **Stacks Project, §11.8 [unverified; cited by GPT].** Splitting fields of central simple algebras (Corollary 5.18).
- **I. Schur (1911); ATLAS.** $H^2(A_7,\mathbb C^\*)=\mathbb Z/6$ and $3.A_7\subset SL_6$.
- **R. A. Wilson et al., *ATLAS of Group Representations* v3 [read: the matrices].** `3A7G1-Ar6B0`, the $\mathbf 6$ of $3.A_7$ over $\mathbb Z[\omega]$; checked in `verify_schur_exact.py`.
- **M. Green, R. Lazarsfeld, Invent. Math. 83 (1986) [unverified wording].** Nonvanishing: $\operatorname{Cliff}\le p$ implies $K_{p,2}\ne0$. With Noether's and Petri's theorems (ACGH Ch. III), the tests of Theorem 7.13.
- **D. Mumford, Invent. Math. 1 (1966).** Theta groups and Mumford classes.
- **M. F. Atiyah, R. Bott, Ann. Math. 88 (1968).** Holomorphic Lefschetz (Ch. 5, 7).
- **E. Kani, J. reine angew. Math. 485 (1997) [unverified]; Birkenhake–Lange, *Complex Abelian Varieties* [unverified section].** Elliptic subcovers via sublattices (§7.9).
- **Cited only in GPT's documents** (`gpt/README.md` §2; not used here). Petrakiev, arXiv:math/0604517 (refined Castelnuovo bounds); Harui, arXiv:1306.5842 (automorphisms of plane curves); Vinokurov, arXiv:2502.03756 (equivariant eigenvalue optimisation); Karpenko–Merkurjev, Invent. Math. 172 (2008) (essential dimension at $p$); Niu–Ulrich, arXiv:1404.5092 (conductor and duality).
- **For §7.11 (not used) [unverified].** M. Baker, specialisation of linear systems from curves to graphs (2008); M. Raynaud (1999) and S. Wewers (2003) on the stable reduction of three-point covers when $p\,\|\,|G|$.

## 5. Spectral gonality and certified eigenvalues

- **Li–Yau.**
  - Hersch, C. R. Acad. Sci. 270 (1970), the balancing lemma.
  - Yang–Yau, Ann. SNS Pisa 7 (1980); Li–Yau, Invent. Math. 69 (1982) [unverified wording]: $\lambda_1\mathrm{Area}\le8\pi d$.
- **Earlier applications**, all using automorphic eigenvalue bounds:
  - Abramovich, IMRN 1996;
  - Ellenberg–Hall–Kowalski, Duke 2012;
  - Cornelissen–Kato–Kool, arXiv:1211.2681, and Amini–Kool (graphs).
- **R. Bryant, Trans. AMS 290 (1985).** Prop. 6.2.
- **Guaranteed lower bounds.**
  - Carstensen–Gedicke, Math. Comp. 83 (2014) [read: Thm 2.1]: $\kappa^2=1/8+j_{1,1}^{-2}$.
  - Liu, Appl. Math. Comput. 267 (2015); Liu–Oishi, SINUM 51 (2013) [unverified wording]: $\lambda_k\ge\lambda_{k,h}/(1+C_h^2\lambda_{k,h})$.
  - Laugesen–Siudeja, JDE 249 (2010).
  - Higham, *Accuracy and Stability*, Thm 10.3.
  - Rump, BIT 46 (2006) [unverified wording].
- **Specific surfaces.**
  - Strohmaier–Uski, CMP 317 (2013): Bolza, $\lambda_1=3.8388872588$.
  - Jenni 1984 (Bolza).
  - Fortier Bourque–Petri, arXiv:2111.14699 (Klein quartic).
  - Lee, arXiv:2311.02632 (genus 7).
  - Kravchuk–Mazáč–Pal, arXiv:2111.12716.
- **Maass forms.**
  - Hejhal's method; Strömberg's vector-valued version.
  - Booker–Strömbergsson–Venkatesh 2006.
  - Seymour-Howell et al., arXiv:2502.01442.
- **Context.**
  - Magee–Naud–Puder; Hide–Magee.
  - Monk–Naud, ICM 2026 survey, arXiv:2601.13988: no deterministic certified computations.

## 6. Novelty search

- **These curves.** No prior work was found on their gonality or spectrum; only existence and census are known.
- **Certified $\lambda_1$.** Prior results exist only in genus $\le7$. None uses guaranteed-lower-bound finite elements on a closed hyperbolic surface, or sign-twisted quotients covering all irreducibles.
- **Follow-ups.** None to Farb–Wolfson 2025 was found.
- **Not relevant.**
  - OpenAI, *Ten advances…* (Aug 2026).
  - *A Differential-Algebraic Solution to Hilbert's 13th Problem* (Jan 2026; not refereed).

**Fetching.** Use `curl -sL -o fw.pdf https://arxiv.org/pdf/2510.22786v1`, because WebFetch to arxiv.org is blocked. Extract text with `pypdf`.
