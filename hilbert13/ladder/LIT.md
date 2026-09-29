# Literature map for the A7 ladder

Status tags: **[read]** full text read in session · **[abstract]** abstract/summary only ·
**[unverified]** cited from memory, check before relying on it.

## Project context

- **B. Farb, J. Wolfson, *Essential dimension relative to branched covers of degree at most n*,
  arXiv:2510.22786 (Oct 2025), 9 pp. [read]** — the backbone. Defines
  $\mathrm{ed}_k(G;\le n)$: the least dimension to which a $G$-variety compresses after a branched cover
  of degree $\le n$. Main Theorem 1.6: if $G$ has no proper subgroup of index $\le n$, contains
  $M\hookrightarrow\mathrm{PSL}_2(k)$ with $|M|>n$, and acts on no curve of genus $\le(n-1)^2$, then
  $\mathrm{ed}_k(G;\le n)>1$. Corollary 1.7(2): $\mathrm{ed}_{\mathbb C}(A_7;\le6)>1$ ("the general septic cannot be
  reduced to a 1-variable algebraic function even allowing an accessory sextic"). Tools used here:
  - Lemma 2.2: if $f_1,\dots,f_j$ of degree $n$ generate $k(C)$ then $g(C)\le(n-1)^2$ (induction on CS);
  - Corollary 2.3: if $|G|\nmid n$ and $G$ acts on no curve of genus $\le(n-1)^2$, no faithful $G$-curve has a
    degree-$n$ function;
  - Lemma 2.4 (norm): a dominant $H\to C$ and a degree-$n$ function on $H$ give one on $C$ (so gonality
    of a quotient is at most gonality of the cover);
  - Lemma 3.1: $\mathrm{ed}(G;\le n)=1$ produces a faithful $G$-curve with a function of degree $\le n$.
  Their $A_7$ proof uses only the Hurwitz bound (genus $\ge31$). With Conder's value 136 (below),
  Corollary 2.3 gives **gonality $\ge13$ for every $A_7$-curve**. Note that hypothesis (1) caps
  $n\le6$ for $A_7$ (because of $A_6$ of index 7), so gonality bounds past 13 serve refined variants.
- **B. Farb, J. Wolfson, *Resolvent degree, Hilbert's 13th Problem and geometry*,
  Enseign. Math. 65 (2019) 303–376, arXiv:1803.04063. [abstract]** — $\mathrm{RD}(A_7)\le3$; Hilbert's
  conjecture $\mathrm{RD}(7)=3$; no $n$ is known with $\mathrm{RD}(n)>1$.
- **Farb–Kisin–Wolfson**: *The essential dimension of congruence covers*, Compositio 157 (2021);
  *Modular functions and resolvent problems*, Math. Ann. 386 (2023); *Essential dimension via
  prismatic cohomology*, Duke 173 (2024). [abstract] — "essential dimension at $p$" methods, which FW
  show cannot give the $A_7$ result (their Remark 1.9).
- Quanta, *Mathematicians Resurrect Hilbert's 13th Problem* (Jan 2021). [abstract]

## The curve

- **Strong symmetric genus of $A_7$ is 136 (Conder).** [abstract: confirmed by search, e.g.
  arXiv:1310.3871] Matches `mingenus.py`: (2,4,7) is the smallest-genus signature generating $A_7$;
  (2,3,7), (2,4,5), (2,4,6), (3,3,4), (2,5,5) do not generate.
- **D. Singerman, *Finitely maximal Fuchsian groups*, J. London Math. Soc. (1972). [unverified]** —
  (2,4,7) is not in the list of non-maximal triangle signatures, so $\mathrm{Aut}(C)=A_7$.
- **K. Takeuchi, arithmetic triangle groups (1977). [unverified]** — (2,4,7) is not arithmetic, so no
  Selberg-type $\lambda_1\ge3/16$ or $1/4$ is available a priori.

## Castelnuovo–Severi and gonality tools

- **H. Stichtenoth, *Algebraic Function Fields and Codes*, Thm 3.11.3 (CS), Cor 3.11.4 (Riemann
  inequality). [via FW]**
- **R. Accola, *On the Castelnuovo–Severi inequality for Riemann surfaces*, Kodai Math. J. 29 (2006)
  299–317. [abstract]** — refinements and equality discussion; worth reading for rung 4.
- **G. Martens, *The gonality of curves on a Hirzebruch surface*, Arch. Math. 67 (1996).
  [unverified]** — needed for Lemma 4.0 (the only $g^1_9$'s on a smooth (9,9) curve are the rulings).
  Check the exact statement.
- **N. Ishii, *Coverings over d-gonal curves*, Tsukuba J. Math. 16 (1992); *Remarks on d-gonal
  curves*, 19 (1995). [via FW, unread]**
- Hodge-index form of CS (equality iff $\Gamma\equiv d_2F_1+d_1F_2$): standard; see e.g. notes by
  Pignatelli, *Surfaces of general type* (CRM 2015). [abstract]

## Spectral route (section 7: used in the proof of both summits)

- **J. Hersch, *Quatre propriétés isopérimétriques de membranes sphériques homogènes*, C. R. Acad. Sci.
  Paris 270 (1970). [unverified wording]** Source of the balancing lemma (a Möbius transformation that
  centres a measure on $S^2$ without atoms). 7.2 proves the degree-$d$ version in full, apart from this
  lemma.
- **P. C. Yang, S.-T. Yau, *Eigenvalues of the Laplacian of compact Riemann surfaces and minimal
  submanifolds*, Ann. Sc. Norm. Sup. Pisa 7 (1980). [unverified wording]** Proves
  $\lambda_1\mathrm{Area}\le8\pi d$ for a degree-$d$ holomorphic map to $\mathbb P^1$.
- **C. Carstensen, J. Gedicke, *Guaranteed lower bounds for eigenvalues*, Math. Comp. 83 (2014)
  2605–2629. [read: Thm 2.1 and its proof]** Two facts are used. First, Theorem 2.1:
  $\|v-I_{NC}v\|\le\kappa H|||v-I_{NC}v|||_{NC}$ with $\kappa^2=1/8+j_{1,1}^{-2}$. Second, its proof: the single-triangle
  estimate $\|f\|_{L^2(T)}\le(\max_{x\in E}|P-x|^2/8+h_T^2/j_{1,1}^2)^{1/2}|f|_{H^1(T)}$ for $\int_Ef=0$.
- **C. Carstensen, J. Gedicke, D. Rim, *Explicit error estimates for Courant, Crouzeix–Raviart and
  Raviart–Thomas finite element methods*, J. Comput. Math. 30 (2012). [via CG14]** Lemma 2.2 is the
  single-triangle estimate above.
- **R. Laugesen, B. Siudeja, *Minimizing Neumann fundamental tones of triangles: an optimal Poincaré
  inequality*, J. Differential Equations 249 (2010). [via CG14]** The $h_T/j_{1,1}$ Poincaré constant on
  triangles.
- **X. Liu, *A framework of verified eigenvalue bounds for self-adjoint differential operators*, Appl.
  Math. Comput. 267 (2015). [unverified wording]** The abstract lower bound
  $\lambda_k\ge\lambda_{k,h}/(1+C_h^2\lambda_{k,h})$. 7.6 re-proves it in the form used here.
- **A. Strohmaier, V. Uski, *An algorithm for the computation of eigenvalues, spectral zeta functions
  and zeta-determinants on hyperbolic surfaces*, Comm. Math. Phys. 317 (2013). [unverified wording]**
  Gives $\lambda_1(\text{Bolza})=3.8388872588\ldots$, used as the external check in `validate_bolza.py`.
- **N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM 2002, Thm 10.3.**
  The componentwise backward error of Cholesky, $|\Delta A|\le\gamma_{n+1}|\hat R^T||\hat R|$. For sparse
  factors, $n$ becomes the largest number of terms in an inner product, which is at most the maximal row
  count.
- **S. M. Rump, *Verification of positive definiteness*, BIT 46 (2006). [unverified wording]** Uses the
  same shift-and-Cholesky idea with a trace bound. 7.7 uses the sharper $\|L\|_1\|L\|_\infty$ bound instead.

- **P. Li, S.-T. Yau, *A new conformal invariant and its applications to the Willmore conjecture and
  the first eigenvalue of compact surfaces*, Invent. Math. 69 (1982). [unverified wording]** —
  $\lambda_1\,\mathrm{Area}\le8\pi d$ for a degree-$d$ conformal map to $S^2$. (Also Yang–Yau 1980.)
- **D. Abramovich, *A linear lower bound on the gonality of modular curves*, IMRN 1996. [abstract]** —
  Li–Yau plus an eigenvalue bound gives gonality $\ge\lambda_1\cdot\mathrm{index}/24$; the template for 6.
- **G. Cornelissen, F. Kato, J. Kool, *A combinatorial Li–Yau inequality and rational points on
  curves*, arXiv:1211.2681, Math. Ann. [abstract]** — graph version; possible route to a combinatorial
  certificate using the (2,4,7) tiling.
- **X. Liu, S. Oishi, *Verified eigenvalue evaluation for the Laplacian over polygonal domains of
  arbitrary shape*, SIAM J. Numer. Anal. 51 (2013). [unverified]** The origin of the Crouzeix–Raviart
  lower-bound method used in 7.6.
- **Magee–Naud–Puder; Hide–Magee (2023)**: spectral gaps of random covers tend to $1/4$.
  [unverified] Context only.

## Checked and not relevant

- OpenAI, *Ten advances in mathematics and theoretical computer science* (Aug 2026) and the repo
  `openai/ten-proofs`: none of the ten results concerns Hilbert 13 (checked the README and grepped the
  Lean sources).
- *A Differential-Algebraic Solution to Hilbert's 13th Problem*, Cambridge Engage preprint
  (Jan 2026): not peer-reviewed; not evaluated.
- No public source found for a GPT result on an A7 subcase; presumably it is the user's own
  collaboration.

## Network notes (for re-fetching)

`WebFetch` was blocked for arxiv.org, openai.com, simonsfoundation.org and simonwillison.net, but
`curl` to arxiv.org works:
`curl -sL -o fw.pdf https://arxiv.org/pdf/2510.22786v1`, then extract text with `pdfminer.six`.
