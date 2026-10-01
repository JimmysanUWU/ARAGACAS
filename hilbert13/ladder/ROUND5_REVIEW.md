# Review of GPT's Round-5 answers (*A7: new bounds and harmonic Hersch*, 1 October 2026)

**Source.** `gpt/A7_Round5_Advances.pdf`, by GPT (Sol 6.1, fast mode). It answers `QUESTIONS_FOR_GPT.md` Round 5 at
`1ef3f6f`.

**Verdict.** Every claim checks out: it is either re-derived on paper (**[P✓]**) or recomputed exactly (**[X]**). The
document found a real bug in our code. Its harmonic Hersch inequality supersedes our topological Hersch. Our
equivariant Riemann–Roch computation (§3) settles the degree-90 test that the document leaves open, which gives
$\mu(A_7)\le90$.

## 1. Claim-by-claim

| § | claim | status |
|---|---|---|
| 1 | Prop. 1.1: for every torsor (also when disconnected after base change), the subgroup upper bound $[G:H]\mu(H)$ holds | [P✓]. Twist $(C,L)$ by an $H$-reduction over a degree-$\le[G:H]$ point, then take a point of a general reduced section's zero set in the free locus. |
| 1 | Citation: Farb–Wolfson, footnote 4 after (1.2) in Problem 1.4, cites FKW Ex. 4.6 and Lemma 4.9 | **verified** against arXiv:2510.22786 (HTML). Vinokurov arXiv:2502.03756 also exists as cited. |
| 2 | Thm 2.2: $\mu(A_7)\ge72$ | [P✓] + [X]. $n\in\{60,66\}$. The only lattice signatures are $(2,6,7)$ and $(3,4,7)$; $\pi(60,8)=220$ forces the section representation to be exactly the standard 6. Its $A_6$-fixed vector has a degree-60 zero divisor, but $A_6$-orbits have size $\ge120$ (resp. 90). |
| 2 | Thm 2.3: the power-sum curve $\{\sum x_j^k=0,\ k\le5\}\subset\mathbb P^6$ is smooth, connected, of genus 481, with linearised $\mathcal O(1)$ of degree 120 | [P✓]. Moment/Vandermonde smoothness, Koszul connectedness, adjunction $K=\mathcal O(8)$. Consistency: $A_7$-signature $(3,7,7)$ ($g=481$) with $N=120$. |
| 2 | Prop. 2.4: $\mu(A_7)\in\{72,84,90,96,108,120\}$ | [X]. 78 and 102 have no lattice signature. For 114: $(3,4,5,7)$, $\pi(114,6)=1221<1354$, and $\prod x_j$ is an invariant section of degree 798, which is not a sum of orbit sizes $\{360,504,630,840,2520\}$. |
| 2 | On a $(2,4,7)$ curve, $\mathrm{Pic}^{A_7}_{\rm lin}=\langle D_2,D_4,D_7\mid2D_2=4D_4=7D_7\rangle\cong\mathbb Z\oplus\mathbb Z/2$. The degree-90 classes are $B$ and $B+T$, and $\mu_C\in\{90,180\}$ | [P✓] (Smith form $1,2$), [X]. **Decided in §3.** |
| 3 | $a(A_5)=2$, $a(A_6)=12$, $a(L_2(7))=4$ | [P✓]. $A_6$: $6\mid n$, $\pi(6,4)=2$. There are no subgroups of order 40 or 45 (normal Sylow 5, normaliser of order 10). The order-36 subgroup is not a Möbius group. |
| 3 | $\mu(A_n)\le(n-2)!$, $a(A_n)\le\min((n-2)!,\,n!/84)$, and the $L_2(q)$ torus-normaliser bounds | [P✓]. Odd dihedral groups lift projectively faithfully to $GL_2$. |
| 4 | Prop. 4.1: correspondence over $B$ with a fixed point and $\mathrm{Alb}(B)=0$ | [P✓] |
| 4 | Prop. 4.2: a quadratic accessory $t^2=q(v)$ keeps full $A_7$ but its smooth rational model has no fixed point | [P✓]. The only invariant line in $6\oplus1$ is the $t$-axis, which is not on the quadric. |
| 5 | **Harmonic Hersch (HH)**: $\lambda_1(A-s)\le4\kappa\sqrt\varepsilon(A-s)+\sqrt{\Lambda(A-s)s(\varepsilon+\lambda_1s)}$ | [P✓]. For holomorphic $x$, $x\times dx={*dx}$ gives $\Theta(x,x,y)=\tfrac12\langle\nabla x,\nabla y\rangle=\tfrac12\lambda_1\|y\|^2$, and $\Theta(x,x,y)=2\Theta(z,y,y)+\Theta(z,y,z)$. This is strictly stronger than our (TH): $4\kappa$ against $6\kappa$, with no $E(z)/2$ term. |
| 5 | Cor. 5.2: $B^*\le3\cdot10^{-4}$ and $\Gamma\le2$ suffice, using only $\kappa^2\le B^*/3$ (no quartic optimisation) | [P✓] + [X] (rational assertions rerun). At the computed constants the allowed inflation is **1.56** (ours was 1.035). |
| 5 | Thm 5.3: trial-space certificate ($\lambda_h\le0.351$, $\rho\le0.005$, $B_h\le2\cdot10^{-4}$, $\Gamma_h\le1.9$) | [P✓], every step. Equal-rank $\|P-P_h\|\le\rho\sqrt b/(b-r)$; $\nu=b-(b-a)\delta^2$; $|c|\le\rho\sqrt{tE_z}$; the trial-space identity $\lambda_ht+c=4\Theta(z,y,y)+2\Theta(z,y,z)$; resolvent monotonicity on the isotypes $10,\overline{10},15,21,35$. [X] the rational assertions (exact root $R=0.206831$). |
| 5 | Complementary-energy upper bound (5.4) for the resolvent form | [P✓] (Galerkin identity plus coercivity on the isotype) |
| 5 | **Bug in `topo_hersch.py`**: the frame matrix should be $L^TSL$, not $LSL^T$ | **confirmed.** Old trace error $3.3\times10^{-5}$, new $10^{-18}$. $\Gamma$: $n=4$: $1.70065\to1.67134$ (exactly as reported); $n=8$: $1.6918\to1.6346$. Fixed. The error was conservative: all earlier conclusions stand. |

## 2. What changes for the gonality track

- **Route.** Harmonic Hersch replaces topological Hersch (`FRAMEWORK_TOPOLOGICAL_HERSCH.md`). The crude, rigorous
  $\kappa^2\le B^*/3$ now suffices, so the "certify the global maximum on $S^{41}$" item disappears.
- **Tolerance.** With the certified eigenvalue inputs, $\sqrt{B^*}$ and $\Gamma$ may exceed their computed values by
  56%, instead of 3.5%.
- **Remaining work.** Certify the four bounds of Thm 5.3 on an exact copy of $14_{(5,2)}$ built from $C/S_5$, by
  equilibrated fluxes and enclosed metric coefficients. This is concrete and well within reach of the existing
  ball-arithmetic machinery.

## 3. New: the degree-90 test is decided, so $\mu(A_7)\le90$ (`equivariant_rr.py`, `equivariant_rr_general.py`)

**Method.** Holomorphic Lefschetz (Atiyah–Bott) on the $(2,4,7)$ curves gives the virtual representation
$\chi_{A_7}(\mathcal O(D))=[H^0]-[H^1]$ for invariant $D=n_2D_2+n_4D_4+n_7D_7$. For $g\ne1$,
$$\sum_q(-1)^q\,\mathrm{tr}(g\mid H^q)=\sum_{p\in\mathrm{Fix}(g)}\frac{a_p^{k_p}}{1-a_p^{-1}} ,$$
where $a_p$ is the rotation of $g$ at $p$ and $k_p=\mathrm{mult}_pD$. Only the classes of the monodromy generators enter.

**Validation.**
- $\chi(K)$ gives $H^0(K)=10+2\cdot\overline{10}+15+21+2\cdot35$ (dimension 136), matching Chevalley–Weil
  ($A^{10}\times E_1^{15}\times E_2^{21}\times S^{35}$).
- The two conventions that are not complex conjugates of this one give non-integral multiplicities, so they are
  rejected.

**Result**, the same on all four $(2,4,7)$ classes:

| class | degree | $\chi_{A_7}$ |
|---|---|---|
| $B=2D_7-D_4$ | 90 | $-\overline{10}-35$ (undecided) |
| $B+T=D_2-3D_4+2D_7$ | 90 | $-6+\mathbf{10}-14_a-14_b-21$ |
| $2B$ | 180 | $6-\overline{10}+14_a+14_b+21$ |

So $H^0(B+T)$ contains a copy of $10$ (or of $\overline{10}$, depending on orientation), and $h^0(B+T)\ge10$. The fixed part
is an invariant divisor of degree $\le90<360$, so it is zero. Hence **$B+T$ is a linearised moving bundle of degree 90 on
every $(2,4,7)$ curve**: $\mu_C(A_7)=90$, and
$$72\le\mu(A_7)\le90,\qquad \mu(A_7)\in\{72,84,90\}.$$
(By minimality, the complete series of $B+T$ is birational onto its image.)

**72 and 84 stay open.**
- The same test on every rigid signature carrying a lattice degree 72 or 84, namely $(2,5,7)$, $(3,5,7)$, $(4,5,7)$, $(5,5,7)$
  and $(3,4,5)$, $(3,5,6)$, $(4,5,6)$, $(5,5,6)$, $(5,6,6)$, $(5,6,7)$, gives $\chi$ with **no positive part** for every curve and every class.
- This is inconclusive, not a refutation.
- The four-point signatures at 84 ($(2,2,3,5)$, $(2,2,5,6)$, $(2,3,3,5)$) were not tested.
