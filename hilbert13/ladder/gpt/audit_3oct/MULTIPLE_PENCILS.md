# Multiple pencils: a mathematical framework

This proof record preserves the general geometric framework re-derived during the 3 October audit. It is independent of the spectral arithmetic. Classical inputs are adjunction, the Hodge index theorem, normalization/conductor comparison, uniform position and Petrakiev's Theorem 2.16(a). It remains subject to independent mathematical and priority review.

All curves are integral over an algebraically closed field of characteristic zero; their genus is the genus of the normalization. A pencil is a morphism from the smooth normalization to P1. Two pencils are distinct when their rational function subfields differ.

## 1. A genus bound on a singular multigraded surface

Let S be an integral hypersurface of multidegree (a,b,c), with a,b,c>0, in (P1)^3. Let B be an integral curve in S not contained in Sing(S). Write v=(m1,m2,m3)^T for its degrees over the three factors, and set

$$M=\begin{pmatrix}0&c&b\\c&0&a\\b&a&0\end{pmatrix}.$$

Then

$$g(B)\le 1+\frac12\left(v^TM^{-1}v+(a-2)m_1+(b-2)m_2+(c-2)m_3\right).$$

**Proof.** Normalize S and take a minimal resolution f:S'->S. Let H_i be the pullback of the ith coordinate divisor, and B' the strict transform of B. The intersection matrix of the H_i is M. Their sum is nef and big. If L=sum_i (M^{-1}v)_i H_i, then B'-L is orthogonal to each H_i, hence to their sum. Hodge index gives B'^2<=L^2=v^TM^{-1}v.

Adjunction for the hypersurface gives the Cartier dualizing class K_S=(a-2,b-2,c-2)|S. The canonical comparison on S' has the form K_S'=f^*K_S+Z. Its nonexceptional coefficients are nonpositive, by the conductor formula. The exceptional intersection matrix is negative definite with nonnegative off-diagonal entries. On a minimal resolution, K_S'.E>=0 for every exceptional irreducible E: adjunction gives 2p_a(E)-2-E^2, and a rational (-1)-curve is absent. Since the conductor contribution intersects E nonpositively, the exceptional coefficients are also nonpositive. Thus Z<=0.

B' is contained in neither conductor nor exceptional support, by the hypothesis on B. Consequently K_S'.B'<=f^*K_S.B'. Arithmetic adjunction on the smooth S', followed by normalization of B', gives 2g(B)-2<=B'^2+K_S'.B'. Substitution proves the formula. □

For equal degrees m the relevant values are

| Surface type | Genus upper bound |
|---|---|
| (2,1,1) | $m^2/2-m+1$ |
| (2,2,1) | $7m^2/16-m/2+1$ |
| (2,2,2) | $3m^2/8+1$ |

The hypothesis B not contained in Sing(S) must not be omitted.

## 2. The graph case, including a curve in the singular locus

An integral surface of type (a,b,1) is the graph of a rational map from P1xP1 given by two forms of bidegree (a,b). Resolve its complete base cluster. If its successive pencil multiplicities are n_i, then sum n_i^2=2ab.

Suppose B has all three coordinate degrees m, and its projection to the first two factors is birational. Let r_i be the multiplicities of the projected curve at the same cluster. Its degree under the third map is

$$m=(a+b)m-\sum_i n_i r_i,$$

so T=sum n_i r_i=(a+b-1)m. Successive blowups and arithmetic adjunction give

$$g(B)\le(m-1)^2-\frac12\sum_i(r_i^2-r_i).$$

As n_i>=1, sum r_i<=T. Cauchy–Schwarz gives sum r_i^2>=T^2/(2ab). Hence

$$g(B)\le(m-1)^2-\frac{T^2}{4ab}+\frac T2.$$

For (a,b)=(2,2) this is 7m^2/16-m/2+1.

If B lies in the singular locus of a (2,2,2) hypersurface, a nonzero first partial derivative gives a hypersurface of type (1,2,2), up to permutation, containing B. Choose the irreducible component containing B. A component with a zero coordinate degree is impossible for m>=8 and pairwise birational projections: its two-factor curve would have to have bidegree (m,m), exceeding (2,2). The remaining component is a graph case as above. Thus the bound for (2,2,1) also handles the singular-locus exception. Factors and zero-coordinate components must be addressed, rather than assuming a derivative is smooth.

## 3. A sharp bound for three pencils

**Theorem.** Let C carry three pencils of the same degree m>=8. Assume every pair gives a birational product map to P1xP1, and the eight products obtained by choosing one section from each pencil are linearly independent. Then

$$g(C)\le\left\lfloor\frac{m^2}{2}-m+1\right\rfloor.$$

The bound is sharp for every even m>=8.

**Proof.** Put A(m)=m^2/2-m+1 and suppose g>A(m). The product map is birational onto an integral curve B in (P1)^3. Its Segre degree is 3m. The eight independent products give a nondegenerate model in P7.

First prove that the product line bundle has exactly eight sections. If there were a ninth, add it to get a birational nondegenerate degree-3m model in P8. Petrakiev, Theorem 2.16(a), gives a surface of degree <=8 when the arithmetic genus exceeds pi_2(3m,8); here arithmetic genus is at least g, and 3m>=24>19. The required numerical comparison is explicit:

$$\lambda=\left\lfloor\frac{3m-1}{9}\right\rfloor,\quad
\epsilon=3m-1-9\lambda,\quad
\mu=\max\left(0,\left\lfloor\frac{\epsilon-4}{2}\right\rfloor\right),$$

$$\pi_2(3m,8)=9\binom{\lambda}{2}+\lambda(\epsilon+2)+\mu\le A(m).$$

For m=3s+r, r=0,1,2, the differences A-pi_2 are respectively s/2, (s+1)/2, and s/2+1.

Project this surface to P7. Its image remains a surface: a curve image of degree <=8 could not contain the degree-3m curve. Its degree is <=8. Any Segre quadric that did not vanish on it would cut a divisor of degree <=16, too small to contain B. Therefore the surface lies in the Segre threefold.

A divisor of type (a,b,c) on that threefold has Segre degree 2(a+b+c). For degree <=8, positive possibilities are (1,1,1) and permutations of (2,1,1). The former contradicts independence of the eight products. For the latter, Section 1 gives g<=A. A zero-coordinate divisor is excluded by pairwise birationality: its two-factor image would have bidegree (m,m), so its degree would be at least 4m. This contradiction proves h0(product)=8.

The P7 model is consequently linearly normal. Let Gamma be a general hyperplane section, a uniform-position set of 3m points in P6. It has h_Gamma(1)=7. If h_Gamma(2)<=18, linear normality gives h_B(2)=8+h_Gamma(2)<=26. The Segre threefold has 27 independent sections in multidegree (2,2,2), all induced by quadrics in P7. Thus B lies on a (2,2,2) hypersurface.

If B is outside its singular locus, Section 1 gives g<=3m^2/8+1<=A for m>=8. If B is inside, Section 2 gives g<=7m^2/16-m/2+1<=A. Reducible factors with zero coordinate degree are excluded as above. Hence h_Gamma(2)>=19.

Uniform-position subadditivity now gives, for j>=0 and j>=1 respectively,

$$h_\Gamma(2j+1)\ge\min(3m,18j+7),\qquad
h_\Gamma(2j)\ge\min(3m,18j+1).$$

The arithmetic genus is at most sum_{i>=1}(3m-h_Gamma(i)), so this also bounds g. Write m=6q+r, 0<=r<6. Summing the displayed deficits gives the following exact polynomials:

| r | Deficit sum |
|---|---|
| 0 | $18q^2-8q+1$ |
| 1 | $18q^2-2q$ |
| 2 | $18q^2+4q$ |
| 3 | $18q^2+10q+2$ |
| 4 | $18q^2+16q+5$ |
| 5 | $18q^2+22q+8$ |

Each is <=floor(A(m)); indeed the difference is 2q except for r=2, where it is 2q+1. Since m>=8, q>=1. This contradicts g>A and proves the theorem.

For sharpness, take a smooth (2,1,1) surface S. It is a degree-4 del Pezzo surface, with -K_S=(0,1,1)|S. A general smooth curve in |-kK_S|, k>=4, has coordinate degrees m=2k and genus 2k^2-2k+1=A(2k). No (1,1,1) section vanishes identically on it: on S its residual class has negative intersection with -K_S. The two birational surface projections restrict birationally; the degree-2 surface projection also restricts birationally on a general curve not invariant under the deck involution. Such curves exist because both invariant and anti-invariant anticanonical section summands are nonzero. These curves satisfy all hypotheses and attain equality. □

The P8 step is essential. Eight chosen products do not by themselves imply completeness or linear normality. The (2,2,2) singular-locus case is also essential.

Primary input: Ivan Petrakiev, *A step in Castelnuovo theory via Gröbner bases*, [arXiv:math/0604517](https://arxiv.org/pdf/math/0604517), Theorem 2.16(a) and the definition of pi_alpha in Section 2.3. The theorem concerns arithmetic genus and does not require linear normality. The remaining construction and surface/singularity comparison are the present re-derivation. Priority of the resulting three-pencil statement is not exhaustively established.

## 4. A field-generation bound

If rational functions of degrees <=d jointly generate a curve's function field, then its genus is <=(d-1)^2.

To prove this, induct on the degree budget. Choose a maximal proper compositum E of a subset of the functions, of index e. Its generators have degree <=d/e, hence g(E)<=(d/e-1)^2 by induction. Adding one further rational field generates the whole field. Castelnuovo–Severi gives

$$g\le e(d/e-1)^2+(e-1)(d-1)
=\frac{d^2}{e}+ed-3d+1.$$

This convex expression is <=(d-1)^2 on 2<=e<=d. A one-generator field is rational. Floors only strengthen the estimate.

## 5. Application to A7 in genus at least 266

**Theorem.** A faithful A7-curve of genus >=266 has gonality >=25. This is an algebraic theorem, independent of the spectral trust base.

Use two finite group inputs from the working chapters: minimum faithful genus 136, and every pencil of degree <=24 has at least 35 conjugate pencil fields (Lemma 6.1). These use the checked generating signatures and subgroup classification.

First, every faithful A7-curve has gonality >=13. Otherwise take the compositum of a degree-<=12 pencil's orbit. Its quotient curve is still faithful: the kernel is normal in simple A7 and its order is bounded by the quotient degree <=12. Section 4 would give genus <=121, contradicting minimum faithful genus 136.

Suppose now gonality d<=24 on C of genus >=266. The orbit compositum is the whole field. If it were proper of index e>=2, its faithful quotient would have a pencil of degree d/e<=12, impossible.

Every pair of distinct conjugate pencils also generates the whole field. If a pair lay in a proper compositum, enlarge it to a maximal proper compositum of orbit fields, of index e. Then e is a proper divisor of d, with 2<=e<d; e=d would make the two pencil fields equal. The preceding Castelnuovo–Severi estimate gives

$$g(C)\le d^2/e+ed-3d+1\le265.$$

For all divisors and d<=24 the last inequality is finite arithmetic; its maximum 265 occurs at (d,e)=(24,2) or (24,12). Thus all pairs are birational. Ordinary Castelnuovo–Severi then forces d>=18.

If some third conjugate gave eight independent products, Section 3 would give g<=A(d)<=A(24)=265. Therefore every third is dependent. Relative to a fixed pair, each is a complete (1,1) pencil on the quadric model. Cancelling a common component would give one of the original two pencil fields, so it has a complete length-two base cluster: two proper points, or a proper point and a first infinitely near point. Distinct pencil fields have distinct clusters.

There are at least 33 such clusters. Regard each as an edge whose vertices carry the quadric curve's multiplicities. For one edge, those multiplicities sum to d, since its moving (1,1) series has degree d rather than 2d. Its genus cost is binom(r,2)+binom(d-r,2). Two disjoint clusters would cost more than the available defect (d-1)^2-g for every 18<=d<=24, g>=266. At d=24, twice the minimum edge cost is 264, while the defect is at most 263; smaller d give a stronger comparison.

Thus the edges pairwise intersect. A family of at least four pairwise-intersecting distinct edges is a star, since the only alternative is a triangle. Let the centre multiplicity be r. The 33 distinct leaves contribute

$$\binom r2+33\binom{d-r}{2}\le(d-1)^2-g.$$

The finite integer check leaves only d=24, r=23. An infinitely near centre would require its fixed proper ancestor as the other point of every length-two cluster, contradicting distinctness of at least 33 edges. Hence the centre is proper.

Blow up this point on the quadric and use (1,1) forms through it. They give a birational plane model of degree 48-23=25. A singular point of multiplicity >=2 on that model would yield a pencil of degree <=23 by projection, contradicting minimal gonality 24. Hence the plane model is smooth, of genus (25-1)(25-2)/2=276.

Finally every faithful A7-curve has genus congruent to 1 modulo 3. In Riemann–Hurwitz each cyclic inertia order is among 2,3,4,5,6,7, and 1260/e is divisible by 3 for each. Genus 276 violates that congruence. This contradiction proves the theorem. □

The global lower bound gamma(A7)>=25 combines this result with the computer-assisted low-genus cases; it is not purely algebraic in all genera.

## 6. Repeated finite arithmetic, without a workspace dependency

The inequalities in Sections 3 and 5 were checked again in the orchestration runtime when preserving this record. Their reproducible pseudocode is:

- Enumerate d=2,...,24 and e|d, 2<=e<d; maximize d*d/e+e*d-3*d+1.
- For d=18,...,24 compare twice min_r[binom(r,2)+binom(d-r,2)] with (d-1)^2-266.
- Enumerate r=0,...,d in the star inequality; retain only (24,23).
- For m>=8 use the explicit pi_2 residue formulas above, and the six explicit Hilbert-deficit polynomials.

These exact finite checks support arithmetic implications. They do not replace the geometric proofs, the classical theorem, the subgroup input, or independent review.
