# Equivariant normalization defects and the degree-60 model

4 October 2026. Baseline: `5a4fd13622239afffc236ff923dcb48a25b8afda`.

**Results [P][X].** On every one of the four `(2,4,7)` covers, the six-section degree-60 morphism is a closed embedding. The involution pencil has base divisor exactly the 18 fixed points, all simple, and moving degree exactly 42; it descends to a degree-21 pencil. These statements use the exact inputs of Theorem 7.5, the minimum faithful genus 136, and Castelnuovo's arithmetic-genus bound. They do not use the spectral lower bound or the numerical equations of Proposition 7.8. They do not establish that 42 or 21 is the gonality.

## 1. Defects of an equivariant normalization

Work over an algebraically closed field of characteristic zero. Let
\[
\nu:C\longrightarrow Y\subset\mathbb P^r
\]
be a finite birational morphism, where `C` is smooth and `Y` is integral and nondegenerate of degree `n`. A finite group `G` acts equivariantly. Put
\[
D=p_a(Y)-g(C),\qquad B=\pi(n,r)-g(C).
\]
The normalization sequence gives
\[
D=\sum_{y\in Y}\delta_y\le B,
\quad
\delta_y=\dim_k\big((\nu_*\mathcal O_C)_y/\mathcal O_{Y,y}\big).
\tag{1}
\]
Here **Castelnuovo bounds the arithmetic genus of the integral image**, not just its normalization. If `n-1=q(r-1)+s`, `0<=s<r-1`, then
\[
\pi(n,r)=(r-1)\binom q2+qs.
\]

**Lemma 1 (local first jets).** Let `b_y` be the number of points above `y`. Choose local parameters `t_i` at them and an affine chart of the ambient projective space. Let `a_y` be the rank of the matrix whose columns are the tangent vectors of the normalization map along these branches. Then
\[
\delta_y\ge 2b_y-1-a_y.
\tag{2}
\]
In particular, if `k_y` of these branches have zero derivative, then
\[
\delta_y\ge b_y-1+k_y.
\tag{3}
\]

*Proof.* The completed normalization is
`S=\prod_{i=1}^{b_y} k[[t_i]]`. Modulo `\prod_i(t_i^2)`, its dimension is `2b_y`. The image of the completed local ring of `Y` has a common constant term on all branches, and its linear coefficients span a space of dimension `a_y`. Thus this image has dimension `1+a_y`. The quotient of `S` by the local ring surjects onto the corresponding quotient of first jets, proving (2). Since `a_y<=b_y-k_y`, (3) follows. No plane-singularity hypothesis is used. QED.

**Theorem 2 (orbit bounds).** Let `R` be the set of points at which `d\nu=0`, and `S` the set of points belonging to a fibre with at least two points. Then
\[
|R|\le D\le B,
\qquad
|S|\le2D\le2B.
\tag{4}
\]
Both sets are finite and `G`-invariant. Consequently:

- If every point orbit on `C` has size greater than `B`, then `\nu` is an immersion.
- If every point orbit has size greater than `2B`, then `\nu` is a closed embedding.
- If `\nu` is an immersion, to prove a closed embedding it suffices to check injectivity on the union of the point orbits of size at most `2B`.

*Proof.* Sum (3) to obtain the first bound. On a fibre with `b_y>=2`, `b_y<=2(b_y-1)<=2\delta_y`, which gives the second. For the last assertion, any collision set must be contained in the stated union, by equivariance and (4). An injective immersion of a smooth proper curve is a closed embedding; equivalently the normalization is then an isomorphism. QED.

The last test is a finite orbit calculation, rather than a check of equations for the whole curve.

There is also a stabilizer form. If `J=G_y` and the point stabilizers on the fibre are `I_1,...,I_t`, one for each `J`-orbit, then
\[
b_y=\sum_{i=1}^t[J:I_i],\qquad
D\ge\sum_{[y]}[G:J]\,(2b_y-1-a_y).
\tag{5}
\]
The sum runs over singular image orbits. Higher jets give further bounds by replacing first jets with `\prod k[t_i]/(t_i^m)` and computing the rank of the restricted local functions. Formula (5) keeps the stabilizers and the tangent geometry in the same inequality.

## 2. The degree-60 embedding

Let `C` be a `(2,4,7)` `A7`-curve and `V` the exceptional six-dimensional section subspace of `L60`. Its morphism is
\[
\varphi:C\longrightarrow\mathbb P(V^*).
\]

**Theorem 3.** This morphism is a closed embedding.

*Proof.*

1. **Degree and birationality.** The section subspace has no base divisor: such a divisor would be invariant, whereas every point orbit has size at least 360, exceeding 60. If the morphism had degree `e>=2` onto its image, the image degree would be `60/e<=30`. The normalized image has a faithful `A7`-action: the kernel is normal, and a trivial image action would force `2520|e`, impossible. Every faithful `A7`-curve has genus at least 136, whereas Castelnuovo gives `g<=pi(30,5)=91`. Hence `e=1`.

2. **The defect bound.** Write `Y=\varphi(C)`. It is an integral nondegenerate degree-60 curve in `P5`, with
   \[
   p_a(Y)\le\pi(60,5)=406,
   \qquad D\le406-136=270.
   \tag{6}
   \]
   The point orbit sizes on `C` are `2520,1260,630,360`, with the last three being the single reduced branch fibres `D2,D4,D7`.

3. **Immersion.** A nonempty ramification set contains an orbit of size at least 360, contradicting `|R|<=270`. Thus `\varphi` is an immersion everywhere, including the 630 four-points.

4. **Only the seven-points could collide.** By (4), `|S|<=540`. Every orbit except `D7` has size at least 630. Thus any nonempty collision set would be contained in `D7`.

5. **The seven-points are separated exactly.** An order-7 lift on `V` or `V*` has characteristic polynomial
   \[
   \Phi_7(T)=T^6+T^5+\cdots+T+1.
   \]
   Each of its six eigenlines has projective stabilizer exactly `C7` in `A7`. This is checked in exact arithmetic in `verify_normalization.py`; the specialization argument is given below. Since `D7` is the transitive orbit `A7/C7`, its equivariant map to the orbit of any such eigenline is a bijection. Consequently the map is injective on `D7`, contradicting a nonempty `S`.

So `S=R=\varnothing`, proving the theorem. QED.

**Exact eigenline certificate.** The script reuses the integral ATLAS matrices and checks every Cayley edge as in `audit_exact.py`. It normalizes an order-7 lift and checks its characteristic polynomial over `Z[omega]`, `omega^2+omega+1=0`. Reduce at `p=43`, under both embeddings `omega=6,36`, and choose a primitive seventh root in `F43`. Every spectral projector
\[
P_j=\frac17\sum_{k=0}^6\zeta_7^{-jk}M^k\qquad(1\le j\le6)
\]
remains rank one and nonzero. The program checks all 2520 group elements against every eigenline, in the representation and its dual. Every stabilizer is exactly the seven-element subgroup. These are 24 exact stabilizer checks.

Any element stabilizing an eigenline in characteristic zero also satisfies the same two-by-two-minor identities after specialization. Therefore the finite-field stabilizer is an **upper bound** for the characteristic-zero stabilizer. The lower bound is its order-7 subgroup itself. This proves equality, rather than inferring it from numerical eigenvectors. The projectors give nonzero specialized columns, so no disappearing eigenvector is assumed.

## 3. The exact involution pencil

Let `tau` be an involution and `R_tau` its reduced fixed divisor, of degree 18. Its order-2 lift has section eigenspaces of dimensions `4,2`, and acts as `+1` on every fixed-point fibre (the exact inputs of Theorem 7.5).

**Corollary 4.** The anti-invariant two-dimensional section space has fixed divisor exactly `R_tau`, with multiplicity one at every point. Its moving part is a base-point-free `g1_42`, pulled back from a base-point-free `g1_21` on `C/tau`.

*Proof.* Its base locus is the inverse image of `P(E_+)` under `\varphi`. If a point maps there, `\varphi(tau p)=\varphi(p)`, so injectivity forces `tau p=p`. Conversely every fixed point lies there by its fibre character.

At a fixed point the tangent action on `C` is `-1`. In the ambient tangent space at a point of `P(E_+)`, the `-1` directions are precisely those towards `E_-`. Since `\varphi` is an immersion, one anti-invariant section has a nonzero first derivative. Thus the common zero has order exactly one. Removing `R_tau` leaves degree `60-18=42`. Ratios of anti-invariant sections are invariant under `tau`, so the morphism factors through the degree-two quotient, leaving degree 21. QED.

This replaces **both** numerical obligations in the old Proposition 7.7: immersion at the four-points and absence of base points off the fixed divisor. It does not establish the entire minimal vanishing sequence at every branch point: immersion proves its first two orders are `(0,1)`, which is all the pencil argument needs.

## 4. Verification and sources

```
python3 gpt/verify_normalization.py
python3 audit_exact.py
```

Saved outputs: `normalization_check.txt`, `normalization_inputs_check.txt`. The minimum faithful genus is independently established by `verify_accessory60.py`; it is the existing algebraic genus input, not a spectral input.

The normalization sequence and local jet estimates above are elementary. For the arithmetic-genus scope of Castelnuovo, see Buczynski–Ilten–Ventura, *Singular Curves of Low Degree and Multifiltrations from Osculating Spaces*, arXiv:1905.11860v3, Introduction, p. 2, and Harris, *A bound on the geometric genus of projective varieties*, Ann. Scuola Norm. Sup. Pisa (4) 8 (1981), 35–68. The former explicitly states that Castelnuovo also bounds the arithmetic genus of integral singular curves.

**Remaining numerical claims.** This proof gives a closed embedding and an exact degree-42 pencil. It does not certify the cubic-and-quartic ideal, its Hilbert function, projective normality, the Hessian identification, quotient Koszul ranks, or exact gonality.

The companion `A7_Quadratic_Systems.md` proves the absence of quadrics and the Klein contact divisors exactly, and identifies the Jacobian/apolar quadratic ideals. The remaining higher restriction ranks and ideal assertions above are separate.
