# 6. Linearized moving degree and accessories

## 6.1 The torsor convention

For a faithful linear representation $V_0$ of a finite group $G$, put $L_0=\mathbb C(V_0)$ and $K=L_0^G$. An accessory is a finite extension $F/K$. Its chosen connected component has monodromy $H=\operatorname{Gal}(L_0F/F)\le G$. A curve compression means a faithful $H$-stable curve field $E\subset L_0F$ with $L_0F=FE$. The full base-changed torsor is induced from this component; faithfulness of a permutation action on disconnected components is insufficient.

Let $a(G)$ be the least degree bound permitting this compression for the generic torsor, and set $\mu(1)=1$. The target may have smaller monodromy after the accessory. Requiring connected full monodromy gives the separate threshold $\mu(G)$.

## 6.2 Divisorial correspondence and induction

**Theorem 6.1 (P).** An irreducible equivariant generically finite degree-$d$ cover $W\dashrightarrow V_0$, with a dominant equivariant map to a faithful $G$-curve $C$, gives a moving linearized bundle of degree at most $d$ on $C$.

*Proof.* The joint-image closure $Z\subset V_0\times C$ is an integral invariant Cartier divisor, of generic degree $e\mid d$ over $V_0$. Its divisor bundle has its canonical linearization. Homotopy invariance gives $\operatorname{Pic}(V_0\times C)=\operatorname{Pic}(C)$, so the generic fibre divisors lie in one degree-$e$ class. Restrict the linearization to the fixed fibre $\{0\}\times C$. This works even if its canonical section restricts to zero. The family moves because $Z$ dominates $C$, so the class has at least two sections. Removing the invariant fixed part preserves linearization and lowers degree. In fact a common point in all generic fibre divisors would give the vertical component $V_0\times\{p\}$, impossible for integral $Z$ dominating $C$; thus the joint-image family itself has no fixed part. $\square$

**Theorem 6.2 (P).** Under this convention,
$$a(G)=\min_{H\le G}[G:H]\mu(H).$$

*Proof.* For an accessory with monodromy $H$, the normalization of $V_0$ in $L_0F$ is an $H$-equivariant cover of degree $e=[L_0F:L_0]$, and $[F:K]=[G:H]e$. Apply Theorem 6.1 to its faithful compressed component. This gives the lower bound.

For the upper bound, the incidence variety for a linearized bundle $L$ of degree $m$ is the kernel of the surjective evaluation map $H^0(L)\otimes\mathcal O_C\to L$. It is irreducible, dominates $C$ and has degree $m$ over the section space. Induce the section representation from $H$ to $G$, add a faithful representation if necessary, and pull back the incidence through its identity-coset projection. Taking $H$-invariants in its field yields accessory degree $[G:H]m$ and the actual induced torsor compression.

There is also an upper bound for every torsor, avoiding a choice of linear model. For a $G$-torsor $T/F$, choose a closed point of $T/H$, of residue degree at most $[G:H]$, giving an $H$-reduction $P$. Twist $C,L$ by $P$. Genuine linearization descends an actual base-point-free degree-$m$ bundle on the twisted curve. Since the base field contains $\mathbb C$, a general section has a reduced zero divisor avoiding the finite nonfree locus. A closed point of that divisor has residue degree at most $m$. Its point on the twisted curve gives the pullback of $C\to C/H$ at a point defined over a field of transcendence degree at most one. Induction gives the required descent. A constant point yields a split torsor. $\square$

## 6.3 Linearized degree and projective genus

**Lemma 6.3 (P).** For quotient inertia orders $e_i$, every linearized bundle has degree in
$$N_G(C)\mathbb Z,\qquad N_G(C)=\frac{|G|}{\operatorname{lcm}(e_i)}.$$

*Proof.* Hilbert 90 supplies an invariant rational section. Its divisor is a sum of point orbits, whose sizes are $|G|$ and $|G|/e_i$. Their greatest common divisor is $N_G(C)$. $\square$

If the multiplier group has exponent $s$, an invariant class satisfies $N_G(C)\mid s\deg L$, since $L^s$ is linearizable. This necessary lattice can be weaker than the exact twisted lattice or descent under a free cyclic subgroup.

**Lemma 6.4 (P).** Let $G$ be nonabelian simple and a global minimizer for $\mu(G)$ have degree $n<|G|$. Its complete series is birational onto a faithful nondegenerate image in $\mathbb P^r$, with $r\ge q(G)-1$, and $g\le\pi(n,r)$. Here $q(G)$ is the least nontrivial genuine representation dimension.

*Proof.* The image-action kernel is normal and has order dividing the map degree, hence is trivial. A nonbirational map of degree $k\ge2$ would pull back a linearized moving hyperplane bundle on the faithful normalized image of degree $n/k$, contradicting minimality. The genuine section representation has a nontrivial constituent, so dimension at least $q(G)$. Apply arithmetic Castelnuovo to the birational image. $\square$

For $A_7$, $q=6$ and the element-order least common multiple is 420, so $6\mid\mu(A_7)$.

## 6.4 The threshold allowing monodromy drop

**Theorem 6.5 (P+X, with the subgroup classification).** $a(A_7)=60$.

*Proof.* If $\mu(A_7)<60$, it is at most 54 and Lemma 6.4 gives $g\le\pi(54,5)=325$. Exact signatures, the seven-sheeted Riemann--Hurwitz test, generating triples and the degree lattice leave only degree 36 on signature $(2,5,7)$ and degree 42 on $(3,4,5)$. Their projective genus bounds are respectively $136<199$ and $190<274$. Thus $\mu(A_7)\ge60$.

The subgroup classification gives all types of index below 60:

| Subgroup | Index | Lower bound for $\mu$ | Product |
|---|---|---|---|
| $A_6$ | 7 | 9 | 63 |
| $\mathrm{PSL}_2(7)$ | 15 | 4 | 60 |
| $S_5$ | 21 | 3 | 63 |
| $(A_4\times C_3):2$ | 35 | 2 | 70 |
| $A_5$ | 42 | 2 | 84 |

For $A_6$, $q=5$, minimum genus ten and $\pi(8,4)=5$ give its bound. For $\mathrm{PSL}_2(7)$, minimum genus three and $\pi(3,2)=1$ exclude degrees at most three. For $S_5$, degree one would require a Möbius embedding; degree two would give a hyperelliptic action with central involution and an impossible faithful quotient action. The order-72 subgroup cannot act on $\mathbb P^1$, because it contains $C_3^2$. Degree one for $A_5$ would require its projective action to lift to a genuine two-dimensional representation, which does not exist. Higher indices give products at least 60 automatically.

The canonical bundle on the Klein quartic has degree four and is genuinely $\mathrm{PSL}_2(7)$-linearized. Its index-15 induction gives degree 60. Apply Theorem 6.2. $\square$

Thus the relative essential-dimension obstruction for a single accessory holds through degree 59 and is sharp at 60. This convention permits monodromy drop.

## 6.5 A degree-90 series

The branch-divisor relations give
$$\operatorname{Pic}_G(C)=\langle D_2,D_4,D_7\mid2D_2=4D_4=7D_7\rangle
\simeq\mathbb Z\oplus\mathbb Z/2.$$
The degree-90 classes are $B=2D_7-D_4$ and $B+T_{\rm tor}$, with $T_{\rm tor}=D_2-2D_4$.

**Proposition 6.6 (P+X).** $B+T_{\rm tor}$ is base-point-free and has at least ten sections, so $\mu(A_7)\le90$.

*Proof.* For a nonidentity element, holomorphic Lefschetz gives
$$\operatorname{tr}(g\mid H^0-H^1)=\sum_{p\in\operatorname{Fix}(g)}\frac{a_p^{k_p}}{1-a_p^{-1}},$$
where $a_p$ is the tangent rotation and $k_p$ the divisor coefficient. Exact cyclotomic arithmetic gives
$$\chi(B+T_{\rm tor})=-6+10-14_a-14_b-21.$$
Its positive multiplicity forces a ten-dimensional section constituent. The fixed divisor is invariant and has degree at most 90, less than the minimum orbit 360. $\square$

The consistency checks are $\chi(B)=-\overline{10}-35$, $\chi(2B)=6-\overline{10}+14_a+14_b+21$ and
$$H^0(K_C)=10+2\overline{10}+15+21+2\cdot35.$$
Orientation reversal interchanges $10$ and $\overline{10}$. A positive Euler-characteristic multiplicity proves a lower bound for sections; negative multiplicities do not prove absence of sections, and the calculation does not assert $h^0(B+T_{\rm tor})=10$.

## 6.6 Excluding degrees below 90

First, degree 66 has no lattice-compatible signature. The remaining degree-60 candidates are $(2,6,7)$ and $(3,4,7)$; Castelnuovo excludes section dimension at least nine, leaving the ordinary standard six. An $A_6$-fixed section would have an $A_6$-invariant zero divisor of degree 60. All $A_6$-orbits on those curves have size at least 120 and 90 respectively, impossible. Hence $\mu\ge72$. Degree 78 likewise has no candidate.

For a candidate minimum $n=72$ or 84, all point orbits exceed $n$, and the signature and Castelnuovo tests give $h^0\le12$. There are no invariant sections, so the genuine section representation is $6,10,\overline{10}$ or $6\oplus6$.

**Lemma 6.7 (orbit semigroup; P).** A degree-$k$ invariant polynomial restricts to zero if $kn$ is not a nonnegative combination of point-orbit sizes.

*Proof.* A nonzero restriction would have an invariant zero divisor of degree $kn$. $\square$

For a standard-six subspace, write $x_1,\ldots,x_7$ with $\sum x_i=0$ and use power sums $p_k$. For $n=72$, the lemma kills $p_2,p_3,p_4,p_6$. Newton identities put the image in the projective curve of ordered roots of $t^7-ut^2-v$, of degree $2\cdot3\cdot4\cdot6=144$. Its monodromy is $S_7$, so the root curve is irreducible and generically reduced; its finite root-configuration map rules out extra components. A curve of degree at most 72 cannot be that image. For $n=84$, it kills $p_2,p_3,p_4,p_7$. Newton gives $p_7=7e_7$, so the image lies in a coordinate hyperplane, and transitivity puts it in every such hyperplane, impossible.

For the ten-dimensional representations, use the half-spin $V_4$ of $2.A_7$, with $\wedge^2V_4=6$ and $\operatorname{Sym}^2V_4=10$ or $\overline{10}$. The curve is a family of symmetric $4\times4$ matrices $Q_p$. The determinant restriction vanishes by the orbit semigroup, so generic rank is at most three.

**Lemma 6.8 (parity; P).** For a lift with eigenvalues $\mu_i$, a fibre eigenvalue $\lambda$ permits only entries with $\mu_i\mu_j=\lambda$. If no $\mu_i^2=\lambda$, every matrix in that eigenspace has even rank.

*Proof.* The eigenvalue blocks pair under the fixed-point-free involution $\mu\mapsto\lambda/\mu$, giving matrices of the form $\left(\begin{smallmatrix}0&A\\A^T&0\end{smallmatrix}\right)$. $\square$

At generic rank two, the second compound matrix has rank one on $\wedge^2V_4$. Taking its inverse Veronese image gives a standard-six map with a linearized bundle $N$ and $N^2=L^2(-B')$. The base divisor $B'$ is invariant and has degree at most $2n$, below every available orbit, so $B'=0$ and $\deg N=n$. This is the excluded standard-six case.

At generic rank three, the adjugate gives a spin kernel map with pullback $\mathcal O(2)=L^3(-B'')$. Exact branch patterns show that at some full branch orbit every singular matrix has rank at most two. Hence every adjugate entry vanishes there. That orbit has size at least 420, exceeding $3n$, impossible.

At rank one, the inverse Veronese is a spin map. Lemma 6.8 excludes it at the order-two branch for $(2,5,7)$ and the order-four branch for $(3,4,5)$. The remaining two $(3,5,6)$ cases yield a degree-42 spin map in $\mathbb P^3$. Its normalized image is faithful. A nonbirational image has degree at most 21 and genus at most $\pi(21,3)=90<136$, so it is birational. Its invariant quadratic ideal is zero because $\operatorname{Sym}^2V_4$ is irreducible. Halphen's no-quadric bound gives $g\le42\cdot39/6+1=274<379$.

**Finite certificate.** Lift spectra are solved exactly by matching their pairwise products to the standard-six spectrum and requiring determinant one. The solutions differ only by sign, or by duality for seven-cycles. A branch eigenspace is a symmetric zero-pattern space. Rank one is possible precisely when the pattern has a diagonal entry; rank three is decided by the nonvanishing of a $3\times3$ minor on the determinant-zero locus, using exact polynomial factorization and radical membership. `verify_mu90_exact.py` closes all 84 curve/class cases, in both orientations. It does not infer matrix rank from floating eigenvectors.

**Theorem 6.9 (P+X).** $\mu(A_7)=90$.

The upper bound is Proposition 6.6; the lattice, standard-six obstruction, exact rank patterns and Halphen argument exclude every smaller multiple of six after Theorem 6.5. No numerical projective equations enter this result.

## 6.7 Amitsur subgroups and the Albanese hypothesis

For a smooth projective $G$-variety $X$, define
$$\operatorname{Am}_G(X)=m_X(\operatorname{Pic}(X)^G)\subset H^2(G,\mathbb C^*).$$
A fixed point makes this subgroup zero by evaluating the multiplier on its fibre.

**Theorem 6.10 (P).** Let $X$ be smooth projective and $C$ a faithful $G$-curve. Suppose
$$\operatorname{Hom}_G(\operatorname{Alb}X,\operatorname{Jac}C)=0.$$
An irreducible equivariant degree-$d$ cover of $X$ dominating $C$ gives invariant bundles $M$ on $X$ and $L$ on $C$ with
$$\deg L=e\mid d,\quad |L|\text{ base-point-free},\quad h^0(L)\ge2,
\quad m_C(L)=-m_X(M).$$

*Proof.* The joint-image closure $Z\subset X\times C$ is an invariant integral Cartier divisor. The canonical, equivariant correspondence sequence
$$0\to\operatorname{Pic}X\oplus\operatorname{Pic}C\to\operatorname{Pic}(X\times C)
\to\operatorname{Hom}(\operatorname{Alb}X,\operatorname{Jac}C)\to0$$
puts its class in the external-product subgroup, by the stated vanishing. Unique factors give invariant $M,L$ with $\mathcal O(Z)=M\boxtimes L$. The fibre divisors have degree $e\mid d$ and move. A common point would be a vertical component of integral $Z$, contradicting dominance of $C$, so they have no fixed part. Additivity of multiplier classes and the canonical linearization of $\mathcal O(Z)$ give the obstruction identity. $\square$

This requires global vanishing of the equivariant Albanese homomorphisms. Constancy for one correspondence gives its own bound, but does not classify all correspondences.

On an open base, $\operatorname{Pic}(B)=0$ alone is insufficient. For the $C_2^2$ action on $\mathbb P^1$ generated by $z\mapsto-z$ and $z\mapsto1/z$, $\mathcal O(1)$ is not linearized, although the invariant open $B=\mathbb P^1\setminus\{0,\infty,\pm1\}$ has trivial Picard group and gives a degree-one correspondence. The lifts anticommute. Similarly the $L_{60}$ incidence gives degree 60 over $\mathbb P(H^0(L_{60}))$ minus an orbit of hyperplanes, a rational Picard-trivial base.

The open obstruction lives in $H^2(G,\mathcal O(B)^*)$. Its constant classes can die through
$$H^1(G,\mathcal O(B)^*/\mathbb C^*)\longrightarrow H^2(G,\mathbb C^*).$$
Constant units together with linearizability of all invariant bundles, or a fixed point alone, restore the relevant conclusion. The punctured quadratic cone $t^2=q(v)$ in the standard-six example is Picard-trivial with constant units: its projective quadric has Picard group generated by $\mathcal O(1)$, and removing the vertex changes neither global functions nor the cone's class group. This specific quadratic example therefore retains the degree-90 lower bound, although the unrestricted open-base criterion is false.

## 6.8 Stable compression and covers

**Theorem 6.11 (Amitsur kernel; P).** For a smooth projective generically free $X$,
$$\operatorname{Am}_G(X)=\ker\left(H^2(G,\mathbb C^*)\to H^2(G,\mathbb C(X)^*)\right)
=\ker\left(H^2(G,\mathbb C^*)\to\operatorname{Br}(\mathbb C(X)^G)\right).$$

*Proof.* Combine the exact sequences of constants, rational functions and principal divisors, and of principal divisors, divisors and line bundles. Hilbert 90 kills the first cohomology of rational functions; the divisor group is a permutation module with zero first cohomology. The resulting kernel is the image of invariant Picard classes. The relative Brauer description gives the second equality. Injectivity of $\operatorname{Br}(K)\to\operatorname{Br}(K(t))$ also proves stable equivariant birational invariance. $\square$

For $A\subset H^2(G,\mathbb C^*)$, let $\mu_A(C)$ be the minimum degree of an invariant moving class with multiplier in $A$. Define $c_X^{\rm st}(C)$ as the minimum degree of an irreducible equivariant cover dominating $C$, whose base is equivariantly birational to $X\times\mathbb P^N$ with trivial action on the extra factor.

**Theorem 6.12 (P).** Under Theorem 6.10's global Albanese hypothesis and generic freeness,
$$c_X^{\rm st}(C)=\mu_{\operatorname{Am}_G(X)}(C).$$

*Proof.* The lower bound is Theorem 6.10 and stable invariance of Albanese and Amitsur subgroup. For the upper bound choose a minimizing $L$; removing an invariant fixed divisor shows it is base-point-free. Choose $M$ on $X$ with opposite multiplier. Then $M\otimes H^0(L)$ is a genuinely linearized vector bundle. Over the free locus it descends to the quotient and has a basis over the generic point. Its projective bundle is therefore equivariantly birational to $X\times\mathbb P^{h^0(L)-1}$ with trivial extra action. Pull back the universal section incidence. This irreducible variety is a projective bundle over $C$ in the section parameter, has degree $\deg L$ over the base and dominates $C$. $\square$

For the specified $(2,4,7)$ target, Chapter 7's twisted lattice and low-degree exclusion give:

| $\operatorname{Am}_{A_7}(X)$ | Possible degree lattice | Stable minimum |
|---|---|---|
| 0 | $90\mathbb Z$ | 90 |
| Order two | $45\mathbb Z$ | 90 |
| Order three | $30\mathbb Z$ | 60 |
| Order six | $15\mathbb Z$ | 60 |

This table retains the lower-bound certificate used to exclude degree 30, and it requires the Albanese hypothesis. The order-two minimum 90 does not make every possible degree divisible by 90.

For a degree-$d$ equivariant generically finite map $Y\dashrightarrow X$ between generically free projective models,
$$\operatorname{Am}_G(X)\subset\operatorname{Am}_G(Y),\qquad
d\operatorname{Am}_G(Y)\subset\operatorname{Am}_G(X).$$
Restriction proves the first inclusion; restriction followed by corestriction is multiplication by $d$, proving the second. Primary parts at primes not dividing $d$ agree.

## 6.9 Brauer splitting and the tower limitation

Twisting an $n$-dimensional projective representation of multiplier $\alpha$ gives a Severi--Brauer variety; its Brauer index divides $n$. The projective degree lists in the Schur sectors have greatest common divisors two, three and six. Over a generic linear base the Amitsur subgroup is zero, so period equals multiplier order. Since period divides index, generic indices are exactly two, three and six. This uses the sector lists as representation-theoretic inputs, not merely their integer gcds.

A maximal subfield of a division algebra of index $r$ gives a splitting extension of degree $r$. Quadratic and cubic steps can therefore remove the two- and three-parts. Their total degree at most six retains full $A_7$ monodromy, since $A_7$ has no proper subgroup of index at most six. These small splitting steps do not themselves provide a curve compression, and their total-space Albanese is not determined. A cubic step followed by degree-60 compression has total degree 180, compatible with $\mu=90$.

For a general joint-image correspondence, the fibre class defines an affine map into the Picard torsor,
$$a_Z:X\dashrightarrow\operatorname{Pic}^e(C),\qquad x\mapsto[Z_x],$$
with equivariant linear part $u_Z:\operatorname{Alb}X\to\operatorname{Jac}C$. It must land in $W_e(C)$ and lift rationally to $\operatorname{Sym}^eC$. If $u_Z=0$, Theorem 6.10's argument gives the lower bound for that correspondence. The stable minimum formula needs all equivariant homomorphisms to vanish. The diagonal on $X=C$ has $e=1$ and $u_Z=\mathrm{id}$, showing that no universal gonality/Amitsur bound holds for overlapping bases.

An origin need not be $G$-fixed; the object is an affine map of Albanese and Picard torsors, not just a linear map. Shared representation constituents in $H^1$ do not alone give a Hodge homomorphism or an effective correspondence. Controlling these maps on actual tower bases remains the precise unresolved problem.
