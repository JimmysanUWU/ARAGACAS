# Curve compressions, Amitsur subgroups, and Albanese variation

**Research note and partial audit • 3 October 2026**

**Repository read:** `JimmysanUWU/ARAGACAS`, branch `claude/continue-previous-qfhm7j`, commit `ea25d8ae94c179cf5315b4e68f1af0f33a7d2bbd`. Chapter numbers below refer to that revision. The guide, seven chapters, earlier-work digest, and relevant script sections were read. This is not a completed certification of all their computational inputs.

## 1. The finding that changes the tower question

Theorem 5.11 has a correct central idea: when the equivariant Albanese-to-Jacobian correspondence vanishes, a family of effective divisors on the target curve stays in a single invariant linear system. Its additional assertion that triviality of the Picard group of an **open** base forces that system to be linearised is false.

The missing datum is the action on invertible regular functions. A nontrivial constant Schur cocycle can become a coboundary with nonconstant invertible coefficients. This can happen even when every line bundle on the open base is trivial.

There is a clean replacement. For a smooth projective base $X$, retain its classical Amitsur subgroup

$$\operatorname{Am}_G(X)=\operatorname{im}\bigl(\operatorname{Pic}(X)^G\longrightarrow H^2(G,\mathbb C^*)\bigr).$$

Under the same Albanese hypothesis, it determines exactly which invariant moving bundles can occur **after stable equivariant birational modification of the base**. Theorem 4 below proves this, including the converse. A finite cover of degree $d$ satisfies

$$\operatorname{Am}_G(X)\subseteq\operatorname{Am}_G(Y),\qquad d\operatorname{Am}_G(Y)\subseteq\operatorname{Am}_G(X).$$

For $A_7$, the associated generic Brauer classes have indices $2$, $3$, and $6$, respectively. In particular, the order-three obstruction can be removed by a cubic extension with full $A_7$ monodromy retained. These facts explain both the usefulness and the limitation of this framework for towers.

**Status.** The counterexample and general theorems below are proved here. The numerical thresholds $60$ and $90$ use the explicitly stated Chapter 5 and Chapter 7 inputs; this note does not independently certify all of those inputs. The Amitsur subgroup itself is established mathematics, not a new invariant [1]. The stable compression formulation is developed here for this project; no literature-priority claim is made.

## 2. Two counterexamples to the open-base criterion

### 2.1 An unconditional example on a projective line

Let $G=C_2\times C_2$ act on $C=\mathbb P^1$ by $s(z)=-z$ and $t(z)=1/z$. Set

$$B=\mathbb P^1\setminus\{0,\infty,1,-1\},\qquad W=B.$$

The identity $W\to B$ is a finite equivariant cover of degree one, and the inclusion $W\to C$ is dominant and equivariant. The smooth projective compactification is $X=\mathbb P^1$, so its Albanese is zero. Moreover,

$$\mathcal O(B)=\mathbb C[z,z^{-1},(z^2-1)^{-1}],\qquad\operatorname{Pic}(B)=0.$$

Nevertheless, $\mathcal O_C(1)$ is not $G$-linearised. Such a linearisation would lift the faithful projective action to a two-dimensional representation of $C_2\times C_2$. Every such representation is a sum of two characters, and its projective image has order at most two. Equivalently, the usual matrix lifts of $s$ and $t$ anticommute. The bundle $\mathcal O_C(2)=K_C^{-1}$ is linearised, so $\mu(G)=2$ while the cover above has degree one.

This disproves Theorem 5.11(2) as written, including its more general premise that every invariant line bundle on $B$ is linearisable. Here all underlying line bundles on $B$ are trivial and admit a linearisation, but that does not eliminate the constant obstruction carried by the target.

### 2.2 The same phenomenon inside the $A_7$ project

Use only the existence assertion of Theorem 7.5: a base-point-free invariant bundle $L=L_{60}$ with a six-dimensional projective space of sections $V\subseteq H^0(C,L)$ and an order-three multiplier. No embedding, equation of the cubic, or exactness of the 42-pencil is needed.

Let $X=\mathbb P(V)$ parametrise section lines. Its incidence variety is

$$I=\{([s],p)\in\mathbb P(V)\times C:s(p)=0\}.$$

The map $I\to X$ is finite of degree $60$: every nonzero section has a zero divisor of that degree. Also, $I\to C$ is the projective bundle of the kernel of the evaluation map $V\otimes\mathcal O_C\to L$. Thus $I$ is irreducible, and it dominates $C$.

The universal divisor bundle is $\mathcal O_X(1)\boxtimes L$. Its two multipliers cancel, so $I$ has a genuine $G$-action. Choose a hyperplane $H\subset X$ and remove its entire $G$-orbit. The resulting nonempty invariant open $B$ is rational and has $\operatorname{Pic}(B)=0$: the boundary contains a hyperplane, whose class generates $\operatorname{Pic}(X)$. Restricting $I$ to $B$ gives an irreducible equivariant degree-60 cover dominating $C$.

Thus $90\mid d$ and $d\ge90$ both fail for **arbitrary rational open bases with trivial Picard group**. This does not contradict the degree-90 theorem over the original linear base. The projective action on this new base has a different Amitsur subgroup.

To see the missing units explicitly, choose lifts $R_g$ with $R_gR_h=\alpha(g,h)R_{gh}$ and a linear form $\ell$ defining $H$. On $B$ the functions

$$u_g([v])=\frac{\ell(R_gv)}{\ell(v)}$$

are invertible and satisfy

$$u_g(h[v])u_h([v])=\alpha(g,h)u_{gh}([v]).$$

The nontrivial constant cocycle $\alpha$ has become a coboundary in the units of $B$.

## 3. The corrected correspondence theorem

Work over $\mathbb C$. All varieties are integral. Write $m_X(M)$ for the obstruction in $H^2(G,\mathbb C^*)$ to linearising an invariant line bundle on a smooth projective $G$-variety $X$.

**Theorem 1 (projective correspondence).** Let $X$ be smooth projective with a generically free $G$-action, and let $C$ be a smooth projective faithful $G$-curve. Assume

$$\operatorname{Hom}_G(\operatorname{Alb}X,\operatorname{Jac}C)=0.$$

Suppose $W\dashrightarrow X$ is an equivariant generically finite map of degree $d$ and $W\dashrightarrow C$ is dominant and equivariant. There are invariant line bundles $M$ on $X$ and $L$ on $C$ such that

$$\deg L=e\mid d,\quad h^0(C,L)\ge2,\quad m_C(L)=-m_X(M).$$

In fact $L$ is base-point-free. In particular $m_C(L)\in\operatorname{Am}_G(X)$.

**Proof.** Let $Z\subset X\times C$ be the closure of the joint image. It is an irreducible invariant divisor dominating both factors. Its generic degree $e$ over $X$ divides $d$, by the intermediate field $\mathbb C(Z)\subseteq\mathbb C(W)$. Since $X\times C$ is smooth, $\mathcal O(Z)$ is an invertible sheaf with its canonical $G$-linearisation.

Use the canonical equivariant exact sequence

$$0\longrightarrow\operatorname{Pic}(X)\oplus\operatorname{Pic}(C)\longrightarrow\operatorname{Pic}(X\times C)\longrightarrow\operatorname{Hom}(\operatorname{Alb}X,\operatorname{Jac}C)\longrightarrow0.$$

The last component of $[Z]$ is equivariant and therefore zero. Hence $\mathcal O(Z)=M\boxtimes L$. The injection of the first two summands is canonical, so their classes are invariant. No equivariant splitting of the whole exact sequence is being asserted or needed.

Over a dense open of $X$, the divisors $Z_x$ belong to $|L|$. They cannot all be the same divisor, since $Z$ dominates $C$, so $h^0(L)\ge2$. If a point $p$ belonged to every $Z_x$, then $X\times\{p\}$ would be an irreducible component of $Z$. That contradicts irreducibility and dominance over $C$. Therefore this family, and hence $|L|$, has no fixed point. Finally the obstruction of $M\boxtimes L$ is $m_X(M)+m_C(L)$ and is zero. This proves the theorem.

**Remark on divisibility.** One need not remove a fixed divisor and risk losing $e\mid d$. Irreducibility of the correspondence already gives base-point-freeness.

**Corollary 2 (safe forms of Theorem 5.11).** The linearised conclusion holds if $\operatorname{Am}_G(X)=0$. It also holds using an invariant open base $B$ if both every invariant line bundle on $B$ is linearisable and the natural map

$$H^2(G,\mathbb C^*)\longrightarrow H^2(G,\mathcal O(B)^*)$$

is injective. Sufficient conditions for this injectivity include $\mathcal O(B)^*=\mathbb C^*$ or a $G$-fixed point of $B$; in the latter case evaluation splits the inclusion of constants. In fact a fixed point alone already linearises $L$, by restricting the canonical linearisation of $\mathcal O(Z)$ to that fibre. Constant units alone still require the stated hypothesis on invariant line bundles.

On an open base, the correct obstruction equation is

$$\delta_B(M)+\iota(m_C(L))=0\quad\hbox{in }H^2(G,\mathcal O(B)^*).$$

Even if $M$ is linearisable, this only says that $m_C(L)$ lies in the kernel of $\iota$. The units exact sequence identifies that kernel as the image of

$$H^1(G,\mathcal O(B)^*/\mathbb C^*)\longrightarrow H^2(G,\mathbb C^*).$$

This is the precise gap in the open-base argument. Related unit-splitting criteria already occur in the universal-torsor literature [1, Proposition 4].

## 4. A sharp stable classification

For a subgroup $A\subseteq H^2(G,\mathbb C^*)$, define

$$\mu_A(C)=\min\{\deg L:[L]\in\operatorname{Pic}(C)^G,\ h^0(L)\ge2,\ m_C(L)\in A\}.$$

A minimising bundle is base-point-free: its fixed divisor is invariant and canonically linearised, so removing it preserves the multiplier and lowers the degree.

Define $c_X^{\mathrm{st}}(C)$ to be the least degree of an irreducible equivariant cover dominating $C$ whose base is equivariantly birational to $X\times\mathbb P^N$ for some $N$, with trivial $G$-action on the extra projective factor. The cover is required to be finite only after passage to an invariant dense open, as in the field-theoretic accessory problem.

**Theorem 3 (Amitsur kernel).** For smooth projective generically free $X$,

$$\operatorname{Am}_G(X)=\ker\bigl(H^2(G,\mathbb C^*)\longrightarrow H^2(G,\mathbb C(X)^*)\bigr).$$

**Proof.** Combine the exact sequences of constants, rational functions and principal divisors, and of principal divisors, divisors and line bundles. Hilbert 90 gives $H^1(G,\mathbb C(X)^*)=0$. The divisor group is a permutation module, so its first cohomology is zero. The resulting kernel is exactly the image of $\operatorname{Pic}(X)^G$. Equivalently, these are the constant multiplier classes whose crossed-product Brauer classes vanish over $\mathbb C(X)^G$.

This description proves birational invariance. Adding variables with trivial $G$-action does not change the kernel, since the Brauer group of a field injects into that of a purely transcendental extension. Thus it also proves the stable invariance needed below.

**Theorem 4 (stable compression formula).** Under the hypotheses on $X,C$ in Theorem 1,

$$\boxed{c_X^{\mathrm{st}}(C)=\mu_{\operatorname{Am}_G(X)}(C).}$$

**Proof of the lower bound.** Stabilisation does not change the Albanese or Amitsur subgroup. Apply Theorem 1 to any proposed cover. It gives a moving bundle of degree $e\mid d$ with multiplier in the indicated subgroup, and hence $d\ge e\ge\mu_{\operatorname{Am}_G(X)}(C)$.

**Proof of the upper bound.** Take a minimising bundle $L$ and put $V=H^0(C,L)$. Choose an invariant bundle $M$ on $X$ with $m_X(M)=-m_C(L)$. Compatible projective actions on $M$ and $V$ make $E=M\otimes V$ a genuine $G$-linearised vector bundle on $X$.

Over the free locus of $X$, the vector bundle $E$ descends to the quotient. At its generic point, that descended vector bundle has a basis. Consequently its projective bundle satisfies

$$\mathbb P_X(E)\dashrightarrow X\times\mathbb P^{\dim V-1}$$

equivariantly and birationally, with trivial action on the last factor. This is the elementary descent proof of the no-name lemma in the form needed here.

As a projective bundle with diagonal action, $\mathbb P_X(E)$ is $X\times\mathbb P(V)$. Pull back the universal incidence divisor of $|L|$. It is $X\times I_L$, irreducible because $I_L\to C$ is a projective bundle, and finite of degree $\deg L$ over $X\times\mathbb P(V)$. It dominates $C$. Transfer this construction through the displayed birational equivalence and shrink the base. It realises the required degree.

**Meaning.** With the Albanese correspondence excluded, the dependence on a base is completely captured by which Schur classes it cancels, after the stated stabilisation. This is stronger than a necessary bound and does not require a fixed point.

## 5. Consequences for a $(2,4,7)$ $A_7$ target

Here use the following upstream inputs explicitly: Corollary 7.2 gives the multiplier-degree congruence; Theorem 7.4 excludes moving invariant classes below $60$; Theorem 7.5 supplies $L_{60}$ with order-three multiplier; Proposition 5.8 supplies the linearised degree-$90$ bundle. For the fixed target, the lattice already forces every nonzero linearised degree to be a multiple of $90$.

Let $A=\operatorname{Am}_{A_7}(X)$, a subgroup of $\mathbb Z/6$, and continue to assume the Albanese hypothesis. Theorem 4 gives the following conditional-on-these-inputs exact table.

| $A$ | Necessary lattice for $e=\deg L$ | $c_X^{\mathrm{st}}(C)$ |
|:---|:---|:---|
| $0$ | $90\mathbb Z$ | $90$ |
| subgroup of order $2$ | $45\mathbb Z$ | $90$ |
| subgroup of order $3$ | $30\mathbb Z$ | $60$ |
| $\mathbb Z/6$ | $15\mathbb Z$ | $60$ |

For the order-two row, degree $45$ is excluded by Theorem 7.4 and degree $90$ is available. For the order-three and full-subgroup rows, degree $60$ is both available and minimal. The table concerns the minimum; it does not assert $90\mid e$ in the order-two row.

Over an ordinary faithful linear base, compactify by $\mathbb P(V\oplus\mathbb C)$. The fixed point implies $A=0$. Thus the original fixed-target assertion $90\mid d$ survives unchanged there. Over the projective base $\mathbb P(H^0(L_{60}))$, the Amitsur subgroup is generated by the order-three multiplier; this explains the incidence cover of degree $60$.

**The particular quadratic accessory in Corollary 5.12(2) survives.** For the nondegenerate standard quadratic form in six variables, the punctured affine cone $t^2=q(v)$ is rational, factorial, and has only constant units. One way to see factoriality is that the base is a smooth projective quadric of dimension five with Picard group generated by $\mathcal O(1)$; the cone class group is its quotient by that generator. Removing the vertex changes neither the class group nor global units. Corollary 2 therefore applies to this particular base. It was the unrestricted open-base criterion, not this specific example, that failed.

## 6. What finite covers can change

**Theorem 5 (degree and Amitsur growth).** If $Y\dashrightarrow X$ is a dominant equivariant generically finite map of degree $d$ between smooth projective generically free $G$-varieties, then

$$\operatorname{Am}_G(X)\subseteq\operatorname{Am}_G(Y),\qquad d\operatorname{Am}_G(Y)\subseteq\operatorname{Am}_G(X).$$

**Proof.** Put $K=\mathbb C(X)^G$ and $F=\mathbb C(Y)^G$. Then $[F:K]=d$, and the generic $G$-torsor of $Y$ is the base change of that of $X$. By Theorem 3, the Amitsur subgroups are kernels of the corresponding multiplier-to-Brauer maps. Restriction gives the first inclusion. If a class becomes zero over $F$, corestriction shows that $d$ times its original Brauer class is zero over $K$, proving the second.

**Corollary 6.** At every prime not dividing $d$, the primary parts of the two Amitsur subgroups agree. Along a connected-monodromy tower beginning at a linear base, a nonzero order-three Amitsur class can first appear only at a step whose degree is divisible by three. A nonzero order-two class can first appear only at an even-degree step.

For an intermediate base still satisfying the Albanese hypothesis, a prefix of degree prime to three therefore retains the lower threshold $90$ in the table above. This is a lower bound on a subsequent compression degree, not a claim that every such degree is divisible by $90$.

## 7. The Brauer obstruction is cheap to remove

There is a decisive limitation on using only the preceding subgroup to attack unrestricted one-variable towers.

**Lemma 7 (representation-degree bound on index).** Let $T$ be a $G$-torsor over a field $K\supset\mathbb C$ and let $\alpha\in H^2(G,\mathbb C^*)$. Write $\beta_T(\alpha)\in\operatorname{Br}(K)$ for the associated Brauer class. For every projective representation of multiplier $\alpha$ and dimension $n$,

$$\operatorname{ind}\beta_T(\alpha)\mid n.$$

**Proof.** Twist its projective space by $T$. The result is the Severi-Brauer variety of a central simple algebra of degree $n$ and Brauer class $\beta_T(\alpha)$, up to the harmless opposite-algebra convention. The index divides the degree. Taking the greatest common divisor over representations proves the assertion.

For $A_7$, Chapter 7 lists the faithful projective degrees in each central character. Their arithmetic gives:

| Multiplier order | Projective irreducible degrees | Greatest common divisor |
|:---|:---|:---|
| $2$ | $4,4,14,14,20,20,36$ | $2$ |
| $3$ | $6,15,15,21,21,24,24$ | $3$ |
| $6$ | $6,6,24,24,36$ | $6$ |

The order-three row can already be bounded using the characteristic-zero dimensions $6$ and $21$ supplied by ATLAS [2]. The remaining projective-degree inputs are the Chapter 7 character table; their construction is an upstream representation-theoretic input, not a consequence of the gcd check. The gcds and the sum of squared degrees $2520$ in each row were checked in integer arithmetic.

**Corollary 8 (generic indices).** For the generic $A_7$-torsor over an ordinary linear base, a multiplier of order $r\in\{2,3,6\}$ has Brauer index exactly $r$.

**Proof.** The base has zero Amitsur subgroup, so its multiplier-to-Brauer map is injective. The period is therefore exactly $r$. Since period divides index, and Lemma 7 bounds the index by the corresponding gcd, equality follows.

**Corollary 9 (splitting without losing $A_7$).** The generic order-two class has a quadratic splitting extension; the generic order-three class has a cubic splitting extension. Both classes can be killed successively by extensions of degrees at most two and three. The resulting extension has degree at most six and preserves connected full $A_7$ monodromy.

**Proof.** A central division algebra of index $r$ has a separable maximal subfield of degree $r$, and that subfield splits it [3]. After splitting the two-part, the three-part still has index at most three, so the second step exists. To check monodromy, let $L/K$ be the original $A_7$ extension and let $F/K$ be any of these extensions of degree at most six. The group $A_7$ has no proper subgroup of index at most six: an action on fewer than seven cosets would, by simplicity, give an impossible injection into $S_6$. Thus $L\cap F=K$, and the Galois extension $L/K$ is linearly disjoint from $F/K$.

Quadratic and cubic equations are solvable by radical towers. Consequently the **Schur-Brauer obstruction alone cannot rule out towers of algebraic functions of one variable**. This does not construct a curve compression after those small extensions, and says nothing about the Albanese of the resulting total spaces. Those are separate geometric questions. Nor does it contradict $\mu(A_7)=90$: a cubic splitting step followed by a degree-60 compression would have total degree $180$, not $60$.

## 8. The surviving obstruction when the Albanese overlaps

For a general correspondence $Z$, record the map

$$a_Z:X\longrightarrow\operatorname{Pic}^e(C),\qquad x\longmapsto[Z_x].$$

It is initially defined on a dense open. Its rational map to the Picard torsor extends over smooth projective $X$, as a map to a torsor under an abelian variety. Its linear part is an equivariant homomorphism

$$u_Z:\operatorname{Alb}X\longrightarrow\operatorname{Jac}C.$$

Choosing an origin writes the map as a translate of $u_Z\circ\operatorname{alb}_X$. Such origins need not be $G$-fixed, so the correct equivariant object is the **affine map between Albanese and Picard torsors**, not just its linear part.

The effective-divisor condition says that $a_Z(X)$ lies in the Brill-Noether locus $W_e(C)$, the locus of degree-$e$ classes with a section. The correspondence also gives a rational lift to $\operatorname{Sym}^e C$, through the Abel map. If $u_Z=0$, this reduces to a single invariant linear system and Theorem 4 applies. If $u_Z\ne0$, it need not do so.

The diagonal $X=C$, $Z=\Delta_C$, has $e=1$ and $u_Z=\mathrm{id}$, even when $C$ has large gonality. Thus no bound depending only on gonality or the Amitsur subgroup can hold for arbitrary overlapping bases. The next precise target is to control the affine maps whose images lie in $W_e(C)$ **and admit the required lift to the symmetric power**, for the actual intermediate bases of an allowed tower. Sharing an abstract $G$-constituent in $H^1$ is only a warning that this can happen; it does not produce the required Hodge homomorphism or effective correspondence.

## 9. Two short repairs from the preliminary audit

**Theorem 7.4, coprimality.** The numerical root-distance test in `twisted_rr.py` is unnecessary once its stated Molien coefficients are accepted. Both $f_{14}$ and $f_{18}$ are invariants of the perfect group $2.A_7$. Their entire gcd is therefore invariant: the group preserves its one-dimensional span, and has no nontrivial characters. A nonconstant gcd of degree at most $14$ can have degree only $8$, $12$, or $14$. Dividing $f_{14}$ in the first two cases would give an invariant of degree $6$ or $2$; dividing $f_{18}$ in the last case would give one of degree $4$. All are absent. Hence the gcd is one.

The 210 lines are distinct as well. In a unitary realisation, an involution lift has two orthogonal two-dimensional eigenspaces with eigenvalues $i$ and $-i$. Either eigenspace determines the other and hence the projective involution. A line shared by two such lifts would force the same involution in $A_7$. Thus the Bezout remainder really is $14\cdot18-210=42$. This repairs the coprimality/line-count subargument; it is not a re-verification of the Molien table or local-character computations.

**Proposition 7.7, step 1.** Parity forces odd vanishing orders for anti-invariant sections; it does not by itself prove the particular orders $\{1,3\}$. An additional jet-separation or embedding argument is required to conclude a simple common zero. A proved immersion at a fixed point would suffice: if both anti-invariant sections vanished to order at least three, then all projective coordinate derivatives would vanish there, since the invariant-coordinate ratios are even. The numerical embedding evidence in Proposition 7.8 is not yet that proof.

Also, the label “universal central extension” in Proposition 7.1 needs its geometric meaning specified. The triangle group $\Delta(2,4,7)$ has abelianisation $C_2$ and is not perfect, so it has no universal central extension in the usual group-theoretic sense. A preimage in the universal covering group of $\mathrm{PSL}_2(\mathbb R)$ is a different, legitimate object. This terminological point does not decide the degree formula or its orientation.

## 10. What remains unaudited

The pivot to the framework above follows the requested priority for structural understanding. The Jacobian-flux/residual/radial-bound certificate, the remaining Chapter 4 arguments, the degree formula of Proposition 7.1, the exact projective-equation and orbit computations in Chapter 7, and the rank-three exclusion in Section 5.4 have **not** received a complete independent audit in this note. Reading their text or saved outputs is not certification. No new value of $\gamma(A_7)$, stable-reduction graph, or nontorsion arithmetic point is claimed.

The concrete chapter edit is to replace Theorem 5.11's open-base linearisation criterion by Theorem 1 and Corollary 2 above, and to state fixed-target divisibility with the base hypothesis attached. The original linear-base results and the particular punctured-quadric corollary are not refuted by this audit. The general stable formula and the splitting limitation supply a more precise starting point for the tower problem.

## References and provenance

1. B. Hassett and Y. Tschinkel, *Torsors and stable equivariant birational geometry*, Nagoya Mathematical Journal 250 (2023), 275-297. Sections 3.5-3.6 give the Amitsur subgroup and a criterion involving the unit sequence of a Picard-trivial open; the no-name lemma is also discussed. [Author PDF](https://www.math.brown.edu/bhassett/papers/gtorsor/gtorsor3.pdf).
2. ATLAS of Finite Group Representations, [alternating group A7 and its covers](https://brauer.maths.qmul.ac.uk/Atlas/v3/alt/A7/). For the full projective-degree lists used here, see also Chapter 7.1 and `twisted_rr_output.txt` at the repository revision identified above.
3. The Stacks Project, [Section 11.8: Splitting fields](https://stacks.math.columbia.edu/tag/074X), especially Lemma 11.8.3, Proposition 11.8.5, and the central-simple-algebra structure theorem.
4. [Repository source at the audited revision](https://github.com/JimmysanUWU/ARAGACAS/tree/ea25d8ae94c179cf5315b4e68f1af0f33a7d2bbd/hilbert13/ladder). The new proofs in this note were derived directly; the gonality certificate and projective character constructions remain upstream inputs where explicitly indicated.
