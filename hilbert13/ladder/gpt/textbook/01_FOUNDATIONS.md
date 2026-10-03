# 1. Function classes and curve invariants

## 1.1 Single superpositions

**Definition.** A single additive superposition of $n$ variables has the form
$$f(x_1,\ldots,x_n)=g\left(\sum_{i=1}^n\varphi_i(x_i)\right).$$
Its expressive power depends on the allowed function class.

**Theorem 1.1 (arbitrary functions; P).** For $n\ge1$, every function $f:\mathbb R^n\to\mathbb R$ is a single additive superposition with arbitrary $g,\varphi_i:\mathbb R\to\mathbb R$.

*Proof.* The $\mathbb Q$-vector spaces $\mathbb R$ and $\mathbb R^n$ both have dimension $\mathfrak c$, the cardinality of the continuum. Indeed, a vector space over a countable field with uncountable cardinality has dimension equal to that cardinality, and $\mathfrak c^n=\mathfrak c$. Choose a $\mathbb Q$-linear isomorphism $e:\mathbb R^n\to\mathbb R$. Set $\varphi_i(t)=e(t\mathbf e_i)$ and $g=f\circ e^{-1}$. Additivity gives $\sum_i\varphi_i(x_i)=e(x)$, proving the identity. The choice of $e$ uses a Hamel basis. $\square$

**Theorem 1.2 (continuous inner functions; P).** There are no functions $g,\varphi,\psi:\mathbb R\to\mathbb R$, with $\varphi$ and $\psi$ continuous, satisfying
$$xy=g(\varphi(x)+\psi(y))\qquad(x,y\in\mathbb R).$$
No regularity hypothesis on $g$ is necessary.

*Proof.* Setting $x=1$ shows that $\psi$ is injective. It is therefore strictly monotone. Thus $\psi(1/2)$ lies in the open interval $J$ between $\psi(0)$ and $\psi(1)$. By continuity of $\varphi$, there exists $t\ne0$ sufficiently close to zero such that $v=\varphi(t)-\varphi(0)+\psi(1/2)\in J$. The intermediate value theorem gives $y'\in[0,1]$ with $\psi(y')=v$. Consequently
$$\frac t2=g(\varphi(t)+\psi(1/2))=g(\varphi(0)+\psi(y'))=0,$$
a contradiction. $\square$

**Proposition 1.3 (two terms; P).** Two continuous superpositions represent multiplication:
$$xy=\frac{(x+y)^2}{4}-\frac{(x-y)^2}{4}.$$
The proof is expansion. These three results, including the finite-dimensional cardinality argument, are the statements formalized in `Hilbert13/Superposition.lean`.

The Kolmogorov--Arnold theorem is a separate classical result: every continuous $f:[0,1]^n\to\mathbb R$ can be expressed as
$$f(x)=\sum_{q=0}^{2n}\Phi_q\left(\sum_{p=1}^n\varphi_{q,p}(x_p)\right)$$
with continuous one-variable functions. It is not formalized by the three propositions above.

## 1.2 Resolvent degree and the algebraic problem

**Definition.** A finite extension $M/K$ of fields containing $\mathbb C$ has resolvent degree at most $d$ if it embeds in the last field of a finite tower
$$K=F_0\subset F_1\subset\cdots\subset F_s,$$
where each step is obtained by base change from a finite extension $E_i'/E_i$ with $E_i\subset F_{i-1}$ and $\operatorname{trdeg}_{\mathbb C}E_i\le d$. The degrees and the number of steps are unrestricted.

For a finite group $G$, its resolvent degree is measured on a generic torsor, for example $\mathbb C(V)/\mathbb C(V)^G$ for a faithful linear representation $V$. The algebraic form of Hilbert's thirteenth problem concerns $\mathrm{RD}(S_7)=3$. The known upper bound is three; even the lower bound $\mathrm{RD}(A_7)>1$ remains open in this project. The composition-factor reduction relating $S_7$ and $A_7$ is an external input from resolvent-degree theory.

An obstruction to a single finite accessory of bounded degree does not exclude all such towers. In particular, the accessory thresholds proved later do not prove a resolvent-degree lower bound greater than one.

## 1.3 Pencils and equivariant line bundles

**Definition.** The gonality $\operatorname{gon}(C)$ is the least degree of a nonconstant morphism $C\to\mathbb P^1$. A pencil consists of a two-dimensional subspace $U\subset H^0(C,L)$. Removing its common zero divisor $B_U$ gives a base-point-free pencil of degree $\deg L-\deg B_U$.

Pencils are identified as maps up to a change of coordinate on $\mathbb P^1$, equivalently by the rational subfield $\mathbb C(f)\subset\mathbb C(C)$. A gonal line bundle has exactly two sections: if $h^0(L)\ge3$, then $h^0(L(-p))\ge2$, and removing its fixed part gives a smaller pencil.

An *invariant* line-bundle class satisfies $g^*L\simeq L$ for every $g\in G$. A *linearization* chooses these isomorphisms compatibly with multiplication in $G$. On a projective connected variety the scalar failure of compatibility defines
$$m(L)\in H^2(G,\mathbb C^*).$$
The class vanishes exactly when a linearization exists. A projective section representation is sufficient for an equivariant map to projective space; a genuine section representation requires a linearization.

Write
$$\gamma(G)=\min_C\operatorname{gon}(C),\qquad
\mu(G)=\min_{C,L}\deg L,$$
where $C$ ranges over faithful $G$-curves and, in the second minimum, $L$ is a linearized base-point-free bundle with $h^0(L)\ge2$. On a specified curve use $\mu_C(G)$. The analogous minimum over invariant moving classes on $C$ is denoted $\widetilde\mu(C)$. These are distinct minimization problems.

## 1.4 Genus inequalities

**Castelnuovo--Severi.** If $C\to C_i$ have degrees $d_i$ and jointly generate $\mathbb C(C)$, then
$$g(C)\le d_1g(C_1)+d_2g(C_2)+(d_1-1)(d_2-1).$$
In particular, a birational pair of degree-$d$ pencils gives a bidegree-$(d,d)$ image in $\mathbb P^1\times\mathbb P^1$ with geometric genus at most $(d-1)^2$.

**Castelnuovo's bound.** For an integral nondegenerate curve $Y\subset\mathbb P^r$, $r\ge2$, of degree $n$, write $n-1=q(r-1)+s$ with $0\le s<r-1$. Then
$$p_a(Y)\le\pi(n,r):=(r-1)\binom q2+qs.$$
This bounds the arithmetic genus, including for singular images. The normalization has genus $p_a(Y)-\sum_y\delta_y$. Consequently a birational product map need not be an embedding, even when its image satisfies a genus bound.

**Riemann--Hurwitz.** A faithful group action with quotient genus $h$ and cyclic inertia orders $e_1,\ldots,e_r$ satisfies
$$2g(C)-2=|G|\left(2h-2+\sum_i(1-1/e_i)\right).$$
For $G=A_7$ the nontrivial element orders are $2,3,4,5,6,7$. Each $1260/e_i$ is divisible by three, so every faithful $A_7$-curve satisfies $g(C)\equiv1\pmod3$.

## 1.5 Notation

Throughout, $G=A_7$, $|G|=2520$, and $C$ is a $(2,4,7)$ cover when no other signature is specified. Its reduced branch fibres are $D_2,D_4,D_7$, with degrees $1260,630,360$. For an involution $\tau$, put $D=C/\langle\tau\rangle$ and $R_\tau=\operatorname{Fix}_C(\tau)$ as a reduced divisor. Use $T_0=D/(C_G(\tau)/\langle\tau\rangle)$ for the genus-three quotient, reserving $T_{\rm tor}=D_2-2D_4$ for a torsion line-bundle class.

The ordinary irreducibles are $1,6,10,\overline{10},14_a,14_b,15,21,35$. The label $14_a$ always means the $(5,2)$ constituent of the permutation representation on two-element subsets. The exceptional projective six of $3.A_7$ is denoted $V$; it is different from the ordinary standard six of $A_7$. Let $S=\operatorname{Sym}V$, the coordinate ring of $\mathbb P(V^*)$.

The four generating-triple indices $0,1,12,14$ are array indices in `triples_data.py`, not intrinsic invariants. Permutations in the programs use letters $0,\ldots,6$ and composition $pq=p\circ q$; mathematical cycle notation here uses $1,\ldots,7$.
