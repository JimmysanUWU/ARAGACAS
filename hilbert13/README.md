# Hilbert's 13th problem: single superpositions

Paper proofs and Lean 4 / Mathlib formalizations of three small results that show why the
choice of function class is the entire content of Hilbert's 13th problem.

All proofs are in [`Hilbert13/Superposition.lean`](Hilbert13/Superposition.lean). They
compile with no `sorry`, and `#print axioms` reports only Lean's standard axioms
(`propext`, `Classical.choice`, `Quot.sound`).

```sh
cd hilbert13
lake exe cache get   # prebuilt Mathlib
lake build
```

## Background

Hilbert (1900) conjectured that the roots of the general septic

$$x^7 + a x^3 + b x^2 + c x + 1 = 0,$$

viewed as a function of $(a, b, c)$, cannot be written as a finite composition of
continuous functions of two variables.

- **Continuous version: false (Kolmogorov 1956, Arnold 1957).** Every continuous
  $f : [0,1]^n \to \mathbb{R}$ can be written
  $$f(x_1,\dots,x_n) = \sum_{q=0}^{2n} \Phi_q\Big(\sum_{p=1}^{n} \varphi_{q,p}(x_p)\Big)$$
  with continuous one-variable $\Phi_q, \varphi_{q,p}$. So addition plus one-variable
  continuous functions already generate everything.
- **Algebraic version: open.** Here the building blocks must be algebraic functions. It is
  measured by the *resolvent degree* $\mathrm{RD}(n)$ (Brauer; Farb–Wolfson). Classically
  $\mathrm{RD}(5) = 1$ (Bring radical) and $\mathrm{RD}(7) \le 3$, and Hilbert's question
  is whether $\mathrm{RD}(7) = 3$. Remarkably, it is not known whether $\mathrm{RD}(n) > 1$
  for *any* $n$.

The results below explain one feature of the Kolmogorov–Arnold formula: why it needs
several outer terms $\Phi_q$ and continuous inner functions, and why dropping regularity
makes the problem trivial.

## Theorem 1: without regularity, one term suffices

**Theorem** (`exists_single_superposition`). Let $n \ge 1$. For every function
$f : \mathbb{R}^n \to \mathbb{R}$ there exist functions $g, \varphi_1, \dots, \varphi_n :
\mathbb{R} \to \mathbb{R}$ with
$$f(x_1, \dots, x_n) = g\big(\varphi_1(x_1) + \dots + \varphi_n(x_n)\big).$$

*Proof.* View $\mathbb{R}$ and $\mathbb{R}^n$ as $\mathbb{Q}$-vector spaces. Since
$|\mathbb{Q}| = \aleph_0 < \mathfrak{c}$, the dimension of each equals its cardinality:
$\dim_{\mathbb{Q}} \mathbb{R} = |\mathbb{R}| = \mathfrak c$, and
$\dim_{\mathbb{Q}} \mathbb{R}^n = \mathfrak c^n = \mathfrak c$. Vector spaces of equal
dimension are isomorphic, so there is a $\mathbb{Q}$-linear bijection
$e : \mathbb{R}^n \to \mathbb{R}$. Put $\varphi_i(t) = e(t\,\mathbf{e}_i)$, where
$\mathbf{e}_i$ is the $i$-th basis vector, and $g = f \circ e^{-1}$. By additivity,
$$\sum_i \varphi_i(x_i) = e\Big(\sum_i x_i \mathbf{e}_i\Big) = e(x),$$
so $g\big(\sum_i \varphi_i(x_i)\big) = f(e^{-1}(e(x))) = f(x)$. $\blacksquare$

The map $e$ comes from a Hamel basis, so it uses the axiom of choice and is wildly
discontinuous. Theorem 2 shows that this is unavoidable.

## Theorem 2: with continuous inner functions, one term fails for $xy$

**Theorem** (`mul_ne_single_superposition`). There are no functions
$g, \varphi, \psi : \mathbb{R} \to \mathbb{R}$ with $\varphi, \psi$ **continuous** such that
$$xy = g\big(\varphi(x) + \psi(y)\big) \quad \text{for all } x, y \in \mathbb{R}.$$
No assumption at all is made on the outer function $g$.

*Proof.* Suppose such $g, \varphi, \psi$ exist.

1. **$g$ vanishes on the line $x = 0$.** For all $y$: $g(\varphi(0) + \psi(y)) = 0 \cdot y = 0$.
2. **$\psi$ is injective.** If $\psi(a) = \psi(b)$ then
   $a = 1 \cdot a = g(\varphi(1) + \psi(a)) = g(\varphi(1) + \psi(b)) = b$.
3. **$\psi(\tfrac12)$ is strictly between $\psi(0)$ and $\psi(1)$.** A continuous injective
   function on $\mathbb{R}$ is strictly monotone. Let $J$ be the open interval with
   endpoints $\psi(0)$ and $\psi(1)$; then $\psi(\tfrac12) \in J$.
4. **Perturb in $x$.** The map $t \mapsto \varphi(t) - \varphi(0) + \psi(\tfrac12)$ is
   continuous and sends $0$ to $\psi(\tfrac12) \in J$. Since $J$ is open, there is some
   $t \neq 0$ with $v := \varphi(t) - \varphi(0) + \psi(\tfrac12) \in J$.
5. **Intermediate value theorem.** $\psi$ is continuous on $[0,1]$ and $v$ lies between
   $\psi(0)$ and $\psi(1)$, so $v = \psi(y')$ for some $y' \in [0,1]$. Rearranging,
   $\varphi(t) + \psi(\tfrac12) = \varphi(0) + \psi(y')$.
6. **Contradiction.** Using step 1,
   $$\tfrac{t}{2} = g\big(\varphi(t) + \psi(\tfrac12)\big) = g\big(\varphi(0) + \psi(y')\big) = 0,$$
   so $t = 0$, contradicting $t \neq 0$. $\blacksquare$

**What is really going on.** Level sets of $s(x,y) = \varphi(x) + \psi(y)$ are the only
thing $g$ can distinguish. The zero set of $xy$ is the union of the two axes, and it
contains the "cross" at the origin. But the level set of $s$ through the vertical axis
$\{0\} \times \mathbb{R}$ also sweeps into nearby vertical lines $x = t$, because
$\psi$ fills an interval and $\varphi$ moves continuously. Nearby points where $xy \neq 0$
would then have to share $g$-values with points where $xy = 0$.

## Theorem 3: two continuous terms suffice for $xy$

**Theorem** (`mul_two_superposition`).
$$xy = \tfrac14 (x + y)^2 - \tfrac14 (x - y)^2.$$

This has the Kolmogorov–Arnold shape $\sum_q \Phi_q(\varphi_q(x) + \psi_q(y))$ with two
terms, where $\Phi_1(t) = t^2/4$, $\Phi_2(t) = -t^2/4$, and all inner functions linear.

## Summary

| Inner functions | Outer function | Terms | Represents every $f$? |
|---|---|---|---|
| arbitrary | arbitrary | 1 | **yes** (Theorem 1) |
| continuous | arbitrary | 1 | **no**, already $xy$ fails (Theorem 2) |
| continuous | continuous | $2n+1$ | **yes** (Kolmogorov–Arnold; not formalized here) |
| algebraic | algebraic | any | **open** (Hilbert 13, algebraic form) |

## References

- A. N. Kolmogorov, *On the representation of continuous functions of many variables by
  superposition of continuous functions of one variable and addition*, Dokl. Akad. Nauk
  SSSR 114 (1957).
- V. I. Arnold, *On functions of three variables*, Dokl. Akad. Nauk SSSR 114 (1957).
- B. Farb, J. Wolfson, *Resolvent degree, Hilbert's 13th Problem and geometry*,
  L'Enseignement Math. 65 (2019), [arXiv:1803.04063](https://arxiv.org/abs/1803.04063).

## The $A_7$ gonality ladder (`ladder/`)

The second project concerns curves with a faithful $A_7$-action. Its main case is the $(2,4,7)$ curves, of genus 136, the minimum.
Full statements, proofs and certificates are in [`ladder/README.md`](ladder/README.md):
- **$\operatorname{gon}(C)\ge25$** for every $(2,4,7)$ curve, so $\operatorname{gon}(C/\langle\tau\rangle)\ge13$ for every involution $\tau$.
  - The route is a certified spectral gap with Li–Yau, $\lambda_1\,\mathrm{Area}\le8\pi\deg$ ([Ch. 2](ladder/2_SPECTRAL.md)).
  - For two classes it adds *harmonic Hersch*: the degree is a cubic form that vanishes on the first eigenspace by representation theory ([Ch. 3](ladder/3_HARMONIC_HERSCH.md)).
- **$\gamma(A_7)\ge25$**: every faithful $A_7$-curve has gonality at least 25 ([Ch. 4](ladder/4_LARGE_GENUS.md)).
- **$\operatorname{gon}(C)\le42$** for the genus-136 curves, from a degree-60 model in $\mathbb P^5$ built from an $A_7$-invariant line bundle whose
  symmetry is only projective (Schur multiplier $\mathbb Z/6$, [Ch. 7](ladder/7_TWISTED.md)). So $25\le\operatorname{gon}(C)\le42$.
- **Accessory irrationalities** ([Ch. 5](ladder/5_ACCESSORY.md)).
  - $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$, and this is sharp: $a(A_7)=60$. The published bound is 6.
  - Keeping connected full $A_7$-monodromy, the exact threshold is $\mu(A_7)=90$.

The logical core of the spectral certificate is checked in Lean in
[`Hilbert13/SpectralCertificate.lean`](Hilbert13/SpectralCertificate.lean):
- the Crouzeix–Raviart eigenvalue lower bound;
- positive definiteness from a perturbed Cholesky factor;
- the final arithmetic.
