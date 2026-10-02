# Chapter 6. Side results, alternative routes, and negative results

These results are not on the critical path of Theorems 3.1, 4.1 and 5.1. They are kept because each one either proves something
independently or explains why a natural route fails.

## 6.1 An algebraic proof of $\operatorname{gon}(D)\ge10$ (Sol/Astra) [P][X]

This is GPT's *Ramification transport* (28 Sep 2026), re-derived. The finite inputs are checked in `review_checks.py`.

*Proof.* Suppose $D=C/\langle\tau\rangle$ has a degree-9 map. Pull it back to $f:C\to\mathbb P^1$ of degree 18, with $f\circ\tau=f$ and $L=f^*\mathcal O(1)$.
1. **Determinant.** Let $\mu\ne\tau$ be an involution commuting with $\tau$.
   - The section $s_0\otimes\mu^*s_1-s_1\otimes\mu^*s_0$ of $L\otimes\mu^*L$ vanishes on $\mathrm{Fix}(\mu)\cup\mathrm{Fix}(\tau\mu)$. These are $18+18$ distinct points, since a common one would have stabiliser $V_4$.
   - It is not identically zero, since otherwise $4\mid18$.
   - Hence $[L]+[\mu^*L]=[R_\mu+R_{\tau\mu}]$.
2. **Average.** With $M=\sum_{h\in C_G(\tau)}h^*L$, summing gives $2M=24[R_\mu+R_{\tau\mu}]$. So $T:=24[R_\tau]+2M=24[R_V]$ for **every**
   Klein four-group $V\ni\tau$. Nothing is divided, so torsion is harmless.
3. **Amalgamation.** $N_G(V)$ preserves $R_V$. Since $\langle N(V_0),N(V_1)\rangle=A_7$ (orders 72 and 24), $T$ is $A_7$-invariant, of degree $24\cdot54=1296$.
4. **Descent.** A subgroup of order 5 acts freely, so an invariant class descends ($H^2(C_5,\mathbb C^\times)=0$). This forces $5\mid1296$, which is false.

Degrees $\le8$ are excluded by the same determinant, because $18>2\cdot8$. $\square$

**Remarks.**
- The proof uses the full $A_7$-action (normalisers that do not centralise $\tau$, and a free $C_5$). So it avoids the audit barrier of Theorem 1.9.
- **Lattice form.** Step 4 is a degree-lattice violation (Lemma 5.4 for invariant classes): $2520/(\exp M(A_7)\cdot28)=15$, and $15\nmid1296$.
  The same mechanism works for $L_2(13)$ with signature $(2,3,7)$, where $13\nmid216$.

## 6.2 Pencil orbits [P]

**Lemma 6.1.** Let $m<60$. A degree-$m$ pencil on a faithful $A_7$-curve has an $A_7$-orbit of size at least 35, and at least 42 if $3\nmid m$.

*Proof.*
1. Let $K$ be the stabiliser of the pencil, and $K_0=\{k:f\circ k=f\}$. Then $|K_0|\mid m$.
2. $K/K_0\subset\mathrm{PGL}_2$ is cyclic, dihedral, $A_4$, $S_4$ or $A_5$.
3. The subgroups of order $>72$ are $A_7,A_6,L_2(7),S_5$. Each forces $K_0\supseteq$ a simple group of order $\ge60$, which is impossible.
4. The only subgroup of order 72 is $(A_4\times3){:}2$. Its polyhedral quotients all need $3\mid m$. $\square$

Lemma 4.3 is the cruder version actually needed. Lemma 6.1 gives at least 33 third pencils for a birational pair.

## 6.3 The conformal route and its ceiling

Li–Yau holds for any conformal metric. One may therefore try to maximise $\lambda_1\cdot\mathrm{Area}$ over $G$-invariant conformal metrics $hg$.

**Proposition 6.2 (non-criticality) [P].** Let $C/G$ be a triangle orbifold. Then no $G$-invariant positive semidefinite form $Q$ on $E_1$ has
$Q(\varphi_1,\dots)\equiv1$. Hence the hyperbolic metric is never critical, and $\Lambda^G>\lambda_1A$ strictly.

*Proof.* Write $Q=\sum\psi_i^2$, and consider $\Phi=(\psi_i):C\to S^N$.
1. $|\nabla\Phi|^2\equiv\lambda_1$, and $\Phi$ is harmonic.
2. The Hopf differential is an invariant holomorphic quadratic differential. Triangle orbifolds carry none, so $\Phi$ is a conformal minimal immersion.
3. Its induced curvature is $-2/\lambda_1<0$. Bryant (Trans. AMS 290, 1985) shows there is no minimal surface of constant negative curvature in $S^n$. $\square$

**Proposition 6.3 (ceiling) [P].** Let $\bar F$ be the $G$-average of $\sum\varphi_i^2$, normalised to mean 1. Then $\lambda_1(hg)\,\mathrm{Area}(hg)\le\lambda_1A/\min\bar F$.

*Proof.* Each $\varphi\in E_1$ is $h\,dA$-orthogonal to the constants, by Schur. Take the mediant of the Rayleigh quotients. $\square$

**Numerics [N]** (`archive/scripts/conformal.py`).
- Classes 0, 1: $\bar F\in[0.998,1.001]$, a gain $\le0.2\%$. But $+2.7\%$ is needed.
- Classes 12, 14: gain $\le0.6\%$.

**Explanation.** The invariant part of $\sum\varphi_i^2$ is nearly constant: it is band-limited near $2\lambda_1$, while the first invariant eigenvalue is about $10.6$. The route is closed, and its replacement is the degree (Chapter 3).

## 6.4 Topological Hersch (precursor of Chapter 3) [P][N]

This was the first inequality to use Lemma 3.3. For balanced holomorphic $x=y+z$ of energy $\lambda_1A+\varepsilon$, it bounds the degree gain by
$6\kappa\sqrt\varepsilon$ plus an $E(z)/2$ term. GPT's harmonic Hersch replaces this with $4\kappa\sqrt\varepsilon$ and no $E(z)$ term, by using
$x\times dx={*dx}$.

The numerics showed a wide margin, but two constants were never certified. Theorem 3.4 makes this obsolete.

Archived: `archive/docs/FRAMEWORK_TOPOLOGICAL_HERSCH.md`, `archive/scripts/topo_hersch*.py`.

**Quadric gap [N].** $\min\frac1A\int(1-|y|)^2$ over $y\in E_1\otimes\mathbb R^3$ with $\|y\|^2=A$ is about $0.016$. This alone gives only $\operatorname{gon}\ge23.7$:
the obstruction is the degree, not the shape.

## 6.5 The equivariant Plücker test (negative) [X]

**The test.** For a $G$-stable series $W\subseteq H^0(L)$ on a three-point $A_7$-curve, the eigenvalues of the stabiliser at a branch point fix the
vanishing orders mod $e$. Plücker's formula then requires
$$(r+1)(d+r(g-1))-\sum_i\tfrac{2520}{e_i}w_{\min}(i)\ \in\ 2520\,\mathbb Z_{\ge0}.$$

**Outcome.** It passes in every case tested, so it excludes nothing beyond Chapter 5 (`archive/scripts/equivariant_plucker.py`).

## 6.6 The $Q_2$ group [X]

$C_{A_7}((16)(23))$ has the following structure (`frontier_checks.py`):
- order 24, with centre $\langle(16)(23)\rangle$;
- element orders $1{:}1,\ 2{:}9,\ 3{:}2,\ 4{:}6,\ 6{:}6$;
- a normal Sylow 3-subgroup and a Sylow 2-subgroup $D_8$.

So it is $C_3\rtimes D_8\cong D_8\times_{C_2}S_3$, not $S_4$.

## 6.7 Open arithmetic: the points $P_E$ and $P_S$ [X][?]

Proposition 1.8 reduces $(\star)$ to the following question: is either of two explicit points of infinite order?
- $P_E\in E=C/L_2(5)$, of genus 1, with $\mathrm{Jac}\sim E_2$. Here $L_2(5)\cong A_5$ acts transitively on 6 letters.
- $P_S\in\mathrm{Jac}(C/(3^2{:}4))$, of genus 2, with $\mathrm{Jac}\sim S$.

Facts:
- Both are images of the Klein difference $D_\delta=\mathrm{Fix}(v_2)+\mathrm{Fix}(v_3)-\mathrm{Fix}(w_2)-\mathrm{Fix}(w_3)$. Its 21- and 35-components are nonzero, and each generates the unique copy.
- A 5-adic argument shows $\delta\ne0$. So if both points are torsion, $\delta$ is torsion of order divisible by 5.
- $E$ is a 6-sheeted cover of $C/A_6\cong\mathbb P^1$. Its map to $C/A_7$ is the degree-7 Belyi map with passport $[2^21^3,\,4\,2\,1,\,7]$. So the question is a finite computation once a model is known.

**Lens.** $\Delta(2,4,7)$ is arithmetic (Takeuchi), but $\ker(\Delta\to A_7)$ is non-congruence. So this is a Manin–Drinfeld-type question with no Hecke operators.
Chapters 2–3 make it unnecessary for gonality. It remains interesting in its own right.
