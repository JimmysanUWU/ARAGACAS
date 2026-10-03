# Chapter 6. Side results, alternative routes, negative results

Each result here either proves something independently of the main line or explains why a natural route fails.

## 6.1 An algebraic proof of $\operatorname{gon}(D)\ge10$ (Astra) [P][X]

**Idea.**
- A nine-pencil on $D$, pulled back to $C$, saturates the determinant against each commuting involution: there are as many forced zeros as the degree allows. Each saturation is an exact Picard identity.
- Averaging over the 24-element centraliser, not over all of $A_7$, removes the choice of pencil without dividing by 5.
- Two Klein four-groups through $\tau$ have normalisers that generate $A_7$, so the averaged class is $A_7$-invariant.
- Its degree, 1296, is not divisible by 5, while a free $C_5$-action forces divisibility by 5.

This is Astra's *ramification transport* (GPT; `gpt/README.md`), re-derived; the finite inputs are in `side_checks.py` (R).

*Proof.* Suppose $D$ has a degree-9 map. Pull it back to $f:C\to\mathbb P^1$ of degree 18, with $f\circ\tau=f$ and $L=f^*\mathcal O(1)$.
1. **Determinant.** Let $\mu\ne\tau$ be an involution commuting with $\tau$.
   - The section $s_0\otimes\mu^*s_1-s_1\otimes\mu^*s_0$ of $L\otimes\mu^*L$ vanishes on $\mathrm{Fix}(\mu)\cup\mathrm{Fix}(\tau\mu)$: 36 distinct points, since no point has stabiliser $V_4$.
   - It is not identically zero, since $4\nmid18$.
   - So $[L]+[\mu^*L]=[R_\mu+R_{\tau\mu}]$.
2. **Average.** With $M=\sum_{h\in C_G(\tau)}h^*L$, we get $T:=24[R_\tau]+2M=24[R_V]$ for every Klein four-group $V\ni\tau$. Nothing is divided, so torsion is harmless.
3. **Amalgamation.** $N_G(V)$ preserves $R_V$, and $\langle N(V_0),N(V_1)\rangle=A_7$ (orders 72 and 24). So $T$ is invariant, of degree $24\cdot54=1296$.
4. **Contradiction.** $1296\notin15\mathbb Z$, the degree lattice of invariant classes (Corollary 7.2). The original argument says the same through a free $C_5$: $5\nmid1296$.

Degrees $d\le8$ fail at step 1. The determinant would have 36 zeros on a bundle of degree $4d\le32$, so it vanishes identically and $f\circ\mu=f$ for every $\mu$. Then $12\mid d$, as in Theorem 1.4. $\square$

The proof uses the full $A_7$-action, so it avoids the audit barrier (Theorem 1.7). The same mechanism works for $L_2(13)$ with signature $(2,3,7)$, where $13\nmid216$: the genus-6 involution quotients of the genus-14 Hurwitz curves have gonality exactly 4 (DAY 2 Cor. 4.4).

## 6.2 Pencil orbits [P]

**Lemma 6.1.** For $m<60$, a degree-$m$ pencil on a faithful $A_7$-curve has an $A_7$-orbit of size $\ge35$, and $\ge42$ if $3\nmid m$.

*Proof.* Let $K$ stabilise the pencil and $K_0=\{k:f\circ k=f\}$. Then:
- $|K_0|$ divides $m$;
- $K/K_0\subset\mathrm{PGL}_2$ is cyclic, dihedral, $A_4$, $S_4$ or $A_5$.

The subgroups of order $>72$ ($A_7,A_6,L_2(7),S_5$) would force a simple $K_0$ of order $\ge60$. The one of order 72, $(A_4\times3){:}2$, needs $3\mid m$. $\square$

$(K,K_0)$ is the pair $(K,N)$ of Proposition 7.11. Lemma 4.3 is the cruder version used in Chapter 4.

## 6.3 The conformal route and its ceiling [P][N]

Li–Yau holds for every conformal metric. So one can try to maximise $\lambda_1\cdot\mathrm{Area}$ over $G$-invariant conformal metrics.

**Proposition 6.2 (non-criticality).** On a triangle orbifold, no $G$-invariant positive form $Q$ on $E_1$ has $Q(\varphi_1,\dots)\equiv1$. So the hyperbolic metric is never critical.

*Proof.* Write $Q=\sum\psi_i^2$, and let $\Phi=(\psi_i):C\to S^N$. Then:
1. $\Phi$ is harmonic with $|\nabla\Phi|^2\equiv\lambda_1$.
2. The Hopf differential is an invariant holomorphic quadratic differential, and a triangle orbifold has none. So $\Phi$ is a conformal minimal immersion, of curvature $-2/\lambda_1<0$.
3. Bryant (1985) shows that no such immersion exists in $S^n$. $\square$

**Proposition 6.3 (ceiling).** Let $\bar F$ be the $G$-average of $\sum\varphi_i^2$, normalised to mean 1. Then $\lambda_1\mathrm{Area}\le\lambda_1A/\min\bar F$ for every invariant conformal metric. (Each $\varphi\in E_1$ is orthogonal to the constants by Schur; take the mediant of Rayleigh quotients.)

**Numerics** (`conformal.py`, in git history at b148338). $\bar F\in[0.998,1.001]$, so the gain is $\le0.2\%$ for classes 0, 1 and $\le0.6\%$ for 12, 14. Classes 0, 1 need $+2.7\%$.

GPT adds that the gain is strict, and that an equivariant maximiser exists, possibly conical (Vinokurov; `gpt/README.md` §2.4). Neither changes the ceiling. So the route is closed. The reason is that $\sum\varphi_i^2$ is band-limited near $2\lambda_1$, while the first invariant eigenvalue is about 10.6. The degree replaces it (Chapter 3).

## 6.4 Negative tests

- **Quadric gap [N].** $\min\frac1A\int(1-|y|)^2\approx0.016$ over $y\in E_1\otimes\mathbb R^3$ with $\|y\|^2=A$. This gives only $\operatorname{gon}\ge23.7$: the obstruction is the degree, not the shape.
- **Equivariant Plücker [X]** (`equivariant_plucker.py`, in git history at b148338). Branch eigenvalues fix the vanishing orders mod $e$, and Plücker then requires
  $$(r+1)(d+r(g-1))-\sum_i\tfrac{2520}{e_i}w_{\min}(i)\in2520\,\mathbb Z_{\ge0}.$$
  It passes in every case tested, so it excludes nothing; Theorem 7.5 uses it as a check.
- **The $Q_2$ group [X]** (`side_checks.py` D1). $C_{A_7}((16)(23))$ has order 24, centre $\langle(16)(23)\rangle$, element orders $1{:}1,2{:}9,3{:}2,4{:}6,6{:}6$, a normal $C_3$ and Sylow-2 subgroup $D_8$. So it is $C_3\rtimes D_8$, not $S_4$.

## 6.5 Open arithmetic: $P_E$ and $P_S$ [X][?]

Proposition 1.6 reduces $(\star)$ to a question: does either of these points have infinite order?
- $P_E$ on the elliptic curve $E=C/L_2(5)$, where $L_2(5)\cong A_5$ is transitive on 6 letters.
- $P_S$ on $\mathrm{Jac}(C/(3^2{:}4))$, a Jacobian of genus 2.

**Facts.**
- Both points are images of the Klein difference $D_\delta=\mathrm{Fix}(v_2)+\mathrm{Fix}(v_3)-\mathrm{Fix}(w_2)-\mathrm{Fix}(w_3)$, whose 21- and 35-components each generate the unique copy.
- A 5-adic argument gives $\delta\ne0$. So if both points are torsion, $\delta$ is torsion of order divisible by 5.
- $E\to C/A_6\cong\mathbb P^1\to C/A_7$ involves the degree-7 Belyi map with passport $[2^21^3,\,4\,2\,1,\,7]$. GPT gives it explicitly over $\mathbb Q(\sqrt{21})$, with resolvent function fields for $E$ and $C/(3^2{:}4)$ (`gpt/README.md` §2.5; the map is checked in `side_checks.py` D6). What remains is Weierstrass models and the images of $D_\delta$.

**Lens.** $\Delta(2,4,7)$ is arithmetic (Takeuchi), but $\ker(\Delta\to A_7)$ is non-congruence. So this is a Manin–Drinfeld-type question with no Hecke operators available. It is not needed for gonality.
