# Chapter 7. Schur-twisted geometry: invariant line bundles, a model in $\mathbb P^5$, and $\operatorname{gon}(C)\le42$

**The idea.** An $A_7$-*invariant* line bundle need not be $A_7$-*linearised*. The obstruction, its *Mumford class*, lies in
$H^2(A_7,\mathbb C^\*)=\mathbb Z/6$, the Schur multiplier. Such a bundle is linearised for the Schur cover $6.A_7$, whose centre $Z\cong\mathbb Z/6$ acts on fibres by a character $\varepsilon$.

Chapter 5 used linearised bundles only, which is right for accessory irrationalities (Theorem 5.2 always produces linearised bundles). For the
geometry of the curve, the twisted classes matter. They are where the low-degree projective models live.

**Main results** on a $(2,4,7)$ curve $C$ (genus 136), with $\tau$ an involution:
1. The degrees of $A_7$-invariant classes form $15\mathbb Z$, against $90\mathbb Z$ for linearised ones. The Mumford class determines the degree mod 90 (Corollary 7.2).
2. Every invariant class of degree $<60$ has no sections beyond the trivial one (Theorem 7.4). The proof uses $\operatorname{gon}\ge25$.
3. There is an invariant class $L_{60}$ of degree 60, with Mumford class of order 3 and $H^0(L_{60})\supseteq\mathbf 6$. Hence a birational $3.A_7$-equivariant model $C\to\mathbb P^5$ of degree 60 (Theorem 7.5).
4. **$\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\langle\tau\rangle)\le21$** (Corollary 7.6). The previous bounds were 56 and 28. So
   $$25\le\operatorname{gon}(C)\le42,\qquad13\le\operatorname{gon}(C/\langle\tau\rangle)\le21 .$$

Status tags as in Chapter 2. Everything is reproduced by `twisted_rr.py` → `twisted_rr_output.txt` (22 s).

## 7.1 The groups

- **$2.A_7\subset SU(4)$** acts on the half-spin space $V_4$, with $\wedge^2V_4=\mathbf 6$. Each $g\in A_7$ is lifted as a product of Clifford reflections $\gamma(e_i-e_j)/\sqrt2$.
- **$3.A_7$** is given by the central lift
  $$\langle x,y,z\mid z^3,\ [x,z],\ [y,z],\ y^5,\ (xy)^7,\ x^3,\ (xyxy^{-1})^2=z,\ (xy^{-2}xy^2)^2\rangle$$
  of $A_7=\langle x,y\mid x^3,y^5,(xy)^7,(xyxy^{-1})^2,(xy^{-2}xy^2)^2\rangle$. Coset enumeration gives order 7560, and the group is perfect.
- **$6.A_7=2.A_7\times_{A_7}3.A_7$** is built with explicit 2-cocycles. It has 40 classes, and its character table is computed by Burnside–Dixon.
  The faithful irreducible degrees are $\{4,4,14,14,20,20,36\}$ (spin), $\{6,15,15,21,21,24,24\}$ (each order-3 character) and $\{6,6,24,24,36\}$ (each order-6 character).

## 7.2 The invariant Picard group [P][X]

Let $x_1,x_2,x_3$ be the branch generators (orders 2, 4, 7, with $x_1x_2x_3=1$), and fix lifts $\hat x_i\in6.A_7$.
- A *local datum* of a twisted class $L$ is the eigenvalue $\lambda_i=e^{2\pi ir_i/e_i}$ of $\hat x_i$ on the fibre over its fixed point.
- Write $\varepsilon(\hat x_1\hat x_2\hat x_3)=e^{2\pi i\rho_0}$.

**Proposition 7.1.**
1. For every character $\varepsilon$ of $Z$ and every local datum with $\lambda_i^{e_i}=\varepsilon(\hat x_i^{e_i})$ there is an $\varepsilon$-twisted invariant line bundle.
2. Its degree is
   $$\deg L=2520\Big(n+\sum_i\frac{r_i}{e_i}-\rho_0\Big),\qquad n\in\mathbb Z .$$

*Proof.*
1. **Setup.** Let $\tilde\Delta=\langle c_1,c_2,c_3\mid c_1^2=c_2^4=c_3^7=c_1c_2c_3\ (=h\text{, central})\rangle$ be the universal central extension of $\Delta(2,4,7)$. Then $K^{w}$ is $\tilde\Delta$-equivariant on $\mathbb H$, with $h$ acting by $e^{-2\pi iw}$ and $c_i$ acting at its fixed point by $e^{-2\pi iw/e_i}$.
2. **Splitting.** An equivariant line bundle $L$ on $\mathbb H$ splits as $K^w\otimes F_\chi$, with $F_\chi$ flat. Here $\chi$ is a character of $\tilde\Delta\times_\Delta\hat\Delta$ with $\chi(h)=e^{2\pi iw}$ and $\chi|_Z=\varepsilon$.
3. **Solving.** The relations $\alpha_i^{e_i}=e^{2\pi iw}\varepsilon(\hat x_i^{e_i})$ and $\prod\alpha_i=e^{2\pi iw}\varepsilon(\hat x_1\hat x_2\hat x_3)$ are solvable exactly when
   $\frac3{28}w\equiv\sum r_i/e_i-\rho_0\pmod1$, where $\frac3{28}=-\chi_{\rm orb}$. So $L$ exists and descends to $C$.
4. **Degree.** By Chern–Weil, $\deg L=270\,w$. $\square$

*Check.* Theta characteristics have $r_i=-\tfrac12$ and $\rho_0=\tfrac12$ (canonical $SL_2$ lifts), which gives degree 135.

**Corollary 7.2.** $\deg\big(\mathrm{Pic}(C)^{A_7}\big)=15\mathbb Z$. By Mumford class:

| Mumford class | trivial | order 2 (spin) | order 3 | order 6 |
|---|---|---|---|---|
| degree mod 90 | 0 | 45 | $\pm30$ | $\pm15$ |

The order-3 twists are not of the form $K^{1/3}$: the Euler class of $\Delta(2,4,7)$ is divisible by 3, because $\tilde\Delta/\langle h^3\rangle\to\Delta$ splits ($c_1h,\ c_2h^2,\ c_3h^2$).

## 7.3 Twisted holomorphic Lefschetz [P][X]

For $\hat g\in6.A_7$ over $g\ne1$:
$$\mathrm{tr}(\hat g\mid H^0)-\mathrm{tr}(\hat g\mid H^1)=\sum_{p\in\mathrm{Fix}(g)}\frac{\lambda(\hat g,p)}{1-a_p^{-1}}.$$
Central elements act by $\varepsilon(z)(\deg L-135)$. The local eigenvalues $\lambda(\hat g,p)$ are transported from the local data by conjugation.

**Validation.** With $\varepsilon=1$ this reproduces `equivariant_rr.py`:
- $\chi(B+T)=-6+10-14_a-14_b-21$;
- $\chi(B)=-\overline{10}-35$;
- $H^0(K)=10+2\cdot\overline{10}+15+21+2\cdot35$.

Both orientation conventions give complex-conjugate tables. Every conclusion below holds in both.

**Forced sections in low degree** (the positive part of $\chi$, so $H^0\supseteq$ it). These are the same on all four classes.

| degree | twist | positive part of $\chi(L)$ |
|---|---|---|
| 15, 105 | order 6 | none |
| 30 | order 3 | none |
| 45, 135 | spin | none |
| **60** | **order 3** | **$\mathbf 6$** (one of the two local data) |
| 75 | order 6 | none |
| 90 | trivial | $10$ (for $B+T$) |
| 120 | order 3 | $6+15$ |
| 165 | order 6 | $6+24$ |

## 7.4 Vanishing below 60

**Lemma 7.3 (local-character base locus) [P].** Let $L$ be twisted with local data $\lambda_i$, and let $f$ be an invariant of degree $k$ of the linear group acting on $H^0(L)$.
Then $f|_C$ is an invariant section of $L^k$. Its order at a point over branch $i$ satisfies $\lambda_i^k=a_i^{o_i}$, so $o_i$ is fixed mod $e_i$.
If the least invariant divisor with these residues, $\sum_i o_i\frac{2520}{e_i}$ (with $0\le o_i<e_i$), exceeds $k\deg L$, then $f|_C\equiv0$.

This sharpens Lemma 5.9, which only uses orbit sizes.

**Theorem 7.4 (no invariant sections below 60) [P][X].** Every $A_7$-invariant class on $C$ of degree $<60$, other than $\mathcal O$, has $h^0=0$.

*Proof.* By Chapters 2–3, $\operatorname{gon}(C)\ge25$, so $\mathrm{Cliff}(C)\ge22$ (Coppens–Martens: $\operatorname{gon}\le\mathrm{Cliff}+3$). $H^0$ of a twisted class is a sum of faithful-type irreducibles.
- **Degree 0.** The classes are $\mathcal O$ and $T$, and $h^0(T)=0$.
- **Degree 15.** $h^0\le1$, below the gonality, and faithful order-6 irreducibles have dimension $\ge6$.
- **Degree 30.** Clifford gives $h^0\le5$, and order-3 irreducibles have dimension $\ge6$.
- **Degree 45 (spin).** Clifford gives $h^0\le12$, so a nonzero $H^0$ contains $V_4$ or $V_4^\*$. This gives $\psi:C\to\mathbb P^3$, birational because $\pi(15,3)=42<136$, with image of degree 45.
  1. *Invariants.* The Molien series of $2.A_7$ on $V_4$ starts $1+t^8+t^{12}+t^{14}+t^{16}+t^{18}+2t^{20}+\dots$
  2. *Forced vanishing.* At degree 45, Lemma 7.3 forces $f_{14}|_C=f_{18}|_C=0$. The least admissible divisors have degree 3150 against 630, and 3330 against 810.
  3. *The lines.* Each involution lifts with eigenvalues $(i,i,-i,-i)$. Since $i^{14}=i^{18}=-1$, both $f_{14}$ and $f_{18}$ vanish on its two fixed lines: 210 distinct lines in all.
  4. *Coprimality.* $f_{14}$ and $f_{18}$ share no factor: a common factor would be an invariant of degree $\le14$, or the product of an $A_6$-orbit of quadrics, and neither divides both. Numerically, their roots on a random line are separated.
  5. *Bézout.* $Z(f_{14})\cap Z(f_{18})$ is a curve of degree 252 containing the 210 lines. So $\psi(C)$, of degree $45>42$, cannot fit. $\square$

In particular: **there is no spinor model** (no $2.A_7$-equivariant map $C\to\mathbb P(V_4)$ of degree 45). The degree-90 series $B+T$ is therefore not a Veronese square.

## 7.5 The degree-60 model in $\mathbb P^5$

**Theorem 7.5 [P][X].** Every $(2,4,7)$ curve carries an invariant class $L_{60}$ of degree 60, with Mumford class of order 3, such that
$$\chi(L_{60})=\mathbf 6-15-2\cdot21-24,\qquad\text{hence }H^0(L_{60})\supseteq\mathbf 6 .$$
The resulting map $\varphi:C\to\mathbb P^5$ has the following properties.
1. **Morphism.** It is base-point-free: the base locus is invariant of degree $\le60<360$.
2. **Birational onto a degree-60 curve.** Otherwise the image would be a faithful $A_7$-curve of degree $\le30$ in $\mathbb P^5$, of genus $\le\pi(30,5)=91<136$.
3. **Equivariant** for $3.A_7$ acting on $\mathbb P^5$ through its exceptional 6-dimensional representation.
4. **The involutions.** For each involution $\tau$, its lift has eigenvalues $(+1)^4(-1)^2$ on $\mathbf 6$, and *all 18 fixed points of $\tau$ have eigenvalue $+1$*. So they lie in the $\mathbb P^3=\mathbb P(E_+)$.
5. **Plücker check.** The vanishing sequences are $(0,1,2,3,4,6)$ at points over the 2- and 4-branches and $(0,1,2,3,4,5)$ at the 7-points. The Plücker total is $6(60+5\cdot135)=4410=1890+2520$, which is consistent.

So the least degree of an invariant class with $h^0\ge2$ on these curves is exactly **60**: Theorems 7.4 and 7.5. Compare $\mu=90$ for linearised classes (Chapter 5).
The value coincides with $a(A_7)=60$. The $\mathbf 6$ restricts irreducibly to both classes of $L_2(7)$, so the model does *not* project $L_2(7)$-equivariantly onto Klein's $\mathbb P^2$; the coincidence stays unexplained.

**Corollary 7.6 (new upper bounds) [P][X].** $\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\langle\tau\rangle)\le21$.

*Proof.*
1. **The pencil.** Let $E_-\subset\mathbf 6\subseteq H^0(L_{60})$ be the 2-dimensional $(-1)$-eigenspace of $\hat\tau$.
2. **Base points.** At a fixed point $p$ with eigenvalue $+1$, an anti-invariant section satisfies $s(p)=-s(p)$, so $s(p)=0$. The pencil $|E_-|$ therefore has the 18 fixed points of $\tau$ as base points, and its moving part has degree $\le60-18=42$.
3. **Descent.** The ratio of two anti-invariant sections is $\tau$-invariant, so the pencil factors through $C/\langle\tau\rangle$ with degree $\le21$. $\square$

**Remark (the $\mu=90$ series) [P].** On $B+T$ the symmetric-matrix family $Q_p\in\mathrm{Sym}^2V_4$ (§5.4) has generic rank 3 or 4:
- **rank 1** would give the excluded spinor model;
- **rank 2** would put a $\mathbf 6$ in $H^0(B)$ or $H^0(B+T)$.

Lemma 7.3 kills $p_2,p_3,p_5,p_6$ (and $p_7$ for $B+T$). Then the image would lie in the root locus of $t^7+ut^3+v$, an irreducible $S_7$-curve of genus 691, or (for $B+T$) in a finite set. Both are impossible.

**A search for better pencils** (`twisted_rr.py` part (8)). For every class with forced sections and degree $\le270$, the scan combines:
- $\tau$-eigenspaces;
- Klein four-group isotypic components;
- fixed-point base loci.

Nothing beats 42.

## 7.6 Open

- **The exact gonality**, in $[25,42]$. Does the pencil $|E_-|$ have base points beyond the 18? These would be $\tau$-symmetric nodes of $\varphi(C)$ on $\mathbb P(E_+)$.
- **Equations.** Equations of the degree-60 curve in $\mathbb P^5$, and the invariant ring of $3.A_7$ on $\mathbf 6$ (degrees $\equiv0\bmod3$).
- **The kernel map.** Is the rank-3 kernel map of the $\mu=90$ series (a spin class of degree 135, $M^2=(B+T)^3$) a genuine $\mathbb P^3$-model?
- **Other signatures.** The same theory on the other faithful $A_7$-curves, and $\tilde\mu(A_7)$ over all of them.
