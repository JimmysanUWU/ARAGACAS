# Chapter 7. Schur-twisted geometry: a model in $\mathbb P^5$ and $\operatorname{gon}(C)\le42$

**The idea.** An $A_7$-*invariant* line bundle need not be $A_7$-*linearised*.
- The obstruction, its *Mumford class*, lies in the Schur multiplier $H^2(A_7,\mathbb C^\*)=\mathbb Z/6$.
- Such a bundle is linearised for the Schur cover $6.A_7$, whose centre $Z\cong\mathbb Z/6$ acts on fibres by a character $\varepsilon$.

Accessory irrationalities only ever produce linearised bundles (Theorem 5.2). The geometry of the curve, however, lives in the twisted classes.

**Main results** for a $(2,4,7)$ curve $C$ (genus 136) and an involution $\tau$:
1. Invariant classes have degrees $15\mathbb Z$; linearised ones have degrees $90\mathbb Z$ (Corollary 7.2).
2. No invariant class of degree $<60$ other than $\mathcal O$ has sections. The proof uses $\operatorname{gon}\ge25$ (Theorem 7.4).
3. A degree-60 class with Mumford class of order 3 gives a birational $3.A_7$-equivariant model $C\to\mathbb P^5$, lying on an explicit invariant cubic fourfold (Theorem 7.5, §7.6).
4. **$\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\langle\tau\rangle)\le21$**, against the previous 56 and 28 (Corollary 7.6). Hence
   $$25\le\operatorname{gon}(C)\le42,\qquad13\le\operatorname{gon}(C/\langle\tau\rangle)\le21,\qquad25\le\gamma(A_7)\le42 .$$

Status tags as in `README.md`. Reproduced by `twisted_rr.py` → `twisted_rr_output.txt` (25 s) and `twisted_survey.py` → `twisted_survey_output.txt` (90 s).

## 7.1 The groups [X]

| group | construction | faithful irreducible degrees |
|---|---|---|
| $2.A_7\subset SU(4)$ | half-spin $V_4$ of the standard 6, so $\wedge^2V_4=\mathbf 6$; lifts are products of Clifford reflections | $4,4,14,14,20,20,36$ |
| $3.A_7$ | the central lift below, by coset enumeration (order 7560, perfect) | $6,15,15,21,21,24,24$ per order-3 character |
| $6.A_7=2.A_7\times_{A_7}3.A_7$ | explicit 2-cocycles; 40 classes; Burnside–Dixon character table | $6,6,24,24,36$ per order-6 character |

The central lift is
$$3.A_7=\langle x,y,z\mid z^3,\ [x,z],\ [y,z],\ x^3,\ y^5,\ (xy)^7,\ (xyxy^{-1})^2=z,\ (xy^{-2}xy^2)^2\rangle .$$
Setting $z=1$ gives $A_7$.

## 7.2 The invariant Picard group [P][X]

**Notation.**
- $x_1,x_2,x_3$ are the branch generators, of orders $e_i=2,4,7$, with $x_1x_2x_3=1$; $\hat x_i\in6.A_7$ are fixed lifts.
- The *local datum* of a twisted class $L$ is the eigenvalue $\lambda_i=e^{2\pi ir_i/e_i}$ of $\hat x_i$ on the fibre over its fixed point.
- $\rho_0$ is defined by $\varepsilon(\hat x_1\hat x_2\hat x_3)=e^{2\pi i\rho_0}$.

**Proposition 7.1.** For every character $\varepsilon$ of $Z$ and every local datum with $\lambda_i^{e_i}=\varepsilon(\hat x_i^{e_i})$, an $\varepsilon$-twisted invariant line bundle exists. Its degree is
$$\deg L=2520\Big(n+\sum_i\frac{r_i}{e_i}-\rho_0\Big),\qquad n\in\mathbb Z .$$

*Proof.*
1. **The universal extension.** Let $\tilde\Delta=\langle c_1,c_2,c_3\mid c_1^2=c_2^4=c_3^7=c_1c_2c_3=:h\rangle$ be the universal central extension of $\Delta(2,4,7)$.
   Then $K^w$ is $\tilde\Delta$-equivariant on $\mathbb H$: $h$ acts by $e^{-2\pi iw}$, and $c_i$ acts at its fixed point by $e^{-2\pi iw/e_i}$.
2. **Splitting.** Every equivariant line bundle on $\mathbb H$ is $K^w\otimes F_\chi$, with $F_\chi$ flat. Here $\chi$ is a character of $\tilde\Delta\times_\Delta\hat\Delta$ with $\chi(h)=e^{2\pi iw}$ and $\chi|_Z=\varepsilon$.
3. **Solving.** The relations $\alpha_i^{e_i}=e^{2\pi iw}\varepsilon(\hat x_i^{e_i})$ and $\prod\alpha_i=e^{2\pi iw}\varepsilon(\hat x_1\hat x_2\hat x_3)$ are solvable exactly when
   $\tfrac3{28}w\equiv\sum r_i/e_i-\rho_0\pmod1$. Here $\tfrac3{28}=-\chi_{\rm orb}$. So $L$ exists, and it descends to $C$.
4. **Degree.** By Chern–Weil, $\deg L=270\,w$. $\square$

*Check.* Theta characteristics ($r_i=-\tfrac12$, $\rho_0=\tfrac12$ for canonical $SL_2$ lifts) get degree 135.

**Corollary 7.2.** $\deg\mathrm{Pic}(C)^{A_7}=15\mathbb Z$, and the Mumford class fixes the degree mod 90:

| Mumford class | trivial | order 2 (spin) | order 3 | order 6 |
|---|---|---|---|---|
| degree mod 90 | 0 | 45 | $\pm30$ | $\pm15$ |

Order-3 twists are not cube roots of $K$. The reason is that $\tilde\Delta/\langle h^3\rangle\to\Delta$ splits (take $c_1h,\ c_2h^2,\ c_3h^2$), so the Euler class of $\Delta$ is divisible by 3.

## 7.3 Twisted holomorphic Lefschetz [P][X]

For $\hat g\in6.A_7$ lying over $g\ne1$:
$$\mathrm{tr}(\hat g\mid H^0)-\mathrm{tr}(\hat g\mid H^1)=\sum_{p\in\mathrm{Fix}(g)}\frac{\lambda(\hat g,p)}{1-a_p^{-1}},$$
where $\lambda(\hat g,p)$ is obtained from the local data by conjugation. A central element $z$ acts by $\varepsilon(z)(\deg L-135)$.

**Checks.**
- With $\varepsilon=1$ this reproduces `equivariant_rr.py`: $\chi(B+T)=-6+10-14_a-14_b-21$, $\chi(B)=-\overline{10}-35$, and $H^0(K)=10+2\cdot\overline{10}+15+21+2\cdot35$.
- The two orientation conventions give complex-conjugate tables. Every conclusion below holds in both.

**Forced sections** (positive part of $\chi$), the same on all four classes:

| degree | twist | $H^0\supseteq$ |
|---|---|---|
| 15, 75, 105 | order 6 | — |
| 30 | order 3 | — |
| 45, 135 | spin | — |
| **60** | **order 3** | **$\mathbf 6$** (for one of the two local data) |
| 90 | trivial | $10$ (for $B+T$) |
| 120 | order 3 | $6+15$ |
| 165 | order 6 | $6+24$ |

## 7.4 Vanishing below degree 60

**Lemma 7.3 (local-character base locus) [P].**
- *Setting.* $L$ is a twisted class with local data $\lambda_i$, and $f$ is an invariant of degree $k$ for the group acting on $H^0(L)$.
- *Orders.* $f|_C$ is an invariant section of $L^k$. Its order $o_i$ at a point over branch $i$ satisfies $\lambda_i^k=a_i^{o_i}$, so it is fixed mod $e_i$.
- *Conclusion.* If $\sum_io_i\cdot2520/e_i$ (with $0\le o_i<e_i$) exceeds $k\deg L$, then $f|_C\equiv0$.

This sharpens Lemma 5.9, which uses orbit sizes only.

**Theorem 7.4 [P][X].** Every $A_7$-invariant class of degree $<60$ other than $\mathcal O$ has $h^0=0$.

*Proof.* Chapters 2–3 give $\operatorname{gon}\ge25$, so $\mathrm{Cliff}\ge22$ (Coppens–Martens: $\operatorname{gon}\le\mathrm{Cliff}+3$). $H^0$ of a twisted class is a sum of faithful irreducibles of the twist.
- **Degree 0.** $T$ has no sections.
- **Degree 15.** $h^0\le1$ (below the gonality), while order-6 irreducibles have dimension $\ge6$.
- **Degree 30.** Clifford gives $h^0\le5$, while order-3 irreducibles have dimension $\ge6$.
- **Degree 45 (spin).** Clifford gives $h^0\le12$, so a nonzero $H^0$ contains $V_4$ or $V_4^\*$. This gives $\psi:C\to\mathbb P^3$ of degree 45, birational because $\pi(15,3)=42<136$.
  1. The invariants of $2.A_7$ on $V_4$ have Hilbert series $1+t^8+t^{12}+t^{14}+t^{16}+t^{18}+2t^{20}+\dots$
  2. Lemma 7.3 forces $f_{14}|_C=f_{18}|_C=0$: the least admissible divisors have degree 3150 against 630, and 3330 against 810.
  3. Involutions lift with eigenvalues $(i,i,-i,-i)$, and $i^{14}=i^{18}=-1$. So $f_{14}$ and $f_{18}$ vanish on the two fixed lines of each involution: 210 distinct lines.
  4. $f_{14}$ and $f_{18}$ are coprime. A common factor would be an invariant of degree $\le14$, or a product of an $A_6$-orbit of quadrics, and neither divides both. Numerically, their roots on a random line are separated.
  5. So $Z(f_{14})\cap Z(f_{18})$ is a curve of degree 252 containing 210 lines, which leaves no room for $\psi(C)$ of degree $45>42$. $\square$

In particular, **there is no spinor model**: no $2.A_7$-equivariant map $C\to\mathbb P(V_4)$ of degree 45.

## 7.5 The degree-60 model and $\operatorname{gon}\le42$

**Theorem 7.5 [P][X].** Every $(2,4,7)$ curve carries an invariant class $L_{60}$ of degree 60, with Mumford class of order 3 and
$$\chi(L_{60})=\mathbf 6-15-2\cdot21-24,\qquad\text{so }H^0(L_{60})\supseteq\mathbf 6 .$$
The map $\varphi:C\to\mathbb P^5$ given by $\mathbf 6$ has these properties:
1. **Base-point-free.** The base locus is invariant of degree $\le60<360$.
2. **Birational onto a curve of degree 60.** Otherwise the image would be a faithful $A_7$-curve of degree $\le30$ in $\mathbb P^5$, of genus $\le\pi(30,5)=91<136$.
3. **Equivariant** for $3.A_7$ acting through its exceptional $\mathbf 6$.
4. **Involutions.** For each involution $\tau$, the lift $\hat\tau$ has eigenvalues $(+1)^4(-1)^2$ on $\mathbf 6$, and *all 18 fixed points of $\tau$ have eigenvalue $+1$*. So they lie in $\mathbb P(E_+)\cong\mathbb P^3$.
5. **Plücker check.** The vanishing sequences are $(0,1,2,3,4,6)$ over the 2- and 4-branches and $(0,1,2,3,4,5)$ over the 7-branch. The total weight $6(60+5\cdot135)=4410=1890+2520$ is consistent.

By Theorems 7.4 and 7.5, **60 is the least degree of an invariant class with $h^0\ge2$**; compare $\mu=90$ for linearised classes. The equality with $a(A_7)=60$ remains unexplained: the $\mathbf 6$ restricts irreducibly to both classes of $L_2(7)$, so there is no $L_2(7)$-projection to Klein's $\mathbb P^2$.

**Corollary 7.6 [P][X].** $\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\langle\tau\rangle)\le21$.

*Proof.*
1. Let $E_-\subset\mathbf 6\subseteq H^0(L_{60})$ be the 2-dimensional $(-1)$-eigenspace of $\hat\tau$.
2. At a fixed point with eigenvalue $+1$, an anti-invariant section satisfies $s(p)=-s(p)$, so it vanishes there. The pencil $|E_-|$ therefore has the 18 fixed points as base points, and its moving part has degree $\le42$.
3. The ratio of two anti-invariant sections is $\tau$-invariant, so the pencil descends to $C/\langle\tau\rangle$ with degree $\le21$. $\square$

**Nothing better from these classes** (`twisted_rr.py` part (8)). Every class with forced sections and degree $\le270$ was tested against involution eigenspaces and Klein four-group isotypic components: none beats 42.

## 7.6 The invariant cubic fourfold [P][X]

- **Invariants.** The invariants of $3.A_7$ on $\mathbf 6$ have Hilbert series $1+t^3+3t^6+5t^9+11t^{12}+18t^{15}+\dots$
- **Forced vanishing.** By Lemma 7.3, every invariant of degree 3, 9 or 15 vanishes on $\varphi(C)$. So **$\varphi(C)$ lies on the unique invariant cubic fourfold $X_3$**.
- **Restrictions.** In degrees 6 and 21 the invariants restrict to sections with divisors $D_7$ and $2D_4$. In degree 42, the pencil $\langle S_6^7,T_{21}^2\rangle$ is the quotient map $C\to C/A_7$.
- **The equation** (`twisted_rr.py` part (10)). The $\mathbf 6$ is the projection of the 21-dimensional $\mathrm{Ind}_{S_5\times\mathbb Z/3}^{3.A_7}(1\otimes\omega)$. In eigencoordinates $x_k$ of an order-7 element $\hat c$ ($\hat c\,x_k=\zeta_7^kx_k$), which its normaliser permutes,
  $$X_3:\ x_1x_2x_4+\beta\,x_3x_5x_6+\gamma\,(x_1^2x_5+x_2^2x_3+x_4^2x_6)+\delta\,(x_1x_3^2+x_2x_6^2+x_4x_5^2)=0,$$
  $$\frac{\gamma^3}\beta=\frac{23-7\sqrt{21}}{16},\qquad\frac{\delta^3}{\beta^2}=\frac{23+7\sqrt{21}}{16}\qquad(\text{the roots of }64x^2-184x-125).$$
  The remaining freedom is to rescale $(x_3,x_5,x_6)$ and to multiply by a cube root of unity.
  The 360 points of $\varphi(C)$ over the 7-branch are coordinate points: the eigenlines of the conjugates of $\hat c$.
- **Uniqueness.** By duality $K-3L_{60}=B+T$, so $H^0(3L_{60})\supseteq6+14_a+14_b+21$, which is 55-dimensional. Since $\mathrm{Sym}^3\mathbf 6=1\oplus(55)$, $X_3$ is the only cubic through $\varphi(C)$, provided $h^0(B+T)=10$.

## 7.7 The $\mu=90$ series revisited [P]

The family $Q_p\in\mathrm{Sym}^2V_4$ of §5.4, attached to $B+T$, has generic rank 3 or 4.
- **Rank 1** would be the spinor model excluded by Theorem 7.4.
- **Rank 2** would put a $\mathbf 6$ into $H^0(B)$ or $H^0(B+T)$. Lemma 7.3 kills $p_2,p_3,p_5,p_6$ (and $p_7$ for $B+T$), so the image would lie in the root curve of $t^7+ut^3+v$ (an irreducible $S_7$-curve of genus 691), or in a finite set. Both are impossible.

## 7.8 The other rigid signatures

`twisted_survey.py` repeats §§7.2–7.5 on all 26 rigid faithful $A_7$-curves of genus $\le335$. It tests every eigenspace of every element with fixed points:

| signature | $g$ | best pencil $\le$ | source |
|---|---|---|---|
| $(2,4,7)$ | 136 | **42** | order-3 class of degree 60; an involution |
| $(3,3,7)$ | 241 | 48 | spin class of degree 60 with $V_4\subseteq H^0$, so a spinor model in $\mathbb P^3$ does exist here; a 3-cycle |
| $(2,7,7)$, one curve | 271 | 84 | spin class of degree 90 |
| $(4,4,4)$ | 316 | 95 | order-6 class of degree 105 |
| $(2,5,7)$, $(2,6,7)$, $(2,7,7)$ | 199–271 | 104–108 | order-3 classes of degree 120 |
| $(3,3,5)$, $(3,3,6)$, $(3,4,4)$, $(3,4,5)$, $(3,4,6)$ | 169–316 | 192–299 | |

So the genus-136 curves realise the bound $\gamma(A_7)\le42$.

## 7.9 Open

- **The exact gonality**, in $[25,42]$. Does $|E_-|$ have base points beyond the 18? Any extra ones are $\tau$-symmetric nodes of $\varphi(C)$ on $\mathbb P(E_+)$.
- **Equations of $\varphi(C)$ beyond $X_3$.** For instance, the sextic and nonic invariants that vanish on it.
- **The kernel map** of the $\mu=90$ series: is it a genuine $\mathbb P^3$-model (a spin class of degree 135 with $M^2=3(B+T)$)?
- **Other curves.** Could another curve beat 42 by a mechanism other than fixed-point base loci?
