# Chapter 7. Schur-twisted geometry: the degree-60 model and $\operatorname{gon}(C)\le42$

**Idea.** An $A_7$-*invariant* line bundle need not be *linearised*.
- **The obstruction.** It is the bundle's Mumford class, in $H^2(A_7,\mathbb C^\*)=\mathbb Z/6$. The bundle is linearised for the Schur cover $6.A_7$, whose centre $Z$ acts on fibres by a character $\varepsilon$.
- **Degrees and sections.** Twisted classes exist in every degree in $15\mathbb Z$, against $90\mathbb Z$ for linearised ones. Holomorphic Lefschetz, twisted by $\varepsilon$, computes their sections.
- **The first twisted class with sections** has degree 60 and Mumford class of order 3. Its six sections embed $C$ in $\mathbb P^5$, on a cubic fourfold.
- **The pencil of degree 42.** An involution's $(-1)$-eigenspace in those six sections is two-dimensional, and its sections vanish at the 18 fixed points. What remains is a pencil of degree 42.
- **Lower bounds feed back.** The lower bound 25 of Chapters 2–3 is what keeps the twisted classes below degree 60 empty (Clifford).

Accessories over bases with a fixed point produce linearised bundles (Theorem 5.2). Over other bases, the twisted classes that can occur are exactly those allowed by the Amitsur subgroup of the base (Theorem 5.15); Example 5.12 realises $L_{60}$ over a rational base. The geometry of $C$ lives in the twisted classes.

**Main results** ($C$ a $(2,4,7)$ curve, $\tau$ an involution).
1. Invariant classes have degrees $15\mathbb Z$; linearised ones $90\mathbb Z$ (Cor. 7.2). None of degree $<60$ has sections (Thm 7.4).
2. A degree-60 class of order 3 gives a $3.A_7$-equivariant embedding $\varphi:C\hookrightarrow\mathbb P^5$:
   - its image lies on the Laza–Zheng $A_7$-cubic fourfold;
   - it is cut out by that cubic and 15 quartics (Thm 7.5, Prop. 7.8).
3. **$\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\tau)\le21$**, from an exact base-point-free pencil (Cor. 7.6, Prop. 7.7). Hence
   $$25\le\operatorname{gon}(C)\le42,\qquad13\le\operatorname{gon}(C/\tau)\le21,\qquad25\le\gamma(A_7)\le42 .$$
4. Each Klein subgroup $L_2(7)$ cuts $\varphi(C)$ by a quadric in exactly $5\times$ its 24 points of order 7 (Prop. 7.9).
5. Elliptic subcovers have degree $\ge60$ (Prop. 7.10).

Scripts: `twisted_rr.py` (25 s), `twisted_survey.py` (90 s), `tau_pencil.py` (1 min), `p5_curve.py` (7 min), `elliptic_subcovers.py` (80 s); outputs `*_output.txt`.

## 7.1 The groups [X]

| group | construction | faithful irreducible degrees |
|---|---|---|
| $2.A_7\subset SU(4)$ | half-spin $V_4$ of the standard 6 ($\wedge^2V_4=\mathbf 6$); Clifford lifts | $4,4,14,14,20,20,36$ |
| $3.A_7$ | coset enumeration of the presentation below (order 7560, perfect) | $6,15,15,21,21,24,24$ per order-3 character |
| $6.A_7=2.A_7\times_{A_7}3.A_7$ | explicit 2-cocycles; 40 classes; Burnside–Dixon | $6,6,24,24,36$ per order-6 character |

$$3.A_7=\langle x,y,z\mid z^3,\ [x,z],\ [y,z],\ x^3,\ y^5,\ (xy)^7,\ (xyxy^{-1})^2=z,\ (xy^{-2}xy^2)^2\rangle\qquad(z=1\text{ gives }A_7).$$

## 7.2 The invariant Picard group [P][X]

Let $x_1x_2x_3=1$ be the branch generators (orders $e_i=2,4,7$), and $\hat x_i\in6.A_7$ fixed lifts. The *local datum* of a twisted class is the eigenvalue $\lambda_i=e^{2\pi ir_i/e_i}$ of $\hat x_i$ on the fibre at its fixed point. Define $\rho_0$ by $\varepsilon(\hat x_1\hat x_2\hat x_3)=e^{2\pi i\rho_0}$.

**Proposition 7.1.** For every $\varepsilon$ and every local datum with $\lambda_i^{e_i}=\varepsilon(\hat x_i^{e_i})$ there is an $\varepsilon$-twisted invariant line bundle, of degree
$$\deg L=2520\Big(n+\sum_i\frac{r_i}{e_i}-\rho_0\Big),\qquad n\in\mathbb Z .$$

*Proof.*
1. $\tilde\Delta=\langle c_i\mid c_1^2=c_2^4=c_3^7=c_1c_2c_3=:h\rangle$ is the centrally extended triangle group (Milnor): the preimage of $\Delta(2,4,7)$ in the universal cover of $\mathrm{PSL}_2(\mathbb R)$, a central extension by $\langle h\rangle\cong\mathbb Z$. It is not a universal central extension in the group-theoretic sense, since $\Delta(2,4,7)$ has abelianisation $C_2$. It acts on $K^w$ over $\mathbb H$: $h$ acts by $e^{-2\pi iw}$, and $c_i$ acts at its fixed point by $e^{-2\pi iw/e_i}$.
2. Every equivariant line bundle on $\mathbb H$ is $K^w\otimes F_\chi$, with $\chi$ a character of $\tilde\Delta\times_\Delta\hat\Delta$ satisfying $\chi(h)=e^{2\pi iw}$ and $\chi|_Z=\varepsilon$.
3. **The local data.** At the fixed point of $c_i$ the generator must act by the local datum, so $\chi(c_i)=\lambda_ie^{2\pi iw/e_i}$.
   - The relation $c_i^{e_i}=h$ becomes $\lambda_i^{e_i}=\varepsilon(\hat x_i^{e_i})$, the stated condition.
   - The relation $c_1c_2c_3=h$ becomes $\prod_i\lambda_i\,e^{2\pi iw\sum1/e_i}=e^{2\pi iw}e^{2\pi i\rho_0}$. Since $1-\sum_i1/e_i=\tfrac3{28}=-\chi_{\rm orb}$, this says $\tfrac3{28}w\equiv\sum_ir_i/e_i-\rho_0\pmod1$.
4. **The degree.** $\deg L=270w$ (Chern–Weil), and $w\in\tfrac{28}3\big(\sum r_i/e_i-\rho_0+\mathbb Z\big)$ gives the formula. $\square$

*Check:* theta characteristics get degree 135.

**Corollary 7.2.** $\deg\mathrm{Pic}(C)^{A_7}=15\mathbb Z$, and the degree mod 90 is fixed by the Mumford class:

| Mumford class | trivial | order 2 | order 3 | order 6 |
|---|---|---|---|---|
| degree mod 90 | 0 | 45 | $\pm30$ | $\pm15$ |

Order-3 twists are not cube roots of $K$: $\tilde\Delta/\langle h^3\rangle\to\Delta$ splits ($c_1h,c_2h^2,c_3h^2$), so the Euler class of $\Delta$ is divisible by 3.

## 7.3 Twisted holomorphic Lefschetz [P][X]

For $\hat g\in6.A_7$ over $g\ne1$, with $a_p$ the rotation of $g$ at $p$:
$$\mathrm{tr}(\hat g\mid H^0-H^1)=\sum_{p\in\mathrm{Fix}(g)}\frac{\lambda(\hat g,p)}{1-a_p^{-1}},\qquad z\in Z\text{ acts by }\varepsilon(z)(\deg L-135).$$

**Checks.**
- With $\varepsilon=1$ it reproduces Proposition 5.8 and Chevalley–Weil.
- The two orientation conventions give complex-conjugate tables.

**Forced sections** (positive part of $\chi$; the same on all four classes):

| degree | 15, 75, 105 | 30 | 45, 135 | **60** | 90 | 120 | 165 |
|---|---|---|---|---|---|---|---|
| twist | order 6 | order 3 | spin | **order 3** | trivial | order 3 | order 6 |
| $H^0\supseteq$ | — | — | — | **$\mathbf 6$** | $10$ ($B+T$) | $6+15$ | $6+24$ |

## 7.4 Nothing below degree 60

**Lemma 7.3 (local-character base locus) [P].** Let $f$ be an invariant of degree $k$ for the group acting on $H^0(L)$. Then $f|_C\in H^0(L^k)^G$, and its order $o_i$ at a branch-$i$ point satisfies $\lambda_i^k=a_i^{o_i}$. If $\sum_io_i\cdot2520/e_i$ (with $0\le o_i<e_i$) exceeds $k\deg L$, then $f|_C\equiv0$. This sharpens Lemma 5.9.

**Theorem 7.4 [P][X].** An invariant class of degree $<60$ other than $\mathcal O$ has $h^0=0$.

*Proof.* $\operatorname{gon}\ge25$ (Chapters 2–3) gives $\mathrm{Cliff}\ge22$, since $\operatorname{gon}\le\mathrm{Cliff}+3$ (Coppens–Martens). $H^0$ is a sum of faithful irreducibles of the twist.
- **Degree 0 ($T$):** no sections.
- **Degree 15:** $h^0\le1<6$.
- **Degree 30:** Clifford gives $h^0\le5<6$.
- **Degree 45 (spin):** $h^0\le12$, so $H^0\supseteq V_4$ or $V_4^\*$. This gives $\psi:C\to\mathbb P^3$ of degree 45, birational since $\pi(15,3)=42<136$.
  1. The invariants of $2.A_7$ on $V_4$ have Hilbert series $1+t^8+t^{12}+t^{14}+t^{16}+t^{18}+2t^{20}+\dots$
  2. Lemma 7.3 kills $f_{14}$ and $f_{18}$ on $C$: their least admissible divisors have degrees 3150 and 3330, against $14\cdot45=630$ and $18\cdot45=810$.
  3. Involutions lift with eigenvalues $(i,i,-i,-i)$, and $i^{14}=i^{18}=-1$, so $f_{14},f_{18}$ vanish on the 210 lines $\mathbb P(E_{\pm i}(\hat\tau))$. These are distinct: in a unitary model $E_{-i}=E_i^\perp$, so either line determines $\hat\tau$ up to sign.
  4. $f_{14},f_{18}$ are coprime. Their gcd is invariant, because $2.A_7$ is perfect and so has no characters. A nonconstant gcd would have degree 8, 12 or 14, and dividing $f_{14}$ or $f_{18}$ by it would give an invariant of degree 6, 2 or 4; there are none. So $Z(f_{14},f_{18})$ is a curve of degree 252 containing 210 lines, leaving degree 42 for $\psi(C)$, which has degree 45. $\square$ (Steps 3–4 in this form are due to GPT, `gpt/A7_Amitsur_Correspondence.pdf` §9.)

So there is **no spinor model** $C\to\mathbb P(V_4)$.

## 7.5 The degree-60 model and the $\tau$-pencil

**Theorem 7.5 [P][X].** Every $(2,4,7)$ curve carries an invariant class $L_{60}$ of degree 60 and Mumford class of order 3 with
$$\chi(L_{60})=\mathbf 6-15-2\cdot21-24,\qquad\text{so }H^0(L_{60})\supseteq\mathbf 6 .$$
The map $\varphi:C\to\mathbb P^5$ given by $\mathbf 6$ has these properties:
1. **Base-point-free.** An invariant base locus would have degree $\le60<360$.
2. **Birational onto a curve of degree 60.** Otherwise the image would have degree $\le30$ and genus $\le\pi(30,5)=91<136$.
3. **Equivariant** for $3.A_7$ on its exceptional $\mathbf 6$.
4. **Involutions.** The lift $\hat\tau$ has eigenvalues $(+1)^4(-1)^2$, and all 18 fixed points of $\tau$ have eigenvalue $+1$, so they lie in $\mathbb P(E_+)\cong\mathbb P^3$.
5. **Plücker [P][N].** Local characters force four even and two odd orders at the 2- and 4-points, with minimal sequence $(0,1,2,3,4,6)$ of weight 1, and six distinct orders mod 7 at the 7-points, with minimal sequence $(0,\dots,5)$. The minimal sequences carry weight 1890 out of the total $6(60+5\cdot135)=4410$. If they hold, as Proposition 7.8 indicates numerically, the remaining 2520 is one free orbit of simple flexes.

So **60 is the least degree of an invariant class with $h^0\ge2$**, against $\mu=90$ for linearised classes.

**Corollary 7.6 [P][X].** $\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\tau)\le21$.

*Proof.* Sections in the 2-dimensional $(-1)$-eigenspace $E_-$ satisfy $s(p)=-s(p)$ at the 18 fixed points. So $|E_-|$ has them as base points, and its moving part has degree $\le42$. Ratios of anti-invariant sections are $\tau$-invariant, so the pencil descends to $C/\tau$. $\square$

**Proposition 7.7 (exactly 42) [P][N].** The base locus of $|E_-|$ is exactly the 18 fixed points, each simple. So the bound comes from a base-point-free $g^1_{42}$, pulled back from a $g^1_{21}$ on $C/\tau$.

*Proof.*
1. **Simple: [P] at 12 points, [N] at 6.** $\hat\tau$ acts on the order-$k$ section at a fixed point by $(-1)^k$ times a constant, so $E_-$ has two odd orders there. Parity alone does not give a simple zero; Plücker (Theorem 7.5(5)) does at the 2-points.
   - **The 12 two-points.** Without an order-1 section the sequence would be at least $(0,2,3,4,5,6)$, of weight 5. That costs $4\cdot1260>2520$.
   - **The 6 four-points.** The same alternative costs exactly $4\cdot630=2520$, so Plücker allows it. The immersion of Proposition 7.8 [N] excludes it. If it held, $|E_-|$ would have multiplicity 3 there, and $\operatorname{gon}(C/\tau)\le15$.
2. **Other base points are double points.** If $x\notin\mathrm{Fix}\,\tau$ and $\varphi(x)\in\mathbb P(E_+)$, then $\varphi(\tau x)=\varphi(x)$.
3. **At most one orbit.** The extra base locus $B'$ is $C(\tau)$-invariant, with $60-18-\deg B'\ge2\operatorname{gon}(C/\tau)\ge26$.
   - Elements of order 4 in $C(\tau)=(D_4\times S_3)\cap A_7$ square to $\tau$.
   - So $C(\tau)$-orbits off $\mathrm{Fix}\,\tau$ have size 12 or 24.
   - Hence $B'$ is empty or one orbit of 12 simple points, each fixed by an involution $\sigma\in C(\tau)\smallsetminus\{\tau\}$.
4. **A plane.** The commuting order-2 lifts of $\sigma,\tau,\sigma\tau$ all have trace 2. So the joint eigenspaces have dimensions $(3,1,1,1)$, and $B'$ maps into the plane $Q_\sigma=\mathbb P(E_{++})$.
5. **Computation** (`tau_pencil.py`, classes 0 and 12). On each of the 8 planes $Q_\sigma$:
   - the cubic $F_3$ and the sextic $G_6$ (§7.6) meet in all 18 of their Bézout points;
   - the other invariants vanishing on $\varphi(C)$ are nonzero at each of them.

   So $B'=\varnothing$. $\square$

Proposition 7.8 confirms this twice: $\varphi$ is an embedding, and a hyperplane of the pencil meets $\varphi(C)$ in $18+42$ points.

**Nothing better from these classes** (`twisted_rr.py` (8)). No class of degree $\le270$ with forced sections beats 42 against involution eigenspaces or Klein four-group isotypic components.

## 7.6 The curve in $\mathbb P^5$

**The cubic [P][X].**
- **The invariants** of $3.A_7$ on $\mathbf 6$ have Hilbert series $1+t^3+3t^6+5t^9+11t^{12}+18t^{15}+33t^{18}+53t^{21}+86t^{24}+\dots$
- **Forced vanishing.** By Lemma 7.3, all invariants of degrees 3, 9 and 15 vanish on $\varphi(C)$. In degrees 6, 12, 18, 21 and 24, restrictions are multiples of the sections with divisors $D_7,2D_7,3D_7,2D_4,4D_7$, which vanish to exact order 1, 2, 3, 2, 4. So **$\varphi(C)$ lies on the unique invariant cubic $X_3$**.
- **The sextic $G_6$.** The sextics whose restriction to the tangent line at a 7-point has a double zero vanish on $\varphi(C)$; they span $\langle F_3^2,G_6\rangle$.
- **The equation** (`twisted_rr.py` (10)). The $\mathbf 6$ is the projection of $\mathrm{Ind}_{S_5\times\mathbb Z/3}^{3.A_7}(1\otimes\omega)$. In eigencoordinates of an order-7 element ($\hat cx_k=\zeta_7^kx_k$),
  $$X_3:\ x_1x_2x_4+\beta\,x_3x_5x_6+\gamma\,(x_1^2x_5+x_2^2x_3+x_4^2x_6)+\delta\,(x_1x_3^2+x_2x_6^2+x_4x_5^2)=0,\qquad\frac{\gamma^3}\beta,\ \frac{\delta^3}{\beta^2}=\frac{23\mp7\sqrt{21}}{16}.$$
  The 360 seven-points of $\varphi(C)$ are coordinate points of these frames.
- **Identification.** $X_3$ is smooth; Koike proves it by reduction mod 2, and numerically $\min|\nabla F|=0.139$ on the unit sphere. It is the second of the two smooth cubic fourfolds with symplectic $A_7$-action (Laza–Zheng), found by Yang–Yu–Zhu via $3.A_7$. Koike gives it over $\mathbb Q$ as
  $$2\textstyle\sum x_i^3+3\sum_{(ij)\in\{12,34,56\}}(x_i^2x_j+x_ix_j^2)+2(x_2x_3x_5+x_1x_4x_5+x_1x_3x_6+x_2x_4x_6)+4(x_1x_3x_5+x_2x_4x_5+x_2x_3x_6+x_1x_4x_6),$$
  with centre $(abab^{-1})^2$, as in our presentation. So **the genus-136 $A_7$-curves live on the Laza–Zheng–Yang–Yu–Zhu cubic fourfold**.
- **Uniqueness [X][N].** $\mathrm{Sym}^3\mathbf 6=1\oplus6\oplus14_a\oplus14_b\oplus21$, multiplicity-free (`twisted_rr.py` (10)). So $X_3$ is the only cubic iff none of the last four constituents vanishes on $\varphi(C)$, which Proposition 7.8 confirms numerically. Since $K-3L_{60}=B+T$, $h^0(3L_{60})=45+h^0(B+T)\ge55$, with equality (3-normality) iff $h^0(B+T)=10$.
- **Lines.** Odd-degree invariants vanish on the $(-1)$-eigenspaces of lifts of order 2 and 4. So $X_3$ contains:
  - the 105 lines $\mathbb P(E_-(\tau))$, the targets of the pencils of Proposition 7.7;
  - the 315 lines $\mathbb P(E_{-1}(\hat h))$, $h$ of order 4, which carry $\varphi(\mathrm{Fix}\,h)$.

**Proposition 7.8 (equations) [N]** (`p5_curve.py`, class 0).
1. **The 4-points.** On each line $\mathbb P(E_{-1}(\hat h))$, $G_6$ has exactly two double zeros, which are singular points of $\{G_6=0\}$. They are $\varphi(\mathrm{Fix}\,h)$.
2. **More equations.** The invariants of degrees 12, 18, 24 vanishing there, and those of degree 21 vanishing at a 7-point, vanish on $\varphi(C)$, by the exact orders above. With $F_3$, $G_6$ and the forced invariants, they cut out a smooth branch through each 7-point. Continuation and the group spread it over $\varphi(C)$.
3. **Hilbert function.** In degrees $2,3,4,5$ the image of $\mathrm{Sym}^k\mathbf 6$ has dimension $21,55,105,165$. So:
   - $\varphi(C)$ lies on no quadric and on one cubic;
   - it lies on 21 quartics, namely $F_3\cdot\mathbf 6$ and 15 more;
   - it is 4- and 5-normal, since $h^0(4L_{60})=105$ and $h^0(5L_{60})=165$.
4. **Embedding.** $5L_{60}$ is very ample (degree $300\ge2g+1$), and all its sections are quintic forms. So $\varphi$ is an embedding.
5. **Cut out.** $X_3$ and the 15 quartics meet a random hyperplane in exactly 60 points, so they cut out $\varphi(C)$. A hyperplane of the $\tau$-pencil meets them in $18+42$ points.

So **$\varphi(C)$ is a smooth curve of degree 60 and genus 136 cut out by the Laza–Zheng cubic and 15 quartics.**

**Proposition 7.9 (Klein quadrics) [P][N].** Let $H\cong L_2(7)$ be one of the 30 Klein subgroups of $A_7$, in two classes of 15.
- $\mathbf 6|_H\cong\mathrm{Sym}^2(3)$, so the $H$-invariant quadrics form a line $\langle Q_H\rangle$: the lift of Klein's quartic through $\mathrm{Sym}^4\subset\mathrm{Sym}^2\mathrm{Sym}^2$.
- Then $\operatorname{div}(Q_H|_C)=5\,O_H$, where $O_H$ is the set of 24 seven-points whose stabiliser lies in $H$.

*Proof.* $Q_H|_C\ne0$, because $\varphi(C)$ lies on no quadric. Its divisor is $H$-invariant of degree 120. The $H$-orbits on $C$ have sizes 24 (only $O_H$), 42, 84 and 168, and $120=24a+42b+84c+168d$ forces $a=5$, $b=c=d=0$. $\square$

*Checks [N]* (`p5_curve.py` (4)):
- $Q_H$ vanishes on $O_H$, to order 5 at the base 7-point (log-log slope 4.99);
- over a class, $\prod_HQ_H/S_6^5$ is constant on $\varphi(C)$ to $10^{-8}$.

**Consequences.**
- **Linear equivalences.** $2L_{60}\sim5O_H$ for all 30 subgroups $H$. Each class partitions $D_7$ into 15 sets $O_H$, and $\prod_HQ_H\equiv c\,S_6^5$ on $C$. This is "$120=5\cdot24$" behind the degree 60.
- **5-torsion.** $O_H-O_K$ is 5-torsion in $\mathrm{Jac}(C)$. It is nonzero for $O_H\ne O_K$: otherwise a degree-$\le24$ pencil would exist, against $\operatorname{gon}\ge25$.

## 7.7 The $\mu=90$ series revisited [P]

The family $Q_p\in\mathrm{Sym}^2V_4$ of §5.4, attached to $B+T$, has generic rank 3 or 4.
- **Rank 1** is the spinor model excluded by Theorem 7.4.
- **Rank 2** would put a $\mathbf 6$ into $H^0(B)$ or $H^0(B+T)$. Lemma 7.3 kills $p_2,p_3,p_5,p_6$ (and $p_7$ for $B+T$). The image would then lie in a finite set, or in the curve of ordered roots of $t^7+ut^3+v$. That curve covers the $s$-line, $s=u^7/v^4$, with degree 5040.
  - Its monodromies are a 7-cycle at $s=0$, an element of type $(4,3)$ at $s=\infty$, and a transposition at the one other zero of the discriminant $-v^2(6912u^7+823543v^4)$.
  - Riemann–Hurwitz gives $2g-2=5040\big(-2+\tfrac67+\tfrac{11}{12}+\tfrac12\big)=1380$, so $g=691>136$, and $C$ cannot map onto it.

## 7.8 The other rigid signatures

`twisted_survey.py` repeats §§7.2–7.5 on all 26 rigid faithful $A_7$-curves of genus $\le335$, against every eigenspace of every element with fixed points:

| signature | $g$ | best pencil $\le$ | source |
|---|---|---|---|
| $(2,4,7)$ | 136 | **42** | order-3 class of degree 60, an involution |
| $(3,3,7)$ | 241 | 48 | spin class of degree 60 with $V_4\subseteq H^0$ (a spinor model exists here), a 3-cycle |
| $(2,7,7)$, one curve | 271 | 84 | spin class of degree 90 |
| $(4,4,4)$ | 316 | 95 | order-6 class of degree 105 |
| $(2,5,7)$, $(2,6,7)$, $(2,7,7)$ | 199–271 | 104–108 | order-3 classes of degree 120 |
| $(3,3,5)$, $(3,3,6)$, $(3,4,4)$, $(3,4,5)$, $(3,4,6)$ | 169–316 | 192–299 | |

## 7.9 Elliptic subcovers [X]

**Proposition 7.10** (`elliptic_subcovers.py`, all four classes). Consider elliptic subcovers $C\to E$ with $H^1(E)=v\otimes M$ inside the $\rho$-isotypic part, $\rho\in\{15,21\}$. They have degree $\ge60$. Equality holds for the Galois quotients $C/A_5$ (two points fixed; $\rho=15$) and $C/L_2(5)$ ($\rho=21$), which are elliptic.

- **The lattice.** $H^1(C,\mathbb Z)$ is computed exactly from the dessin, with its cup product and $A_7$-action:
  - the 5040 triangles form a $\Delta$-complex over $0<1<\infty$;
  - a tree–cotree basis gives integral cocycles;
  - the Alexander–Whitney product is unimodular, and the character is $3(10+\overline{10})+2\cdot15+2\cdot21+4\cdot35$.
- **The planes.** $\rho=15,21$ have multiplicity 2, with multiplicity space $M$ of rank 2. For $v\in V_\rho(\mathbb Q)$, $\Lambda_v=(v\otimes M)\cap H^1(C,\mathbb Z)$ is a sub-Hodge lattice, i.e. $f^*H^1(E_v)$. Its degree is $n(v)=|\langle e_1,e_2\rangle|$ on a basis. If $\mathrm{End}(E_\rho)=\mathbb Z$, these are all elliptic subcovers in that part.
- **The search.** $n(v)=|\Phi(v)|/g(v)$, where:
  - $\Phi$ is the cup-product form on translates of a $K$-fixed vector ($K=A_5$, $L_2(5)$);
  - $g(v)$ comes from the finite glue, of squarefree exponent.

  Fincke–Pohst finds no $v$ with $|\Phi(v)|<60\,g(v)$.

So elliptic subcovers give only $\operatorname{gon}\le120$. Conversely, $\operatorname{gon}\ge25$ forces degree $\ge13$ for *every* elliptic subcover. That includes any arising from complex multiplication, from $E_{15}\sim E_{21}$, or from the parts $10\oplus\overline{10}$ and $35$.

## 7.10 Open

- **The exact gonality** in $[25,42]$. The $\tau$-pencil is exactly 42, elliptic subcovers give $\ge60$, and symmetric projections of $\varphi(C)$ give nothing better. A smaller pencil would have to be non-symmetric.
- **Lower bounds by degeneration.** Baker's specialisation lemma gives $\operatorname{gon}(C)\ge\operatorname{dgon}(\Gamma)$ for the dual graph $\Gamma$ of a stable reduction. At $p=7\,\|\,|A_7|$, Raynaud–Wewers describe the stable reduction of three-point covers. Could $\Gamma$, which carries an $A_7$-action, have divisorial gonality $>25$?
- **The ideal of $\varphi(C)$.** Is it generated by $X_3$ and the 15 quartics? Is $\varphi(C)$ projectively normal, i.e. $h^0(2L_{60})=21$ and $h^0(B+T)=10$? How does $\varphi(C)$ meet the 420 lines and sit in the hyperkähler Fano variety of $X_3$?
- **The kernel map** of the $\mu=90$ series: is it a $\mathbb P^3$-model, a spin class of degree 135 with $M^2=3(B+T)$?
- **Immersion at the 4-points.** Prove that $\varphi$ is an immersion at the points with stabiliser $C_4$, i.e. that the sequence there is $(0,1,2,3,4,6)$ rather than $(0,2,3,4,5,6)$. This would make Proposition 7.7 a paper proof. If it failed, $\operatorname{gon}(C/\tau)\le15$.
