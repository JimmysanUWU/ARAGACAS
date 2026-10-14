# Chapter 7. Schur-twisted geometry: the degree-60 model and $\operatorname{gon}(C)\le42$

**Idea.** An $A_7$-*invariant* line bundle need not be *linearised*.
- **The obstruction.** It is the bundle's Mumford class, in $H^2(A_7,\mathbb C^\*)=\mathbb Z/6$. The bundle is linearised for the Schur cover $6.A_7$, whose centre $Z$ acts on fibres by a character $\varepsilon$.
- **Degrees and sections.** Twisted classes exist in every degree in $15\mathbb Z$, against $90\mathbb Z$ for linearised ones. Holomorphic Lefschetz, twisted by $\varepsilon$, computes their sections.
- **The first twisted class with sections** has degree 60 and Mumford class of order 3. Its six sections embed $C$ as a degree-60 curve in $\mathbb P^5$, on a cubic fourfold. The normalization-defect bound and exact eigenline stabilizers prove the embedding (§7.12).
- **The pencil of degree 42.** An involution's $(-1)$-eigenspace in those six sections is two-dimensional, and its sections vanish at the 18 fixed points. What remains is a pencil of degree 42.
- **Lower bounds feed back.** The lower bound 25 of Chapters 2–3 is what keeps the twisted classes below degree 60 empty (Clifford).

Accessories over a linear base produce linearised bundles (Theorem 5.2). For smooth projective generically free bases with global equivariant Albanese vanishing, the stable formula admits exactly the twisted classes allowed by the Amitsur subgroup (Theorem 5.15); Example 5.12 realises $L_{60}$ over a rational base. The geometry of $C$ lives in the twisted classes.

**Main results** ($C$ a $(2,4,7)$ curve, $\tau$ an involution).
1. Invariant classes have degrees $15\mathbb Z$; linearised ones $90\mathbb Z$ (Cor. 7.2). Apart from the trivial bundle, none of degree $<60$ has sections (Thm 7.4; its Clifford step depends on the certified lower bound).
2. A degree-60 class of order 3 gives a $3.A_7$-equivariant closed embedding $\varphi:C\to\mathbb P^5$ [P][X]:
   - its image lies on the Laza–Zheng $A_7$-cubic fourfold;
   - it lies on no quadric [P][X] (§7.12); being cut out by the cubic and 15 quartics remains [N] (Prop. 7.8).
3. **$\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\tau)\le21$**, (Cor. 7.6 [P][X]). The pencil is base-point-free of degree exactly 42 by Prop. 7.7 [P][X]. Hence
   $$25\le\operatorname{gon}(C)\le42,\qquad13\le\operatorname{gon}(C/\tau)\le21,\qquad25\le\gamma(A_7)\le42 .$$
4. Each Klein subgroup $L_2(7)$ cuts $\varphi(C)$ by a quadric in exactly $5\times$ its 24 points of order 7 (Prop. 7.9).
   - That quadric is the Plücker quadric of the Weil representation of $SL_2(7)$.
   - The 24 points are the Veronese images of the flexes of Klein's quartic.
   - $6L_{60}\sim D_7$ exactly; identifying the cutting sextic with $\mathrm{Hess}(X_3)$ remains [N].
5. **Every elliptic subcover has degree $\ge60$** (Thm 7.21 [P][X]). Pure 15/21 planes are priced exactly (Prop. 7.10). Frobenius at $p=37,43,67$ shows that no other planes exist: $A$ and $S$ have no elliptic factors, and $E_{15}$, $E_{21}$ have no CM and are not isogenous (Prop. 7.20).
6. **Pencils by symmetry type.** The kernel $N$ of a gonal pencil's class stabiliser on the pencil costs degree $|N|\operatorname{gon}(C/N)$. This cost is $\ge42$ unless $N$ is trivial or one of five small groups (Prop. 7.11, Thm 7.13, Cor. 7.14 [P][N]).
   - The kernel $C_7$ is impossible [P][X][C] (Thm 7.22).
   - Free symmetries of the class divide the degree (Lemma 7.17). Symmetries acting on the pencil restrict the balanced energy to allowed isotypes (Thm 7.18).
7. **Orbifold descent** (§7.14 [P][X]).
   - A symmetric pencil descends to $C/K\to\mathbb P^1/(K/N)$, and inertia forces local degrees modulo rotation orders (Lemma 7.24). For example, faithful Klein symmetry forces odd degree.
   - In the record stratum, pencils are 3-torsion relations $(1-\iota)\Delta\sim\pm\theta$ on the genus-10 curve $C/U$ (Prop. 7.25).
   - Frobenius freezes the free point (Lemma 7.26). **No pencil of degree $<42$ has the symmetry of the $\tau$-pencil** (Thm 7.27).
   - The two curves are Galois conjugate, so every algebraic statement holds on all four classes.

Scripts: `twisted_rr.py` (25 s), `twisted_survey.py` (90 s), `tau_pencil.py` (1 min), `p5_curve.py` (7 min), `elliptic_subcovers.py` (80 s), `quotient_gonality.py` (9 min per class), `audit_exact.py` (8 s), `equivariant_spectral.py` (3 min; `--strata` 10 min), `frobenius_lpolys.py` (about 40 min for $p=37,43,67$), `orbifold_congruence.py` (3 min), `tau_stratum.py` (Sage; about 90 min at $p=43$); outputs `*_output.txt`.

## 7.1 The groups [X]

| group | construction | faithful irreducible degrees |
|---|---|---|
| $2.A_7\subset SU(4)$ | half-spin $V_4$ of the standard 6 ($\wedge^2V_4=\mathbf 6$); Clifford lifts | $4,4,14,14,20,20,36$ |
| $3.A_7$ | coset enumeration of the presentation below (order 7560, perfect) | $6,15,15,21,21,24,24$ per order-3 character |
| $6.A_7=2.A_7\times_{A_7}3.A_7$ | explicit 2-cocycles; 40 classes; Burnside–Dixon | $6,6,24,24,36$ per order-6 character |

$$3.A_7=\langle x,y,z\mid z^3,\ [x,z],\ [y,z],\ x^3,\ y^5,\ (xy)^7,\ (xyxy^{-1})^2=z,\ (xy^{-2}xy^2)^2\rangle\qquad(z=1\text{ gives }A_7).$$

**Exact check [X]** (`audit_exact.py`). The ATLAS matrices `3A7G1-Ar6B0` (over $\mathbb Z[\omega]$) and `2A7G1-Ar4aB0` (over $\mathbb Z[b_7]$) satisfy the relations up to scalars, and every Cayley-graph edge is central. So they give the $\mathbf 6$ and $V_4$ exactly. The same script computes the Molien series of §§7.4 and 7.6 exactly, and the complete central-3 character sector $6,15,15,21,21,24,24$.

## 7.2 The invariant Picard group [P][X]

Let $x_1x_2x_3=1$ be the branch generators (orders $e_i=2,4,7$), and $\hat x_i\in6.A_7$ fixed lifts. The *local datum* of a twisted class is the eigenvalue $\lambda_i=e^{2\pi ir_i/e_i}$ of $\hat x_i$ on the fibre at its fixed point. Define $\rho_0$ by $\varepsilon(\hat x_1\hat x_2\hat x_3)=e^{2\pi i\rho_0}$.

**Proposition 7.1.** For every $\varepsilon$ and every local datum with $\lambda_i^{e_i}=\varepsilon(\hat x_i^{e_i})$ there is an $\varepsilon$-twisted invariant line bundle, of degree
$$\deg L=2520\Big(n+\sum_i\frac{r_i}{e_i}-\rho_0\Big),\qquad n\in\mathbb Z .$$

*Proof.*
1. $\tilde\Delta=\langle c_i\mid c_1^2=c_2^4=c_3^7=c_1c_2c_3=:h\rangle$ is the centrally extended triangle group (Milnor): the preimage of $\Delta(2,4,7)$ in the universal cover of $\mathrm{PSL}_2(\mathbb R)$, a central extension by $\langle h\rangle\cong\mathbb Z$. It acts on $K^w$ over $\mathbb H$: $h$ acts by $e^{-2\pi iw}$, and $c_i$ acts at its fixed point by $e^{-2\pi iw/e_i}$.
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

**Theorem 7.4 [P][X][C].** An invariant class of degree $<60$ other than $\mathcal O$ has $h^0=0$.

*Proof.* $\operatorname{gon}\ge25$ (Chapters 2–3) gives $\mathrm{Cliff}\ge22$, since $\operatorname{gon}\le\mathrm{Cliff}+3$ (Coppens–Martens). $H^0$ is a sum of faithful irreducibles of the twist.
- **Degree 0 ($T$):** no sections.
- **Degree 15:** $h^0\le1<6$.
- **Degree 30:** Clifford gives $h^0\le5<6$.
- **Degree 45 (spin):** $h^0\le12$, so $H^0\supseteq V_4$ or $V_4^\*$. This gives $\psi:C\to\mathbb P^3$ of degree 45. Its normalized image has a faithful $A_7$-action; a nonbirational image would have degree $\le15$ and genus $\le\pi(15,3)=42<136$, so the map is birational.
  1. The invariants of $2.A_7$ on $V_4$ have Hilbert series $1+t^8+t^{12}+t^{14}+t^{16}+t^{18}+2t^{20}+\dots$
  2. Lemma 7.3 kills $f_{14}$ and $f_{18}$ on $C$: their least admissible divisors have degrees 3150 and 3330, against $14\cdot45=630$ and $18\cdot45=810$.
  3. Involutions lift with eigenvalues $(i,i,-i,-i)$, and $i^{14}=i^{18}=-1$, so $f_{14},f_{18}$ vanish on the 210 lines $\mathbb P(E_{\pm i}(\hat\tau))$. These are distinct: in a unitary model $E_{-i}=E_i^\perp$, so either line determines $\hat\tau$ up to sign.
  4. $f_{14},f_{18}$ are coprime. Their gcd is invariant, because $2.A_7$ is perfect and so has no characters. A nonconstant gcd would have degree 8, 12 or 14, and dividing $f_{14}$ or $f_{18}$ by it would give an invariant of degree 6, 2 or 4; there are none. So $Z(f_{14},f_{18})$ is a curve of degree 252 containing 210 lines, leaving degree 42 for $\psi(C)$, which has degree 45. $\square$ (Steps 3–4: GPT, `gpt/A7_Amitsur_Correspondence.pdf` §9.)

So there is **no degree-45 spinor model** $C\to\mathbb P(V_4)$. The local data of step 2 are checked exactly in `audit_exact.py`.

## 7.5 The degree-60 model and the $\tau$-pencil

**Theorem 7.5 [P][X]** (exact: `audit_exact.py`). Every $(2,4,7)$ curve carries an invariant class $L_{60}$ of degree 60 and Mumford class of order 3 with
$$\chi(L_{60})=\mathbf 6-15-2\cdot21-24,\qquad\text{so }H^0(L_{60})\supseteq\mathbf 6 .$$
The map $\varphi:C\to\mathbb P^5$ given by $\mathbf 6$ has these properties:
1. **Base-point-free.** An invariant base locus would have degree $\le60<360$.
2. **A closed embedding of degree 60.** First it is birational: a nonbirational image would have degree $\le30$ and normalized genus $\le\pi(30,5)=91$, whereas its $A_7$-action is faithful and requires genus $\ge136$. (A trivial image action would force the map degree to be at least 2520.) The normalization-defect argument of Proposition 7.15 then proves the embedding exactly.
3. **Equivariant** for $3.A_7$ on its exceptional $\mathbf 6$.
4. **Involutions.** The lift $\hat\tau$ has eigenvalues $(+1)^4(-1)^2$, and all 18 fixed points of $\tau$ have eigenvalue $+1$, so they lie in $\mathbb P(E_+)\cong\mathbb P^3$.
5. **Plücker [P][N].** Local characters force four even and two odd orders at the 2- and 4-points, with minimal sequence $(0,1,2,3,4,6)$ of weight 1, and six distinct orders mod 7 at the 7-points, with minimal sequence $(0,\dots,5)$. The minimal sequences carry weight 1890 out of the total $6(60+5\cdot135)=4410$. If they hold, as Proposition 7.8 indicates numerically, the remaining 2520 is one free orbit of simple flexes.

So **60 is the least degree of an invariant class with $h^0\ge2$** [P][X][C], against $\mu=90$ for linearised classes. The exact six-dimensional subsystem does not assert $h^0(L_{60})=6$.

**Corollary 7.6 [P][X].** $\operatorname{gon}(C)\le42$ and $\operatorname{gon}(C/\tau)\le21$.

*Proof.* Sections in the 2-dimensional $(-1)$-eigenspace $E_-$ satisfy $s(p)=-s(p)$ at the 18 fixed points. So $|E_-|$ has them as base points, and its moving part has degree $\le42$. Ratios of anti-invariant sections are $\tau$-invariant, so the pencil descends to $C/\tau$. $\square$

The inputs are certified in exact arithmetic over $\mathbb Q(\zeta_{84})$, on all four classes and both orientations (`audit_exact.py`):
- the multiplicity $\langle\chi(L_{60}),\chi_{\mathbf 6}\rangle=1$, with $\langle\chi(L_{60}),\chi(L_{60})\rangle=7$, at the local datum $r=(0,2,6)$;
- $\operatorname{tr}\hat\tau=2$;
- eigenvalue $+1$ at all 18 fixed points (since $\lambda_1=1$ and $\lambda_2^2=1$).

**Proposition 7.7 (exactly 42) [P][X].** The base locus of $|E_-|$ is exactly the 18 fixed points, each simple. So the bound comes from a base-point-free $g^1_{42}$, pulled back from a $g^1_{21}$ on $C/\tau$.

*Proof.* By Proposition 7.15, $\varphi$ is an embedding. A base point maps into $\mathbb P(E_+)$ and therefore satisfies $\varphi(\tau p)=\varphi(p)$, forcing $\tau p=p$. Conversely all 18 fixed points lie there by their fibre characters. At a fixed point, the tangent action is $-1$. In the ambient tangent space at $\mathbb P(E_+)$, the $-1$ directions are precisely those towards $E_-$. Immersion therefore gives an anti-invariant section with a nonzero first derivative, making the common zero simple. Removing the fixed divisor leaves degree $60-18=42$; ratios of the sections are $\tau$-invariant, giving quotient degree 21. $\square$

The numerical plane-intersection check in `tau_pencil.py` is now supplementary. Neither it nor the lower bound 25 is needed for this proposition.

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

**Proposition 7.8 (equations) [N], with no quadrics and embedding now [P][X]** (`p5_curve.py`, class 0; §7.12 for the exact statements).
1. **The 4-points.** On each line $\mathbb P(E_{-1}(\hat h))$, $G_6$ has exactly two double zeros, which are singular points of $\{G_6=0\}$. They are $\varphi(\mathrm{Fix}\,h)$.
2. **More equations.** The invariants of degrees 12, 18, 24 vanishing there, and those of degree 21 vanishing at a 7-point, vanish on $\varphi(C)$, by the exact orders above. With $F_3$, $G_6$ and the forced invariants, they cut out a smooth branch through each 7-point. Continuation and the group spread it over $\varphi(C)$.
3. **Hilbert function.** In degrees $2,3,4,5$ the image of $\mathrm{Sym}^k\mathbf 6$ has dimension $21,55,105,165$. So:
   - $\varphi(C)$ lies on no quadric [P][X] and on one cubic [N for uniqueness among all cubics];
   - it lies on 21 quartics, namely $F_3\cdot\mathbf 6$ and 15 more;
   - the reported ranks imply 4- and 5-normality [N]. Riemann--Roch gives $h^0(4L_{60})=105$ using Theorem 7.4 to kill the invariant degree-30 residual, so this equality inherits [C]. The degree-300 bundle is nonspecial and $h^0(5L_{60})=165$ unconditionally.
4. **Embedding [P][X].** Proposition 7.15 proves this independently of the numerical quintic restriction rank. The latter provides a supplementary check using very ampleness of $5L_{60}$.
5. **Cut out.** $X_3$ and the 15 quartics meet a random hyperplane in exactly 60 points, so they cut out $\varphi(C)$. A hyperplane of the $\tau$-pencil meets them in $18+42$ points.

So $\varphi(C)$ is a smooth curve of degree 60 and genus 136 [P][X]; being cut out by the Laza–Zheng cubic and 15 quartics remains [N].

**Proposition 7.9 (Klein quadrics) [P][X].** Let $H\cong L_2(7)$ be one of the 30 Klein subgroups of $A_7$, in two classes of 15.
- $\mathbf 6|_H\cong\mathrm{Sym}^2(3)$, so the $H$-invariant quadrics form a line $\langle Q_H\rangle$: the lift of Klein's quartic through $\mathrm{Sym}^4\subset\mathrm{Sym}^2\mathrm{Sym}^2$.
- Then $\operatorname{div}(Q_H|_C)=5\,O_H$, where $O_H$ is the set of 24 seven-points whose stabiliser lies in $H$.

*Proof.* Proposition 7.16 proves $Q_H|_C\ne0$ exactly, since $\varphi(C)$ lies on no quadric. Its divisor is $H$-invariant of degree 120. The $H$-orbits on $C$ have sizes 24 (only $O_H$), 42, 84 and 168, and $120=24a+42b+84c+168d$ forces $a=5$, $b=c=d=0$. $\square$

*Checks [N]* (`p5_curve.py` (4)):
- $Q_H$ vanishes on $O_H$, to order 5 at the base 7-point (log-log slope 4.99);
- over a class, $\prod_HQ_H/S_6^5$ is constant on $\varphi(C)$ to $10^{-8}$.

**Consequences.** The nonzero restriction is now exact (Prop. 7.16): both irreducible quadratic systems have empty projective base locus. The Hessian identification below remains [N].
- **Linear equivalences.** $2L_{60}\sim5O_H$ for all 30 subgroups $H$. Each class partitions $D_7$ into 15 sets $O_H$, and $\prod_HQ_H\equiv c\,S_6^5$ on $C$. This is "$120=5\cdot24$" behind the degree 60.
- **5-torsion.** $O_H-O_K$ is 5-torsion in $\mathrm{Jac}(C)$ [P][X]. Its nonvanishing for $O_H\ne O_K$ additionally uses [C]: otherwise a degree-$\le24$ pencil would exist, against $\operatorname{gon}\ge25$.
- **A sixth root of the 7-points.** [P] for $6L_{60}\sim D_7$; the last equivalence uses Prop. 7.9. $K_C\sim D_2+3D_4-8D_7$ and $K_C-3L_{60}=B+T=D_2-3D_4+2D_7$ (§7.6, Prop. 5.8) give $3L_{60}\sim2D_4-3D_7$. With $4D_4\sim7D_7$:
  $$6L_{60}\sim D_7\sim15\,O_H .$$
  The invariant sextic that cuts $D_7$ may be taken to be the Hessian of $X_3$, which is nonzero on $\varphi(C)$ [N]. So $\varphi(C)\cap\mathrm{Hess}(X_3)=D_7$. This parallels Klein's quartic $X$, whose 24 flexes are $X\cap\mathrm{Hess}(X)\sim6K_X$.
- **Plücker and Klein's quartic [P].** Let $W_H$ be the 4-dimensional Weil representation of $SL_2(7)$. Then $\wedge^2W_H\cong\mathbf 6|_H\cong\mathrm{Sym}^2(3_H)$, since $L_2(7)$ has a single irreducible of degree 6. So $\mathbb P^5$ carries two $H$-structures at once.
  - $Q_H$ is the Plücker quadric $G(2,W_H)$, since that quadric is $SL(W_H)$-invariant and the invariant quadric is unique. So $\varphi(C)$ meets the Grassmannian of lines of $\mathbb P(W_H)$ only in $O_H$, with contact 5 at each point.
  - The Veronese surface $V_H=v_2\mathbb P(3_H)$ meets $Q_H$ in $v_2(X_H)$, where $X_H$ is Klein's quartic. This is the construction of $Q_H$.
  - **$\varphi(O_H)=v_2(\mathrm{Flex}\,X_H)$** on one of the two Veronese surfaces of $H$ (those of $3_H$ and $\bar3_H$):
    1. A Sylow 7-subgroup $\langle c\rangle\le H$ has six eigenlines on $\mathbf 6$, and the three on a given Veronese surface are $v_2$ of its fixed points on $\mathbb P(3_H)$, which are flexes of $X_H$.
    2. The three points of $C$ fixed by $c$ are permuted by $N(\langle c\rangle)=7{:}3$, so the rotations of $c$ at them are $r,2r,4r$.
    3. The local datum of Theorem 7.5(5) (orders $0,\dots,5$ mod 7) sends the point with rotation $r$ to the eigenline of $\zeta^{\pm r}$. So the three points land on the eigenlines with exponents $\pm r\{1,2,4\}$, i.e. on one Veronese triangle.

    Numerically $\varphi(p)$ lies on $V_H$ to $7\cdot10^{-16}$, for both Klein subgroups containing $c$. So each 7-point of $C$ is a flex of two Klein quartics, one for each class.

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

So elliptic subcovers give only $\operatorname{gon}\le120$. Theorem 7.21 (§7.13) extends the bound 60 to *every* elliptic subcover. Frobenius shows that complex multiplication, $E_{15}\sim E_{21}$, and elliptic factors of the parts $10\oplus\overline{10}$ and $35$ do not occur.

## 7.10 Pencils by symmetry type [P][N]

**Idea.**
- **The stabiliser and its kernel.** The class of a pencil has a stabiliser $K\le A_7$, and $K$ acts on the pencil $\mathbb P^1$ through a finite subgroup of $PGL_2$.
- **Where symmetry costs degree.** The kernel $N$ of that action is the part that costs degree: the pencil descends to $C/N$.
- **Pricing the quotients.** So the question becomes the gonality of the quotients $C/K$. It is computed from holomorphic differentials on $C/K$, syzygies, and Castelnuovo–Severi through the subgroup lattice.
- **The answer.** Unless $N$ is trivial or one of six small groups, the cost is at least 42.

**Proposition 7.11 (class stabiliser) [P].** Let $f:C\to\mathbb P^1$ have degree $d=\operatorname{gon}(C)$, and $M=f^\*\mathcal O(1)$. Let $K=\{g:g^\*M\cong M\}$, and let $N\trianglelefteq K$ be the kernel of the action of $K$ on $|M|\cong\mathbb P^1$.
1. $K/N$ is a finite subgroup of $PGL_2(\mathbb C)$: cyclic, dihedral, $A_4$, $S_4$ or $A_5$.
2. $f$ factors through $C/N$, so $d\ge|N|\operatorname{gon}(C/N)$.
3. The orbit of $[M]$ consists of $[A_7:K]$ pencils of degree $d$.

*Proof.* $h^0(M)=2$, since otherwise $M(-x)$ would be a pencil of degree $d-1$. For $k\in K$ an isomorphism $k^\*M\cong M$ is unique up to scalar, so $K$ acts on $\mathbb P(H^0(M))$. An element of $N$ acts on $H^0(M)$ by a scalar, so it fixes $f=s_1/s_2$. $\square$

$K$ is the stabiliser of the pencil field $\mathbb C(f)$, the group $H$ of Lemma 4.3. That lemma uses parts 1–2 for $d\le24$; here they are applied for every $d<42$.

**Example (the record pencil) [P][X].** For the $\tau$-pencil of Proposition 7.7, $K=C(\tau)$, $N=\langle\tau\rangle$ and $K/N\cong C_2\times S_3\cong D_6$ (the group $H$ of Chapter 1).
- $g$ fixes the class $[L_{60}-\mathrm{Fix}\,\tau]$ iff $\mathrm{Fix}(g\tau g^{-1})\sim\mathrm{Fix}\,\tau$. For $g\tau g^{-1}\ne\tau$ the two sets are disjoint (stabilisers are cyclic), so their difference would be the divisor of a function of degree 18, against $\operatorname{gon}\ge25$. Hence $K=C(\tau)$.
- On $E_-(\tau)$, only $1$ and $\tau$ act by scalars [X], so $N=\langle\tau\rangle$.
- So $C_2\times S_3$ acts faithfully on the pencil, as a dihedral group of $\mathbb P^1$, as Proposition 7.11 requires.

**Lemma 7.12 (Castelnuovo–Severi on the subgroup lattice) [P].** Let $S<T\le A_7$ with $m=[T:S]$, and let $f$ be a pencil of degree $d$ on $C/S$. Then either
- $f$ factors through $C/\langle S,j\rangle$ for some $j\in T\smallsetminus S$, or
- $g(C/S)\le m\,g(C/T)+(m-1)(d-1)$.

*Proof.* If $(f,\pi):C/S\to\mathbb P^1\times C/T$ is birational onto its image, this is the Castelnuovo–Severi inequality. Otherwise it factors through a cover $C/S\to Z$ of degree $\ge2$ with $Z\to C/T$. By Galois theory for $C\to C/T$, $Z=C/K''$ with $S\subsetneq K''\le T$. Then $f$ factors through $C/\langle S,j\rangle$ for any $j\in K''\smallsetminus S$. $\square$

**Theorem 7.13 (quotient pencils) [P][N]** (`quotient_gonality.py`, classes 0 and 12). For every $K\le A_7$ with $|K|\le60$, $\operatorname{gon}(C/K)$ is at least the value in the table below. In particular $|K|\operatorname{gon}(C/K)\ge42$ unless $K$ is conjugate to $1$, $C_2$, $C_3$, $C_4$, $V_4$, $S_3$ or $C_7$.

**The inputs [N].**
- **Differentials on $C/K$.** Holomorphic differentials on $C$ are vector-valued weight-2 forms for $\Delta(2,4,7)$, one solve per irreducible in $H^0(K_C)$, with gaps $>10^8$. Pairing them with $K$-fixed vectors gives the differentials on $C/K$, sampled over the whole curve.
- **Noether.** $\mathrm{Sym}^2$ has rank $2g-1$ if $C/K$ is hyperelliptic and $3g-3$ otherwise.
- **Petri.** The quadrics through the canonical curve have Jacobian rank $g-3$ if it is trigonal or a plane quintic ($g=6$, which does not occur here), and $g-2$ otherwise.
- **Green–Lazarsfeld nonvanishing.** $\operatorname{Cliff}\le p$ forces $K_{p,2}\ne0$. So $K_{p,2}=0$ gives $\operatorname{gon}\ge p+3$.

These give lower bounds for ten subgroups: $\operatorname{gon}\ge3$ for $3^2{:}2$; $\ge4$ for $5{:}4$, $D_{12}$, $3{:}4$, two $A_4$ and $C_6\times C_2$; $\ge5$ for $D_{10}$ and $C_3^2$ ($K_{2,2}=0$); and $\ge6$ for $D_8$ ($K_{3,2}=0$).

**The propagation [P].** Start from $|K|\operatorname{gon}(C/K)\ge\operatorname{gon}(C)\ge25$ and propagate with Lemma 7.12 over the 36 classes of subgroups of order $\le120$ until nothing changes. The result:

| $\lvert K\rvert$ | $K$ | $g(C/K)$ | $\operatorname{gon}\ge$ | $\lvert K\rvert\operatorname{gon}\ge$ |
|---|---|---|---|---|
| 2 | $C_2$ | 64 | 13 | 26 |
| 3 | $C_3$ (two classes) | 46 | 9 | 27 |
| 4 | $C_4$ | 32 | 9 | 36 |
| 4 | $V_4$ (two classes) | 28 | 7 | 28 |
| 5 | $C_5$ | 28 | 9 | 45 |
| 6 | $C_6$ | 22 | 8 | 48 |
| 6 | $S_3$, orbits $3+2$ / $3+3$ | 19 | 6 / 5 | 36 / 30 |
| 7 | $C_7$ | 19 | 4 | 28 |
| 8 | $D_8$ | 12 | 6 | 48 |
| 9 | $C_3^2$ | 16 | 6 | 54 |
| 10 | $D_{10}$ | 10 | 5 | 50 |
| 12 | $A_4$ (four classes), $3{:}4$, $D_{12}$, $C_6\times C_2$ | 7–11 | 4–6 | $\ge48$ |
| 18, 20, 21 | $3^2{:}2$, $5{:}4$, $7{:}3$ | 4, 5, 7 | 3, 4, 2 | 54, 80, 42 |
| 24–60 | $S_4$ (four), $C_3\rtimes D_8$, $3^2{:}4$, $3\times A_4$, $A_5$ (two) | 1–4 | 2–3 | $\ge48$ |

The script prints every row. Class 12 gives the same table.

**Corollary 7.14 [P][N].** If $\operatorname{gon}(C)\le41$, every pencil computing it has kernel $N$ conjugate to one of
$$1,\quad C_2,\quad C_3\ (\text{two classes}),\quad C_4,\quad V_4\ (\text{two classes}),\quad S_3\ (\text{two classes}),\quad C_7 .$$
Moreover $K/N\subset PGL_2$. For $N=C_7$ this forces $K\le7{:}3$. For $N=1$, the stabiliser $K$ itself acts faithfully on the pencil, and the pencil descends to one of the same degree on $C/K$ whose fibres over the branch values of $\mathbb P^1\to\mathbb P^1/K$ are divisible by the branching orders, away from points with nontrivial stabiliser.

Theorem 7.22 removes $C_7$ from this list (§7.13). So a pencil beating 42 either is pulled back from $C/N$ for one of the five remaining nontrivial groups, or every symmetry of its class acts faithfully on it. The record pencil lives in the stratum $N=C_2$, where 42 is optimal iff $\operatorname{gon}(C/\tau)=21$.

## 7.11 Open

- **The exact gonality** in $[25,42]$.
  - The $\tau$-pencil has degree exactly 42 [P][X] (Prop. 7.7). Every elliptic subcover gives $\ge60$ (Thm 7.21), and symmetric projections of $\varphi(C)$ give nothing better.
  - By Corollary 7.14 and Theorem 7.22, a smaller pencil has kernel $1$ or one of five small groups. $C_7$ is closed (Thm 7.22).
  - Closing $S_3$ needs $\operatorname{gon}(C/S_3)\ge7$, e.g. $K_{4,2}=0$ on the genus-19 quotients; only $d=30,36$ remain (Cor. 7.23).
  - $C_2$ needs $\operatorname{gon}(C/\tau)=21$. Its top stratum $K=C(\tau)$ is closed: the $\tau$-pencil is optimal there (Thm 7.27). What remains is $K\subsetneq C(\tau)$, with the degrees of Lemma 7.24.
  - The orbifold descent of §7.14 makes every stratum with small $g(C/K)$ a finite computation on an explicit curve. The three $S_4$-strata of genus 3 ($N=1$, $d\in\{27,33,39\}$; $N=C_3$, $K/N\cong D_8$; $N=V_4$, $K/N\cong S_3$) are next.
  - The $[P][C]$ list of Corollary 7.23(2) shows which numerical quotient bounds of Theorem 7.13 matter: those for $C_5$, $C_6$, $D_8$, $D_{10}$, $D_{12}$ and $5{:}4$.
- **Lower bounds by degeneration.** Baker's specialisation lemma gives $\operatorname{gon}(C)\ge\operatorname{dgon}(\Gamma)$ for the dual graph $\Gamma$ of a stable reduction. At $p=7\,\|\,|A_7|$, Raynaud–Wewers describe the stable reduction of three-point covers. Could $\Gamma$, which carries an $A_7$-action, have divisorial gonality $>25$?
- **The ideal of $\varphi(C)$.** Is it generated by $X_3$ and the 15 quartics? Is $\varphi(C)$ projectively normal, i.e. $h^0(2L_{60})=21$ and $h^0(B+T)=10$? How does $\varphi(C)$ meet the 420 lines and sit in the hyperkähler Fano variety of $X_3$?
- **The kernel map** of the $\mu=90$ series: is it a $\mathbb P^3$-model, a spin class of degree 135 with $M^2=3(B+T)$?
- **Full branch vanishing sequences.** Immersion at all 4-points and absence of extra pencil base points are proved in §7.12. This gives first orders $(0,1)$, not the entire minimal vanishing sequence; the remaining orders and free flex orbit of Theorem 7.5(5) still require proof.

## 7.12 Normalization defects and the two quadratic systems [P][X]

Proofs: `gpt/A7_Normalization_Defects.md`, `gpt/A7_Quadratic_Systems.md`. Exact checks: `gpt/verify_normalization.py`, `gpt/verify_quadrics.py`. These replace the embedding and no-quadric numerical inputs above.

**Proposition 7.15 (embedding).** The morphism $\varphi$ of Theorem 7.5 is a closed embedding on all four classes.

*Proof.* For $Y=\varphi(C)$, Castelnuovo bounds the arithmetic genus by $\pi(60,5)=406$. Hence the total normalization defect is $D=p_a(Y)-136\le270$. Locally, a fibre with $b$ branches, $k$ of them nonimmersed, contributes at least $b-1+k$ to the defect. Thus the nonimmersion locus has at most $D$ points, and the collision locus at most $2D$ points. They are $A_7$-invariant. Every point orbit has size at least 360, so nonimmersion is impossible; every orbit except $D_7$ has size at least 630, so collisions could occur only within $D_7$. Each order-7 eigenline has stabilizer exactly $C_7$, so the map on the transitive orbit $D_7=A_7/C_7$ is injective. All collisions are excluded. $\square$

*The eigenline stabilizers, by hand.* A stabilizer of an eigenline $\ell$ of $\hat c$ contains $c$, so it is $C_7$, $7{:}3$, one of the two $L_2(7)$, or $A_7$ (the transitive subgroups of $A_7$ containing a 7-cycle).
- An element of $7{:}3\smallsetminus C_7$ conjugates $c$ to $c^2$ or $c^4$. It therefore maps the $\zeta^k$-eigenline to the $\zeta^{2k}$- or $\zeta^{4k}$-eigenline. The six eigenvalues of $\hat c$ are distinct, so no eigenline is fixed.
- $\mathbf 6$ is irreducible on both Klein subgroups and on $A_7$, so neither fixes a line.

The exact ATLAS specialization at $p=43$ (`gpt/verify_normalization.py`) confirms this.

**Proposition 7.16 (quadrics and apolarity).** Write $\operatorname{Sym}^2\mathbf6=Q_6\oplus Q_{15}$, with $Q_6\cong\mathbf6^*$.
1. Both quadratic systems have empty projective base locus. Thus every nonempty $A_7$-invariant closed subset of the exceptional $\mathbb P^5$, including $\varphi(C)$, lies on no quadric.
2. $Q_6$ is the Jacobian system of the invariant cubic. Its quotient ring has Hilbert series $(1+t)^6$.
3. $Q_{15}$ generates the apolar ideal of the dual invariant cubic, with quotient Hilbert function $(1,6,6,1)$.

*Proof.* The exact character projector $P_6=\frac6{2520}\sum_g\operatorname{tr}(M_g)\operatorname{Sym}^2(M_g)$ and $1-P_6$ give the two systems. Modular integer rank checks at $p=43$, for both embeddings of $\mathbb Z[\omega]$ and both dual conventions, prove $Q_6S_5=S_7$ (rank 792) and $Q_{15}S_2=S_4$ (rank 126). Full rank after specialization proves full rank in characteristic zero. Each system therefore has empty base locus; any nonzero invariant quadratic ideal would contain one of them.

Differentiating the invariant cubic identifies its six partials with $Q_6$, a regular sequence by the empty base locus. Differentiating the dual cubic twice identifies its apolar quadrics with $Q_{15}$. The degree-three multiplication rank is exactly 55: specialization bounds it below, and the invariant cubic constituent is absent from $Q_{15}\otimes\mathbf6$. Degree four vanishes in the quotient. Hence the apolar ideal is generated by these quadrics, with Hilbert function $(1,6,6,1)$. $\square$

The same script checks that the six restricts irreducibly to each Klein subgroup and has exactly one invariant quadric. Proposition 7.9 is consequently exact. The identification of the Hessian section, the proposed quartic ideal, and projective normality remain separate.

## 7.13 Symmetric pencils, elliptic subcovers and the kernel $C_7$ [P][X][C]

**Idea.**
- **Symmetries of a pencil's class cost degree.** A free cyclic symmetry linearises the pencil bundle, which then descends; so its order divides the degree. A symmetry that acts on the pencil can be built into Hersch's balancing, which confines the balanced map to the isotypes that the symmetry allows.
- **Elliptic subcovers are expensive.** Frobenius at $p=37,43,67$ shows that the abelian varieties $A$ and $S$ of $\mathrm{Jac}(C)\sim A^{10}E_{15}^{15}E_{21}^{21}S^{35}$ have no elliptic factors, and that $E_{15}$, $E_{21}$ have no CM and are not isogenous. So every elliptic subcover lies in a pure plane of Proposition 7.10, and costs at least 60.
- **The kernel $C_7$ is impossible.** A pencil of degree $\le41$ pulled back from $C/C_7$ would make $C/C_7$ bielliptic, with an elliptic subcover of degree 14.

Throughout, $f:C\to\mathbb P^1$ is a gonal pencil of degree $d$, $M=f^*\mathcal O(1)$, and $K\supseteq N$ are its class stabiliser and kernel (Proposition 7.11), with $\rho:K/N\hookrightarrow PGL_2$.

**Lemma 7.17 (free divisibility) [P].** Let $F\le K$ act freely on $C$, with $\rho(F)$ cyclic. Then $|F|$ divides $d$. In particular:
- $3\mid d$ if $K$ contains an element of order 3;
- $9\mid d$ if $9\mid|K|$;
- $5\mid d$ if $5\mid|K|$.

*Proof.* The Mumford class of $M$ on $F$ is the class of the projective representation $\rho|_F$ on $H^0(M)\cong\mathbb C^2$. That representation factors through the cyclic group $\rho(F)$, so it lifts to $GL_2$, and $M$ is $F$-linearisable. A linearised bundle descends along the étale quotient $C\to C/F$, so $|F|\mid\deg M$.

For the special cases: elements of orders 3 and 5 act freely (Corollary 1.2). A subgroup of order 9 is $C_3^2$, all of whose elements have order 3. Its image in $PGL_2$ is abelian of odd order, hence cyclic. $\square$

**Theorem 7.18 (equivariant Li–Yau) [P].** Let $R:K\to SO(3)$ be $\rho$ conjugated into $SO(3)$. Let $\Lambda$ be the least nonzero eigenvalue of $C$ whose eigenspace contains an irreducible $\sigma$ with $\mathrm{Hom}_K(R,\sigma)\ne0$. Then $d\ge\Lambda\,\mathrm{Area}(C)/8\pi=67.5\,\Lambda$.

*Proof.*
1. **Invariant measure.** $\mu=f_*dA$ has no atoms. It is $\rho(K)$-invariant, since $f\circ k=\rho(k)\circ f$ and $k$ is an isometry.
2. **Equivariant balancing.** The Douady–Earle barycentre of $\mu$ is unique and conformally natural, so $\rho(K)$ fixes it. A Möbius $\gamma$ moving it to $0$ balances $f$ and conjugates $\rho(K)$ into $\mathrm{Stab}(0)=SO(3)$. Then $x=\gamma\circ f$ has $\int x=0$ and $x\circ k=R(k)\,x$.
3. **Isotypes.** $\ell\mapsto\ell(x)$ is a $K$-map from $\mathbb R^{3*}$ to $L^2(C)$. Its $\sigma$-isotypic component vanishes unless $\mathrm{Hom}_K(R,\sigma)\ne0$.
4. **Li–Yau.** Hence $8\pi d=E(x)\ge\Lambda\|x\|^2=\Lambda\,\mathrm{Area}$. $\square$

**Computing the multiplicities [X].**
- $m_\sigma(K,N)=\langle\chi_R,\sigma|_K\rangle$ depends only on $N$. The reason: $\chi_R(k)=1+2\cos\theta$ depends only on the order of $kN$, except at orders 5 and 7, where the angle is not determined.
- At those two orders the mean of $1+2\cos\theta$ over a cyclic group is $\tfrac12$ (order 5) and $\tfrac23$ (order 7), also on each class of 7-cycles. $A_7$-characters are constant there ($10+\overline{10}$ taken together).
- `equivariant_spectral.py` lists $m_\sigma$ for all 40 classes of subgroups $K$ and all admissible $N$.

**Corollary 7.19 [P][C]** (classes 0, 1; spectral inputs of §2.5 and $Q_2$). The only nonzero eigenvalue below 0.55998 is $\lambda_1$, with $E_1\cong14_a$. Every nonzero eigenvalue below 0.689027 has its eigenspace in the $14_a$-, $14_b$- and $21$-isotypes.
- **$K\cong A_5$ acting faithfully** (either class): $m_{14_a}=0$ gives $d\ge38$, and $15\mid d$. So $d\ge45$.
- **$K=(A_4\times3){:}2$ with $N=C_3$ ($K/N\cong S_4$):** $m_\sigma=0$ for $\sigma=\mathbf6,14_a,14_b,21$. So $d\ge\lceil67.5\cdot0.689027\rceil=47$, and $9\mid d$. So $d\ge54$.
- For small $K$, $m_{14_a}>0$ and the theorem gives no more than Chapter 3.

**Proposition 7.20 (Frobenius on the factors) [X]** (`frobenius_lpolys.py` → `frobenius_lpolys_output.txt`). On both curves (both signs of $\sqrt{21}$ in the septic Belyi model of §6.5):
1. $A$ and $S$ have no elliptic isogeny factor over $\mathbb C$;
2. $E_{15}$ and $E_{21}$ have no complex multiplication;
3. $E_{15}\not\sim E_{21}$.

*Method.*
1. **Arithmetic monodromy.** Over $\mathbb F_p$ with $21$ and $-7$ squares, $\mathrm{disc}(\beta-T)=-7r^{-6}T^2(T-1)^4$ is a square. So the Galois closure of $\beta$ has arithmetic monodromy $A_7$, and each $C/H$ is defined over $\mathbb F_p$.
2. **Unramified fibres.** The fibre over $T$ contributes the permutation character of $\mathrm{Frob}_T$, read off from the factorisation pattern of $\beta(z)-T$.
3. **Ramified fibres.** Over $T=0,1,\infty$ the inertia generator and a Frobenius lift act on the Puiseux branches $z=r+c\,(T-T_b)^{1/e}+\cdots$. The fibre contributes the Frobenius-stable inertia orbits on $H\backslash A_7$.
4. **The L-polynomials.**
   - $E_{15}=C/A_5$, $E_{21}=C/L_2(5)$: from counts over $\mathbb F_p$.
   - $S=\mathrm{Jac}(C/(3^2{:}4))$ (genus 2): from counts over $\mathbb F_p,\mathbb F_{p^2}$.
   - $A$: since $\mathrm{Jac}(C/(7{:}3))\sim A\times E_{15}\times E_{21}\times S$ (genus 7), from counts over $\mathbb F_{p^k}$, $k\le3$, and the functional equation.

*Checks.*
- The genus-0 quotients $C/A_6$, $C/S_5$, $C/L_2(7)$ have exactly $p^k+1$ points.
- All roots lie on $|x|=\sqrt p$.
- $L_A$ splits over $\mathbb Q(\sqrt{-7})$ into two cubics, as $L(10)L(\overline{10})$ must.

*Proof.*
- **(1) No elliptic factors.** At $p=37$, the Weil polynomials of $S$ (degree 4) and $A$ (degree 6) have only $\mathbb Q$-irreducible factors of degree $\ge3$. No ratio of two distinct roots is a root of unity: the resultant $\mathrm{Res}_x(P(x),P(yx))/(y-1)^{2g}$ has no cyclotomic factor. So $[\mathbb Q(\alpha^k):\mathbb Q]\ge3$ for every root $\alpha$ and every $k$. By Honda–Tate, no base change of the reduction has an elliptic factor. An elliptic factor over $\overline{\mathbb Q}$ would survive good reduction.
- **(2) No CM.** A CM curve has Frobenius field $K$ at every ordinary prime and is supersingular only at primes not split in $K$ (Deuring).
  - Curve $+$, $E_{15}$: ordinary at 37 and 43 with field $\mathbb Q(\sqrt{-3})$, but supersingular at 67, which splits in $\mathbb Q(\sqrt{-3})$.
  - Curve $-$, $E_{15}$: supersingular at 37, which splits in the field $\mathbb Q(\sqrt{-34})$ of $p=43$.
  - $E_{21}$ on both curves: different fields at 37 and 43.
- **(3) Not isogenous.** At $p=37$: on curve $+$ the fields are $\mathbb Q(\sqrt{-3})$ and $\mathbb Q(\sqrt{-11})$; on curve $-$, $E_{15}$ is supersingular and $E_{21}$ ordinary. $\square$

| curve | $p$ | $a_p(E_{15})$ | $a_p(E_{21})$ | $L_S$ | $L_A$ |
|---|---|---|---|---|---|
| $+$ | 37 | $-10$ | $-7$ | $1+9T+48T^2+\dots$ | $1-2T-5T^2+228T^3+\dots$ |
| $+$ | 43 | $-8$ | $-11$ | $1+3T+78T^2+\dots$ | $1-14T+70T^2-240T^3+\dots$ |
| $+$ | 67 | $0$ | $-12$ | $1-12T+158T^2+\dots$ | $1+10T-30T^2-760T^3+\dots$ |
| $-$ | 37 | $0$ | $-12$ | $1+8T+50T^2+\dots$ | $1+6T+62T^2+402T^3+\dots$ |
| $-$ | 43 | $6$ | $-6$ | $1+8T+50T^2+\dots$ | $1+0T-7T^2+432T^3+\dots$ |
| $-$ | 67 | $2$ | $-7$ | $1-3T+80T^2+\dots$ (splits mod 67) | $1-8T+81T^2-640T^3+\dots$ |

$E_{15}$ on curve $-$ is also non-CM by its two ordinary fields $\mathbb Q(\sqrt{-34})$ ($p=43$) and $\mathbb Q(\sqrt{-66})$ ($p=67$). The reduction of $S$ splits at 67, which reduction may do. The conclusion uses $p=37,43$, where $L_S$ is irreducible. Not every prime is good. At $p=79$ the reduction of $A$ on curve $+$ acquires an elliptic factor: $L_A$ has a root in $\mathbb Q(\sqrt{-7})$, and two roots differ by the factor $-1$. Reduction can only add endomorphisms, so one good prime suffices.

**Theorem 7.21 (elliptic subcovers) [P][X].** Every nonconstant map from $C$ to an elliptic curve has degree $\ge60$. Equality holds for the Galois quotients $C\to C/A_5$ and $C\to C/L_2(5)$.

*Proof.* Let $H_Z\subset H^1(C,\mathbb Q)$ be the pullback of $H^1(Z)$, a simple rank-2 Hodge structure.
1. **Isotypic projections.** The projections to the isotypic parts $10+\overline{10}$, $15$, $21$, $35$ are Hodge morphisms, so each is injective or zero.
2. **Ruling out $A$, $S$ and mixing.** A nonzero projection to the $10+\overline{10}$- or $35$-part would give an elliptic factor of $A$ or $S$. Nonzero projections to both the 15- and the 21-part would give $E_{15}\sim Z\sim E_{21}$. Proposition 7.20 excludes all three.
3. **A pure plane.** So $H_Z$ lies in $V_\rho\otimes M_\rho$ for one $\rho\in\{15,21\}$. Since $\mathrm{End}(E_\rho)=\mathbb Z$, the sub-Hodge structures of $V_\rho\otimes M_\rho$ isomorphic to $M_\rho$ are exactly the $v\otimes M_\rho$ with $v\in V_\rho(\mathbb Q)$.
4. **The degree.** The pullback of $H^1(Z,\mathbb Z)$ is a sublattice of $\Lambda_v$. The degree of the map is its cup product, a multiple of $n(v)\ge60$ (Proposition 7.10). $\square$

This replaces the earlier "$\ge13$ for every elliptic subcover", which came from $\operatorname{gon}\ge25$.

**Theorem 7.22 (no kernel $C_7$) [P][X][C].** If $\operatorname{gon}(C)\le41$, no gonal pencil has kernel $C_7$.

*Proof.* Suppose $N=\langle c\rangle\cong C_7$, and write $f=\bar f\circ\pi$ with $\pi:C\to Y=C/C_7$, $g(Y)=19$, and $\deg\bar f=e$.
1. **The degree and the stabiliser.** From $25\le d\le41$ and $7\mid d$, we get $d\in\{28,35\}$. $N\trianglelefteq K$ gives $K\le N_{A_7}(C_7)=7{:}3$. If $K=7{:}3$, Lemma 7.17 gives $3\mid d$, which is impossible. So $K=C_7$.
2. **A second pencil on $Y$.** Let $s\in7{:}3$ have order 3. It descends to $\bar s\in\mathrm{Aut}(Y)$. Since $s\notin K$, the pencil fields $\mathbb C(\bar f)$ and $\mathbb C(\bar f\circ\bar s)$ are distinct.
3. **$e=5$ is impossible.** The compositum then has index 1, and Castelnuovo–Severi gives $19\le16$.
4. **$e=4$ makes $Y$ bielliptic.**
   - Index 1 would give $19\le9$, so the index is 2. The compositum is the function field of a double quotient $Z$ of $Y$, carrying two distinct degree-2 pencils, so $g(Z)\le1$.
   - $g(Z)=0$ would make $Y$ hyperelliptic and $\operatorname{gon}(C)\le14$. So $Z$ is elliptic, and $C\to Z$ has degree 14.
5. **Theorem 7.21 gives a contradiction**, since $14<60$. $\square$

*A direct finish, without Theorem 7.21.* Only (1) and (3) of Proposition 7.20 are needed.
- **The symmetry acts on $Z$.** The bielliptic involution $\iota$ of $Y$ is unique: two would give $19\le5$ by Castelnuovo–Severi. So $\bar s$ commutes with $\iota$.
- **It acts by a translation.** $\bar s$ acts freely on $Y$, because every element of $7{:}3\smallsetminus C_7$ has order 3 and acts freely on $C$. So the induced automorphism of $Z$ has order 3 and no fixed point. It is a translation, and acts trivially on $H^1(Z)$.
- **$H_Z$ is $(7{:}3)$-invariant.** Each of $10$, $\overline{10}$, $15$, $21$, $35$ has a one-dimensional $(7{:}3)$-fixed space. So $H^1(C/(7{:}3))$ is the sum of the $(7{:}3)$-invariants of the isotypic parts, of ranks 6, 2, 2, 4, and $H_Z$ lies in it.
- **The rank-2 parts are pure planes.** The 15- and 21-parts there have rank 2, so each is a pure plane $v_0\otimes M_\rho$ of degree $\ge60$ by Proposition 7.10. The other parts are excluded as in Theorem 7.21.

**Corollary 7.23 (the kernels that remain).**
1. **[P][N]** In Corollary 7.14, $C_7$ is removed: $N\in\{1,C_2,C_3,C_4,V_4,S_3\}$. These are six isomorphism types and nine conjugacy classes, since $C_3$, $V_4$ and $S_3$ have two classes each.
2. **[P][C]** Without the numerical inputs of Theorem 7.13, using only $\operatorname{gon}(C)\ge25$, Lemma 7.12, Lemma 7.17, Theorem 7.18 and Castelnuovo–Severi on the conjugate pencils of $C/N$ (`equivariant_spectral.py --strata`), the possible kernels of a gonal pencil of degree $\le41$ are:
   - $1,C_2,C_3,C_4,V_4$ (all degrees allowed by $|N|\operatorname{gon}(C/N)$);
   - $C_5$ ($d=30,35,40$); $C_6$ ($d=36$); $S_3$ ($d=36$ for orbits $3+2$; $d=30,36$ for orbits $3+3$); $D_8$ ($d=32,40$); $D_{10}$ ($d=30,40$, only with $K=5{:}4$); $D_{12}$ ($d=36$, only with $K/N=C_2$); $5{:}4$ ($d=40$).
   - Theorem 7.21 also removes the elliptic-subcover branches of the enumeration. This excludes $D_{10}$ with $K=N$ at $d=40$, in addition to $C_7$.

*Classes 12, 14.* There the spectrum is less favourable. $\lambda_1$ (on $14_a$) is followed by an eigenvalue near $0.38$ on the $21$-isotype [N] (§2.5), so Corollary 7.19 cannot be proved there directly. But the two curves are Galois conjugate (§7.14, Remark), so its algebraic conclusions ($d\ge45$ for faithful $A_5$, $d\ge54$ for $(A_4\times3){:}2$ with $N=C_3$) hold on all four classes, as do Lemma 7.17, Theorems 7.21–7.22 and Corollary 7.23.

## 7.14 Orbifold descent of symmetric pencils; the $\tau$-pencil is optimal in its stratum [P][X]

**Idea.**
- **Descend to $C/K$.** A pencil with class stabiliser $K$ and kernel $N$ becomes a map $h:C/K\to\mathbb P^1/(K/N)$, and $C/N\to C/K$ is the pull-back of the Galois cover $\mathbb P^1\to\mathbb P^1/(K/N)$. Every inertia point therefore forces its local degree modulo the order of a rotation. Lemma 7.17 is the case without inertia.
- **Dihedral norms.** In the stratum of the record pencil, $w=f^6$ lives on the genus-10 curve $X_{\rm rot}=C/U$, and its divisor is $3(\Delta-\iota\Delta)$ for an effective $\Delta$ of degree $d/6$. So a pencil of the stratum is a solution of $(1-\iota)\Delta\sim\pm\theta$ for a fixed 3-torsion class $\theta$. The $\tau$-pencil is the case $\deg\Delta=7$.
- **Frobenius freezes the free point.** Modulo a good prime, Frobenius fixes $\theta$, and Castelnuovo–Severi on $X_{\rm rot}\to T$ forbids small functions. So the one free point of a degree-30 candidate must be rational over the residue field. A one-parameter search becomes a finite one.

### Orbifold congruences

Let $f$ be a gonal pencil of degree $d$ with class stabiliser $K$ and kernel $N$ (Proposition 7.11). Let $\bar K=K/N$ act on $\mathbb P^1$ through $\rho$, $Y=C/N$, $Z=C/K$ and $k=d/|N|$. Then $f$ induces $h:Z\to\mathbb P^1/\bar K$ of degree $k$. Over a branch value $b$ of $\mathbb P^1\to\mathbb P^1/\bar K$, let $r_b$ be the rotation by $2\pi/n_b$ about a point above $b$. For $P\in Z$, let $\gamma_P\in\bar K$ be the canonical generator of inertia of $Y\to Z$ at a point above $P$, which acts on the tangent space by $e^{2\pi i/|\gamma_P|}$.

**Lemma 7.24 (orbifold congruences) [P].** If $h(P)=b$ with local degree $e$, then $\rho(\gamma_P)$ is conjugate to $r_b^{\,e}$. Hence, for every branch value $b$,
$$k=\sum_{h(P)=b}e_P,\qquad e_P\equiv e(\gamma_P)\pmod{n_b},$$
where $e(\gamma)$ is the exponent with $\rho(\gamma)\sim r_b^{e(\gamma)}$, and $e(1)=0$.

*Proof.* The map $Y\to Z\times_{\mathbb P^1/\bar K}\mathbb P^1$, $y\mapsto(\bar y,f_Y(y))$, is $\bar K$-equivariant and generically injective, because $\rho$ is faithful. So $Y$ is the normalised pull-back. Locally, $u^m=t$ on $Y$, $v^n=s$ on $\mathbb P^1$ and $s=t^e$. Then $f_Y$ is $v=c\,u^{me/n}$, and $u\mapsto e^{2\pi i/m}u$ maps to $v\mapsto e^{2\pi ie/n}v$. $\square$

**Computing it [X]** (`orbifold_congruence.py` → `orbifold_congruence_output.txt`).
- The points of $C$ with nontrivial stabiliser are the cosets $g\langle c\rangle$ for $c$ of order 2, 4, 7, with inertia generator $gcg^{-1}$. This $A_7$-set is the same on all four classes.
- For each stratum $(K,N)$ of Corollary 7.23(2), the script lists the inertia classes of $Y\to Z$ and their admissible pairs (branch value, residue): the rotation axes of $\rho(\gamma)$ in the cyclic, dihedral, $A_4$, $S_4$ and $A_5$ cases. A dynamic programme then keeps the $k$ that satisfy all the congruences and $k\ge\sum$ (least residues).

New exclusions, all on all four classes:
- **Faithful Klein symmetry forces odd degree.** If $N=1$ and $K\supseteq V_4$, then $d$ is odd. Each of the three involutions is a half-turn with its own axis and fixes 18 points, giving 9 points of $Z$ with odd local degree over its branch value. The same holds for faithful $D_8$, and for $N=\langle\tau\rangle$ with $K/N\cong V_4$ ($d\equiv2\pmod4$).
- **Faithful $A_4$, $D_6$, $S_4$** allow only $d\in\{27,33,39\}$. The same holds for $N=C_3$ with $K/N\cong V_4$ or $D_8$. For $N=C_3$ with $K/N\cong A_4$, only $d=27$ remains.
- **Kernel $C_5$** survives only with $K=C_5$. For $K=D_{10}$ and $K=5{:}4$, the 18 (resp. $8+2$) inertia points need $k\ge9$, while $k\le8$. The same mechanism removes $N=C_4$ with $K=D_8$ at $d=32$. It also removes two of the four strata $K\supsetneq N=C_6$: those with $g(C/K)=7$ and $3$.
- **The record stratum.** For $K=C(\tau)$, $N=\langle\tau\rangle$, $\bar K=C_2\times S_3\cong D_6$. $D\to T=C/C(\tau)$ (genus 3) has 13 inertia points:
  - nine of type $R_a$, the reflections $[(12)(56)]$;
  - one of type $R_b$, the reflections $[(1324)(56)]$;
  - three of type $z$, the central half-turn $[(13)(24)]$.

  Over the $C_6$-poles, the three $z$-points need $e\equiv3\pmod6$, and every other point needs $e\equiv0\pmod6$. So $k\equiv3\pmod6$, i.e. $d\equiv6\pmod{12}$. Below 42 only $d=30$ survives. The $\tau$-pencil itself has $k=21$.

So the kernel list of Corollary 7.23(2) is unchanged, but many of its degrees are removed. The output lists every stratum.

### The record stratum

Fix $\tau=(12)(34)$, $z=(13)(24)$, $s=(12)(56)$, and let
- $U=\langle\tau,z,(567)\rangle\cong V_4\times C_3$, the preimage of the rotation group $C_6\subset D_6$;
- $T=C/C(\tau)$, of genus 3;
- $X_{\rm rot}=C/U$, of genus 10, with involution $\iota$ (induced by $s$) and $X_{\rm rot}/\iota=T$;
- $\theta\in\mathrm{Pic}^0(X_{\rm rot})[3]$, the class of the étale cyclic cover $C/\langle\tau,z\rangle\to X_{\rm rot}$ (Galois group $U/\langle\tau,z\rangle\cong C_3$).

**Proposition 7.25 (dihedral norms) [P].** Let $f$ be a pencil of degree $d$ with $(K,N)=(C(\tau),\langle\tau\rangle)$, with coordinates in which $D_6$ is generated by $x\mapsto\zeta_6x$ and $x\mapsto1/x$.
1. $w=f^6$ descends to $X_{\rm rot}$, with $\iota^*w=w^{-1}$ and $\operatorname{div}w=3(\Delta-\iota\Delta)$.
2. Here $\Delta\ge0$ has degree $d/6$, and $\Delta=P^0+2E$, where $P^0$ contains one preimage of each $z$-point.
3. $(1-\iota)\Delta\sim\pm\theta$.
4. On $T$ the pencil is $y=\tfrac14(w+2+w^{-1})=\operatorname{Tr}(g)^2/4\operatorname{Nm}(g)$ for any $g$ with $w=g/\iota g$, and
   $$y^*(\infty)=\textstyle\sum_je_jP_j+6(\cdots),\qquad y^*(0)=A+2E_a,\qquad y^*(1)=B+2E_b,$$
   where $A$ is the nine $R_a$-points and $B$ the $R_b$-point.

*Proof.*
- **The descent.** Rotations fix $x^6$ and reflections invert it, so $w$ descends with $\iota^*w=w^{-1}$.
- **Orders of $w$.** Over a point of $D$ with $f=0$, $C\to X_{\rm rot}$ has ramification 2 at the $z$-points and 1 elsewhere. So $\operatorname{ord}w=3\operatorname{ord}f$ (odd) at the $z$-points and $6\operatorname{ord}f$ elsewhere; similarly at the poles.
- **The torsion class.** $w$ is a cube on $C/\langle\tau,z\rangle$, namely $f^2$. This means exactly that $(1-\iota)\Delta$ is the class of that cyclic cover, up to sign.
- **The pencil on $T$.** The last formula is the identity $\tfrac14(w+2+w^{-1})=(g+\iota g)^2/4g\,\iota g$. The fibres are Lemma 7.24, with the labelling $b_a=0$, $b_b=1$, $b_3=\infty$. $\square$

So $42=6\cdot7$: **the $\tau$-pencil is an effective $\Delta$ of degree 7 with $(1-\iota)\Delta$ equal to a 3-torsion point** of the Prym $P=\ker(\mathrm{Nm})^0$, which has dimension $10-3=7$ and $P\sim A\times E_{15}\times E_{21}\times S$. On $T$ the naive parameter count for $k=21$ is 40 against 45 conditions, so the $\tau$-pencil is not expected from dimensions. Its existence comes from the twisted class $L_{60}$.

**Lemma 7.26 (Frobenius rigidity) [P].** Let $X\to X/\iota$ be a double cover of curves over $\mathbb F_q$ such that every function of degree $\le2m$ on $X$ is pulled back from $X/\iota$. Let $P^0\in\mathrm{Div}(X)$ and $\theta\in\mathrm{Pic}(X)$ be $\mathbb F_q$-rational. If $Q\in X(\overline{\mathbb F}_q)$ satisfies $(1-\iota)(P^0+mQ)\sim\theta$, then $Q\in X(\mathbb F_q)$ or $\iota Q=Q$.

*Proof.* Applying Frobenius $\phi$ and subtracting gives a function $F$ with $\operatorname{div}F=m(\phi Q+\iota Q)-m(Q+\iota\phi Q)$, of degree $\le2m$. By hypothesis $F$ is $\iota$-invariant. But $\iota$ negates this divisor, so $\operatorname{div}F=0$, and $\phi Q+\iota Q=Q+\iota\phi Q$. $\square$

For $X_{\rm rot}\to T$ and $m=2$ the hypothesis is Castelnuovo–Severi: a function of degree $n\le4$ not from $T$ would give $10\le2\cdot3+(n-1)\le9$. It holds in every characteristic.

**Theorem 7.27 (the $\tau$-pencil is optimal in its stratum) [P][X].** On all four classes, every pencil with class stabiliser $C(\tau)$ and kernel $\langle\tau\rangle$ has degree $\ge42$. The $\tau$-pencil attains 42.

*Proof.*
1. **Reduction to $d=30$.** By Lemma 7.24, $d=30$, $k=15$, and $y^*(\infty)$ is either $3P_1+3P_2+3P_3+6Q$ (free type) or $3P_i+3P_j+9P_l$ (fixed type).
2. **Reduction mod 43.** Reduce at a prime above 43 of a field of definition, with good reduction. All conditions are closed:
   - $y\in L(F)$;
   - $\operatorname{div}(y)+F-A$ and $\operatorname{div}(y-\lambda)+F-B$ effective and even.

   So a solution in characteristic 0 gives one over $\overline{\mathbb F}_{43}$. Here $3\mid42$, so $\phi\theta=\theta$.
3. **The explicit models [X]** (`tau_stratum.py`, from the septic of §6.5).
   - $W=C/(S_4\times S_3)^+$ (3-sets of roots, genus 0) is parametrised by Riemann–Roch.
   - $T$ is the cubic resolvent of the complementary quartic over $W$, and $X_{\rm rot}=T(\sqrt{\operatorname{disc}A})$.
   - The 13 inertia points are identified from the specialisations of the 3-set and the pairing. The odd part of $\operatorname{div}(\operatorname{disc}A)$ is exactly $A+B$, which independently confirms the identification and gives $g(X_{\rm rot})=10$.
   - The preimages of $P_1,P_2,P_3$ in $X_{\rm rot}$ are $\mathbb F_{43}$-rational.
4. **The test [X].** For given $F,A',B'$, $y$ ranges over $L(F-A')\cap(L(F-B')+\mathbb F_p)$. The norms of $y$ and $y-\lambda$ to the $w$-line of $W$, with their fixed parts removed, must be squares of binary forms of degrees $\tfrac12(\deg F-\deg A')$ and $\tfrac12(\deg F-\deg B')$. Gröbner bases decide this over $\overline{\mathbb F}_p$. An empty ideal excludes the configuration.

   Positive control: a planted dihedral norm $y=a^2/(a^2-b^2(w-c))$ is found, and the same $F$ with a wrong $B'$ is excluded.
5. **Fixed types.** All three are excluded.
6. **Free type.** By Lemma 7.26, the reduction of the point of $\Delta$ above $Q$ is $\mathbb F_{43}$-rational or a ramification point.
   - **Rational:** the 53 rational points of $T$ outside $P\cup A\cup B$ are tested one by one. Those inside $P$ are the fixed types. All 53 are excluded (`tau_stratum_output.txt`, 87 min).
   - **Ramified:** then $(1-\iota)P^0\sim\pm\theta$. This gives a function $\tfrac14(w_0+2+w_0^{-1})$ on $T$ with poles $3(P_1+P_2+P_3)$, sending the ten ramification points to $0$ and $1$ in a split $A'\sqcup B'$ of odd sizes. A Frobenius-swapped split is impossible, because $A\cup B$ has rational points. All 16 Galois-stable splits are excluded.
7. **The other curve.** The curves of classes $\{0,1\}$ and $\{12,14\}$ are the Galois closures of the two septics, which are conjugate under $\sqrt{21}\mapsto-\sqrt{21}$. So the curves are $\mathrm{Aut}(\mathbb C)$-conjugate, and existence of a pencil of given degree and symmetry type is invariant. $\square$

**Remark (Galois conjugacy).** Step 7 applies to every algebraic statement. Gonality, the strata of Corollary 7.23 and the elliptic subcovers are the same on both curves; only spectral data differ (§2.5). In particular, a proof of $\operatorname{gon}\ge26$ on classes 0, 1 would give it on classes 12, 14, where the harmonic Hersch route is closed.

**Consequences.**
- **Why 42 is hard to beat by symmetry.** A pencil below 42 cannot have the full symmetry of the record pencil.
- **What remains in the $C_2$ stratum.** Kernel $\langle\tau\rangle$ with $K\subsetneq C(\tau)$ has $K/N\in\{1,C_2,V_4,C_3,S_3,C_6\}$, with the degrees of `orbifold_congruence_output.txt`.
- **The method transfers.** The other strata with $g(C/K)=3$ ($K\cong S_4$, three classes in $(S_4\times S_3)^+$) have quotients that are trigonal over $W$ and can be treated the same way.
