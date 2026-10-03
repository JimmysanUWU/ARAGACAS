# Chapter 1. The $(2,4,7)$ curve and its involution quotients

**Idea.** Everything starts from fixed points.
- An involution $\tau$ fixes 18 points of $C$. Riemann–Hurwitz turns this into the genus of $D=C/\tau$ and of every quotient of $D$.
- Castelnuovo–Severi then compares a low-degree pencil on $D$ with those quotients. Below degree 9 the pencil is forced to be invariant under all of $H=C_G(\tau)/\langle\tau\rangle$, which is impossible.
- At degree 9 the comparison becomes an equality. The pencil's $H$-orbit must then take one of three shapes, and each forces the same linear equivalence $(\star)$ on the genus-3 curve $D/H$.
- Whether $(\star)$ holds is an arithmetic question. An explicit "audit curve" shows that the $H$-action on $D$ alone cannot decide it, so going further needs the metric (Chapter 2) or the whole of $A_7$ (§6.1).

**Setting.** $C$ is a smooth connected curve with a faithful $G=A_7$-action and $C/G\cong\mathbb P^1$, branched at three points with inertia orders $2,4,7$. Equivalently, $C$ is given by a generating triple $(a,b,c)$ of $A_7$ with orders $2,4,7$ and $abc=1$. Throughout, $\tau$ is an involution and $D=C/\langle\tau\rangle$.

## 1.1 The group and the curves

- **Triples [X]** (`curve_checks.py` §1). With $a=(12)(34)$ fixed there are 96 choices of $b$, and 10080 triples in all.
  - They form 4 classes under $A_7$-conjugation, labelled 0, 1, 12, 14 (indices into `triples_data.py`).
  - They form 2 classes under $S_7$: $\{0,1\}$ and $\{12,14\}$. Within an $S_7$-orbit the curves are isomorphic, and the actions differ by an outer automorphism.
- **Conjugacy data [P].**

| elements | cycle type | count | centraliser |
|---|---|---|---|
| involutions | $2^21^3$ | 105 | 24 |
| order 4 | $(4,2,1)$ | 630 | 4 |
| order 7 | two classes | 360 each | 7 |

- **Genus [P].** $2g-2=2520\big(1-\frac12-\frac14-\frac17\big)=270$, so $g(C)=136$. This is the minimum genus of a faithful $A_7$-curve (exhaustive search, `verify_accessory60.py`).

## 1.2 Fixed points and the quotient $D$

**Lemma 1.1 [P].** Over a branch point with stabiliser $\langle x\rangle$, an element $h$ fixes $|C_G(h)|\cdot|h^G\cap\langle x\rangle|/|x|$ points. (The points are the cosets $g\langle x\rangle$, and $h$ fixes $g\langle x\rangle$ iff $g^{-1}hg\in\langle x\rangle$. There are $|C_G(h)|\cdot|h^G\cap\langle x\rangle|$ such $g$.)

**Corollary 1.2 [P][X]** (`curve_checks.py` §2). Fixed points by order of the element:

| order | 2 | 4 | 7 | 3, 5, 6 |
|---|---|---|---|---|
| fixed points | $24\cdot\frac12+24\cdot\frac14=18$ | 2 | 3 | 0 |

Hence $270=2(2g(D)-2)+18$, and $g(D)=64$. Point stabilisers are cyclic, so distinct commuting involutions have disjoint fixed sets.

**The action on $D$.** $H=C_G(\tau)/\langle\tau\rangle\cong C_2\times S_3$ acts on $D$. For $\tau=(12)(34)$ it is $\{(x,y)\in V_4\times S_3:\operatorname{sgn}x=\operatorname{sgn}y\}$, where $V_4$ is the quotient of the dihedral group on $\{1,2,3,4\}$ by $\tau$.

**Proposition 1.3 [P][X]** (`curve_checks.py` §3). Use $\#\mathrm{Fix}_{C/V}(\bar w)=\frac1{|V|}\sum_{g\in wV}\#\mathrm{Fix}_C(g)$.
- **Good involutions of $H$.** There are 4: the images of the involutions $\mu\ne\tau$ commuting with $\tau$. Each has 18 fixed points on $D$, and $g(D/\bar\mu)=28$. One of them, the image of $(13)(24)$, is central.
- **Bad involutions.** There are 3: the images of the order-4 elements with $\mu^2=\tau$. Each has 2 fixed points on $D$, and $g(D/\bar\mu)=32$.
- **Orders 3 and 6** act freely.
- **Quotient genera.** Over the 16 subgroups of $H$, $g(D/K)$ takes the values $64,28,32,22,12,10,11,7,3$. In particular $g(D/C_3)=22$, $g(D/S_3)=7$ for the $S_3$ containing the non-central good involutions, and $g(D/H)=3$.

## 1.3 The first algebraic barrier: $\operatorname{gon}(D)\ge9$

**Castelnuovo–Severi.** If $f_i:X\to Y_i$ have degrees $d_i$ and do not both factor through a common nontrivial map, then $g(X)\le d_1g(Y_1)+d_2g(Y_2)+(d_1-1)(d_2-1)$.

**Theorem 1.4 [P].** $\operatorname{gon}(D)\ge9$.

*Proof.*
1. Take $f:D\to\mathbb P^1$ of degree $d\le8$ and a good involution $v$.
2. If $f$ did not factor through $D\to D/v$ (degree 2, genus 28), Castelnuovo–Severi would give $64\le56+(d-1)$. So $f\circ v=f$.
3. The four good involutions generate $H$. So $f$ factors through $D/H$, and $12\mid d$. Contradiction. $\square$

## 1.4 Degree 9: what $D$ alone cannot decide

At $d=9$, Castelnuovo–Severi against each good involution is an equality ($64=56+8$), so it no longer forces $f\circ v=f$. Moreover $f\circ v=f$ is now impossible, because 9 is odd. Everything therefore happens at the level of pencils.

Write $f_h=f\circ h$, $[f_h]$ for the pencil $\mathbb C(f_h)\subseteq\mathbb C(D)$, and $F=\mathbb C(f_h:h\in H)$. Let $T=D/H$, of genus 3.

**Theorem 1.5 (structure) [P].** Suppose $\operatorname{gon}(D)=9$, realised by $f$.
1. **Three shapes.** The index $[\mathbb C(D):F]$ is 9, 3 or 1, and accordingly:
   - **(a)** the pencil is $H$-invariant, $H$ acts faithfully on its target as the dihedral group $D_6$, and $D$ is the normalised pullback of $\mathbb P^1\to\mathbb P^1/D_6$ along a degree-9 map $\varphi:T\to\mathbb P^1$;
   - **(b)** $F=\mathbb C(Y)$ with $g(Y)=4$, $H$ acts faithfully on $Y$ with signature $(2,2,2,2,2)$, and $D$ is the normalised pullback of $Y\to Y/H\cong\mathbb P^1$ along a degree-3 map $\varphi:T\to\mathbb P^1$;
   - **(c)** $D$ is a smooth $(9,9)$-curve in $\mathbb P^1\times\mathbb P^1$ whose rulings are pencils $[f_h]$, and the good involutions exchange the rulings.
2. **Branch data.** $D\to T$ has 13 branch points, all with monodromy of order 2:
   - $t_1,t_2,t_3$ from the central good involution;
   - $g_1,\dots,g_9$ from the three other good involutions;
   - one point $t_b$ from the bad involutions.
3. **The relation.** In every case
$$(\star)\qquad g_1+\dots+g_9\ \sim\ 3(t_1+t_2+t_3)\quad\text{on }T .$$

*Proof.*
1. **Trichotomy.** For $[f_h]\ne[f_{h'}]$, the field $\mathbb C(f_h,f_{h'})$ has index dividing 9 and less than 9, so 1 or 3.
   - *Index 1.* Castelnuovo–Severi for two maps of degree 9 gives $64\le8^2$, an equality. So $(f_h,f_{h'})$ embeds $D$ as a smooth $(9,9)$-curve, on which the only $g^1_9$'s are the rulings (Martens). This is (c).
   - *Index 3.* Then $\mathbb C(f_h,f_{h'})=\mathbb C(Y)$, where $Y$ carries two degree-3 maps generating its field, so $g(Y)\le4$. For any $k$, $\mathbb C(Y,f_k)$ has index 1 or 3. Index 1 is excluded by Castelnuovo–Severi ($64\le3\cdot4+2\cdot8=28$), so every $f_k\in\mathbb C(Y)$ and $F=\mathbb C(Y)$. This is (b).
   - If all the $[f_h]$ coincide, $F=\mathbb C(f)$. This is (a).
2. **Case (a).** Write $f\circ h=\rho(h)\circ f$ with $\rho:H\to\mathrm{PGL}_2$. Then $|\ker\rho|$ divides 9, so $\ker\rho$ is 1 or $C_3$.
   - Suppose it is $C_3$. Then $f=f_3\circ q$ with $q:D\to D/C_3$ (genus 22) and $\deg f_3=3$. A non-central good involution $w$ acts on $D/C_3$ with quotient $D/S_3$ of genus 7. Castelnuovo–Severi gives $22\le14+2$ unless $f_3$ factors through that quotient, which is impossible since 3 is odd.
   - So $\rho$ is faithful and $\rho(H)\cong D_6$. Since $H$ acts faithfully on $\mathbb C(f)$, Galois theory gives $\mathbb C(f)\cdot\mathbb C(T)=\mathbb C(D)$. That is the pullback statement, with $\deg\varphi=[\mathbb C(D):\mathbb C(f)]=9$.
3. **Case (b).** $H$ acts faithfully on $Y$, since a kernel $C_3$ would make $Y=D/C_3$, of genus $22>4$. The same Galois argument gives the pullback, with $\deg\varphi=3$. The signature is forced in step 6.
4. **Case (c).** A good involution preserving both rulings acts as $(\alpha,\beta)$.
   - If one factor is the identity, it fixes a projection, which is impossible for degree 9.
   - Otherwise it has at most 4 fixed points, not 18.

   So the good involutions exchange the rulings.
5. **Branch data.** The points of $D$ with nontrivial stabiliser in $H$ are:
   - the 18 fixed points of the central good involution, in orbits of size 6, giving $t_1,t_2,t_3$;
   - $3\cdot18$ points of the others, giving $g_1,\dots,g_9$;
   - $3\cdot2$ points of the bad involutions, giving $t_b$.

   Check: $126=12\cdot4+6\cdot13$.
6. **The relation $(\star)$.** In a pullback the monodromy over $t$ is the monodromy over $\varphi(t)$ raised to the power $e_\varphi(t)$, and $D\to T$ has no monodromy of order 3 or 6. In each case the fibres of $\varphi$ are linearly equivalent.
   - *(a)* $\mathbb P^1\to\mathbb P^1/D_6$ is branched at three points: $y_6$ (rotation $s$, with $s^3$ the central involution), $y_1$ (good reflections) and $y_2$ (bad reflections).
     - Over $y_6$ every ramification index is $\equiv0\pmod3$, and the monodromy is central exactly when $e\equiv3\pmod6$. Three such points in a fibre of size 9 force $\varphi^\*(y_6)=3(t_1+t_2+t_3)$.
     - Good monodromy needs odd $e$ over $y_1$, and nine such points force $\varphi^\*(y_1)=\sum g_j$.
   - *(b)* Each good branch point of $Y\to Y/H$ yields at most 3 of the $g_j$, so there are at least three of them. There is also at least one central and one bad branch point.
     - Riemann–Hurwitz gives $2g(Y)-2\ge12(2h-2)+5\cdot6$. With $g(Y)\le4$ this forces $h=0$, $g(Y)=4$ and exactly these five branch points.
     - Then $\varphi^\*(y_{\rm central})=t_1+t_2+t_3$, and each good fibre is three of the $g_j$.
   - *(c)* Let $\Delta$ and $\Gamma_1,\Gamma_2,\Gamma_3$ be the fixed curves of the central and the non-central good involutions. They are the diagonal and three graphs, all $(1,1)$-curves.
     - $\Psi=G_1G_2G_3/\Delta^3$ has $\operatorname{div}\Psi=\sum\Gamma_i-3\Delta$ and $\Psi\circ h=\pm\Psi$, so $\Psi^2$ descends to a function $\psi$ on $T$.
     - $D$ meets each of these curves transversally, exactly in the 18 fixed points ($18=(9,9)\cdot(1,1)$). Hence $\pi^\*\operatorname{div}\psi=\operatorname{div}(\Psi^2|_D)=\pi^\*\big(\sum g_j-3\sum t_i\big)$. $\square$

**The Jacobian [X]** (`curve_checks.py` §4, Chevalley–Weil). $H^1(C)=3(10+\overline{10})+2\cdot15+2\cdot21+4\cdot35$, so
$$\mathrm{Jac}(C)\sim A^{10}E_{15}^{15}E_{21}^{21}S^{35},\qquad\mathrm{Jac}(D)\sim A^4E_{15}^7E_{21}^{11}S^{17},\qquad\mathrm{Jac}(T)\sim E_{21}\times S .$$
Here:
- $A$ is a threefold with $\mathbb Q(\sqrt{-7})$-multiplication;
- $E_{15}=C/A_5$ ($A_5$ fixing two points) and $E_{21}=C/L_2(5)$ are elliptic curves, named by the isotype they carry (§7.9);
- $S\sim\mathrm{Jac}(C/(3^2{:}4))$ is an abelian surface.

The irreducibles $1,6,14_a,14_b$ do not occur in $H^1(C)$.

**Proposition 1.6 [P][X]** (`curve_checks.py` §5, `side_checks.py` D4).
- The $\mathbb Q[G]$-span of the involution fixed divisors is a quotient of $\mathbb Q[G/C(\tau)]=1+6+2\cdot14_a+14_b+21+35$. Components in $1+6+14_a+14_b$ are torsion in $\mathrm{Jac}(C)$.
- The lift of $(\star)$ to $C$ has nonzero 21- and 35-components.
- So $(\star)$ fails as soon as either of two explicit points has infinite order: $P_E\in E_{21}$ or $P_S\in\mathrm{Jac}(C/(3^2{:}4))$ (§6.5).

**Theorem 1.7 (audit curve) [P][X]** (`model99.py`). A smooth $(9,9)$-curve $D'\subset\mathbb P^1\times\mathbb P^1$ of genus 64 carries a faithful $C_2\times S_3$-action with exactly the fixed-point data of $D$, and $\operatorname{gon}(D')=9$.

*Proof.*
1. **The curve.** $D'=\{F=0\}$ has integer coefficients $c_{ij}$ supported on $i+j\equiv0\pmod3$, symmetric under $i\leftrightarrow j$ and $(i,j)\mapsto(9-i,9-j)$. The group $S_3=\langle z\mapsto\omega z,\ z\mapsto1/z\rangle$ acts diagonally, and the swap is the central involution.
2. **Checks modulo $p=1000003$**, which are open conditions and so hold over $\mathbb Q$:
   - $D'$ is smooth in the chart $x,y\ne\infty$. The involution $(x,y)\mapsto(1/x,1/y)$ moves every other point of $D'$ into that chart, because the corner coefficients $c_{09},c_{90}\ne0$ keep $(0,\infty)$ and $(\infty,0)$ off $D'$.
   - $F(x,x)$ and $F(x,1/x)$ are squarefree of degree 18. So the swap and the good involutions have exactly 18 fixed points.
   - $F(1,1),F(-1,-1)\ne0=F(1,-1)=F(-1,1)$, so each bad involution fixes exactly 2 points.
   - The corner coefficients are nonzero, so elements of order 3 and 6 act freely.
3. **Gonality.** A ruling gives $\operatorname{gon}\le9$. The proof of Theorem 1.4 uses only this fixed-point data, so $\operatorname{gon}\ge9$. $\square$

So the fixed points, the quotient genera and the $H$-character of $H^0(K)$ (which holomorphic Lefschetz computes from the fixed points) cannot decide degree 9. One needs more of $C$:
- the hyperbolic metric (Chapter 2), or
- the full $A_7$-action (§6.1).

## 1.5 The ladder

| rung | statement | where |
|---|---|---|
| 1 | $g(C)=136$; 18 fixed points; $g(D)=64$ | §1.2 |
| 2 | $\operatorname{gon}(D)\ge9$; structure at 9; reduction to $(\star)$ | §§1.3–1.4 |
| 3 | $\operatorname{gon}(D)\ge10$ by algebra | §6.1 |
| 4 | $\operatorname{gon}(C)\ge25$, $\operatorname{gon}(D)\ge13$ | Chapters 2–3 |
| upper | $\operatorname{gon}(C)\le42$, $\operatorname{gon}(D)\le21$, from a degree-60 twisted class | Cor. 7.6, Prop. 7.7 |
