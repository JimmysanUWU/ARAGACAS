# Chapter 1. The $(2,4,7)$ curve and its involution quotients

**Setting.** $C$ is a smooth connected curve with a faithful $G=A_7$-action and $C/G\cong\mathbb P^1$ branched at three points with inertia orders $2,4,7$. Equivalently it is given by a generating triple $(a,b,c)$ of $A_7$ with orders $2,4,7$ and $abc=1$. Throughout, $\tau$ is an involution and $D=C/\langle\tau\rangle$.

## 1.1 The group and the curves

- **Triples [X]** (`curve_checks.py` §1). With $a=(12)(34)$ fixed there are 96 choices of $b$, and 10080 triples in all.
  - These form 4 classes under $A_7$-conjugation, labelled 0, 1, 12, 14 (`triples_data.py`).
  - They form 2 classes under $S_7$: $\{0,1\}$ and $\{12,14\}$. Within an $S_7$-orbit the curves are isomorphic, and the actions differ by an outer automorphism.
- **Conjugacy data [P].**

| elements | cycle type | count | centraliser |
|---|---|---|---|
| involutions | $2^21^3$ | 105 | 24 |
| order 4 | $(4,2,1)$ | 630 | 4 |
| order 7 | two classes | 360 each | 7 |

- **Genus [P].** $2g-2=2520(1-\frac12-\frac14-\frac17)=270$, so $g(C)=136$. This is the minimum genus of a faithful $A_7$-curve (exhaustive search, `verify_accessory60.py`).

## 1.2 Fixed points and the quotient $D$

**Lemma 1.1 [P].** Over a branch point with stabiliser $\langle x\rangle$, an element $h$ fixes $|C_G(h)|\cdot|h^G\cap\langle x\rangle|/|x|$ points.

**Corollary 1.2 [P][X]** (`curve_checks.py` §2). Fixed points by order of the element:

| order | 2 | 4 | 7 | 3, 5, 6 |
|---|---|---|---|---|
| fixed points | $24\cdot\frac12+24\cdot\frac14=18$ | 2 | 3 | 0 |

Hence $270=2(2g(D)-2)+18$ and $g(D)=64$.

**The action on $D$.** $H=C_G(\tau)/\langle\tau\rangle\cong C_2\times S_3$ acts on $D$. For $\tau=(12)(34)$ it is $\{(x,y)\in V_4\times S_3:\operatorname{sgn}x=\operatorname{sgn}y\}$.

**Proposition 1.3 [P][X]** (`curve_checks.py` §3). Use $\#\mathrm{Fix}_{C/V}(\bar w)=\frac1{|V|}\sum_{g\in wV}\#\mathrm{Fix}_C(g)$.
- **Good involutions of $H$.** There are 4, the images of involutions $\mu\ne\tau$ commuting with $\tau$. Each has 18 fixed points on $D$, and $g(D/\bar\mu)=28$.
- **Bad involutions.** There are 3, the images of order-4 elements with $\mu^2=\tau$. Each has 2 fixed points, and $g(D/\bar\mu)=32$.
- **Orders 3 and 6** act freely.
- **Quotient genera.** Over the 16 subgroups of $H$, $g(D/K)$ ranges over $64,28,32,22,12,10,11,7,3$. In particular $g(D/H)=3$.

## 1.3 The first algebraic barrier: $\operatorname{gon}(D)\ge9$

**Castelnuovo–Severi.** If $f_i:X\to Y_i$ have degrees $d_i$ and no common factorisation, then $g(X)\le d_1g(Y_1)+d_2g(Y_2)+(d_1-1)(d_2-1)$.

**Theorem 1.4 [P].** $\operatorname{gon}(D)\ge9$.

*Proof.*
1. Take $f$ of degree $d\le8$ and a good involution $v$.
2. Castelnuovo–Severi against $D\to D/v$ would give $64\le56+(d-1)$. So $f\circ v=f$.
3. The four good involutions generate $H$. So $f$ factors through $D/H$, and $12\mid d$. Contradiction. $\square$

## 1.4 Degree 9: what $D$ alone cannot decide

At $d=9$, Castelnuovo–Severi is an equality ($64=56+8$) against each good involution. Then $f\circ v=M_v\circ f$ with $M_v$ Möbius.

**Theorem 1.5 (structure) [P].** Suppose $\operatorname{gon}(D)=9$, realised by $f$. Write $T=D/H$ (genus 3). The cover $D\to T$ has 13 branch points: $t_1,t_2,t_3$ from the central good involution, $g_1,\dots,g_9$ from the others, and one from the bad ones.

The translates $f\circ h$ generate a field of index 9, 3 or 1, and accordingly:
- **(a)** $D$ is the pullback of $\mathbb P^1\to\mathbb P^1/D_6$ along a degree-9 map $T\to\mathbb P^1$;
- **(b)** $D$ is the pullback of $Y\to Y/H$, with $g(Y)\le4$, along a degree-3 map;
- **(c)** $D$ is a smooth $(9,9)$-curve whose rulings are translates of $f$ (Martens 1996).

In each case the monodromy forces the linear equivalence
$$(\star)\qquad g_1+\dots+g_9\ \sim\ 3(t_1+t_2+t_3)\quad\text{on }T .$$

**The Jacobian [X]** (`curve_checks.py` §4, Chevalley–Weil). $H^1(C)=3(10+\overline{10})+2\cdot15+2\cdot21+4\cdot35$, so
$$\mathrm{Jac}(C)\sim A^{10}E_1^{15}E_2^{21}S^{35},\qquad\mathrm{Jac}(D)\sim A^4E_1^7E_2^{11}S^{17},\qquad\mathrm{Jac}(T)\sim E_2\times S .$$
Here:
- $A$ is a threefold with $\mathbb Q(\sqrt{-7})$-multiplication;
- $E_1=C/A_5$ ($A_5$ fixing two points) and $E_2=C/L_2(5)$ are elliptic curves (§7.9);
- $S\sim\mathrm{Jac}(C/(3^2{:}4))$ is an abelian surface.

**Proposition 1.6 [P][X]** (`curve_checks.py` §5, `side_checks.py` D4).
- The $\mathbb Q[G]$-span of the involution fixed divisors is a quotient of $\mathbb Q[G/C(\tau)]=1+6+14_a+2\cdot14_b+21+35$.
- The lift of $(\star)$ has nonzero 21- and 35-components.
- So $(\star)$ fails as soon as either of two explicit points has infinite order: $P_E\in E_2$ or $P_S\in\mathrm{Jac}(C/(3^2{:}4))$ (§6.5).

**Theorem 1.7 (audit curve) [P][X]** (`model99.py`). A smooth $(9,9)$-curve $D'\subset\mathbb P^1\times\mathbb P^1$ of genus 64 carries a faithful $C_2\times S_3$-action with exactly the fixed-point data of $D$, and has $\operatorname{gon}(D')=9$.

So the fixed points, the quotient genera and the $H$-character of $H^0(K)$ cannot decide degree 9. One needs more of $C$:
- the hyperbolic metric (Chapter 2), or
- the full $A_7$-action (§6.1).

## 1.5 The ladder

| rung | statement | where |
|---|---|---|
| 1 | $g(C)=136$; 18 fixed points; $g(D)=64$ | §1.2 |
| 2 | $\operatorname{gon}(D)\ge9$; structure at 9; reduction to $(\star)$ | §§1.3–1.4 |
| 3 | $\operatorname{gon}(D)\ge10$ by algebra | §6.1 |
| 4 | $\operatorname{gon}(C)\ge25$, $\operatorname{gon}(D)\ge13$ | Chapters 2–3 |
| upper | $\operatorname{gon}(C)\le42$, $\operatorname{gon}(D)\le21$: an exact base-point-free pencil on a degree-60 model in $\mathbb P^5$ | Prop. 7.7 |

Before Chapter 7 the best upper bound was the quotient $C\to C/P\to\mathbb P^1$, with $P\ni\tau$ of order 8 and $g(C/P)=12$. It gives only $56$ and $28$.
