# GUIDE

How to read this repository, for a reader (human or AI) starting cold. Live status and open problems: [`HANDOFF.md`](HANDOFF.md).

## 1. The problem and the results

Hilbert's 13th problem, in its algebraic form, asks whether the general septic needs algebraic functions of three variables, i.e. whether $\mathrm{RD}(7)=3$; it is open. It is approached here through two invariants of $A_7$:
- **gonality**: the least degree of a map to $\mathbb P^1$ from a curve with a faithful $A_7$-action;
- **accessory degree**: the least degree of a finite extension after which the generic $A_7$-torsor descends to a curve.

Let $C$ be a $(2,4,7)$ $A_7$-curve (genus 136, the minimum) and $\tau$ an involution.
1. $25\le\operatorname{gon}(C)\le42$ and $13\le\operatorname{gon}(C/\tau)\le21$.
2. $25\le\gamma(A_7)\le42$, where $\gamma$ is the least gonality of a faithful $A_7$-curve.
3. $a(A_7)=60$, i.e. $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$, sharp (the published bound was 6).
4. $\mu(A_7)=90$ when the full monodromy must stay connected. For a fixed $(2,4,7)$ target and a smooth projective generically free base with $\mathrm{Hom}_{A_7}(\mathrm{Alb}\,X,\mathrm{Jac}\,C)=0$, the stable threshold is 60 or 90 according to the base's Amitsur subgroup.

## 2. The mathematics in one page

Each chapter opens with a paragraph headed **Idea**; read it first.

**Fixed points and Castelnuovo–Severi** (Ch. 1, §6.1).
- An involution fixes 18 points, which fixes the genus of every quotient of $D=C/\tau$.
- Castelnuovo–Severi forces a pencil on $D$ of degree $\le8$ to be invariant under a group that is too large, so $\operatorname{gon}(D)\ge9$.
- At degree 9 the pencil's orbit has three possible shapes, each forcing a linear equivalence on a genus-3 curve.
- Ramification transport (§6.1) uses all of $A_7$ to exclude degree 9.
- *Limit:* the $H$-action on $D$ alone cannot decide degree 9 (the audit curve).

**Degree versus energy** (Ch. 2–3).
- A balanced map of degree $d$ has energy $8\pi d$, so $\lambda_1$ bounds $d$ (Li–Yau). Certifying $\lambda_1$ gives 24.
- The degree is also a cubic form, which vanishes on the first eigenspace because $\wedge^314_{(5,2)}$ has no invariants. Harmonic Hersch makes this quantitative: 25.
- *Limit:* the spectral route stops at 25 or 26.

**Schur obstructions** (Ch. 5, 7).
- An invariant line bundle need not carry a $G$-action. The obstruction lies in $H^2(A_7,\mathbb C^\*)=\mathbb Z/6$. Invariant classes have degrees $15\mathbb Z$, linearised ones $90\mathbb Z$.
- **Upper bound 42.** The first twisted class with sections has degree 60. Its six sections map $C$ birationally into $\mathbb P^5$, and an involution's $(-1)$-eigenspace is a pencil of degree $\le42$.
- **Embedding and exact pencil.** The arithmetic-genus defect is at most $\pi(60,5)-136=270$, too small for a nonimmersed orbit. Exact order-7 eigenline stabilizers exclude the remaining possible collisions. So the map is an embedding and the displayed pencil has degree exactly 42 (§7.12).
- **The two quadratic systems.** They are the cubic's Jacobian system and its dual apolar system. Both have empty projective base locus, so the curve lies on no quadric; the Klein quadrics' fivefold contact divisors are exact (§7.12).
- **Accessories.** Over a linear base a compression is one linearised moving series, so accessory degrees are degrees of linearised series induced from subgroups: $a=15\cdot4=60$ (via the Klein quartic of $L_2(7)$), and $\mu=90$ with full monodromy. A general base cancels exactly its Amitsur subgroup's worth of obstruction.
- *Limit:* in towers the Schur–Brauer obstruction is cheap to remove; what remains is the Albanese part.

**Pencil orbits and symmetry** (Ch. 4, §7.10).
- A pencil of degree $\le24$ has at least 15 conjugates. Pairs give $\mathbb P^1\times\mathbb P^1$ models and triples $\mathbb P^7$ models, which bounds the genus of an $A_7$-curve of gonality $\le24$ by 397.
- A gonal pencil's class stabiliser acts on it through $PGL_2$. The kernel $N$ costs degree $|N|\operatorname{gon}(C/N)$, and quotient gonalities put this cost at $\ge42$ unless $N$ is trivial or one of six small groups.

## 3. Reading order

| step | file | what you get |
|---|---|---|
| 1 | this guide | the mathematics in one page, notation, dependencies |
| 2 | [`hilbert13/ladder/README.md`](hilbert13/ladder/README.md) | every result with its tag and chapter; scripts |
| 3 | [`1_CURVE.md`](hilbert13/ladder/1_CURVE.md) | the group, the curves, $D=C/\tau$; the first bounds; degree 9 |
| 4 | [`2_SPECTRAL.md`](hilbert13/ladder/2_SPECTRAL.md), [`3_HARMONIC_HERSCH.md`](hilbert13/ladder/3_HARMONIC_HERSCH.md) | the lower bound 25 |
| 5 | [`7_TWISTED.md`](hilbert13/ladder/7_TWISTED.md) | the upper bound 42; the curve in $\mathbb P^5$; pencils by symmetry type |
| 6 | [`4_LARGE_GENUS.md`](hilbert13/ladder/4_LARGE_GENUS.md) | from one curve to all faithful $A_7$-curves |
| 7 | [`5_ACCESSORY.md`](hilbert13/ladder/5_ACCESSORY.md) | $a=60$, $\mu=90$, Amitsur subgroups, towers |
| 8 | [`6_SIDE_RESULTS.md`](hilbert13/ladder/6_SIDE_RESULTS.md) | independent proofs, dead ends, the arithmetic question |
| 9 | [`gpt/README.md`](hilbert13/ladder/gpt/README.md) | what GPT proved, and where it is used |
| — | [`LITERATURE.md`](hilbert13/ladder/LITERATURE.md), [`QUESTIONS_FOR_GPT.md`](hilbert13/ladder/QUESTIONS_FOR_GPT.md) | sources; open questions |

Chapter 5 can be read straight after Chapter 1. The Lean project [`hilbert13/`](hilbert13/README.md) is self-contained.

## 4. Notation

| symbol | meaning |
|---|---|
| $G=A_7$; $(a,b,c)$ | a generating triple of orders $2,4,7$ with $abc=1$. In code $a=(01)(23)$ on letters $0..6$; in the text $(12)(34)$ on $1..7$ |
| classes 0, 1, 12, 14 | the four $A_7$-classes of triples (indices into `triples_data.triples`); $\{0,1\}$ and $\{12,14\}$ are $S_7$-orbits |
| $C$, $D$, $H$, $T$ | the genus-136 curve; $D=C/\langle\tau\rangle$ (genus 64); $H=C_G(\tau)/\langle\tau\rangle\cong C_2\times S_3$; $T=D/H$ (genus 3) |
| good, bad involutions | the involutions of $H$ with 18, resp. 2, fixed points on $D$ (§1.2) |
| irreducibles | $1,6,10,\overline{10},14_a,14_b,15,21,35$, with $14_a=14_{(5,2)}$; orientation conventions swap $10\leftrightarrow\overline{10}$ |
| $E_1$, $\lambda_1$ | first eigenspace and eigenvalue of the Laplacian; classes 0, 1: $E_1\cong14_a$, $\lambda_1\approx0.34627$ |
| $Q_0,Q_1,Q_2$ | the three sign-twisted quotient problems that see every irreducible (§2.3) |
| $D_2,D_4,D_7$ | the reduced fibres over the branch points, of degrees $1260,630,360$ |
| $B$, $T$ | $B=2D_7-D_4$ (degree 90); $T=D_2-2D_4$ (2-torsion); $B+T$ is linearised with $h^0\ge10$ |
| invariant vs linearised | a class fixed by $G$, versus one carrying a $G$-action; the obstruction (Mumford class) lies in $\mathbb Z/6$ |
| $L_{60}$, $\varphi$, $X_3$ | the degree-60 class of Mumford order 3; $\varphi:C\to\mathbb P^5$ given by its $\mathbf 6$; the Laza–Zheng cubic fourfold containing $\varphi(C)$ |
| $E_\pm(\tau)$ | the eigenspaces of the order-2 lift $\hat\tau$ on the $\mathbf 6$; $\lvert E_-\rvert$ is the $\tau$-pencil |
| $\gamma$, $\mu$, $a$, $\tilde\mu(C)$ | least gonality of a faithful $G$-curve; least degree of a linearised moving bundle; least accessory degree; least degree of an invariant class with $h^0\ge2$ |
| $\mathrm{Am}_G(X)$, $\mu_A(C)$, $c^{\rm st}_X(C)$ | the Amitsur subgroup of a base; least degree of a moving class with obstruction in $A$; least stable compression degree. $c^{\rm st}_X=\mu_{\mathrm{Am}(X)}$ (§5.5) |
| $E_{15}=C/A_5$, $E_{21}=C/L_2(5)$ | the elliptic quotients carrying 15 and 21 (GPT's "$E$" is $E_{21}$) |

## 5. Dependencies and tags

```
Ch.1 group data, fixed points ───────────────┬──────────────────────────┐
  │                                          │                          │
Ch.2 certified λ1 ≥ 0.34089 → gon ≥ 24       Ch.5 a(A7)=60, μ ≤ 90      §6.1 gon(D) ≥ 10; §4.5 g ≥ 266 ⇒ gon ≥ 25 (algebra only)
  │   (E1 = 14a, λ' ≥ 0.55998)               │
  ▼                                          Ch.5 μ ≠ 72, 84 (exact) → μ = 90
Ch.3 harmonic Hersch → gon ≥ 25              │
  ├──► Ch.4 γ(A7) ≥ 25 (+ Ch.2 for 10 other signatures, + the window §4.4)
  └──► Ch.7 Thm 7.4 (uses gon ≥ 25) ──► Cor 5.16 (60 or 90, by the Amitsur subgroup)
       Ch.7 Cor 7.6 gon ≤ 42 (exact; independent of every lower bound)
       Ch.7 §7.12 normalization defects + exact eigenlines → embedding → Prop 7.7 exactly 42 [P][X]
       Ch.7 §7.12 quadratic systems → no quadrics → Prop 7.9 Klein divisors [P][X]
       Ch.7 Cor 7.14 a pencil below 42 has a small kernel ([P] given [N] quotient data)
```

Every statement carries a tag:
- **[P]** paper proof in the chapter;
- **[X]** exact integer, rational, finite-field or algebraic-number computation; floating-point character sums rounded to integers retain [N] status unless independently replaced by an exact check;
- **[C]** computer-assisted proof (ball arithmetic, or floating point with a-priori error bounds);
- **[L]** checked in Lean;
- **[N]** numerical only;
- **[G]** GPT's proof, not re-derived here.

The lower bound 25 is [P][C]. The upper bound, the displayed pencil's exact degree 42, the embedding, and the Klein contact divisors are [P][X]. Exact degree 42 does not assert gonality 42. $a=60$ and $\mu=90$ are [P][X].

## 6. Verification

```sh
pip install numpy scipy sympy python-flint cvxopt pypdf
cd hilbert13/ladder && python3 cover.py        # 1 s smoke test
```

[`hilbert13/ladder/README.md`](hilbert13/ladder/README.md) lists every script, its claim and its command. Every saved output reproduces. The finite-element jobs need about 7 GB each; run them one at a time. Lean: `HANDOFF.md` §Setup.

## 7. Pitfalls

- **Invariant is not linearised.** Most of Chapter 7 lives in the gap between $15\mathbb Z$ and $90\mathbb Z$.
- **Open bases.** On a non-projective base, invertible functions can cancel a Schur obstruction, so $\mathrm{Pic}(B)=0$ does not force linearisation (Example 5.12). Apply Theorem 5.11 on a projective model.
- **Parity is not order.** A sign argument fixes vanishing orders mod 2, not their values. Proposition 7.7 uses immersion, proved by normalization defects, to obtain simple zeros.
- **"Algebraic" bounds.** $\operatorname{gon}\ge25$ for $g\ge266$ (Theorem 4.9) and $\operatorname{gon}(D)\ge10$ (§6.1) are purely algebraic. Below genus 266 the spectral input is needed.
- **Accessory conventions.** $a=60$ lets the monodromy drop after the accessory; with connected full monodromy the answer is $\mu=90$. The bound $\le59$ holds in both.
- **[N] claims** (the proposed ideal and remaining restriction ranks in Prop. 7.8, the Hessian section, and the inputs of Thm 7.13) are numerics, not proofs. Embedding, the exact displayed pencil, and the Klein contact divisors now have independent exact proofs.
- **GPT material** is input, not authority; `gpt/README.md` records what has been re-derived.
- **Labels.** Letters are $0..6$ in code and $1..7$ in the text; class labels are list indices, not invariants.
