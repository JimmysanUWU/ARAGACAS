# GUIDE: how to read this repository

This guide is for a reader, human or AI, who starts with no context. It covers:
- what the mathematics is, in one page (§2);
- what to read, and in what order (§3);
- what each object is (§4);
- which results depend on which (§5);
- how to check them (§6).

The current status, open problems and working rules are in [`HANDOFF.md`](HANDOFF.md).

## 1. What the repository is

It is research on Hilbert's 13th problem through curves with an $A_7$-action. The algebraic form of the problem asks whether the general septic needs algebraic functions of three variables, i.e. whether $\mathrm{RD}(7)=3$. That remains open. Here it is approached through two invariants of $A_7$:
- **gonality**: the least degree of a map to $\mathbb P^1$ from a curve with a faithful $A_7$-action;
- **accessory degree**: the least degree of a single finite extension after which the generic $A_7$-torsor descends to a curve.

**Headline results.** Let $C$ be any $(2,4,7)$ $A_7$-curve (genus 136, the minimum).
1. $25\le\operatorname{gon}(C)\le42$, and $13\le\operatorname{gon}(C/\tau)\le21$ for every involution $\tau$.
2. $25\le\gamma(A_7)\le42$, where $\gamma$ is the least gonality of any faithful $A_7$-curve.
3. $a(A_7)=60$, i.e. $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$, and this is sharp. The published bound was 6.
4. $\mu(A_7)=90$: the threshold when the full $A_7$-monodromy must stay connected. Over other bases the threshold is 60 or 90, according to the base's Amitsur subgroup.

## 2. The mathematics in one page

Four mechanisms do all the work. Each chapter opens with a paragraph headed **Idea**; read it first.

**Fixed points and Castelnuovo–Severi** (Ch. 1, §6.1).
- An involution fixes 18 points. That fixes the genus of every quotient of $D=C/\tau$.
- Castelnuovo–Severi forces a pencil of degree $\le8$ to be invariant under a group that is too large, so $\operatorname{gon}(D)\ge9$.
- At degree 9 the pencil's orbit has three possible shapes, each forcing a linear equivalence $(\star)$ on a genus-3 curve.
- Ramification transport (§6.1) goes one step further, using all of $A_7$:
  - saturated determinants give exact Picard identities;
  - averaging over a centraliser removes the choice of pencil;
  - amalgamating normalisers makes the class invariant;
  - a free $C_5$ forbids its degree, 1296.
- *Limit:* the $H$-action on $D$ alone cannot decide degree 9 (the audit curve).

**Degree versus energy** (Ch. 2–3).
- A balanced map of degree $d$ has energy $8\pi d$, so $\lambda_1$ bounds $d$ (Li–Yau). Certifying $\lambda_1$ on the genus-136 surface gives 24.
- The degree is also a cubic form, and it vanishes on the first eigenspace because $\wedge^3 14_{(5,2)}$ has no invariants. Harmonic Hersch makes this quantitative and gives 25.
- *Limit:* the spectral route stops at 25 or 26 here.

**Schur obstructions** (Ch. 5, 7).
- An invariant line bundle need not carry a $G$-action. The obstruction lives in $H^2(A_7,\mathbb C^*)=\mathbb Z/6$. On a $(2,4,7)$ curve, invariant classes have degrees $15\mathbb Z$ and linearised classes $90\mathbb Z$.
- **Upper bound 42.** The first twisted class with sections has degree 60. It embeds $C$ in $\mathbb P^5$, and an involution's $(-1)$-eigenspace gives a pencil of degree 42.
- **Pricing symmetry (§7.10).** A pencil's class stabiliser acts on the pencil through $PGL_2$, and the kernel $N$ of that action costs degree $|N|\operatorname{gon}(C/N)$. Gonalities of quotients come from differentials on $C/K$, syzygies, and Castelnuovo–Severi through the subgroup lattice. They put the cost at $\ge42$ unless $N$ is trivial or one of six small groups.
- **Accessories.** Over a linear base, a compression is a single linearised moving series, so accessory degrees are degrees of linearised series induced from subgroups:
  - $a=15\cdot4=60$, through the Klein quartic of $L_2(7)$;
  - $\mu=90$ with full monodromy.
- **Other bases.** The base can cancel exactly its Amitsur subgroup's worth of obstruction, which gives 60 or 90.
- *Limit:* in towers the Schur–Brauer obstruction is cheap to remove (a quadratic then a cubic step). What remains is the Albanese part.

**Castelnuovo theory of pencil orbits** (Ch. 4).
- A low-degree pencil has at least 15 conjugates. Pairs give $\mathbb P^1\times\mathbb P^1$ models; triples give $\mathbb P^7$ models or singularity costs.
- This bounds the genus of any faithful $A_7$-curve of gonality $\le24$ by 397. The certificates of Chapters 2–3 then cover every such genus.

## 3. Reading order

| step | file | time | what you get |
|---|---|---|---|
| 1 | this guide | 15 min | the mathematics in one page, notation, the dependency map |
| 2 | [`hilbert13/ladder/README.md`](hilbert13/ladder/README.md) | 5 min | every result with its status tag and chapter; script table |
| 3 | [`1_CURVE.md`](hilbert13/ladder/1_CURVE.md) | 20 min | the group, the curves, fixed points, $D=C/\tau$; the first algebraic bounds; the degree-9 structure theorem |
| 4 | [`2_SPECTRAL.md`](hilbert13/ladder/2_SPECTRAL.md) → [`3_HARMONIC_HERSCH.md`](hilbert13/ladder/3_HARMONIC_HERSCH.md) | 40 min | the lower bound 25: certified $\lambda_1$, then the cubic-form argument |
| 5 | [`7_TWISTED.md`](hilbert13/ladder/7_TWISTED.md) | 45 min | the upper bound 42: Schur-twisted bundles, the curve in $\mathbb P^5$; which pencils could beat 42 (§7.10) |
| 6 | [`4_LARGE_GENUS.md`](hilbert13/ladder/4_LARGE_GENUS.md) | 15 min | from one curve to all faithful $A_7$-curves |
| 7 | [`5_ACCESSORY.md`](hilbert13/ladder/5_ACCESSORY.md) | 40 min | $a=60$, $\mu=90$, Amitsur subgroups, and what survives in towers |
| 8 | [`6_SIDE_RESULTS.md`](hilbert13/ladder/6_SIDE_RESULTS.md) | 15 min | independent proofs, dead ends, the arithmetic question |
| 9 | [`gpt/README.md`](hilbert13/ladder/gpt/README.md) | 15 min | what Sol and Astra (GPT) proved, and where it is used |
| — | [`LITERATURE.md`](hilbert13/ladder/LITERATURE.md), [`QUESTIONS_FOR_GPT.md`](hilbert13/ladder/QUESTIONS_FOR_GPT.md) | as needed | sources with read/unverified tags; open questions |

Chapter 5 can be read straight after Chapter 1. The Lean project in [`hilbert13/`](hilbert13/README.md) is self-contained: three theorems on single superpositions, plus the logical core of the spectral certificate.

## 4. Notation

| symbol | meaning |
|---|---|
| $G=A_7$; $(a,b,c)$ | a generating triple of orders $2,4,7$ with $abc=1$. In code $a=(01)(23)$ on letters $0..6$; in the text $(12)(34)$ on $1..7$ |
| classes 0, 1, 12, 14 | the four $A_7$-classes of triples, as indices into `triples_data.triples`. $\{0,1\}$ and $\{12,14\}$ are $S_7$-orbits |
| $C$, $D$, $H$, $T$ | the genus-136 curve; $D=C/\langle\tau\rangle$ (genus 64); $H=C_G(\tau)/\langle\tau\rangle\cong C_2\times S_3$; $T=D/H$ (genus 3) |
| good, bad involutions | the involutions of $H$ with 18 fixed points on $D$, and those with 2 (§1.2) |
| irreducibles | $1,6,10,\overline{10},14_a,14_b,15,21,35$, with $14_a=14_{(5,2)}$. Orientation conventions swap $10\leftrightarrow\overline{10}$ |
| $E_1$, $\lambda_1$ | first eigenspace and eigenvalue of the hyperbolic Laplacian. Classes 0, 1: $E_1\cong14_a$ and $\lambda_1\approx0.34627$ |
| $Q_0,Q_1,Q_2$ | the three sign-twisted quotient problems that together see every irreducible (§2.3) |
| $D_2,D_4,D_7$ | the reduced fibres over the branch points, of degrees $1260,630,360$ |
| $B$, $T$ | $B=2D_7-D_4$ (degree 90); $T=D_2-2D_4$, a 2-torsion class. $B+T$ is the linearised degree-90 series with $h^0\ge10$ |
| invariant vs linearised | a class fixed by $G$, versus one carrying a $G$-action. The obstruction (Mumford class) lies in $H^2(A_7,\mathbb C^*)=\mathbb Z/6$ |
| $L_{60}$, $\varphi$, $X_3$ | the degree-60 invariant class, with Mumford class of order 3; its embedding $\varphi:C\hookrightarrow\mathbb P^5$; the Laza–Zheng cubic fourfold containing $\varphi(C)$ |
| $E_\pm(\tau)$ | the eigenspaces of a lift $\hat\tau$ on the $\mathbf 6$. $\lvert E_-\rvert$ is the $\tau$-pencil of degree 42 |
| $\gamma(G)$, $\mu(G)$, $a(G)$, $\tilde\mu(C)$ | least gonality of a faithful $G$-curve; least degree of a linearised moving bundle; least accessory degree; least degree of an invariant class with $h^0\ge2$ |
| $\mathrm{Am}_G(X)$, $\mu_A(C)$, $c^{\rm st}_X(C)$ | the Amitsur subgroup of a base; the least degree of a moving invariant class with obstruction in $A$; the least stable compression degree over $X$. They satisfy $c^{\rm st}_X=\mu_{\mathrm{Am}(X)}$ (§5.5) |
| $E_{15}=C/A_5$, $E_{21}=C/L_2(5)$ | the elliptic quotients carrying 15 and 21; not to be confused with the eigenspace $E_1$. GPT's "$E$" is $E_{21}$ |

## 5. How the results depend on each other

```
Ch.1 group data, fixed points ───────────────┬──────────────────────────┐
  │                                          │                          │
Ch.2 certified λ1 ≥ 0.34089 → gon ≥ 24       Ch.5 a(A7)=60, μ ≤ 90      §6.1 gon(D) ≥ 10
  │   (and E1 = 14a, λ' ≥ 0.55998)              (algebra + ATLAS)           (algebra only; Astra)
  ▼                                          │
Ch.3 harmonic Hersch → gon ≥ 25 (all four    Ch.5 μ ≠ 72, 84 (exact script) → μ = 90
  │   classes)                               │
  ├──► Ch.4 γ(A7) ≥ 25 (+ Ch.2 for 10 other signatures, + the window §4.4)
  │                                          │
  └──► Ch.7 Thm 7.4 (uses gon ≥ 25 via Clifford) ──► Cor 5.16 (60 or 90, by the Amitsur subgroup)
       Ch.7 Cor 7.6 gon ≤ 42 (independent of every lower bound)
       Ch.7 Prop 7.7 exactly 42 (uses gon(C/τ) ≥ 13; simple base points [P] at 12 of 18, [N] at 6)
```

**What is proved, and how.** Every statement carries a tag:
- **[P]** a paper proof, written out in the chapter;
- **[X]** an exact finite computation;
- **[C]** a computer-assisted proof, in ball arithmetic or floating point with a-priori error bounds;
- **[L]** checked in Lean;
- **[N]** numerical only;
- **[G]** GPT's proof, not re-derived here.

The lower bound 25 is [P]+[C]. The upper bound 42 is [P]+[X], plus [N] for exactness. $a=60$ and $\mu=90$ are [P]+[X]. Nothing rests on an [N] or [G] claim except where the tag says so.

**Rigor sweep (3 October 2026).** Every [P] statement was checked to have its proof in the chapter.
- **Restored or added:**
  - the proof of the degree-9 structure theorem (Thm 1.5);
  - the proofs of the audit curve (Thm 1.7) and of Lemma 5.5;
  - the degree formula of Prop. 7.1;
  - the genus 691 used in §7.7.
- **Corrected:** Theorem 5.11 (open bases) and Prop. 7.7 step 1 (parity).

## 6. How to verify

```sh
pip install numpy scipy sympy python-flint cvxopt pypdf
cd hilbert13/ladder
python3 cover.py                 # 1 s smoke test
```

- `hilbert13/ladder/README.md` lists every script, its claim, and the command that reproduces its saved `*_output.txt`.
- Times there are for an otherwise idle 4-core machine.
- The finite-element jobs need about 7 GB each, so run them one at a time.
- Lean: see `HANDOFF.md` §Setup. `lake build` in `hilbert13/` takes a few minutes with the Mathlib cache.

**Last full check (2 October 2026).** Every script was re-run, and every output was compared with the saved one:
- 17 of the 19 output files are identical apart from timings, and so are the main certificate and the `certify_th.py` runs;
- the two signature certificates agree to $4\cdot10^{-16}$, floating-point noise far below every margin;
- the window run stops at $(4,4,7)$, as recorded in §4.4;
- the Lean project builds with no `sorry`.

Since then only `twisted_rr.py` has changed (it adds the $\mathrm{Sym}^3\mathbf 6$ decomposition), and its output was regenerated.

## 7. Pitfalls

- **Invariant is not linearised.** Degrees of invariant classes form $15\mathbb Z$; linearised ones form $90\mathbb Z$. Most of Chapter 7 lives in the gap between them.
- **Open bases.** On a non-projective base, invertible functions can cancel a Schur obstruction, so $\mathrm{Pic}(B)=0$ does not force linearisation (Example 5.12). Apply Theorem 5.11 on a projective model.
- **Parity is not order.** A sign argument fixes vanishing orders mod 2, not their values. Plücker weights decide more (Prop. 7.7, step 1).
- **"Algebraic" bounds.** $\operatorname{gon}\ge23$ in Theorem 4.1 still uses the spectral input for genus $\le335$. Only $\operatorname{gon}(D)\ge10$ (§6.1) is purely algebraic.
- **Accessory conventions.** $a(A_7)=60$ allows the monodromy to drop after the accessory. With connected full monodromy the answer is $\mu=90$. The bound $\le59$ holds in both conventions.
- **[N] claims** (Props. 7.7–7.9, the inputs of Thm 7.13) are high-precision numerics, not proofs. The bound $\operatorname{gon}\le42$ itself does not depend on them.
- **GPT material** is input, not authority. Use a GPT claim only after re-deriving it; `gpt/README.md` records which ones have been.
- **Code labels.** Letters are $0..6$ in code and $1..7$ in the text. Class labels are list indices, not invariants.
- **Killing jobs.** `pkill -f pattern` also kills the calling shell. Use `pgrep -f "[p]attern" | xargs -r kill`.

## 8. Where the frontier is

- **The exact gonality**, somewhere in $[25,42]$, is open; every symmetric construction stops at 42. A pencil below 42 must have kernel $N\in\{1,C_2,C_3,C_4,V_4,S_3,C_7\}$ (Cor. 7.14). The most promising new inputs are:
  - stable reduction at $p=7$ with graph gonality (§7.11);
  - immersion at the 630 four-points (§7.11). It would make Prop. 7.7 a paper proof, and a failure would give $\operatorname{gon}(C/\tau)\le15$.
- **The tower problem**, $\mathrm{RD}(A_7)>1$, is reduced in §5.5 to correspondences whose Albanese part $u_Z$ is nonzero. The Schur–Brauer part is cheap (Cor. 5.18).

The prioritised list is in `HANDOFF.md`.
