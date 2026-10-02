# GUIDE: how to read this repository

This guide is written for a reader, human or AI, who starts with no context. It covers what to read, in what order, what each object is, which results depend on which, and how to check them. The current status, open problems and working rules are in [`HANDOFF.md`](HANDOFF.md).

## 1. What the repository is

It is research on Hilbert's 13th problem through curves with an $A_7$-action. The algebraic form of the problem asks whether the general septic needs algebraic functions of three variables, i.e. whether $\mathrm{RD}(7)=3$. That remains open; here it is approached through two invariants of $A_7$:

- **gonality**: the least degree of a map to $\mathbb P^1$ from a curve with a faithful $A_7$-action;
- **accessory degree**: the least degree of a single finite extension after which the generic $A_7$-torsor descends to a curve.

**Headline results.** Let $C$ be any $(2,4,7)$ $A_7$-curve (genus 136, the minimum).
1. $25\le\operatorname{gon}(C)\le42$, and $13\le\operatorname{gon}(C/\tau)\le21$ for every involution $\tau$.
2. $25\le\gamma(A_7)\le42$, where $\gamma$ is the least gonality of any faithful $A_7$-curve.
3. $a(A_7)=60$, i.e. $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$, and this is sharp. The published bound was 6.
4. $\mu(A_7)=90$: the same threshold when the full $A_7$-monodromy must stay connected.

## 2. Reading order

| step | file | time | what you get |
|---|---|---|---|
| 1 | this guide | 10 min | orientation, notation, the dependency map |
| 2 | [`hilbert13/ladder/README.md`](hilbert13/ladder/README.md) | 5 min | every result with its status tag and chapter; script table |
| 3 | [`1_CURVE.md`](hilbert13/ladder/1_CURVE.md) | 15 min | the group, the curves, fixed points, $D=C/\tau$; the first algebraic bounds |
| 4 | [`2_SPECTRAL.md`](hilbert13/ladder/2_SPECTRAL.md) → [`3_HARMONIC_HERSCH.md`](hilbert13/ladder/3_HARMONIC_HERSCH.md) | 40 min | the lower bound 25: certified $\lambda_1$, then the cubic-form argument |
| 5 | [`7_TWISTED.md`](hilbert13/ladder/7_TWISTED.md) | 40 min | the upper bound 42: Schur-twisted bundles, the curve in $\mathbb P^5$ |
| 6 | [`4_LARGE_GENUS.md`](hilbert13/ladder/4_LARGE_GENUS.md) | 15 min | from one curve to all faithful $A_7$-curves |
| 7 | [`5_ACCESSORY.md`](hilbert13/ladder/5_ACCESSORY.md) | 30 min | $a=60$, $\mu=90$, and what survives in towers |
| 8 | [`6_SIDE_RESULTS.md`](hilbert13/ladder/6_SIDE_RESULTS.md) | 15 min | independent proofs, dead ends, the arithmetic question |
| 9 | [`gpt/README.md`](hilbert13/ladder/gpt/README.md) | 15 min | what the parallel GPT collaboration proved and where it is used |
| — | [`LITERATURE.md`](hilbert13/ladder/LITERATURE.md), [`QUESTIONS_FOR_GPT.md`](hilbert13/ladder/QUESTIONS_FOR_GPT.md) | as needed | sources with read/unverified tags; open questions |

Chapter 5 can be read straight after Chapter 1. The Lean project in [`hilbert13/`](hilbert13/README.md) is self-contained: three theorems on single superpositions, plus the logical core of the spectral certificate.

## 3. Notation

| symbol | meaning |
|---|---|
| $G=A_7$; $(a,b,c)$ | a generating triple of orders $2,4,7$ with $abc=1$; in code $a=(01)(23)$ on letters $0..6$, in the text $(12)(34)$ on $1..7$ |
| classes 0, 1, 12, 14 | the four $A_7$-classes of triples, as indices into `triples_data.triples`; $\{0,1\}$ and $\{12,14\}$ are $S_7$-orbits |
| $C$, $D$, $H$ | the genus-136 curve; $D=C/\langle\tau\rangle$ (genus 64); $H=C_G(\tau)/\langle\tau\rangle\cong C_2\times S_3$ |
| irreducibles | $1,6,10,\overline{10},14_a,14_b,15,21,35$, with $14_a=14_{(5,2)}$; orientation conventions swap $10\leftrightarrow\overline{10}$ |
| $E_1$, $\lambda_1$ | first eigenspace and eigenvalue of the hyperbolic Laplacian. Classes 0, 1: $E_1\cong14_a$, $\lambda_1\approx0.34627$ |
| $Q_0,Q_1,Q_2$ | the three sign-twisted quotient problems that together see every irreducible (§2.3) |
| $D_2,D_4,D_7$ | the reduced fibres over the branch points, of degrees $1260,630,360$ |
| $B$, $T$ | $B=2D_7-D_4$ (degree 90); $T=D_2-2D_4$, a 2-torsion class. $B+T$ is the linearised degree-90 series with $h^0\ge10$ |
| invariant vs linearised | a class fixed by $G$, versus one carrying a $G$-action. The obstruction (Mumford class) lies in $H^2(A_7,\mathbb C^*)=\mathbb Z/6$ |
| $L_{60}$, $\varphi$, $X_3$ | the degree-60 invariant class (Mumford class of order 3); its embedding $\varphi:C\hookrightarrow\mathbb P^5$; the Laza–Zheng cubic fourfold containing $\varphi(C)$ |
| $E_\pm(\tau)$ | the eigenspaces of a lift $\hat\tau$ on the $\mathbf 6$; $\lvert E_-\rvert$ is the $\tau$-pencil of degree 42 |
| $\gamma(G)$, $\mu(G)$, $a(G)$, $\tilde\mu(C)$ | least gonality of a faithful $G$-curve; least degree of a linearised moving bundle; least accessory degree; least degree of an invariant class with $h^0\ge2$ |
| $\mathrm{Am}_G(X)$, $\mu_A(C)$ | the Amitsur subgroup of a base (obstructions of its invariant classes); the least degree of a moving invariant class on $C$ with obstruction in $A$ (§5.5) |
| $E_{15}=C/A_5$, $E_{21}=C/L_2(5)$ | the elliptic quotients carrying 15 and 21 (not to be confused with the eigenspace $E_1$). GPT's "$E$" is $E_{21}$ |

## 4. How the results depend on each other

```
Ch.1 group data, fixed points ───────────────┬──────────────────────────┐
  │                                          │                          │
Ch.2 certified λ1 ≥ 0.34089 → gon ≥ 24       Ch.5 a(A7)=60, μ ≤ 90      §6.1 gon(D) ≥ 10
  │   (and E1 = 14a, λ' ≥ 0.55998)              (algebra + ATLAS)           (algebra only)
  ▼                                          │
Ch.3 harmonic Hersch → gon ≥ 25 (all four    Ch.5 μ ≠ 72, 84 (exact script) → μ = 90
  │   classes)                               │
  ├──► Ch.4 γ(A7) ≥ 25 (+ Ch.2 for 10 other signatures, + the window §4.4)
  │                                          │
  └──► Ch.7 Thm 7.4 (uses gon ≥ 25 via Clifford) ──► Cor 5.16 (60 or 90, by the Amitsur subgroup)
       Ch.7 Cor 7.6 gon ≤ 42 (independent of every lower bound)
       Ch.7 Prop 7.7 exactly 42 (uses gon(C/τ) ≥ 13)
```

**What is proved, and how.** Every statement carries a tag:
- **[P]** a paper proof in the chapter;
- **[X]** an exact finite computation;
- **[C]** a computer-assisted proof, in ball arithmetic or floating point with a-priori error bounds;
- **[L]** checked in Lean;
- **[N]** numerical only;
- **[G]** GPT's proof, not re-derived here.

The lower bound 25 is [P]+[C]. The upper bound 42 is [P]+[X], plus [N] for exactness. $a=60$ and $\mu=90$ are [P]+[X]. Nothing rests on an [N] or [G] claim except where the tag says so.

## 5. How to verify

```sh
pip install numpy scipy sympy python-flint cvxopt pypdf
cd hilbert13/ladder
python3 cover.py                 # 1 s smoke test
```

- `hilbert13/ladder/README.md` lists every script, its claim, and the command that reproduces its saved `*_output.txt`.
- Times there are for an otherwise idle 4-core machine.
- The finite-element jobs need about 7 GB each, so run them one at a time.
- Lean: see `HANDOFF.md` §Setup; `lake build` in `hilbert13/` takes a few minutes with the Mathlib cache.

**Last full check (2 October 2026).** Every script was re-run, and every output was compared with the saved one:
- 17 of the 19 output files are identical apart from timings, and so are the main certificate and the `certify_th.py` runs;
- the two signature certificates agree to $4\cdot10^{-16}$, floating-point noise far below every margin;
- the window run stops at $(4,4,7)$, as recorded in §4.4;
- the Lean project builds with no `sorry`.

## 6. Pitfalls

- **Invariant is not linearised.** Degrees of invariant classes form $15\mathbb Z$; linearised ones form $90\mathbb Z$. Most of Chapter 7 lives in the gap between them.
- **"Algebraic" bounds.** $\operatorname{gon}\ge23$ in Theorem 4.1 still uses the spectral input for genus $\le335$; only $\operatorname{gon}(D)\ge10$ (§6.1) is purely algebraic.
- **Open bases.** On a non-projective base, invertible functions can cancel a Schur obstruction, so $\mathrm{Pic}(B)=0$ does not force linearisation (Example 5.12). Apply Theorem 5.11 on a projective model.
- **Accessory conventions.** $a(A_7)=60$ allows the monodromy to drop after the accessory. With connected full monodromy the answer is $\mu=90$. The bound $\le59$ holds in both conventions.
- **[N] claims** (Props. 7.7–7.9) are high-precision numerics, not proofs. The bound $\operatorname{gon}\le42$ itself does not depend on them.
- **GPT material** is input, not authority. Use a GPT claim only after re-deriving it; `gpt/README.md` records which ones have been.
- **Code labels.** Letters are $0..6$ in code and $1..7$ in the text. Class labels are list indices, not invariants.
- **Killing jobs.** `pkill -f pattern` also kills the calling shell. Use `pgrep -f "[p]attern" | xargs -r kill`.

## 7. Where the frontier is

The exact gonality, somewhere in $[25,42]$, is open; every symmetric construction stops at 42. The most promising new input is stable reduction at $p=7$ with graph gonality (§7.10). The unrestricted tower problem, $\mathrm{RD}(A_7)>1$, is untouched beyond Theorem 5.11. The prioritised list is in `HANDOFF.md`.
