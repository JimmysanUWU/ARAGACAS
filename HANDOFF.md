# HANDOFF — Hilbert 13 / A7 gonality ladder

Written at the end of the first session (2026-09-29) so the next session can continue without it.
**Read this first**, then `hilbert13/ladder/NOTES.md` (all proofs) and `hilbert13/ladder/LIT.md`.

## 0. Restore

The work is 6 commits on branch `claude/hopeful-allen-my45u6`, on top of base commit `299078f`
("Create README.md"), which exists in both `x5ilky/ARAGACAS` and the fork `JimmysanUWU/ARAGACAS`.

```sh
# in a clone of JimmysanUWU/ARAGACAS, with aragacas-ladder.bundle in the repo root:
git fetch aragacas-ladder.bundle claude/hopeful-allen-my45u6:claude/hopeful-allen-my45u6
git checkout claude/hopeful-allen-my45u6
git push -u origin claude/hopeful-allen-my45u6      # then open a draft PR
```
Fallback without git: `aragacas-transfer.tar.gz` holds the same files as a plain snapshot.

Why a new session: pushing to `x5ilky/ARAGACAS` failed (403, Claude GitHub App not installed on
Ethan's account), and the fork could not be attached to the old session (same repo name).

## 1. Environment setup (about 5 minutes)

```sh
# Lean 4 + Mathlib (needs release.lean-lang.org, github.com and the Mathlib cache hosts allowed)
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o elan-init.sh
bash elan-init.sh -y --default-toolchain leanprover/lean4:stable
export PATH=$HOME/.elan/bin:$PATH
cd hilbert13 && lake exe cache get && lake build     # the project pins Lean/Mathlib v4.34.1

# Python (group theory, character table, FEM spectrum)
pip install numpy scipy sympy python-flint pdfminer.six
```

## 2. Files

| path | what |
|---|---|
| `hilbert13/Hilbert13/Superposition.lean` | first mini-project (Lean, no `sorry`): single superpositions $g(\varphi(x)+\psi(y))$ — every function with arbitrary inner functions; $xy$ impossible with continuous inner functions; two terms suffice |
| `hilbert13/README.md` | paper proofs for that mini-project |
| `hilbert13/ladder/NOTES.md` | **main notes**: rungs 1–6, statuses [P]/[C]/[?]/[X] |
| `hilbert13/ladder/LIT.md` | literature map |
| `hilbert13/ladder/a7.py` | brute-force A7 toolkit (permutations of 0..6, `mul(p,q)=p∘q`) |
| `triples.py`, `triples_data.py`, `classreps.py` | (2,4,7) generating triples: 96 with $a=(01)(23)$ fixed, 4 $A_7$-classes = indices 0, 1, 12, 14 |
| `rung1.py` | fixed points computed from the branch orbits |
| `quot.py` | the group $H=C(\tau)/\langle\tau\rangle$ on $D$: fixed points and quotient genera |
| `mingenus.py` | (2,4,7) is the minimal-genus generating signature of $A_7$ (genus 136) |
| `chartab.py` | $A_7$ character table (Dixon method, checked orthonormal) |
| `jacobian.py` | Chevalley–Weil decomposition of $H^1(C)$ and $\mathrm{Jac}(D)$ |
| `hdecomp.py` | $H$-isotypic structure of $\rho^\tau$ for each irreducible $\rho$ |
| `cusp.py`, `star.py` | isotypic projections of divisors supported on the branch fibres |
| `model99.py` | the audit curve $D'$ (explicit smooth $H$-invariant (9,9) curve with gonality 9) |
| `spectrum.py` | P1-FEM Laplace spectrum of $C$: `python3 spectrum.py <n> <triple-index>` |

## 3. The task (the user's ladder, verbatim in substance)

Let $C$ be a smooth connected complex curve with a faithful $A_7$-action, $C/A_7\cong\mathbb P^1$, with exactly
three branch points of inertia orders 2, 4, 7. Fix an involution $\tau\in A_7$, $D=C/\langle\tau\rangle$;
$\operatorname{gon}(X)$ = least degree of a nonconstant map $X\to\mathbb P^1$. Work independently; if a rung resists,
state the exact implication you cannot justify; invent frameworks freely.
1. Compute $g(C)$, the fixed points of an involution (from the branch orbits), and $g(D)$.
2. For an involution $\mu\ne\tau$ commuting with $\tau$: how it descends to $D$ and its fixed points there.
3. Prove or refute $\operatorname{gon}(D)\ge9$; distinguish $f\circ v=f$ from $f\circ v=M\circ f$.
4. Assume a degree-9 map on $D$: extract every rigorous consequence (fixed points, other degree-9
   maps obtained by automorphisms); find a contradiction or a precise description.
5. **First summit:** settle $\operatorname{gon}(C/\langle\tau\rangle)\ge10$ for every such $C$ and $\tau$ (all generating triples).
6. **Further summit:** investigate $\operatorname{gon}(C)\ge17$ with a fresh audit — what the quotient result
   contributes, what is independent, prove or disprove parts.

The user's standing instruction: "Aim for QED of the ladder — or your own invention of a mathematical
structure / novel framework." They also want Lean formalization where feasible. The user and GPT
say they have an audited proof of the first summit ("rung 4 is where it reaches the door to our
proof"); we have not seen it.

## 4. Results so far (details and proofs in NOTES.md)

- **Rung 1 [P][C]:** $g(C)=136$; an involution fixes $12+6=18$ points; $g(D)=64$. Order 4: 2 fixed
  points; order 7: 3; orders 3, 5, 6: none. $A_7$ has 4 conjugacy classes of generating triples, 2 up to
  $S_7$; (2,4,7) curves are the minimal-genus $A_7$-curves (Conder: 136).
- **Rung 2 [P][C]:** $H=C(\tau)/\langle\tau\rangle\cong C_2\times S_3$ acts on $D$. There are 4 **good** involutions
  (18 fixed points, $g(D/v)=28$; one central), 3 **bad** involutions (lifting to order 4; 2 fixed points,
  $g=32$), and order 3 and 6 elements act freely ($g(D/C_3)=22$). $g(D/H)=3$.
- **Rung 3 [P]:** $\operatorname{gon}(D)\ge9$. For degree $\le8$, CS forces $f\circ v=f$ for every good $v$; they generate $H$,
  so $12\mid d$.
- **Rung 4 [P] (modulo Lemma 4.0, a standard fact on (9,9) curves, citation to verify):** at degree 9
  every good involution gives CS **equality**, and $D=\{s_0x^2+s_1xy+s_2y^2=0\}\subset E_v\times\mathbb P^1$. Pencil
  trichotomy: (a) $H$-invariant pencil (then $H\to\mathrm{PGL}_2$ is faithful, dihedral of order 12);
  (b) a factor through a genus-4 $H$-curve $Y$; (c) $D$ is a smooth (9,9) curve. Each case exhibits the
  $H$-cover $D\to T=C/C(\tau)$ (genus 3) as a **pullback**.
- **Rung 5: not settled.** Theorem 4.5 [P]: $\operatorname{gon}(D)\le9\Rightarrow(\star)$: $g_1+\dots+g_9\sim3(t_1+t_2+t_3)$ on $T$, where the
  $t_i$ are images of the central good involution's fixed points and the $g_j$ of the other good ones.
  - $\mathrm{Jac}(C)\sim A^{10}E_1^{15}E_2^{21}S^{35}$, $\mathrm{Jac}(D)\sim A^4E_1^7E_2^{11}S^{17}$, $\mathrm{Jac}(T)\sim E_2\times S$
    [C]; the irreducibles $6,14_a,14_b$ are absent from $H^1(C)$. The divisor behind $(\star)$ has nonzero
    21- and 35-components, so $(\star)$ is a genuine arithmetic relation (Mordell–Weil, torsion) —
    open.
  - **Audit Theorem 5.1 [P][C]:** an explicit smooth $H$-curve $D'$ of genus 64 with *identical*
    fixed-point data has gonality 9. So no argument using only the $H$-action on $D$ can prove the
    first summit; a correct proof must use more of $C$.
- **Rung 6:** $\operatorname{gon}(C)\ge13$ rigorously (FW Lemma 2.2: $A_7$ is simple, so translates of $f$ generate $k(C)$,
  and $136>11^2$). **Spectral framework [P]:** by Li–Yau, $\operatorname{gon}(C)\ge67.5\,\lambda_1(C)$ and
  $\operatorname{gon}(C/\langle\tau\rangle)\ge33.75\,\lambda_1(C)$. So $\lambda_1(C)>0.2667$ proves **both summits**.
  **Numerics [C, not rigorous]:**

  | triple class | n=3 | n=5 | multiplicity |
  |---|---|---|---|
  | 0 (same spectrum as 1) | 0.34560 | 0.34604 | 14 |
  | 12 (same spectrum as 14) | 0.35847 | — | 14 |

  Richardson extrapolation from $n=3,5$ gives $\lambda_1\approx0.3463$ for class 0. If accurate:
  $\operatorname{gon}(C)\ge24$ and $\operatorname{gon}(D)\ge12$ for every such curve.

## 5. Facts derived but not yet in NOTES.md [P]

- $C/A_6\cong\mathbb P^1$ (genus 0): $C$ is the Galois closure of a **degree-7 Shabat polynomial** $p(x)$ with
  passport $[2^21^3,\ 4\,2\,1,\ 7]$ (7-cycle at $\infty$). Some polynomials with this passport have
  monodromy $\mathrm{PSL}_3(\mathbb F_2)$ instead; keep only the $A_7$ ones.
- $T=C/C(\tau)$ maps with **degree 3** to $C/N(V_I)\cong\mathbb P^1_w$, where $N(V_I)=(S_4\times S_3)\cap A_7$ (order 72;
  $C/N(V_I)$ has genus 0) and $w$ = sum of 4 roots of $p(x)=t$; the triple cover is the resolvent cubic
  of those 4 roots. Two fibres: $t_0+t_1+t_2$ over $w_a$ and $t_b+2t_3$ over $w_b$ ($t_0$ = image of the
  $\tau$-fixed points over the order-2 branch point, which is also a branch point of $C\to D$). So
  $t_0+t_1+t_2\sim t_b+2t_3$ on $T$.

## 6. Recommended next steps (priority order)

1. **Converge $\lambda_1$.** Run `spectrum.py` at $n=8$ and $n=12$ for classes 0 and 12, one run at a time
   (running several in one shell call got a process OOM-killed), then extrapolate. $n=5$ has 62,730
   nodes and takes about 50 s.
2. **Certify $\lambda_1(C)>0.2667$** — this would prove both summits. Idea: decompose $L^2(C)$ into
   $A_7$-isotypic pieces (twisted Laplacians on the (2,4,7) orbifold, rep dimension $\le35$) and apply a
   guaranteed-lower-bound FEM (Liu–Oishi Crouzeix–Raviart type) adapted to the hyperbolic metric,
   with interval arithmetic. The margin is about 30%. Also double-check the Li–Yau constant
   ($8\pi d$) against the original paper.
3. **Decide $(\star)$ explicitly:** compute $p(x)$ exactly, build $T$ via the resolvent cubic over
   $\mathbb P^1_w$ (needs the degree-35 genus-0 map $\mathbb P^1_w\to\mathbb P^1_t$), locate $t_i, g_j$, and test linear
   equivalence (numerical Abel–Jacobi, or Riemann–Roch mod $p$). If $(\star)$ fails for all classes,
   the first summit is proved algebraically.
4. **Lean:** Mathlib has no Castelnuovo–Severi or Riemann–Hurwitz for curves, so formalize the finite
   skeleton (genus and fixed-point arithmetic, $A_7$ counts by `decide`, the case exclusions in
   4.2–4.5).
5. Verify the citation for Lemma 4.0 (Martens 1996) or write out its proof.

## 7. Infrastructure notes

- The auto-mode safety classifier intermittently returned "no verdict"; retrying once usually worked.
- The `Claude_Code_Remote` connector (`send_later`) needed re-authorization in claude.ai settings.
- Network: arxiv via `curl` works; `WebFetch` is blocked for arxiv.org and openai.com.
