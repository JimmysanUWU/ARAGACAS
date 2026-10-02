# Review: *Ramification transport* (synthesis, 28 Sep 2026) and *A ladder to the septic* (verification note, 29 Sep 2026)

Reviewer: Claude, on the ARAGACAS branch. Finite checks: `review_checks.py`. The documents credit Sol and Astra
for the underlying work. This review treats the arguments, not the authors.

## Verdict at a glance

| Claim | Where | Verdict |
|---|---|---|
| $\operatorname{gon}(C/\langle\tau\rangle)\ge10$ for every $(2,4,7)$ $A_7$-curve and every involution | Synthesis Thm 2.3, Ladder §3 | **Correct.** Complete proof; every step checked below. |
| $\operatorname{gon}(D)\ge9$ | Synthesis (2.4) | **Correct.** Same as our rung 3. |
| $\mathrm{Aut}(C)=A_7$; no map of degree $\le18$ on $C$ has a double first factor | Cor. 2.4 | **Correct.** |
| Ramification transport (Thm 3.1), amalgamation (Thm 3.2) | §3.1–3.2 | **Correct** as stated. |
| $A_5$ genus-6 boundary example | §3.3 | **Correct**, checked by computer. It shows the connectivity hypothesis cannot be dropped. |
| Effective-defect identity (3.4)–(3.5) | §3.4 | Correct as identities. Correctly labelled as not yet an obstruction. |
| Transfer lemma (Lemma 4.1) | §4.1 | **Correct.** The generic-fibre choice also repairs a small gap in Farb–Wolfson Lemma 2.4 (a norm can be constant). |
| Genus gap and Table 1 | §4.2 | **Correct and complete.** All nine signatures are realised by generating triples. |
| $\mathrm{ed}_{\mathbb C}(A_7;\le17)>1$ | Prop 4.2 | **Correct**, given Karpenko–Merkurjev. Apparently new beyond Farb–Wolfson's $\le6$. |
| Lemma 5.1 (orbit field, common-field bound $B_m(e)$) | §5.1 | **Correct**; table re-derived. |
| Lemma 5.2 (three-pencil bound) | §5.2 | **Correct**; all numbers re-derived. |
| Lemma 5.3 (singularity budget) | §5.2 | Argument sound; table entries spot-checked (136/15, 136/16, 169/16–18, 199/18, 211/18). |
| Lemma 5.4 (Segre gap) | §5.3 | **Correct.** |
| Lemma 5.5 (two exceptional independent triples) | §5.3 | **Not verified.** It depends on the ranges of Petrakiev's refined Castelnuovo theorems. |
| Thm 6.1 ($\operatorname{gon}\ge19$ off genus 136; $\ge17$ at 136) | §6 | Rows 199–274 rest on verified or classical inputs. Row 169 at $m=18$ needs Lemma 5.5. Row 136 is superseded by our $\operatorname{gon}(C)\ge24$. |
| Torsion/norm model | §7 | Consistent with our independent rung-4/5 analysis. Not needed for any headline result. |
| Finite certificates | §8, Ladder p.1 | **All reproduced** (see below). |

**No mathematical error was found in the main line.** The minor issues are listed in §6 of this review.

## 1. The central theorem: $\operatorname{gon}(C/\langle\tau\rangle)\ge10$

**The argument, re-derived.** Let $f:C\to\mathbb P^1$ be the pull-back of a hypothetical degree-9 map on $D=C/\langle\tau\rangle$.
It has degree 18, is $\tau$-invariant, and has $L=f^*\mathcal O(1)$.

1. Take an involution $\mu\ne\tau$ commuting with $\tau$. The section
   $s_0\otimes\mu^*s_1-s_1\otimes\mu^*s_0\in H^0(L\otimes\mu^*L)$, of degree 36, vanishes at $p$ iff $f(p)=f(\mu p)$. That
   happens on $\mathrm{Fix}(\mu)$, and also on $\mathrm{Fix}(\tau\mu)$ because $f\circ\tau=f$. These are 18 + 18 distinct points: a common
   fixed point would have a non-cyclic stabiliser $V_4$.
2. The section is not identically zero. Otherwise $f$ would be invariant under $\langle\tau,\mu\rangle\cong V_4$, and $4\mid18$ would
   follow. So its divisor is exactly $R_\mu+R_{\tau\mu}$, and $[L]+[\mu^*L]=[R_\mu+R_{\tau\mu}]$.
3. Each translate $h^*L$ with $h\in C_G(\tau)$ comes from the $\tau$-invariant map $f\circ h$ and satisfies the same identity.
   With $M=\sum_{h\in C_G(\tau)}h^*L$ we have $\mu^*M=M$, and summing gives $2M=24[R_\mu+R_{\tau\mu}]$.
4. Hence $T:=24[R_\tau]+2M=24[R_V]$ for **every** Klein four-group $V\ni\tau$, where $R_V=\sum_{1\ne v\in V}R_v$. Nothing
   is divided, so torsion in $\mathrm{Pic}$ is harmless.
5. $N_G(V)$ permutes the involutions of $V$, so it preserves the divisor $R_V$. Hence $T$ is fixed by
   $\langle N(V_0),N(V_1)\rangle$.
6. Checked: $|N(V_0)|=72$, $|N(V_1)|=24$, and together they generate all 2520 elements. So $T$ is
   $A_7$-invariant, of degree $24\cdot54=1296$.
7. An order-5 subgroup acts freely, since there is no order-5 inertia. An invariant class for a freely
   acting cyclic group can be linearised ($H^2(C_n,\mathbb C^\*)=0$) and descended. So $5\mid1296$, which is false.

Degrees $\le8$ are excluded by (2.4). This uses the same determinant with $18>2\cdot8$ fixed points, and the fact
that $\bar z,\bar\sigma,\overline{\gamma\sigma}$ generate $H$ (checked).

**Assessment.** The proof is correct, short, and purely algebraic. The finite inputs are small enough to check
by hand, and we also checked them by computer.

**Relation to our audit theorem (NOTES §5.1).** That theorem says no argument using only the $H$-equivariant data
of $D$ can prove $\operatorname{gon}\ge10$. This proof uses the full $A_7$-action on $C$: normalisers of Klein four-groups
that do not centralise $\tau$, and a freely acting $C_5$. So it lies exactly outside the audit barrier, which is
where the audit said a proof would have to be. Our spectral proof also lies outside it, but for a different
reason: it uses the hyperbolic metric.

**Comparison with our result.**

| | this proof | our proof (§7) |
|---|---|---|
| bound on $\operatorname{gon}(C/\langle\tau\rangle)$ | 10 | 12 |
| method | algebraic | computer-assisted |
| verification | fully human-verifiable | trusted floating-point model |

The two proofs are logically independent. Each confirms the other's conclusion at the level $\ge10$.

## 2. The transport language (Thms 3.1, 3.2; §3.3–3.4)

- **Theorem 3.1.** $G$-equivariance comes from the centralizer sum. Connectivity of the Klein incidence graph
  makes all $T_\tau$ equal. Degree $3cr$ follows. For $A_7$ the graph was checked: 105 vertices, degree 8, 420
  edges, 140 triangles, connected. This gives a second proof of Theorem 2.3.
- **Theorem 3.2** is correct as stated. The hypotheses $m\equiv2\pmod4$ and saturation $\deg R_{\mu_i}+\deg R_{\nu_i}=2m$ are exactly
  what the argument uses.
- **§3.3 ($A_5$, signature $(2,2,2,3)$, genus 6).** Checked by computer: the branch cycles multiply to 1 and
  generate $A_5$. Each involution fixes 6 points, the quotients have genus 2 and hence degree-3 maps, and
  $5\nmid72$. So the theorem genuinely fails without connectivity. The component-sum repair (degree $3crb$) is
  correct.
- **§3.4.** Correctly presented as necessary identities only. Beyond saturation no obstruction is claimed.

## 3. The accessory barrier $\mathrm{ed}_{\mathbb C}(A_7;\le17)>1$ (Prop 4.2)

**Framework.** Farb–Wolfson (1.1) quantify over faithful $G$-varieties $\tilde X$ with $\tilde X\dashrightarrow X$ of degree $\le n$. Put
$K'=k(\tilde X)^G$. Then $k(X)\cap K'\subseteq k(X)^G=K$, so $k(\tilde X)=k(X)K'$. Farb–Wolfson covers are therefore
exactly the synthesis's accessories $F_{acc}=K'$, of degree $d$.

**The case split.**
- **$d$ prime to 2 or to 3.** Compression to a curve would give $\mathrm{ed}_p(A_7)\le1$. But $\mathrm{ed}_2(A_7)=\mathrm{ed}(D_8)=2$ and
  $\mathrm{ed}_3(A_7)=\mathrm{ed}(C_3^2)=2$ (Karpenko–Merkurjev).
- **$d\in\{6,12\}$.** Connectivity holds because $A_7$ has no subgroup of index 2, 3, 4, 6 or 12 (no subgroup of
  order 210). Then gonality $\ge13$ (Farb–Wolfson Lemma 2.2 with genus $\ge136$) and the transfer lemma with
  $M=A_5$, $|M|=60$, finish.

**Consistency with Farb–Wolfson Remark 1.9.** That remark says prime-local methods *alone* stop at 1. This
proof mixes them with the geometric argument degree by degree, so there is no conflict.

**A side remark on that paper.** Remark 1.9 states the Sylow-2 subgroup of $A_7$ is $\mathbb Z/2\times\mathbb Z/2$. It is
$D_8$, of order 8. This does not affect their theorem, and the synthesis uses the correct $D_8$.

**Assessment.** Correct, and apparently new: our novelty search found no improvement on Farb–Wolfson's
$\le6$. The Ladder rightly avoids a priority claim.

## 4. Pencil geometry and Theorem 6.1

**Verified:**
- **Lemma 5.1.** The orbit field is faithful, and an index $\ge2$ would give genus $\le64$. The $B_{18}(e)$ table
  is $145, 109, 73, 82$. For $m\le17$ all pairs are birational.
- **Lemma 5.2.** Values $108$, $127$, $192$, $217$ for $m=13,14,17,18$, and $\pi_0(51,7)=184$, $\pi_0(54,7)=208$.
- **Lemma 5.4.** The divisor degree on the Segre threefold is $2(a+b+c)$; the case analysis is correct.
- **Lemma 5.3** (budget). The reasoning holds: pairwise-intersecting centre pairs form a star or a triangle,
  and multiplicities 0 or 1 are unaffordable. Spot-checked entries:
  - $(15,136)$: $\delta=60$, at most 1 third;
  - $(16,136)$: $\delta=89$, at most 3 thirds;
  - $(16,169)$: $\delta=56$, at most 1;
  - $(17,169)$: $\delta=87$, at most 2;
  - $(18,169)$: $\delta=120$, star $15;3^5$, at most 5;
  - $(18,199)$ and $(18,211)$: at most 1.
- **Table 1** is complete for genus $\le289$. Every non-listed triangle fails the seven-sheet test
  (e.g. $(2,5,6)$ and $(2,6,6)$). Every listed signature is realised by generating triples of $A_7$ (counts in
  `review_checks.py`).

**Dependency of Theorem 6.1, row by row:**

| genus | signatures | what the exclusion of degree $\le18$ uses |
|---|---|---|
| 136 | $(2,4,7)$ | Lemma 5.5 at $m=16$, budget. **Superseded:** our $\operatorname{gon}(C)\ge24$. |
| 169 | $(3,3,5)$ | orbit field ($m\le13$); Lemma 5.2 ($m=14,15$); $\pi_0(48,7)=161$ + budget ($m=16$); Eisenbud–Harris $\pi_1(51,7)=161$ + Lemma 5.4 + budget ($m=17$); **Lemma 5.5** + budget ($m=18$). |
| 199 | $(2,5,7)$ | orbit field ($m\le15$); Lemma 5.2 ($m=16,17$); $\pi_1(54,7)=182$ (Eisenbud–Harris) + Lemma 5.4 + budget ($m=18$). |
| 211 | $(3,3,6)$, $(3,4,4)$ | orbit field ($m\le15$); Lemma 5.2 ($m=16,17$); $\pi_0(54,7)=208$ + budget ($m=18$). Elementary. |
| 241, 271, 274 | $(2,6,7)$, $(3,3,7)$, $(2,7,7)$, $(3,4,5)$ | orbit field + Lemma 5.2 only. Elementary. |

So the single unverified input is Lemma 5.5 (Petrakiev ranges), needed only for $(3,3,5)$ at $m=18$. The Ladder
itself asks for a referee pass there, correctly.

## 5. The torsion/norm model (§7)

This is independently consistent with our first-session analysis:
- Genera $E=22$, $Q=28$, $F=10$, $B=T=3$, $Y=4$, $X=0$.
- The trigonal map $B\to X$ has special fibres $P_1+P_2+P_\tau$ and $2P_3+P_t$. We had $t_0+t_1+t_2$ and $t_b+2t_3$.
- The derived condition $3\Delta_z\sim\Delta_\sigma$ is the analogue of our $(\star)$ ($\sum g_j\sim3\sum t_i$, NOTES 4.5).

The moving-case correction (three divisors $W_\pm,W_0$) and the root-versus-square cautions are sound. Nothing
here is needed for the theorems, as the synthesis says.

## 6. Minor issues and suggestions

1. **Lemma 2.1.** As stated it covers only $\mathrm{Fix}(v)$. The application (2.5) uses the two-involution form:
   zeros on $\mathrm{Fix}(\mu)\cup\mathrm{Fix}(\tau\mu)$ for $\tau$-invariant $f$. It is worth stating that form explicitly.
2. **Orbit sizes ("at least eight").** Orbits have size $[A_7:\mathrm{Stab}]$, and after 7 the next subgroup index is
   15. So a minimal pencil's orbit has $\ge15$ members, giving at least 13 third pencils rather than 6. This
   strengthens every budget argument in Thm 6.1 at no cost.
3. **Table 1 and marked covers.** State explicitly that each row has several marked classes (our counts:
   18, 96, 18, 72, 96, 90, 144, 90 generating pairs per fixed $x$-class). Any certificate must cover every
   class.
4. **Lemma 5.5.** Needs the Petrakiev ranges checked line by line, as the Ladder notes. Alternatively bypass it
   (§7 below).
5. **Ladder §4** ("leave only … degrees 17, 18"). This is now resolved: our Theorem 7.1 gives
   $\operatorname{gon}(C)\ge24$, so $W^1_{17}(C)=W^1_{18}(C)=\emptyset$. **Rungs G17 and G18 are proved.**

## 7. What the two bodies of work give together

1. **First summit, proved twice independently:** algebraically, $\ge10$ (Sol/Astra), and spectrally, $\ge12$
   (ours).
2. **Rungs G17 and G18 are closed** by our $\operatorname{gon}(C)\ge24$.
3. **Accessory barrier to degree 23 (conditional).** A degree-18 accessory is connected (their index
   argument). Farb–Wolfson transfer with $A_5$ then requires a faithful $A_7$-curve of gonality $\le18$. By Lemma
   5.1 such a curve has genus $\le289$, so it lies in Table 1. Genus 136 is excluded by our theorem, the other
   rows by their Theorem 6.1. For $d=19,\dots,23$, each $d$ is prime to 2 or to 3, and Karpenko–Merkurjev applies.
   Hence
   $$\mathrm{ed}_{\mathbb C}(A_7;\le23)>1,$$
   conditional only on Lemma 5.5 for the $(3,3,5)$ row (and on our certificate). This is the natural limit
   of the combined method. Degree 24 would need $\operatorname{gon}\ge25$ in genus 136, and Hersch/Yang–Yau cannot give
   that for classes 0, 1, since $67.5\cdot0.3463<24$.
4. **Removing the condition.** Run our certified spectral method on the $(3,3,5)$ curves, or on all eight
   Table-1 signatures to bypass §5–6 entirely. The required gaps are small: $\lambda_1>36/(g-1)$, i.e.
   $0.214$ for genus 169 and $0.13$ for genus 274. The code needs only a general $(p,q,r)$ reference triangle
   and a cover of irreducibles per signature.
5. **What remains open,** as both documents correctly say:
   - $\mathrm{RD}(A_7)>1$ (tower-stable obstructions, T1 and T2);
   - the accessory barrier beyond degree 23.

## 8. Overall assessment

**Rigour.** The synthesis is careful. It states quantifiers and marks statuses honestly. It withdraws the
earlier moving-case claim and documents its failed routes.

**Main theorem.** Theorem 2.3 is correct. It is the cleanest known proof of the first summit and needs no
computer.

**New tools.** The transport and amalgamation theorems are genuinely new, reusable tools: local averaging over
a centraliser, then normaliser amplification, then descent under a freely acting cyclic group. The
connectivity counterexample shows the authors tested the tool's limits.

**Barrier result.** $\mathrm{ed}_{\mathbb C}(A_7;\le17)>1$ is correct and extends the published bound.

**Remaining risk.** It sits entirely in the refined-Castelnuovo chain (Lemma 5.5), which the authors themselves
flag.

Combined with the spectral certificate, the two projects give a coherent package:
- the first summit, proved twice independently;
- both geometric rungs G17 and G18;
- one unverified lemma away from $\mathrm{ed}_{\mathbb C}(A_7;\le23)>1$, a step our spectral method could also bypass.
