# Independent audit and continuation, 4 October

Scope:
- every numbered statement of Chapters 1–7;
- the GPT material they rely on: `gpt/audit_3oct/MULTIPLE_PENCILS.md`, `gpt/A7_Normalization_Defects.md`, `gpt/A7_Quadratic_Systems.md` and the Amitsur note;
- the scripts behind them.

Each proof was re-derived by hand, and each computation was rerun and compared with its saved output. The new mathematics of this round is in §2.5 (classes 12, 14) and §7.13.

## 1. Verdict

No error affecting a stated result was found. The tags are accurate: [P] proofs are complete given their cited inputs, and [N] items are labelled as such. Six expository gaps were closed (§3). The trust base is unchanged (`README.md`, Trust base), plus one new classical input in §7.13: the Douady–Earle conformal barycentre.

## 2. Reproduction

| script | result |
|---|---|
| `cover`, `curve_checks`, `verify_exact`, `gonality_large_genus`, `verify_accessory60`, `equivariant_rr`, `model99`, `audit_exact`, `verify_mu90_exact`, `side_checks`, `elliptic_subcovers`, `gpt/verify_normalization`, `gpt/verify_quadrics`, `gpt/textbook/verify_textbook_inputs` | byte-identical to the saved outputs |
| `twisted_rr`, `twisted_survey`, `tau_pencil` | identical apart from timing lines |
| `run_certificate.py` (Ch. 2) | identical certified values: $Q_0\ge0.997833$, $Q_1\ge0.340893$, $Q_2\ge0.689027$ (classes 0, 1), $\ge0.762988$ (12, 14) |
| `hh_certify.py` classes 0, 1 at $192\times128$ (Ch. 3) | identical: left $0.330681>0.271290$ (class 0), $0.329591>0.279466$ (class 1); the same $a$, $\rho$, $\Gamma_h$, $B_h$, $\varphi_\sigma$ |
| `certify_th.py` count, `s5`, `a6`, classes 0, 1 | identical (one negative pivot; $U=0.36318$; $\ge0.99783$) |
| `hejhal_solve.py 12 90` | $\lambda=0.3596713548916$, against the recorded $0.359671354895$ |

Not rerun, for time:
- the Lean build (toolchain not installed in this session);
- `certify_signatures.py` (window), `validate_bolza.py`, `p5_curve.py` [N], `quotient_gonality.py` [N].

## 3. Gaps closed

1. **§6.1, step 2 (averaging).** $T=24[R_V]$ for *every* $V\ni\tau$ needs step 1 applied to each $\tau$-invariant pencil $f\circ h$, $h\in C_G(\tau)$. Applied to $f$ alone, the average is $8\sum_jR_{V_1^{(j)}}$ for the second Klein type. The text now says so.
2. **Theorem 4.9, step 1.** The bound $d^2/e+ed-3d+1$ covers chains of any length only through the maximal-compositum argument of `MULTIPLE_PENCILS.md` §4. This is now explicit.
3. **Theorem 4.8, step 2.** The estimate is $h_B(2)\le h^0(\mathcal O_B(1))+h_\Gamma(2)$, which is where linear normality enters. This is now explicit.
4. **Proposition 7.15.** The eigenline stabilisers are now proved by hand:
   - the subgroups containing a 7-cycle are $C_7$, $7{:}3$, $L_2(7)$ (two classes) and $A_7$;
   - $7{:}3$ moves every eigenline of $\hat c$;
   - $\mathbf 6$ is irreducible on $L_2(7)$ and on $A_7$.

   The exact check at $p=43$ is now confirmatory.
5. **Lemma 3.5.** The sign in $\int w\{u,v\}=\pm\int\nabla w\cdot J$ depends on the $\nabla^\perp$ convention. Only the bound is used. No change.
6. **Checked, no change.**
   - The Crouzeix–Raviart constant $\kappa_e^2=h^2(1/16+1/j_{1,1}^2)$ exceeds Carstensen–Gedicke's $h^2(1/48+1/j_{1,1}^2)$, which is the safe direction.
   - Petrakiev's Theorem 2.16(a) (arXiv:math/0604517) needs $n\ge8$ and $d\ge2n+3=19$, and its $\pi_2$ is the one used. The differences $A(m)-\pi_2(3m,8)=s/2,(s+1)/2,s/2+1$ and the six deficit polynomials hold for all $8\le m<2000$ (independent recomputation).

## 4. New results (§2.5, §7.13)

**Lemma 7.17 [P].** Free cyclic symmetries of a gonal pencil's class divide its degree: $3\mid d$, $9\mid d$ or $5\mid d$ according to $K$.

**Theorem 7.18 [P].** Equivariant Li–Yau.
- Douady–Earle balancing is conformally natural, so the balanced map is $K$-equivariant. Its energy lives on the isotypes $\sigma$ with $\mathrm{Hom}_K(R,\sigma)\ne0$.
- Multiplicities for all 40 subgroup classes: `equivariant_spectral.py` [X].
- Consequence: an $A_5$-symmetric gonal pencil has $d\ge45$, and the $(A_4\times3){:}2$ stratum with kernel $C_3$ has $d\ge54$.

**Proposition 7.20 [X].** Frobenius at $p=37,43,67$ on both curves (`frobenius_lpolys.py`):
- $A$ and $S$ have no elliptic isogeny factor;
- $E_{15}$ and $E_{21}$ have no CM and are not isogenous.

The point counts use the septic Belyi model and the local Galois action at the branch points. The checks are the genus-0 quotients, the Weil bounds, and $L_A=L(10)L(\overline{10})$ over $\mathbb Q(\sqrt{-7})$.

**Theorem 7.21 [P][X].** Every elliptic subcover of $C$ has degree $\ge60$. Proposition 7.10 now covers all elliptic subcovers, not only pure planes.

**Theorem 7.22 [P][X][C].** No gonal pencil of degree $\le41$ has kernel $C_7$. This settles half of open problem 3(a) of `HANDOFF.md`, without the $K_{3,2}$ computation it proposed.

**Corollary 7.23.**
1. The kernel list of Corollary 7.14 becomes $1,C_2,C_3,C_4,V_4,S_3$ [P][N].
2. Without any numerical input [P][C], the list is $1,C_2,C_3,C_4,V_4,C_5,C_6,S_3,D_8,D_{10},D_{12},5{:}4$, with the admissible degrees (`equivariant_spectral.py --strata`).

**Classes 12, 14 (§2.5).**
- [C] $14_a$ has an eigenvalue in $[0.355695,0.37722]$, and $Q_1$ has at most two eigenvalues below $0.55998$.
- [N] The second lies on the $21$-isotype near $0.38$.
- [P] Because $\langle21,\wedge^214_a\rangle=1$, the harmonic Hersch route to $\operatorname{gon}\ge26$ is closed on that curve.

## 5. Unchanged

- Exact gonality in $[25,42]$.
- The kernels $S_3$ and $C_2$, and the faithful stratum.
- $\operatorname{gon}\ge26$ on classes 0, 1.
- Towers.
- $P_E$, $P_S$.

These remain as in `HANDOFF.md`.
