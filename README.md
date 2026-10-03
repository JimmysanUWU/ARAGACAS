# Hilbert 13 and the $A_7$ gonality ladder

Hilbert's 13th problem through curves with an $A_7$-action. The [textbook edition](hilbert13/ladder/gpt/textbook/README.md) organizes the whole repository into definitions, statements, proofs and verification records. Start with [`GUIDE.md`](GUIDE.md) for the research summary. Lean checks specified logical cores; the finite and analytical inputs have separate certificates.

**Results.** $C$ is a $(2,4,7)$ $A_7$-curve (genus 136, the minimum) and $\tau$ an involution.
- $25\le\operatorname{gon}(C)\le42$ and $13\le\operatorname{gon}(C/\tau)\le21$.
  - Lower: a certified spectral gap and a cubic-form refinement of Hersch's argument.
  - Upper: a degree-60 Schur-twisted line bundle, certified exactly. Its sections map $C$ to the Laza–Zheng $A_7$-cubic fourfold in $\mathbb P^5$.
- $25\le\gamma(A_7)\le42$, where $\gamma(A_7)$ is the least gonality of a faithful $A_7$-curve.
- $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$, sharp: $a(A_7)=60$ (the published bound was 6). With connected full monodromy the threshold is $\mu(A_7)=90$.

| path | contents |
|---|---|
| [`GUIDE.md`](GUIDE.md) | the mathematics in one page; reading order, notation, dependencies, pitfalls |
| [`HANDOFF.md`](HANDOFF.md) | status, open problems, working rules |
| [`textbook/`](hilbert13/ladder/gpt/textbook/README.md) | systematic edition of all mathematical sources, with explicit hypotheses and evidence status |
| [`hilbert13/ladder/`](hilbert13/ladder/README.md) | seven chapters, scripts with saved outputs, literature, GPT's documents |
| [`hilbert13/`](hilbert13/README.md) | Lean 4: single superpositions, and the core of the spectral certificate |

Nothing here is refereed. Every claim is tagged by how it is established: paper proof, exact computation, computer-assisted proof, Lean, or numerical.
