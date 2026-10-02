# ARAGACAS
------
:rock_boom:

## Hilbert 13 and the $A_7$ gonality ladder

This is research on Hilbert's 13th problem through curves with an $A_7$-action, written as a short textbook with machine-checked certificates. **New readers, human or AI, should start with [`GUIDE.md`](GUIDE.md).**

**Main results.** $C$ is a $(2,4,7)$ $A_7$-curve, of genus 136 (the minimum), and $\tau$ an involution.
- $25\le\operatorname{gon}(C)\le42$, and $13\le\operatorname{gon}(C/\tau)\le21$.
  - The lower bound combines a certified spectral gap with a cubic-form refinement of Hersch's argument.
  - The upper bound comes from a degree-60 model in $\mathbb P^5$ on the Laza–Zheng $A_7$-cubic fourfold.
- $25\le\gamma(A_7)\le42$, where $\gamma(A_7)$ is the least gonality of any faithful $A_7$-curve.
- **Accessory irrationalities.**
  - $\mathrm{ed}_{\mathbb C}(A_7;\le59)>1$, sharp: $a(A_7)=60$. The published bound was 6.
  - With connected full monodromy the threshold is exactly $\mu(A_7)=90$.

**Map.**

| path | contents |
|---|---|
| [`GUIDE.md`](GUIDE.md) | how to read the repository: order, notation, dependencies, verification, pitfalls |
| [`HANDOFF.md`](HANDOFF.md) | current status, open problems, working rules |
| [`hilbert13/ladder/`](hilbert13/ladder/README.md) | seven chapters, scripts with saved outputs, the literature, GPT's documents |
| [`hilbert13/`](hilbert13/README.md) | Lean 4 / Mathlib: single superpositions, and the core of the spectral certificate |

Nothing here is refereed. Every claim is tagged by how it is established: paper proof, exact computation, computer-assisted proof, Lean, or numerical.
