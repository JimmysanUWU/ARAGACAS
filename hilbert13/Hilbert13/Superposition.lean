/-
Copyright (c) 2026 ARAGACAS contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: ARAGACAS contributors
-/
import Mathlib

/-!
# Single superpositions and Hilbert's 13th problem

Hilbert's 13th problem asks whether functions of several variables can be built from
functions of fewer variables by composition. This file records why the choice of function
class is the whole problem:

* `exists_single_superposition`: with **no regularity**, every `f : ℝⁿ → ℝ` is
  `g (φ₁ x₁ + ⋯ + φₙ xₙ)`, via a `ℚ`-linear isomorphism `ℝⁿ ≃ ℝ` (Hamel basis).
* `mul_ne_single_superposition`: with **continuous inner functions** this fails already for
  `x * y`, even if the outer function `g` is arbitrary.
* `mul_two_superposition`: two continuous terms suffice for `x * y`, consistent with the
  Kolmogorov–Arnold theorem, which needs several outer terms.
-/

open Cardinal

/-- `#(Fin n → ℝ) = 𝔠` for `n ≥ 1`. -/
lemma mk_fin_real {n : ℕ} (hn : 1 ≤ n) : #(Fin n → ℝ) = 𝔠 := by
  rw [← Cardinal.power_nat_eq (c := 𝔠) Cardinal.aleph0_le_continuum hn]
  simp [Cardinal.mk_real]

lemma rank_rat_pi {n : ℕ} (hn : 1 ≤ n) : Module.rank ℚ (Fin n → ℝ) = 𝔠 := by
  rw [Module.Free.rank_eq_mk_of_infinite_lt ℚ (Fin n → ℝ), mk_fin_real hn]
  simp [mk_fin_real hn, Cardinal.aleph0_lt_continuum]

lemma rank_rat_real : Module.rank ℚ ℝ = 𝔠 := by
  rw [Module.Free.rank_eq_mk_of_infinite_lt ℚ ℝ, Cardinal.mk_real]
  simp [Cardinal.mk_real, Cardinal.aleph0_lt_continuum]

/-- A `ℚ`-linear isomorphism `ℝⁿ ≃ ℝ` (Hamel bases). -/
noncomputable def hamelEquiv {n : ℕ} (hn : 1 ≤ n) : (Fin n → ℝ) ≃ₗ[ℚ] ℝ :=
  (nonempty_linearEquiv_of_rank_eq ((rank_rat_pi hn).trans rank_rat_real.symm)).some

/-- **Without continuity, Hilbert 13 is trivial.** Every function of `n ≥ 1` real
variables is `g (φ₁ x₁ + ⋯ + φₙ xₙ)` for one-variable functions `g, φᵢ`. -/
theorem exists_single_superposition {n : ℕ} (hn : 1 ≤ n) (f : (Fin n → ℝ) → ℝ) :
    ∃ (g : ℝ → ℝ) (φ : Fin n → ℝ → ℝ), ∀ x, f x = g (∑ i, φ i (x i)) := by
  let e := hamelEquiv hn
  refine ⟨fun t => f (e.symm t), fun i t => e (Pi.single i t), fun x => ?_⟩
  simp only [← map_sum, Finset.univ_sum_single, LinearEquiv.symm_apply_apply]

/-- Two-variable form of `exists_single_superposition`. -/
theorem exists_single_superposition₂ (f : ℝ → ℝ → ℝ) :
    ∃ g φ ψ : ℝ → ℝ, ∀ x y, f x y = g (φ x + ψ y) := by
  obtain ⟨g, φ, h⟩ := exists_single_superposition (n := 2) (by norm_num)
    (fun v => f (v 0) (v 1))
  refine ⟨g, φ 0, φ 1, fun x y => ?_⟩
  simpa [Fin.sum_univ_two] using h ![x, y]

/-- **With continuous inner functions it fails**, even for `x * y` and even allowing an
arbitrary (discontinuous) outer function `g`. -/
theorem mul_ne_single_superposition :
    ¬ ∃ g φ ψ : ℝ → ℝ, Continuous φ ∧ Continuous ψ ∧ ∀ x y, x * y = g (φ x + ψ y) := by
  rintro ⟨g, φ, ψ, hφ, hψ, h⟩
  -- `ψ` is injective: `y = g (φ 1 + ψ y)`.
  have hinj : Function.Injective ψ := fun a b hab => by
    simpa [hab] using (h 1 a).trans ((congrArg g (by rw [hab])).trans (h 1 b).symm)
  -- so `ψ (1/2)` lies strictly between `ψ 0` and `ψ 1`.
  have hmid : ψ (1 / 2) ∈ Set.Ioo (min (ψ 0) (ψ 1)) (max (ψ 0) (ψ 1)) := by
    rcases hψ.strictMono_of_inj hinj with hm | ha
    · have h1 := hm (show (0 : ℝ) < 1 / 2 by norm_num)
      have h2 := hm (show (1 / 2 : ℝ) < 1 by norm_num)
      constructor
      · exact min_lt_of_left_lt h1
      · exact lt_max_of_lt_right h2
    · have h1 := ha (show (0 : ℝ) < 1 / 2 by norm_num)
      have h2 := ha (show (1 / 2 : ℝ) < 1 by norm_num)
      constructor
      · exact min_lt_of_right_lt h2
      · exact lt_max_of_lt_left h1
  -- Perturb: for some `t ≠ 0`, `φ t - φ 0 + ψ (1/2)` is still in that interval.
  have hcont : Continuous fun t => φ t - φ 0 + ψ (1 / 2) := by fun_prop
  have hev : ∀ᶠ t in nhdsWithin (0 : ℝ) {0}ᶜ,
      φ t - φ 0 + ψ (1 / 2) ∈ Set.Ioo (min (ψ 0) (ψ 1)) (max (ψ 0) (ψ 1)) := by
    apply eventually_nhdsWithin_of_eventually_nhds
    apply hcont.continuousAt.preimage_mem_nhds
    simpa using isOpen_Ioo.mem_nhds hmid
  obtain ⟨t, ht, ht0⟩ := (hev.and self_mem_nhdsWithin).exists
  -- By the IVT it equals `ψ y'` for some `y'`.
  obtain ⟨y', -, hy'⟩ := intermediate_value_uIcc (a := (0 : ℝ)) (b := 1) hψ.continuousOn
    (Set.Ioo_subset_Icc_self ht)
  -- Then `t / 2 = g (φ t + ψ (1/2)) = g (φ 0 + ψ y') = 0`.
  have key : φ t + ψ (1 / 2) = φ 0 + ψ y' := by linarith
  have := (h t (1 / 2)).trans (key ▸ (h 0 y').symm)
  simp only [zero_mul, mul_eq_zero] at this
  rcases this with h0 | h0
  · exact ht0 h0
  · norm_num at h0

/-- Contrast: two continuous terms suffice for `x * y`. -/
theorem mul_two_superposition (x y : ℝ) :
    x * y = (fun t => t ^ 2 / 4) (x + y) + (fun t => -(t ^ 2) / 4) (x + -y) := by
  ring
