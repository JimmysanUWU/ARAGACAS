/-
Copyright (c) 2026 ARAGACAS contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: ARAGACAS contributors
-/
import Mathlib

/-!
# The logical core of the spectral certificate (ladder/2_SPECTRAL.md, §2.4)

The proof of `gon(C) ≥ 23` and `gon(C/⟨τ⟩) ≥ 12` for the (2,4,7) `A₇`-curves combines analysis
(Hersch's inequality), a finite element computation and floating-point linear algebra.  This file
formalizes the three purely logical steps that turn the computed numbers into the theorem:

1. `cr_lower_bound` — the abstract guaranteed lower bound for the lowest eigenvalue
   (Liu 2015; Carstensen–Gedicke 2014), used in NOTES 7.6: an `a`-orthogonal projection `P` onto a
   discrete set `Vh` with `‖u - P u‖² ≤ C² a(u - P u, u - P u)` and a discrete lower bound
   `λ_h ‖v‖² ≤ a(v, v)` on `Vh` give `λ_h / (1 + C² λ_h) ‖u‖² ≤ a(u, u)` for every `u`.
   No eigenfunction, compactness or min-max theory is needed for the lowest eigenvalue.
2. `posDef_of_shifted_factor(_perm)` — NOTES 7.7: if `L Lᵀ = B - c I + Δ` with the componentwise
   backward error bound `|Δᵢⱼ| ≤ γ ∑ₖ |Lᵢₖ| |Lⱼₖ|` of floating-point Cholesky, row sums of `|L|`
   at most `R`, column sums at most `Cc`, and `γ R Cc < c`, then `B` is positive definite
   (also after the symmetric permutation chosen by the sparse factorization).
3. `gon_C_ge_23`, `gon_C_ge_24`, `gon_D_ge_12` — the arithmetic of NOTES 7.2 from the certified
   `λ₁(C) ≥ 0.33335` (resp. `≥ 0.34089`), and the Riemann–Hurwitz numbers of rung 1.
-/

namespace Hilbert13.Spectral

section CR

variable {W : Type*} [NormedAddCommGroup W] [Module ℝ W]

/-- **Guaranteed lower bound** (abstract form of Liu's framework, lowest eigenvalue).
`W` carries the mass norm, `a` is the (broken) energy form, `S` the continuous space, `Vh` the
discrete space, `P` the Crouzeix–Raviart interpolation. -/
theorem cr_lower_bound (a : LinearMap.BilinForm ℝ W) (ha_symm : ∀ x y, a x y = a y x)
    (ha_nonneg : ∀ x, 0 ≤ a x x) (S Vh : Set W) (P : W → W) (hP : ∀ u ∈ S, P u ∈ Vh)
    (horth : ∀ u ∈ S, ∀ v ∈ Vh, a (u - P u) v = 0) {C lh : ℝ} (hC0 : 0 ≤ C)
    (hC : ∀ u ∈ S, ‖u - P u‖ ^ 2 ≤ C ^ 2 * a (u - P u) (u - P u))
    (hlh : 0 ≤ lh) (hmin : ∀ v ∈ Vh, lh * ‖v‖ ^ 2 ≤ a v v) {u : W} (hu : u ∈ S) :
    lh / (1 + C ^ 2 * lh) * ‖u‖ ^ 2 ≤ a u u := by
  have hvVh : P u ∈ Vh := hP u hu
  -- the energy splits orthogonally: a(u,u) = a(Pu,Pu) + a(e,e), e = u - Pu
  have hsplit : a u u = a (P u) (P u) + a (u - P u) (u - P u) := by
    have h1 : a (u - P u) (P u) = 0 := horth u hu (P u) hvVh
    have h2 : a (P u) (u - P u) = 0 := by rw [ha_symm]; exact h1
    have huv : u = P u + (u - P u) := by abel
    conv_lhs => rw [huv]
    simp only [map_add, LinearMap.add_apply, h1, h2]
    ring
  set t := Real.sqrt (a (u - P u) (u - P u)) with ht_def
  have ht0 : 0 ≤ t := Real.sqrt_nonneg _
  have ht : t ^ 2 = a (u - P u) (u - P u) := Real.sq_sqrt (ha_nonneg _)
  have hCt : 0 ≤ C * t := mul_nonneg hC0 ht0
  have hCe : ‖u - P u‖ ≤ C * t := by
    have h := hC u hu
    rw [← ht, ← mul_pow] at h
    exact (pow_le_pow_iff_left₀ (norm_nonneg _) hCt two_ne_zero).mp h
  have htri : ‖u‖ - ‖u - P u‖ ≤ ‖P u‖ := by
    have h := norm_sub_norm_le u (u - P u)
    rwa [sub_sub_cancel] at h
  have hmv := hmin (P u) hvVh
  have hden : 0 < 1 + C ^ 2 * lh := by positivity
  have hav : 0 ≤ a (P u) (P u) := ha_nonneg _
  rw [div_mul_eq_mul_div, div_le_iff₀ hden, hsplit, ← ht]
  by_cases hcase : C * t ≤ ‖u‖
  · have h1 : ‖u‖ - C * t ≤ ‖P u‖ := by linarith
    have h2 : (‖u‖ - C * t) ^ 2 ≤ ‖P u‖ ^ 2 := pow_le_pow_left₀ (by linarith) h1 2
    have h3 : lh * (‖u‖ - C * t) ^ 2 ≤ a (P u) (P u) :=
      le_trans (mul_le_mul_of_nonneg_left h2 hlh) hmv
    have h4 := mul_le_mul_of_nonneg_right h3 hden.le
    nlinarith [sq_nonneg ((1 + C ^ 2 * lh) * t - lh * C * ‖u‖)]
  · push Not at hcase
    have hu2 : ‖u‖ ^ 2 ≤ (C * t) ^ 2 := pow_le_pow_left₀ (norm_nonneg _) hcase.le 2
    have h5 := mul_le_mul_of_nonneg_left hu2 hlh
    have h6 : 0 ≤ a (P u) (P u) * (1 + C ^ 2 * lh) := mul_nonneg hav hden.le
    nlinarith [sq_nonneg t]

end CR

section Cholesky

open Matrix

/-- weighted Cauchy–Schwarz: `(∑ wᵢ yᵢ)² ≤ (∑ wᵢ)(∑ wᵢ yᵢ²)` for `w ≥ 0`. -/
lemma weighted_cs {ι : Type*} (s : Finset ι) (w y : ι → ℝ) (hw : ∀ i, 0 ≤ w i) :
    (∑ i ∈ s, w i * y i) ^ 2 ≤ (∑ i ∈ s, w i) * ∑ i ∈ s, w i * y i ^ 2 := by
  have h := Finset.sum_mul_sq_le_sq_mul_sq s (fun i => Real.sqrt (w i))
    (fun i => Real.sqrt (w i) * y i)
  have e1 : ∀ i, Real.sqrt (w i) * (Real.sqrt (w i) * y i) = w i * y i := fun i => by
    rw [← mul_assoc, Real.mul_self_sqrt (hw i)]
  have e2 : ∀ i, Real.sqrt (w i) ^ 2 = w i := fun i => Real.sq_sqrt (hw i)
  have e3 : ∀ i, (Real.sqrt (w i) * y i) ^ 2 = w i * y i ^ 2 := fun i => by
    rw [mul_pow, e2]
  simp only [e1, e2, e3] at h
  exact h

/-- the Schur-test bound `|x|ᵀ |L| |L|ᵀ |x| ≤ R Cc ‖x‖²`. -/
lemma schur_bound {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) {R Cc : ℝ}
    (hrow : ∀ i, ∑ k, |L i k| ≤ R) (hcol : ∀ k, ∑ i, |L i k| ≤ Cc) (x : Fin n → ℝ) :
    ∑ i, ∑ j, |x i| * |x j| * ∑ k, |L i k| * |L j k| ≤ R * Cc * ∑ i, x i ^ 2 := by
  have hCc : ∀ k, 0 ≤ Cc := fun k => le_trans (Finset.sum_nonneg fun i _ => abs_nonneg _) (hcol k)
  -- rewrite as ∑_k (∑_i |L i k| |x i|)^2
  have hre : ∑ i, ∑ j, |x i| * |x j| * ∑ k, |L i k| * |L j k|
      = ∑ k, (∑ i, |L i k| * |x i|) ^ 2 := by
    have h : ∀ k, (∑ i, |L i k| * |x i|) ^ 2
        = ∑ i, ∑ j, (|L i k| * |x i|) * (|L j k| * |x j|) := by
      intro k; rw [sq, Finset.sum_mul_sum]
    simp_rw [h, Finset.mul_sum]
    conv_rhs => rw [Finset.sum_comm]; enter [2, i]; rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ =>
      Finset.sum_congr rfl fun k _ => ?_
    ring
  rw [hre]
  calc ∑ k, (∑ i, |L i k| * |x i|) ^ 2
      ≤ ∑ k, (∑ i, |L i k|) * ∑ i, |L i k| * |x i| ^ 2 :=
        Finset.sum_le_sum fun k _ => weighted_cs _ _ _ fun i => abs_nonneg _
    _ ≤ ∑ k, Cc * ∑ i, |L i k| * |x i| ^ 2 := by
        refine Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_right (hcol k) ?_
        exact Finset.sum_nonneg fun i _ => mul_nonneg (abs_nonneg _) (sq_nonneg _)
    _ = Cc * ∑ i, x i ^ 2 * ∑ k, |L i k| := by
        rw [← Finset.mul_sum, Finset.sum_comm]
        congr 1
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [Finset.mul_sum]
        refine Finset.sum_congr rfl fun k _ => ?_
        rw [sq_abs]; ring
    _ ≤ Cc * ∑ i, x i ^ 2 * R := by
        rcases Nat.eq_zero_or_pos n with hn | hn
        · subst hn; simp
        · refine mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun i _ =>
            mul_le_mul_of_nonneg_left (hrow i) (sq_nonneg _)) (hCc ⟨0, hn⟩)
    _ = R * Cc * ∑ i, x i ^ 2 := by rw [← Finset.sum_mul]; ring

/-- **Verified positive definiteness from a perturbed Cholesky factor** (NOTES 7.7). -/
theorem posDef_of_shifted_factor {n : ℕ} (B L Δ : Matrix (Fin n) (Fin n) ℝ) {c γ R Cc : ℝ}
    (hfac : L * Lᵀ = B - c • (1 : Matrix (Fin n) (Fin n) ℝ) + Δ)
    (hΔ : ∀ i j, |Δ i j| ≤ γ * ∑ k, |L i k| * |L j k|) (hγ : 0 ≤ γ)
    (hrow : ∀ i, ∑ k, |L i k| ≤ R) (hcol : ∀ k, ∑ i, |L i k| ≤ Cc) (hc : γ * R * Cc < c)
    (x : Fin n → ℝ) (hx : x ≠ 0) : 0 < x ⬝ᵥ (B *ᵥ x) := by
  have hB : B = L * Lᵀ + c • (1 : Matrix (Fin n) (Fin n) ℝ) - Δ := by rw [hfac]; abel
  -- the three quadratic forms
  have hLL : 0 ≤ x ⬝ᵥ ((L * Lᵀ) *ᵥ x) := by
    rw [← mulVec_mulVec, dotProduct_mulVec, ← mulVec_transpose]
    simp only [dotProduct]
    exact Finset.sum_nonneg fun i _ => mul_self_nonneg _
  have hI : x ⬝ᵥ ((c • (1 : Matrix (Fin n) (Fin n) ℝ)) *ᵥ x) = c * ∑ i, x i ^ 2 := by
    rw [smul_mulVec, one_mulVec, dotProduct_smul, smul_eq_mul]
    simp [dotProduct, sq]
  have hΔx : |x ⬝ᵥ (Δ *ᵥ x)| ≤ γ * R * Cc * ∑ i, x i ^ 2 := by
    have e : x ⬝ᵥ (Δ *ᵥ x) = ∑ i, ∑ j, x i * Δ i j * x j := by
      simp [dotProduct, mulVec, Finset.mul_sum, mul_assoc]
    rw [e]
    calc |∑ i, ∑ j, x i * Δ i j * x j| ≤ ∑ i, ∑ j, |x i * Δ i j * x j| :=
          (Finset.abs_sum_le_sum_abs _ _).trans
            (Finset.sum_le_sum fun i _ => Finset.abs_sum_le_sum_abs _ _)
      _ ≤ ∑ i, ∑ j, |x i| * |x j| * (γ * ∑ k, |L i k| * |L j k|) := by
          refine Finset.sum_le_sum fun i _ => Finset.sum_le_sum fun j _ => ?_
          rw [abs_mul, abs_mul]
          have := hΔ i j
          nlinarith [abs_nonneg (x i), abs_nonneg (x j), abs_nonneg (Δ i j),
            mul_nonneg (abs_nonneg (x i)) (abs_nonneg (x j))]
      _ = γ * ∑ i, ∑ j, |x i| * |x j| * ∑ k, |L i k| * |L j k| := by
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun i _ => ?_
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun j _ => ?_
          ring
      _ ≤ γ * (R * Cc * ∑ i, x i ^ 2) :=
          mul_le_mul_of_nonneg_left (schur_bound L hrow hcol x) hγ
      _ = γ * R * Cc * ∑ i, x i ^ 2 := by ring
  have hxx : 0 < ∑ i, x i ^ 2 := by
    obtain ⟨i, hi⟩ := Function.ne_iff.mp hx
    rw [Pi.zero_apply] at hi
    exact lt_of_lt_of_le (lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hi)))
      (Finset.single_le_sum (fun j _ => sq_nonneg (x j)) (Finset.mem_univ i))
  rw [hB, sub_mulVec, add_mulVec, dotProduct_sub, dotProduct_add, hI]
  have h1 := le_abs_self (x ⬝ᵥ (Δ *ᵥ x))
  have h2 := mul_lt_mul_of_pos_right hc hxx
  linarith

/-- The same with the fill-reducing permutation used by sparse Cholesky: the factor is computed for
`Q (B - c I) Qᵀ`, i.e. for `B.submatrix σ σ`. -/
theorem posDef_of_shifted_factor_perm {n : ℕ} (B L Δ : Matrix (Fin n) (Fin n) ℝ)
    (σ : Equiv.Perm (Fin n)) {c γ R Cc : ℝ}
    (hfac : L * Lᵀ = B.submatrix σ σ - c • (1 : Matrix (Fin n) (Fin n) ℝ) + Δ)
    (hΔ : ∀ i j, |Δ i j| ≤ γ * ∑ k, |L i k| * |L j k|) (hγ : 0 ≤ γ)
    (hrow : ∀ i, ∑ k, |L i k| ≤ R) (hcol : ∀ k, ∑ i, |L i k| ≤ Cc) (hc : γ * R * Cc < c)
    (x : Fin n → ℝ) (hx : x ≠ 0) : 0 < x ⬝ᵥ (B *ᵥ x) := by
  have hy : x ∘ σ ≠ 0 := by
    intro h; apply hx; funext i
    simpa using congrFun h (σ.symm i)
  have h := posDef_of_shifted_factor (B.submatrix σ σ) L Δ hfac hΔ hγ hrow hcol hc (x ∘ σ) hy
  have e1 : (x ∘ σ) ∘ σ.symm = x := by funext i; simp
  rw [submatrix_mulVec_equiv, e1] at h
  have e2 : (x ∘ σ) ⬝ᵥ ((B *ᵥ x) ∘ σ) = x ⬝ᵥ (B *ᵥ x) := by
    simp only [dotProduct, Function.comp]
    exact Equiv.sum_comp σ (fun i => x i * (B *ᵥ x) i)
  rwa [e2] at h

end Cholesky

section Arithmetic

/-- Rung 1: Riemann–Hurwitz for `C → C/A₇` with signature `(2,4,7)` gives genus 136. -/
theorem genus_C : (2 * 136 - 2 : ℚ) = 2520 * (-2 + (1 - 1 / 2) + (1 - 1 / 4) + (1 - 1 / 7)) := by
  norm_num

/-- Rung 1: Riemann–Hurwitz for `C → D = C/⟨τ⟩` with 18 fixed points gives `g(D) = 64`. -/
theorem genus_D : (2 * 136 - 2 : ℤ) = 2 * (2 * 64 - 2) + 18 := by norm_num

/-- NOTES 7.2: with `Area(C) = 540π`, Hersch's inequality `λ₁ Area ≤ 8π d` and the certified
`λ₁(C) ≥ 0.33335` force every map `C → ℙ¹` to have degree at least 23. -/
theorem gon_C_ge_23 {lam : ℝ} (hlam : (0.33335 : ℝ) ≤ lam) {d : ℕ}
    (hd : lam * (540 * Real.pi) ≤ 8 * Real.pi * d) : 23 ≤ d := by
  have hpi := Real.pi_pos
  have h1 : (0.33335 : ℝ) * 540 ≤ 8 * d := by nlinarith
  have h2 : (22 : ℝ) < d := by nlinarith
  exact_mod_cast (show (22 : ℕ) < d by exact_mod_cast h2)

/-- NOTES 7.2: a degree-`d` map `D → ℙ¹` gives a degree-`2d` map `C → ℙ¹`, so `d ≥ 12`. -/
theorem gon_D_ge_12 {lam : ℝ} (hlam : (0.33335 : ℝ) ≤ lam) {d : ℕ}
    (hd : lam * (540 * Real.pi) ≤ 8 * Real.pi * (2 * d)) : 12 ≤ d := by
  have hpi := Real.pi_pos
  have h1 : (0.33335 : ℝ) * 540 ≤ 16 * d := by nlinarith
  have h2 : (11 : ℝ) < d := by nlinarith
  exact_mod_cast (show (11 : ℕ) < d by exact_mod_cast h2)

/-- NOTES 7.8: with the sharper certified value `λ₁(C) ≥ 0.34089` (mesh `n = 96` for `Q₁`),
every map `C → ℙ¹` has degree at least 24. -/
theorem gon_C_ge_24 {lam : ℝ} (hlam : (0.34089 : ℝ) ≤ lam) {d : ℕ}
    (hd : lam * (540 * Real.pi) ≤ 8 * Real.pi * d) : 24 ≤ d := by
  have hpi := Real.pi_pos
  have h1 : (0.34089 : ℝ) * 540 ≤ 8 * d := by nlinarith
  have h2 : (23 : ℝ) < d := by nlinarith
  exact_mod_cast (show (23 : ℕ) < d by exact_mod_cast h2)

end Arithmetic

end Hilbert13.Spectral
