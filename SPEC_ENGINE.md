# Vector-level engine: spec (structures, 28/09 00:10). Generic in V; ℤ¹⁷ reference implementation in z17/stageB.py

**Input.** V (Gram matrix), y1, y2, K = V ∩ y^⊥ (basis and Gram), and a skeleton: classes g_2..g_7 ∈ K/2K, with class(1,d) = g_d
and class(c,d) = g_c + g_d.

## 0. Parity functional (depends only on the skeleton)

For x ∈ V, put p(x) := (0, p_2, …, p_7) ∈ F₂⁷ with p_e = (x, w_{1e}) mod 2, where w_{1e} is any vector of class g_e.

Canonical form mod 𝟏: if bit 1 is set, complement all 7 bits. So ε_1 ≡ (0,1,1,1,1,1,1) and ε_c for c ≥ 2 is itself.

## 1. Cell domains

For each cell c < d, the domain is the set of ordered signed pairs (w⁺, w⁻) of norm-10 vectors in class(c,d) with:
- w⁺ ≠ ±w⁻ and (w⁺, w⁻) ∈ {0, ±4};
- **M1 (orientation):** p((w⁺ + w⁻)/2) ≡ ε_c (mod 𝟏).
  - Then b = (w⁺ − w⁻)/2 automatically satisfies p(b) = ε_d.
  - Both orders (w⁺, w⁻) and (w⁺, −w⁻) must be tried; this condition selects the right one.
- **L1 (refined):** each of w⁺ and w⁻ must avoid a norm-1 vector in U. That means: no v ∈ V with |w + 2v|² = 2 and
  p(v) ∈ {0, ε_c, ε_d, ε_c + ε_d} (mod 𝟏).
  - Cheap sufficient part: w ≡ r mod 4V for a root r (v ∈ 2V, p = 0).
- **Sign normalisation** (the h_c flips act by a_cd → −a_cd for d > c and b_dc → −b_dc for d < c; a single representative
  survives per orbit): a_12, b_12 and b_1d for d ≥ 3 have a positive first non-zero coordinate (in any fixed basis).
- An empty domain means the skeleton is dead (reason "empty cell").

## 2. Pairwise rule (M2 and the Gram condition)

For cells X ≠ Y and s, t ∈ {+,−}, require

  ((w_s^X, w_t^Y) + r_s^X · r_t^Y) / 4 ∈ {0, ±1}

with r⁺ = e_c + e_d and r⁻ = e_c − e_d (cell c < d).

## 3. Search

- **Order:** the cells of vertex 1 first, (1,2), …, (1,7); then (2,3), …, (2,7); and so on (lexicographic). Star c closes at the
  last cell containing c.
- **Forward checking:** after each assignment, filter every unassigned domain by rule 2. An empty domain kills the node.
- **Star closure:** Σ_{d>c} (w⁺_cd + w⁻_cd) + Σ_{d<c} (w⁺_dc − w⁻_dc) = 0.
- **Frame:** the residual R = 28 P_W − Σ (w⁺w⁺ᵀ + w⁻w⁻ᵀ) must stay PSD (check the minimum eigenvalue, tolerance 1e-7).
  - In K-coordinates, use 28 G_K⁻¹ − Σ c cᵀ.
  - At the leaf, R = 0 exactly.
- **Verdicts:**
  - EXCLUDED, with reason empty-cell, root (1 node) or deeper(n nodes);
  - FOUND, with the solution;
  - LIMIT, when the node cap is hit.
  - FOUND and LIMIT stop everything (policy).

## 4. Certificates for pass 2

- **Empty cell:** the cell and class, plus a recount of its norm-10 vectors. Recheckable independently.
- **Root / deeper:** the search is exhaustive. The record is (skeleton, node count, cell order); re-running must reproduce it.
- **Exact arithmetic:** all tests are integer except the PSD check. That check only prunes, so a floating-point error there
  could lose a solution but never create a false one. Pass 2 re-runs with an exact LDLᵀ over ℚ at pruned nodes.

## 5. ℤ¹⁷ specifics (for reference)

- Classes are odd supports.
- Weight-10 classes: ±1 on S.
- Weight-6 classes: ±1 on S plus one ±2 off S.
- Weight-2 classes: ±1 on S plus two ±2 off S. The shape (3,1) is removed by L1.
- **Pair 0 (PROVEN, from L3):** a weight-10 class disjoint from the 3 free coordinates has (w, w′) ≡ 2 mod 4 for all w, w′, so it
  is never usable.
