# V = ℤ¹⁷ subcase: search design (structures, 27/09 night)

Notation as in RANK17_HANDOFF.md.
- K = y^⊥ ∩ ℤ¹⁷, of rank 15.
- For each cell {c,d} (c < d) there are two vectors w⁺_cd and w⁻_cd in K of norm 10: Same and Opp, oriented as in §4.

## Necessary conditions used

1. **Glue parity** (handoff §5–6). There are π_1,…,π_7 ∈ F₂¹⁷ such that w_cd ≡ π_c + π_d (mod 2), for both w⁺ and w⁻.
   - The π are defined up to adding one common vector (complementing a column). Normalise π_1 = 0.
   - (a) **Norm 10 parity.** A norm-10 integer vector has an odd support S with |S| ≡ 2 (mod 4). Hence wt(π_c + π_d) ∈ {2,6,10} for
     all c ≠ d.
     - With π_1 = 0: wt(π_c) ∈ {2,6,10} and |π_c ∩ π_d| is odd.
     - So the Gram matrix of π_2..π_7 mod 2 is I + J (6×6, invertible): a plus-type E₇ embedding.
   - (b) **Units (min U ≥ 2).** Column i of (π_1..π_7) is non-constant. With k_i = min(#1, #0) we need Y(e_i) := y₁[i]² + y₂[i]² ≤
     4, 2, 1 for k_i = 1, 2, 3 (parity-energy inequality). With π_1 = 0, the supports of π_2..π_7 cover all 17 coordinates.
   - (c) **⊥ y parity.** w ⊥ y_j forces |S ∩ O(y_j)| even, where O(y_j) is the set of odd coordinates of y_j. So every π_c ⊥ 𝟏_{O(y_j)}
     (mod 2), for j = 1, 2.
2. **Coordinate budgets** (diagonal of Σ ww^T = 28 P_K): Σ_w w_i² = 28 − 4Y(e_i).
   - The w with odd i-th entry number 2k_i(7−k_i) ∈ {12, 20, 24}.
   - Off-diagonal: Σ_w w_i w_j = −4(y₁[i]y₁[j] + y₂[i]y₂[j]).
3. **Gram.** (w_u, w_v) = 4s − r_u·r_v with s ∈ {0, ±1} for u ≠ v.
   - Disjoint cells: {0, ±4}.
   - The same cell: (w⁺, w⁻) = 4t, t ∈ {0, ±1}.
   - Cells sharing a block: values ≡ −r_u·r_v (mod 4), in [−5, 5].
4. **Stars.** For each c: Σ_{d>c}(w⁺_cd + w⁻_cd) + Σ_{d<c}(w⁺_dc − w⁻_dc) = 0.
5. **Tight frame.** Σ w wᵀ = 28 (I − y₁y₁ᵀ/7 − y₂y₂ᵀ/7).

## Stages

- **Stage 1: π-enumeration.** Enumerate the 6-word systems (π_2..π_7) satisfying 1(a)(b)(c) and pairwise distance ≤ 10, up to:
  - the coordinate symmetries of the y-pair (signed permutations fixing {±y₁, ±y₂});
  - S₇ on blocks.

  Search: DFS over words of weight 2, 6 or 10, with pairwise intersection parity and distance pruning, and canonical ordering.
- **Stage 2: per π.** Build the candidate lists per cell: norm-10 vectors of K with odd support S_cd (at most a few hundred). Then run
  an LP relaxation of 2–5 (variables per cell and candidate, orientation as Same/Opp) and branch-and-bound with exact certificates,
  as for the node LP.
- **Output:** for each of the 8 y-pairs, either "no π" or "all π excluded (certificates)". Otherwise, surviving structures for inspection.
