# ℤ¹⁷ findings (structures, night 27→28/09)

Labels: PROVEN / DERIVED / PATTERN (sampled). Nothing here counts as an exclusion until PASS 2 finishes, except:
- **Pair 7: COMPUTED, excluded at the skeleton level.** Its exhaustive stage A (PE + usability + incremental M1 + coverage) has 0 leaves (36 s).
  C1's independent engine confirms it (28/09 ~05:30): no skeleton survives.
- Pairs 0–6 are MOD4_FEASIBLE in C1's engine, so the integer search is needed. That is z17pipe19, PASS 2, running.

## Kill mechanisms (PATTERN, pairs 0, 1, 2, 4, 5)

1. **M1 (orientation).** Cell domains need p((w⁺+w⁻)/2) = ε_c. This kills about 30% of usable skeletons outright: some cell ends up empty.
2. **Lift layer alone kills.** Turn off the stars and the frame/PSD, and keep only the pairwise rules (G)/M2 plus the M1/L1 domains.
   - The relaxed exact search still kills every resolved skeleton (32/32 random ones and the 4 heavy leads).
   - Two cases are open at a 500k-node cap.
   - Nothing was FOUND.
3. **Arc consistency on the pairwise rules**, with no search, wipes out a domain in 21 of 23 random skeletons.
   - Only 2/23 have a dead cell pair, i.e. an all-zero compatibility matrix.
   - Minimal wipeout cores have 3–6 cells.
   - The smallest core found is 3 pairwise disjoint cells (1,2), (3,4), (5,6) with class weights (2, 10, 6). Every X entry has a Y-partner or a Z-partner, never both. The reason is twos-parity mod 4, together with the ⊥y balance and L3.
4. **Purely ℤ₄ (PATTERN, 12 skeletons, pairs 0 and 4).** Replace the integer cross-cell rule (w,w′)+r·r ∈ {0,±4} by the mod-4 rule alone,
   (w,w′)+r·r ≡ 0 mod 4, and keep the cell domains as they are (the integer pair rule, M1, L1, the shapes).
   - AC wipes out exactly the same skeletons: 11/11 wipeouts agree, and the 1 surviving skeleton survives both versions.
   - So the |(u,v)| ≤ 1 bound plays no role in the kill. The obstruction is ℤ₄ self-orthogonality between the shaped lifts of the 21 residue classes.
   - This is the candidate for a uniform lemma: a statement about self-dual ℤ₄ codes of type 4⁶2¹² whose residue is the graph of ρ.
   - **CORRECTION (28/09 ~05:10): not universal.** The ℤ¹⁷ pair-6 lead 040c0408… is mod-4 feasible (and LP-feasible). It dies only in
     the integer search.
     C1 finds hundreds of thousands of mod-4 survivors on A11E6⁺, D8²⁺ and E7²⁺. There the kill is integer-level (star 1).
     So the mod-4 layer kills most skeletons (99% of pair 0), but NOT all. A uniform proof needs the integer layer too.
5. **Contrast (C2's data):**
   - E₈⊕ℤ⁹: AC kills 15/15, like ℤ¹⁷.
   - D₁₂⁺⊕ℤ⁵: AC kills 0/12. The per-skeleton LP kills are driven by the frame rows ("one pair per cell + frame" is infeasible in 21/23), not by the stars.

## Engine speed-ups (all verdict-neutral, checked)

- Incremental M1 in stage A (55% of leaves pruned) and S₇ canonical dedupe (79–87% of leaves are duplicates).
- Star look-ahead with 1, 2 and 3 open cells. The 3-open version uses pair sum-sets of the last three cells.
- Bitset tier:
  - lazy compatibility rows;
  - full arc consistency first;
  - then the search, ordered dynamically by smallest domain first. This gives 5–500× fewer nodes on heavy cases.
- Keep the first tier in fixed order. Dynamic ordering there is 10× slower, because it breaks the first-cell kill.

## Exact mod-4 test (Z4) — the current main filter (DERIVED, checked)

- If the frame codewords q_u = (r_u | w_u) are pairwise orthogonal mod 4, they lie in one self-dual ℤ₄ code Q with residue D.
- Such codes correspond to φ with
  - φ(d+d′) = φ(d) + φ(d′) + d∗d′, and
  - φ(d)(d′) + φ(d′)(d) = |d∩d′|/2.

  On the basis b_c = d_{1c}, this is a 6×6 matrix M with M + Mᵀ fixed: 21 bits.
- Each cell gives a signature set Σ_cd ⊂ F₂⁶, where σ(q) = ((q − q̃)/2 mod D^⊥) and M1 is σ(q⁺) = σ(q⁻). The skeleton is mod-4
  feasible ⟺ some M satisfies φ(d_cd) ∈ Σ_cd for all 21 cells.
- In C (z17pipe11), this runs before the search. It kills 99.2% of pair-0 and about 92% of pair-4 distinct skeletons, in microseconds.
  It matches z4test.py on 80 leaves.
- Every previously known heavy skeleton (pairs 0 and 4) has 0 feasible M.
- Skeletons that pass Z4 have 350–1,900 feasible M. Each M shrinks the total domain by about 30×.

## Symmetry C: coordinate sign flips (28/09 ~05:00)

- Every coordinate i with y1_i = y2_i = 0 gives a symmetry f_i → −f_i. Pair 0 has 3 such coordinates, pair 4 has 6, pair 6 has 9.
- Breaking it has two parts:
  - (a) the h-flip normalisation looks only at coordinates outside that set, which makes it C-invariant;
  - (b) a lex-leader rule along the search path: the first non-zero entry at each such coordinate (w⁺ before w⁻) must be positive.
- Result:
  - pair-4 heavy cases need up to 63× fewer nodes, with the same verdicts;
  - the pair-6 LIMIT lead (040c0408…, Z4-feasible and LP-feasible, 5·10⁷ nodes) is now EXCLUDED in 72,757 nodes.
- In a general V the same symmetry comes from reflections in unit vectors e ∈ K. They are not in W(R_K).

## Rejected experiments (verdict-neutral, but slower)

- Dynamic ordering in the first (plain) tier: 10× slower.
- Partial Z4 inside stage A: pair 7 goes from 45 s to more than 300 s, because rebuilding partial domains at every deep node is too costly.
- Splitting by feasible M, with a fixed-order search per M: 23× slower than AC plus dynamic ordering on the full domains.

## Scale

- Pair 0 has 13 top-level branches and pair 4 has 23. After 1–3 minutes each is still in its first branch, with 10–15k distinct skeletons.
- Estimate: 10⁵–10⁶ distinct skeletons per pair.
- PASS 2 at about 3–15 distinct skeletons per second per thread means hours to days per pair. More threads, or a lemma built from the AC cores, would change that.
