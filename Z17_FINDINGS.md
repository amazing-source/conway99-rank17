# ℤ¹⁷ findings (structures, night 27→28/09)

Labels: PROVEN / DERIVED / PATTERN (sampled). Nothing here counts as an exclusion until PASS 2 finishes.

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
4. **Contrast (C2's data):**
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

## Scale

- Pair 0 has 13 top-level branches and pair 4 has 23. After 1–3 minutes each is still in its first branch, with 10–15k distinct skeletons.
- Estimate: 10⁵–10⁶ distinct skeletons per pair.
- PASS 2 at about 3–15 distinct skeletons per second per thread means hours to days per pair. More threads, or a lemma built from the AC cores, would change that.
