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

- Partial Z4 as a pruning test (28/09): at depths 2–5 it killed about 1% of nodes on a 16-walk sample, so it is useless.
  - WITHDRAWN (28/09 ~09:10): the claim that "without L1 even the full test kills ~1%, so L1 does the work". It rested on those
    16 walks and on a test that is not equivalent to z4_ok.

## Correction (28/09 ~09:10): the Z4 kill rate on ℤ¹⁷

- The "97% of distinct skeletons killed by Z4" (pair 0) came from the first ~1% of the tree in DFS order, so it is biased.
- Unbiased Knuth estimates:
  - pair 0: 77% of leaves pass Z4 (2048 walks, se ~12%);
  - pairs 1 and 2: 46–87% pass, depending on the sample.
- So on ℤ¹⁷ most leaves need the integer stage.

## Sol's short-dual-class lemma (DC, PROVEN; proof read by structures) on ℤ¹⁷

- Implemented as a stage-A prefix filter (`dc` flag, z17pipe22+). The forbidden classes are {[7κ] : 7|κ|² ∈ {2,4,6,8,10,14}}:
  41–84 per pair, i.e. all weight-2 classes plus O1\{i} and O2\{i}.
- The Knuth leaf estimates are unchanged within ±5% for every pair: on ℤ¹⁷ DC cuts about 1×.
  (C1's types: about 9× on A11E6⁺ #0, 2× on E7²⁺ #20, ~100× on D8²⁺ #26.)
- Class-energy extension (`ce N` flag, leaf level). For every κ with [7κ] in the span:
  - the congruence (κ,w_u) ≡ −t·r_u (mod 4);
  - the budget Σ s_u² = (n + 3|t|²)/4, using per-cell minima over the real domains, with λ = 0.
  It killed 0 of the sampled Z4 survivors. It needs Sol's λ (the star sums), and it only bites at slack ≤ 1, which is rare on ℤ¹⁷.
- Remark: PE is exactly Sol's inequality (IB) applied to x ∈ U_H with λ = 1_P, and DC is IB applied to the half-lifts.
  Together they exhaust the single-vector IB arguments on U.

## Scale (Knuth estimates, 28/09 ~06:40, PATTERN)

- Mode `knuth NW seed` of z17pipe (z17pipe20.exe): random root-to-leaf walks of stage A, one child drawn uniformly per level.
  The product of the branching factors is an unbiased estimate of the leaf count.
- 8192 walks per pair, a few ms each. Counts are stage-A leaves including S₇ duplicates (observed dup rate 80–99%).

| pair | leaves (se) | per-depth nodes, depths 1–5 |
|---|---|---|
| 0 | 1.7e8 (9e6) | 13, 164, 5.1e3, 2.8e5, 1.5e7 |
| 1 | 1.6e10 (1.5e9) | 46, 2.4e3, 2.1e5, 2.6e7, 1.3e9 |
| 2 | 7.2e9 (1.3e9) | 63, 4.2e3, 4.5e5, 6.0e7, 1.4e9 |
| 3 | 2.3e6 (8e5) | 29, 763, 1.7e4, 5.2e5, 2.2e6 |
| 4 | 5.4e7 (9e6) | 23, 554, 1.9e4, 9.8e5, 1.8e7 |
| 5 | 3.6e6 (1.1e6) | 42, 1.6e3, 2.8e4, 6.1e5, 6.7e6 |
| 6 | 4.4e5 (1e5) | 29, 782, 1.9e4, 5.8e5, 2.4e6 |
| 7 | 0 | 32, 723, 737, 2.8e3, 1.2e4 (matches COMPUTED) |

- Throughput of z17pipe19 in real runs: pair 0 ~250 leaves/s, pair 6 ~50 leaves/s (99% duplicates, canon-bound).
  - Of the distinct skeletons in pair 0, 97% die in the exact Z4 test.
  - The integer solver spends ~0.6 s on each Z4 survivor; pairs 4–6 average ~1 s per skeleton.
- So on one thread:
  - pairs 6 and 3: hours;
  - pair 5: about a day;
  - pairs 4 and 0: days to a week;
  - pairs 1 and 2: years.
- Better dedupe and a φ-tree solver are worth ~2–10×. That is not enough for pairs 1 and 2: they need a lemma that cuts the skeleton tree.
- The old estimate "10⁵–10⁶ distinct skeletons per pair" was wrong for pairs 0, 1, 2 and 4.
- Knuth walks with full leaf processing (`knuth 512 seed solve`, z17pipe21, no dedupe; se ~30–70%):

  | pair | leaves | leaves passing Z4 (share) | processing s (before dedupe) |
  |---|---|---|---|
  | 1 | 2.7e10 | 1.2e10 (46%) | 1.3e10 |
  | 2 | 4.7e9 | 2.3e9 (48%) | 1.2e9 |
  | 3 | 4.2e5 | 5.1e4 (12%) | 9.3e4 |
  | 5 | 7.6e5 | 7.4e5 (97%) | 7.2e5 |

  - The Z4 kill rate depends strongly on the pair: 97% for pair 0, but only about half for pairs 1 and 2.
  - No sampled leaf was FOUND or LIMIT.

## Sol's graph-origin star rules (28/09 ~09:40–10:20)

- Rules (PROVEN by Sol from λ = 1, μ = 2; signs checked by C1):
  - (R1) in a star's local orientation the u-Gram is Q − P, two perfect matchings: each u has at most one +1 and one −1 there,
    in equal numbers.
  - (R2) an outside row has at most 2 of each sign, again equal.
  - (R3) every u has exactly 10 neighbours.
- In the stage-B solver (`sbdfs_gs.exe ... star`): the pair-3 LIMIT skeleton (59 M nodes in z17pipe) is EXCLUDED in 61,505 nodes
  (12 s). Every branch dies within the first four anchor cells.
- As a stage-A prefix test on the anchor star (z17pipe `astar k`), with the star-equation look-ahead:
  - Look-ahead rule (PROVEN): with m anchors missing, (T,v) ∈ {−m−2, −m, −m+2} for every assigned label v, with the matching
    count conditions.
  - Pair 0: no Knuth walk reaches a leaf. Every complete anchor star dies.
  - But the depth-5 node count only halves (8e6 vs 1.4e7 over 16 walks; the smaller 8-walk figures were noise), and each test
    costs ~0.14 s.
  - So it is a strong LEAF filter but not a tree cut. Dropped.
- V-free model `n5_signed.py` (CP-SAT):
  - contents: S1 matchings, S2 outside rows, S3 degree 10; optionally the n = 5 cut vector with Na = 4a (Sol's boxed lemma) and
    N² = N + 12I − G_r;
  - without N², both variants are UNKNOWN at 3–5 min, and C1 points out that this relaxation is too weak to be informative;
  - a 30-min, 6-worker run of n = 5 + N² is in progress.
