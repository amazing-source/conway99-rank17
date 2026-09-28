Follow-up to your short-dual-class answer (same notation).

WHAT HAPPENED WITH DC
- Your proof was read and checked step by step: IB, the lifting argument, the table, and the n = 14 case. DC is now a prefix
  filter in both engines.
- Measured effect (random-walk estimates of the skeleton-tree leaves, before/after):
  - V = ℤ¹⁷: about 1×. There are 41–84 forbidden classes per y-pair: all weight-2 classes, plus O1\{i} and O2\{i}. The trees are
    dominated by skeletons with only weight-6/10 cells, and their spans rarely meet a short-dual class.
  - C1's types: about 9× on A11E6⁺ #0, 2× on E7²⁺ #20 and ~100× on D8²⁺ #26. Real, but far from the 10⁴ we need.

CORRECTIONS TO THE DATA WE SENT YOU
- "Mod 4 kills 97% of ℤ¹⁷ pair-0 skeletons" was biased: it came from the first 1% of the tree in DFS order. Unbiased sampling
  gives about 23% for pair 0 and 15–50% for the two hardest pairs. So on ℤ¹⁷ nearly every leaf needs the integer stage.
- "Without the refined L1 it kills about 1%, so min U ≥ 2 does most of the work" is WITHDRAWN (bad measurement).

WHAT WE TRIED ON TOP OF DC
- For every κ ∈ K* with [7κ] in the span: the congruences (κ,w_u) ≡ −t·r_u (mod 4) and the exact budget
  Σ_u s_u² = (n + 3‖t‖²)/4, using per-cell minima over the actual cell domains.
  - With λ = 0 it kills nothing: the star sums Rᵀs = 3t are what make IB bite.
  - It only has force at slack ≤ 1 (n within 4 of your thresholds), which is rare on ℤ¹⁷.
- Remark (please confirm): PE is IB applied to x ∈ U_H with λ = 1_P, and DC is IB applied to the half-lifts P_W x = κ/2.
  So together they exhaust the single-vector IB bounds on U.

QUESTIONS
1. Multi-vector IB. Take several short x_1..x_m ∈ U of PE and DC type, with small budgets. Their pairing vectors
   s^(i) = ((x_i,u))_u ∈ ℤ⁴² have:
   - Gram Σ_u s^(i)_u s^(j)_u = 7(x_i,x_j) − (t_i,t_j) − Σ_k (x_i,y_k)(x_j,y_k);
   - star sums Rᵀs^(i) = 3t_i.
   Does realizing this Gram by integer vectors give obstructions that can be computed from the (partial) skeleton alone? One
   example: two slack-0 structures whose forced supports cannot have the required overlap.
2. ℤ¹⁷ with unit vectors. K = K₀ ⊥ ℤ^z, where z = 3–10 is the number of y-free coordinates. The two hardest y-pairs:
   - (A) y1 = 1⁷0¹⁰ and y2 = (1,−1,0⁵,1⁵,0⁵), with z = 5;
   - (B) y1 = 1⁷0¹⁰ and y2 = (1,1,−1,−1,0³,1³,0⁷), with z = 7.
   Each has ~10¹⁰ leaves, 50–85% pass mod 4, and all die at the integer level. For a y-free coordinate j, the column
   c_j = ((w_u)_j)_u ∈ ℤ⁴² has:
   - entries in {0,±1,±2}, odd exactly on the cut edges of a bipartition of the 7 vertices (sizes p_j ∈ {1,2,3});
   - n₂(j) = (14 − p_j(7−p_j))/2 entries equal to ±2;
   - norm 28 and zero star sums;
   - N·c_j = 4c_j;
   - pairwise orthogonality with every other column.
   Is there a lemma, on these columns or on their interaction with K₀, that kills (A) and (B), or cuts their trees by 10⁴?
3. More generally: anything that cuts the ℤ¹⁷ trees by 10⁴ or more. Label every claim PROVEN / DERIVED / CONJECTURE.
