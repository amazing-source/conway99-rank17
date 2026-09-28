Follow-up to your graph-star spectral bound and your n = 5/6 reduction (same notation). A uniform finish line.

DATA (exact, two independent scans)
- Your min K* ≥ 5/7 kills 492 of the 1,272 y-pairs.
- For the other 780, the exact minimum of n = 7|κ|² over K* \ {0} is: 5 for 438 pairs, 6 for 239, 7 for 20, 8 for 83.
  It is never above 8.
- Hence:

  min K* ≥ 9/7, i.e. no κ ∈ K* with 7κ² ∈ {5,6,7,8}, excludes ALL 1,272 cases: all nine V and the involution.

WHAT IS ALREADY KNOWN
- n ≤ 4: impossible (your theorem).
- n = 5: only k = 2 remains (your boxed statement: a = ±1 on the 20 labels of a 2|5 cut, N_CC a_C = 4a_C,
  N_{C̄C} a_C = 0).
- n = 6: your patterns J3, J2, J1 (1a, 1b) and J0.
  - (1b) is dead by your star bound: its star at the cut vertex has energy 11 + 9 = 20 > ⌊120/7⌋ = 17.
  - J0 (six ±1): after switching by b, each support label has exactly 4 positive neighbours among the other 5 and one
    non-neighbour. So the support is an octahedron K_{2,2,2}.
  - Also |S ∩ star(z)| ∈ {0,2,4}.
  - If it is 4, the induced graph is the 4-cycle made of two antipodal pairs, with ε = σ_z b of opposite signs inside each pair.
- n = 7 (|a|² = 28): the star bound ⌊20n/7⌋ = 20 kills the k = 1 pattern with two ±3. What remains:
  - k = 0: b = a/2 has seven ±1 (max-coordinate argument);
  - k = 1: 12×±1 plus 4×±2 off the cut; or 11×±1 plus one ±3 plus 2×±2; or 12×±1 plus one ±4;
  - k = 2: 20×±1 plus 2×±2; or 19×±1 plus one ±3;
  - k = 3: 24×±1 plus one ±2.
- n = 8 (|a|² = 32, star bound 22): similar short lists (k = 0: eight ±1).
- Every case satisfies N a = 4a, Rᵀa = 0, |a|² = 4n, parity = a cut of size k, and the star bound. N has your star structure:
  - each star is Q − P;
  - outside rows are B − A with |A| = |B| = 2;
  - degree 10;
  - N² = N + 12I − G_r.
  We have a CP-SAT model of exactly this (V-free). It is UNKNOWN after minutes, so we need proofs.

QUESTIONS
1. Prove min K* ≥ 9/7: rule out every (n, k, pattern) above, using the graph structure (λ = 1, μ = 2, P/Q matchings,
   outside-row bounds), not only the linear identities. Your n = 5 k = 2 target is the largest single piece, 438 pairs.
2. If a full proof is out of reach, prove as many rungs as possible, in order n = 5, 6, 7, 8. Each rung kills the pairs whose
   minimum is that n (438 / 239 / 20 / 83).
3. If some rung is genuinely possible, describe the local configuration explicitly (the support of a, the signs, and the
   N-entries around it) so that we can test it against the lattices.
Label every claim PROVEN / DERIVED / CONJECTURE.
