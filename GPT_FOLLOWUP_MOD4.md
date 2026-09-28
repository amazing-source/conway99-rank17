CAVEAT (added 28/09 ~05:10): the "purely mod 4" pattern is NOT universal. Some skeletons (ℤ¹⁷ pair 6, and C1's A11E6⁺ / D8²⁺ / E7²⁺)
pass the exact mod-4 test and die only on the integer rules |(u,v)| ≤ 1, the stars and the frame. So ask for a mod-4 lemma covering
MOST skeletons, or for a lemma that combines mod 4 with the integer bounds.

FOLLOW-UP (same notation as the previous prompt): the lift-layer kill is purely mod 4.

NEW DATA (V = ℤ¹⁷, PATTERN)
- Keep each cell's domain exactly as before: ordered pairs (w⁺,w⁻) of norm-10 vectors in K with odd support S_cd,
  (w⁺,w⁻) ∈ {0,±4}, M1: p((w⁺+w⁻)/2) = ε_c, and L1.
- Between different cells, use ONLY the mod-4 rule (w_X, w_Y) + r_X·r_Y ≡ 0 (mod 4), i.e. ℤ₄ self-orthogonality of the codewords r_u + w_u.
- Then plain arc consistency (no search, no stars, no frame) empties some cell domain for exactly the same skeletons as the full
  integer rule |(u,v)| ≤ 1. Tested on 12 random skeletons of pairs 0 and 4: 11 wiped out under both rules, 1 survived both.
- Minimal wipeout cores have 3–6 cells. The smallest is three pairwise disjoint cells {1,2}, {3,4}, {5,6} with class weights 2, 10, 6.
  - Every choice (w⁺,w⁻) for the weight-2 cell has a mod-4-compatible partner in the weight-10 cell or in the weight-6 cell, never both.
  - The mechanism: twos-parity. When S_X ∩ S_Y = ∅, (w_X,w_Y) ≡ 2(|T_X ∩ S_Y| + |S_X ∩ T_Y|) mod 4, where T = the positions of the ±2 entries.
    This combines with the ⊥y balance and with L3.

QUESTION
In ℤ₄ language: residue D = {(ξ, ρ(ξ)) : ξ ∈ E} of dimension 6, doubly even, length 24.
- Can a self-orthogonal ℤ₄ code with this residue contain, in each of the 21 residue classes d_cd, two codewords
  q± = (e_c ± e_d | w±) whose f-parts have Euclidean weight exactly 10, in the shapes (1¹⁰), (2,1⁶), (2²,1²)?
- Two further requirements on those codewords:
  - (w⁺,w⁻) ∈ {0,±4} and the orientation M1 (the latter is automatic inside one code);
  - the integer conditions w ⊥ y1, y2.
Prove impossibility (for ℤ¹⁷ and the 8 y-pairs, or a sub-case), or describe the structure of a solution.
A natural route might be the quadratic refinement φ with λ(d) = d̃ + 2φ(d), φ(d+d′) = φ(d) + φ(d′) + d∗d′ mod D^⊥,
together with the Euclidean-weight-10 constraint on the f-parts.
Label everything PROVEN / DERIVED / CONJECTURE.
