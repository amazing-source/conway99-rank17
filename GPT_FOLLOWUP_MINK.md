Follow-up to your graph-star spectral bound (same notation).

RESULT OF YOUR THEOREM
- min(K*) ≥ 5/7 was read and checked. Two independent exact scans agree: 492 of the 1,272 y-pairs have a dual vector with
  7κ² ≤ 4 and are excluded. Per type: ℤ¹⁷ 5, E8² 7, D16⁺ 10, E8 32, D12⁺ 78, E7²⁺ 82, D8²⁺ 107, A15⁺ 80, A11E6⁺ 91.
- After this and the earlier computational exclusions, 754 y-pairs are open.
- Among the 780 pairs not killed by the theorem, min 7κ² is:
  - 5 for 438 pairs;
  - 6 for 239 pairs;
  - > 6 for 103 pairs.
- Your §7.2 rule (norm 5/7 ⇒ cut size exactly 2) is now a prefix filter. It cuts the two big ℤ¹⁷ trees only 2–5×.

THE HIGH-VALUE QUESTION
Push the bound further.
- min(K*) ≥ 6/7 would kill the 438 pairs with norm-5/7 vectors outright.
- min(K*) ≥ 1 would also kill the 239 with norm 6/7, leaving about 103 pairs.

1. n = 7κ² = 5, k = 2 (the only case left by §7.2). The vector a = Wκ is ±1 on the 20 labels of the 10 cut cells (c,x) and
   (c',x), x ∉ {c,c'}, and 0 elsewhere, including cell (c,c'). It satisfies Na = 4a and aᵀa = 20.
   - Each cut u needs (Na)_u = 4a_u.
   - Within its own star it gets at most a_{Q(u)} − a_{P(u)} (|·| ≤ 2).
   - From the other cut star it gets at most 4 entries, by the outside-row bound.
   - Every uncut u needs (Na)_u = 0.
   Can the P/Q matching structure of the three stars at c, c' and x, plus the outside-row bounds, λ = 1 and μ = 2, rule this out?
   For example: count the ±4 requirement over the 20 cut labels against the matchings, or use the star-x Gram bound on the
   restriction of a to star(x) (entries ±1 on the 4 labels of cells (c,x) and (c',x)).
2. n = 6 (κ² = 6/7, aᵀa = 24): k ∈ {0,1,2,3}. Can each k be ruled out as in §4/§7.2?
   - k = 3 forces a = ±1 on exactly the 24 cut labels.
   - k = 0 gives b = a/2 with b² = 6 and Nb = 4b.
   - k = 1 and k = 2 have small even remainders.
3. If a complete proof is out of reach, give the strongest PROVEN restriction on the cut size k (or on the pattern of a) for n = 5
   and n = 6 that we can test on partial skeletons.
Label every claim PROVEN / DERIVED / CONJECTURE.
