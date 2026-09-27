# conway99-rank17 — rank-17 switch for the involution case (team copy)

- `RANK17_HANDOFF.md` and the helper files come from the user's GPT Pro handoff of 27/09/2026 (unchanged, SHA256SUMS.txt).
- Core re-derived independently by "structures" (laptop) and C2 (lab PC):
  - U ∩ H = 2Z^7 + Z·1;
  - the switch U′ = U_H + ⊕Z(h_c/2) = Z^7 ⊥ V (V unimodular, rank 17, 9 types);
  - w_u ∈ K = V ∩ y^⊥ (det 49; ⟨y1,y2⟩ is primitive because a²+b² is anisotropic mod 7);
  - the star sums as 7 linear equations in K, and Σ ww^T = 28 I_15.
- Lanes (protocol point 3):
  - structures: rank-17 search, starting with V = Z^17 and its 8 y-pair types;
  - C2: an independent second implementation of the Z^17 subcase, plus the mod-2 glue filter;
  - C1: y-pair types for the 8 other V's.
- Status: necessary-condition reduction verified. Nothing excluded yet through this route.
