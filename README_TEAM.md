# conway99-rank17 — rank-17 switch for the involution case (team copy)

- `RANK17_HANDOFF.md` and the helper files come from the user's GPT Pro handoff of 27/09/2026 (unchanged, SHA256SUMS.txt).
- The core was re-derived independently by structures (laptop) and C2 (lab PC):
  - U ∩ H = 2ℤ⁷ + ℤ𝟏;
  - U′ = U_H + ⊕ℤ(h_c/2) = ℤ⁷ ⊥ V, with V unimodular of rank 17, in 9 types;
  - w_u ∈ K = V ∩ y^⊥, det 49;
  - the 7 star equations;
  - Σ ww^T = 28 I₁₅.

## Scope (user, 27/09 night)

- Finish the 9 lattices V: 1272 y-pair cases in total.
  - 8 in ℤ¹⁷.
  - 1264 in the other 8 (C1, conway99-c1 notes/146).
- It is a necessary-condition search: every case must be killed.
- Policy (C2 relaying the user): PASS 1 samples fast; PASS 2 re-runs exhaustively with certificates.
- Verdict labels: PROVEN / COMPUTED / PATTERN. They are never mixed, and a survivor stops everything.

## Lanes

- **structures:** V = ℤ¹⁷, all 8 pairs. Engine in `z17/`. Spec in `SPEC_ENGINE.md`.
- **C2:** an independent ℤ¹⁷ implementation, plus E₈⊕ℤ⁹ (53) and D₁₂⁺⊕ℤ⁵ (158) with their own general-V engine.
- **C1:** the other 6 V's. Skeleton test, then stage B with `z17/sbdfs_g.exe` (general V, Gram input).

## Lemmas in use (all PROVEN unless marked)

- **L1:** min U ≥ 2 forces min(w_u + 4V) ≥ 10.
  - In ℤ¹⁷ there is no (3,1)-shaped w.
  - Refined form: there is no v with |w_u + 2v|² = 2 and p(v) ∈ ⟨ε_c, ε_d, 𝟏⟩.
- **L2 (g-form):** p(a_cd) = ε_c and p(b_cd) = ε_d. This is **M1**.
- **L3 (ℤ¹⁷, pair 0):** every a_cd and b_cd has exactly one ±1 on the 3 free coordinates.
  - Consequence: a weight-10 class that misses them is unusable.
- **L4 (ℤ¹⁷):** U = ½(Q + 4ℤ²⁴), where Q is a self-dual ℤ₄ code of type 4⁶2¹².
- **M2:** (w_X, w_Y) ≡ −r_X·r_Y mod 4, i.e. self-orthogonality of the ℤ₄ lift.
- **PE:** 7|x|² − Y(x) ≥ k(7−k)/2 for all x ∈ V, applied to every κ ∈ K* with |κ|² < 6/7.

## ℤ¹⁷ engine (`z17/`)

- `stageA.py` / `run_pair.py`:
  - π-row DFS with partition refinement on the y-coordinate classes;
  - S₇ breaking: N_12 minimal, then N_13 ≤ … ≤ N_17;
  - usability prune;
  - all PE items pruned during the DFS;
  - `--dedupe` gives the full S₇ canonical form.
- `stageB.py`: cell domains.
  - M1 is keyed on (w·π_e) mod 4.
  - Refined L1 is vectorised.
  - Sign normalisation for the 2⁷ h-flips.
- `sbdfs.c`, `sbdfs3.c` (bitsets + 2-open star look-ahead), `sbdfs_g.c` (general V): exact DFS.
  - All reproduce the Python node counts exactly.
  - `sbdfs_g` was also tested in a random basis.

## Status (28/09 ~04:40)

### ℤ¹⁷

- **Pair 7: COMPUTED.** Exhaustive stage A has 0 skeletons.
- **Pairs 0 and 6:** exhaustive PASS 2 running on `z17pipe17`, 0 survivors.
  - Pipeline: stage A + incremental M1 + S₇ dedupe, then exact Z4, then the Z4 code-consistency filter, then plain DFS, then AC + the per-code AC split + dynamic-ordering bitset DFS.
- **Queue:** pairs 1 and 2 after pair 0; pairs 3 and 5 after pair 6; pair 4 after pair 7's stats.
- **Independent check:** C1 is running all 8 pairs through its engine (`c1_recherche`: `rank17_z4.c`).

### Other V's

- **C1's engine** (PASS 2, exhaustive):
  - canonical skeleton enumeration under W(R_K), with the pointwise stabiliser chain;
  - the exact Z4 test on partial and full skeletons.
  - A11E6⁺ #31 (10¹⁰ raw skeletons) is decided in 4 s. E8²⊕ℤ is fully excluded.
  - All 1,053 pairs of C1's 6 V's are expected in about 2 h.
- **C2 (E₈⊕ℤ⁹, D₁₂⁺⊕ℤ⁵):** PATTERN-dead so far, using LP + AC + MITM.
  - D₁₂⁺ is AC-resistant: there the frame rows are needed.
  - The lab PC has no compiler.

### Key lemma material

See `Z17_FINDINGS.md`: the lift-layer kill is purely ℤ₄, and the exact Z4 test uses 21 bits.
