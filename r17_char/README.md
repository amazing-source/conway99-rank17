# Rank-17 characteristic test: 228 of the 751 open (V, y-pair) classes closed (C2, 29/09)

**Source list:** C1 `data/rank17/open_after_sol.json` + `ypairs_*.json` at conway99-c1 888e731 (754 open; 3 are Z17).

## Theorem (Astra 29/09, `ASTRA_PROOF_29_09.md`; re-proved by C2)
In the switched frame, the 42 residual vectors a_cd, b_cd lie in K = V ∩ ⟨y1,y2⟩⊥, all have ODD norm (3, 5, 7), and
Σ_x x xᵀ = 14 P_K (handoff §4). So Σ_x x xᵀ + 2y1y1ᵀ + 2y2y2ᵀ = 14 I on V. For any characteristic χ of V, every (χ,x) is odd:

    (χ,y1)² + (χ,y2)² ≤ 7|χ|² − 21.

(Height 7 against a norm-9 χ kills y1 whatever the partner: 49 + 1 > 42.)

## Completeness of the test (C2)
The bound says min Char(K*) ≥ 3. A violating class lifts, after subtracting even multiples of y1, y2 (χ − 2y stays characteristic),
to χ with odd heights |h| ≤ 7, hence |χ|² < 3 + 98/7 = 17. Characteristic norms are ≡ 1 mod 8, so **χ of norm 1 or 9 decide the test
exactly**.

## Result (`r17_char_scan.py`, exact integer arithmetic; `rank17_char_kills_2909.json` = one witness per class)
| V | open | closed | partner-independent |
|---|---|---|---|
| E8 | 21 | 5 | 3 |
| D12p | 79 | 20 | 12 |
| E72p | 108 | 30 | 16 |
| A15p | 122 | 34 | 17 |
| D82p | 208 | 65 | 30 |
| A11E6p | 213 | 74 | 36 |
| **total** | **751** | **228** | 114 |

Z17 (3 classes) is untouched: Z17 has no characteristic vector of norm ≤ 9. Astra's 9 named cases (E8:5,9,11,14,19; A11E6p:8,9,10,102)
are all included; A11E6p:100 survives (as Astra reports).

Each witness: V, id, chi_coef (integer coefficients in the V basis = int_basis(core) ⊕ Z^m, as in rank17_cores.py + ypairs coordinates),
chi_norm, heights, bound. The script re-checks: characteristic congruence, norms, Gram 7I of the pair, integrality of y coefficients.
It needs `frame_dfs.py` (C2 lab: `short`, `solve_mod2`, `one_per_pair`) and conway99-c1 `code/rank17_cores.py`.
Conditional on the established switch (RANK17_HANDOFF §3–4); scope = these rank-17 classes, not new 156-envelope exclusions.
