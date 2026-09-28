# Involution case: the main proof obligation, its dependency graph, and the open branches

*28/09/2026, ~13:00. Written at the user's request: reorganise the work around the main implication chain instead of
collecting side contradictions.*

**DECISION (28/09 ~13:30; C1, with the user's delegation; structures agrees).**
- The finish line is **R6 ∧ R8**. R5, the old target (B), is dropped as a main target; it stays only as a base for J2.
- Lanes:
  - **C1 = R6**: J3, then J2, then J1a;
  - **structures = R8**: k = 0 first, then k = 3, 2, 1.
- Each closed cell gets a note and a one-line [RESULTAT].

**Sources read.** This repo: RANK17_HANDOFF, README_TEAM, Z17_FINDINGS, GPT_FOLLOWUP_MINK/MIN97, out_kstar_all.json and
out_kstar_ladder.json. conway99-c1: notes/139, 146–150 and data/rank17/open_after_sol.json. involution_exploratoire:
ROUTE_C_PROOF.md. involution_lab_c1: audit_refutations.md. c2_experimental_ro: SYNTHESE.md.

**Labels.**

| label | meaning |
|---|---|
| PROVEN | written proof, read by at least two sessions |
| PROVEN-ext | proof by an external model (Sol/GPT); the team read and checked it, but the text is in no repo |
| COMPUTED | exhaustive exact computation |
| DERIVED-here | new in this note; needs a second reader |
| PATTERN | observed on samples only |
| OPEN | not settled |

---

## A. The formal obligation

**Theorem C2 (the target).**
- *Setting.* Let Γ be a graph on 99 vertices. Its adjacency matrix A is symmetric, 0/1, with zero diagonal, and satisfies
  A² + A − 12I = 2J; that is, Γ is an srg(99,14,1,2).
- *Assumption.* Let σ be a permutation of the vertices with P_σ A = A P_σ, σ² = 1 and σ ≠ 1.
- *Claim.* No such pair exists.

- **Object assumed to exist:** the pair X = (Γ, σ).
- **Conditions it must satisfy:** the srg(99,14,1,2) axioms, plus σ ∈ Aut(Γ) of order 2. Nothing else.
- **What must be shown impossible:** the existence of X.
- **Why it matters.** The umbrella result is "if the 99-graph exists it is asymmetric". It follows by combining three inputs:
  - C2 (this obligation);
  - C3: no automorphism of order 3 (Lean, conway99-c3-experimental);
  - the literature: |Aut| ∈ {1, 2, 7, 11} by hand, {1, 2} with Behbahani–Lam (c3_experimental JOURNAL).
- **Not the target:**
  - nonexistence of Γ itself;
  - the trivial-Aut branch (the C1 branch).

---

## B. The implication chain and the dependency graph

### B.1 The shortest known chain

Every arrow below is a consequence that holds for EVERY X.

| step | statement | status |
|---|---|---|
| F1 | σ fixes exactly one vertex x and swaps the two ends of each of the 7 lines through x | PROVEN: Makhnev–Minakova 2004, Th. 1 (order 2); team proof in ROUTE_C_PROOF, Lemma 1.1 |
| F2 | The σ-odd orbit matrix S (49×49: 7 block orbits, 42 root orbits) satisfies S² + S = 12I, with spectrum 3²⁷(−4)²². In block form S = [[−I, C], [Cᵀ, T]], with CT = 0 and T² + T = 12I − CᵀC | PROVEN (ROUTE_C_PROOF, Lemma 1.2) |
| F3 | The frame lattice N (rank 22) is spanned by 42 u (norm 3) and 7 h (norm 4), with Gram 3I − S and Σvvᵀ = 7I | PROVEN (ROUTE_C_PROOF, Th. 2.1) |
| F4 | min N* ≥ 8/7: no integer (−4)-eigenvector b of S has \|b\|² ≤ 7 | PROVEN (handoff §10; re-read in notes/139) |
| F5 | The 7-part of disc(N) is non-split (B′), with d even | PROVEN, two proofs (C2's B′; handoff §9) |
| F6 | There is an envelope U ⊇ N ⊕ ⟨y₁, y₂⟩: odd unimodular, rank 24, min U ≥ 2 | PROVEN (verified independently by C1 and structures, notes/139) |
| F7 | Switch: U′ = ℤ⁷ ⊥ V, with V one of the 9 unimodular lattices of rank 17 and y₁, y₂ ∈ V. Then K = V ∩ y⊥ has rank 15 and det 49, w_u ∈ K, and ρ: E₇ → K/2K is injective | PROVEN (handoff Th. 3.1, Prop. 4.1 and 5.1; re-derived by structures and C2); the 9 lattices come from Conway–Sloane |
| F8 | (V, y₁, y₂) is one of 1,272 classes. The unit bound Y(e) ≤ 4 comes from PE and min U ≥ 2 | COMPUTED: the 8 ℤ¹⁷ classes by hand (handoff Cor. 6.2); the other 8 V in notes/146, cross-checked by two enumerators |
| F9 | Dictionary: take κ ∈ K* with 7\|κ\|² = n. Then a = ((κ, w_u))_u ∈ ℤ⁴² satisfies Ta = −4a, Ca = 0, \|a\|² = 4n, and a mod 2 is a cut of K₇ of size k. PE gives n ≥ k(7−k)/2. In graph terms, h = 2(κ, f_p) is an integer σ-odd (−4)-eigenvector of Γ that vanishes on x ∪ N(x) | PROVEN (PE: handoff Th. 6.1); DERIVED (notes/150 §1) |
| F10 | Rungs R_n for n ≤ 4: no such a exists | PROVEN-ext (Sol: "read and checked", GPT_FOLLOWUP_MINK) |
| F11 | 492 classes have a κ with n ≤ 4, so they are excluded. The other 780 have min n ∈ {5, 6, 7, 8} | COMPUTED (two independent scans, explicit witnesses) |
| F12 | 26 of those 780 are excluded by the exhaustive skeleton engine (Z4 plus integer DFS). 754 remain open | COMPUTED (notes/147, data/rank17/z4_*.jsonl) |
| F13 | **Every one of the 754 open classes has a K* vector with n = 6 or n = 8** | COMPUTED here (§D.3) |
| **gap** | **R6 and R8** | **OPEN** |

Here **R_n** is the statement: "no a ∈ ℤ⁴² satisfies Ta = −4a, Ca = 0, |a|² = 4n, with a mod 2 a cut". In team notation
this reads Na = 4a and Rᵀa = 0, with N = −T and R = Cᵀ. By F9, R_n forbids every κ ∈ K* with 7|κ|² = n, in all nine V at
once. It is **V-free**: a statement about Γ alone.

**Where the chain stops:** at R6 ∧ R8. Nothing else is missing mathematically (see §E3 for the rigour items).

### B.2 Dependency graph (✓ = established, ✗ = missing arrow)

```
THEOREM C2  (no involution)
└── ⊥  ⟸  [every X lands in one of the 1,272 classes]  ∧  [every class is closed]
    │
    ├── every X lands in a class  ✓
    │     F1 → F2 → F3 → F4, F5 → F6 → F7 → F8
    │     F1 ⟸ Makhnev–Minakova (lit.)   F4 ⟸ path counting (§10)   F7 ⟸ Conway–Sloane (lit.)   F8 COMPUTED
    │
    └── every class is closed:
          ├── 492 classes: K* has κ with n ≤ 4   ✓   F11 (COMPUTED) + F9 + R≤4 (PROVEN-ext)
          ├──  26 classes: engine                ✓   M1, M2, L1–L4, DC, Z4 (PROVEN) + integer DFS (COMPUTED)
          └── 754 classes: K* has κ with n ∈ {6, 8}   F13 (COMPUTED) + F9 + R6 + R8
                 R6 (n = 6) = J0 ✓ (F4) + J1b ✓ (star bound) + J1a ✗ + J2 ✗ + J3 ✗
                 R8 (n = 8) = k0 ✗ + k1 ✗ + k2 ✗ + k3 ✗
```

Missing arrows (✗):
- the open cells of R6 and R8;
- rigour arrows: PROVEN-ext → PROVEN for Sol's inputs, and a completeness certificate for F8.

**R5 and R7 are not on the graph.** Given R6 ∧ R8 they close nothing more (§D.3).

---

## C. Existing results and their relevance to Theorem C2

Relevance classes: **1** = directly advances the chain; **2** = eliminates a rigorously defined class; **3** = could become
relevant if an extra implication is proved; **4** = currently disconnected.

| result | status | rel. | how it connects / what is missing |
|---|---|---|---|
| F1–F9, F11, F13 (§B.1) | PROVEN / COMPUTED | 1 | the chain itself |
| R≤4: min K* ≥ 5/7 (Sol) | PROVEN-ext | 1 | closes 492 classes; the proof text must be archived (§E3) |
| **min N* ≥ 8/7 also closes every k = 0 cell with n ≤ 7** | DERIVED-here (one line, below) | 1 | closes J0 of R6 and k = 0 of R7; both are still listed as open in GPT_FOLLOWUP_MIN97 |
| Sol's star bound (star energy ≤ ⌊20n/7⌋) and star rules R1–R3 | PROVEN-ext (C1 checked the signs) | 1 | closes J1b; will carry much of R6/R8; needs archiving |
| §7.2: n = 5 ⇒ k = 2 | PROVEN-ext | 2 | concerns R5 only |
| 59 engine exclusions (26 beyond K*) | COMPUTED | 1–2 | 26 of them are needed until R6 ∧ R8 (+ R5 for 2 classes) replace them |
| R5 partial (notes/150): classes A–D; class A has β ∈ {12, 14, 16}; E_out gives β ∈ [8, 16] | DERIVED | 2 | eliminates sub-cells of R5. R5 is redundant, but the machinery transfers to J2 (R6) and to k = 2 of R8 |
| DC lemma: [7κ] ∉ span(skeleton) for n ∈ {2, 4, 6, 8, 10, 14} | PROVEN | 3 | the mod-2 shadow of R6/R8; matters only if the engine route resumes |
| Z4 layer (self-dual ℤ₄ codes 4⁶2¹²) | PATTERN; not universal | 3 | a uniform ℤ₄ lemma would only help the engine route |
| IB / class-energy inequalities (Sol) | PROVEN-ext | 3 | a tool for R6/R8 |
| Route (C): 34/156 envelopes excluded | PROVEN / COMPUTED mix | 2 → 4 | superseded. Each U maps to a (V, y) class, and no transfer lemma "all U over a class dead ⇒ class dead" is used |
| Block/pivot programme 26–27/09: all-sib dead (C2 Th. D, S39); Σ-motifs with ≥ 3 sib pairs; anchored family; g⁻ = 16 dead; Σ = ∅ one-block 25–33/79 (DRAT) | PROVEN / COMPUTED | 2 → 3 | these live in the D-structure partition (sib / one / dis). They are relevant if the R6/R8 proofs split by the D-structure of the cut stars (notes/149 §2: Q_c, P_c are perfect matchings and Q_c ∩ P_c = D-edges) |
| g⁻ ∈ {18, 20}, q:7 rigidity (notes/129–132) | PROVEN | 3 | relevant only if a rung proof must separate K = L (g = 20) from [K:L] = 7 (g = 18) |
| L_β (S40), even sector (notes/135) | DERIVED | 3 | relevant if R6/R8 need the σ-even half of Γ |
| Target (A): consistency of the V-free N axioms (structures) | OPEN | 3 | if the axioms are inconsistent, C2 follows directly. If a model carries an n ∈ {6, 8} witness, R6/R8 must use the σ-even part |
| n5_local_solution.json: a solution of the two-star relaxation carrying an n = 5 witness | COMPUTED (CP-SAT, before the ban) | 3 | not an elimination but a test object: any R5/J2 argument that only uses the two cut stars fails on it |
| C1 branch (trivial Aut): L1–L26, twisting, P ≤ 1151, 7-adic k bounds, Hodge, Krein, order-8 positivity | PROVEN / COMPUTED | 4 | about Γ without σ |
| Pass-1 PATTERN verdicts, "every sampled skeleton dies in star 1", Knuth estimates, engine speed-ups | PATTERN / measurements | 4 | not eliminations |
| SAT / CP-SAT UNKNOWN runs (n5_signed, alldis probe, β-cubes, 79-type cubes) | no verdict | 4 | carry no information |
| "No C₅-symmetric class-A model" (notes/150 §9) | DERIVED | 4 | assumes a symmetry the problem does not impose |
| Kills from diameter-3 witnesses or partial configurations (audit_refutations, 11 over-scoped labels) | re-scoped | 4 | outside the chain |
| "min K* ≥ 9/7 kills all 1,272" | correct but too strong | – | R5 and R7 are not needed (§D.3) |

**The k = 0 remark (DERIVED-here).**
- For κ ∈ K* we have (κ, h_c) = 0 and (κ, w_u) = 2(κ, u).
- So k = 0 (all a_u even) ⟺ (κ, u) ∈ ℤ for every u ⟺ κ ∈ N*, since the u generate N (handoff §8).
- Then F4 gives n = 7|κ|² ≥ 8. Hence k = 0 is impossible for n ≤ 7.
- In particular J0 (n = 6, b = a/2 with six ±1) and the case "k = 0, seven ±1" at n = 7 are closed. §10's counting gives
  Σd² + 7e ≤ m(m − 7) < 0 for m = 6 in one line.

---

## D. Open branches

### D.1 Exhaustive partition by (V, y) class: S = S₁ ∪ … ∪ S₁₂₇₂

| V | classes | closed by R≤4 | closed by engine only | **open** | min n over the open classes | open classes with n = 6 / n = 8 |
|---|---|---|---|---|---|---|
| ℤ¹⁷ | 8 | 5 | 0 | **3** | 5: 2, 6: 1 | 3 / 1 |
| E₈ ⊕ ℤ⁹ | 53 | 32 | 0 | **21** | 5: 18, 6: 3 | 21 / 13 |
| D₁₂⁺ ⊕ ℤ⁵ | 158 | 78 | 1 | **79** | 5: 63, 6: 16 | 78 / 62 |
| E₇²⁺ ⊕ ℤ³ | 191 | 82 | 1 | **108** | 5: 75, 6: 32, 7: 1 | 101 / 95 |
| A₁₅⁺ ⊕ ℤ² | 203 | 80 | 1 | **122** | 5: 75, 6: 45, 7: 2 | 105 / 115 |
| E₈² ⊕ ℤ | 12 | 7 | 5 | **0** | – | – |
| D₁₆⁺ ⊕ ℤ | 18 | 10 | 8 | **0** | – | – |
| D₈²⁺ ⊕ ℤ | 319 | 107 | 4 | **208** | 5: 111, 6: 80, 7: 17 | 151 / 194 |
| A₁₁E₆⁺ | 310 | 91 | 6 | **213** | 5: 68, 6: 62, 8: 83 | 96 / 196 |
| **total** | **1,272** | **492** | **26** | **754** | 5: 412, 6: 239, 7: 20, 8: 83 | 555 / 676 |

- The 754 open classes carry 287 distinct dual theta series (n ≤ 14), so they hold at least 287 distinct lattices K.
- The same K can arise from several (V, y).

### D.2 V-free partition of a dual witness by (n, k, pattern)

Throughout: |a|² = 4n; the odd labels of a are exactly the 2k(7 − k) labels of the cut.

| n | k | cell | status |
|---|---|---|---|
| ≤ 4 | all | – | PROVEN-ext (Sol) |
| 5 | 0 / 1 / 2 | k = 0 (F4); k = 1 (§7.2); k = 2: ±1 on the 20 labels of a 2\|5 cut, classes A–D | closed / closed / **OPEN, but not needed** |
| 6 | 0 | J0: b = a/2 with six ±1 | **closed by F4** (DERIVED-here) |
| 6 | 1 | J1a: 12 × ±1 on star(c) + three ±2 off it; J1b: 11 × ±1 + ±3 + ±2 | **J1a OPEN**; J1b closed (star energy 20 > 17) |
| 6 | 2 | J2: 20 × ±1 on a 2\|5 cut + one ±2 off it | **OPEN** |
| 6 | 3 | J3: ±1 exactly on the 24 labels of a 3\|4 cut; 8 sign classes (n6_k3_classes.json) | **OPEN** |
| 7 | 0–3 | k = 0 closed by F4; k = 1 with two ±3 closed by the star bound; rest listed in MIN97 | OPEN, **not needed** |
| 8 | 0 | eight ±1 on root orbits (a ±2 is impossible by §10 step 1) | **OPEN** (first reduction below) |
| 8 | 1 | 12 × ±1 + five ±2; 12 × ±1 + ±4 + ±2; 11 × ±1 + ±3 + three ±2 | **OPEN**; 10 × ±1 + two ±3 + ±2 closed (star energy 28 > 22) |
| 8 | 2 | 20 × ±1 + three ±2; 19 × ±1 + ±3 + ±2 | **OPEN** |
| 8 | 3 | 24 × ±1 + two ±2; 23 × ±1 + ±3 | **OPEN** |

- The n = 8 pattern lists are norm/parity bookkeeping (DERIVED-here); check them against Sol's list.
- The star bound prunes further. It applies to every star, not only the cut star.

**First reduction of R8, k = 0 (DERIVED-here).** Let P be the positive support, with |P| = 8.
1. §10's inequality Σd² + 7e ≤ m(m − 7) = 8 forces e = 0, so P is independent.
2. Each p ∈ P has exactly 4 neighbours in σP.
3. The P–σP pattern B is symmetric, 4-regular, with zero diagonal (p ≁ σp off N(x)).
4. B² = 2I + 2J − E, where E counts the common neighbours outside P ∪ σP. E is the adjacency matrix of a 2-regular
   multigraph (row sums 2).

Spectra alone do not kill this: B = C₈(1,2) with E = 2C₄ satisfies it. The kill has to use the root-layer (D₇) geometry:
- Ca = 0 at the blocks;
- equal numbers of neighbours in P and σP at every vertex outside P ∪ σP;
- λ = 1.

### D.3 Which rungs are needed (COMPUTED here from out_kstar_ladder.json)

"Rung n covers a class" means that the class's K* contains a vector with 7|κ|² = n, not necessarily its minimum.

| rung set | open classes covered (of 754) | of the 780 (not using the engine) |
|---|---|---|
| R5 | 412 | 438 |
| R6 | 555 | 569 |
| R7 | 274 | 274 |
| R8 | 676 | 696 |
| **R6 ∧ R8** | **754** | 778 |
| R5 ∧ R6 ∧ R8 | 754 | 780 |
| R5 ∧ R8 / R5 ∧ R6 / R5 ∧ R6 ∧ R7 | 716 / 651 / 671 | – |

1. **Every cover drawn from rungs 5–8 contains both R6 and R8.**
   - 24 open classes have no vector with n ∈ {5, 7, 8} but do have one with n = 6.
   - The 83 A₁₁E₆⁺ classes with min 8 have nothing below n = 8.
2. **R6 ∧ R8 is sufficient** for the 754.
   - R5 and R7 add nothing once R6 and R8 hold.
   - The earlier accounting "R5 = 438 pairs, the largest piece" counts minima only. Rung n actually kills every class with
     any n-vector.
3. The two classes that escape R6 ∧ R8 among the 780 are A₁₁E₆⁺ #31 and D₁₆⁺ #7.
   - They have identical dual theta series, so they are probably one K.
   - Both are engine-excluded.
   - An engine-free chain needs R5 (or R9) for them.
4. Minimal covers using rungs 5–10: {6, 8}, {6, 9}, {7, 10}, {8, 10}, {9, 10}.
   - {6, 9} is the fallback if R8 stalls.
   - Every open class has vectors with n = 12, 13 and 14. So the rungs are not lattice facts; they can only hold because of Γ.

---

## E. The three most valuable next directions

**E1 — Prove R6 (open cells J1a, J2, J3).**
- *If it succeeds:* it removes the 555 open classes whose K* has an n = 6 vector. R6 lies in every rung cover (24 classes
  have no other rung ≤ 8).
- *How to start:*
  - J3 first: a = ±1 exactly on a 3|4 cut, with 8 sign classes already enumerated. This is the same shape as C1's R5 analysis
    (classes, two-star identities, E_out, β-bounds), so that machinery transfers directly.
  - Then J2: the 2|5 cut of R5 plus one ±2. Here notes/150 applies almost verbatim, with a modified Ca = 0 balance.
  - J1a: star(c) carries twelve ±1 and three ±2 sit off it. The star bound and IB are the natural tools.

**E2 — Prove R8, starting with k = 0.**
- *If it succeeds:* it removes 676 of the 754 open classes by itself. Together with R6 it closes all 754, and with F1–F12
  that proves Theorem C2.
  - This proof still rests on the engine for 26 classes.
  - Adding R5 for A₁₁E₆⁺ #31 and D₁₆⁺ #7 makes it engine-free.
- *How to start:* k = 0 ("N* ∩ W has no vector of norm 8/7") is a hand-sized extension of an existing proof (§10). The
  reduction in §D.2 leaves one independent 8-set with a 4-regular B, to be killed by the D₇ geometry. Then treat k = 3, 2, 1
  with the star bound.

**E3 — Turn the PROVEN-ext and COMPUTED arrows the finish depends on into proofs.**
- *If it succeeds:* the open set does not shrink. It makes R6 ∧ R8 an actual proof rather than a proof modulo unchecked
  inputs.
- *Items:*
  - (a) Put Sol's proofs into this repo, each checked line by line by two sessions: R≤4 (carries 492 classes), §7.2, the star
    bound, R1–R3 and IB. R1–R3 and the star bound will also carry R6/R8.
  - (b) An independent completeness check of the 1,272-class list (F8), e.g. a mass-formula check on the unfiltered y-pair
    enumeration.
  - (c) Decide whether the 26 engine exclusions stand as COMPUTED or get replaced by rungs.
- This is not new mathematics, but without it the chain is not a proof.

**Deliberately not in the top three:**
- R5 as the main target (redundant given R6 ∧ R8; keep it only as the pilot for J2).
- R7 (redundant).
- Target (A) by hand construction (low expected yield without solvers; it is a diagnostic, see F).
- Multi-vector per-class arguments: these are V-dependent. They are the fallback if R8 stalls on the 83 min-8 A₁₁E₆⁺ classes.

---

## F. Evidence that the conjecture may be false or under-constrained

1. **Nothing resembles an actual X.**
   - No complete frame appeared in any search.
   - All 59 exhaustively searched classes died.
   - Every sampled skeleton died.
   - None of this is proof.
2. **Every relaxation tried is satisfiable.** So the obstruction, if it exists, is global and combinatorial. Examples:
   - the two-star local system with an n = 5 witness (n5_local_solution.json);
   - the B-level system of R8/k = 0 (C₈(1,2) / 2C₄);
   - 1,213 classes that are mod-4 feasible;
   - the whole genus / 7-adic layer, which is consistent: "no genus obstruction will kill C2" (C2 N09);
   - all moment identities, which are automatic (notes/149–150).

   Any proof of R6/R8 must use integrality and the exact λ = 1, μ = 2 beyond the cut stars.
3. **Uniformity in m is forbidden.**
   - Paley(9) (m = 2) and BvLS (m = 11) have involutions with one fixed point.
   - Both are all-sib, and all-sib is dead at m = 7. So a proof may use a non-sib orbit, but not an argument valid for every m.
   - C1 found no n = 5 analogue in BvLS (dimension 0).
   - The same exact linear-algebra check has NOT been run for the J/k cells of R6 and R8. It is cheap and should come before
     any investment: a cell whose analogue exists in BvLS needs an m = 7-specific argument.
4. **Difficulty risk.**
   - Restricted to the 754 classes, R6 ∧ R8 is equivalent to Theorem C2 there. It is easier only if small-support
     eigenvectors are locally rigid.
   - The evidence so far: R≤4 and F4 fell by hand. R5 (20-label support) resisted a morning, with every identity automatic.
     R6/J3 has 24 labels and R8 up to 32.
5. **The rank-15 frame problem is badly under-constrained on its own.**
   - It drops the σ-even half of Γ. Knuth estimates give 10⁴–10¹² skeleton leaves per class.
   - That is why the engine cannot finish, and why the rungs, which bring λ and μ back through Sol's star rules, are the right
     level.
6. **Bookkeeping drift.**
   - J0 and k = 0 of R7 are still listed as open although F4 closes them.
   - "min K* ≥ 9/7" asks for more than is needed.
   - None of this affects correctness, but it wastes effort.

---

### Reproduce §D.3

Plain arithmetic on the stored exact scans; no solver.

```
python - <<'EOF'
import json, itertools
d = json.load(open('conway99-rank17/out_kstar_ladder.json'))
o = set(json.load(open('c1_recherche/data/rank17/open_after_sol.json')))
lad = {'%s:%d' % (r['V'], r['id']): r for r in d}
cov = lambda keys, rs: {k for k in keys if any(lad[k]['counts'].get(str(n), 0) for n in rs)}
for rs in [(5,), (6,), (7,), (8,), (6, 8), (5, 6, 8)]:
    print(rs, len(cov(o, rs)), len(cov(set(lad), rs)))
EOF
```
