# evenframe — the even-sector class route at the hard core (structures, 01/10/2026)

Status: **work in progress, nothing closed yet.** The reduction and Lemma A are hand derivations that need a second
reader. The solver is cross-checked only in its first version (see Validation).

## 1. The idea
The new session's hard core is (g, k) ∈ {(20,26), (20,28)}, where the mod-7 two-star search has no grip. At k = 28 the
even sector carries no 7-adic information at all, but for the same reason its integral lattice is as small as it can be.

- L_β = ℤ⟨β_U⟩ is the even frame lattice of S40 (c99_structures_lab/NOTE_S40): 42 vectors of norm 5 in ℝ¹⁵, Gram
  6I − 2B − WWᵀ + J with B = X + 2D (Layer Q), tight frame Σββᵀ = 14·I, star sums Σ_{U∋c} β_U = 0.
- det L_β = 2⁸·7^{15−g′} with g′ = k − 13 (S40 §2, re-derived 01/10). So **k = 28 ⟺ det L_β = 2⁸**.
- L_β is 2-elementary and every vector of L_β* has integral norm (M = ℤ⁴² ∩ V₋₄(B) is even). So the discriminant form is
  alternating, it has a Lagrangian, and **L_β lies with index 16 in a unimodular lattice U of rank 15, with 2U ⊂ L_β**.
- There are five such U: ℤ¹⁵, E₈⊕ℤ⁷, D₁₂⁺⊕ℤ³, E₇²⁺⊕ℤ, A₁₅⁺.
- At k = 26 the same construction gives a rank-15 lattice of determinant 49, i.e. the genus of the K's of the rank-17
  route (the 1,272 classes). Not started.

## 2. The mod-2 skeleton (DERIVED)
- I := 2L*/2U ⊂ U/2U is 4-dimensional, totally isotropic and doubly even (norms ≡ 0 mod 4). C := I^⊥ = L/2U.
- Let u₀ be the characteristic class of U. W₆ := (C ∩ u₀^⊥)/I with q(v) = |v|²/2 mod 2 is an O⁺₆(2) space: its 28
  non-singular vectors are the pairs of 8 "points".
- The Gram mod 2 is J + WWᵀ, so the two orbits of cell {c,d} have class u₀ + ω_cd + I, where the ω are the 21 pairs of
  7 of the 8 points. A skeleton is (I, omitted point); the labelling of the 7 lines is then free (S₇ symmetry of Layer Q).
- Roots of L have classes in W₆ (never in I), so each root is a pair of points and its reflection swaps them. Omitted
  points joined by roots give equivalent skeletons.

## 3. Lemma A (hand proof): L_β has no vector of norm 1
For e ∈ L_β with |e|² = 1 put a_U = (e, β_U) ∈ ℤ. The tight frame gives Σa_U² = 14. Mod 2, a_U ≡ 1 + (w, ω_p) where
w ∈ W₆ is the class of e minus χ. Counting the odd entries: w = 0 gives 42; w a pair gives 22 or 30; w a 4-set gives 18.
All exceed 14. ∎
(The argument only uses the Gram mod 2 and the frame identity, so it holds for every g′.)

## 4. What the computation gives so far
- `evenpipe.py`: the five lattices, shells, classes mod 2, and all I up to root transvections (lexmin flags).
  Counts of I-representatives: ℤ¹⁵ 397, E₈ℤ⁷ 239, D₁₂ℤ³ 116, E₇²ℤ 115, A₁₅⁺ 221. Every one has W₆ of plus type.
- `Linv.py`: the 1,088 representatives fall into **31 invariant classes of L** (shell sizes and root system).
  18 of them contain a norm-1 vector and are dead by Lemma A. The 13 others, by root system of L:
  ∅, A₁⁴, A₁⁸, A₁²A₂², A₁⁶A₃, A₃², A₁⁴A₃², A₃+rank 4, A₁⁵+rank 5, D₄², A₇, D₄D₅, D₈.
- The D₈ class is L = √2E₇* ⊕ D₈ and is dead by a dimension count: every norm-5 vector has E₇-part of norm 3, and
  42·3 ≠ 14·7.
- `evenframe.c`: exact frame search per skeleton. Skeletons whose L splits die at the root on such dimension counts
  ("budgets": Σ_U |P_F β_U|² = 14·dim F for any subspace F). In the partial runs: ℤ¹⁵ 130 of 143 skeletons excluded,
  E₈ℤ⁷ 85 of 102; the rest hit the node cap.
- **Open: the indecomposable L's.** Example A₁₅⁺ skeleton 0 (L ⊃ A₇, all 28 pairs are roots): 6·10⁷ nodes without
  leaving the first depth-4 branch, even with the S₇×±1 lexmin test. Plain search will not close these.

## 5. Next step (not done)
Two-level search: split each orbit's candidates by their projection on the span of the roots of L ("types"), solve
the tight-frame equation on that subspace first, then lift. For L ⊃ A₇ the types of cell {c,d} are, with τ = 4t ∈ ℤ⁸:
±(3,3,−1⁶), ±4(e_c − e_d), ±(2,2,6@z,−2⁵); the count n₃ of the last kind is 8 or 12 (hand computation from the diagonal
of Στ τᵀ = 224I − 28J). A first naive type-level search (`a7types.c`) did not finish in 10 minutes; it needs the
arc-consistency propagation.

## 6. Validation
- `validate.py`: independent Python count of valid partial assignments against the C depth profile. MATCH on ℤ¹⁵
  skeletons 3, 8 (depth 3) and 40 (depth 2), for the solver with pairwise rules, line counts, PSD, sister order and
  count-aware forward checking.
- **Not yet re-validated:** budgets, the common-neighbour bound, reflection lexmin, group lexmin, fail-first ordering.
  All are meant as necessary conditions, and verdicts were unchanged where both versions finished, but no EXCLUDED
  verdict from the later versions should be quoted before that check is done.
- No positive control exists (no Layer-Q solution is known), and L_β exists only at m = 7, so BvLS/Paley say nothing here.
- Two skeletons returned a solver error (ℤ¹⁵ 384, E₈ℤ⁷ 257); not investigated.

## 7. Dead end recorded
The 42 A₇-orbit ovoids of PG(3,2) (weight-5 words of the [15,11] Hamming code) looked like a candidate frame in ℤ¹⁵.
Their intersection parities are 5 odd / 36 even per ovoid, against 21 / 20 required. Not a frame (`pg32.py`, `ovoid42.py`).

## 8. Hand mathematics, 01/10 afternoon (DERIVED, no computation; needs a second reader)
Notation: N := √2·L_β* (even, rank 15, det 2⁷, defined at k = 28), so L_β = √2·N*; x_U := β_U/√2 ∈ N*, norm 5/2.

**Lemma B.** For every root n of N and every orbit U, (n, x_U) ∈ {0, ±1}, and exactly 14 of the 42 values are ±1.
- Proof. The frame gives Σ_U (n, x_U)² = 7|n|² = 14, and Cauchy–Schwarz gives |(n, x_U)| ≤ 2. If (n, x_U) = ±2 then
  x_U ∓ n is a vector of N* of norm 1/2, i.e. √2 times it is a norm-1 vector of L_β, against Lemma A. ∎

**Dictionary between the two decompositions (at k = 28).** With Z_p := √2(β_{p+} + β_{p−})/2 ∈ N and ρ_Q := β_U + β_DU ∈ L_β:
- cell p is sibling ⟺ Z_p = 0; sisters X-adjacent ⟺ Z_p is a root of N; otherwise |Z_p|² = 4;
- a D-square is sibling ⟺ ρ_Q = 0; one-block (anchored or twisted) ⟺ ρ_Q is a root of L_β; disjoint ⟺ |ρ_Q|² = 4.
So a lattice class with no roots in L_β has only sibling or disjoint D-squares. This is the first map between C2's
D-type tree and a lattice classification.

**The rootless class is a problem about PG(3,2).** Here L_β = Λ(Hamming [15,11,3]) ⊂ ℤ¹⁵ and N ⊃ A₁¹⁵.
- The norm-5 vectors of L_β are the 168 projective frames ("ovoids") of PG(3,2) with signs.
- Each ovoid is an ovoid of the generalized quadrangle W(2) of exactly one of the 28 symplectic forms on F₂⁴; a form has
  6 ovoids, two of them meet in exactly one point, and each point lies on 2 of them.
- The 28 forms are the pairs of the 8 points; the seven lines are seven forms σ_x with Σσ_x = 0 and all σ_x + σ_y
  non-degenerate; cell {x,y} has form σ_x + σ_y. Explicit model: F₂⁴ = F₂ ⊕ F₈,
  σ_x((t,u),(t′,u′)) = Tr(ζ^{−x}(tu′ + t′u)) + Tr(ζ^{−3x}(uu′² + u²u′)), with Pfaffian Tr(ab²).
- Consequences for a frame in this class:
  1. both orbits of a cell are the two ovoids through one point i_p, so the unsigned frame is a map
     p ↦ i_p from the 21 cells to the 15 points; the signs at i_p are opposite;
  2. no cell is sibling and no sisters are adjacent (κ_p = 0 for all p);
  3. every D-square is of disjoint type: the two ovoids meet in 3 points with all three signs opposite;
  4. star relation mod 2: Σ_{p∋c} ω_p(i_p, ·) = 0 for every line c;
  5. point regularity: 2·#{p : i_p = x} + #{p : ω_p(i_p, x) = 1} = 14 for every point x;
  6. counting 5 in two ways: the number of pairs of cells with the same point equals the number of pairs of cells with
     the same plane i_p^⊥.
- Not yet shown: that no such map with signs exists.
