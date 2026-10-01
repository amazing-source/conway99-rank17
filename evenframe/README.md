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
  1. **(corrected 01/10 evening)** the two orbits of a cell are either two different ovoids of its form, meeting in one
     point i_p with opposite signs there (κ_p = 0), or the SAME ovoid with 3, 4 or 5 signs changed (κ_p = 0, 1, 2).
     The first version of this list assumed the two ovoids were always different; that was wrong.
  2. so sibling cells and adjacent sisters are possible; they are exactly the same-ovoid cells with 5 or 4 sign changes;
  3. every D-square is sibling or of disjoint type (two ovoids meeting in 3 points with all three signs opposite);
     one-block D-squares cannot occur;
  4. star relation mod 2, over the cells with two different ovoids: Σ_{p∋c} ω_p(i_p, ·) = 0 for every line c;
  5. point regularity: every point lies in exactly 14 of the 42 ovoids.
  (Items 4–6 of the first version were stated for different-ovoid cells only and are superseded by 4–5.)
- Not yet shown: that no such map with signs exists.


## 9. Hand mathematics, 01/10 evening (DERIVED; the CHECK lines are small exact verifications of the hand claims)

**Lemma C (minuscule projections).** Let R be an irreducible component of the root system of N. For every orbit U the
projection of x_U on span(R) is zero or a minuscule weight of R, and these projections form a tight frame:
Σ_U π_R(x_U) π_R(x_U)ᵀ = 7·Id on span(R).
- Proof. By Lemma B, (n, x_U) ∈ {0, ±1} for every root n; a weight with all root pairings in {0, ±1} is zero or minuscule.
  The frame identity restricts to any subspace. ∎
- Minuscule norms: A₁ 1/2; A₃ 3/4, 1; A₅ 5/6, 4/3, 3/2; A₇ 7/8, 3/2, 15/8, 2; D_m 1 (vector), m/4 (spinor); E₇ 3/2.
- For D_m ⊂ ℤ^m: every axis carries 7 − n_s/4 vector weights, n_s ≡ 0 (mod 4), and the spinor signs form m mutually
  orthogonal ±1 columns of length n_s.

**The 13 classes without norm-1 vectors (lattice data, `RN.py`):**

| # | roots of L_β | roots of N | rank of R(N) |
|---|---|---|---|
| 1 | none | A₁¹⁵ | 15 |
| 2 | A₁⁴ | A₁A₃³ | 10 |
| 3 | A₁⁸ | A₁¹¹D₄ | 15 |
| 4 | A₁²A₂² | A₁²A₃A₅ | 10 |
| 5 | A₁⁶A₃ | A₃²A₅ | 11 |
| 6 | A₃² | A₁A₃D₅ | 9 |
| 7 | A₁⁴A₃² | A₁⁷D₄² | 15 |
| 8 | rank 3 + rank 4 | A₁³A₇ | 10 |
| 9 | A₁⁵ + rank 5 | A₅D₅ | 10 |
| 10 | D₄² | A₁⁹D₆ | 15 |
| 11 | A₇ | A₁D₇ | 8 |
| 12 | D₄D₅ | A₁⁵D₄D₆ | 15 |
| 13 | D₈ | A₁⁸E₇ | 15 |

- **Class 13 is dead by arithmetic:** the E₇ budget needs (3/2)·n = 7·7, and 98/3 is not an integer. (Same kill as the
  dimension count of §4, now a one-line proof.)
- No other class dies on the norm budgets alone (checked by hand for A₁, A₃, A₅, A₇, D₄, D₅, D₆, D₇).

**Rootless class (#1), verified facts (CHECK, `z7_H.py`):** the seven forms σ_x in the F₂ ⊕ F₈ model have zero sum,
non-degenerate pairwise sums and degenerate triple sums; every cell form has 6 ovoids, two of them meet in one point;
ovoids of cells sharing a line meet in 0 or 2 points (1 : 5), of disjoint cells in 1 or 3 points (4 : 2); for an ovoid O
and an adjacent form, one ovoid is disjoint from O and the other five cut O in the edges of a pentagon.

**No 7-fold symmetric frame in the rootless class (CHECK, `z7_H.py`).** With the element of order 7 acting on the lines,
no choice of six signed base ovoids is even pairwise compatible with its rotations. (An invariant frame can be taken
invariant under the pure permutation, since H¹(ℤ₇, F₂¹⁵) = 0.)

**Not done.** Whether the rootless class has a frame at all is open. A star-level count did not finish (Python version
stopped by the time limit, C version cut at 15 minutes), so I do not know if even one star exists there.
