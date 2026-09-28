# Targets A / B — structures' notes (28/09, afternoon, by hand; no solver, no enumeration)

## 1. The complete axiom system (DERIVED)
A commutes with σ, so A² + A − 12I = 2J splits into an antisymmetric and a symmetric part. With the 42 labels u,
representatives p_u, N_uw = A(p_u,σp_w) − A(p_u,p_w) and M_uw = A(p_u,σp_w) + A(p_u,p_w):

- **(I) antisymmetric:** N² = N + 12I − RRᵀ, and NR = 0 (inner–outer block).
- **(II) symmetric:** M² + M + |R||R|ᵀ = 12I + 4J, and |R|ᵀM = 4J − 2|R|ᵀ (inner–outer block).
- **Entries:** M = |N| + 2D, where D is a 0/1 symmetric matrix disjoint from supp N.
  - The diagonals of (I) and (II) force deg N = 10, p ≁ σp, and D a perfect matching (the D-partners).

**(N, D) satisfying (I) + (II) ⟺ a graph with the involution.** Everything else (x, inner lines, inner–inner blocks) is
automatic.

So the N-axioms of Target A are exactly (I) plus the star consequences, and **the only missing graph axiom is (II)**:
- the symmetric part;
- equivalently, the separate λ/μ counts for (p_u, p_w) and (p_u, σp_w), not just their difference.

"D-consistency" (C1) is a piece of (II).

The rank-17 engine (SPEC_ENGINE §1–3) searches for exactly an N-frame: norm-10 vectors z_u with z_u·z_v = 4N_uv − r_u·r_v
and tight frame 28 P_W. So an explicit model of the N-axioms = a frame found by that kind of search in some rank-15 lattice.

## 2. Symmetric ansätze for an N-model — all dead (DERIVED)
| group (action on the 7 inner lines) | why it dies |
|---|---|
| S7 | invariant {0,±1} N is block-diagonal on the ± labels; spectrum 5/−2 instead of 4/−3 |
| A7, and anything containing L3(2) | via L3(2) below |
| L3(2), untwisted | star mates forced onto Fano lines; degree ≤ 8 |
| L3(2), twisted by sgn of S4 (the only other class, H¹ = F₂) | see below |
| S5, A5, F20 on the 5-side of the cut (graph level) | star-0 mates forced to cell-mates; star-1 row of (0,q,+) then lies in {o±}, which is zero, but needs 4 entries |
| C5, D5 in class A (C1) | β = 0, but class A needs β ≥ 11 |

Twisted L3(2), in detail:
- All 42 labels form one orbit. τ_k fixes each label up to sign −1, so every star block is 0.
- The outside rows then give 4a + 4b = 4 at star 3, while 2a + 2e = 4 and 2b + 2e = 4 at stars 4 and 5. So a = b and
  a + b = 1, a contradiction.

Point-transitive groups still open: C7, D7, F21, F42. Too many unknowns to settle by hand.

## 3. Target B at graph level (DERIVED)
Notation:
- cut lines 0, 1; T_c = the 20 cut vertices on line c;
- O = the 4 vertices of cell {0,1}; Z = the 40 vertices of cells inside Q = {2..6};
- X_c^± = the outer neighbours of c^±.

Facts:
- **Degree splits.** For v ∈ T_0: v has 2 − ι neighbours in T_0, 4 − ι in T_1, ι in O, and 6 + ι in Z.
  - For z ∈ Z: n_{T_0}(z) = n_{T_1}(z) = 4 − n_O(z).
- **J is balanced and takes three values.** Counting T–Z edges both ways gives J := #(T_0–O edges) = #(T_1–O edges).
  Also J = 8 − 2e(O), and e(O) ∈ {0, 2, 4} since M_{o₊o₋} ∈ {0, 1, 2}.
  - So J ∈ {8, 4, 0} and C1's o₁ = o₂ = J/2 ∈ {4, 2, 0}.
- **A quadratic identity.** From |A1_{T_0}|² and A² = 12I − A + 2J: Σ_O (2 − n_O)² + Σ_Z (4 − n_O)² = 360 + 11J.
- **J = 8.** Σ_Z n_O(n_O − 1) = 16 exactly equals the μ-bound, so every O-pair's outer common neighbours lie in Z.
  Two sub-configurations: k₂ = 8, or k₃ = 2 with k₂ = 2.
- **J = 0 (rigid).**
  - O is a 4-cycle and every z ∈ Z has exactly one O-neighbour, so Z = ⊔ Z_o with |Z_o| = 10.
  - Each T vertex has 6 T-, 6 Z- and 0 O-neighbours.
  - For each o ∈ X_c^α, z ↦ (its unique neighbour in T_c ∩ X_c^α) is a bijection Z_o → T_c ∩ X_c^α.
- **Dictionary with C1.** For v ∈ T_0, β_U = −h(v)·(h(m_0) + h(q_0)) ∈ {0, 2, 1} for types A, B, C, so
  β = 2·#B + #C (labels).
  - In class A, Σ_{z ∈ Z⁺⁺} (T_0⁺-nbrs − T_0⁻-nbrs) = β − 20 = 2σ₊. This reproduces C1's identity.
- Every identity above that is linear in h is automatic (as C1 found). What can still bite is (II): D-partners, O-adjacency,
  and the separate λ/μ counts.

## 4. Edge eigenvectors (DERIVED)
- For every edge uv with third triangle vertex w, y_uv = −4(e_u − e_v) + 1_{N(u)∖{v,w}} − 1_{N(v)∖{u,w}} lies in E₋₄.
  Its support has size 26.
- The block vector b_c is y for the edge c⁺c⁻, and ⟨h, y_uv⟩ = 7(h(v) − h(u)).

## 5. R8, k = 0 (structures' lane after the 28/09 re-split; DERIVED, by hand)
Setting (C1's reduction):
- P = 8 outer vertices with h' = +1 and σP with h' = −1; P is independent.
- B_ij = A(p_i, σp_j) is symmetric, 4-regular with zero diagonal, and B² = 2I + 2J − E. Here E_ij = the number of common
  neighbours of p_i, p_j outside P∪σP, so E is 2-regular.
- The roots ρ_i (the signed inner neighbours of p_i) sum to 0, and t_z ≤ 2 by the star bound.

**5.1 Only two B's.** (B²)_ij ≤ 2 for i ≠ j.
- For B_ij = 0 this forces λ = 0 in the complement.
- For B_ij = 1 it forces μ ≤ 2 in the complement.
- So the complement is a triangle-free cubic graph on 8 vertices, i.e. the cube Q₃ or the Wagner graph V₈. Hence:
  - **B₁ = K₄□K₂** (complement of Q₃): E = 2·(the K₂-matching), i.e. p_i and its matching mate p_i* share 2 outside
    neighbours;
  - **B₂ = C₈(1,2)** (complement of V₈): E = C₈(2), two 4-cycles with multiplicity 1.

**5.2 Connectors.** Let t(y) = #P-nbrs = #σP-nbrs of an outside vertex y. Σt(t−1) = 16, and E forbids t ≥ 3.
- Every t = 2 vertex y has P-nbrs forming an E-edge e and σP-nbrs σ(f) for another E-edge f ≠ e.
  - f = e would put y in the D-orbits of two different labels, which contradicts D being a matching.
  - Call y a connector e → f. σ maps e → f to f → e, and each E-edge is the P-side of exactly E_e connectors.
- **B₁.**
  - The connector graph on the 4 pairs M_a is a 4-cycle.
  - Each connector is the apex of the two same-layer H-edges between its pairs.
  - Matching H-edges p_iσp_i* and diagonal-pair H-edges have t = 1 apexes.
  - A t_z = 2 line is an inner connector pair on one C₄-edge. Two such lines at the same M_a would force ρ_i = ρ_i*, so
    w := #(t_z = 2) ∈ {1, 2}, and for w = 2 the two lines sit on opposite C₄-edges.
- **B₂.** The connector map is a fixed-point-free involution φ on the 8 E-edges.
  - An allowed partner of {0,2} is one of {4,6}, {3,5}, {5,7}. This is forced by N(y) being a matching: B_ij + B_ij' ≤ 1.
- **In both cases:** for c = c_z^s, the outside vertices of X_c have Σt = 16 − 6t_z.

**5.3 w = 1 (DERIVED; both B's).** Let w := #{z : t_z = 2}.
- **Setup.** t ≤ 2 everywhere, since E has no triangles. Σt(t−1) = 16, so exactly 8 vertices have t = 2, of which 2w are the
  inner c_z^± with t_z = 2.
- **A t_z = 0 line kills w ≥ 2.** Take a line ℓ0 with t_z = 0 (one exists iff w ≥ 2).
  - Every support vertex p has μ(p, c_ℓ0^±) = 2, with both common neighbours in X_ℓ0^±.
  - So Σ_{y∈X_ℓ0^+} t(y) = 16 over 12 outer vertices. That needs ≥ 4 outer t = 2 vertices in X_ℓ0^+, and likewise in the
    disjoint X_ℓ0^−.
  - That is 8 ≤ 8 − 2w, a contradiction.
- **Hence w = 1.**
  - Exactly one line z* carries 4 support roots (2 of each sign); the other 6 lines carry 2 each (opposite signs).
  - There are 6 outer connectors and exactly 10 outer t = 0 vertices.
- **Refinement.**
  - At c_z*^± the 8 outside vertices of X have #t=0 = 4 + #t=2.
  - At each t = 1 inner vertex the 10 outside vertices have #t=0 = #t=2.
  - With r = #(t=0 vertices off star z*), r ∈ {0, 2}, and each X_z*^± holds exactly 1 − r/2 outer connectors.
- **B₁ specifics.**
  - The z*-line hosts one C₄-edge {M1, M2}. The outer connectors 1→4, 4→1, 2→3 and 3→2 avoid star z*.
  - Π (the pairs sharing a t=1 line) excludes M1–M2 pairs and same-layer C₄-adjacent pairs.
  - The only double edges of the root multigraph are M1 or M2 at z*.

## 6. R8, k = 0 is impossible (DERIVED, structures, 28/09 evening). Proof.
**Setting.**
- b ∈ {0,±1}⁴² with 8 nonzeros (a ±2 is impossible: its row needs absolute sum 8 > 4), Nb = 4b, Rᵀb = 0.
- Lift to the graph: h' = +1 on P (8 outer vertices), −1 on σP, 0 elsewhere; Ah' = −4h'.
- For y ∉ S := P∪σP, put t(y) = #P-nbrs = #σP-nbrs. Let t_z = the number of support roots with sign + at line z.

**G1.** P is independent and each p has 4 neighbours in σP (§10 count: Σd² + 7e ≤ 8 ⟹ e = 0). The slack is
Σ_{y∉S} t(t−1) = 16.

**G2. B ∈ {K₄□K₂, C₈(1,2)}; t ≤ 2; "connectors".**
- B ∈ {K₄□K₂, C₈(1,2)} is §5.1. Hence t ≤ 2, and exactly 8 vertices have t = 2.
- Each t = 2 vertex has P-nbrs an E-edge g and σP-nbrs σ(g') for an E-edge g' ≠ g; g' = g would give D two partners.
- σ maps g → g' to g' → g.

**G3. w = 1.** Proved in §5.3: exactly one line z* has t = 2 and the other six have t = 1.

**G4. Counts.**
- For y ∉ S: n₁' + 2n₂' = 16 − 5t and n₀' = 3t − 2 + n₂', where n_j' = #nbrs with t = j. This is λ/μ counting against P.
- There are 10 outer t = 0 vertices, and T0 = {x} ∪ those 10.
- Double counting T0–T2 edges gives e(T0) = 5 + e(T2) = 8 − r + e_oo.
  - r = #(outer t=0 vertices off star z*).
  - e_oo = #edges among the 6 outer connectors.
  - n₂(X_z*^±) = 1 − r/2 by σ, so r ∈ {0, 2}.

**G5. Cauchy–Schwarz in E₋₄.**
- r̂ := P₋₄1_S/2 ∈ E₋₄ has |r̂|² = 4/9 and ⟨1_T0, r̂⟩ = 11·8/63. Also |P₋₄1_T0|² = (418 − 18e(T0))/63, where
  P₋₄ = (27I − 9A + J)/63.
- So e(T0) ≤ 7, hence **r = 2 and e_oo ≤ 1**: no outer connector lies in star z*.

**G6.** Each outer connector y has 2 + 2 nbrs in X_z*^±. Of these, v(y) := |g∩(e∪f)| + |g'∩(e∪f)| are support, where e and f
are the ± pairs at z*. Σv = 8, so there are 16 connector edges into the outside part of star z*.

**G7. Only two t=0 vertices lie off star z*.** Write Y = the 8 outer t=0 vertices in star z*, τ, στ = the 2 off it, and
c(·) = #outer-connector nbrs.
- Σ_T0o n₂' = 2e(T0o) + 20 = 32 + 2e_oo, while the left side is 8 + Σ_Y c + 2c(τ).
- c(τ) ≤ 4: two full σ-pairs of connectors at τ would be 4 common nbrs of τ and στ.
- Hence e_oo = 0, c(τ) = 4, and every connector edge into star z* lands in Y.

**G8. v(y) ≤ 2, and at most two connectors have v = 2.**
- Every outer connector has n₀' = 4, with 4 − v(y) of those in Y. So it has exactly v(y) nbrs in {τ, στ}.
- The connectors with v = 2 are common nbrs of τ and στ, so there are at most 2 of them.

**G9. Cases.**
- **K₄□K₂.** The hosted C₄-edge {M_a, M_b} makes the outer connectors a→d, d→a, b→c, c→b have v = 2. That is four ✗.
- **C₈(1,2), e = {0,2}, f = {4,6}.** The edges {2,4} and {6,0} each carry v ≥ 2 on both connectors ✗.
- **C₈(1,2), e = {0,2}, f = {3,5}** (f = {5,7} is the mirror image).
  - (2,1,1) leaves four matchings φ, each containing an overlapping pair at root 6 or 7. That pair's two connectors would be
    the D-orbit of u₆ or u₇.
  - But X_z*^+ has 4 t=1 vertices matching the P-sources {p0,p1,p6,p7} with the σP-sources {σp4,σp5,σp6,σp7}.
    - p0 and σp5 have outside P_z*-partners w0 and w5, so λ forces w0 ↔ (0,4) and w5 ↔ (1,5).
    - Then (6,7),(7,6) would give λ(p7,σp6) = 2 ✗.
    - So (6,6),(7,7): D(u6) and D(u7) are t=1 vertices, contradicting the overlap ✗.  ∎

Consequence: the k = 0 cell of R8 is closed. Remaining R8 cells, after C1's row-capacity kills:
- k = 1: 12 × ±1 + five ±2;
- k = 2: 20 × ±1 + three ±2;
- k = 3: 24 × ±1 + two ±2.
