# PARADIGM_BREAK — representations of the involution case (01/10/2026)

Setting: Γ = srg(99,14,1,2), σ an involution fixing exactly x. Units = the 14 neighbours of x (7 lines of mates),
roots = the 84 other vertices = non-mate pairs {a,b} of units, orbits = the 42 σ-pairs of roots.
N = natural odd layer (N_UV = a(u,v) − a(u,σv)), B = X + 2D even layer, X = |N|.
Status words: DERIVED (by hand), CHECKED (machine, on a control), COLLAPSED (same content as existing work), OPEN.
None of the records below implies anything new about the hard core (g,k) ∈ {(20,26),(20,28)}.

---

## 1. Wedge flow
- PRIMITIVES: the 14 units; for each unit u a perfect matching μ_u of the 12 units off u's line.
- DEFINITION: μ_u(v) = w iff the roots {u,v}, {u,w} are adjacent (their triangle has apex u).
  Flow Φ(u,v) = (v, μ_v(u)) on the 168 ordered non-collinear pairs of units.
- MAKES SIMPLE: flatness at x is one sentence (μ_u = mate map for all u); the edges inside the far set split into a
  2-factor (wedge edges, orbits of Φ) and 140 triangles (a configuration 84₅,140₃).
- LAWS: (i) Φ is a permutation, its cycles have length ≥ 4 (a 3-cycle puts an edge in two triangles);
  (ii) flat ⟺ every cycle has length 4 ⟺ 21 quadrangles; (iii) μ_{−u} = σ μ_u σ.
- COUNTEREXAMPLES: the μ_u alone do not determine Γ (the 140 triangles are extra data).
- RELATION: Φ-cycles on orbits = the μ-edges of the seven stars; "adjacent sisters" = Φ-steps that land on a mate.
- STATUS: COLLAPSED. This is C1's shadow connection (notes/10 F1, notes/99).

## 2. Point reflections
- PRIMITIVES: for every vertex v of any (1,2)-graph, the involution τ_v: fixes v, swaps mates in N(v), sends the far
  vertex {a,b} to {a′,b′}.
- DEFINITION: x ▷ y := τ_x(y).
- MAKES SIMPLE: C2 ⟺ some τ_x is an automorphism (σ = τ_x; a non-trivially permuted line would fix a far vertex).
- LAWS: (i) x▷x = x, x▷(x▷y) = y; (ii) y ~ x ⟹ x▷y = y▷x = third point of the line;
  (iii) in the flat controls every τ_v is an automorphism (inversion through v).
- COUNTEREXAMPLES: self-distributivity fails in general; it holds iff τ_x ∈ Aut.
- RELATION: none needed; lattice language starts after σ is given.
- STATUS: COLLAPSED (commutation of reflections = grid closure, C1 notes/104).

## 3. Two pairings on 42 orbits
- PRIMITIVES: the 42 antipodal root pairs of D₇; pairing P (sisters, fixed); pairing D (σ-invariant quadrangles).
- DEFINITION: D pairs U with the orbit formed by the two common neighbours of u and σu. The 21 invariant
  quadrangles partition the 84 roots. Type of a D-pair: sister / one common line / disjoint; counts d₀+d₁+d₂ = 21.
  Tips of a quadrangle {±r,±r′}: third points of the triangles on r r′ and r(−r′).
- MAKES SIMPLE: D is intrinsic (no equation needed to define it); signed triangle counts.
- LAWS: (i) diag N³ = −2(n−2) = −10, so t⁺ − t⁻ = −5 at every orbit; (ii) t⁺(U) = 5 − (2−s_U) − tip(U);
  (iii) T⁻ = 140 − d₁ − 2d₂, T⁺ = 70 − d₁ − 2d₂; (iv) the number of σ-antipodal root hexagons is 140 whatever D is.
  (i) and the flat case of (iii) CHECKED on srg(243,22,1,2) (bvls_check.py: T⁺ 330, T⁻ 660).
- COUNTEREXAMPLES: parity counting of invariant 6-cycles gives nothing (the count is constant: 42+42+140).
- RELATION: the three D-types are ρ_Q = 0 / root of L_β / norm 4 in the even-frame dictionary.
- STATUS: COLLAPSED (trace identities of N; same content as the layer equations).

## 4. Star as three involutions
- PRIMITIVES: the 12 flags of a line c (roots containing +e_c); ν (other sister), μ_c (wedge), ψ_c (cross-mate).
- DEFINITION: a 3-coloured cubic multigraph on 12 flags; the star block of N is P_μ − P_ψ.
- MAKES SIMPLE: the star block as a signed sum of two matchings; coupling ψ_c fixes e_d ⟺ μ_d(e_c) = −e_c.
- LAWS: (i) ψ_c is a fixed-point-free involution; (ii) νψ_c is the holonomy of the fixed triangle {x,e_c,−e_c};
  its fixed points are the adjacent sister pairs on line c.
- COUNTEREXAMPLES: det(5 − 4ν + 2μ − 2ψ) is not a function of the cycle data of the three involutions
  (gem_det.py: 1,003 keys, only 489 with a single value; [det ≡ 0 mod 7] constant on 738 only).
- RELATION: star blocks of N; C1's line holonomy ρ_ℓ.
- STATUS: COLLAPSED, and the one prediction tested was false.

## 5. Triangle chart
- PRIMITIVES: a triangle T; the 36 vertices with one neighbour on T; the 60 with none.
- DEFINITION: the 36 carry a cubic graph (one mate edge, two cross edges); f(T) = number of prisms on T =
  fixed points of the holonomy of T.
- MAKES SIMPLE: counts of triangles by number of edges to T.
- LAWS: n₃ = f, n₂ = 6(n−1) − 3f, n₁ = 6(n−1)(n−3) + 3f, n₀ = n(2n²+1)/3 − 1 − 9(n−1) − 6(n−1)(n−3) − f.
  For n = 7: n₀ = 32 − f. CHECKED on n = 11 (300, 540, 0, 20). The linear data (E₁+2E₂+3E₃ = C² − 5C − 18) does
  not see f.
- COUNTEREXAMPLES: interlacing on T ∪ 36 is vacuous for every holonomy.
- RELATION: at a fixed triangle the chart is star c plus the D₆ roots.
- STATUS: COLLAPSED (C1: prisms, P = 1386 − Σ|W|/6, notes/100 L15; Hodge form notes/10 F21).

## 6. Lorentzian similitude
- PRIMITIVES: M = 2A + I (entries 0,1,2), M² = 49I + 8J.
- DEFINITION: S = [[M, 2·1],[4·1ᵀ, 29]], H = diag(I₉₉, −2). Then S H Sᵀ = 49 H.
- MAKES SIMPLE: Γ as a multiplier-49 self-similitude of the lattice I₉₉ ⊕ ⟨−2⟩.
- LAWS: (i) the identity above (hand); (ii) rank₇(S) is the 7-rank of A − 3I.
- COUNTEREXAMPLES: rational existence is unobstructed (multiplicities are integral).
- RELATION: it is the 7-rank language again; (g,k) are its σ-components.
- STATUS: COLLAPSED (p-rank).

## 7. Schur squares of the σ-blocks
- PRIMITIVES: even and odd blocks of the idempotents on roots: 63E⁺₋₄ = 27I − 9B + 2J, 63E⁻₋₄ = 27I − 9N.
- DEFINITION: entrywise squares and products of these blocks, each positive semidefinite separately.
- MAKES SIMPLE: inequalities on X and D alone.
- LAWS: 729I + 81X ⪰ 0, 837I + 4J + 45X + 252D ⪰ 0, 1168I − 8J − 41X − 280D ⪰ 0.
- COUNTEREXAMPLES: all three are far from tight (λ_min(X) ≥ −9 against the trivial −10).
- RELATION: Krein-type; entrywise it is only X = N∘N.
- STATUS: COLLAPSED and vacuous.

## 8. Curved (1,2)-geometry (counterfactual world)
- PRIMITIVES: any connected graph with every edge in one triangle and two common neighbours at distance 2;
  "flat" = wedge matching equals the mate map everywhere.
- DEFINITION: in a Cayley graph on order-3 subgroups, the braiding r(d,e) = the other factorisation of de;
  ρ_d(e) = first letter of r(d,e); wedge matching μ_d = ρ_d ∘ inverse ∘ ρ_d⁻¹; triangle holonomy = ρ_d⁻³.
- MAKES SIMPLE: the question "does curvature exist anywhere?" independent of 99.
- LAWS: (i) flat ⟺ ρ_d commutes with inversion; (ii) infinite curved examples exist by free completion
  (adding a new vertex never breaks the ≤ conditions), so flatness is not a local theorem.
- COUNTEREXAMPLES / DATA: exhaustive over connection sets (curved_cayley.py, curved2.py; A₆ and the three groups
  of order 243 capped): Heis(27), Heis×3, three groups of order 243, p:3 for p = 7,13,19,31,37,43, A₄, 2⁴:3, 5²:3,
  7²:3, A₅, L₃(2), A₆, SL(2,3). Every hit is abelian and flat. No curved example found.
- RELATION: none to the lattice language. Analogy: among graphs with (λ,μ) = (0,2) the curved diameter-3
  examples are the biplanes that are not cube quotients; Γ would be the analogous object for (1,2).
- STATUS: OPEN question with data. No implication for C2: Γ is not vertex-transitive.
  Not in C1's notes or the invention registry (checked by grep: no Cayley search, no biplane analogy).

## 9. Parallelism (second mutation; needs a hypothesis)
- PRIMITIVES: a partition of the 231 triangles into 7 classes of 33 pairwise disjoint triangles, each oriented;
  the seven permutations T₁..T₇ of the 99 vertices (T_i = rotate every class-i triangle).
- DEFINITION: Γ = Schreier graph of ⟨T₁..T₇⟩ on 99 points. Braiding field r_x(a,b) = (c,d) iff T_bT_a x = T_dT_c x.
- MAKES SIMPLE: existence becomes a question about a permutation representation; flat ⟺ the T_i commute.
- LAWS: (i) a class of 33 disjoint triangles meets the Hoffman bound 231·3/(18+3) = 33 in the triangle graph
  (spectrum 18, 7⁵⁴, 0⁴⁴, (−3)¹³²); (ii) if σ preserves the classes and reverses the orientations, then
  ⟨σ, T_i⟩ is a quotient of the Coxeter group of the star K₁,₇ (centre σ, leaves σT_i, leaf pairs free),
  Γ is a degree-99 permutation representation of it in which σ has one fixed point.
- COUNTEREXAMPLES: nothing forces a parallelism to exist in Γ, nor σ to respect one. Both are assumptions.
- RELATION: none to the lattice language; the 7 classes restrict to the 7 lines at x.
- STATUS: OPEN, hypothesis-dependent. Global and not spectral, which is why it is kept.

## 10. Earlier attempts (same day, before this file)
- Roots of N as balanced 7|7 sets of orbits, A_k components as designs with r − λ = 7: COLLAPSED (tight sets, C1 107).
- Star defect Δ_S = 14G_S − G_S² = FFᵀ: COLLAPSED (C1 shell completion).

---

## Why they collapse
1. Every primitive in 1–5 is a function of adjacency inside the ball of radius 2 at x, and that ball is the whole
   graph. They are coordinate changes on (N, B, D). A coordinate change cannot add a constraint.
2. The local data is genuinely free: f(T), the wedge matchings and the D-types are not fixed by the parameters, and
   every linear count cancels them (n₀ + n₁ + n₂ + n₃, E₁+2E₂+3E₃, the hexagon count 140).
3. The only global facts are the eigenvalue multiplicities. Every object that sees them (6, 7) is a rank or a
   positive-semidefinite statement, which is the existing language.
4. Record 8 is the only one that leaves the graph itself. It asks about the class of objects Γ would belong to.

## Where the next mutation has to go
A new primitive has to be global and not spectral. Two candidates, neither started:
- a resolution of the 231 triangles into 7 classes of 33 (global direction labels without flatness);
- a proof or a counterexample for "a (1,2)-Cayley graph on order-3 subgroups has abelian ⟨D⟩".
