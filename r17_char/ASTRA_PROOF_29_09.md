# Rank-17 characteristic-root closures

**Date:** 29 September 2026.

**Result:** ten explicitly specified normal-vector orbits and nine explicitly named catalogue pairs cannot support the full switched C2 frame. The nine catalogue pairs occur in the pinned `open_after_sol.json` at C1 commit `888e731149a9935e189b5af8774d51fc9df91447`. This does not assert that they remain unclosed by every other route today. No whole rank-24 envelope count is increased here, and the full involution case is not closed.

All norms in this note are squared Euclidean lengths.

## 1. The exact premise and normalization

Use the established rank-17 switch. An integral unimodular lattice V contains orthogonal vectors y1,y2 of norm 7. In

    K = V intersect span(y1,y2)^perp

there are 42 lattice vectors x, the a_cd,b_cd, with cell norm pairs (5,5), (3,7), or (7,3), and

    sum_x x x^T = 14 P_K.                         (1)

Every one of these 42 vectors has odd norm. In particular (1) implies the following operator identity on the entire 17-dimensional ambient space:

    sum_x x x^T + 2 y1 y1^T + 2 y2 y2^T = 14 I. (2)

This identity uses the **a,b** vectors, not the **w** vectors. The latter have tight constant 28. No rank-22 characteristic bound is transferred to K here.

For completeness, if w+,w- are the two norm-10 residual vectors of a cell, put a=(w++w-)/2 and b=(w+-w-)/2. Their original u-inner product t is 0 or +-1 and their D7 columns are orthogonal. Thus w+ dot w-=4t, so a^2=5+2t and b^2=5-2t. The contribution aa^T+bb^T is half the two w outer products. This gives the odd norms and (1) from the established 28-tight w frame. Lattice membership of a,b is part of the switch (subtract the integral unit H-projection from u++-u- in the switched lattice).

## 2. PROVEN: characteristic energy

A characteristic covector chi of V satisfies

    (chi,v) == v^2 mod 2 for every v in V.

V is unimodular, so such a covector belongs to V. For each of the 42 x, (chi,x) is an odd integer, hence its square is at least 1. Apply (2) to chi:

    14 chi^2 = sum_x (chi,x)^2 + 2[(chi,y1)^2+(chi,y2)^2].

Consequently

    (chi,y1)^2+(chi,y2)^2 <= 7 chi^2 - 21.       (3)

For chi^2=9 the right-hand side is 42. The two heights are odd. A norm-seven normal of height 7 therefore contradicts (3), irrespective of the partner: 7^2+1^2=50>42. A pair of heights 5,5 likewise contradicts it.

For all nine explicit pair certificates in this package the available x-energy is exactly

    14*9 - 2*(49+1) = 26,

or equivalently 14*9 - 2*(25+25)=26. The 42 odd pairings require at least 42. Every certificate has strict gap 16.

## 3. PROVEN: characteristic-root normal-orbit theorem

Suppose y in V has norm 7, r in V has norm 2, (r,y)=0, and chi=y+r is characteristic. Then y cannot be either normal in a switched C2 frame.

**Proof.** Every residual x is perpendicular to y, so

    (r,x) = (chi,x)

is odd. Thus sum_x(r,x)^2 >=42. On the other hand (1) gives

    sum_x(r,x)^2 =14 ||P_K r||^2 <=14 r^2=28,

a contradiction. The partner need never be chosen or enumerated. QED.

Equivalently, chi^2=9 and (chi,y)=7, so section 2 applies. Isometries and sign changes transport both the root and characteristic witness. Thus each explicit y below closes its entire Aut(V)-orbit.

## 4. Ten explicit normal-orbit closures

The listed chi is characteristic, y^2=7, chi^2=9, and r=chi-y is a norm-two vector orthogonal to y. Coordinates are Euclidean; commas between components denote orthogonal sums. Complete vectors are in certificates.json.

### 4.1 E8 + Z9

    chi=(0;1^9), y=(0;1^7,0^2).

E8 is even and the Z9 coordinates of chi are odd. The root r is the sum of the last two unit vectors.

### 4.2 D12+ + Z5: two orbits

Use D12+ = D12 + Z(1/2)^12 and

    chi=(2e1;1^5).

It pairs evenly with D12 root generators, oddly with the norm-three spinor glue, and oddly with the unit generators. The two y are

    (e1+e2;1^5),
    (2e1;1,1,1,0,0).

Their root witnesses are (e1-e2;0) and (0;0,0,0,1,1), respectively.

### 4.3 (E7+E7)+ + Z3: two orbits

Realize E7 in E8 perpendicular to (0^6,1,1), and put

    omega=(1/2,1/2,1/2,1/2,1/2,1/2,0,0), omega^2=3/2.

Use the unimodular glue (omega,omega). Then

    chi=(2omega,0;1,1,1)

pairs evenly with roots and oddly with the norm-three glue and units. Take

    y=(e1+e2+e3+e4,0;1,1,1),
    y=(2omega,0;1,0,0).

The roots chi-y are (e5+e6,0;0) and (0,0;0,1,1).

### 4.4 A15+ + Z2: two orbits

In the A15 sum-zero hyperplane of R16 use the glue (1/4)^12,(-3/4)^4. Put

    c=(7/4,7/4,(-1/4)^14), c^2=7,
    chi=(c;1,1).

Its pairing with every root is even; its pairing with the norm-three glue is 1. Take

    y=(c-(e1-e3);1,1),
    y=(c;0,0).

The roots are (e1-e3;0,0) and (0;1,1).

### 4.5 (D8+D8)+ + Z: one orbit

Use the two glues (s,v),(v,s), where s=(1/2)^8 and v=e1. Set

    chi=(2e1,2e1;1),
    y=(e1+e2,2e1;1).

The glue pairings of chi are both 3, hence odd; all root pairings are even. The root chi-y is (e1-e2,0;0).

### 4.6 A11E6+: two orbits

Use the model in section 5. Write L_i for the vector in R12 with entry -11/6 in position i and 1/6 elsewhere. One normal is

    y0=(L_11;
        1/2,1/2,-3/2,-1/2,1/2,-1/6,-1/6,1/6).

Take

    r0=(0;
        1/2,1/2,1/2,-1/2,1/2,1/2,1/2,-1/2).

Here r0 is an E6 root. The resulting chi is

    chi0=(L_11;1,1,-1,-1,1,1/3,1/3,-1/3).

The norm and characteristic checks are proved by the generator test in section 5.

For the other orbit take

    chi=(-3/2,-3/2,-3/2,(1/2)^9;0^8),
    r=e4-e1 in A11,
    y=chi-r.

Root pairings are even, and the glue pairing of chi is -1, hence odd. Its norm is 9 and (chi,r)=2, as required.

## 5. Hand proof for the A11E6p:8, :9, :10 bottlenecks

### The model

The root lattice is A11 perpendicular E6, in R12 perpendicular R8. A11 is the sum-zero integer lattice. E6 is the subsystem of E8 perpendicular to

    (0^6,1,1), (0^5,1,-1,0).

The unimodular lattice is generated by these roots and the glue

    g=((1/6)^10,(-5/6)^2;
       1,0,0,0,0,1/3,1/3,-1/3),    g^2=3.

The root-lattice determinant is 12*3=36 and the glue has order 6, so the integral overlattice has determinant 1. The exact model audit also verifies this directly by HNF.

The E6 roots are precisely:
- +-e_i+-e_j for 1<=i<j<=5;
- (epsilon_1,...,epsilon_5,t,t,-t)/2, where all signs are +-1 and t*product(epsilon_i)=-1.

### A uniform characteristic-vector check

For any i, any signs s1,...,s5 with product +1, the vector

    chi=(L_i;s1,s2,s3,s4,s5,1/3,1/3,-1/3)     (4)

is characteristic of norm 9.

A11 root pairings are 0 or +-2. Pairings with integer E6 roots are even. On a half-root, twice the pairing is

    s1 epsilon1+...+s5 epsilon5+t.

These six signs have product -1, so their sum is 0 or +-4, making the pairing even. Finally, pairing with the glue is s1 or s1+2, both odd. Since the roots have even norm and g has odd norm, (4) is characteristic on the entire integral lattice. Its norm is 11/3+16/3=9.

### The named partners

The relevant normal in each pair is:

| Catalogue pair | Normal y2 | Root r in the E6 block |
|---|---|---|
| A11E6p:8 | (L_11; -1,0,1,0,1,1/3,1/3,-1/3) | (0,1,0,-1,0,0,0,0) |
| A11E6p:9 | (L_3; 0,0,-1,-1,1,1/3,1/3,-1/3) | (1,1,0,0,0,0,0,0) |
| A11E6p:10 | (L_2; -1,0,0,-1,1,1/3,1/3,-1/3) | (0,1,1,0,0,0,0,0) |

In every row r^2=2 and (y2,r)=0. Each chi=y2+r has the form (4). Hence the entire orbit of that normal is impossible, regardless of the other normal. In particular all three named pair classes are closed.

The common catalogue y1 has A11 coordinates (1,1,7,1,1,1,-5,-5,1,1,-5,1)/6 and E6 coordinates (1/2,1/2,-3/2,-1/2,1/2,-1/6,-1/6,1/6). Its characteristic heights for the three witnesses are respectively 1,1,-1; y2's height is 7. All pair certificates therefore have Y_chi=50.

This is a genuine infinite-under-isometries normal-orbit statement, not a check of a truncated partner list.

## 6. Further named catalogue closures

The same 42-versus-26 contradiction closes:

| Pair | Characteristic heights |
|---|---|
| E8:5 | (7,-1) |
| E8:9 | (7,-1) |
| E8:11 | (7,-1) |
| E8:14 | (5,-5) |
| E8:19 | (5,-5) |
| A11E6p:102 | (-5,5) |

For all E8 cases chi=(0^8;1^9). For A11E6p:102 it has zero E6 component and A11 entries -3/2 at positions 3,6,10 and +1/2 elsewhere. Roots pair evenly and the glue pairs to -1, so this too is characteristic of norm 9. Full pair coordinates are retained in certificates.json.

Combined with section 5, these are nine named catalogue entries, all present in the pinned open list. This is not a claim that exactly nine current classes are newly excluded: no full scan of all pair files was performed this round, and another route may already imply some exclusions.

## 7. Limits of the result

This argument has not closed V=Z17, all of A11E6+, or the involution case. In particular, explicit pair A11E6p:100 attains maximum 34 over a constructed 1088-vector norm-nine characteristic shell and survives this witness family. That is a limited witness-family non-rejection, not a full-frame candidate and not a proof the overall route cannot close.

The proof here derives a rank-15 characteristic-covector lower bound 3 from the odd a,b frame. It does not import the rank-22 bound 46/7 and does not confuse the 14-tight a,b frame with the 28-tight w frame.

Nor does death after assigning nine **complete rows** prove that a nine-vertex principal submatrix is intrinsically impossible. For selected i,j,

    (T^2)ij=sum_{ell=1}^{42} Tiell Tellj

uses all outside columns. A depth-nine rejection can depend on assignments and completion bounds across the entire 42-label universe. A universal nine-row obstruction still needs proof.

## 8. Exact verification and source boundary

`verify.py` uses only Python fractions and integers. It verifies generator integrality, rank and span, characteristic congruences, every witness norm, every named pair's Gram 7 I2, the root-completion identities, and strict energy gaps. It also rejects 19 deliberately corrupted certificates.

`audit_models.py` independently uses exact SymPy HNF, determinants and rational matrix operations to verify that all six lattice models are positive definite, integral, odd and unimodular, and that every supplied normal, partner, characteristic vector and root actually belongs to its model. It verifies 57 vector memberships. It is not an optimization solver or a frame search.

These checks do not independently reprove the graph-to-switch implication or the completeness of the team's pair enumeration; neither is needed to prove nonexistence for the specifically displayed subcases conditional on the established switch.

Sources read:
- amazing-source/conway99-rank17, RANK17_HANDOFF.md at 7ac7971d3a2ea333b0b997762a77a0a7bb8db1e8, section 4.
- amazing-source/conway99-c1, code/rank17_cores.py at 888e731149a9935e189b5af8774d51fc9df91447.
- At that same C1 commit: data/rank17/ypairs_E8.json, ypairs_A11E6p.json, open_after_sol.json; notes/154_2x5_ENGINE_AND_ROUTE_CHANGE.md.
- The standard characteristic-vector definition agrees with Elkies, Lattices and codes with long shadows, arXiv:math/9906086. No modular-form theorem is required for the new contradiction.

## Verdict

CLOSED for the specified normal orbits and named rank-17 pair classes. The full C2 involution case remains unclosed.
