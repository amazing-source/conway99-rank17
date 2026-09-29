# Characteristic triangle parity: Type 0 and Type II

**Date:** 29 September 2026.  
**Input:** `amazing-source/conway99-rank17`, commit
`d421b842dc4eafa501c704c8d4822353b36f5aa9`,
`GPT_PROMPT_TYPE0_TYPE2.md`, and the same commit's `RANK17_HANDOFF.md`.

## Status and scope

PROVEN below: the characteristic triangle congruence, its cut normal form,
the missing norm congruence, completeness of the characteristic-parity tests
on the switched frame lattice, and the stated norm-37 exclusions.

COMPUTED: 3 Type-0 and 74 Type-II anchor orbits, with two independent complete
integer enumerations. These are **anchors Phi, not matrices N**. No complete
N-matrix extension was searched, no counted open K-class is newly claimed
excluded, and no bound of one day for the extension problem is claimed.

In particular the earlier suggestion that a genuine characteristic anchor
could occupy any one of 54 two-graph classes was unnecessarily weak. It can
occupy only the class of the all-Opp transversal.

## 1. Notation and hypotheses

Index the 21 cells by 0 <= i < j <= 6. Their labels are Same (+) then Opp (-),
with rows of R equal to e_i+e_j and e_i-e_j, respectively. Thus R has size 42x7.
The signed label matrix is N=Gram(u)-3I, not the negative of this convention.
We use

    NR=0,  N^2=N+12I-RR^T,  diag N=0,
    N_ij in {0,+1,-1},  exactly ten nonzeros per row.

The residual w-vectors are rows of W and satisfy

    ||w_u||^2=10,
    (w_u,w_v)=4N_uv-r_u.r_v for u != v,
    W W^T=12I+4N-RR^T,
    W^T W=28 I_15.

They belong to K, and a_ij=(w_ij^++w_ij^-)/2 and
b_ij=(w_ij^+-w_ij^-)/2 belong to K and have odd norms 5+2t and 5-2t.
Let c in Char(K*) and define Phi_u=(c,w_u)/2. It follows that

    Phi integral, N Phi=4 Phi, R^T Phi=0,
    Phi_ij^+ and Phi_ij^- have opposite parities,
    ||Phi||^2=7||c||^2.

We also use the actual graph-origin half-triangle vectors. Their membership
in K is already established in Proposition 5.1 of the handoff and is
re-derived below.

## 2. PROVEN: the half-triangle congruence

For i<j<k define

    v_ijk=(w_ij^+ + w_ik^+ + w_jk^+)/2.

Then v_ijk belongs to K and has even norm. Consequently

    Phi_ij^+ + Phi_ik^+ + Phi_jk^+ = 0 (mod 2).       (TP)

### Proof of membership and evenness

The three original Same u-vectors have total H-projection e_i+e_j+e_k,
where e_l=h_l/2. Their sum has even pairings with every h_l and lies in U_H.
Subtracting its integral H-projection in the switched lattice leaves v_ijk
in V and orthogonal to the normals, so v_ijk lies in K.

The three Same r-vectors have pairwise product 1. Therefore

    ||v_ijk||^2
      = [30 + 2(4(N_ij,ik + N_ij,jk + N_ik,jk)-3)]/4
      = 6 + 2(N_ij,ik + N_ij,jk + N_ik,jk),

which is even. Characteristicness says (c,v_ijk) is even; this pairing is
exactly the left side of (TP). This proves the result. QED.

There are only 15 independent conditions: take the triangles 0ij with
1 <= i < j <= 6. All 35 triangle conditions follow from those.

## 3. PROVEN: all-Opp gauge, not 54 possible two-graph classes

Set beta_ij=Phi_ij^+ mod 2. Equation (TP) implies

    beta_ij=b_i+b_j  (mod 2),

where b_0=0 and b_i=beta_0i. Change signs of the h-coordinates by
eta_i=(-1)^b_i, permuting the Same/Opp labels in each cell when necessary,
and reorient each u so the first nonzero r-coordinate is again positive.
This is a permissible signed-permutation change of the labelled frame and
preserves the matrix equations, norms and star rules. It gives

    Phi_ij^+ even,  Phi_ij^- odd,  for every cell.     (G)

Indeed changing one endpoint sign swaps the two label parities, precisely
when beta_ij=1. Overall sign of Phi can also be reversed.

The same statement is a single cut-equivalence class of transversals. A
nontrivial two-graph is impossible for a genuine characteristic anchor.
This assertion is stronger than merely saying all characteristic anchors
have mutually switching-equivalent transversals.

### A strictness control

Assign, at vertex i in Z/7, positive half-edges toward i+1,i+2,i+5 and
negative half-edges toward the other three vertices. Select the corresponding
signed D7 root in each cell. The resulting Phi has norm 21, one odd entry
per cell and R^T Phi=0. But its Same coordinates sum to 1 on triangle 012.
It therefore cannot be a characteristic anchor. No claim is made that an N
matrix realizes this deliberately invalid control.

## 4. PROVEN: the norm is 5 modulo 8

In gauge (G) write Phi^+=2t and Phi^-=o with o odd. Summing all seven
coordinate equations R^T Phi=0 gives

    2 sum_cells Phi^+ = 0, hence sum_cells t=0.

It follows that

    ||Phi||^2 = sum o^2 + 4 sum t^2
              = 21 + 4 sum t (mod 8)
              = 5 (mod 8).                           (NC)

Thus 25 and 33 are impossible. Since every cell contributes at least one,
the only possible norms <=37 are 21,29,37.

An independent lattice-side check is available. Every characteristic c of
K* lifts to a characteristic chi of the unimodular V. If a=(chi,y1) and
b=(chi,y2), then a,b are odd, chi^2=17 (mod 8), and

    7c^2=7chi^2-a^2-b^2=5 (mod 8).

The elementary proof (NC) needs no modular-form theorem or shell enumeration.

**Scope caution.** The minimum m* as literally defined in the prompt ranges
over arbitrary integral eigenvectors with one odd coordinate per cell. It
omits (TP). We have not proved (TP) for every vector in that larger set.
Results here apply to genuine characteristic vectors, or to the enlarged
V-free axiom system with (TP) explicitly included. Enlarging a necessary-
condition system in that way is sound for every physical K-class.

## 5. PROVEN/COMPUTED: all Type-0 anchors reduce to three tournaments

At norm 21, (G) forces Phi^+=0 and Phi^-=+1 or -1 in every cell. Give ij the
arrow i->j when Phi_ij^-=+1. Coordinate balance is

    outdeg(i)-indeg(i)=0.

Thus the tournament is regular, with outdegree three at every vertex.
Conversely every such tournament gives an anchor satisfying integrality,
parity, balance and norm 21. It does not automatically have an N-extension.

Residual symmetries of gauge (G) are vertex permutations and global sign
reversal, the latter reversing all arrows. Exact enumeration gives

    labelled regular tournaments: 2640,
    orbits under S7 and reversal: 3,
    orbit sizes: 720, 1680, 240.

Cells are lexicographic, bit 1 means i->j for i<j. Representatives:

    17976, 22066, 85298.

The Python enumerator fills edges recursively and prunes only impossible
outdegrees. The independent C++ checker scans all 2^21 tournaments and
finds orbits by BFS under generators. They agree on the representative
masks and all orbit sizes, not just the orbit count.

## 6. PROVEN/COMPUTED: all Type-II anchors reduce to 74 patterns

At norm 29 the excess over 21 is eight. There are two possible magnitude
profiles:

1. one Opp entry has magnitude three, all other Opp entries magnitude one,
   all Same entries zero (the already excluded Type I);
2. two Same entries have magnitude two, all Opp entries magnitude one.

For (1), the magnitude-three row requires 12 but its ten neighbors have
magnitude at most one. This is impossible.

For (2), sum Phi^+=0 forces the two Same defects to have opposite signs:
+2 and -2. There are two incidence types.

### Intersecting defects

Normalize the +2 defect to cell 01 and the -2 defect to cell 02. The Opp
entries are a tournament with outdegree vector

    (3,2,4,3,3,3,3).

The remaining group has order 48: S4 on vertices 3,4,5,6 and the operation
which swaps 1,2 while reversing every arrow (global Phi sign).

There are 1720 labelled tournaments and 42 orbits. Orbit-size histogram:

    30 of size 48, 11 of size 24, 1 of size 16.

### Disjoint defects

Normalize the +2 defect to 01 and the -2 defect to 23. The outdegree vector is

    (2,2,4,4,3,3,3).

The remaining group has order 48: S2 on each defect pair, S3 on vertices
4,5,6, and interchange of the two defect pairs together with reversal.

There are 1112 labelled tournaments and 32 orbits. Orbit-size histogram:

    15 of size 48, 16 of size 24, 1 of size 8.

Thus Type II has 42+32=74 characteristic anchor orbits. The accompanying
`anchors_42.json` supplies all 74 literal vectors, plus the three Type-0
vectors, in the stated canonical R convention.

### Completeness argument

Every genuine Type-II anchor can first be put in gauge (G). Permuting the
seven coordinates puts its ordered positive/negative defect edges in one
of the displayed forms. The exact balance equations force the displayed
score vector. Every remaining tournament is generated by the exhaustive
score enumeration, and quotienting by precisely the normal-form stabilizer
retains a representative. This proves completeness at the anchor level.
It neither bounds nor solves the remaining N-matrix completion problem.

## 7. PROVEN: magnitude profiles through norm 37

In (G) let d=sum (o_e^2-1)/8. Then

    ||Phi||^2=21+8d+4 sum t_e^2,  sum t_e=0.

Below or at 37, o_e is +/-1 or +/-3. The only nonzero even entries are +/-2:
a single +/-4 cannot be balanced within the remaining norm budget.

The complete surviving magnitude possibilities are:

| norm | Opp entries | Same entries |
|---|---|---|
| 21 | 21 unit entries | zero |
| 29 | 21 unit entries | one +2, one -2 |
| 37 | 21 unit entries | two +2, two -2 |
| 37 | one magnitude 3 and 20 units | one +2, one -2 |
| 37 | two magnitude 3 and 19 units | zero |

The norm-29 lone magnitude-three profile is impossible by row capacity.
For the norm-37 one-three/two-two profile, the three's row must be adjacent
to both two's and eight units, with all ten eigen-equation contributions
aligned. This is forced because 12=2+2+8 is its maximum possible row sum.
The three's cell may or may not be one of the two Same-defect cells; both
cell-incidence variants are retained in the table.

## 8. PROVEN: with two magnitude-threes at norm 37 they form a directed path

Let i,j be the two magnitude-three labels. Switch signs of label vectors
for the proof so all 21 nonzero Phi entries are positive. Write r_i,r_j
for the correspondingly oriented D7 roots. All these labels are Opp in the
previous canonical gauge, so distinct labels have root products in {-1,0,1}.

Each three's eigen-equation is 12. With ten nonzero N-entries, equality
forces N_ij=1 and nine unit neighbors, all with coefficient +1. It has no
neighbors in the 21 zero labels. Let A,B be the two nine-element unit
neighbor sets, let C=A intersect B, and let Z be the unit labels in neither.
There are 19 unit labels, so |Z|=1+|C|.

The quadratic matrix equation gives

    |C|=(N^2)_ij=1-r_i.r_j.

NR=0 on the two defect rows gives

    sum_A r=-r_j,  sum_B r=-r_i.

The anchor balance gives sum_units r=-3r_i-3r_j. Subtracting yields

    sum_Z r - sum_C r = -2(r_i+r_j).                (8.1)

If r_i.r_j=1, then |C|=0 and |Z|=1. The left side has norm two,
whereas the right side has norm 24. Impossible.

If r_i.r_j=0, then |C|=1 and |Z|=2. The left side is a signed sum of three
roots from distinct cells. Their pairwise products are at most one, so
its norm is at most 6+6=12. The right side has norm 16. Impossible.

Thus r_i.r_j=-1: the two oriented Opp roots must be successive edges of
a directed path on three vertices. A disjoint pair, or a same-direction
fork at a vertex, is completely excluded. The directed-path case is not
claimed excluded.

## 9. PROVEN: what is the full characteristic-parity content?

Let L be the rank-15 lattice obtained by switching the original frame
lattice N_frame, before any odd-index 7-gluing. An element of L is the
projection of sum m_u u_u with R^T m even, and equals (1/2)sum m_u w_u.

Over F2, R^T is the incidence matrix of the complete graph with each edge
doubled (Same and Opp). Its kernel is generated by the 21 double-edge
cycles and 15 based triangles. L is therefore generated integrally by

    w_u,
    (w_ij^+ + w_ij^-)/2,
    (w_0i^+ + w_0j^+ + w_ij^+)/2.

To see integral generation, subtract representatives of the mod-2 cycle
basis from m; the remainder has even coefficients and is generated by the
2e_u vectors. Their residual images are precisely w_u.

For an integral functional c on L, characteristicness on these generators
is exactly:

1. Phi integral (w has even norm);
2. one odd Phi entry per cell (the double-edge generator has odd norm);
3. (TP) on the 15 based triangles (the half-triangle has even norm).

Because norm parity is additive on an integral lattice, these conditions
are sufficient as well as necessary for c to be characteristic on L.
If L is contained in K with odd index and c belongs to K*, the condition
extends to K: for odd d with dv in L,

    d(c,v)=(c,dv)=||dv||^2=d^2||v||^2 (mod 2).

This is a completeness statement about characteristic parity, not about
all lift/Gram conditions on a putative frame. It does not discard M1/M2,
R1-R3, the eigen-equations, or candidate-lattice membership.

The missing relation can be written modulo four on Wc:

    (c,w_ij^+ + w_ik^+ + w_jk^+) = 0 (mod 4).

There is no assertion here that arbitrary mod-4 conditions on Phi itself
are automatic. For multiple characteristic vectors one common N,R must be
used; they cannot all be normalized independently without transporting
the same unknown matrix consistently.

## 10. Deliverable and limitation

The package closes the nontrivial characteristic-transversal switching
classes, the m=25 and m=33 possibilities, and the stated norm-37 defect
geometries. It gives a complete explicit list of 77 anchors for the two
central cases, independently checked.

It does **not** prove no Type 0 or no Type II. It does not justify changing
the ledger by 76 or 378, or by any other number. The scalar prompt's R5
status is not treated here as an independently checked exclusion.

The next matrix search may start directly with these 77 R,Phi inputs and
must retain every necessary graph condition. No claim is made that doing
so finishes in under one day, or that the star-matching count is reduced
by the ratio of anchor counts.

## 11. COMPUTED: an explicitly symmetry-restricted construction attempt

As a separate experiment I tried an N invariant under the order-21 affine
group of the Paley regular tournament. This is an **additional symmetry
assumption**, not something implied by the anchor.

Let the vertices be F7 and take differences 1,2,4. Both the Same and Opp
label sets are regular orbits of the group x -> a*x+b, a in {1,2,4}. Choose
r^+_(a,b)=e_b+e_(a+b), r^-_(a,b)=e_b-e_(a+b). The Type-0 anchor is zero on
the Same orbit and one on the Opp orbit.

Every invariant symmetric N is described by two inverse-symmetric functions
A,B on this group, zero at the identity, and a function C for the cross
block. The group has 10 nonidentity inverse pairs. NR=0 implies the Same
row sum is zero. N Phi=4 Phi implies the Opp row sum is four and the cross
row/column sums are zero. Ten nonzeros per row leave only these two modes:

| mode | Same intra-orbit inverse pairs | Opp inverse pairs | cross entries |
|---|---|---|---|
| 0 | 1 positive, 1 negative | 2 positive | 3 positive, 3 negative |
| 1 | 2 positive, 2 negative | 3 positive, 1 negative | 1 positive, 1 negative |

The optional exact script enumerates both modes, uses NR=0 and the R1
per-star row capacities, then checks N^2=N+12I-RR^T. It found:

| mode | cross patterns | NR-feasible triples | row-cap feasible | N^2 solutions |
|---|---:|---:|---:|---:|
| 0 | 1,085,280 | 4,208 | 481 | 0 |
| 1 | 420 | 219 | 69 | 0 |

All arithmetic is integer. The script and log are included; unlike the
anchor census this optional experiment has not received a second independent
implementation. It does not exclude a general N over the Paley anchor:
the matrix need not inherit the tournament's automorphisms. The argument
has not been promoted to a Type-0 closure or a new ledger exclusion.
