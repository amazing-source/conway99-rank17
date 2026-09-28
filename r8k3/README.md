# R8, k = 3: witness sign classes (C2, 28/09)

`r8k3_classes.json`: every orbit-level pattern h of type n = 8, k = 3, up to symmetry. The pattern is:
- ±1 on the 24 cut orbits of the 3|4 cut I = {1,2,3} | J = {4,5,6,7};
- ±2 on exactly two non-cut orbits;
- 0 elsewhere;
- block balance at the 14 inner vertices.

**Symmetry group** (order 36,864): S₃(I) × S₄(J) × per-block flips × global sign.

**Counts.** 612 placements of the ±2 pair give **8 placement types**; in all there are **209 classes**. Per type:

| type | ±2 placement | orbit | classes |
|---|---|---|---|
| 0 | (1,2,0), (1,2,1) | 12 | 12 |
| 1 | (1,2,0), (1,3,0) | 24 | 15 |
| 2 | (1,2,0), (4,5,0) | 288 | 78 |
| 3 | (1,2,0), (2,3,1) | 24 | 18 |
| 4 | (4,5,0), (4,5,1) | 24 | 11 |
| 5 | (4,5,0), (4,6,0) | 96 | 15 |
| 6 | (4,5,0), (6,7,0) | 48 | 25 |
| 7 | (4,5,0), (5,7,1) | 96 | 35 |

**Conventions.** Those of orbit_sat:
- orbit U = (i, j, e) with i < j;
- representative {(i,0), (j,e)}, i.e. the D₇ root e_i + (−1)^e e_j;
- h_U is the value on the representative.

**Scope.** This is only the witness-pattern layer: exact enumeration, no solver, no graph model. The script `r8k3_classes.py` regenerates it:
- `--placements` lists the 8 placement types;
- `--place t` enumerates type t.
