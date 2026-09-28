# R6 data checks (28/09 evening, C1). Exact short vectors; no solver. Scripts are in this folder.

Inputs:
- the 754 open classes: conway99-c1 data/rank17/open_after_sol.json;
- the K Grams: conway99-c1 data/rank17/K_*.json.

## n6_config.py: n=6 witnesses in the 555 open R6 classes
- Number of n=6 vectors (up to sign) per class:
  | vectors | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 14 |
  |---|---|---|---|---|---|---|---|---|---|
  | classes | 220 | 129 | 92 | 48 | 28 | 22 | 12 | 2 | 2 |
- Pairwise |7(k,k')| between n=6 witnesses lies in {0, 1, 2} only.
- There is NO linked pair (|7(k,k')| = 3, which would make k ∓ k' another n=6 witness), in all 555 classes.

## n6_joint.py / n6_triples.py: short triples
A short triple is k, k', k ∓ k' with all three norms n in 5..8. By linearity it forces three graph eigenvector patterns with
a1 − a2 = a3.
- 395 of the 754 open classes contain a short triple.
- By type:
  | (n, n', n'') | classes |
  |---|---|
  | (6, 8, 8) | 302 |
  | (6, 6, 8) | 157 |
  | (5, 5, 6) | 88 |

  A class can contain several types.

## Why this matters
- The single-witness cells (J1a types (iii) and (v)) pass every local test tried. A triple theorem uses a strictly stronger
  hypothesis: three patterns on ONE graph, coupled linearly, with cut(a3) = cut(a1) + cut(a2).
- "No graph realises a (6,8,8) triple" would close 302 classes directly, without the full R6 or R8.
