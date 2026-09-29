# Type 0 / Type II: characteristic parity and exact anchor normal forms

Read `PROOF.md` first. This responds to `GPT_PROMPT_TYPE0_TYPE2.md` at rank17
commit `d421b842dc4eafa501c704c8d4822353b36f5aa9`.

## Results

* The Same entries of a genuine characteristic anchor satisfy every triangle
  parity equation. They form a cut, so coordinate switching makes every Same
  entry even and every Opp entry odd. Only one two-graph class is possible.
* With R^T Phi=0, this forces ||Phi||^2 = 5 mod 8. Norms 25 and 33 are impossible.
* There are exactly 3 Type-0 anchor orbits and 74 Type-II anchor orbits under
  the specified signed-coordinate symmetries. Literal inputs are in
  `anchors_42.json` (canonical R and each 42-coordinate Phi).
* The complete magnitude-pattern table through norm 37 and a proof excluding
  disjoint/fork geometry for the two-three norm-37 case are in `PROOF.md`.

No N-matrix extension has been searched or claimed impossible. No new open
K-class count is claimed. These are not full frame candidates.

## Reproduce

```sh
python enumerate_anchors.py --output anchors.json
g++ -O2 -std=c++17 independent_verify.cpp -o independent_verify
./independent_verify > independent_counts.json
python check_anchors.py
```

Python uses score-constrained recursive generation and all stabilizer
permutations. C++ independently scans all 2^21 tournaments and uses group-
generator orbit closure. Both agree on every representative and orbit size.
The check also verifies every coordinate sign switch of all 77 exported
anchors, and a deliberately invalid balanced Type-0 parity control.

No external Python dependencies, no optimizer/solver, no network access,
no writes outside this directory. The C++ binary is not included.

## Engine input

Vertices are 0..6, cells lexicographic, Same at index 2e and Opp at 2e+1.
Use N=Gram(u)-3I, so N Phi=4 Phi and N^2=N+12I-RR^T.
Each `anchors_42.json` record includes its type, norm, and Phi.
Transport N,R,Phi together under sign switches. In a multi-characteristic
search, fixing one anchor's gauge does not independently fix all others.

## Sources

* Team prompt and handoff at the pinned commit above.
* The optional independent check using characteristic norms in unimodular V
  uses the standard norm-congruence fact also stated by N. D. Elkies,
  "Lattices and codes with long shadows", arXiv:math/9906086.

The primary proof of the 5 mod 8 result in this package is elementary and
does not depend on that external theorem.

## Optional restricted construction experiment

`explore_paley_invariant.py` searches an EXTRA affine-symmetry ansatz for N
over the Paley Type-0 tournament. It uses NumPy integer arithmetic, not an
optimization solver. It found no such invariant N. The ansatz is not a
complete reduction of Type 0 and must not be assumed for other N extensions.
Its record is `paley_invariant_experiment.json`; it has only one implementation.
The main anchor enumeration and proofs do not depend on this experiment.
