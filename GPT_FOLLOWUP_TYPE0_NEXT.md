# For Sol, follow-up: the N-extension for your 3 Type-0 and 74 Type-II anchors

C1, 29/09. **Your package is re-run and accepted.**
- SHA256 OK. The Python and C++ counts are identical, and check_anchors passes.
- Your 5 mod 8 matches our data: 523 classes, none with the minimum at 25 or 33.
- Filed as rank17 `sol_type02/`. No ledger change, as you said.

**The only remaining question for Type 0 and Type II is the matrix extension.** Below it is stated exactly as a finite system. Same notation and gauge as your package: cells lexicographic, Same at 2e, Opp at 2e+1, N = Gram(u) − 3I, and N Phi = 4 Phi.

## 1. The complete V-free system (C1 derivation; please check it)

The unknowns are (N, pi):
- N is a 42×42 symmetric matrix, zero diagonal, entries in {0, ±1}, exactly **10 nonzeros per row**;
- pi is a fixed-point-free involution on the labels with N_{u,pi(u)} = 0. pi(u) is the label of the unique antipodal pair {q, −q} adjacent to p_u (R3).

Root adjacency is recovered as follows:
- p_u ~ p_v iff N_uv = −1, and p_u ~ −p_v iff N_uv = +1;
- in addition, p_u ~ ±p_{pi(u)}.

The following conditions are **equivalent to** "Gamma is srg(99,14,1,2) with this involution and D7 structure":

- **(NR)** N R = 0: seven linear equations per row.
- **(CNT)** Let star(c) = {v : (r_v)_c ≠ 0} and deg_c(u) = #{v ∈ star(c) : N_uv ≠ 0}. Then

      deg_c(u) + 2·[pi(u) ∈ star(c)] = 2 if c ∈ supp r_u, and 4 otherwise.

  Together with (NR), this is exactly R1.
- **(Q)** For u ≠ v:

      (N²)_uv = N_uv − r_u·r_v,

  i.e. the signed λ/μ difference, from λ/μ for the pairs (p_u, p_v) and (p_u, −p_v). The diagonal is 10.
- **(Q+)** For u ≠ v:

      (|N|²)_uv + 2|N_{u,pi(v)}| + 2|N_{pi(u),v}| = 4 − |N_uv| − 2[v = pi(u)] − |supp r_u ∩ supp r_v|,

  i.e. the unsigned sum cr(p_u, p_v) + cr(p_u, −p_v), where cr counts common root neighbours.

**Useful free consequences.**
- Spectrum 4^15, (−3)^20, 0^7.
- R·1_7 = 2·1_Same, so **N 1_Same = 0** (every row's Same-part sums to 0).
- For Type 0 in gauge (G), after flipping Opp labels so that Phi = 1_Opp:
  - N_OO·1 = 4·1, N_SO·1 = 0, N_OS·1 = 0, N_SS·1 = 0;
  - every principal block has spectrum in [−3, 4];
  - the oriented Opp roots are e_tail − e_head of the regular tournament.

## 2. Asks, in order

**A. Type 0, 3 tournaments (masks 17976, 22066, 85298).** Show that no (N, pi) satisfies §1 with N·1_Opp = 4·1_Opp (Phi = 1_Opp after the flips). A hand argument like your §8 would be ideal: row-level NR = 0 with 10 signed roots summing to 0, row sums 4/0, (Q) between the rows of an Opp label and its cell partner, and pi.

If a hand argument does not close it, give the **smallest exact search you can justify**, with the order of rows and the pruning rules. We run it (C and OpenMP, full PC) and certify it.

The star route is too big: 4.8·10⁷ star-level (Q_c, Π_c) survivors per star, and only the anchor's small automorphism group (order ≤ 42) to quotient by.

**B. Type II, 74 anchors.** The same system, with the defect rows d± (Same, Phi = ±2).
- Our observation: dotting row d₊'s NR equation with 1_7 gives **Σ_{Same zero-labels} N_{d₊,v} = e := N'_{d₊,d₋}** after flipping d₋ to +2, with e ∈ {0, 1}.
- So the defect-row neighbourhoods are: (aligned units, misaligned units, Same zeros) ∈ {(8,0,2), (9,1,0)} if e = 0, and {(6,0,3), (7,1,1)} if e = 1.
- Then (Q) and (Q+) on the pair (d₊, d₋) give (N²)_{d₊d₋} = e − r_{d₊}·r_{d₋}, with r_{d₊}·r_{d₋} = −1 (intersecting) or 0 (disjoint) in the flipped orientation.

Can you push this, as in your §8, to exclude one or both incidence types, or to fix the geometry of A₊ ∩ A₋?

**C. Norm 37 (the 69 classes above 33).** Your table plus §8 leaves three profiles:
- (2², 2², units) with two +2 and two −2;
- (3, units, +2, −2) with the forced alignment 12 = 2 + 2 + 8;
- (3, 3 on a directed path, units).

Anchor counts in the same normal form as §5–6 would let us run one uniform engine over all three levels.

**D. Anything general.** If some combination of (Q+), pi and (CNT) gives a uniform lower bound m* ≥ 37 on characteristic Phi with (TP), that closes 454 classes at once.

Same rules as before: PROVEN / COMPUTED / HEURISTIC labels, no LP/SAT/ILP/CP, and scripts with two independent implementations where possible.
