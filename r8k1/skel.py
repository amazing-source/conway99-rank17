"""R8 k=1 enumeration, step A-B: W-skeletons.
Five W labels = five distinct D6 labels (cells inside B0 = lines 0..5, c = 6), each with the root rho_k of its h=+2
representative.  Constraint: net_d = sum_k (rho_k)_d in {0,+1,-1} for every line d (balance at d^+-).
Symmetry used: vertex permutations of the cell multigraph (full S6 via canonical shapes) and sign flips of lines
(normalised: first incident edge endpoint of each used line has sign +1)."""
import itertools
LINES = 6
def shapes():
    """multigraphs with 5 edges on {0..5}, multiplicity <= 2, up to isomorphism (brute-force canonical form)."""
    pairs = list(itertools.combinations(range(LINES), 2))
    seen = set(); out = []
    for combo in itertools.combinations_with_replacement(range(len(pairs)), 5):
        cnt = {}
        for p in combo: cnt[p] = cnt.get(p, 0) + 1
        if max(cnt.values()) > 2: continue
        edges = [pairs[p] for p in combo]
        best = None
        for perm in itertools.permutations(range(LINES)):
            e2 = tuple(sorted(tuple(sorted((perm[a], perm[b]))) for a, b in edges))
            if best is None or e2 < best: best = e2
        if best in seen: continue
        seen.add(best); out.append(best)
    return out
def signed_skeletons(shape):
    """assign (alpha_i, alpha_j) to each edge; parallel edges must have different products; normalise line flips."""
    E = list(shape)
    res = []
    for signs in itertools.product([(1, 1), (1, -1), (-1, 1), (-1, -1)], repeat=5):
        ok = True
        # parallel edges distinct labels: product alpha_i*alpha_j must differ
        for x in range(5):
            for y in range(x+1, 5):
                if E[x] == E[y] and signs[x][0]*signs[x][1] == signs[y][0]*signs[y][1]: ok = False
        if not ok: continue
        # normalisation: for each line, first incident edge (lowest index) has alpha = +1 at that line
        first = {}
        for x, (i, j) in enumerate(E):
            for d, al in ((i, signs[x][0]), (j, signs[x][1])):
                if d not in first: first[d] = al
        if any(v != 1 for v in first.values()): continue
        net = [0]*LINES
        for x, (i, j) in enumerate(E):
            net[i] += signs[x][0]; net[j] += signs[x][1]
        if any(abs(n) > 1 for n in net): continue
        rhos = []
        for x, (i, j) in enumerate(E):
            r = [0]*7; r[i] = signs[x][0]; r[j] = signs[x][1]; rhos.append(tuple(r))
        res.append(tuple(rhos))
    return res
if __name__ == '__main__':
    S = shapes(); print('shapes', len(S))
    tot = 0
    for sh in S:
        sk = signed_skeletons(sh); tot += len(sk)
    print('signed skeletons (flip-normalised, shape automorphisms not removed):', tot)
