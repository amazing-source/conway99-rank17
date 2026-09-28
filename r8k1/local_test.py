import numpy as np, itertools
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')[:800]
def star_adj(N):
    # label-level star graph: u ~ w iff N_SS[u,w] != 0 (pi or Q partner; D-edges cancel -> treat via adjacency lookup)
    return {u: [w for w in range(12) if w != u and N[u, w] != 0] for u in range(12)}
from collections import Counter
stats = Counter()
for ln in lines:
    rho, adj, h = load(ln); M, N, R = MU(rho, adj)
    reps = [15 + u for u in range(12)]
    # star adjacency incl. D edges (use vertex adjacency: rep_u ~ rep_w or rep_u ~ sigma rep_w)
    SA = {u: [w for w in range(12) if w != u and (((adj[reps[u]] >> reps[w]) & 1) or ((adj[reps[u]] >> sig(reps[w])) & 1))] for u in range(12)}
    Wlabs = {k: [u for u in range(12) if N[u, 12 + k] != 0 or ((adj[39 + k] >> reps[u]) & 1) or ((adj[39 + k] >> sig(reps[u])) & 1)] for k in range(5)}
    res = {}
    for rad in (0, 1, 2):
        for nW in (1, 2):
            kill = False
            for Ws in itertools.combinations(range(5), nW):
                S = set(u for k in Ws for u in Wlabs[k])
                for _ in range(rad):
                    S = S | set(w for u in S for w in SA[u])
                idx = sorted(S) + [12 + k for k in Ws]
                if np.linalg.eigvalsh(M[np.ix_(idx, idx)])[0] < -1e-9: kill = True; break
            stats[(rad, nW, kill)] += 1
for rad in (0, 1, 2):
    for nW in (1, 2):
        k = stats[(rad, nW, True)]; t = k + stats[(rad, nW, False)]
        print('radius %d around %d W(s): kill fraction %.3f' % (rad, nW, k / t))
