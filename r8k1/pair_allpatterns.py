"""For leaf pairs (W_k, W_l) killed at radius 0: is the kill independent of the star-edge pattern inside S0?
Enumerate all (partial pi-matching, partial Q-matching) on S0 (labels adjacent to W_k or W_l) and test PSD."""
import numpy as np, itertools, random
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')
random.seed(3)
def partial_matchings(items):
    items = list(items)
    if not items: yield []; return
    a = items[0]; rest = items[1:]
    for pm in partial_matchings(rest): yield pm          # a unmatched
    for i, b in enumerate(rest):
        for pm in partial_matchings(rest[:i] + rest[i+1:]): yield [(a, b)] + pm
tested = 0; indep = 0; sizes = []
for ln in random.sample(lines, 60):
    rho, adj, h = load(ln); M, N, R = MU(rho, adj)
    reps = [15 + u for u in range(12)]
    Wl = {k: [u for u in range(12) if ((adj[39 + k] >> reps[u]) & 1) or ((adj[39 + k] >> sig(reps[u])) & 1)] for k in range(5)}
    for k, l in itertools.combinations(range(5), 2):
        S0 = sorted(set(Wl[k]) | set(Wl[l]))
        idx = S0 + [12 + k, 12 + l]
        base = M[np.ix_(idx, idx)]
        if np.linalg.eigvalsh(base)[0] >= -1e-9: continue      # only pairs killed in the true leaf
        tested += 1; sizes.append(len(S0))
        # true star block entries come from 4N_SS; replace by all patterns
        n0 = len(S0)
        RR = (np.array([root_vec(reps[u], rho) for u in S0], float))
        RRt = RR @ RR.T
        allfail = True
        pms = list(partial_matchings(range(n0)))
        for pp in pms:
            Npi = np.zeros((n0, n0))
            for a, b in pp: Npi[a, b] = Npi[b, a] = 1
            for pq in pms:
                NS = Npi.copy()
                for a, b in pq: NS[a, b] -= 1; NS[b, a] -= 1
                B = base.copy(); B[:n0, :n0] = 4 * NS + 12 * np.eye(n0) - RRt
                if np.linalg.eigvalsh(B)[0] >= -1e-9: allfail = False; break
            if not allfail: break
        indep += allfail
    if tested >= 40: break
print('pairs tested', tested, 'S0 sizes', sorted(set(sizes)), 'killed for ALL star patterns inside S0:', indep)
