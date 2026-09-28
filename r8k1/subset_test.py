import numpy as np, sys, itertools, random
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')[:400]
rng = random.Random(1)
for sz in (6, 8, 9, 10, 11, 12):
    fails = 0; tot = 0; best_fail = 0
    for ln in lines:
        rho, adj, h = load(ln)
        M, N, R = MU(rho, adj)
        # best subset of size sz (search over random subsets) including all W
        found = False
        subs = [tuple(sorted(rng.sample(range(12), sz))) for _ in range(40)] if sz < 12 else [tuple(range(12))]
        for S in subs:
            idx = list(S) + list(range(12, 17))
            if np.linalg.eigvalsh(M[np.ix_(idx, idx)])[0] < -1e-9: found = True; break
        tot += 1; fails += found
    print('star subset size', sz, '(+ all 5 W): fraction of leaves with SOME failing subset among 40 random:', fails / tot)
# subsets with fewer W
for nw in (2, 3, 4):
    fails = 0
    for ln in lines:
        rho, adj, h = load(ln); M, N, R = MU(rho, adj)
        found = False
        for Ws in itertools.combinations(range(5), nw):
            idx = list(range(12)) + [12 + k for k in Ws]
            if np.linalg.eigvalsh(M[np.ix_(idx, idx)])[0] < -1e-9: found = True; break
        fails += found
    print('all 12 star +', nw, 'W: fraction failing for some W-subset:', fails / len(lines))
