import numpy as np, itertools
from collections import Counter
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')[:800]
dist = Counter(); dist1 = Counter()
for ln in lines:
    rho, adj, h = load(ln); M, N, R = MU(rho, adj)
    nf = 0
    for k, l in itertools.combinations(range(5), 2):
        idx = list(range(12)) + [12 + k, 12 + l]
        if np.linalg.eigvalsh(M[np.ix_(idx, idx)])[0] < -1e-9: nf += 1
    n1 = sum(1 for k in range(5) if np.linalg.eigvalsh(M[np.ix_(list(range(12)) + [12 + k], list(range(12)) + [12 + k])])[0] < -1e-9)
    dist[nf] += 1; dist1[n1] += 1
print('number of failing pairs (12 star + 2 W) per leaf:', sorted(dist.items()))
print('number of failing singles (12 star + 1 W) per leaf:', sorted(dist1.items()))
