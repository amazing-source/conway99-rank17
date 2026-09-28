import numpy as np, random
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')[:1000]
from collections import Counter
first_fail = Counter(); first_fail_1W = Counter(); first_fail_2W = Counter()
for ln in lines:
    rho, adj, h = load(ln); M, N, R = MU(rho, adj)
    # order star labels: those adjacent (N != 0) to W0 first, then to W1, ... then the rest
    order = []
    for k in range(5):
        for u in range(12):
            if N[u, 12 + k] != 0 and u not in order: order.append(u)
    for u in range(12):
        if u not in order: order.append(u)
    ff = None
    for s in range(1, 13):
        idx = order[:s] + list(range(12, 17))
        if np.linalg.eigvalsh(M[np.ix_(idx, idx)])[0] < -1e-9: ff = s; break
    first_fail[ff] += 1
    # all star + single W / pairs of W
    one = any(np.linalg.eigvalsh(M[np.ix_(list(range(12)) + [12 + k], list(range(12)) + [12 + k])])[0] < -1e-9 for k in range(5))
    first_fail_1W[one] += 1
print('first failing prefix size (star labels, ordered by W-adjacency) with all 5 W:', sorted(first_fail.items(), key=lambda t: (t[0] is None, t[0] or 0)))
print('all 12 star + 1 W fails:', dict(first_fail_1W))
