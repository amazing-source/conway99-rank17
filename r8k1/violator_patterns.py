"""Collect minimal violating label subsets and canonicalise them by their integer Gram (M restricted) up to permutation."""
import numpy as np, random, itertools
from collections import Counter
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')
random.seed(11)
def canon(Mi):
    n = len(Mi); best = None
    for p in itertools.permutations(range(n)):
        t = tuple(int(round(Mi[p[i], p[j]])) for i in range(n) for j in range(n))
        if best is None or t < best: best = t
    return (n, best)
pats = Counter(); per_leaf_min = []
for ln in random.sample(lines, 250):
    rho, adj, h = load(ln); M, N, R = MU(rho, adj)
    found = set()
    for trial in range(12):
        S = list(range(17)); random.shuffle(S); cur = list(range(17))
        for x in S:
            t = [y for y in cur if y != x]
            if np.linalg.eigvalsh(M[np.ix_(t, t)])[0] < -1e-9: cur = t
        if len(cur) <= 7:
            found.add(canon(M[np.ix_(cur, cur)]))
    for f in found: pats[f] += 1
    per_leaf_min.append(len(found))
print('distinct minimal violating Gram patterns (size <= 7):', len(pats))
by = Counter(k[0] for k in pats)
print('by size:', dict(by))
print('leaves with at least one small (<=7) violator found:', sum(1 for x in per_leaf_min if x > 0), '/', len(per_leaf_min))
for (n, t), c in pats.most_common(6):
    print('count', c, 'size', n); print(np.array(t).reshape(n, n))
