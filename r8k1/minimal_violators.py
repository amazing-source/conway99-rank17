"""Greedy minimal label subsets U' (star labels + W labels) with M_{U'} not PSD, over dumped leaves.
Also include inner lines implicitly (M_U already has them Schur-complemented)."""
import numpy as np, random
from collections import Counter
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')
random.seed(5)
sizes = Counter(); comp = Counter(); examples = []
for ln in random.sample(lines, 300):
    rho, adj, h = load(ln); M, N, R = MU(rho, adj)
    best = None
    for trial in range(8):
        S = list(range(17)); random.shuffle(S)
        cur = list(range(17))
        for x in S:
            t = [y for y in cur if y != x]
            if np.linalg.eigvalsh(M[np.ix_(t, t)])[0] < -1e-9: cur = t
        if best is None or len(cur) < len(best): best = cur
    sizes[len(best)] += 1
    nW = sum(1 for x in best if x >= 12)
    comp[(len(best) - nW, nW)] += 1
    if len(examples) < 5: examples.append((best, M[np.ix_(best, best)]))
print('minimal violating subset sizes:', sorted(sizes.items()))
print('(star labels, W labels) composition:', sorted(comp.items()))
for b, m in examples[:3]:
    print('example', b); print(np.round(m).astype(int))
