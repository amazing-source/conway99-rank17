import numpy as np
from collections import Counter
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')
cyc = Counter(); gam = Counter(); dedge = Counter(); samecell = Counter()
for ln in lines:
    rho, adj, h = load(ln)
    stars = list(range(15, 39))
    A = lambda a, b: (adj[a] >> b) & 1
    # vertex-level star graph components
    seen = set(); comps = []
    for v in stars:
        if v in seen: continue
        comp = []; stack = [v]
        while stack:
            x = stack.pop()
            if x in seen: continue
            seen.add(x); comp.append(x)
            for y in stars:
                if A(x, y) and y not in seen: stack.append(y)
        comps.append(len(comp))
    cyc[tuple(sorted(comps))] += 1
    g = tuple(sorted(sum(1 for w in range(39, 49) if A(v, w)) for v in stars))
    gam[g] += 1
    nd = sum(1 for v in range(15, 27) for w in range(15, 27) if v < w and A(v, w) and A(v, sig(w)))
    dedge[nd] += 1
print('vertex-level star cycle length multisets (top):', cyc.most_common(6))
print('gamma multisets (top):', gam.most_common(4))
print('D-edges inside star:', dict(dedge))
