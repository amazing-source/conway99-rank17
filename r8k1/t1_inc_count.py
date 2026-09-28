"""Count T1-cyclic incidence structures (G+, G-): labeled 5-edge graphs on the six -1 roots / six +1 roots of O_c+,
degrees in {1,2,3}, line exclusions, mu budgets (aligned <= 2, non-aligned <= 1 common star nbrs), mu(y, v_k) <= 2."""
import itertools, sys
from collections import Counter
eps = list(map(int, sys.argv[1:7])) if len(sys.argv) > 6 else [1, 1, 1, -1, -1, -1]
C = 6
# O_c+ roots: (j, beta) meaning e_c + beta e_j ; h = eps_j
minus = [(j, b) for j in range(6) for b in (1, -1) if eps[j] == -1]
plus = [(j, b) for j in range(6) for b in (1, -1) if eps[j] == 1]
def excl_plus(k):   # G+ vertices (O_c+ -1 roots) excluded for v_k: (k,+), (k+1,-)
    return {(k, 1), ((k + 1) % 5, -1)}
def excl_minus(k):  # G- vertices (O_c+ +1 roots, = sigma of v_k's O_c- picks): excluded (k,-), (k+1,+)
    return {(k, -1), ((k + 1) % 5, 1)}
def graphs(verts, excl):
    E = list(itertools.combinations(verts, 2))
    out = []
    def rec(k, used, deg):
        if k == 5:
            if all(1 <= deg.get(v, 0) <= 3 for v in verts): out.append(tuple(used))
            return
        for e in E:
            if e in used: continue
            if e[0] in excl(k) or e[1] in excl(k): continue
            if deg.get(e[0], 0) >= 3 or deg.get(e[1], 0) >= 3: continue
            deg[e[0]] = deg.get(e[0], 0) + 1; deg[e[1]] = deg.get(e[1], 0) + 1
            rec(k + 1, used + [e], deg)
            deg[e[0]] -= 1; deg[e[1]] -= 1
    rec(0, [], {})
    return out
Gp = graphs(minus, excl_plus); Gm = graphs(plus, excl_minus)
print('eps', eps, '|G+|', len(Gp), '|G-|', len(Gm))
def share(g, k, l): return len(set(g[k]) & set(g[l]))
aligned = {(k, (k + 1) % 5) for k in range(5)} | {((k + 1) % 5, k) for k in range(5)}
# precompute sharing vectors
def svec(g): return tuple(share(g, k, l) for k in range(5) for l in range(k + 1, 5))
pairs = [(k, l) for k in range(5) for l in range(k + 1, 5)]
bud = [2 if (k, l) in aligned else 1 for (k, l) in pairs]
cp = Counter(svec(g) for g in Gp); cm = Counter(svec(g) for g in Gm)
tot = 0
for sa, na in cp.items():
    for sb, nb in cm.items():
        if all(sa[i] + sb[i] <= bud[i] for i in range(10)): tot += na * nb
print('pairs satisfying mu budgets:', tot)
