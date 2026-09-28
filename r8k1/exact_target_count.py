"""For a few leaf star structures: count W-incidence structures (5 unlabeled W's, each = (A+ pair, A-' pair) of O_c+ roots)
with EXACT column targets p_y - q_y = tau_y, the 10 pairs pairwise distinct and not Q-edges, per-W h-sum in {-4,-2,0,2,4}."""
import itertools, random
from collections import Counter
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')
random.seed(2)
for ln in random.sample(lines, 4):
    rho, adj, h = load(ln)
    A = lambda a, b: (adj[a] >> b) & 1
    Op = list(range(15, 27))
    hv = {v: h[v - 15] for v in range(15, 39)}
    Q = {y: [z for z in Op if z != y and A(y, z)][0] for y in Op}
    P = {y: [z for z in range(27, 39) if A(y, z)][0] for y in Op}
    tau = {y: (-4 * hv[y] - hv[Q[y]] - hv[P[y]]) // 2 for y in Op}   # p_y - q_y
    pairs = [p for p in itertools.combinations(Op, 2) if Q[p[0]] != p[1]]
    # W row: (plus pair, minus pair) with h-sum = h(plus) - h(minus) in {-4,..,4} (always true); count solutions by DFS
    need = dict(tau)
    sols = 0
    rows = [(a, b) for a in pairs for b in pairs if a != b]
    def rec(k, start, used, need):
        global sols
        if k == 5:
            if all(v == 0 for v in need.values()): sols += 1
            return
        # capacity prune: each remaining W changes each coordinate by at most +1 (plus) and -1 (minus)
        rem = 5 - k
        if sum(max(v, 0) for v in need.values()) > 2 * rem or sum(max(-v, 0) for v in need.values()) > 2 * rem: return
        for i in range(start, len(rows)):
            a, b = rows[i]
            if a in used or b in used: continue
            for y in a: need[y] -= 1
            for y in b: need[y] += 1
            rec(k + 1, i + 1, used | {a, b}, need)
            for y in a: need[y] += 1
            for y in b: need[y] -= 1
    rec(0, 0, frozenset(), need)
    print('targets', sorted(Counter(tau.values()).items()), 'unlabeled W-incidence structures with exact targets:', sols)
