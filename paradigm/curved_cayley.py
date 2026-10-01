"""Counterfactual experiment: does a CURVED (1,2)-graph exist anywhere?
(1,2)-graph: every edge in exactly one triangle, every two vertices at distance 2 have exactly 2 common neighbours.
Flat = locally Hamming (wedge matching at every vertex = the mate map).  All known examples (H(n,3), quotients, BvLS) are flat.
Here: Cayley graphs on class-2 exponent-3 groups G = V x Z, (v,z)(w,y) = (v+w, z+y+beta(v,w)), connection set = n subgroups
of order 3.  Exhaustive over all connection sets (increasing subgroup index), tests at the identity.
For every (1,2) hit: is <D> abelian, and is the graph flat (mu_d(e) = e^-1 for all wedges)?"""
import itertools, sys, collections
def make(m, k, mats):
    def mul(a, b):
        v = tuple((x + y) % 3 for x, y in zip(a[0], b[0]))
        z = tuple((a[1][j] + b[1][j] + sum(a[0][p] * mats[j][p][q] * b[0][q] for p in range(m) for q in range(m))) % 3
                  for j in range(k))
        return (v, z)
    els = [(v, z) for v in itertools.product(range(3), repeat=m) for z in itertools.product(range(3), repeat=k)]
    idx = {e: i for i, e in enumerate(els)}
    N = len(els)
    M = [[idx[mul(a, b)] for b in els] for a in els]
    e0 = idx[((0,) * m, (0,) * k)]
    inv = [next(j for j in range(N) if M[i][j] == e0) for i in range(N)]
    return N, M, inv, e0
def E(m, i, j):
    A = [[0] * m for _ in range(m)]; A[i][j] = 1; return A
def add(A, B): return [[x + y for x, y in zip(r, s)] for r, s in zip(A, B)]
GROUPS = {
 'Heis27':   (2, 1, [E(2, 0, 1)]),
 'Heis27x3': (3, 1, [E(3, 0, 1)]),
 '3^(1+4)':  (4, 1, [add(E(4, 0, 1), E(4, 2, 3))]),
 'Heisx9':   (4, 1, [E(4, 0, 1)]),
 'V3Z2':     (3, 2, [E(3, 0, 1), E(3, 0, 2)]),
}
def run(name, nmax, cap):
    m, k, mats = GROUPS[name]
    N, M, inv, e0 = make(m, k, mats)
    subs = sorted({min(g, inv[g]) for g in range(N) if g != e0})
    res = collections.Counter(); ex = {}; nodes = [0]
    def flat(D):
        Ds = set(D)
        for d in D:
            for d2 in D:
                if d2 in (d, inv[d]): continue
                g = M[d][d2]; h = M[d][inv[d2]]
                e = [a for a in D if a != d and M[inv[a]][g] in Ds]
                e2 = [a for a in D if a != d and M[inv[a]][h] in Ds]
                if len(e) != 1 or len(e2) != 1 or e2[0] != inv[e[0]]: return False
        return True
    def rec(start, D, cnt):
        nodes[0] += 1
        if nodes[0] > cap: return
        n = len(D) // 2
        if n >= 2:
            Ds = set(D)
            if all((c == 1 if g in Ds else c == 2) for g, c in cnt.items() if g != e0):
                ab = all(M[a][b] == M[b][a] for a in D for b in D)
                # order of <D>
                seen = {e0}; fr = [e0]
                while fr:
                    nf = []
                    for x in fr:
                        for d in D:
                            y = M[x][d]
                            if y not in seen: seen.add(y); nf.append(y)
                    fr = nf
                key = (n, len(seen), ab, flat(D))
                res[key] += 1
                ex.setdefault(key, list(D))
        if n == nmax: return
        for i in range(start, len(subs)):
            d = subs[i]; di = inv[d]
            if d in D: continue
            D2 = D + [d, di]; Ds = set(D2); c2 = dict(cnt); ok = True
            for a in D2:
                for b in (d, di):
                    for g in ({M[a][b], M[b][a]} if a not in (d, di) else {M[a][b]}):
                        pass
            # recount from scratch (small)
            c2 = collections.Counter(M[a][b] for a in D2 for b in D2)
            for g, c in c2.items():
                if g == e0: continue
                if (g in Ds and c > 1) or c > 2: ok = False; break
            if ok: rec(i + 1, D2, c2)
    rec(0, [], collections.Counter())
    print(name, 'order', N, 'subgroups', len(subs), 'nodes', nodes[0], '(capped)' if nodes[0] > cap else '')
    for key in sorted(res): print('   n', key[0], '|<D>|', key[1], 'abelian', key[2], 'flat', key[3], 'count', res[key], 'e.g.', ex[key][:6])
for name, nmax, cap in [('Heis27', 3, 10**6), ('Heis27x3', 6, 3 * 10**6), ('V3Z2', 5, 3 * 10**6), ('3^(1+4)', 4, 3 * 10**6), ('Heisx9', 4, 3 * 10**6)]:
    run(name, nmax, cap); sys.stdout.flush()
