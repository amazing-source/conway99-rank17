"""Same question as curved_cayley.py on groups that are far from abelian: V:3 (Frobenius-like), A5, PSL(2,7), A6, SL(2,3), A4xA4.
A hit = connection set of n subgroups of order 3 giving a (1,2)-graph; report |<D>|, abelian?, flat?"""
import itertools, sys, collections
def table_from_perms(gens):
    n = len(gens[0]); e = tuple(range(n)); els = [e]; idx = {e: 0}; fr = [e]
    while fr:
        nf = []
        for x in fr:
            for g in gens:
                y = tuple(x[i] for i in g)
                if y not in idx: idx[y] = len(els); els.append(y); nf.append(y)
        fr = nf
    M = [[idx[tuple(a[i] for i in b)] for b in els] for a in els]
    return M
def affine3(p, m, T):
    """V = F_p^m, T an m x m matrix of order 3; elements (v, i)"""
    def app(i, v):
        for _ in range(i): v = tuple(sum(T[r][c] * v[c] for c in range(m)) % p for r in range(m))
        return v
    els = [(v, i) for v in itertools.product(range(p), repeat=m) for i in range(3)]
    idx = {e: k for k, e in enumerate(els)}
    M = [[idx[(tuple((x + y) % p for x, y in zip(a[0], app(a[1], b[0]))), (a[1] + b[1]) % 3)] for b in els] for a in els]
    return M
def cyc(n, *cs):
    p = list(range(n))
    for c in cs:
        for i in range(len(c)): p[c[i]] = c[(i + 1) % len(c)]
    return tuple(p)
def run(name, M, cap):
    N = len(M); e0 = next(i for i in range(N) if all(M[i][j] == j for j in range(N)))
    inv = [next(j for j in range(N) if M[i][j] == e0) for i in range(N)]
    subs = sorted({min(g, inv[g]) for g in range(N) if g != e0 and M[M[g][g]][g] == e0})
    nmax = int(((N - 1) / 2) ** 0.5)
    res = collections.Counter(); ex = {}; nodes = [0]; deep = [0]
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
        n = len(D) // 2; deep[0] = max(deep[0], n)
        if n >= 2:
            Ds = set(D)
            if all((c == 1 if g in Ds else c == 2) for g, c in cnt.items() if g != e0):
                ab = all(M[a][b] == M[b][a] for a in D for b in D)
                seen = {e0}; fr = [e0]
                while fr:
                    nf = []
                    for x in fr:
                        for d in D:
                            y = M[x][d]
                            if y not in seen: seen.add(y); nf.append(y)
                    fr = nf
                key = (n, len(seen), ab, flat(D)); res[key] += 1; ex.setdefault(key, list(D))
        if n == nmax: return
        for i in range(start, len(subs)):
            d = subs[i]; di = inv[d]
            D2 = D + [d, di]; Ds = set(D2)
            c2 = collections.Counter(M[a][b] for a in D2 for b in D2); ok = True
            for g, c in c2.items():
                if g == e0: continue
                if (g in Ds and c > 1) or c > 2: ok = False; break
            if ok: rec(i + 1, D2, c2)
    # first subgroup: only one per start index (no symmetry reduction, capped)
    rec(0, [], collections.Counter())
    print(name, 'order', N, 'order-3 subgroups', len(subs), 'nmax', nmax, 'nodes', nodes[0], '(capped)' if nodes[0] > cap else '', 'max depth', deep[0])
    for key in sorted(res): print('   n', key[0], '|<D>|', key[1], 'abelian', key[2], 'flat', key[3], 'count', res[key], 'e.g.', ex[key])
    sys.stdout.flush()
cap = 1500000
for p, w in [(7, 2), (13, 3), (19, 7), (31, 5), (37, 10), (43, 6)]:
    run('%d:3' % p, affine3(p, 1, [[w]]), cap)
run('A4', affine3(2, 2, [[0, 1], [1, 1]]), cap)
run('2^4:3', affine3(2, 4, [[0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1]]), cap)
run('5^2:3', affine3(5, 2, [[0, 4], [1, 4]]), cap)
run('7^2:3', affine3(7, 2, [[2, 0], [0, 4]]), cap)
run('A5', table_from_perms([cyc(5, (0, 1, 2)), cyc(5, (2, 3, 4))]), cap)
run('PSL(2,7)', table_from_perms([cyc(7, (0, 1, 2, 3, 4, 5, 6)), cyc(7, (1, 2, 4), (3, 6, 5)), cyc(7, (0, 1), (2, 5))]), cap) if False else None
run('L3(2)', table_from_perms([cyc(7, (0, 1, 2, 3, 4, 5, 6)), cyc(7, (1, 2), (3, 5))]), cap)
run('A6', table_from_perms([cyc(6, (0, 1, 2)), cyc(6, (1, 2, 3)), cyc(6, (2, 3, 4)), cyc(6, (3, 4, 5))]), cap)
run('SL(2,3)', table_from_perms([cyc(8, (0, 1, 2, 3), (4, 5, 6, 7)), cyc(8, (0, 4, 2, 6), (1, 7, 3, 5)), cyc(8, (1, 4, 7), (3, 6, 5))]), cap)
