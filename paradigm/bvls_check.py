"""Check of hand-derived laws on the flat control srg(243,22,1,2) (n = 11), sigma = -1.
 (a) triangle chart: for a triangle T, f = #prisms on T, n0 = #triangles with no edge to T; law n0 = c(n) - f.
 (b) signed triangles of the odd layer N (natural sign): t+ - t- = -(n-2), t- = 2(n-2) - (2-s_U) - tip(U)."""
import itertools, numpy as np
g = [2, 0, 1, 2, 1, 1]            # g(x) = x^5 + x^4 + 2x^3 + x^2 + 2, coefficients low -> high
def xpow(i):
    v = [0] * 5; v[0] = 1
    for _ in range(i):
        c = v[4]; v = [0] + v[:4]
        v = [(v[j] - c * g[j]) % 3 for j in range(5)]
    return tuple(v)
H = [xpow(i) for i in range(11)]
V = list(itertools.product(range(3), repeat=5)); idx = {v: i for i, v in enumerate(V)}
D = H + [tuple((-a) % 3 for a in h) for h in H]
add = lambda a, b: tuple((x + y) % 3 for x, y in zip(a, b))
A = np.zeros((243, 243), dtype=np.int64)
for v in V:
    for d in D: A[idx[v], idx[add(v, d)]] = 1
A2 = A @ A
assert (A2 == 22 * np.eye(243, dtype=np.int64) + A + 2 * (1 - np.eye(243, dtype=np.int64) - A)).all(); print('srg(243,22,1,2) OK')
n = 11
nb = [set(np.nonzero(A[i])[0]) for i in range(243)]
tris = sorted({tuple(sorted((i, j, k))) for i in range(243) for j in nb[i] for k in nb[i] & nb[j]})
print('triangles', len(tris), 'expected', n * 243 // 3)
T = tris[0]; Ts = set(T)
f = n0 = 0; cnt = {0: 0, 1: 0, 2: 0, 3: 0}
for S in tris:
    if Ts & set(S): continue
    e = sum(1 for a in T for b in S if A[a, b]); cnt[e] += 1
print('edges-to-T distribution', cnt, ' law: n3=f, n2=6(n-1)-3f, n1=6(n-1)(n-3)+3f, n0=', n * (2 * n * n + 1) // 3 - 1 - 9 * (n - 1) - 6 * (n - 1) * (n - 3) - cnt[3])
# sigma = -1, x = 0
neg = lambda v: tuple((-a) % 3 for a in v)
x = idx[(0,) * 5]; units = nb[x]
roots = [i for i in range(243) if i != x and i not in units]
orb = sorted({min(i, idx[neg(V[i])]) for i in roots}); m = len(orb); pos = {u: k for k, u in enumerate(orb)}
N = np.zeros((m, m), dtype=np.int64)
for u in orb:
    for v in orb:
        if u != v: N[pos[u], pos[v]] = A[u, v] - A[u, idx[neg(V[v])]]
N3 = N @ N @ N
print('orbits', m, 'diag N^3 values', set(np.diag(N3)), ' law: -2(n-2) =', -2 * (n - 2))
tp = tm = 0
for a, b, c in itertools.combinations(range(m), 3):
    p = N[a, b] * N[b, c] * N[a, c]
    if p > 0: tp += 1
    elif p < 0: tm += 1
print('T+ =', tp, 'T- =', tm, ' law (flat, D = sisters): T+ = m(n-2)/3 =', m * (n - 2) // 3, ', T- = 2m(n-2)/3 =', 2 * m * (n - 2) // 3)
