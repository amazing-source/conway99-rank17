import sys, time, collections, numpy as np
sys.path.insert(0, '.')
import wfe
tag, pid = sys.argv[1], int(sys.argv[2]); NS = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
t0 = time.time()
case = wfe.Wcase(tag, pid)
nw = len(case.wn)
pairs = [(t, t2) for t in range(nw) for t2 in range(t, nw)]
PT = np.array(pairs)
target = np.array([int(4 * (case.wx[t] @ case.D @ case.wx[t2])) // 7 for t, t2 in pairs])
W = case.W.astype(np.int64); cls = case.cls
order = np.argsort(cls, kind='stable'); W = W[order]; cls = cls[order]
b = np.flatnonzero(np.diff(cls)) + 1; st = np.concatenate([[0], b]); en = np.concatenate([b, [len(cls)]])
GK = case.GK.astype(np.int64); G2 = GK % 2; pw = 1 << np.arange(15)
TAB = {}   # class -> (masks array, mn (nm, npairs), mx)
npair_tot = 0
for s, e in zip(st, en):
    C = int(cls[s])
    if C not in case.usable: continue
    Wc = W[s:e]; vals = Wc @ case.wx.T
    ip = Wc @ GK @ Wc.T
    I, J = np.nonzero(np.isin(np.abs(ip), (0, 4)) & (np.arange(len(Wc))[:, None] < np.arange(len(Wc))[None, :]))
    if len(I) == 0: continue
    npair_tot += len(I)
    A = vals[I]; B = vals[J]
    con = A[:, PT[:, 0]] * A[:, PT[:, 1]] + B[:, PT[:, 0]] * B[:, PT[:, 1]]
    m1 = (((Wc[I] + Wc[J]) // 2 @ G2) % 2) @ pw
    m2 = (((Wc[I] - Wc[J]) // 2 @ G2) % 2) @ pw
    M = np.concatenate([m1, m2]); CON = np.concatenate([con, con])
    um, inv = np.unique(M, return_inverse=True)
    mn = np.full((len(um), len(pairs)), 10**6); mx = np.full((len(um), len(pairs)), -10**6)
    np.minimum.at(mn, inv, CON); np.maximum.at(mx, inv, CON)
    TAB[C] = (um, mn, mx)
print('tables', len(TAB), 'pairs', npair_tot, 'avg masks/class', np.mean([len(v[0]) for v in TAB.values()]), round(time.time() - t0, 1), 's', flush=True)
S, ex, nn = case.skeletons(max_count=NS, limit=10**7)
res = collections.Counter()
t1 = time.time()
for s in S:
    g = [0] + list(s)           # g[0] = line 1
    G = np.array(s, dtype=np.int64)
    lo = np.zeros(len(pairs), dtype=np.int64); hi = np.zeros(len(pairs), dtype=np.int64); dead = None
    for c in range(7):
        for d in range(c + 1, 7):
            C = g[c] ^ g[d]
            um, mn, mx = TAB[C]
            par = (np.bitwise_count(um[:, None] & G[None, :]) & 1)
            eps = np.ones(6, dtype=np.int64) if c == 0 else (np.arange(2, 8) == c + 1).astype(np.int64)
            ok = (par == eps).all(1)
            if not ok.any():
                dead = 'M1-empty'; break
            lo += mn[ok].min(0); hi += mx[ok].max(0)
        if dead: break
    if dead:
        res[dead] += 1; continue
    bad = (lo > target) | (hi < target)
    res['budget kill' if bad.any() else 'pass'] += 1
print(dict(res), 'per skeleton ms', round(1000 * (time.time() - t1) / max(1, len(S)), 2))
# ---- kill statistics per witness pair
killers = collections.Counter(); per_s = []
for s in S:
    g = [0] + list(s); G = np.array(s, dtype=np.int64)
    lo = np.zeros(len(pairs), dtype=np.int64); hi = np.zeros(len(pairs), dtype=np.int64)
    for c in range(7):
        for d in range(c + 1, 7):
            um, mn, mx = TAB[g[c] ^ g[d]]
            par = (np.bitwise_count(um[:, None] & G[None, :]) & 1)
            eps = np.ones(6, dtype=np.int64) if c == 0 else (np.arange(2, 8) == c + 1).astype(np.int64)
            ok = (par == eps).all(1)
            lo += mn[ok].min(0); hi += mx[ok].max(0)
    bad = np.flatnonzero((lo > target) | (hi < target))
    per_s.append(set(bad.tolist()))
    for b_ in bad: killers[b_] += 1
print('kills per skeleton (min/mean)', min(len(x) for x in per_s), np.mean([len(x) for x in per_s]))
top = killers.most_common(25)
print('top pairs', [(pairs[k], (case.wn[pairs[k][0]], case.wn[pairs[k][1]]), v) for k, v in top])
# greedy cover
unc = list(range(len(per_s))); chosen = []
while unc:
    cnt = collections.Counter(b_ for i in unc for b_ in per_s[i])
    k, v = cnt.most_common(1)[0]; chosen.append(k)
    unc = [i for i in unc if k not in per_s[i]]
print('greedy cover size', len(chosen), [(pairs[k], (case.wn[pairs[k][0]], case.wn[pairs[k][1]])) for k in chosen])
diag = [i for i, (t, t2) in enumerate(pairs) if t == t2]
print('skeletons killed by a diagonal (Parseval) pair', sum(1 for x in per_s if x & set(diag)))
