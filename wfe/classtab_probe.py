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
# per usable class: valid pairs (|ip| in {0,4}) -> min/max contributions per witness pair
W = case.W.astype(np.int64); cls = case.cls
order = np.argsort(cls, kind='stable'); W = W[order]; cls = cls[order]
b = np.flatnonzero(np.diff(cls)) + 1; st = np.concatenate([[0], b]); en = np.concatenate([b, [len(cls)]])
GK = case.GK.astype(np.int64)
tabmn, tabmx = {}, {}
for s, e in zip(st, en):
    C = int(cls[s])
    if C not in case.usable: continue
    Wc = W[s:e]; vals = Wc @ case.wx.T
    ip = Wc @ GK @ Wc.T
    I, J = np.nonzero(np.isin(np.abs(ip), (0, 4)) & (np.arange(len(Wc))[:, None] < np.arange(len(Wc))[None, :]))
    if len(I) == 0: continue
    # contributions for each pair of vectors: a_t(I)a_t2(I) + a_t(J)a_t2(J)
    A = vals[I]; B = vals[J]
    con = A[:, PT[:, 0]] * A[:, PT[:, 1]] + B[:, PT[:, 0]] * B[:, PT[:, 1]]
    tabmn[C] = con.min(0); tabmx[C] = con.max(0)
print('tables', len(tabmn), round(time.time() - t0, 1), 's', flush=True)
S, ex, nn = case.skeletons(max_count=NS, limit=10**7)
kill_tab = 0; res = collections.Counter()
for s in S:
    g = [0] + list(s)
    cl = [g[d] for d in range(1, 7)] + [g[c] ^ g[d] for c in range(1, 7) for d in range(c + 1, 7)]
    if any(c not in tabmn for c in cl):
        res['tab: class without pair'] += 1; continue
    lo = sum(tabmn[c] for c in cl); hi = sum(tabmx[c] for c in cl)
    bad = (lo > target) | (hi < target)
    res['tab kill' if bad.any() else 'tab pass'] += 1
print(dict(res), round(time.time() - t0, 1), 's')
