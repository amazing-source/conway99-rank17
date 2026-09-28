import sys, time, itertools, collections, numpy as np
sys.path.insert(0, '.')
import wfe
tag, pid = sys.argv[1], int(sys.argv[2])
t0 = time.time()
case = wfe.Wcase(tag, pid)
wn = case.wn; nw = len(wn)
# F2 rank of witnesses mod 2
M = (case.wx % 2).astype(np.uint8)
def f2rank(M):
    M = M.copy() % 2; r = 0
    for c in range(M.shape[1]):
        piv = [i for i in range(r, M.shape[0]) if M[i, c]]
        if not piv: continue
        M[[r, piv[0]]] = M[[piv[0], r]]
        for i in range(M.shape[0]):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
    return r
print('witness norms', wn, 'F2 rank', f2rank(M), flush=True)
# valid vectors in usable classes
U = set(case.usable.keys())
inU = np.isin(case.cls, list(U))
V = case.W[inU].astype(np.int64) @ case.wx.T         # values (nvec, nw)
print('valid vectors in usable classes', len(V), flush=True)
par = V & 1
# per pair (t,t2): nonzero product sets per type
def m_of(n, k):
    return (4 * n - 2 * k * (7 - k)) // 4
cuts = [p for p in range(64)]                       # p over lines 2..7 (bits 0..5), line 1 = 0
def cutvec(p): return [0] + [(p >> i) & 1 for i in range(6)]
CELLS = [(c, d) for c in range(7) for d in range(c + 1, 7)]
allowed_cuts = {}
for n in set(wn):
    allowed_cuts[n] = [p for p in range(64) if min(bin(p).count('1'), 7 - bin(p).count('1')) in wfe.K_ALLOWED[n]]
print({n: len(v) for n, v in allowed_cuts.items()})
def types(p, q):
    P, Q = cutvec(p), cutvec(q)
    cnt = collections.Counter()
    for c, d in CELLS:
        cnt[(P[c] ^ P[d], Q[c] ^ Q[d])] += 2
    return cnt
TY = {}
for p in range(64):
    for q in range(64):
        TY[(p, q)] = types(p, q)
surv = {}
tot_pairs = 0; tot_surv = 0
for t in range(nw):
    for t2 in range(t + 1, nw):
        n1, n2 = wn[t], wn[t2]
        target = int(4 * (case.wx[t] @ case.D @ case.wx[t2])) // 7
        prod = V[:, t] * V[:, t2]
        Ps = {}
        for ty in [(1, 1), (1, 0), (0, 1), (0, 0)]:
            sel = (par[:, t] == ty[0]) & (par[:, t2] == ty[1])
            pr = prod[sel]
            Ps[ty] = sorted(set(pr[pr != 0].tolist()))
        ok = set()
        for p in allowed_cuts[n1]:
            k1 = min(bin(p).count('1'), 7 - bin(p).count('1'))
            for q in allowed_cuts[n2]:
                k2 = min(bin(q).count('1'), 7 - bin(q).count('1'))
                cnt = TY[(p, q)]
                n11, n10, n01, n00 = cnt[(1, 1)], cnt[(1, 0)], cnt[(0, 1)], cnt[(0, 0)]
                if 7 in (n1, n2):
                    # conservative: every label anything in its type's product set (plus 0 unless type 11)
                    lo = hi = 0
                    for ty, nn in cnt.items():
                        vals = Ps[ty] + ([] if ty == (1, 1) else [0])
                        if not vals: lo = 10**9; break
                        lo += nn * min(vals); hi += nn * max(vals)
                    if lo <= target <= hi: ok.add((p, q))
                    continue
                m1, m2 = m_of(n1, k1), m_of(n2, k2)
                P11 = Ps[(1, 1)]
                if n11 and not P11: continue
                lo = n11 * min(P11) if n11 else 0; hi = n11 * max(P11) if n11 else 0
                feas = True
                # type 10: x10 labels (odd t, |a_t2| = 2), range [max(0, m2 - n00), min(m2, n10)]
                for (nn, mm, Pz) in ((n10, m2, Ps[(1, 0)]), (n01, m1, Ps[(0, 1)])):
                    xlo, xhi = max(0, mm - n00), min(mm, nn)
                    if xlo > xhi: feas = False; break
                    if not Pz:
                        if xlo > 0: feas = False; break
                        continue
                    lo += min(xlo * min(Pz), xhi * min(Pz)); hi += max(xlo * max(Pz), xhi * max(Pz))
                if not feas: continue
                P00 = Ps[(0, 0)]
                zhi = min(m1, m2, n00)
                if P00:
                    lo += min(0, zhi * min(P00)); hi += max(0, zhi * max(P00))
                if lo <= target <= hi: ok.add((p, q))
        tot_pairs += len(allowed_cuts[n1]) * len(allowed_cuts[n2]); tot_surv += len(ok)
        surv[(t, t2)] = ok
print('cut pairs total', tot_pairs, 'surviving', tot_surv, round(time.time() - t0, 1), 's', flush=True)
# per witness, cuts that have support with every other witness
dom = {t: set(allowed_cuts[wn[t]]) for t in range(nw)}
changed = True
while changed:
    changed = False
    for (t, t2), ok in surv.items():
        a = {p for p in dom[t] if any((p, q) in ok for q in dom[t2])}
        b = {q for q in dom[t2] if any((p, q) in ok for p in dom[t])}
        if a != dom[t] or b != dom[t2]:
            dom[t], dom[t2] = a, b; changed = True
print('arc-consistent domain sizes', [len(dom[t]) for t in range(nw)], flush=True)
import pickle
pickle.dump({'surv': surv, 'dom': dom, 'wn': wn, 'wx': case.wx}, open('cutcsp_%s_%d.pkl' % (tag, pid), 'wb'))
