import sys, json, numpy as np, collections
sys.path.insert(0, '.')
import wfe
from rank17_vsearch import CELLS
case = wfe.Wcase('E8', 5)
S, ex, nn = case.skeletons(max_count=400, limit=10**6)
shown = 0
for s in S:
    r = case.vector_level(s)
    if r.get('why', '').startswith('budget'):
        t, t2 = map(int, r['why'].split()[1].split(','))
        print('skel', s, r['why'], 'n', case.wn[t], case.wn[t2], 'sizes1', r['sizes1'])
        # recompute contributions
        doms = case.domains(s)
        g = [0] + list(s)
        p = np.array([[bin(gc & xm).count('1') & 1 for gc in g] for xm in case.xmask])
        o = p.sum(1); k = np.minimum(o, 7 - o)
        print('  k', k[t], k[t2], 'p_t', p[t].tolist(), 'p_t2', p[t2].tolist())
        target = int(4 * (case.wx[t] @ case.D @ case.wx[t2])) // 7
        lo = hi = 0
        for ci, (c, d) in enumerate(CELLS):
            Pp, Mn = doms[ci]
            vp = (Pp.astype(np.int64) @ case.T) // 49; vm = (Mn.astype(np.int64) @ case.T) // 49
            cut = p[:, c - 1] ^ p[:, d - 1]
            keep = np.ones(len(Pp), bool)
            for tt in range(len(case.wn)):
                ev, od = wfe.label_allowed(case.wn[tt], k[tt]); al = list(od if cut[tt] else ev)
                keep &= np.isin(np.abs(vp[:, tt]), al) & np.isin(np.abs(vm[:, tt]), al)
                if case.wn[tt] == 6 and k[tt] == 2 and not cut[tt]:
                    keep &= ~((np.abs(vp[:, tt]) == 2) & (np.abs(vm[:, tt]) == 2))
            cc = vp[keep, t] * vp[keep, t2] + vm[keep, t] * vm[keep, t2]
            lo += cc.min(); hi += cc.max()
            print('   cell', (c, d), 'cut', cut[t], cut[t2], 'n', keep.sum(), 'contrib', dict(collections.Counter(cc.tolist())),
                  'vals_t', sorted(set(np.abs(vp[keep, t]).tolist() + np.abs(vm[keep, t]).tolist())),
                  'vals_t2', sorted(set(np.abs(vp[keep, t2]).tolist() + np.abs(vm[keep, t2]).tolist())))
        print('  target', target, 'lo', lo, 'hi', hi)
        shown += 1
        if shown >= 2: break
