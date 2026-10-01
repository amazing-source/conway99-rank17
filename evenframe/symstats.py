import numpy as np, sys, pickle, collections, itertools
import evenpipe as ep, skel15, one
nm = sys.argv[1]; idxs = [int(x) for x in sys.argv[2].split(',')]
L = ep.Lat(nm)
good = pickle.load(open('skel_%s.pkl' % nm, 'rb'))
scls = ep.cls(L.S5); rcls = ep.cls(L.S2)
par = lambda a, b: np.bitwise_count(a & b) & 1
for gi in idxs:
    Ib, omit, _ = good[gi]
    info = skel15.analyse(L, list(Ib)); I = np.array(info['I'])
    canon = lambda x: int(np.min(I ^ int(x)))
    pairof = {canon(v): k for k, v in info['pair'].items()}
    inC = np.ones(len(L.S2), dtype=bool)
    for i in Ib: inC &= par(rcls, int(L.FM[i])) == 0
    R = L.S2[inC]; Rc = rcls[inC]
    prs = collections.Counter(pairof.get(canon(c), None) for c in Rc)
    print(nm, gi, 'I', Ib, 'omit', omit, 'roots in L', len(R), 'pairs realised', len(prs), sorted(prs.items(), key=str)[:28])
    # components of the root graph on the 8 points
    comp = list(range(8))
    def f(a):
        while comp[a] != a: a = comp[a]
        return a
    for p in prs:
        if p is None: continue
        a, b = f(p[0]), f(p[1])
        if a != b: comp[b] = a
    print('   point components', collections.Counter(f(a) for a in range(8)).values())
