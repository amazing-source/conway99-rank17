import numpy as np, sys, pickle, collections
import evenpipe as ep, skel15, one
par = lambda a, b: np.bitwise_count(a & b) & 1
tot = collections.Counter(); ex = {}
for nm in ['Z15', 'E8Z7', 'D12Z3', 'E7E7Z', 'A15p']:
    L = ep.Lat(nm); L.S3 = ep.shell(L.G, 3); L.S4 = ep.shell(L.G, 4)
    Is = pickle.load(open('I_%s.pkl' % nm, 'rb'))
    cl = {k: ep.cls(S) if len(S) else np.zeros(0, dtype=np.int64) for k, S in ((1, L.S1), (2, L.S2), (3, L.S3), (4, L.S4), (5, L.S5))}
    for Ib in Is:
        inv = []
        for k in (1, 2, 3, 4, 5):
            ok = np.ones(len(cl[k]), dtype=bool)
            for i in Ib: ok &= par(cl[k], int(L.FM[i])) == 0
            inv.append(int(ok.sum()))
        # root system components (ranks) of L
        rcomp = sorted(len(b) for b in one.root_components(L, Ib))
        key = (tuple(inv), tuple(rcomp))
        tot[key] += 1; ex.setdefault(key, (nm, Ib))
print('distinct invariant classes of L:', len(tot))
for key, c in sorted(tot.items()):
    print(key[0], 'root comps', key[1], 'count', c, 'e.g.', ex[key])
