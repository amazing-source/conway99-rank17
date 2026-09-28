import sys, time, collections, numpy as np
sys.path.insert(0, '.')
import wfe
import rank17_skeleton as sk
from rank17_ypairs import short_vectors_chunked
tag, pid = sys.argv[1], int(sys.argv[2])
case = wfe.Wcase(tag, pid)
GK = case.GK
X, nrm = short_vectors_chunked(GK.tolist(), 2)
R = X[nrm == 2].astype(np.int64)
print('roots of K (up to sign):', len(R))
rc = sorted(set(int(((r & 1) << np.arange(15)).sum()) for r in R))
print('distinct root classes mod 2:', len(rc))
orb = sk.root_orbits(GK, list(case.usable.keys()))
cnt = collections.Counter(orb[x] for x in case.usable)
print('usable classes', len(case.usable), 'H-orbits on usable', len(cnt), 'sizes', sorted(cnt.values(), reverse=True)[:30])
# root system type: Gram of roots, connected components of the Dynkin-ish graph via a simple basis is complex; just count
G = R @ GK @ R.T
print('roots inner product histogram', collections.Counter(np.abs(G[np.triu_indices(len(R), 1)]).tolist()))
