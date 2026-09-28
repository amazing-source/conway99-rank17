import sys, time, numpy as np
sys.path.insert(0, '.')
import wfe
t0 = time.time()
case = wfe.Wcase('E8', 5)
print('init', round(time.time() - t0, 1), case.stats, flush=True)
t1 = time.time()
S, ex, nn = case.skeletons(max_count=10**7, limit=5*10**7)
print('skeletons', len(S), 'exhaustive', ex, 'nodes', nn, 's', round(time.time() - t1, 1), flush=True)
import collections
print('masks sizes', sorted(collections.Counter(len(m) for m in case.masks.values()).items())[:20])
import pickle
pickle.dump(S, open('E8_5_skel.pkl', 'wb'))
