import pickle, itertools, collections, numpy as np
d = pickle.load(open('pg32.pkl', 'rb'))
orb, A7, idx = d['orb'], d['A7'], d['idx']
O42 = [o for o in d['w5orbits'] if len(o) == 42][0]
def act(perm, plane): return frozenset(frozenset(perm[x] for x in l) for l in plane)
def perm15(p): return [idx[act(p, orb[i])] for i in range(15)]
def appw(P, w):
    r = 0
    for i in range(15):
        if w >> i & 1: r |= 1 << P[i]
    return r
# stabilizer of the first ovoid in A7 and its orbits on the 7 points
w0 = O42[0]
stab = [p for p in A7 if appw(perm15(p), w0) == w0]
print('stabilizer order', len(stab))
orbs7 = []
seen = set()
for x in range(7):
    if x in seen: continue
    o = {p[x] for p in stab}; seen |= o; orbs7.append(sorted(o))
print('stabilizer orbits on 7 points', orbs7)
# label each ovoid by the image of the stabilizer structure
fixed = [o for o in orbs7 if len(o) == 1]
print('fixed points', fixed)
# label via transport: find for each ovoid a group element mapping w0 to it
lab = {}
for p in A7:
    w = appw(perm15(p), w0)
    key = tuple(p[o[0]] for o in orbs7 if len(o) == 1)
    lab.setdefault(w, set()).add(key)
print('labels per ovoid', collections.Counter(len(v) for v in lab.values()))
labels = {w: sorted(v)[0] for w, v in lab.items()}
print('distinct labels', len(set(labels.values())), list(labels.values())[:6])
# intersection sizes by label relation
M = np.array([[bin(a & b).count('1') for b in O42] for a in O42])
rel = collections.defaultdict(collections.Counter)
L = [labels[w] for w in O42]
for i in range(42):
    for j in range(42):
        if i == j: continue
        a, b = L[i], L[j]
        if len(a) == 2:
            if set(a) == set(b): r = 'sister(reversed)'
            elif a[0] == b[0]: r = 'same first'
            elif a[1] == b[1]: r = 'same second'
            elif a[0] == b[1] or a[1] == b[0]: r = 'head-tail'
            else: r = 'disjoint'
        else:
            r = 'same pt' if a == b else 'diff pt'
        rel[r][M[i, j]] += 1
for r, c in rel.items(): print(r, dict(c))
# point coverage
cov = collections.Counter()
for w in O42:
    for i in range(15):
        if w >> i & 1: cov[i] += 1
print('coverage per point', set(cov.values()))
pickle.dump({'O42': O42, 'labels': labels}, open('ov42.pkl', 'wb'))
