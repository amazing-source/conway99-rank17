import itertools, collections, numpy as np
pts7 = range(7)
# all Fano planes on 7 points
triples = list(itertools.combinations(pts7, 3))
def is_fano(lines):
    seen = set()
    for l in lines:
        for p in itertools.combinations(l, 2):
            if p in seen: return False
            seen.add(p)
    return len(seen) == 21
base = [(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]
assert is_fano(base)
def act(perm, plane): return frozenset(frozenset(perm[x] for x in l) for l in plane)
planes = set()
for perm in itertools.permutations(range(7)):
    planes.add(act(perm, base))
planes = list(planes); print('fano planes', len(planes))
def parity(p):
    p = list(p); s = 0
    for i in range(7):
        for j in range(i+1, 7):
            if p[i] > p[j]: s ^= 1
    return s
# A7 orbit of base
A7 = [p for p in itertools.permutations(range(7)) if parity(p) == 0]
orb = list({act(p, base) for p in A7}); print('A7 orbit', len(orb))
idx = {pl: i for i, pl in enumerate(orb)}
# lines of PG(3,2): for each triple, the planes containing it
lines = []
for t in triples:
    L = [i for i, pl in enumerate(orb) if frozenset(t) in pl]
    lines.append(tuple(L))
print('line sizes', collections.Counter(len(l) for l in lines))
# Hamming code = span of lines over F2
vecs = []
for L in lines:
    v = 0
    for i in L: v |= 1 << i
    vecs.append(v)
basis = []
for v in vecs:
    for b in basis: v = min(v, v ^ b)
    if v: basis.append(v)
print('dim span lines', len(basis))
words = {0}
for b in basis: words |= {w ^ b for w in words}
wd = collections.Counter(bin(w).count('1') for w in words); print('weights', sorted(wd.items()))
w5 = [w for w in words if bin(w).count('1') == 5]
# A7 action on the 15 planes
def perm15(p): return [idx[act(p, orb[i])] for i in range(15)]
gens = [(1,2,0,3,4,5,6), (0,1,3,4,5,6,2)]   # 3-cycle and 5-cycle generate A7
G15 = [perm15(g) for g in gens]
def appw(P, w):
    r = 0
    for i in range(15):
        if w >> i & 1: r |= 1 << P[i]
    return r
seen = set(); orbits = []
for w in w5:
    if w in seen: continue
    O = {w}; fr = [w]
    while fr:
        nf = []
        for x in fr:
            for P in G15:
                y = appw(P, x)
                if y not in O: O.add(y); nf.append(y)
        fr = nf
    seen |= O; orbits.append(sorted(O))
print('A7 orbits on weight-5 words', [len(o) for o in orbits])
import pickle
pickle.dump({'orb': orb, 'lines': lines, 'w5orbits': orbits, 'A7': A7, 'idx': idx}, open('pg32.pkl', 'wb'))
