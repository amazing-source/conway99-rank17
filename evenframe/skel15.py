"""Per I: quotient W6, type, the 8 'points', skeletons (omitted point), per-cell class sets and domain sizes."""
import numpy as np, itertools, sys, pickle, collections
import evenpipe as ep

def span(vs):
    s = {0}
    for v in vs: s |= {x ^ v for x in s}
    return s

def analyse(L, Ib):
    FM, N4, u0 = L.FM, L.N4, L.u0
    I = sorted(span(Ib))
    allc = np.arange(1 << 15)
    par = lambda a, b: np.bitwise_count(a & b) & 1
    inC = np.ones(1 << 15, dtype=bool)
    for v in Ib: inC &= par(allc, int(FM[v])) == 0
    even = inC & (par(allc, int(FM[u0])) == 0)             # C cap u0^perp
    # representatives of (C cap u0^perp)/I : min of coset
    Iarr = np.array(I, dtype=np.int64)
    E = np.flatnonzero(even)
    rep = np.minimum.reduce([E ^ i for i in Iarr])
    reps = np.unique(rep)
    if len(reps) != 64: return None
    q = (N4[reps] // 2) % 2                                  # q = |v|^2/2 mod 2 (classes are even)
    assert (N4[reps] % 2 == 0).all()
    nons = reps[q == 1]
    if len(nons) != 28: return {'type': 'minus' if len(nons) == 36 else 'deg%d' % len(nons)}
    # triangular graph T(8) on the 28 nonsingular classes: adjacent iff b = 1
    def canon(x): return int(np.min(Iarr ^ x))
    adj = {int(a): {int(b) for b in nons if b != a and L.b(int(a), int(b))} for a in nons}
    # points = maximal cliques of size 7
    pts = []
    for a in nons:
        a = int(a)
        for b in adj[a]:
            common = adj[a] & adj[b]
            # a = {x,y}, b = {x,z}: common nbrs = pairs containing x (5 more) + {y,z}
            for c in common:
                cl = {a, b, c} | (adj[a] & adj[b] & adj[c])
                if len(cl) == 7 and all(v2 in adj[v1] for v1 in cl for v2 in cl if v1 != v2):
                    pts.append(frozenset(cl))
    pts = list(set(pts))
    if len(pts) != 8: return {'type': 'plus?', 'pts': len(pts)}
    # pair class of points i,j = the unique nonsingular class in both cliques
    pair = {}
    for i, j in itertools.combinations(range(8), 2):
        c = pts[i] & pts[j]; assert len(c) == 1
        pair[(i, j)] = next(iter(c))
    return {'type': 'plus', 'I': I, 'pair': pair}

def cell_classes(L, info, omit):
    """labelling: the 7 points other than `omit`, in increasing order, are lines 0..6; cell (c,d) -> 16 classes"""
    pts = [i for i in range(8) if i != omit]
    Iarr = np.array(info['I'], dtype=np.int64)
    out = {}
    for c in range(7):
        for d in range(c + 1, 7):
            w = info['pair'][(pts[c], pts[d])]
            out[(c, d)] = sorted(int(x) for x in (Iarr ^ w ^ L.u0))
    return out

if __name__ == '__main__':
    nm = sys.argv[1]
    L = ep.Lat(nm)
    Is = pickle.load(open('I_%s.pkl' % nm, 'rb'))
    scls = ep.cls(L.S5)
    cnt = np.bincount(scls, minlength=1 << 15)
    types = collections.Counter(); good = []
    for Ib in Is:
        info = analyse(L, list(Ib))
        if info is None: types['notrank'] += 1; continue
        types[info['type']] += 1
        if info['type'] != 'plus': continue
        for omit in range(8):
            cc = cell_classes(L, info, omit)
            sizes = [int(cnt[v].sum()) for v in cc.values()]
            # every cell needs two sister vectors; norm-5 classes must be nonempty
            if min(sizes) >= 2:
                good.append((Ib, omit, sizes))
    print(nm, 'I types', dict(types), 'skeletons with all cells nonempty:', len(good))
    for g in good[:10]: print('  I', g[0], 'omit', g[1], 'sizes', sorted(g[2]))
    pickle.dump(good, open('skel_%s.pkl' % nm, 'wb'))
