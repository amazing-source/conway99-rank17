"""Lattice data (not a search): root system of N = sqrt2 L* for one representative of each surviving class of L.
Roots of N = v/sqrt2 for v in U with |v|^2 = 4 and class in I."""
import numpy as np, sympy, collections
import evenpipe as ep, skel15, one
reps = [('Z15', (255, 3855, 13107, 21845)), ('E8Z7', (769, 3074, 12308, 21860)), ('Z15', (15, 1011, 7285, 25782)),
        ('E8Z7', (769, 1284, 6152, 26934)), ('E8Z7', (769, 1284, 6152, 24606)), ('E8Z7', (769, 3073, 12290, 21812)),
        ('Z15', (15, 240, 3891, 28981)), ('E8Z7', (769, 1284, 2310, 28941)), ('E8Z7', (769, 1284, 6152, 24584)),
        ('Z15', (15, 51, 2005, 30805)), ('E8Z7', (769, 3073, 12289, 21772)), ('Z15', (15, 51, 960, 31765)),
        ('Z15', (15, 51, 85, 32640))]
def typ(rank, nroots):
    for name, f in (('A', lambda r: r * (r + 1)), ('D', lambda r: 2 * r * (r - 1))):
        if f(rank) == nroots and not (name == 'D' and rank < 4): return '%s%d' % (name, rank)
    return {(6, 72): 'E6', (7, 126): 'E7', (8, 240): 'E8'}.get((rank, nroots), '?(%d,%d)' % (rank, nroots))
cache = {}
for nm, Ib in reps:
    if nm not in cache:
        L = ep.Lat(nm); L.S4 = ep.shell(L.G, 4); cache[nm] = L
    L = cache[nm]
    I = skel15.span(list(Ib))
    c4 = ep.cls(L.S4)
    R = L.S4[np.isin(c4, list(I))]
    R = R[one.firstpos(R)]
    G = R @ L.G @ R.T // 2                      # Gram of the roots of N (norm 2)
    parent = list(range(len(R)))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a in range(len(R)):
        for b in np.flatnonzero(G[a, a + 1:] != 0) + a + 1:
            ra, rb = find(a), find(int(b))
            if ra != rb: parent[rb] = ra
    comps = collections.defaultdict(list)
    for a in range(len(R)): comps[find(a)].append(a)
    desc = []
    for idx in comps.values():
        rk = sympy.Matrix(R[idx].tolist()).rank()
        desc.append(typ(rk, 2 * len(idx)))
    tot = sum(int(d[1:]) for d in desc if d[0] in 'ADE')
    rl = sorted(len(b) for b in one.root_components(L, Ib))
    print('R(L) comps', rl, '| R(N) =', ' + '.join(sorted(desc)), '| rank', tot, '| roots of N', 2 * len(R))
