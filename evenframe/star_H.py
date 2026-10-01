"""CHECK (small exact enumeration): does the rootless class admit even ONE star?
Star(0) = the 6 cells {0,d}; 12 signed ovoids, 2 per cell, with
  sisters: G in {-1,-3,-5};  different cells of the star: G in {0,-2,-4};  each vector: sum of B over the star = 2;
  sum of the 12 vectors = 0.
Sign changes at the 15 points are symmetries: the first vector is taken all +."""
import itertools, numpy as np, sys
exec(open('z7_H.py').read().split('# ---- Z7 action')[0].split('pat = collections.Counter()')[0])
star = [(0, d) for d in range(1, 7)]
cand = {}
for c in star:
    L = []
    for O in OV[c]:
        O = sorted(O)
        for signs in itertools.product((1, -1), repeat=5):
            v = np.zeros(15, dtype=np.int64); v[O] = signs; L.append(v)
    cand[c] = np.array(L)
# pairs per cell: (i<j) with G in {-1,-3,-5}
pairs = {}
for c in star:
    X = cand[c]; G = X @ X.T
    I, J = np.nonzero(np.isin(G, (-1, -3, -5)) & (np.arange(len(X))[:, None] < np.arange(len(X))[None, :]))
    pairs[c] = list(zip(I.tolist(), J.tolist()))
print('ordered sister pairs per cell', [len(pairs[c]) for c in star])
count = 0; nodes = 0
sols = []
def rec(k, chosen, vecs, total):
    """chosen: list of (cell, i, j); vecs: list of vectors"""
    global count, nodes
    nodes += 1
    if k == 6:
        if not total.any():
            V = np.array(vecs); G = V @ V.T
            B = np.zeros_like(G)
            for a in range(12):
                for b in range(12):
                    if a == b: continue
                    same = (a // 2 == b // 2)
                    B[a, b] = (-1 - G[a, b]) // 2 if same else (-G[a, b]) // 2
            if (B.sum(axis=1) == 2).all():
                count += 1
                if len(sols) < 3: sols.append(V.copy())
        return
    c = star[k]; X = cand[c]
    for (i, j) in pairs[c]:
        if k == 0 and not (X[i] >= 0).all(): continue          # sign normalisation: first vector all +
        v, w = X[i], X[j]
        ok = True
        for u in vecs:
            if int(u @ v) not in (0, -2, -4) or int(u @ w) not in (0, -2, -4): ok = False; break
        if not ok: continue
        t2 = total + v + w
        # remaining cells can change each coordinate by at most 2 per cell
        if (np.abs(t2) > 2 * (5 - k)).any(): continue
        rec(k + 1, chosen + [(c, i, j)], vecs + [v, w], t2)
rec(0, [], [], np.zeros(15, dtype=np.int64))
print('nodes', nodes, 'stars (first vector all +):', count)
if sols: np.save('starH_example.npy', sols[0]); print(sols[0])
