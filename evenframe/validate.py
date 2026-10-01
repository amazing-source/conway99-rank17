"""Independent Python count of valid partial assignments (first d orbits) vs evenframe.c depth profile."""
import numpy as np, sys, pickle, struct, subprocess, os, itertools
import evenpipe as ep, skel15
nm, gi, depth = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
L = ep.Lat(nm)
good = pickle.load(open('skel_%s.pkl' % nm, 'rb'))
Ib, omit, sizes = good[gi]
info = skel15.analyse(L, list(Ib)); cc = skel15.cell_classes(L, info, omit)
scls = ep.cls(L.S5)
cells = [(c, d) for c in range(7) for d in range(c + 1, 7)]
def firstpos(X):
    nz = X != 0; f = X[np.arange(len(X)), nz.argmax(axis=1)]; return f > 0
dom = []
for p in cells:
    X = L.S5[np.isin(scls, cc[p])]
    order = np.argsort(scls[np.isin(scls, cc[p])], kind='stable')   # same order as runframes: sorted by class
    # replicate runframes ordering exactly: vstack over classes in cc order of S sorted by class
    X = np.vstack([L.S5[scls == v] for v in cc[p]])
    P = X[firstpos(X)]; dom.append(np.vstack([P, -P]))
path = 'val_%s.bin' % nm
with open(path, 'wb') as f:
    f.write(struct.pack('<i', 15)); f.write(L.G.astype('<i4').tobytes()); f.write(struct.pack('<ii', 1000, 1))
    for X in dom: f.write(struct.pack('<i', len(X))); f.write(X.astype(np.int8).tobytes())
out = subprocess.run(['./evenframe2.exe', path], capture_output=True, text=True, env=dict(os.environ, MAXDEPTH=str(depth))).stdout
cprof = list(map(int, out.split('|')[1].split()))[:depth + 1]
# ---- independent reference
G = L.G
Lc = np.linalg.cholesky(G.astype(float))
orb = [(x, s) for x in range(21) for s in (0, 1)]
def ov(x, y): return len(set(cells[x]) & set(cells[y]))
def qof(s, ip):
    b = 1 - s - ip
    return {0: 0, 2: 1, 4: 2}.get(b, -1)
IP = {}
def ipm(x, y):
    if (x, y) not in IP: IP[(x, y)] = dom[x] @ G @ dom[y].T
    return IP[(x, y)]
counts = [0] * (depth + 1)
def valid(assign):
    k = len(assign)
    # pairwise + D + line counts
    Q = {}
    for i in range(k):
        for j in range(i):
            q = qof(ov(orb[i][0], orb[j][0]), int(ipm(orb[i][0], orb[j][0])[assign[i], assign[j]]))
            if q < 0: return False
            Q[(i, j)] = Q[(j, i)] = q
    for i in range(k):
        if sum(1 for j in range(k) if j != i and Q[(i, j)] == 2) > 1: return False
        for c in range(7):
            tgt = 2 if c in cells[orb[i][0]] else 4
            cur = sum(Q[(i, j)] for j in range(k) if j != i and c in cells[orb[j][0]])
            rem = sum(1 for j in range(k, 42) if c in cells[orb[j][0]])
            if cur > tgt or cur + 2 * rem < tgt: return False
    # sister order
    for i in range(k):
        if orb[i][1] == 1 and assign[i] <= assign[i - 1]: return False
    # PSD
    R = 14 * np.eye(15)
    for i in range(k):
        p = Lc.T @ dom[orb[i][0]][assign[i]]; R -= np.outer(p, p)
    if np.linalg.eigvalsh(R).min() < -1e-7: return False
    return True
def fc_ok(assign):
    """mirror of the count-aware forward checking: a future candidate is live iff it could be chosen now w.r.t. counts"""
    k = len(assign)
    Q = {}
    for i in range(k):
        for j in range(i):
            Q[(i, j)] = Q[(j, i)] = qof(ov(orb[i][0], orb[j][0]), int(ipm(orb[i][0], orb[j][0])[assign[i], assign[j]]))
    lineQ = [[sum(Q[(i, j)] for j in range(k) if j != i and c in cells[orb[j][0]]) for c in range(7)] for i in range(k)]
    Dc = [any(Q[(i, j)] == 2 for j in range(k) if j != i) for i in range(k)]
    for j in range(k, 42):
        x = orb[j][0]; nq = len(dom[x]); ok = np.ones(nq, dtype=bool)
        qs = np.zeros((nq, k), dtype=np.int64)
        for i in range(k):
            s = ov(orb[i][0], x); ips = ipm(x, orb[i][0])[:, assign[i]]
            b = 1 - s - ips
            ok &= np.isin(b, (0, 2, 4)); qs[:, i] = b // 2
            if orb[j][1] == 1 and orb[i] == (x, 0): ok &= np.arange(nq) > assign[i]
        ok &= (qs == 2).sum(axis=1) <= 1
        for i in range(k):
            if Dc[i]: ok &= qs[:, i] != 2
            for c in cells[x]:
                tgtU = 2 if c in cells[orb[i][0]] else 4
                ok &= lineQ[i][c] + qs[:, i] <= tgtU
        for c in range(7):
            tgt = 2 if c in cells[x] else 4
            cur = sum(qs[:, i] for i in range(k) if c in cells[orb[i][0]]) if k else np.zeros(nq, dtype=np.int64)
            rem = sum(1 for jj in range(k, 42) if jj != j and c in cells[orb[jj][0]])
            ok &= (cur <= tgt) & (cur + 2 * rem >= tgt)
        if not ok.any(): return False
    return True
def rec(assign):
    k = len(assign); counts[k] += 1
    if k == depth: return
    x = orb[k][0]
    for t in range(len(dom[x])):
        if k == 0 and t >= len(dom[x]) // 2: continue
        a = assign + [t]
        if valid(a) and fc_ok(a): rec(a)
rec([])
print('C     ', cprof); print('Python', counts); print('MATCH' if cprof == counts else 'MISMATCH')
