import numpy as np, sys
from collections import Counter
C = 6
def sig(v):
    if v == 0: return 0
    if v < 15: return v + 1 if (v - 1) % 2 == 0 else v - 1
    if v < 39: w = v - 15; al = w // 12; r = w % 12; return 15 + (1 - al)*12 + (r ^ 1)
    return v + 5 if v < 44 else v - 5
def star_root(v):
    w = v - 15; al = 1 if w < 12 else -1; w %= 12; return al, w // 2, (1 if w % 2 == 0 else -1)
def load(line):
    x = line.split()
    rho = [list(map(int, x[6*k:6*k+6])) + [0] for k in range(5)]
    adj = [int(t, 16) for t in x[30:79]]
    h = list(map(int, x[79:103]))
    return rho, adj, h
def root_vec(v, rho):
    r = [0]*7
    if 15 <= v < 39:
        al, j, be = star_root(v); r[C] = al; r[j] = be
    elif 39 <= v < 44: r = rho[v-39][:]
    elif v >= 44: r = [-t for t in rho[v-44]]
    return r
def MU(rho, adj):
    reps = [15 + u for u in range(12)] + [39 + k for k in range(5)]   # star label u=(j,s): rep e_c + s e_j = vertex 15+u
    A = lambda a, b: (adj[a] >> b) & 1
    N = np.zeros((17, 17)); R = np.array([root_vec(p, rho) for p in reps], float)
    for i, p in enumerate(reps):
        for j, q in enumerate(reps):
            if i != j: N[i, j] = A(p, sig(q)) - A(p, q)
    return 4*N + 12*np.eye(17) - R @ R.T, N, R
lines = open(sys.argv[1]).read().strip().split('\n')
mins = []; wparts = []; rk = Counter()
th = 4*np.pi/5
for ln in lines:
    rho, adj, h = load(ln)
    M, N, R = MU(rho, adj)
    w, V = np.linalg.eigh(M)
    mins.append(w[0]); rk[int(np.sum(np.abs(w) < 1e-8))] += 1
    v = V[:, 0]
    wparts.append(np.linalg.norm(v[12:])**2)
print('n', len(lines), 'M_U min eig: max %.3f min %.3f mean %.3f' % (max(mins), min(mins), np.mean(mins)))
print('nullity distribution', dict(rk))
print('W-part weight of min eigvec: min %.3f mean %.3f' % (min(wparts), np.mean(wparts)))

print('--- Schur analysis ---')
seen = Counter(); out = []
for ln in lines[:3000]:
    rho, adj, h = load(ln)
    M, N, R = MU(rho, adj)
    MSS, MSW, MWW = M[:12,:12], M[:12,12:], M[12:,12:]
    F = MSW.T @ np.linalg.pinv(MSS) @ MSW
    Sig = MWW - F
    key = tuple(np.round(np.linalg.eigvalsh(MWW), 3))
    seen[key] += 1
    crude = MWW - MSW.T @ MSW / 20
    out.append((np.linalg.eigvalsh(Sig)[0], np.linalg.eigvalsh(crude)[0], np.linalg.eigvalsh(MSS)[-1]))
o = np.array(out)
print('M_WW spectra seen:', dict(seen))
print('Sigma min eig: max %.3f | crude T0 min eig: max %.3f frac<0 %.3f | lambda_max(M_SS) range %.2f..%.2f' % (o[:,0].max(), o[:,1].max(), np.mean(o[:,1] < -1e-9), o[:,2].min(), o[:,2].max()))

print('--- single-direction response test (x = -m/mu), y over M_WW low eigvecs and Sigma min eigvec ---')
res = []
for ln in lines[:3000]:
    rho, adj, h = load(ln)
    M, N, R = MU(rho, adj)
    MSS, MSW, MWW = M[:12,:12], M[:12,12:], M[12:,12:]
    wv, Vw = np.linalg.eigh(MWW)
    ys = [Vw[:, 0], Vw[:, 1]] + [np.cos(t)*Vw[:, 0] + np.sin(t)*Vw[:, 1] for t in np.linspace(0, np.pi, 13)]
    best = 1e9
    for y in ys:
        m = MSW @ y
        A2 = m @ MSS @ m; B = m @ m
        val = y @ MWW @ y - B*B/A2
        best = min(best, val)
    res.append(best)
r = np.array(res)
print('single-direction (true N_SS) kills fraction %.3f, max value %.3f' % (np.mean(r < -1e-9), r.max()))
