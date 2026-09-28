import numpy as np, itertools
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')[:600]
from collections import Counter
res = {s: 0 for s in range(1, 6)}; res0 = {s: 0 for s in range(1, 6)}
for ln in lines:
    rho, adj, h = load(ln); M, N, R = MU(rho, adj)
    # Q matrix on star labels: Q edges are rep-rep adjacencies within O_c+ : N_SS entry -1 from A(p_u,p_w)
    reps = [15 + u for u in range(12)]
    PiQ = np.zeros((12, 12))
    for i in range(12):
        for j in range(12):
            if i != j and (adj[reps[i]] >> reps[j]) & 1: PiQ[i, j] = 1
    L = np.zeros((12, 12))
    for u in range(12): L[u, u ^ 1] = 1
    KQ = 15*np.eye(12) + L - 4*PiQ - np.ones((12, 12))
    K0 = 20*np.eye(12) - np.ones((12, 12))
    for s in range(1, 6):
        killQ = kill0 = False
        for Ws in itertools.combinations(range(5), s):
            idx = list(range(12)) + [12 + k for k in Ws]
            B = M[np.ix_(idx, idx)].copy(); B[:12, :12] = KQ
            if np.linalg.eigvalsh(B)[0] < -1e-9: killQ = True
            B0 = M[np.ix_(idx, idx)].copy(); B0[:12, :12] = K0
            if np.linalg.eigvalsh(B0)[0] < -1e-9: kill0 = True
        res[s] += killQ; res0[s] += kill0
print('fraction killed by K_Q bound using the best subset of s W labels:', {s: round(res[s]/len(lines), 3) for s in res})
print('fraction killed by K_0 bound:', {s: round(res0[s]/len(lines), 3) for s in res0})
res = {s: 0 for s in range(1, 6)}
for ln in lines:
    rho, adj, h = load(ln); M, N, R = MU(rho, adj)
    reps = [15 + u for u in range(12)]
    Pipi = np.zeros((12, 12))
    for i in range(12):
        for j in range(12):
            if i != j and (adj[reps[i]] >> sig(reps[j])) & 1: Pipi[i, j] = 1
    L = np.zeros((12, 12))
    for u in range(12): L[u, u ^ 1] = 1
    Kt = 15*np.eye(12) + L + 4*Pipi - np.ones((12, 12))
    for s in range(1, 6):
        kill = False
        for Ws in itertools.combinations(range(5), s):
            idx = list(range(12)) + [12 + k for k in Ws]
            B = M[np.ix_(idx, idx)].copy(); B[:12, :12] = Kt
            if np.linalg.eigvalsh(B)[0] < -1e-9: kill = True; break
        res[s] += kill
print('fraction killed by K_tau (pi known, Q unknown):', {s: round(res[s]/len(lines), 3) for s in res})
