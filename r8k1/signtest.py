"""Q,P-independent sign-constrained bounds on dumped leaves:
 even: F0(f) = 27|f|^2 - 9 f^T A0 f + (sum f)^2 with star-star edges removed, f even, f >= 0 on star roots  (form <= F0)
 odd : f odd, f >= 0 on O_c+ (so <= 0 on O_c-): Q edges drop, P edges bounded by 9 sum_star f^2."""
import numpy as np
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')[:600]
NVv = 49
def orbits():
    ev = {}; od = {}; rep = {}
    ne = no = 0
    for v in range(NVv):
        s = sig(v)
        if v == 0: ev[0] = ne; ne += 1; continue
        if v < s: ev[v] = ev[s] = ne; ne += 1; od[v] = od[s] = no; no += 1; rep[v] = 1; rep[s] = -1
    return ev, od, rep, ne, no
EV, OD, REP, NE, NO = orbits()
def forms(adj, drop_star=True):
    A = np.zeros((NVv, NVv))
    for a in range(NVv):
        for b in range(NVv):
            if (adj[a] >> b) & 1: A[a, b] = 1
    if drop_star:
        for a in range(15, 39):
            for b in range(15, 39): A[a, b] = 0
    P = 27*np.eye(NVv) - 9*A + np.ones((NVv, NVv))
    Be = np.zeros((NVv, NE)); Bo = np.zeros((NVv, NO))
    for v in range(NVv):
        Be[v, EV[v]] = 1
        if v != 0: Bo[v, OD[v]] = REP[v]
    return Be.T @ P @ Be, Bo.T @ P @ Bo
star_e = sorted({EV[v] for v in range(15, 39)}); star_o_plus = sorted({OD[v] for v in range(15, 27)})
def min_signed(Fm, starset, extra_diag=0.0, iters=300):
    """heuristic min of z^T F z / |z|^2 (orbit coords, weight ignored) with z >= 0 on starset (orbit coords where rep sign fixed)"""
    F = Fm.copy()
    for i in starset: F[i, i] += extra_diag
    n = F.shape[0]; w, V = np.linalg.eigh(F); z = V[:, 0].copy()
    best = 1e9
    for sgn in (1, -1):
        z = sgn * V[:, 0].copy()
        for it in range(iters):
            z[starset] = np.maximum(z[starset], 0)
            nz = np.linalg.norm(z)
            if nz == 0: break
            z /= nz
            val = z @ F @ z
            best = min(best, val)
            z = z - 0.05 * (F @ z - val * z)
    return best
ev_res = []; od_res = []
for ln in lines:
    rho, adj, h = load(ln)
    Fe, Fo = forms(adj, True)
    ev_res.append(min_signed(Fe, star_e, 0.0))
    # odd: orbit rep of star orbits in O_c+ is the O_c+ vertex (v < sig(v)) -> REP=+1 ; z>=0 there; add 9*2 on star orbit diag (P edges, both vertices of orbit)
    od_res.append(min_signed(Fo, star_o_plus, 18.0))
e = np.array(ev_res); o = np.array(od_res)
print('even sign-constrained (f>=0 on star, no correction): kills %.3f (max %.3f)' % (np.mean(e < -1e-9), e.max()))
print('odd  sign-constrained (f>=0 on O_c+, P-correction):  kills %.3f (max %.3f)' % (np.mean(o < -1e-9), o.max()))
print('either:', np.mean((e < -1e-9) | (o < -1e-9)))
