"""Interlacing on Y0 = {x} u N(x) u W2 only (fully determined by skeleton + W-W state): odd and even parts of
27I - 9A + J (PSD, rank) and 44I + 11A - 2J (PSD).  Run over all (skeleton, ww) with delta in {2,3}."""
import numpy as np, itertools
from skel import shapes, signed_skeletons
import model
from step_c import ww_states
sk = [s for sh in shapes() for s in signed_skeletons(sh)]
def forms(rhos, st):
    # vertices: x, 14 inner, 10 W2
    V = ['x'] + [('i', t, s) for t in range(7) for s in (1, -1)] + [('v', k) for k in range(5)] + [('s', k) for k in range(5)]
    idx = {v: i for i, v in enumerate(V)}; n = len(V)
    A = np.zeros((n, n))
    def add(a, b): A[idx[a], idx[b]] = A[idx[b], idx[a]] = 1
    for t in range(7): add('x', ('i', t, 1)); add('x', ('i', t, -1)); add(('i', t, 1), ('i', t, -1))
    for k in range(5):
        for t in range(6):
            if rhos[k][t]: add(('v', k), ('i', t, rhos[k][t])); add(('s', k), ('i', t, -rhos[k][t]))
    for (k, l), s in st.items():
        if s in (1, 3): add(('v', k), ('s', l)); add(('v', l), ('s', k))
        if s in (2, 3): add(('v', k), ('v', l)); add(('s', k), ('s', l))
    sigma = {}
    for v in V:
        if v == 'x': sigma[v] = v
        elif v[0] == 'i': sigma[v] = ('i', v[1], -v[2])
        elif v[0] == 'v': sigma[v] = ('s', v[1])
        else: sigma[v] = ('v', v[1])
    Sg = np.zeros((n, n))
    for v in V: Sg[idx[sigma[v]], idx[v]] = 1
    w, U = np.linalg.eigh((np.eye(n) + Sg) / 2); Be = U[:, w > .5]
    w, U = np.linalg.eigh((np.eye(n) - Sg) / 2); Bo = U[:, w > .5]
    P = 27*np.eye(n) - 9*A + np.ones((n, n)); Q = 44*np.eye(n) + 11*A - 2*np.ones((n, n))
    return [np.linalg.eigvalsh(B.T @ M @ B).min() for M in (P, Q) for B in (Be, Bo)]
from collections import Counter
st_all = 0; kill = Counter()
for s_i, rhos in enumerate(sk):
    for st in ww_states(rhos):
        d = [sum(1 for (k, l), v in st.items() if k == m or l == m if v == 1) - sum(1 for (k, l), v in st.items() if (k == m or l == m) and v == 2) for m in range(5)]
        al = [sum(1 for (k, l), v in st.items() if (m in (k, l)) and v == 1) - sum(1 for (k, l), v in st.items() if (m in (k, l)) and v == 2) for m in range(5)]
        if not all(2 <= x <= 3 for x in al): continue
        st_all += 1
        f = forms(rhos, st)
        kill[tuple(x < -1e-9 for x in f)] += 1
print('W-W states with delta in {2,3}:', st_all)
print('(P4 even fails, P4 odd fails, P3 even fails, P3 odd fails):', dict(kill))
