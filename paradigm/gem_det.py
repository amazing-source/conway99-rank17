"""Experiment (law discovery, not a decision run): a star is three fixed-point-free involutions (nu, mu, psi) on 12 flags.
For every pair (mu, psi) with mu taken up to the symmetries of nu, compute
  - the integer det of L = 5 - 4 nu + 2 mu - 2 psi (its vanishing mod 7 is the odd rank drop of the star block),
  - the 'topology' of the 3-coloured cubic graph: components, bicoloured cycle counts, bipartiteness, Euler char.
Question: is det L (or det L mod 7) a function of the topological data?"""
import itertools, collections, sys
import numpy as np
from fractions import Fraction
def matchings(pts):
    if not pts: yield []; return
    a = pts[0]
    for i in range(1, len(pts)):
        b = pts[i]; rest = pts[1:i] + pts[i+1:]
        for m in matchings(rest): yield [(a, b)] + m
P12 = list(range(12))
nu = [i ^ 1 for i in range(12)]                    # flags 2j, 2j+1 are antipodal
allm = []
for m in matchings(P12):
    inv = [0] * 12
    for a, b in m: inv[a] = b; inv[b] = a
    allm.append(tuple(inv))
print('matchings', len(allm))
def cycles(a, b):
    """cycle lengths (in flags) of the 2-coloured graph a, b"""
    seen = [False] * 12; out = []
    for s in range(12):
        if seen[s]: continue
        n = 0; x = s
        while not seen[x]:
            seen[x] = True; x = a[x]; seen[x] = True; x = b[x]; n += 2
        out.append(n)
    return tuple(sorted(out))
def comps(a, b, c):
    par = list(range(12))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for inv in (a, b, c):
        for i in range(12): par[f(i)] = f(inv[i])
    return len({f(i) for i in range(12)})
def bip(a, b, c):
    col = [-1] * 12
    for s in range(12):
        if col[s] >= 0: continue
        col[s] = 0; st = [s]
        while st:
            x = st.pop()
            for inv in (a, b, c):
                y = inv[x]
                if col[y] < 0: col[y] = 1 - col[x]; st.append(y)
                elif col[y] == col[x] and y != x: return False
    return True
def Pm(inv):
    M = np.zeros((12, 12), dtype=np.int64)
    for i in range(12): M[i, inv[i]] = 1
    return M
Pnu = Pm(nu)
def det_int(M):
    from sympy import Matrix
    return int(Matrix(M.tolist()).det())
# mu up to the group of nu: one representative per cycle type of (nu, mu)
reps = {}
for m in allm:
    ct = cycles(nu, m)
    if ct not in reps: reps[ct] = m
print('mu classes (cycle types of nu.mu):', sorted(reps))
table = collections.defaultdict(collections.Counter)
tot = 0
for ct, mu in reps.items():
    Pmu = Pm(mu)
    for psi in allm:
        L = 5 * np.eye(12, dtype=np.int64) - 4 * Pnu + 2 * Pmu - 2 * Pm(psi)
        d = int(round(np.linalg.det(L.astype(float))))
        key = (ct, cycles(nu, psi), cycles(mu, psi), comps(nu, mu, psi), bip(nu, mu, psi))
        table[key][d] += 1; tot += 1
print('pairs', tot, 'topological keys', len(table))
single = sum(1 for k, v in table.items() if len(v) == 1)
print('keys with a single det value:', single)
mod7 = sum(1 for k, v in table.items() if len({d % 7 == 0 for d in v}) == 1)
print('keys on which [det = 0 mod 7] is constant:', mod7)
import pickle; pickle.dump(dict(table), open('gem_det.pkl', 'wb'))
# show a few keys with several values
shown = 0
for k, v in sorted(table.items(), key=lambda kv: -len(kv[1])):
    print(k, dict(list(v.items())[:8])); shown += 1
    if shown >= 6: break
