"""V-free feasibility of a short dual witness (structures, 28/09) — generalises n5_signed.py.

A dual vector kappa in K* with n = 7|kappa|^2 gives a = W kappa in Z^42 with
  |a|^2 = 4n,  N a = 4a,  R^T a = 0,  a_u odd exactly on the labels of the cells cut by a vertex set P (|P| = k, 0 <= k <= 3),
  and Sol's star bound  sum_{u in star(c)} a_u^2 <= floor(20 n / 7)  for every vertex c.
N is a graph-origin signed matrix: (S1) each star = Q_z - P_z (perfect matchings, local orientation (r_u)_z); (S2) every outside
row into a star = B - A with |A| = |B| = 2; (S3) degree 10; (N2, optional) N^2 = N + 12 I - G_r.
If the model is infeasible for (n, k), no K* vector of that norm and cut size exists in any graph-origin configuration.
Usage: python eigvec_signed.py N K [--n2] [--amax M] [--time T] [--workers W]"""
import sys, itertools, time, json
from ortools.sat.python import cp_model

n = int(sys.argv[1]); k = int(sys.argv[2])
use_n2 = '--n2' in sys.argv
amax = int(sys.argv[sys.argv.index('--amax') + 1]) if '--amax' in sys.argv else 3
T = float(sys.argv[sys.argv.index('--time') + 1]) if '--time' in sys.argv else 300.0
W = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 2
P = set(range(k))                                            # WLOG the cut set is {0..k-1}

cells = [(i, j) for i in range(7) for j in range(i + 1, 7)]
labels = [(c, s) for c in cells for s in (1, -1)]
L = len(labels)
R = []
for (i, j), s in labels:
    v = [0] * 7; v[i] = 1; v[j] = s; R.append(v)
star = {z: [u for u in range(L) if R[u][z] != 0] for z in range(7)}
cutlab = [((labels[u][0][0] in P) != (labels[u][0][1] in P)) for u in range(L)]

m = cp_model.CpModel()
pos, neg = {}, {}
for u in range(L):
    for v in range(u + 1, L):
        pos[u, v] = m.NewBoolVar(''); neg[u, v] = m.NewBoolVar('')
        m.AddBoolOr([pos[u, v].Not(), neg[u, v].Not()])
Nv = {}
for (u, v) in pos:
    x = m.NewIntVar(-1, 1, ''); m.Add(x == pos[u, v] - neg[u, v]); Nv[u, v] = Nv[v, u] = x


def NE(u, v):
    a, b = min(u, v), max(u, v)
    return pos[a, b] - neg[a, b]


for z in range(7):
    S = star[z]
    Pm = {e: m.NewBoolVar('') for e in itertools.combinations(S, 2)}
    Qm = {e: m.NewBoolVar('') for e in itertools.combinations(S, 2)}
    for u in S:
        m.Add(sum(Pm[e] for e in Pm if u in e) == 1)
        m.Add(sum(Qm[e] for e in Qm if u in e) == 1)
    for (u, v) in Pm:
        m.Add(R[u][z] * R[v][z] * NE(u, v) == Qm[u, v] - Pm[u, v])
for z in range(7):
    S = star[z]
    for u in range(L):
        if R[u][z] != 0:
            continue
        A = {v: m.NewBoolVar('') for v in S}; B = {v: m.NewBoolVar('') for v in S}
        m.Add(sum(A.values()) == 2); m.Add(sum(B.values()) == 2)
        for v in S:
            m.Add(R[v][z] * NE(u, v) == B[v] - A[v])
for u in range(L):
    m.Add(sum(pos[min(u, v), max(u, v)] + neg[min(u, v), max(u, v)] for v in range(L) if v != u) == 10)

# the witness a: parity by the cut, |a|^2 = 4n, star bound, R^T a = 0, N a = 4 a
a = {}
sq = {}
for u in range(L):
    if cutlab[u]:
        a[u] = m.NewIntVarFromDomain(cp_model.Domain.FromValues([x for x in range(-amax, amax + 1) if x % 2]), f'a{u}')
    else:
        a[u] = m.NewIntVarFromDomain(cp_model.Domain.FromValues([x for x in range(-amax + (amax % 2), amax + 1, 2)]), f'a{u}')
    sq[u] = m.NewIntVar(0, amax * amax, ''); m.AddMultiplicationEquality(sq[u], [a[u], a[u]])
m.Add(sum(sq.values()) == 4 * n)
for c in range(7):
    m.Add(sum(sq[u] for u in star[c]) <= (20 * n) // 7)
for z in range(7):
    m.Add(sum(R[u][z] * a[u] for u in range(L)) == 0)
for u in range(L):
    terms = []
    for v in range(L):
        if v == u:
            continue
        t = m.NewIntVar(-amax, amax, ''); m.AddMultiplicationEquality(t, [Nv[u, v], a[v]]); terms.append(t)
    m.Add(sum(terms) == 4 * a[u])
# symmetry: a -> -a (fix the sign of the first nonzero-forced label if any)
firstcut = next((u for u in range(L) if cutlab[u]), None)
if firstcut is not None:
    m.Add(a[firstcut] >= 1)
if use_n2:
    for u in range(L):
        for v in range(u + 1, L):
            terms = []
            for w in range(L):
                if w == u or w == v:
                    continue
                t = m.NewIntVar(-1, 1, ''); m.AddMultiplicationEquality(t, [Nv[u, w], Nv[w, v]]); terms.append(t)
            rr = sum(R[u][q] * R[v][q] for q in range(7))
            m.Add(sum(terms) == NE(u, v) - rr)
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = T
solver.parameters.num_workers = W
t0 = time.time()
st = solver.Solve(m)
print(f'n={n} k={k} n2={use_n2} amax={amax}: {solver.StatusName(st)} ({time.time() - t0:.1f}s)', flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    sol = {'a': [solver.Value(a[u]) for u in range(L)],
           'N': [[(solver.Value(Nv[u, v]) if u != v else 0) for v in range(L)] for u in range(L)]}
    json.dump(sol, open(f'eigvec_n{n}_k{k}_solution.json', 'w'))
    print('solution saved')
