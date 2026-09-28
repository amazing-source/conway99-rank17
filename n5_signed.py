"""V-free test of Sol's boxed n = 5 statement (structures, 28/09).

Claim to test: there is no graph-origin signed matrix N on the 42 labels u = (cell {i<j}, s = +/-), r_u = e_i + s e_j, with
  (S1) for every vertex z: N restricted to star(z), in the local orientation sig_z(u) = (r_u)_z, equals Q_z - P_z for two perfect
       matchings P_z, Q_z of the 12 labels of star(z)                                  [Sol, PROVEN from lambda = 1, mu = 2]
  (S2) for every vertex z and label u outside star(z): its row into star(z), sig_z(v) N_uv, equals B - A for 2-subsets A, B of
       star(z) (A = neighbours in X, B = neighbours in sigma X)                        [Sol, PROVEN]
  (S3) every label has exactly 10 nonzero entries                                      [diagonal of N^2 = N + 12I - G_r]
  (N2) optional: N^2 = N + 12 I - G_r entrywise                                        [frame identity]
carrying a vector a with a = +-1 on the 20 labels of the cut {0,1} | {2..6} and 0 elsewhere, N a = 4 a     [Sol n = 5, k = 2]
If the model is infeasible, every norm-5/7 dual vector is impossible, so min K* >= 6/7.
Usage: python n5_signed.py [--n2] [--time T] [--workers W]"""
import sys, itertools, time
from ortools.sat.python import cp_model

use_n2 = '--n2' in sys.argv
no_a = '--noa' in sys.argv
T = float(sys.argv[sys.argv.index('--time') + 1]) if '--time' in sys.argv else 300.0
W = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 2
cut = (0, 1)

cells = [(i, j) for i in range(7) for j in range(i + 1, 7)]
labels = [(c, s) for c in cells for s in (1, -1)]           # r = e_i + s e_j
L = len(labels)


def r(u):
    (i, j), s = labels[u]
    v = [0] * 7; v[i] = 1; v[j] = s
    return v


R = [r(u) for u in range(L)]
star = {z: [u for u in range(L) if R[u][z] != 0] for z in range(7)}
m = cp_model.CpModel()
pos, neg = {}, {}
for u in range(L):
    for v in range(u + 1, L):
        pos[u, v] = m.NewBoolVar(f'p{u}_{v}'); neg[u, v] = m.NewBoolVar(f'm{u}_{v}')
        m.AddBoolOr([pos[u, v].Not(), neg[u, v].Not()])


def Nexpr(u, v):
    a, b = min(u, v), max(u, v)
    return pos[a, b] - neg[a, b]


def absN(u, v):
    a, b = min(u, v), max(u, v)
    return pos[a, b] + neg[a, b]


# (S1) stars: perfect matchings P_z, Q_z; sig sig N = Q - P on star pairs
for z in range(7):
    S = star[z]
    P = {e: m.NewBoolVar(f'P{z}_{e}') for e in itertools.combinations(S, 2)}
    Q = {e: m.NewBoolVar(f'Q{z}_{e}') for e in itertools.combinations(S, 2)}
    for u in S:
        m.Add(sum(P[e] for e in P if u in e) == 1)
        m.Add(sum(Q[e] for e in Q if u in e) == 1)
    for (u, v) in P:
        sg = R[u][z] * R[v][z]
        m.Add(sg * Nexpr(u, v) == Q[u, v] - P[u, v])
# (S2) outside rows: sig_z(v) N_uv = B_v - A_v with |A| = |B| = 2
for z in range(7):
    S = star[z]
    for u in range(L):
        if R[u][z] != 0:
            continue
        A = {v: m.NewBoolVar(f'A{z}_{u}_{v}') for v in S}
        B = {v: m.NewBoolVar(f'B{z}_{u}_{v}') for v in S}
        m.Add(sum(A.values()) == 2); m.Add(sum(B.values()) == 2)
        for v in S:
            m.Add(R[v][z] * Nexpr(u, v) == B[v] - A[v])
# (S3) degree 10
for u in range(L):
    m.Add(sum(absN(u, v) for v in range(L) if v != u) == 10)
# the cut vector a
C = [u for u in range(L) if (R[u][cut[0]] != 0) != (R[u][cut[1]] != 0)]
assert len(C) == 20
sgn = {u: m.NewBoolVar(f's{u}') for u in C}                  # a_u = 2 s_u - 1


def prod(u, v):
    """N_uv * a_v as a linear expression via AND variables (v in C)"""
    a, b = min(u, v), max(u, v)
    pp = m.NewBoolVar(''); mm = m.NewBoolVar('')
    m.AddMultiplicationEquality(pp, [pos[a, b], sgn[v]])
    m.AddMultiplicationEquality(mm, [neg[a, b], sgn[v]])
    return 2 * pp - 2 * mm - pos[a, b] + neg[a, b]


for u in (range(L) if not no_a else []):
    terms = [prod(u, v) for v in C if v != u]
    if u in C:
        m.Add(sum(terms) == 4 * (2 * sgn[u] - 1))
    else:
        m.Add(sum(terms) == 0)
# R^T a = 0 (implied, cheap)
for z in (range(7) if not no_a else []):
    m.Add(sum(R[u][z] * (2 * sgn[u] - 1) for u in C) == 0)
# symmetry: fix one sign (a -> -a)
m.Add(sgn[C[0]] == 1)
# (N2) optional frame identity N^2 = N + 12 I - G_r
if use_n2:
    for u in range(L):
        for v in range(u + 1, L):
            terms = []
            for w in range(L):
                if w == u or w == v:
                    continue
                # N_uw N_wv in {-1,0,1}: t = (p_uw p_wv + m_uw m_wv) - (p_uw m_wv + m_uw p_wv)
                a1, b1 = min(u, w), max(u, w); a2, b2 = min(w, v), max(w, v)
                pp = m.NewBoolVar(''); mm = m.NewBoolVar(''); pm = m.NewBoolVar(''); mp = m.NewBoolVar('')
                m.AddMultiplicationEquality(pp, [pos[a1, b1], pos[a2, b2]])
                m.AddMultiplicationEquality(mm, [neg[a1, b1], neg[a2, b2]])
                m.AddMultiplicationEquality(pm, [pos[a1, b1], neg[a2, b2]])
                m.AddMultiplicationEquality(mp, [neg[a1, b1], pos[a2, b2]])
                terms.append(pp + mm - pm - mp)
            rr = sum(R[u][k] * R[v][k] for k in range(7))
            m.Add(sum(terms) == Nexpr(u, v) - rr)
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = T
solver.parameters.num_workers = W
t0 = time.time()
st = solver.Solve(m)
print('n2' if use_n2 else 'no-n2', 'status', solver.StatusName(st), 'time %.1fs' % (time.time() - t0), flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    a = {u: 2 * solver.Value(sgn[u]) - 1 for u in C}
    print('a on C:', [a[u] for u in C])
    Nm = [[0] * L for _ in range(L)]
    for (u, v) in pos:
        Nm[u][v] = Nm[v][u] = solver.Value(pos[u, v]) - solver.Value(neg[u, v])
    import json
    json.dump({'N': Nm, 'a': a, 'labels': [[list(c), s] for c, s in labels]}, open('n5_solution.json', 'w'))
    print('saved n5_solution.json')
