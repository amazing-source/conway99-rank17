"""n = 5, k = 2 witness, restricted to the two cut stars (structures, 28/09).

Cut P = {0,1}. Labels: C0 = cells {0,x} (x = 2..6, 10 labels), C1 = cells {1,x} (10), O = cell {0,1} (2), Z = cells {x,y} in {2..6}
(20). star(0) = C0 + O, star(1) = C1 + O. Every entry of N touching star(0) u star(1) is a variable, so the 22 labels there have
their full rows. Constraints, each PROVEN:
  - stars 0 and 1: local orientation sig_z(u) = (r_u)_z, N = Q_z - P_z with perfect matchings                  [Sol, S1]
  - outside rows: every label outside star(z) (z = 0,1) has sig_z(t) N_wt = B - A over star(z), |A| = |B| = 2      [Sol, S2]
  - degree exactly 10 for the 22 labels of star(0) u star(1)                                               [N^2 diagonal]
  - witness a = +-1 on C = C0 u C1, 0 elsewhere; (N a)_u = 4 a_u on C, 0 on O u Z; R^T a = 0             [n = 5, k = 2]
  - optional --starx: the star structure S1 also at the five outer vertices x (their 12 labels: cells {0,x},{1,x},{x,y})
This is a RELAXATION of the full graph-origin system (it drops constraints), so INFEASIBLE proves there is no n = 5 dual
witness at all; FEASIBLE only means the relaxation is too weak.
Usage: python n5_local.py [--starx] [--time T] [--workers W]"""
import sys, itertools, time, json
from ortools.sat.python import cp_model

starx = '--starx' in sys.argv
n2core = '--n2core' in sys.argv
s2all = '--s2all' in sys.argv
T = float(sys.argv[sys.argv.index('--time') + 1]) if '--time' in sys.argv else 600.0
Wk = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 4

cells = [(i, j) for i in range(7) for j in range(i + 1, 7)]
labels = [(c, s) for c in cells for s in (1, -1)]
L = len(labels)
R = []
for (i, j), s in labels:
    v = [0] * 7; v[i] = 1; v[j] = s; R.append(v)
star = {z: [u for u in range(L) if R[u][z] != 0] for z in range(7)}
C = [u for u in range(L) if (R[u][0] != 0) != (R[u][1] != 0)]
Cset = set(C)
core = sorted(set(star[0]) | set(star[1]))                 # 22 labels with complete rows
coreset = set(core)

m = cp_model.CpModel()
pos, neg = {}, {}


def pair(u, v):
    a, b = min(u, v), max(u, v)
    if (a, b) not in pos:
        pos[a, b] = m.NewBoolVar(''); neg[a, b] = m.NewBoolVar('')
        m.AddBoolOr([pos[a, b].Not(), neg[a, b].Not()])
    return a, b


def NE(u, v):
    a, b = pair(u, v)
    return pos[a, b] - neg[a, b]


# variables: every pair touching the core
for u in core:
    for v in range(L):
        if v != u:
            pair(u, v)


def star_structure(z):
    S = star[z]
    Pm = {e: m.NewBoolVar('') for e in itertools.combinations(S, 2)}
    Qm = {e: m.NewBoolVar('') for e in itertools.combinations(S, 2)}
    for u in S:
        m.Add(sum(Pm[e] for e in Pm if u in e) == 1)
        m.Add(sum(Qm[e] for e in Qm if u in e) == 1)
    for (u, v) in Pm:
        m.Add(R[u][z] * R[v][z] * NE(u, v) == Qm[u, v] - Pm[u, v])


for z in (0, 1):
    star_structure(z)
    for w in range(L):
        if R[w][z] != 0:
            continue
        A = {t: m.NewBoolVar('') for t in star[z]}; B = {t: m.NewBoolVar('') for t in star[z]}
        m.Add(sum(A.values()) == 2); m.Add(sum(B.values()) == 2)
        for t in star[z]:
            m.Add(R[t][z] * NE(w, t) == B[t] - A[t])
if starx:
    for z in range(2, 7):
        star_structure(z)
if s2all:                                                   # outside rows into the five outer stars as well
    for z in range(2, 7):
        for w in range(L):
            if R[w][z] != 0:
                continue
            A = {t: m.NewBoolVar('') for t in star[z]}; B = {t: m.NewBoolVar('') for t in star[z]}
            m.Add(sum(A.values()) == 2); m.Add(sum(B.values()) == 2)
            for t in star[z]:
                m.Add(R[t][z] * NE(w, t) == B[t] - A[t])
if n2core:                                                  # N^2 = N + 12I - G_r on every pair of core labels (full rows known)
    Nv = {}
    def NV(u, v):
        a, b = pair(u, v)
        if (a, b) not in Nv:
            x = m.NewIntVar(-1, 1, ''); m.Add(x == pos[a, b] - neg[a, b]); Nv[a, b] = x
        return Nv[a, b]
    for i, u in enumerate(core):
        for v in core[i + 1:]:
            terms = []
            for w in range(L):
                if w == u or w == v:
                    continue
                t = m.NewIntVar(-1, 1, ''); m.AddMultiplicationEquality(t, [NV(u, w), NV(w, v)]); terms.append(t)
            rr = sum(R[u][q] * R[v][q] for q in range(7))
            m.Add(sum(terms) == NE(u, v) - rr)
for u in core:
    a_, b_ = None, None
    m.Add(sum(pos[min(u, v), max(u, v)] + neg[min(u, v), max(u, v)] for v in range(L) if v != u) == 10)
sgn = {u: m.NewBoolVar(f's{u}') for u in C}


def prodNa(u, v):
    a, b = pair(u, v)
    pp = m.NewBoolVar(''); mm = m.NewBoolVar('')
    m.AddMultiplicationEquality(pp, [pos[a, b], sgn[v]])
    m.AddMultiplicationEquality(mm, [neg[a, b], sgn[v]])
    return 2 * pp - 2 * mm - pos[a, b] + neg[a, b]         # N_uv * (2 s_v - 1)


for u in range(L):
    terms = [prodNa(u, v) for v in C if v != u]
    m.Add(sum(terms) == (4 * (2 * sgn[u] - 1) if u in Cset else 0))
for z in range(7):
    m.Add(sum(R[u][z] * (2 * sgn[u] - 1) for u in C) == 0)
m.Add(sgn[C[0]] == 1)
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = T
solver.parameters.num_workers = Wk
t0 = time.time()
st = solver.Solve(m)
print(f'n5_local starx={starx} s2all={s2all} n2core={n2core}: {solver.StatusName(st)} ({time.time() - t0:.1f}s), vars {len(pos)} pairs', flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    sol = {'a': {u: 2 * solver.Value(sgn[u]) - 1 for u in C},
           'N': {f'{u},{v}': solver.Value(pos[u, v]) - solver.Value(neg[u, v]) for (u, v) in pos}}
    json.dump(sol, open(f'n5_local_solution{"_x" if starx else ""}{"_s2" if s2all else ""}{"_n2" if n2core else ""}.json', 'w'))
    print('solution saved')
