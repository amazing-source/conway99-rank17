"""Rank-17 characteristic scan (Astra 29/09, re-proved by C2). For the switched frame: sum_x x x^T + 2 y1 y1^T + 2 y2 y2^T = 14 I on V, the
42 residual x have odd norm, so for every characteristic chi of V: (chi,y1)^2 + (chi,y2)^2 <= 7|chi|^2 - 21.
Completeness of the shell: a violating class of Char(K*) (norm < 3) lifts, after subtracting even multiples of y1, y2, to chi with odd
heights |h| <= 7, hence |chi|^2 < 3 + 98/7 = 17, i.e. |chi|^2 in {1, 9}. So the characteristic vectors of norm <= 9 decide the test exactly.
All arithmetic below is in integer coefficients (exact). Usage: python r17_char_scan.py"""
import sys, os, json, collections, numpy as np
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import frame_dfs as FD
C1 = r"C:\Users\RESEARCHER\AppData\Local\Temp\claude\C--Users-RESEARCHER-Desktop\19a7c411-3446-4f59-9793-6c42892ed6cd\scratchpad\ext\conway99-c1"
sys.path.insert(0, os.path.join(C1, "code")); import rank17_cores as RC
opened = collections.defaultdict(set)
for x in json.load(open(os.path.join(C1, "data", "rank17", "open_after_sol.json"))): v, i = x.split(":"); opened[v].add(int(i))
cores = RC.cores(); out = {}; CERT = []
for Vn in ("E8", "D12p", "E72p", "A15p", "D82p", "A11E6p"):
    d = json.load(open(os.path.join(C1, "data", "rank17", "ypairs_%s.json" % Vn))); m = d["units"]; gens, dim = cores[d["core"]]
    Bc = [[F(x) for x in row] for row in RC.int_basis(gens)]                 # exact core basis rows (ambient dim)
    Bx = [row + [F(0)] * m for row in Bc] + [[F(0)] * dim + [F(int(i == j)) for j in range(m)] for i in range(m)]
    n = len(Bx); assert n == 17
    Gq = [[sum(a * b for a, b in zip(Bx[i], Bx[j])) for j in range(n)] for i in range(n)]
    assert all(x.denominator == 1 for r in Gq for x in r); G = np.array([[int(x) for x in r] for r in Gq], dtype=np.int64)
    c = FD.solve_mod2(G, np.diag(G)); Z = FD.short(G, 3, c); Z = FD.one_per_pair(Z)
    nz = np.einsum("ij,jk,ik->i", Z, G, Z); Z = Z[(nz == 1) | (nz == 9)]; nz = nz[(nz == 1) | (nz == 9)]
    assert (((Z @ G) - np.diag(G)[None, :]) % 2 == 0).all()                  # characteristic, exactly
    Bf = np.array([[float(x) for x in r] for r in Bx]); Binv = np.linalg.pinv(Bf)
    def coef(amb):
        a = [F(x) for x in amb]; cf = np.rint(np.array([float(x) for x in a]) @ Binv).astype(np.int64)
        assert all(sum(int(cf[i]) * Bx[i][j] for i in range(n)) == a[j] for j in range(len(a))), "y not in V"
        return cf
    kill = 0; kill_single = 0; rows = []
    for p in d["pairs"]:
        if p["id"] not in opened[Vn]: continue
        y1 = coef(p["y1_core"] + [str(u) for u in p["y1_units"]]); y2 = coef(p["y2_core"] + [str(u) for u in p["y2_units"]])
        assert y1 @ G @ y1 == 7 and y2 @ G @ y2 == 7 and y1 @ G @ y2 == 0
        h1 = Z @ G @ y1; h2 = Z @ G @ y2; bound = 7 * nz - 21
        bad = np.flatnonzero(h1 ** 2 + h2 ** 2 > bound); single = np.flatnonzero((h1 ** 2 > bound - 1) | (h2 ** 2 > bound - 1))
        if len(bad): kill += 1; rows.append((p["id"], int(nz[bad[0]]), int(h1[bad[0]]), int(h2[bad[0]]))); CERT.append({"V": Vn, "id": p["id"], "chi_coef": [int(t) for t in Z[bad[0]]], "chi_norm": int(nz[bad[0]]), "heights": [int(h1[bad[0]]), int(h2[bad[0]])], "bound": int(bound[bad[0]]), "partner_independent": bool(len(single))})
        kill_single += len(single) > 0
    out[Vn] = (len(opened[Vn]), kill, kill_single)
    print("%-7s shell(norm<=9) %4d  open %3d  KILLED %3d  (partner-independent %3d)  e.g. %s" % (Vn, len(Z), len(opened[Vn]), kill, kill_single, rows[:3]), flush=True)
print("TOTAL open (non-Z17):", sum(v[0] for v in out.values()), " killed:", sum(v[1] for v in out.values()))
json.dump({"source": "C1 888e731 open_after_sol.json + ypairs_*.json; test: (chi,y1)^2+(chi,y2)^2 > 7|chi|^2-21, chi characteristic of V, |chi|^2 in {1,9} (complete)", "killed": CERT}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "rank17_char_kills_2909.json"), "w"), indent=0)
