"""Per-skeleton LP relaxation for V = Z^17 (structures, 28/09; C2's formulation).
Variables x_{X,t} in [0,1] for every domain entry t = (w+, w-) of every cell X (domains as stageB.Skeleton.build).
  (1) one pair per cell:        sum_t x_{X,t} = 1                                    (21 rows)
  (2) stars (g-form):           for each vertex v: sum over cells X at v of sum_t x_{X,t} * contrib_v(t) = 0   (7 x 17 rows)
                                 contrib_v = w+ + w- if v = c (X = (c,d)), else w+ - w-
  (3) tight frame:              sum_{X,t} x_{X,t} (w+ w+^T + w- w-^T) = 28 I - 4 (y1 y1^T + y2 y2^T)          (153 rows)
Infeasible => the skeleton is dead.  Exact certificate: integer Farkas vector lam with lam.b - sum_j max(0,(A^T lam)_j) > 0
(0 <= x <= 1), verified in Python integers (same format as routeC_check/node_lp.py and exact_helpers.box_farkas_gap).
Usage: python skel_lp.py PAIR HEXCOLS"""
import sys, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack, identity, vstack
import stageB as B


def build_lp(sk):
    z = sk.z
    cols_blocks = []
    rows_b = []
    iu = [(i, j) for i in range(17) for j in range(i, 17)]
    nvar = sum(len(sk.dom[X][0]) for X in B.CELLS)
    A_one = np.zeros((21, nvar), dtype=np.int64)
    A_star = np.zeros((7 * 17, nvar), dtype=np.int64)
    A_fr = np.zeros((len(iu), nvar), dtype=np.int64)
    off = 0
    I_, J_ = np.array([a for a, b in iu]), np.array([b for a, b in iu])
    for k, X in enumerate(B.CELLS):
        P, M = sk.dom[X]
        n = len(P)
        c, d = X
        A_one[k, off:off + n] = 1
        cp = P + M      # contribution at vertex c
        cm = P - M      # contribution at vertex d
        A_star[17 * c:17 * c + 17, off:off + n] = cp.T
        A_star[17 * d:17 * d + 17, off:off + n] = cm.T
        A_fr[:, off:off + n] = (P[:, I_] * P[:, J_] + M[:, I_] * M[:, J_]).T
        off += n
    F = z.F
    b_one = np.ones(21, dtype=np.int64)
    b_star = np.zeros(7 * 17, dtype=np.int64)
    b_fr = np.array([F[i, j] for i, j in iu], dtype=np.int64)
    A = np.vstack([A_one, A_star, A_fr])
    b = np.concatenate([b_one, b_star, b_fr])
    return A, b


def lp_verdict(A, b):
    m, n = A.shape
    Af = hstack([csr_matrix(A.astype(float)), identity(m, format='csr'), -identity(m, format='csr')]).tocsr()
    cst = np.concatenate([np.zeros(n), np.ones(2 * m)])
    res = linprog(cst, A_eq=Af, b_eq=b.astype(float), bounds=[(0, 1)] * n + [(0, None)] * (2 * m), method='highs-ds')
    if res.status != 0:
        return 'LPERR', res.status, None
    if res.fun < 1e-7:
        return 'FEASIBLE', res.fun, None
    lam = np.asarray(res.eqlin.marginals, dtype=float)
    for D in (1, 2, 4, 8, 16, 64, 256, 1024, 4096, 65536, 10 ** 6, 10 ** 8):
        li = np.rint(lam * D).astype(object)
        red = A.T.astype(object) @ li
        val = int(b.astype(object) @ li) - int(sum(v for v in red if v > 0))
        if val > 0:
            return 'INFEASIBLE_CERT', res.fun, (D, val)
    return 'INFEASIBLE_UNCERT', res.fun, None


if __name__ == '__main__':
    pid = int(sys.argv[1]); hx = sys.argv[2]
    z = B.Z17(pid)
    cols, pis = B.parse_raw(hx)
    sk = B.Skeleton(z, pis).build()
    t = time.time()
    A, b = build_lp(sk)
    v = lp_verdict(A, b)
    print('pair', pid, hx, 'vars', A.shape[1], 'rows', A.shape[0], 'LP:', v[0], 'phase1 %.3g' % v[1] if v[1] is not None else '', 'cert', v[2], '%.1fs' % (time.time() - t))
