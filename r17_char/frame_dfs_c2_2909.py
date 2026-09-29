"""Finish line (C): generic frame search in an envelope U with the pair (y1, y2) FIXED. This is C2's generalisation of o24_dfs.py
(structures' request, 27/09), for the pairs the frame LP cannot exclude and as an independent reader of C1's h-stage.

Lattice: U given by rational generators (ambient rows) and an inner-product scale s, so (u, v) = s u.v. Basis: HNF + LLL.
Coefficient vectors are integers and the Gram matrix G is integral and unimodular. Everything below is exact (integers), except
Fincke-Pohst, which enumerates in floats and re-checks every norm exactly.
- New roots: the norm-2 vectors of the shadow chi/2 + U, as 2 beta = 2x + c (c characteristic). One "point" per +- pair.
- Frame identity at a new root beta (C1 notes/140-141, all proven):
  - |(h, beta)| <= 1 and w(beta) = sum_c (h_c, beta)^2 is in {1, 3};
  - w(beta) = 1 exactly on the beta hit at 3/2 by an odd frame vector.
  - RIGID lattices (no norm-3 vector hits a new root at 3/2), e.g. O24 and 10, 44, 87, 115, 132: target w = 1 on the y-hits
    and 3 elsewhere. Only the rigid case is implemented; a non-rigid lattice stops with an error.
- Candidates (one per +- pair):
  - h: norm 4, orthogonal to y1 and y2, |(h, beta)| <= 1 for every new root;
  - u (words): norm 3, orthogonal to y1 and y2.
  - V(h, u) = (h, u) is in {0, +-1, +-2}; "Good" means no |V| = 2.
- Stage 1 (DFS over the h, points = new-root pairs, need = target), as in o24_dfs.py:
  - pairwise orthogonal;
  - (a) Same/Opp sets non-empty and orthogonal to the other chosen h (C1's CELL);
  - (b) star0 >= 2 (7 - |P|) and #common-zero words >= (7 - |P|)(6 - |P|);
  - no two h through a target-1 point; need <= left.
- Stage 2 at the leaves:
  - Lemma Z: sum h in 2U;
  - L7 local, for each c: sum over the 12 words of c (one Same and one Opp per partner), oriented so that (u, h_c) = +1, equals
    3 h_c (exact coefficients; pruning by the inner products with the new roots);
  - then, in Python: L7 joint (42 distinct words), T^2 + T = 12 I - C^T C, and sum over the 51 frame vectors of x x^T = 7 G^-1.
Output: {nodes, counts, survivors, verdict}.
Usage: python frame_dfs.py --o24 <class> (control against o24_dfs)  |  python frame_dfs.py <U.json> <pairs.json> <k> [<k> ...]
       python frame_dfs.py --batch <U.json> <pairs.json> <out.jsonl> [from] [to]   (resumable, one JSON line per pair)"""
import sys, os, json, time, itertools
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"): os.environ.setdefault(_v, "1")
from fractions import Fraction as F
import numpy as np
import numba as nb
LAB = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(LAB, "external_runs", "c1_96a94fe", "code"))
from envelope_trap import to_int_basis, lll                             # standard HNF + exact LLL (not the h-stage bitsets)
U1 = np.uint64(1)


# ---------------------------------------------------------------- lattice
def setup_lattice(gens, scale=F(1)):
    B = lll(to_int_basis([[F(x) for x in g] for g in gens]))
    assert len(B) == 24
    Gq = [[scale * sum(a * b for a, b in zip(u, v)) for v in B] for u in B]
    assert all(x.denominator == 1 for row in Gq for x in row), "not integral"
    G = np.array([[int(x) for x in row] for row in Gq], dtype=np.int64)
    assert round(abs(np.linalg.det(G.astype(float)))) == 1, "not unimodular"
    return B, G


def solve_mod2(G, rhs):
    n = len(G); A = [[int(G[i][j]) & 1 for j in range(n)] + [int(rhs[i]) & 1] for i in range(n)]; row = 0; piv = []
    for c in range(n):
        p = next((i for i in range(row, n) if A[i][c]), None)
        if p is None: continue
        A[row], A[p] = A[p], A[row]
        for i in range(n):
            if i != row and A[i][c]: A[i] = [a ^ b for a, b in zip(A[i], A[row])]
        piv.append(c); row += 1
    assert row == n
    x = [0] * n
    for i, c in enumerate(piv): x[c] = A[i][n]
    return np.array(x, dtype=np.int64)


@nb.njit(cache=True)
def fp_enum(G, Lu, d, N4, shift2, maxout):
    """all integer x with (2x + shift2)^T G (2x + shift2) <= N4 (exact check), where G = Lu^T diag(d) Lu (Lu unit upper
    triangular, floats). Returns the rows 2x + shift2."""
    n = G.shape[0]; out = np.zeros((maxout, n), np.int64); cnt = 0
    x = np.zeros(n, np.int64); z = np.zeros(n)                         # z = x + shift2 / 2
    rem = np.zeros(n + 1); cen = np.zeros(n); hi = np.zeros(n, np.int64)
    s = np.empty(n)
    for i in range(n): s[i] = shift2[i] / 2.0
    bound = N4 / 4.0 + 1e-7
    i = n - 1; rem[n] = bound
    c = 0.0
    for j in range(i + 1, n): c -= Lu[i, j] * z[j]
    cen[i] = c
    r = np.sqrt(max(rem[i + 1], 0.0) / d[i])
    x[i] = np.int64(np.ceil(cen[i] - s[i] - r)); hi[i] = np.int64(np.floor(cen[i] - s[i] + r))
    while True:
        if x[i] > hi[i]:
            i += 1
            if i >= n: break
            x[i] += 1
            continue
        z[i] = x[i] + s[i]
        t = z[i] - cen[i]; rr = rem[i + 1] - d[i] * t * t
        if rr < -1e-9:
            x[i] += 1
            continue
        if i == 0:
            v = np.empty(n, np.int64)
            for j in range(n): v[j] = 2 * x[j] + shift2[j]
            q = 0
            for a in range(n):
                if v[a] != 0:
                    for b in range(n): q += v[a] * G[a, b] * v[b]
            if q <= N4:
                if cnt < maxout:
                    for j in range(n): out[cnt, j] = v[j]
                cnt += 1
            x[i] += 1
            continue
        rem[i] = rr; i -= 1
        c = 0.0
        for j in range(i + 1, n): c -= Lu[i, j] * z[j]
        cen[i] = c
        r = np.sqrt(max(rem[i + 1], 0.0) / d[i])
        x[i] = np.int64(np.ceil(cen[i] - s[i] - r)); hi[i] = np.int64(np.floor(cen[i] - s[i] + r))
    return out[:min(cnt, maxout)], cnt


def ldl(G):
    n = len(G); Gf = G.astype(float); Lu = np.eye(n); d = np.zeros(n)
    Lc = np.linalg.cholesky(Gf)                                         # G = Lc Lc^T, Lc lower
    for i in range(n): d[i] = Lc[i, i] ** 2
    for i in range(n):
        for j in range(i + 1, n): Lu[i, j] = Lc[j, i] / Lc[i, i]
    return Lu, d


def short(G, N, shift2=None, maxout=3_000_000):
    """vectors (as 2x + shift2 rows) of norm <= N in U (shift2 = 0) or in the shadow (shift2 = c)"""
    Lu, d = ldl(G); sh = np.zeros(len(G), np.int64) if shift2 is None else np.asarray(shift2, np.int64)
    out, cnt = fp_enum(G, Lu, d, 4 * N, sh, maxout)
    assert cnt <= maxout, "maxout too small"
    return out


def one_per_pair(V):
    """keep the rows whose first non-zero entry is positive"""
    keep = []
    for v in V:
        nz = np.flatnonzero(v)
        if len(nz) and v[nz[0]] > 0: keep.append(v)
    return np.array(keep, dtype=np.int64)


def lattice_data(gens, scale=F(1)):
    t0 = time.time(); B, G = setup_lattice(gens, scale)
    V2 = short(G, 4)                                                    # rows = 2x (norm <= 4)
    X = one_per_pair(V2 // 2); nrm = np.einsum("ij,jk,ik->i", X, G, X)
    assert nrm.min() >= 2, "norm-1 vectors"
    R, X3, X4 = X[nrm == 2], X[nrm == 3], X[nrm == 4]
    c = solve_mod2(G, np.diag(G))
    S2 = one_per_pair(short(G, 2, c))                                    # rows = 2 beta, one per +- pair
    assert all(int(v @ G @ v) == 8 for v in S2)
    r = 2 * len(R)
    assert 2 * len(X3) == 4096, "norm-3 count %d" % (2 * len(X3))
    assert 2 * len(S2) == r + 48, "shadow roots %d, r = %d" % (2 * len(S2), r)
    hits3 = np.abs(X3 @ G @ S2.T) == 3                                 # 2 (x, beta) = +-3
    return {"B": B, "G": G, "R": R, "X3": X3, "X4": X4, "S2": S2, "r": r, "rigid": not hits3.any(), "secs": round(time.time() - t0, 1)}


# ---------------------------------------------------------------- one pair
def pair_setup(L, y1, y2):
    G, X3, X4, S2 = L["G"], L["X3"], L["X4"], L["S2"]
    for y in (y1, y2): assert int(y @ G @ y) == 7
    assert int(y1 @ G @ y2) == 0
    Y = np.array([y1, y2]); YS = Y @ G @ S2.T                             # 2 (y, beta)
    if np.abs(YS).max() > 3: return None, "a y has |(y, beta)| > 3/2 at a new root (impossible in a frame)"
    yhit = (np.abs(YS) == 3).any(axis=0)
    if ((np.abs(YS) == 3).sum(axis=0) > 1).any(): return None, "both y hit the same new root (impossible)"
    H = X4[((X4 @ G @ Y.T) == 0).all(axis=1)]
    HS = H @ G @ S2.T                                                  # 2 (h, beta) in {0, +-2}
    H = H[np.abs(HS).max(axis=1) <= 2]; HS = H @ G @ S2.T
    W = X3[((X3 @ G @ Y.T) == 0).all(axis=1)]
    WS = W @ G @ S2.T                                                  # 2 (w, beta): +-1, or +-3 at a 3/2-hit
    V = H @ G @ W.T                                                    # (h, w)
    assert V.size == 0 or np.abs(V).max() <= 2                        # same check as set(np.unique(V)) <= {0,+-1,+-2}, cheaper
    whit = np.abs(WS) == 3                                             # 3/2-hits of the candidate u
    # GENERAL mode (C1 notes/142): a new root met once by the h is hit at 3/2 by exactly one odd frame vector (a y or a u), and a u
    # is a norm-3 vector orthogonal to y. So w = 1 on the y-hits, w in {1,3} on the other roots hit by some candidate u
    # ("flex"), w = 3 elsewhere. Rigid lattices have no flex point.
    flex = (~yhit) & whit.any(axis=0)
    target = np.where(yhit, 1, np.where(flex, 0, 3)).astype(np.int64)   # 0 marks a flex point
    return {"H": H, "HS": HS, "HSfull": HS, "W": W, "WS": WS, "V": V, "target": target, "yhit": yhit, "flex": flex, "whit": whit, "Y": Y}, None


def pair_setup_fast(L, y1, y2):
    """Same output as pair_setup (C2 29/09): lattice-level caches (X4 G S2^T for all norm-4 lines, admissibility, X3 G S2^T) and
    float64 BLAS products (exact: every entry is a small integer), instead of per-pair int64 matmuls."""
    G, X3, X4, S2 = L["G"], L["X3"], L["X4"], L["S2"]
    if "_X4S" not in L:
        Gf = G.astype(np.float64); L["_X4Gf"] = X4.astype(np.float64) @ Gf; L["_X3Gf"] = X3.astype(np.float64) @ Gf
        L["_X4S"] = np.rint(L["_X4Gf"] @ S2.T.astype(np.float64)).astype(np.int64)
        L["_X4adm"] = np.abs(L["_X4S"]).max(axis=1) <= 2
        L["_X3S"] = np.rint(L["_X3Gf"] @ S2.T.astype(np.float64)).astype(np.int64)
    for y in (y1, y2): assert int(y @ G @ y) == 7
    assert int(y1 @ G @ y2) == 0
    Y = np.array([y1, y2]); YS = Y @ G @ S2.T
    if np.abs(YS).max() > 3: return None, "a y has |(y, beta)| > 3/2 at a new root (impossible in a frame)"
    yhit = (np.abs(YS) == 3).any(axis=0)
    if ((np.abs(YS) == 3).sum(axis=0) > 1).any(): return None, "both y hit the same new root (impossible)"
    Yf = Y.T.astype(np.float64)
    ih = np.flatnonzero((np.abs(L["_X4Gf"] @ Yf) < 0.5).all(axis=1) & L["_X4adm"])
    H = X4[ih]; HS = L["_X4S"][ih]
    iw = np.flatnonzero((np.abs(L["_X3Gf"] @ Yf) < 0.5).all(axis=1))
    W = X3[iw]; WS = L["_X3S"][iw]
    V = np.rint(L["_X4Gf"][ih] @ W.T.astype(np.float64)).astype(np.int64)
    assert V.size == 0 or np.abs(V).max() <= 2
    whit = np.abs(WS) == 3
    flex = (~yhit) & whit.any(axis=0)
    target = np.where(yhit, 1, np.where(flex, 0, 3)).astype(np.int64)
    return {"H": H, "HS": HS, "HSfull": HS, "W": W, "WS": WS, "V": V, "target": target, "yhit": yhit, "flex": flex, "whit": whit, "Y": Y}, None


def pack64(Bm):
    Bm = np.asarray(Bm, dtype=bool); nw = (Bm.shape[1] + 63) // 64
    P8 = np.packbits(Bm, axis=1, bitorder="little")
    P8 = np.concatenate([P8, np.zeros((Bm.shape[0], nw * 8 - P8.shape[1]), dtype=np.uint8)], axis=1)
    return np.ascontiguousarray(P8).view("<u8").astype(np.uint64).reshape(Bm.shape[0], nw)


@nb.njit(cache=True, inline="always")
def popc(x):
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return np.int64((x * np.uint64(0x0101010101010101)) >> np.uint64(56))


@nb.njit(cache=True)
def leaf_l7(P, Hc, Wc, V, WS, HS, Z, Pp, Pn, allW, slots, slen, part, pos, npt):
    """0: Lemma Z fails; 1: some local star sum fails; 2: all 7 pass. Hc, Wc: coefficient vectors; WS, HS: 2 (., beta)."""
    nWw = Z.shape[1]; n = Hc.shape[1]
    for i in range(n):
        s = 0
        for x in range(7): s += Hc[P[x], i]
        if s % 2 != 0: return 0
    for x in range(7):
        c = P[x]; ns = 0
        for y in range(7):
            if y == x: continue
            d = P[y]
            for kd in range(2):
                L = 0
                for w in range(nWw):
                    zo = allW[w]
                    for f in range(7):
                        if f != x and f != y: zo &= Z[P[f], w]
                    if kd == 0: s_ = ((Pp[c, w] & Pp[d, w]) | (Pn[c, w] & Pn[d, w])) & zo
                    else: s_ = ((Pp[c, w] & Pn[d, w]) | (Pn[c, w] & Pp[d, w])) & zo
                    while s_ != 0:
                        low = s_ & (~s_ + U1); slots[ns, L] = w * 64 + popc(low - U1); L += 1; s_ ^= low
                slen[ns] = L; ns += 1
        order = np.argsort(slen[:12], kind="mergesort")
        for i in range(npt): part[0, i] = 0
        lvl = 0; entering = True; found = False
        while True:
            if entering:
                r = 12 - lvl; okb = True
                for q in range(npt):                                     # target 3 * 2 (h_c, beta_q); each word adds +-1
                    dl = 3 * HS[c, q] - part[lvl, q]
                    if dl > r or dl < -r or (dl - r) % 2 != 0:
                        okb = False; break
                if okb and lvl == 12:
                    # exact coefficients: sum of the 12 oriented words = 3 h_c
                    ok2 = True
                    for i in range(n):
                        s = 0
                        for l2 in range(12):
                            t = slots[order[l2], pos[l2]]; sg = 1 if V[c, t] > 0 else -1
                            s += sg * Wc[t, i]
                        if s != 3 * Hc[c, i]:
                            ok2 = False; break
                    if ok2:
                        found = True; break
                    okb = False
                if not okb:
                    lvl -= 1
                    if lvl < 0: break
                    entering = False; continue
                pos[lvl] = -1; entering = False
            sl = order[lvl]; pos[lvl] += 1
            if pos[lvl] >= slen[sl]:
                lvl -= 1
                if lvl < 0: break
                continue
            t = slots[sl, pos[lvl]]; sg = 1 if V[c, t] > 0 else -1
            for q in range(npt): part[lvl + 1, q] = part[lvl, q] + sg * WS[t, q]
            lvl += 1; entering = True
        if not found: return 1
    return 2


@nb.njit(cache=True)
def dfs_generic(compat, through, supp, suppn, Z, Pp, Pn, Bd, NZ, target, allM, allW, node_limit, cap, Hc, Wc, V, WS, HS, p0, fmask):
    """o24_dfs kernel 2 on general points (no support tables: only need <= left), with leaf_l7 at the leaves"""
    nMw = compat.shape[1]; nWw = Z.shape[1]; npt = target.shape[0]; nW = Wc.shape[0]
    cand = np.zeros((8, nMw), np.uint64); Cs = np.zeros((8, nMw), np.uint64)
    nzw = np.zeros((8, nMw), np.int64); nnz = np.zeros(8, np.int64)
    need = np.zeros((8, npt), np.int64)
    ge1 = np.zeros((8, nWw), np.uint64); ge2 = np.zeros((8, nWw), np.uint64); good = np.zeros((8, nWw), np.uint64)
    P = np.zeros(8, np.int64); curj = np.zeros(8, np.int64); nodes = np.zeros(8, np.int64)
    leaves = np.zeros((cap, 7), np.int64); nleaves = 0; overflow = False; counts = np.zeros(3, np.int64)
    slots = np.zeros((12, nW), np.int64); slen = np.zeros(12, np.int64); part = np.zeros((13, npt), np.int64); pos = np.zeros(13, np.int64)
    for w in range(nMw):
        cand[0, w] = allM[w]; nzw[0, w] = w
    nnz[0] = nMw
    for w in range(nWw): good[0, w] = allW[w]
    for i in range(npt): need[0, i] = target[i]
    total = 0; stop = False; d = 0; action = 0
    while True:
        if action == 0:
            nodes[d] += 1; total += 1
            if node_limit > 0 and total >= node_limit: stop = True
            if stop:
                action = 2; continue
            if d == 7:
                counts[0] += 1
                res = leaf_l7(P, Hc, Wc, V, WS, HS, Z, Pp, Pn, allW, slots, slen, part, pos, npt)
                if res >= 1: counts[1] += 1
                if res == 2:
                    counts[2] += 1
                    if nleaves < cap:
                        for x in range(7): leaves[nleaves, x] = P[x]
                        nleaves += 1
                    else: overflow = True
                action = 2; continue
            left = 7 - d; best = -1; bkey = 1 << 30; pruned = False
            for i in range(npt):
                r = need[d, i]
                if r > left:
                    pruned = True; break
                if r > 0:
                    n = 0
                    for t in range(nnz[d]):
                        w = nzw[d, t]; n += popc(cand[d, w] & through[i, w])
                    if n < r:
                        pruned = True; break
                    if n - r < bkey:
                        best = i; bkey = n - r
            if pruned:
                action = 2; continue
            for t in range(nnz[d]):
                w = nzw[d, t]; Cs[d, w] = cand[d, w] & through[best, w]
            curj[d] = 0
            action = 1; continue
        elif action == 1:
            j = -1
            while curj[d] < nnz[d]:
                w = nzw[d, curj[d]]; x = Cs[d, w]
                if x != 0:
                    low = x & (~x + U1)
                    j = w * 64 + popc(low - U1); Cs[d, w] = x ^ low
                    break
                curj[d] += 1
            if j < 0:
                action = 2; continue
            P[d] = j
            e1 = d + 1
            for w in range(nWw):
                good[e1, w] = good[d, w] & ~Bd[j, w]
                ge1[e1, w] = ge1[d, w] | NZ[j, w]
                ge2[e1, w] = ge2[d, w] | (ge1[d, w] & NZ[j, w])
            left = 7 - e1; ok = True
            for x in range(e1):
                if not ok: break
                for y in range(x + 1, e1):
                    c = P[x]; e = P[y]; same = np.uint64(0); opp = np.uint64(0)
                    for w in range(nWw):
                        zo = good[e1, w]
                        for f in range(e1):
                            if f != x and f != y: zo &= Z[P[f], w]
                        same |= ((Pp[c, w] & Pp[e, w]) | (Pn[c, w] & Pn[e, w])) & zo
                        opp |= ((Pp[c, w] & Pn[e, w]) | (Pn[c, w] & Pp[e, w])) & zo
                    if same == 0 or opp == 0:
                        ok = False; break
            if ok and left > 0:
                n = 0
                for w in range(nWw): n += popc(good[e1, w] & ~ge1[e1, w])
                if n < left * (left - 1): ok = False
                else:
                    for x in range(e1):
                        c = P[x]; n = 0
                        for w in range(nWw): n += popc(NZ[c, w] & good[e1, w] & ~ge2[e1, w])
                        if n < 2 * left:
                            ok = False; break
            if ok:
                for i in range(npt): need[e1, i] = need[d, i]
                cnt = 0
                for t in range(nnz[d]):
                    w = nzw[d, t]; v = cand[d, w] & compat[j, w]
                    if v != 0:
                        cand[e1, w] = v; nzw[e1, cnt] = w; cnt += 1
                nnz[e1] = cnt
                for q in range(suppn[j]):
                    x = supp[j, q]; need[e1, x] -= 1
                    if need[e1, x] == 0:
                        for t in range(cnt):
                            w = nzw[e1, t]; cand[e1, w] &= ~through[x, w]
                d = e1; action = 0; continue
            cand[d, j >> 6] &= ~(U1 << np.uint64(j & 63))
            action = 1; continue
        else:
            if d == 0: break
            d -= 1; j = P[d]
            cand[d, j >> 6] &= ~(U1 << np.uint64(j & 63))
            action = 2 if stop else 1
    return nodes, leaves[:nleaves], not stop, overflow, counts


@nb.njit(cache=True)
def leaf_general(P, cov, isflex, whitb, Hc, Wc, V, WS, HS, Z, Pp, Pn, allW, slots, slen, part, pos, npt, allowed):
    """GENERAL mode leaf: F1 = flex points met once. A candidate u is allowed iff all its 3/2-hits lie in F1 (elsewhere w = 3, or the
    point is a y-hit). Then Lemma Z and the local star sums with the allowed words (a word adds +-1, or +-3 at its own hits, to
    2 (., beta); at most one word of the frame hits a given point). Returns 0 / 1 / 2 as leaf_l7, and fills `allowed`."""
    nWw = Z.shape[1]; n = Hc.shape[1]; nW = Wc.shape[0]; nPw = whitb.shape[1]
    F1 = np.zeros(nPw, np.uint64)
    for p in range(npt):
        if isflex[p] and cov[p] == 1: F1[p >> 6] |= U1 << np.uint64(p & 63)
    for w in range(nWw): allowed[w] = np.uint64(0)
    for t in range(nW):
        ok = True
        for q in range(nPw):
            if whitb[t, q] & ~F1[q]:
                ok = False; break
        if ok: allowed[t >> 6] |= U1 << np.uint64(t & 63)
    for w in range(nWw): allowed[w] &= allW[w]
    for i in range(n):
        s = 0
        for x in range(7): s += Hc[P[x], i]
        if s % 2 != 0: return 0
    for x in range(7):
        c = P[x]; ns = 0
        for y in range(7):
            if y == x: continue
            d = P[y]
            for kd in range(2):
                L = 0
                for w in range(nWw):
                    zo = allowed[w]
                    for f in range(7):
                        if f != x and f != y: zo &= Z[P[f], w]
                    if kd == 0: s_ = ((Pp[c, w] & Pp[d, w]) | (Pn[c, w] & Pn[d, w])) & zo
                    else: s_ = ((Pp[c, w] & Pn[d, w]) | (Pn[c, w] & Pp[d, w])) & zo
                    while s_ != 0:
                        low = s_ & (~s_ + U1); slots[ns, L] = w * 64 + popc(low - U1); L += 1; s_ ^= low
                slen[ns] = L; ns += 1
        order = np.argsort(slen[:12], kind="mergesort")
        for i in range(npt): part[0, i] = 0
        lvl = 0; entering = True; found = False
        while True:
            if entering:
                r = 12 - lvl; okb = True
                for q in range(npt):
                    dl = 3 * HS[c, q] - part[lvl, q]
                    fl = 0
                    if (F1[q >> 6] >> np.uint64(q & 63)) & U1: fl = 2
                    if dl > r + fl or dl < -r - fl or (dl - r) % 2 != 0:
                        okb = False; break
                if okb and lvl == 12:
                    ok2 = True
                    for i in range(n):
                        s = 0
                        for l2 in range(12):
                            t = slots[order[l2], pos[l2]]; sg = 1 if V[c, t] > 0 else -1
                            s += sg * Wc[t, i]
                        if s != 3 * Hc[c, i]:
                            ok2 = False; break
                    if ok2:
                        found = True; break
                    okb = False
                if not okb:
                    lvl -= 1
                    if lvl < 0: break
                    entering = False; continue
                pos[lvl] = -1; entering = False
            sl = order[lvl]; pos[lvl] += 1
            if pos[lvl] >= slen[sl]:
                lvl -= 1
                if lvl < 0: break
                continue
            t = slots[sl, pos[lvl]]; sg = 1 if V[c, t] > 0 else -1
            for q in range(npt): part[lvl + 1, q] = part[lvl, q] + sg * WS[t, q]
            lvl += 1; entering = True
        if not found: return 1
    return 2


@nb.njit(cache=True)
def dfs_general(compat, through, supp, suppn, Z, Pp, Pn, Bd, NZ, target, isflex, whitb, allM, allW, node_limit, cap, Hc, Wc, V, WS, HS, p0, fmask):
    """GENERAL mode (flex points): target[p] = 1 or 3 exactly, or 0 for a flex point (final coverage 1 or 3). Pending points:
    exact points with need > 0, flex points with coverage 0 or 2. Branch on a pending point (smallest-index member through it);
    with no pending point, on the smallest-index remaining member (all candidates). Saturated: exact need 0, flex coverage 3."""
    nMw = compat.shape[1]; nWw = Z.shape[1]; npt = target.shape[0]; nW = Wc.shape[0]
    cand = np.zeros((8, nMw), np.uint64); Cs = np.zeros((8, nMw), np.uint64)
    nzw = np.zeros((8, nMw), np.int64); nnz = np.zeros(8, np.int64)
    need = np.zeros((8, npt), np.int64); cov = np.zeros((8, npt), np.int64)
    ge1 = np.zeros((8, nWw), np.uint64); ge2 = np.zeros((8, nWw), np.uint64); good = np.zeros((8, nWw), np.uint64)
    P = np.zeros(8, np.int64); curj = np.zeros(8, np.int64); nodes = np.zeros(8, np.int64)
    leaves = np.zeros((cap, 7), np.int64); leafcov = np.zeros((cap, npt), np.int64); nleaves = 0; overflow = False
    counts = np.zeros(3, np.int64)
    slots = np.zeros((12, nW), np.int64); slen = np.zeros(12, np.int64); part = np.zeros((13, npt), np.int64); pos = np.zeros(13, np.int64)
    allowed = np.zeros(nWw, np.uint64)
    for w in range(nMw):
        cand[0, w] = allM[w]; nzw[0, w] = w
    nnz[0] = nMw
    for w in range(nWw): good[0, w] = allW[w]
    for i in range(npt): need[0, i] = target[i]
    total = 0; stop = False; d = 0; action = 0
    while True:
        if action == 0:
            nodes[d] += 1; total += 1
            if node_limit > 0 and total >= node_limit: stop = True
            if stop:
                action = 2; continue
            if d == 7:
                valid = True
                for i in range(npt):
                    if isflex[i]:
                        if cov[d, i] != 1 and cov[d, i] != 3:
                            valid = False; break
                    elif need[d, i] != 0:
                        valid = False; break
                if valid:
                    counts[0] += 1
                    res = leaf_general(P, cov[d], isflex, whitb, Hc, Wc, V, WS, HS, Z, Pp, Pn, allW, slots, slen, part, pos, npt, allowed)
                    if res >= 1: counts[1] += 1
                    if res == 2:
                        counts[2] += 1
                        if nleaves < cap:
                            for x in range(7): leaves[nleaves, x] = P[x]
                            for i in range(npt): leafcov[nleaves, i] = cov[d, i]
                            nleaves += 1
                        else: overflow = True
                action = 2; continue
            left = 7 - d; best = -1; bkey = 1 << 30; pruned = False
            for i in range(npt):
                if isflex[i]:
                    c = cov[d, i]
                    r = 1 if (c == 0 or c == 2) else 0
                else:
                    r = need[d, i]
                if r > left:
                    pruned = True; break
                if r > 0:
                    n = 0
                    for t in range(nnz[d]):
                        w = nzw[d, t]; n += popc(cand[d, w] & through[i, w])
                    if n < r:
                        pruned = True; break
                    if n - r < bkey:
                        best = i; bkey = n - r
            if pruned:
                action = 2; continue
            if d == 0 and p0 >= 0: best = p0                             # C2 28/09: symmetry breaking (forced first point)
            if best >= 0:
                for t in range(nnz[d]):
                    w = nzw[d, t]; Cs[d, w] = cand[d, w] & through[best, w]
                    if d == 0 and p0 >= 0: Cs[d, w] &= fmask[w]
            else:                                                        # no pending point: any remaining member
                for t in range(nnz[d]):
                    w = nzw[d, t]; Cs[d, w] = cand[d, w]
            curj[d] = 0
            action = 1; continue
        elif action == 1:
            j = -1
            while curj[d] < nnz[d]:
                w = nzw[d, curj[d]]; x = Cs[d, w]
                if x != 0:
                    low = x & (~x + U1)
                    j = w * 64 + popc(low - U1); Cs[d, w] = x ^ low
                    break
                curj[d] += 1
            if j < 0:
                action = 2; continue
            P[d] = j
            e1 = d + 1
            for w in range(nWw):
                good[e1, w] = good[d, w] & ~Bd[j, w]
                ge1[e1, w] = ge1[d, w] | NZ[j, w]
                ge2[e1, w] = ge2[d, w] | (ge1[d, w] & NZ[j, w])
            left = 7 - e1; ok = True
            for x in range(e1):
                if not ok: break
                for y in range(x + 1, e1):
                    c = P[x]; e = P[y]; same = np.uint64(0); opp = np.uint64(0)
                    for w in range(nWw):
                        zo = good[e1, w]
                        for f in range(e1):
                            if f != x and f != y: zo &= Z[P[f], w]
                        same |= ((Pp[c, w] & Pp[e, w]) | (Pn[c, w] & Pn[e, w])) & zo
                        opp |= ((Pp[c, w] & Pn[e, w]) | (Pn[c, w] & Pp[e, w])) & zo
                    if same == 0 or opp == 0:
                        ok = False; break
            if ok and left > 0:
                n = 0
                for w in range(nWw): n += popc(good[e1, w] & ~ge1[e1, w])
                if n < left * (left - 1): ok = False
                else:
                    for x in range(e1):
                        c = P[x]; n = 0
                        for w in range(nWw): n += popc(NZ[c, w] & good[e1, w] & ~ge2[e1, w])
                        if n < 2 * left:
                            ok = False; break
            if ok:
                for i in range(npt):
                    need[e1, i] = need[d, i]; cov[e1, i] = cov[d, i]
                cnt = 0
                for t in range(nnz[d]):
                    w = nzw[d, t]; v = cand[d, w] & compat[j, w]
                    if v != 0:
                        cand[e1, w] = v; nzw[e1, cnt] = w; cnt += 1
                nnz[e1] = cnt
                for q in range(suppn[j]):
                    x = supp[j, q]; cov[e1, x] += 1
                    sat = False
                    if isflex[x]:
                        sat = cov[e1, x] == 3
                    else:
                        need[e1, x] -= 1; sat = need[e1, x] == 0
                    if sat:
                        for t in range(cnt):
                            w = nzw[e1, t]; cand[e1, w] &= ~through[x, w]
                d = e1; action = 0; continue
            cand[d, j >> 6] &= ~(U1 << np.uint64(j & 63))
            action = 1; continue
        else:
            if d == 0: break
            d -= 1; j = P[d]
            cand[d, j >> 6] &= ~(U1 << np.uint64(j & 63))
            action = 2 if stop else 1
    return nodes, leaves[:nleaves], leafcov[:nleaves], not stop, overflow, counts


@nb.njit(cache=True)
def dfs_general_u(compat, through, supp, suppn, Z, Pp, Pn, Bd, NZ, target, isflex, whitb, allM, allW, node_limit, cap, Hc, Wc, V, WS, HS, hitb, p0, fmask):
    """GENERAL mode + C1 LEMMA U+ (notes/145): a frame word hits (at 3/2) ONLY new roots that end at w = 1 and are not y-hits, and has the
    D7 pattern on the 7 h (exactly two (u,h) = +-1). At a node a word is VALID iff: no |(u,h)| = 2 and <= 2 nonzero on the chosen h,
    no 3/2-hit on a y-hit point, none on a flex point already covered >= 2 (monotone). The cells (Same/Opp) and the counting checks
    use the valid words only; a flex point at w <= 2 without a valid hitter must reach 3; at the leaf a w = 1 flex point needs a
    hitter with exactly 2 nonzero.
    GENERAL mode (flex points): target[p] = 1 or 3 exactly, or 0 for a flex point (final coverage 1 or 3). Pending points:
    exact points with need > 0, flex points with coverage 0 or 2. Branch on a pending point (smallest-index member through it);
    with no pending point, on the smallest-index remaining member (all candidates). Saturated: exact need 0, flex coverage 3."""
    nMw = compat.shape[1]; nWw = Z.shape[1]; npt = target.shape[0]; nW = Wc.shape[0]
    cand = np.zeros((8, nMw), np.uint64); Cs = np.zeros((8, nMw), np.uint64)
    nzw = np.zeros((8, nMw), np.int64); nnz = np.zeros(8, np.int64)
    need = np.zeros((8, npt), np.int64); cov = np.zeros((8, npt), np.int64)
    ge1 = np.zeros((8, nWw), np.uint64); ge2 = np.zeros((8, nWw), np.uint64); good = np.zeros((8, nWw), np.uint64)
    ge3 = np.zeros((8, nWw), np.uint64); bad = np.zeros((8, nWw), np.uint64); vm = np.zeros(nWw, np.uint64)
    P = np.zeros(8, np.int64); curj = np.zeros(8, np.int64); nodes = np.zeros(8, np.int64)
    leaves = np.zeros((cap, 7), np.int64); leafcov = np.zeros((cap, npt), np.int64); nleaves = 0; overflow = False
    counts = np.zeros(3, np.int64)
    slots = np.zeros((12, nW), np.int64); slen = np.zeros(12, np.int64); part = np.zeros((13, npt), np.int64); pos = np.zeros(13, np.int64)
    allowed = np.zeros(nWw, np.uint64)
    for w in range(nMw):
        cand[0, w] = allM[w]; nzw[0, w] = w
    nnz[0] = nMw
    for w in range(nWw): good[0, w] = allW[w]
    for i in range(npt):
        need[0, i] = target[i]
        if target[i] == 1:
            for w in range(nWw): bad[0, w] |= hitb[i, w]
    total = 0; stop = False; d = 0; action = 0
    while True:
        if action == 0:
            nodes[d] += 1; total += 1
            if node_limit > 0 and total >= node_limit: stop = True
            if stop:
                action = 2; continue
            if d == 7:
                valid = True
                for i in range(npt):
                    if isflex[i]:
                        if cov[d, i] != 1 and cov[d, i] != 3:
                            valid = False; break
                        if cov[d, i] == 1:                                   # Lemma U at the leaf: a hitter with exactly 2 nonzero
                            hv = False
                            for w in range(nWw):
                                if good[d, w] & ge2[d, w] & ~ge3[d, w] & ~bad[d, w] & hitb[i, w]:
                                    hv = True; break
                            if not hv:
                                valid = False; break
                    elif need[d, i] != 0:
                        valid = False; break
                if valid:
                    counts[0] += 1
                    res = leaf_general(P, cov[d], isflex, whitb, Hc, Wc, V, WS, HS, Z, Pp, Pn, allW, slots, slen, part, pos, npt, allowed)
                    if res >= 1: counts[1] += 1
                    if res == 2:
                        counts[2] += 1
                        if nleaves < cap:
                            for x in range(7): leaves[nleaves, x] = P[x]
                            for i in range(npt): leafcov[nleaves, i] = cov[d, i]
                            nleaves += 1
                        else: overflow = True
                action = 2; continue
            left = 7 - d; best = -1; bkey = 1 << 30; pruned = False
            for i in range(npt):
                if isflex[i]:
                    c = cov[d, i]
                    if c >= 3: r = 0
                    else:
                        hv = False
                        for w in range(nWw):
                            if good[d, w] & ~ge3[d, w] & ~bad[d, w] & hitb[i, w]:
                                hv = True; break
                        if hv: r = 1 if (c == 0 or c == 2) else 0
                        else: r = 3 - c                                      # Lemma U: w = 1 impossible, must reach 3
                else:
                    r = need[d, i]
                if r > left:
                    pruned = True; break
                if r > 0:
                    n = 0
                    for t in range(nnz[d]):
                        w = nzw[d, t]; n += popc(cand[d, w] & through[i, w])
                    if n < r:
                        pruned = True; break
                    if n - r < bkey:
                        best = i; bkey = n - r
            if pruned:
                action = 2; continue
            if d == 0 and p0 >= 0: best = p0                             # C2 28/09: symmetry breaking (forced first point)
            if best >= 0:
                for t in range(nnz[d]):
                    w = nzw[d, t]; Cs[d, w] = cand[d, w] & through[best, w]
                    if d == 0 and p0 >= 0: Cs[d, w] &= fmask[w]
            else:                                                        # no pending point: any remaining member
                for t in range(nnz[d]):
                    w = nzw[d, t]; Cs[d, w] = cand[d, w]
            curj[d] = 0
            action = 1; continue
        elif action == 1:
            j = -1
            while curj[d] < nnz[d]:
                w = nzw[d, curj[d]]; x = Cs[d, w]
                if x != 0:
                    low = x & (~x + U1)
                    j = w * 64 + popc(low - U1); Cs[d, w] = x ^ low
                    break
                curj[d] += 1
            if j < 0:
                action = 2; continue
            P[d] = j
            e1 = d + 1
            for w in range(nWw):
                good[e1, w] = good[d, w] & ~Bd[j, w]
                ge1[e1, w] = ge1[d, w] | NZ[j, w]
                ge2[e1, w] = ge2[d, w] | (ge1[d, w] & NZ[j, w])
                ge3[e1, w] = ge3[d, w] | (ge2[d, w] & NZ[j, w])
                bad[e1, w] = bad[d, w]
            for q in range(suppn[j]):                                        # flex points reaching coverage >= 2: their hitters die
                x = supp[j, q]
                if isflex[x] and cov[d, x] + 1 >= 2:
                    for w in range(nWw): bad[e1, w] |= hitb[x, w]
            for w in range(nWw): vm[w] = good[e1, w] & ~ge3[e1, w] & ~bad[e1, w]
            left = 7 - e1; ok = True
            for x in range(e1):
                if not ok: break
                for y in range(x + 1, e1):
                    c = P[x]; e = P[y]; same = np.uint64(0); opp = np.uint64(0)
                    for w in range(nWw):
                        zo = vm[w]
                        for f in range(e1):
                            if f != x and f != y: zo &= Z[P[f], w]
                        same |= ((Pp[c, w] & Pp[e, w]) | (Pn[c, w] & Pn[e, w])) & zo
                        opp |= ((Pp[c, w] & Pn[e, w]) | (Pn[c, w] & Pp[e, w])) & zo
                    if same == 0 or opp == 0:
                        ok = False; break
            if ok and left > 0:
                n = 0
                for w in range(nWw): n += popc(vm[w] & ~ge1[e1, w])
                if n < left * (left - 1): ok = False
                else:
                    for x in range(e1):
                        c = P[x]; n = 0
                        for w in range(nWw): n += popc(NZ[c, w] & vm[w] & ~ge2[e1, w])
                        if n < 2 * left:
                            ok = False; break
            if ok:
                for i in range(npt):
                    need[e1, i] = need[d, i]; cov[e1, i] = cov[d, i]
                cnt = 0
                for t in range(nnz[d]):
                    w = nzw[d, t]; v = cand[d, w] & compat[j, w]
                    if v != 0:
                        cand[e1, w] = v; nzw[e1, cnt] = w; cnt += 1
                nnz[e1] = cnt
                for q in range(suppn[j]):
                    x = supp[j, q]; cov[e1, x] += 1
                    sat = False
                    if isflex[x]:
                        sat = cov[e1, x] == 3
                    else:
                        need[e1, x] -= 1; sat = need[e1, x] == 0
                    if sat:
                        for t in range(cnt):
                            w = nzw[e1, t]; cand[e1, w] &= ~through[x, w]
                d = e1; action = 0; continue
            cand[d, j >> 6] &= ~(U1 << np.uint64(j & 63))
            action = 1; continue
        else:
            if d == 0: break
            d -= 1; j = P[d]
            cand[d, j >> 6] &= ~(U1 << np.uint64(j & 63))
            action = 2 if stop else 1
    return nodes, leaves[:nleaves], leafcov[:nleaves], not stop, overflow, counts



@nb.njit(cache=True)
def compat_cached(ORTH1, pos, inv, posmask, Ppb, Pnb, SAb):
    """Same output as compat_nb, but orthogonality read from the per-y1 cache ORTH1 (bits over H1); only orthogonal pairs visited."""
    nH = pos.shape[0]; nw = (nH + 63) // 64; out = np.zeros((nH, nw), np.uint64); w1 = ORTH1.shape[1]
    for i in range(nH):
        r = ORTH1[pos[i]]
        for w in range(w1):
            x = r[w] & posmask[w]
            while x != 0:
                low = x & (~x + U1); b = w * 64 + popc(low - U1); x ^= low
                j = inv[b]
                if j <= i: continue
                same = False; opp = False
                for t in range(Ppb.shape[1]):
                    if (Ppb[i, t] & Ppb[j, t]) | (Pnb[i, t] & Pnb[j, t]): same = True
                    if (Ppb[i, t] & Pnb[j, t]) | (Pnb[i, t] & Ppb[j, t]): opp = True
                if not (same and opp): continue
                ov = False
                for t in range(SAb.shape[1]):
                    if SAb[i, t] & SAb[j, t]: ov = True; break
                if ov: continue
                out[i, j >> 6] |= U1 << np.uint64(j & 63); out[j, i >> 6] |= U1 << np.uint64(i & 63)
    return out


def y1_cache(L, y1):
    """Per-y1 shared setup (C2 29/09): H1 = admissible norm-4 vectors orthogonal to y1 (X4 order kept) and their orthogonality bits."""
    k = y1.tobytes()
    if L.get("_y1c", (None,))[0] == k: return L["_y1c"][1]
    if "_A4" not in L:
        X4 = L["X4"]; L["_A4"] = X4[np.abs(X4 @ L["G"] @ L["S2"].T).max(axis=1) <= 2]
    A4 = L["_A4"]; H1 = A4[(A4 @ L["G"] @ y1) == 0]; n1 = len(H1); w1 = (n1 + 63) // 64
    Hf = H1.astype(np.float32); HG = (H1 @ L["G"]).astype(np.float32); ORTH1 = np.zeros((n1, w1), np.uint64)
    for s0 in range(0, n1, 2048):
        ORTH1[s0:s0 + 2048] = pack64(np.abs(HG[s0:s0 + 2048] @ Hf.T) < 0.5)
    idx = {H1[i].tobytes(): i for i in range(n1)}
    c = (H1, ORTH1, idx); L["_y1c"] = (k, c); return c

@nb.njit(cache=True)
def compat_fast(HG, H, Ppb, Pnb, through, t1supp, t1n):
    """Same output as compat_nb (C2 29/09): the no-common-target-1-point test is done first with bitsets (union of the candidates
    through the target-1 points of h_i), so the dot product and the Same/Opp cells are only evaluated on the survivors."""
    n = H.shape[0]; nw = (n + 63) // 64; nWw = Ppb.shape[1]
    out = np.zeros((n, nw), np.uint64); acc = np.zeros(nw, np.uint64)
    for i in range(n):
        for w in range(nw): acc[w] = 0
        for q in range(t1n[i]):
            p = t1supp[i, q]
            for w in range(nw): acc[w] |= through[p, w]
        for w in range(i >> 6, nw):
            x = ~acc[w]
            if w == (i >> 6):
                sh = (i & 63) + 1
                x = x & ~((U1 << np.uint64(sh)) - U1) if sh < 64 else np.uint64(0)
            while x != 0:
                low = x & (~x + U1); j = w * 64 + popc(low - U1); x ^= low
                if j >= n: break
                dot = 0
                for q in range(H.shape[1]): dot += HG[i, q] * H[j, q]
                if dot != 0: continue
                same = np.uint64(0); opp = np.uint64(0)
                for t in range(nWw):
                    same |= (Ppb[i, t] & Ppb[j, t]) | (Pnb[i, t] & Pnb[j, t])
                    opp |= (Ppb[i, t] & Pnb[j, t]) | (Pnb[i, t] & Ppb[j, t])
                if same != 0 and opp != 0:
                    out[i, j >> 6] |= U1 << np.uint64(j & 63); out[j, i >> 6] |= U1 << np.uint64(i & 63)
    return out


@nb.njit(cache=True)
def pack_V(V):
    """The five bitsets of V in one pass (C2 29/09): V==0, V==1, V==-1, |V|==2, |V|==1 (rows = h, bits = u); same as pack64(...)."""
    nH, nW = V.shape; nw = (nW + 63) // 64
    Z = np.zeros((nH, nw), np.uint64); Pp = np.zeros((nH, nw), np.uint64); Pn = np.zeros((nH, nw), np.uint64)
    B2 = np.zeros((nH, nw), np.uint64); NZ = np.zeros((nH, nw), np.uint64)
    for i in range(nH):
        for j in range(nW):
            v = V[i, j]; bit = U1 << np.uint64(j & 63); w = j >> 6
            if v == 0: Z[i, w] |= bit
            elif v == 1: Pp[i, w] |= bit; NZ[i, w] |= bit
            elif v == -1: Pn[i, w] |= bit; NZ[i, w] |= bit
            elif v == 2 or v == -2: B2[i, w] |= bit
    return Z, Pp, Pn, B2, NZ


def sym_reps(H, Cidx, Rs, G):
    """Union-find of the candidates through beta_0 under the reflections in Rs, vectorised (C2 29/09). Same result as the loop:
    the representative of a component is its minimum index. Hash matches are re-checked exactly (no false edge)."""
    m = len(Cidx)
    if m == 0 or len(Rs) == 0: return np.array(sorted(Cidx.tolist()), dtype=np.int64)
    C = H[Cidx].astype(np.int64); RG = Rs @ G; M = C @ RG.T
    wts = np.random.default_rng(12345).integers(1, 2 ** 62, size=C.shape[1], dtype=np.int64).view(np.uint64)
    hsh = lambda X: (X.view(np.uint64) * wts).sum(axis=1)
    kC = np.concatenate([hsh(C), hsh(-C)]); owner = np.concatenate([np.arange(m), np.arange(m)]); sgn = np.concatenate([np.ones(m, np.int64), -np.ones(m, np.int64)])
    order = np.argsort(kC, kind="stable"); kS = kC[order]
    ii, rr = np.nonzero(M)
    img = C[ii] - M[ii, rr][:, None] * Rs[rr]
    hi = hsh(np.ascontiguousarray(img)); pos = np.minimum(np.searchsorted(kS, hi), len(kS) - 1)
    hit = kS[pos] == hi                                             # C has one vector per +- line, so a hash has one owner
    o = order[pos[hit]]; j = owner[o]
    ok = (img[hit] == sgn[o][:, None] * C[j]).all(axis=1)           # exact re-check: a collision can only drop an edge
    src = ii[hit][ok]; dst = j[ok]
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    ncomp, lab = connected_components(coo_matrix((np.ones(len(src)), (src, dst)), shape=(m, m)), directed=False)
    mins = np.full(ncomp, np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(mins, lab, Cidx.astype(np.int64))
    return np.sort(mins)


@nb.njit(cache=True)
def compat_nb(HG, H, Ppb, Pnb, SAb):
    """compat[i] bit j: (h_i, h_j) = 0, Same(i,j) and Opp(i,j) non-empty, no common target-1 point (same as the float version)"""
    n = H.shape[0]; nw = (n + 63) // 64; nWw = Ppb.shape[1]; nAw = SAb.shape[1]
    out = np.zeros((n, nw), np.uint64)
    for i in range(n):
        for j in range(i + 1, n):
            dot = 0
            for q in range(H.shape[1]): dot += HG[i, q] * H[j, q]
            if dot != 0: continue
            clash = False
            for w in range(nAw):
                if SAb[i, w] & SAb[j, w]:
                    clash = True; break
            if clash: continue
            same = np.uint64(0); opp = np.uint64(0)
            for w in range(nWw):
                same |= (Ppb[i, w] & Ppb[j, w]) | (Pnb[i, w] & Pnb[j, w])
                opp |= (Ppb[i, w] & Pnb[j, w]) | (Pnb[i, w] & Ppb[j, w])
            if same != 0 and opp != 0:
                out[i, j >> 6] |= U1 << np.uint64(j & 63); out[j, i >> 6] |= U1 << np.uint64(i & 63)
    return out


@nb.njit(cache=True)
def compat_masked(ORTH, Ppb, Pnb, SAb):
    """same as compat_nb, with the orthogonality test precomputed (ORTH, from one float32 BLAS product: exact for these small integers)"""
    n = ORTH.shape[0]; nw = (n + 63) // 64; nWw = Ppb.shape[1]; nAw = SAb.shape[1]
    out = np.zeros((n, nw), np.uint64)
    for i in range(n):
        for j in range(i + 1, n):
            if not ORTH[i, j]: continue
            clash = False
            for w in range(nAw):
                if SAb[i, w] & SAb[j, w]:
                    clash = True; break
            if clash: continue
            same = np.uint64(0); opp = np.uint64(0)
            for w in range(nWw):
                same |= (Ppb[i, w] & Ppb[j, w]) | (Pnb[i, w] & Pnb[j, w])
                opp |= (Ppb[i, w] & Pnb[j, w]) | (Pnb[i, w] & Ppb[j, w])
            if same != 0 and opp != 0:
                out[i, j >> 6] |= U1 << np.uint64(j & 63); out[j, i >> 6] |= U1 << np.uint64(i & 63)
    return out


def char_trap(L, y1, y2):
    """Astra 28/09 (PROVEN, checked by C2): z characteristic of norm 16 (z = 2s, s a norm-4 shadow vector). Frame identity at z:
    sum_f (f,z)^2 = 112; the 42 u pair oddly with z (>= 42), h >= 0, so (z,y1)^2 + (z,y2)^2 <= 70, and a sum of two odd squares
    <= 70 is <= 66. Returns an explicit witness z (lattice coefficients) or None."""
    if "_ZG" not in L:
        c = solve_mod2(L["G"], np.diag(L["G"])); Z = one_per_pair(short(L["G"], 4, c))
        L["_Z"] = Z[np.einsum("ij,jk,ik->i", Z, L["G"], Z) == 16]; L["_ZG"] = L["_Z"] @ L["G"]; L["_A1"] = (None, None)
    if "_ZGf" not in L: L["_ZGf"] = L["_ZG"].astype(np.float64)       # exact: small integers, float64 BLAS
    k1 = y1.tobytes()
    if L["_A1"][0] != k1: L["_A1"] = (k1, np.rint(L["_ZGf"] @ y1.astype(np.float64)).astype(np.int64) ** 2)
    Yz = L["_A1"][1] + np.rint(L["_ZGf"] @ y2.astype(np.float64)).astype(np.int64) ** 2
    L["_Yz_last"] = (k1, y2.tobytes(), Yz)
    bad = np.flatnonzero(Yz > 66)
    return None if len(bad) == 0 else L["_Z"][bad[0]]


def search(L, y1, y2, node_limit=0, cap=2_000_000):
    t0 = time.time()
    if os.environ.get("FD_CHAR", "1") == "1":
        z = char_trap(L, y1, y2)
        if z is not None:
            return {"verdict": "EXCLUDED (char-vector trap)", "complete": True, "witness_z": [int(v) for v in z],
                    "secs": {"setup": round(time.time() - t0, 3), "dfs": 0.0, "stage2": 0.0}}
    ps, why = (pair_setup_fast if os.environ.get("FD_FASTSETUP", "1") == "1" else pair_setup)(L, y1, y2)
    if ps is None: return {"verdict": "EXCLUDED (%s)" % why}
    nb_ = int(os.environ.get("FD_BND", "1000"))   # 29/09: all boundary points (= the twin envelope's new roots): nodes x3 fewer
    if nb_ > 0 and "_ZG" in L:
        # Boundary-shadow lemma (Astra/5.6 29/09, checked by C2): z char. norm 16 with Y_z = 58 behaves like a non-y-hit new root:
        # every |(z,h)/2| <= 1, coverage in {1,3}, and 1 only if an actual u has |(z,u)| = 3. Added as extra points (2(h,s) = (h,z)).
        c_ = L.get("_Yz_last")
        if c_ is not None and c_[0] == y1.tobytes() and c_[1] == y2.tobytes(): Yz = c_[2]
        else: Yz = (L["_ZG"] @ y1) ** 2 + (L["_ZG"] @ y2) ** 2
        zi = np.flatnonzero(Yz == 58)
        if len(zi):
            HZa = np.rint(ps["H"].astype(np.float64) @ L["_ZGf"][zi].T).astype(np.int64) if "_ZGf" in L else ps["H"] @ L["_ZG"][zi].T
            cnt_ = (HZa != 0).sum(0); o_ = np.argsort(cnt_)[:nb_]; zi = zi[o_]
            ZGs = L["_ZG"][zi]; HZ = HZa[:, o_]; okh = (np.abs(HZ) <= 2).all(1)
            WZ = ps["W"] @ ZGs.T; wh = np.abs(WZ) == 3; fl = wh.any(0)
            ps["H"] = ps["H"][okh]; ps["HS"] = np.hstack([ps["HS"][okh], HZ[okh]]); ps["V"] = ps["V"][okh]
            ps["WS"] = np.hstack([ps["WS"], WZ]); ps["whit"] = np.hstack([ps["whit"], wh])
            ps["target"] = np.concatenate([ps["target"], np.where(fl, 0, 3)]); ps["flex"] = np.concatenate([ps["flex"], fl])
            ps["yhit"] = np.concatenate([ps["yhit"], np.zeros(len(zi), bool)])
    if os.environ.get("FD_WPRUNE", "0") == "1" and "_ZG" in L:
        # C2 29/09 (PROVEN, energy identity at a norm-16 characteristic z): 4 sum_i a_i^2 + sum_u b_u^2 = 112 - Y_z with
        # sum_i a_i^2 >= 1 (odd) and every b_u odd  =>  sum_u (b_u^2 - 1) <= 66 - Y_z and 4 a_i^2 <= 70 - Y_z.
        # Y_z in {50,58}: no frame u has |(z,u)| = 5 (cost 24 > 16).  Y_z = 58: |(z,h)| <= 2;  Y_z = 50: |(z,h)| <= 4.
        Yz_ = (L["_ZG"] @ y1) ** 2 + (L["_ZG"] @ y2) ** 2
        Z50 = L["_ZG"][Yz_ >= 50]
        if len(Z50):
            wdel = (np.abs(ps["W"] @ Z50.T) >= 5).any(axis=1)
            HZ_ = np.abs(ps["H"] @ Z50.T); lim_ = np.where(Yz_[Yz_ >= 50] == 58, 2, 4)
            hdel = (HZ_ > lim_[None, :]).any(axis=1)
            if wdel.any():
                kw = ~wdel
                ps["W"] = ps["W"][kw]; ps["WS"] = ps["WS"][kw]; ps["V"] = ps["V"][:, kw]; ps["whit"] = ps["whit"][kw]
                fl = (~ps["yhit"]) & ps["whit"].any(axis=0)
                ps["flex"] = fl; ps["target"] = np.where(ps["yhit"], 1, np.where(fl, 0, 3)).astype(np.int64)
            if hdel.any():
                kh = ~hdel; ps["H"] = ps["H"][kh]; ps["HS"] = ps["HS"][kh]; ps["V"] = ps["V"][kh]
            ps["_wprune"] = (int(wdel.sum()), int(hdel.sum()))
    if os.environ.get("FD_PREBUDGET") == "1" and len(L["R"]):          # C2 28/09: exact per-pair root-class budget ILP (pair_budget.py)
        import pair_budget as PB
        if "_pb_ctx" not in L: L["_pb_ctx"] = PB.setup(None, L)
        v_ = PB.pair_test(L["_pb_ctx"], y1, y2, ps)
        if v_.startswith("EXCLUDED"):
            return {"verdict": v_, "complete": True, "secs": {"setup": round(time.time() - t0, 2), "dfs": 0.0, "stage2": 0.0}}
    _tm = {"pre": time.time() - t0}; _t = time.time()
    H, HS, W, WS, V, target = ps["H"], ps["HS"], ps["W"], ps["WS"], ps["V"], ps["target"]; G = L["G"]
    nH, nW, npt = len(H), len(W), len(target)
    SUP = (HS != 0)
    HG = H @ G
    Zp_, Pp_, Pn_, B2p_, NZp_ = pack_V(np.ascontiguousarray(V.astype(np.int64))); THR_ = pack64(SUP.T)
    if os.environ.get("FD_FASTCOMPAT") == "1":
        Hf = H.astype(np.float32); ORTH = (Hf @ (G.astype(np.float32) @ Hf.T)) == 0
        CBb = compat_masked(ORTH, pack64(V == 1), pack64(V == -1), pack64(SUP & (target == 1)[None, :]))
    elif os.environ.get("FD_SHARED", "0") == "1":   # 29/09: exact but no gain (too few search pairs per y1 after the trap)
        H1, ORTH1, idx1 = y1_cache(L, y1); pos = np.array([idx1[h.tobytes()] for h in H], dtype=np.int64)
        inv = -np.ones(len(H1), dtype=np.int64); inv[pos] = np.arange(len(H)); pm = pack64(inv[None, :] >= 0)[0]
        CBb = compat_cached(ORTH1, pos, inv, pm, pack64(V == 1), pack64(V == -1), pack64(SUP & (target == 1)[None, :]))
        if os.environ.get("FD_CHECK_SHARED"):
            assert np.array_equal(CBb, compat_nb(HG, H, pack64(V == 1), pack64(V == -1), pack64(SUP & (target == 1)[None, :]))), "shared compat differs"
    elif os.environ.get("FD_FASTCOMPAT2", "1") == "1":
        T1 = SUP & (target == 1)[None, :]; t1n = T1.sum(axis=1).astype(np.int64)
        t1s = np.zeros((nH, max(1, int(t1n.max()) if nH else 1)), dtype=np.int64)
        rr, cc = np.nonzero(T1); t1s[rr, np.arange(len(rr)) - np.repeat(np.cumsum(t1n) - t1n, t1n)] = cc
        CBb = compat_fast(HG, H, Pp_, Pn_, THR_, t1s, t1n)
        if os.environ.get("FD_CHECK_COMPAT2"):
            assert np.array_equal(CBb, compat_nb(HG, H, pack64(V == 1), pack64(V == -1), pack64(SUP & (target == 1)[None, :]))), "fast compat differs"
    else:
        CBb = compat_nb(HG, H, pack64(V == 1), pack64(V == -1), pack64(SUP & (target == 1)[None, :]))
    if os.environ.get("FD_CHECK_COMPAT"):                               # float version, for the cross-check
        Pf = (V == 1).astype(np.float32); Nf = (V == -1).astype(np.float32); SA = (SUP & (target == 1)[None, :]).astype(np.float32)
        CB = np.zeros((nH, nH), dtype=bool)
        for s in range(0, nH, 512):
            e = min(nH, s + 512)
            CB[s:e] = ((HG[s:e] @ H.T) == 0) & ((Pf[s:e] @ Pf.T + Nf[s:e] @ Nf.T) > 0.5) & ((Pf[s:e] @ Nf.T + Nf[s:e] @ Pf.T) > 0.5) \
                & ((SA[s:e] @ SA.T) < 0.5)
        assert np.array_equal(pack64(CB), CBb), "compat: numba and float versions differ"
    _tm["compat"] = time.time() - _t; _t = time.time()
    supn = SUP.sum(axis=1); supp = np.full((nH, int(supn.max()) if nH else 1), -1, dtype=np.int64)
    rr_, cc_ = np.nonzero(SUP); supp[rr_, np.arange(len(rr_)) - np.repeat(np.cumsum(supn) - supn, supn)] = cc_
    u = {"compat": CBb, "through": THR_, "Z": Zp_, "Pp": Pp_, "Pn": Pn_, "Bd": B2p_, "NZ": NZp_}
    if os.environ.get("FD_EC", "1") == "1" and nW:
        # Endpoint confinement (29/09, PROVEN, checked by C2): if h covers a point hit at 3/2 by u, then (u,h) != 0 in a frame.
        T_ = (SUP.astype(np.float32) @ ps["whit"].T.astype(np.float32)) > 0.5          # nH x nW: u hits a point covered by h
        if os.environ.get("FD_ECSIGN", "0") != "1":
            u["Bd"] = B2p_ | (pack64(T_) & ~NZp_)                           # = pack64(|V|==2 | (T_ & |V|!=1)); padding bits stay 0
        badm = None if os.environ.get("FD_ECSIGN", "0") != "1" else (np.abs(V) == 2) | (T_ & (np.abs(V) != 1))
        if os.environ.get("FD_ECSIGN", "0") == "1":
            # signed EC (same proof): if h covers p and u hits p at 3/2 then (u,h) = sign(2(s,u)) * (h,s) at EVERY such p
            Hp = (HS // 2).astype(np.float32); Wp = (np.sign(ps["WS"]) * (np.abs(ps["WS"]) == 3)).astype(np.float32)
            Ssum = Hp @ Wp.T; Ncnt = np.abs(Hp) @ np.abs(Wp).T
            badm |= (Ncnt > 0.5) & (np.abs(Ssum - V * Ncnt) > 0.5)
            u["Bd"] = pack64(badm)
    allM = pack64(np.ones((1, nH), dtype=bool))[0]; allW = pack64(np.ones((1, nW), dtype=bool))[0]
    _tm["packs"] = time.time() - _t; _t = time.time()
    t1 = time.time()
    general = bool(ps["flex"].any()) or bool(os.environ.get("FD_FORCE_GENERAL"))   # forcing it on a rigid lattice = control
    p0, fmask = -1, allM
    if os.environ.get("FD_SYM") == "1" and general and (target == 1).any():
        # C2 28/09: symmetry breaking. G_ = reflections in the roots of U orthogonal to y1, y2 and to beta_0, where beta_0 is a
        # target-1 point (covered EXACTLY once by every solution). G_ fixes y1, y2, beta_0 and preserves every constraint, so it
        # suffices to branch on one representative per G_-orbit of the candidates through beta_0.
        cntp = SUP.sum(axis=0); cand_pts = np.flatnonzero(target == 1)
        p0 = int(cand_pts[np.argmin(cntp[cand_pts])])
        Y = ps["Y"]; RT = L["R"]; b0 = L["S2"][p0]
        Rs = RT[((RT @ G @ Y.T) == 0).all(axis=1) & ((RT @ G @ b0) == 0)]
        Cidx = np.flatnonzero(SUP[:, p0])
        if os.environ.get("FD_SYMOLD") != "1":
            reps = sym_reps(H, Cidx, Rs, G); fm = np.zeros(nH, dtype=bool); fm[reps] = True; fmask = pack64(fm[None, :])[0]
    if os.environ.get("FD_SYM") == "1" and general and (target == 1).any() and os.environ.get("FD_SYMOLD") == "1":
        keyi = {}
        for i in Cidx: keyi[H[i].tobytes()] = i; keyi[(-H[i]).tobytes()] = i
        par = {i: i for i in Cidx}
        def fnd(i):
            while par[i] != i: par[i] = par[par[i]]; i = par[i]
            return i
        RG = Rs @ G
        for i in Cidx:
            for r_, v_ in zip(Rs, H[i] @ RG.T):
                if v_ == 0: continue
                j = keyi.get((H[i] - v_ * r_).tobytes())
                if j is not None: a_, b_ = fnd(i), fnd(j); par[max(a_, b_)] = min(a_, b_)
        reps = sorted({fnd(i) for i in Cidx}); fm = np.zeros(nH, dtype=bool); fm[reps] = True; fmask = pack64(fm[None, :])[0]
    _tm["sym"] = time.time() - _t
    if general: cap = min(cap, int(os.environ.get("FD_CAP_GENERAL", "50000")))    # leafcov = cap x points int64 (2e6 x 136 = 2.2 GB!)
    if not general:
        nodes, leaves, complete, overflow, counts = dfs_generic(u["compat"], u["through"], supp, supn.astype(np.int64), u["Z"], u["Pp"], u["Pn"],
                                                                u["Bd"], u["NZ"], target, allM, allW, node_limit, cap, H, W,
                                                                V.astype(np.int64), WS.astype(np.int64), HS.astype(np.int64), -1, allM)
        lcov = [None] * len(leaves)
    else:
        if os.environ.get("FD_LEMMA_U"):                                  # C1 notes/145 U+ (strong), per-node bitset version (C2)
            hitb = pack64(ps["whit"].T)                                       # points x words
            nodes, leaves, lcov, complete, overflow, counts = dfs_general_u(u["compat"], u["through"], supp, supn.astype(np.int64), u["Z"], u["Pp"],
                                                                          u["Pn"], u["Bd"], u["NZ"], target, ps["flex"], pack64(ps["whit"]), allM,
                                                                          allW, node_limit, cap, H, W, V.astype(np.int64), WS.astype(np.int64),
                                                                          HS.astype(np.int64), hitb, p0, fmask)
        else:
            nodes, leaves, lcov, complete, overflow, counts = dfs_general(u["compat"], u["through"], supp, supn.astype(np.int64), u["Z"], u["Pp"],
                                                                          u["Pn"], u["Bd"], u["NZ"], target, ps["flex"], pack64(ps["whit"]), allM,
                                                                          allW, node_limit, cap, H, W, V.astype(np.int64), WS.astype(np.int64),
                                                                          HS.astype(np.int64), p0, fmask)
    assert not overflow
    t2 = time.time()
    joint_ok, full = [], []
    for Pl, cv in zip(leaves, lcov):
        F1 = None if cv is None else (ps["flex"] & (np.asarray(cv) == 1))
        sols = joint(list(Pl), H, W, V, ps, G, F1=F1)
        joint_ok.append(len(sols))
        full += [full_check(list(Pl), s, H, W, V, ps, G) for s in sols]
    out = {"nH": nH, "nW": nW, "points": npt, "yhits": int(ps["yhit"].sum()), "flex": int(ps["flex"].sum()), "nodes": [int(x) for x in nodes], "complete": bool(complete),
           "counts": {"leaves": int(counts[0]), "lemmaZ": int(counts[1]), "local_star_sums": int(counts[2]),
                      "joint_star_sums": sum(1 for x in joint_ok if x), "frame": sum(1 for f in full if f == 3)},
           "tm": _tm, "secs": {"setup": round(t1 - t0, 1), "dfs": round(t2 - t1, 1), "stage2": round(time.time() - t2, 1)}}
    c = out["counts"]
    out["verdict"] = ("INCOMPLETE" if not complete else "FRAME FOUND" if c["frame"] else "EXCLUDED (stage 1 empty: no 7 h)" if not c["leaves"]
                      else "EXCLUDED (Lemma Z)" if not c["lemmaZ"] else "EXCLUDED (no local star sums)" if not c["local_star_sums"]
                      else "EXCLUDED (no joint star sums)" if not c["joint_star_sums"] else "EXCLUDED (joint star sums, but (F) fails)")
    return out


def joint(P, H, W, V, ps, G, cap=1000, F1=None):
    """all choices of the 42 words (Same and Opp per pair of blocks, distinct) with every star sum = 3 h_c (exact).
    GENERAL mode (F1 = boolean array of the flex points met once): only words whose 3/2-hits lie in F1, and the chosen words
    must hit every point of F1 exactly once.
    Pruning (27/09, as in o24_dfs.joint and leaf_l7): for each block x, the partial star sum read on the new roots, 2(.,beta),
    must still reach 3 * 2(h_x, beta) with the remaining words (each adds an odd number: +-1, or +-3 at an F1 point it hits)."""
    Pw = [(x, y) for x, y in itertools.combinations(range(7), 2)]
    others = lambda x, y: [P[f] for f in range(7) if f not in (x, y)]
    allowed = np.ones(len(W), dtype=bool) if F1 is None else ~(ps["whit"] & ~F1[None, :]).any(axis=1)
    WS = W @ G @ ps["S2"].T if "S2" in ps else None
    slots = []
    for x, y in Pw:
        c, d = P[x], P[y]; zero = np.all(V[others(x, y)] == 0, axis=0) & allowed
        same = np.flatnonzero(zero & (np.abs(V[c]) == 1) & (V[c] == V[d])); opp = np.flatnonzero(zero & (np.abs(V[c]) == 1) & (V[c] == -V[d]))
        slots += [(x, y, same), (x, y, opp)]
    slots.sort(key=lambda s: len(s[2]))
    tgt = [3 * H[c] for c in P]; part = [np.zeros(len(H[0]), dtype=np.int64) for _ in range(7)]; rem = [12] * 7
    WSb = ps["WS"]; HSb = ps["HS_P"] if "HS_P" in ps else None
    tgt_ws = [3 * ps["HSfull"][c] for c in P]
    part_ws = [np.zeros(WSb.shape[1], dtype=np.int64) for _ in range(7)]
    fl = np.zeros(WSb.shape[1], dtype=np.int64) if F1 is None else 2 * F1.astype(np.int64)
    sols, choice = [], [None] * len(slots)

    def feas(x):
        dl = tgt_ws[x] - part_ws[x]; r = rem[x]
        if (dl > r + fl).any() or (dl < -r - fl).any(): return False
        return not ((dl - r) % 2).any()

    def go(i):
        if len(sols) >= cap: return
        if i == len(slots):
            if all(np.array_equal(part[x], tgt[x]) for x in range(7)):
                if F1 is not None:
                    hc = ps["whit"][choice].sum(axis=0)
                    if not np.array_equal(hc, F1.astype(np.int64)): return
                sols.append(list(choice))
            return
        x, y, ts = slots[i]
        for t in ts:
            if t in choice[:i]: continue
            sx = 1 if V[P[x], t] > 0 else -1; sy = 1 if V[P[y], t] > 0 else -1
            part[x] += sx * W[t]; part[y] += sy * W[t]; part_ws[x] += sx * WSb[t]; part_ws[y] += sy * WSb[t]
            rem[x] -= 1; rem[y] -= 1; choice[i] = t
            if feas(x) and feas(y) and (rem[x] or np.array_equal(part[x], tgt[x])) and (rem[y] or np.array_equal(part[y], tgt[y])): go(i + 1)
            part[x] -= sx * W[t]; part[y] -= sy * W[t]; part_ws[x] -= sx * WSb[t]; part_ws[y] -= sy * WSb[t]
            rem[x] += 1; rem[y] += 1; choice[i] = None
    go(0)
    return [(slots, s) for s in sols]


def full_check(P, sol, H, W, V, ps, G):
    """1: T^2 + T fails; 2: T ok, frame identity fails; 3: FRAME"""
    slots, s = sol
    U, Cm = [], np.zeros((7, 42), dtype=np.int64)
    for o, ((x, y, ts), t) in enumerate(zip(slots, s)):
        u = W[t] * (1 if V[P[x], t] > 0 else -1); U.append(u)       # (u, h_x) = +1 -> C_xo = -1
        Cm[x, o] = -1; Cm[y, o] = -int(u @ G @ H[P[y]])
    Um = np.array(U, dtype=np.int64); GU = Um @ G @ Um.T
    T = -(GU - np.diag(np.diag(GU)))
    if not np.array_equal(T @ T + T, 12 * np.eye(42, dtype=np.int64) - Cm.T @ Cm): return 1
    Hm = H[P]; Y = ps["Y"]
    S = Hm.T @ Hm + Um.T @ Um + Y.T @ Y
    Gi = np.rint(np.linalg.inv(G.astype(float))).astype(np.int64); assert np.array_equal(G @ Gi, np.eye(24, dtype=np.int64))
    return 3 if np.array_equal(S, 7 * Gi) else 2


# ---------------------------------------------------------------- O24 control (C1's Golay coordinates, scale 1/8)
def o24_gens():
    sys.path.insert(0, os.path.join(LAB, "scripts")); from o24_dfs import DAT
    gens = [[2 if (b >> i) & 1 else 0 for i in range(24)] for b in DAT["golay_basis"]]
    for i in range(1, 24): gens.append([4 if j in (0, i) else 0 for j in range(24)])
    gens.append([4] + [-4] + [0] * 22); gens.append([1] * 24)
    return gens, DAT


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "--o24":
        gens, DAT = o24_gens(); L = lattice_data(gens, F(1, 8))
        print("O24: roots %d, norm-3 lines %d, norm-4 lines %d, new-root pairs %d, rigid %s (%.1fs)" % (L["r"], len(L["X3"]), len(L["X4"]), len(L["S2"]), L["rigid"], L["secs"]), flush=True)
        Binv = np.linalg.inv(np.array([[float(x) for x in row] for row in L["B"]]))
        lim = int(os.environ.get("NODE_LIMIT", "0"))
        for k in map(int, a[1:]):
            cl = DAT["orbits"][k]
            y1, y2 = (np.rint(np.array(cl[key], dtype=float) @ Binv).astype(np.int64) for key in ("a1", "a2"))
            for y, key in ((y1, "a1"), (y2, "a2")):
                assert np.allclose(y @ np.array([[float(x) for x in row] for row in L["B"]]), cl[key])
            res = search(L, y1, y2, lim)
            print("O24 class %d via frame_dfs: %s" % (k, json.dumps(res)), flush=True)
    elif a and a[0] == "--batch":                                      # --batch U.json|- pairs.json out.jsonl [from] [to]: resumable
        # U given by U/NNN.json ("gens", standard inner product), or "-" to take it from the pairs file itself
        # ("basis" rows in ambient coordinates and "scale" s, inner product a.b / s, as in structures' Golay-frame files)
        P_ = json.load(open(a[2])); pairs = P_["pairs"]; outp = a[3]
        if a[1] == "-": gens, sc = P_["basis"], F(1, int(P_.get("scale", 1)))
        else: gens, sc = json.load(open(a[1]))["gens"], F(1)
        k0 = int(a[4]) if len(a) > 4 else 0; k1 = min(int(a[5]) if len(a) > 5 else len(pairs), len(pairs))
        lim = int(os.environ.get("FD_NODE_LIMIT", "300000000"))            # per pair; a pair hitting it is reported INCOMPLETE
        done = set()
        if os.path.exists(outp):
            for l in open(outp): done.add(json.loads(l)["k"])
        L = lattice_data(gens, sc); Bq = np.array([[F(x) for x in row] for row in L["B"]], dtype=object)
        Bm = np.array([[float(x) for x in row] for row in L["B"]]); Binv = np.linalg.inv(Bm)
        print("U %s: r %d, rigid %s, new-root pairs %d (%.1fs); pairs %d..%d, %d done, node limit %d" % (a[1], L["r"], L["rigid"], len(L["S2"]), L["secs"], k0, k1, len(done), lim), flush=True)
        for k in range(k0, k1):
            if k in done: continue
            yy = []
            for c in pairs[k]:
                v = np.array([float(F(x)) for x in c]); cf = np.rint(v @ Binv).astype(np.int64)
                assert np.array_equal(cf @ Bq, np.array([F(x) for x in c], dtype=object)), "y not in U"
                yy.append(cf)
            t_ = time.time(); r = search(L, yy[0], yy[1], node_limit=lim); r["k"] = k; r["secs_total"] = round(time.time() - t_, 2)
            open(outp, "a").write(json.dumps(r) + chr(10))
        print("batch done", flush=True)
    else:
        rec = json.load(open(a[0])); pairs = json.load(open(a[1]))["pairs"]
        L = lattice_data(rec["gens"])
        print("U %s: r %d, rigid %s, new-root pairs %d (%.1fs)" % (a[0], L["r"], L["rigid"], len(L["S2"]), L["secs"]), flush=True)
        Binv = [[F(x) for x in row] for row in np.linalg.inv(np.array([[float(x) for x in row] for row in L["B"]]))]
        for k in map(int, a[2:]):
            yy = []
            for c in pairs[k]:
                v = np.array([float(F(x)) for x in c]); cf = np.rint(v @ np.array(Binv, dtype=float)).astype(np.int64)
                assert np.allclose(cf @ np.array([[float(x) for x in row] for row in L["B"]]), v); yy.append(cf)
            print("pair %d: %s" % (k, json.dumps(search(L, yy[0], yy[1]))), flush=True)
