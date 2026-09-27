"""Stage B for V = Z^17 (structures, 27/09 night): given a binary skeleton (pi_1 = 0, pi_2..pi_7), decide whether the 21 cells can be
filled.  Exact search; every pruning rule below is a necessary condition (RANK17_HANDOFF sections 4-6 + lemmas L1-L2).

Cell (c,d), c<d (0-based blocks): an ordered signed pair (w+, w-) of norm-10 vectors of K = Z^17 cap y^perp with odd support
pi_c ^ pi_d, such that
  - (w+, w-) in {0, +-4}, w+ != +-w-;
  - a = (w+ + w-)/2 has p(a) = eps_c (then b = (w+ - w-)/2 has p(b) = eps_d automatically);   [G_cd = e_c + a in U_H]
  - L1 (refined): no norm-1 vector in u + U_H, i.e. for w of shape (2,2,1,1) (odd {i,j}, twos {k,l}):
    (col_k ^ col_l ^ <col_i, col_j>) avoids {0, eps_c, eps_d, eps_c+eps_d} mod 1;  shape (3,1) is excluded outright.
Between cells X, Y: (u_s, v_t) = ((w_s, w_t) + r_s.r_t)/4 in {0, +-1} for s,t in {+,-}  (r+ = e_c + e_d, r- = e_c - e_d).
Stars: for each c, sum_{d>c} 2 a_cd + sum_{d<c} 2 b_dc = 0.   Frame: sum (w+ w+^T + w- w-^T) = 28 I - 4 (y1 y1^T + y2 y2^T).
Sign normalisation (flips of h_c act as a_cd -> -a_cd for all d>c, b_dc -> -b_dc for all d<c): a_01 and b_01 and b_0d (d>=2)
have a positive first non-zero coordinate.  Every orbit keeps exactly one representative.
"""
import itertools, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from stageA import load_pair

CELLS = [(c, d) for c in range(7) for d in range(c + 1, 7)]
CIDX = {cd: k for k, cd in enumerate(CELLS)}


def canon7(x):
    return x ^ 127 if x & 1 else x


CANON7 = np.array([x ^ 127 if x & 1 else x for x in range(128)], dtype=np.int64)


def rvec(c, d, s):
    r = np.zeros(7, np.int64); r[c] = 1; r[d] = s
    return r


class Z17:
    def __init__(self, pid):
        self.pid = pid
        self.y1, self.y2 = load_pair(pid)
        self.Y = np.array([self.y1, self.y2], dtype=np.int64)
        self.Yi = np.array([self.y1[i] ** 2 + self.y2[i] ** 2 for i in range(17)])
        self.F = 28 * np.eye(17, dtype=np.int64) - 4 * (np.outer(self.y1, self.y1) + np.outer(self.y2, self.y2))
        self.cache = {}
        self.pow2 = (1 << np.arange(17)).astype(np.int64)

    def class_vectors(self, S):
        """All norm-10 w in K with odd support S (bitmask); shape (3,1) excluded (L1)."""
        if S in self.cache:
            return self.cache[S]
        idx = [i for i in range(17) if S >> i & 1]
        rest = [i for i in range(17) if not S >> i & 1]
        n = len(idx)
        signs = np.array(list(itertools.product((1, -1), repeat=n)), dtype=np.int64)
        blocks = []
        if n == 10:
            W = np.zeros((len(signs), 17), np.int64); W[:, idx] = signs; blocks.append(W)
        elif n == 6:
            for k in rest:
                for s2 in (2, -2):
                    W = np.zeros((len(signs), 17), np.int64); W[:, idx] = signs; W[:, k] = s2; blocks.append(W)
        elif n == 2:
            for k, l in itertools.combinations(rest, 2):
                for sk in (2, -2):
                    for sl in (2, -2):
                        W = np.zeros((len(signs), 17), np.int64); W[:, idx] = signs; W[:, k] = sk; W[:, l] = sl
                        blocks.append(W)
        else:
            raise ValueError(S)
        W = np.vstack(blocks)
        W = W[(W @ self.Y.T == 0).all(axis=1)]
        self.cache[S] = W
        return W


def pvals(masks, pis):
    """p(x) as canonical 7-bit functional for an array of bitmasks x (numpy int64)."""
    out = np.zeros(len(masks), np.int64)
    for e, pe in enumerate(pis):
        par = np.zeros(len(masks), np.int64)
        v = masks & pe
        # popcount parity
        while True:
            nz = v != 0
            if not nz.any():
                break
            par ^= (v & 1)
            v = v >> 1
        out |= par << e
    return np.where(out & 1, out ^ 127, out)


def first_pos(A):
    """True where the first non-zero coordinate of each row is positive."""
    nz = A != 0
    j = nz.argmax(axis=1)
    return A[np.arange(len(A)), j] > 0


class Skeleton:
    def __init__(self, z, pis):
        self.z = z
        self.pis = pis                      # 7 masks, pis[0] = 0
        self.cols = [sum(((pis[e] >> i) & 1) << e for e in range(7)) for i in range(17)]
        self.dom = {}
        self.info = {}
        self.PI = np.array([[(pis[e] >> i) & 1 for i in range(17)] for e in range(7)], dtype=np.int64)
        self.P4 = (4 ** np.arange(7)).astype(np.int64)

    def build(self):
        z, pis = self.z, self.pis
        for (c, d) in CELLS:
            S = pis[c] ^ pis[d]
            W = z.class_vectors(S)
            odd = [i for i in range(17) if S >> i & 1]
            # L1 refined for N = 2 cells
            if len(odd) == 2:
                i, j = odd
                forb = {0, canon7(1 << c), canon7(1 << d), canon7((1 << c) ^ (1 << d))}
                span = {0, self.cols[i], self.cols[j], self.cols[i] ^ self.cols[j]}
                if len(W):
                    tw = np.nonzero(np.abs(W) == 2)[1].reshape(len(W), 2)
                    colsa = np.array(self.cols, dtype=np.int64)
                    base = colsa[tw[:, 0]] ^ colsa[tw[:, 1]]
                    bad = np.zeros(len(W), bool)
                    for s_ in span:
                        bad |= np.isin(CANON7[base ^ s_], list(forb))
                    W = W[~bad]
            n = len(W)
            if n < 2:
                self.dom[(c, d)] = (np.zeros((0, 17), np.int64), np.zeros((0, 17), np.int64))
                self.info[(c, d)] = (len(odd), n, 0)
                continue
            # M1 keyed: t_e(w) = (w . pi_e) mod 4; p(a)_e = ((t_e(w+) + t_e(w-)) / 2) mod 2 must equal target_e, where the
            # target is eps_c in canonical form (bit of block 1 is 0 because pi_1 = 0): block 1 -> all other bits, else bit c.
            tgt = np.array([(e != 0) if c == 0 else (e == c) for e in range(7)], dtype=np.int64)
            T = (W @ self.PI.T) % 4                       # (n, 7)
            key = (T * self.P4).sum(axis=1)
            need = (((2 * tgt - T) % 4) * self.P4).sum(axis=1)
            groups = {}
            for idx, k in enumerate(key.tolist()):
                groups.setdefault(k, []).append(idx)
            Is, Js = [], []
            for k, I0 in groups.items():
                I0 = np.array(I0)
                nk = int(need[I0[0]])
                J0 = groups.get(nk)
                if not J0:
                    continue
                J0 = np.array(J0)
                G = W[I0] @ W[J0].T
                a_, b_ = np.nonzero((G == 0) | (G == 4) | (G == -4))
                Is.append(I0[a_]); Js.append(J0[b_])
            I = np.concatenate(Is) if Is else np.zeros(0, np.int64)
            J = np.concatenate(Js) if Js else np.zeros(0, np.int64)
            Wp, Wm = W[I], W[J]
            a, b = (Wp + Wm) // 2, (Wp - Wm) // 2
            keep = np.ones(len(I), bool)
            if c == 0:
                keep &= first_pos(b)
                if d == 1:
                    keep &= first_pos(a)
            self.dom[(c, d)] = (Wp[keep], Wm[keep])
            self.info[(c, d)] = (len(odd), n, int(keep.sum()))
        return self

    def sizes(self):
        return [len(self.dom[cd][0]) for cd in CELLS]


def rr(X, Y):
    (c, d), (e, f) = X, Y
    out = {}
    for s in (1, -1):
        for t in (1, -1):
            out[(s, t)] = int(rvec(c, d, s) @ rvec(e, f, t))
    return out


RR = {(X, Y): rr(X, Y) for X in CELLS for Y in CELLS if X != Y}


def compat(Xw, Yd, X, Y):
    """Xw = (w+, w-) chosen for cell X; Yd = (P+, P-) domain arrays of cell Y.  Boolean mask over Y's domain."""
    r = RR[(X, Y)]
    ok = np.ones(len(Yd[0]), bool)
    for s, ws in ((1, Xw[0]), (-1, Xw[1])):
        for t, P in ((1, Yd[0]), (-1, Yd[1])):
            g = P @ ws + r[(s, t)]
            ok &= (g == 0) | (g == 4) | (g == -4)
    return ok


def solve(sk, node_limit=10 ** 7, stats=None):
    """DFS over cells in CELLS order with forward checking; star check when a vertex closes; PSD residual.
    Returns 'EXCLUDED' (no filling), 'FOUND' (a filling; returned in stats['sol']), or 'LIMIT'."""
    z = sk.z
    F = z.F.astype(float)
    order = CELLS[:]
    doms = {cd: sk.dom[cd] for cd in CELLS}
    if stats is None:
        stats = {}
    stats.setdefault('nodes', 0)
    last_at = {}
    for k, (c, d) in enumerate(order):
        last_at[c] = k; last_at[d] = k
    closes = {}
    for v, k in last_at.items():
        closes.setdefault(k, []).append(v)
    chosen = {}

    def star_ok(v):
        s = np.zeros(17, np.int64)
        for (c, d), (wp, wm) in chosen.items():
            if c == v:
                s += wp + wm
            elif d == v:
                s += wp - wm
        return not s.any()

    stats.setdefault('depth', [0] * (len(order) + 1))

    def rec(k, doms, R):
        stats['nodes'] += 1
        stats['depth'][k] += 1
        if stats['nodes'] > node_limit:
            return 'LIMIT'
        if k == len(order):
            if np.abs(R).max() < 1e-9:
                stats['sol'] = {str(cd): [chosen[cd][0].tolist(), chosen[cd][1].tolist()] for cd in CELLS}
                return 'FOUND'
            return None
        X = order[k]
        P, M = doms[X]
        res = None
        for t in range(len(P)):
            wp, wm = P[t], M[t]
            R2 = R - np.outer(wp, wp) - np.outer(wm, wm)
            if np.diag(R2).min() < -1e-9:
                continue
            if np.linalg.eigvalsh(R2).min() < -1e-7:
                continue
            chosen[X] = (wp, wm)
            if all(star_ok(v) for v in closes.get(k, [])):
                nd = {}
                dead = False
                for Y in order[k + 1:]:
                    m = compat((wp, wm), doms[Y], X, Y)
                    if not m.any():
                        dead = True; break
                    nd[Y] = (doms[Y][0][m], doms[Y][1][m])
                if not dead:
                    # star look-ahead: a vertex with exactly one open cell fixes that cell's a (or b)
                    rest = order[k + 1:]
                    for v in range(7):
                        open_v = [Y for Y in rest if v in Y]
                        if len(open_v) != 1:
                            continue
                        Y = open_v[0]
                        s = np.zeros(17, np.int64)
                        for (c, d), (xp, xm) in chosen.items():
                            if c == v:
                                s += xp + xm
                            elif d == v:
                                s += xp - xm
                        P2, M2 = nd[Y]
                        contrib = (P2 + M2) if Y[0] == v else (P2 - M2)
                        m = (contrib == -s).all(axis=1)
                        if not m.any():
                            dead = True; break
                        nd[Y] = (P2[m], M2[m])
                if not dead:
                    r = rec(k + 1, nd, R2)
                    if r in ('FOUND', 'LIMIT'):
                        del chosen[X]
                        return r
            del chosen[X]
        return res

    r = rec(0, doms, F)
    return r if r else 'EXCLUDED'


def parse_raw(line):
    cols = [int(line[2 * i:2 * i + 2], 16) for i in range(17)]
    pis = [sum(((cols[i] >> e) & 1) << i for i in range(17)) for e in range(7)]
    return cols, pis


if __name__ == '__main__':
    pid = int(sys.argv[1]); fn = sys.argv[2]
    nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    z = Z17(pid)
    t0 = time.time()
    for n, line in enumerate(open(fn)):
        if n >= nmax:
            break
        cols, pis = parse_raw(line.strip())
        sk = Skeleton(z, pis).build()
        sz = sk.sizes()
        st = {}
        v = solve(sk, 10 ** 6, st) if min(sz) > 0 else 'EXCLUDED(empty cell)'
        print(n, 'N', [bin(pis[c] ^ pis[d]).count('1') for c, d in CELLS], 'dom', sz, v, st.get('nodes'), '%.1fs' % (time.time() - t0),
              flush=True)
