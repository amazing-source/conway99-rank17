r"""wfe.py -- witness-filtered per-class frame search (structures, 28/09 night; filter validated by C1).

Setting: a class = (V, y-pair), K = V cap y^perp (rank 15, det 49), frame = 21 cells (w+, w-) of norm-10 vectors of K with
sum_u w_u w_u^T = 28 P_W (tight frame).  For every kappa in K* (dual coordinates x, 7|kappa|^2 = n):
    a_u = (kappa, w_u) in Z,  sum_u a_u^2 = 28|kappa|^2 = 4n  (Parseval),  |a_u| <= sqrt(10 n / 7)  (Cauchy-Schwarz),
    a_u = p_c + p_d mod 2 on the labels of cell (c,d), p_c = (kappa, g_c) mod 2 (g_1 = 0): the cut of the skeleton,
    k = min(o, 7 - o), o = #{d : p_d = 1}; 2k(7-k) odd labels.
Rules (C1-validated from PROVEN inputs; n = 7 uses Parseval + C-S + min N* only):
    n=5: k = 2; odd labels +-1, even labels 0.
    n=6: k in {2,3}; odd +-1; even 0 (k=3) or exactly one +-2 over the frame (k=2).
    n=7: k in {1,2,3}; odd in {+-1,+-3} (+-1 only if k=3), even in {0,+-2}; Parseval.
    n=8: k in {1,2,3} (k=0 is R8 k=0); odd +-1; even in {0,+-2}; Parseval (5/3/2 labels at +-2).
Pipeline per class:
    1. norm-10 shell of K (C1's Case), per-vector filter (k-independent allowed |a|) -> filtered shell;
    2. usable classes + M1 masks on the filtered shell (C1's usable_classes; basic L1 inside);
    3. skeleton search (C1's search with H-reduction, M1, unit rule) + witness o-rules;
    4. per skeleton: C1's domains (M1, refined L1, sign normalisation), exact-k value filter (parity asserted), and the
       linear budgets (Parseval per witness, cross sums 28(kappa,kappa') per witness pair) -> sbdfs_w.exe
       (pairwise rule, star closure, PSD residual, budgets).
Usage: python wfe.py <tag> <pair ids comma-separated> [max_skeletons] [solver node cap] [dyn 0/1] [skeleton node limit]"""
import json, os, sys, time, struct, subprocess
import numpy as np

sys.path.insert(0, 'C:/home/researcher/Desktop/c1_recherche/code')
import rank17_engine as eng
import rank17_skeleton as sk
from rank17_ypairs import short_vectors_chunked
from rank17_vsearch import CELLS
from envelope_trap import inverse_exact

sys.path.insert(0, 'C:/home/researcher/Desktop/conway99-rank17/z17')
import stageB

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, 'sbdfs_w.exe')
WORK = 'C:/home/researcher/AppData/Local/Temp/claude/runs/wfe'
NMAX = 8

VEC_ALLOWED = {5: ({0}, {1}), 6: ({0, 2}, {1}), 7: ({0, 2}, {1, 3}), 8: ({0, 2}, {1})}   # (even, odd), any k
K_ALLOWED = {5: {2}, 6: {2, 3}, 7: {1, 2, 3}, 8: {1, 2, 3}}


def o_allowed(n):
    return frozenset(o for o in range(7) if min(o, 7 - o) in K_ALLOWED[n])


def label_allowed(n, k):
    """allowed |a| on (even labels, odd labels) given n and the cut size k of the skeleton"""
    if n == 5:
        return {0}, {1}
    if n == 6:
        return ({0} if k == 3 else {0, 2}), {1}
    if n == 7:
        return {0, 2}, ({1} if k == 3 else {1, 3})
    return {0, 2}, {1}


def witnesses(GK, nmax=NMAX):
    """kappa in K* with 7|kappa|^2 <= nmax, up to sign: (D = 49 GK^-1, list of (x, n))"""
    Gi = inverse_exact([[int(x) for x in row] for row in GK])
    assert all((x * 49).denominator == 1 for row in Gi for x in row)
    D = np.array([[int(x * 49) for x in row] for row in Gi], dtype=np.int64)
    X, nrm = short_vectors_chunked(D.tolist(), 7 * nmax)
    out = []
    for x, v in zip(X, nrm):
        v = int(v)
        if v == 0:
            continue
        assert v % 7 == 0
        out.append((np.array(x, dtype=np.int64), v // 7))
    return D, out


class Wcase(eng.Case):
    def __init__(self, tag, pid, nmax=NMAX):
        t0 = time.time()
        super().__init__(tag, pid, with_masks=False)
        self.tag = tag
        self.D, wit = witnesses(self.GK, nmax)
        self.wx = np.array([x for x, _ in wit], dtype=np.int64).reshape(-1, 15)
        self.wn = [n for _, n in wit]
        assert not self.wn or min(self.wn) >= 5, 'dual vector with 7|kappa|^2 < 5: handled elsewhere'
        self.xmask = [int(((x % 2) << np.arange(15)).sum()) for x in self.wx]
        W = self.W.astype(np.int64)
        vals = W @ self.wx.T
        ok = np.ones(len(W), dtype=bool)
        for t, n in enumerate(self.wn):
            ev, od = VEC_ALLOWED[n]
            a = np.abs(vals[:, t])
            ok &= np.where(a % 2 == 1, np.isin(a, list(od)), np.isin(a, list(ev)))
        self.stats = {'shell': int(len(W)), 'shell_kept': int(ok.sum()), 'witnesses': dict(
            (str(n), self.wn.count(n)) for n in sorted(set(self.wn)))}
        self.usable_before = len(self.usable)
        self.W = self.W[ok]
        self.cls = self.cls[ok]
        self.usable, self.info, self.masks = sk.usable_classes(self.GK, self.Kv, self.rootkeys4, W10=self.W,
                                                               with_masks=True)
        self.stats['usable_before'] = self.usable_before
        self.stats['usable_after'] = len(self.usable)
        # (kappa, w) = w_V . T / 49 for w in V coordinates
        self.T = self.Gv @ self.Kv.T @ self.D @ self.wx.T
        self.stats['init_s'] = round(time.time() - t0, 1)

    # ---------------------------------------------------------------- skeletons (C1's search + witness o-rules)
    def skeletons(self, max_count=10 ** 6, limit=10 ** 8, use_orbits=True):
        usable = self.usable
        if len(usable) < 6:
            return [], True, 0
        orbits = sk.root_orbits(self.GK, list(usable.keys())) if use_orbits else None
        if orbits is not None:
            U = np.array(sorted(usable, key=lambda x: (orbits[x], x)), dtype=np.int64)
        else:
            U = np.array(sorted(usable), dtype=np.int64)
        n = len(U)
        umask = np.zeros(1 << 15, dtype=bool)
        umask[U] = True
        bits = ((U[:, None] >> np.arange(15)) & 1).astype(np.int64)
        Gm2 = (self.GK % 2).astype(np.int64)
        Fm = (((bits @ Gm2) % 2) << np.arange(15)).sum(axis=1).astype(np.int16)
        P = (np.bitwise_count(Fm[:, None] & U.astype(np.int16)[None, :]) & 1).astype(bool)
        nbr = []
        for i in range(n):
            okr = P[i].copy()
            okr[:i + 1] = False
            okr &= umask[U ^ U[i]]
            nbr.append(int.from_bytes(np.packbits(okr, bitorder='little').tobytes(), 'little'))
        nw = len(self.wn)
        FK = [tuple(int(bin(int(u) & xm).count('1') & 1) for xm in self.xmask) for u in U]
        # feasibility of (o so far, r further elements) per witness norm
        FEAS = {m: {(o, r) for o in range(7) for r in range(7) if any(o <= q <= o + r for q in o_allowed(m))}
                for m in set(self.wn)}
        WN = self.wn
        FU = [tuple(int(((bits[i] @ d) % 2)) for d in self.units) for i in range(n)]
        nu = len(self.units)
        masks = self.masks
        Ul = U.tolist()
        out = []
        nodes = [0]

        def m1_full(G):
            Ga = np.array(G, dtype=np.int64)
            for i in range(6):
                M = masks[G[i]]
                Par = np.bitwise_count(M[:, None] & Ga[None, :]) & 1
                onlyi = (Par[:, i] == 1) & (Par.sum(axis=1) == 1)
                if not (Par.all(axis=1) | onlyi).any():
                    return False
                for i2 in range(i + 1, 6):
                    M = masks[G[i] ^ G[i2]]
                    Par = np.bitwise_count(M[:, None] & Ga[None, :]) & 1
                    one = Par.sum(axis=1) == 1
                    if not (one & ((Par[:, i] == 1) | (Par[:, i2] == 1))).any():
                        return False
            return True

        def rec(chosen, span, cand, o, ou):
            nodes[0] += 1
            if nodes[0] > limit:
                raise TimeoutError
            depth = len(chosen)
            need = 5 - depth
            c = cand
            while c:
                low = c & -c
                i = low.bit_length() - 1
                c ^= low
                x = Ul[i]
                if x in span:
                    continue
                if depth == 0 and orbits is not None and orbits[x] != x:
                    continue
                fk = FK[i]
                o2 = tuple(a + b for a, b in zip(o, fk))
                r = 5 - depth
                bad = False
                for t in range(nw):
                    if (o2[t], r) not in FEAS[WN[t]]:
                        bad = True
                        break
                if bad:
                    continue
                ou2 = tuple(a + b for a, b in zip(ou, FU[i])) if nu else ou
                if depth + 1 == 6:
                    if nu and min(ou2) < 1:
                        continue
                    G = [Ul[t] for t in chosen] + [x]
                    if m1_full(G):
                        out.append(G)
                        if len(out) >= max_count:
                            return True
                    continue
                c2 = cand & nbr[i]
                if c2.bit_count() < need:
                    continue
                if rec(chosen + [i], span | {s_ ^ x for s_ in span}, c2, o2, ou2):
                    return True
            return False

        try:
            stopped = rec([], {0}, (1 << n) - 1, tuple([0] * nw), tuple([0] * nu))
            exhaustive = not stopped
        except TimeoutError:
            exhaustive = False
        return out, exhaustive, nodes[0]

    # ---------------------------------------------------------------- vector level
    def vector_level(self, S, cap=10 ** 7, dyn=0, keep_files=False):
        doms = self.domains(S)
        g = [0] + [int(v) for v in S]
        p = np.array([[bin(gc & xm).count('1') & 1 for gc in g] for xm in self.xmask], dtype=np.int64).reshape(-1, 7)
        o = p.sum(axis=1)
        k = np.minimum(o, 7 - o)
        nwit = len(self.wn)
        for t in range(nwit):
            assert k[t] in K_ALLOWED[self.wn[t]], 'skeleton violates the o-rules'
        newdoms, vals = [], []
        sizes0, sizes1 = [], []
        for ci, (c, d) in enumerate(CELLS):
            Pp, Mn = doms[ci]
            sizes0.append(len(Pp))
            if len(Pp) == 0:
                return {'verdict': 'EXCLUDED', 'why': 'empty %s (C1 domains)' % str((c, d)), 'nodes': 0,
                        'sizes0': sizes0}
            vp = Pp.astype(np.int64) @ self.T
            vm = Mn.astype(np.int64) @ self.T
            assert (vp % 49 == 0).all() and (vm % 49 == 0).all()
            vp //= 49
            vm //= 49
            cut = p[:, c - 1] ^ p[:, d - 1]
            assert ((vp - cut[None, :]) % 2 == 0).all() and ((vm - cut[None, :]) % 2 == 0).all(), 'cut parity'
            keep = np.ones(len(Pp), dtype=bool)
            for t in range(nwit):
                ev, od = label_allowed(self.wn[t], k[t])
                al = list(od if cut[t] else ev)
                keep &= np.isin(np.abs(vp[:, t]), al) & np.isin(np.abs(vm[:, t]), al)
                if self.wn[t] == 6 and k[t] == 2 and not cut[t]:
                    keep &= ~((np.abs(vp[:, t]) == 2) & (np.abs(vm[:, t]) == 2))
            sizes1.append(int(keep.sum()))
            if not keep.any():
                return {'verdict': 'EXCLUDED', 'why': 'empty %s (value filter)' % str((c, d)), 'nodes': 0,
                        'sizes0': sizes0, 'sizes1': sizes1}
            newdoms.append((Pp[keep], Mn[keep]))
            vals.append((vp[keep], vm[keep]))
        # budgets: Parseval per witness, cross sums per witness pair; drop constant ones (checked exactly here)
        cons = []
        for t in range(nwit):
            for t2 in range(t, nwit):
                target = int(4 * (self.wx[t] @ self.D @ self.wx[t2]))
                assert target % 7 == 0
                target //= 7
                contrib = [vp[:, t] * vp[:, t2] + vm[:, t] * vm[:, t2] for vp, vm in vals]
                lo = sum(int(cc.min()) for cc in contrib)
                hi = sum(int(cc.max()) for cc in contrib)
                if lo > target or hi < target:
                    return {'verdict': 'EXCLUDED', 'why': 'budget %d,%d unreachable' % (t, t2), 'nodes': 0,
                            'sizes0': sizes0, 'sizes1': sizes1}
                if lo == hi:
                    continue
                assert all(np.abs(cc).max() <= 127 for cc in contrib)
                cons.append((target, contrib))
        os.makedirs(WORK, exist_ok=True)
        path = os.path.join(WORK, 'wfe_%s_%d_%d.bin' % (self.tag, self.pid, os.getpid()))
        stageB.export_g(path, cap, self.Gv, self.M, newdoms)
        with open(path, 'ab') as f:
            f.write(struct.pack('<i', len(cons)))
            for target, contrib in cons:
                f.write(struct.pack('<i', target))
                for cc in contrib:
                    f.write(np.ascontiguousarray(cc, dtype=np.int8).tobytes())
        res = subprocess.run([EXE, path, str(dyn)], capture_output=True, text=True).stdout.split('\n')
        if not keep_files:
            os.remove(path)
        head = res[0].split('|')[0].split()
        r = {'verdict': head[0], 'nodes': int(head[1]), 'depth': [int(x) for x in head[2:]], 'cons': len(cons),
             'solver_tail': res[0].split('|')[1].strip() if '|' in res[0] else '', 'sizes0': sizes0, 'sizes1': sizes1}
        if head[0] == 'FOUND':
            r['sol'] = [list(map(int, l.split())) for l in res[1:] if l.strip()]
        return r


def run_class(tag, pid, max_skel=10 ** 6, cap=10 ** 7, dyn=0, skel_limit=10 ** 8, log=None):
    t0 = time.time()
    case = Wcase(tag, pid)
    out = {'tag': tag, 'pair': pid, 'V': case.V, **case.stats}
    ts = time.time()
    S, exhaustive, snodes = case.skeletons(max_count=max_skel, limit=skel_limit)
    out.update({'skeletons': len(S), 'skeleton_search_exhaustive': exhaustive, 'skeleton_nodes': snodes,
                'skeleton_s': round(time.time() - ts, 1)})
    kinds, total_nodes, worst = {}, 0, None
    for j, s in enumerate(S):
        r = case.vector_level(s, cap=cap, dyn=dyn)
        key = r['verdict'] + (' ' + r['why'].split()[0] if 'why' in r else '')
        kinds[key] = kinds.get(key, 0) + 1
        total_nodes += r['nodes']
        if worst is None or r['nodes'] > worst[1]:
            worst = (s, r['nodes'])
        if r['verdict'] in ('FOUND', 'LIMIT'):
            out['stop'] = {'skeleton': s, 'result': r}
            break
        if log and j % 100 == 0:
            print('  skel %d/%d %s nodes %d (%.0fs)' % (j, len(S), key, r['nodes'], time.time() - t0), file=log,
                  flush=True)
    out['vector_level'] = kinds
    out['solver_nodes'] = total_nodes
    out['worst'] = worst
    if 'stop' in out:
        out['verdict'] = out['stop']['result']['verdict']
    elif exhaustive:
        out['verdict'] = 'EXCLUDED'
    else:
        out['verdict'] = 'PARTIAL'
    out['s'] = round(time.time() - t0, 1)
    return out


if __name__ == '__main__':
    tag = sys.argv[1]
    ids = [int(x) for x in sys.argv[2].split(',')]
    ms = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 6
    cap = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 7
    dyn = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    sl = int(sys.argv[6]) if len(sys.argv) > 6 else 10 ** 8
    outf = os.path.join(HERE, 'results_%s.jsonl' % tag)
    for pid in ids:
        o = run_class(tag, pid, ms, cap, dyn, sl, log=sys.stderr)
        print(json.dumps(o), flush=True)
        with open(outf, 'a') as f:
            f.write(json.dumps(o) + '\n')
