r"""wfe2.py -- per-class driver: witness budgets at the skeleton leaves (skel_w.exe), survivors to the vector level (wfe.py).
Steps for a class (tag, pair):
  1. Wcase (wfe.py): witnesses (7|kappa|^2 <= 8), filtered shell, usable classes + M1 masks.
  2. Pair selection: random skeletons (random DFS orders) are checked with M1-aware tables over ALL witness pairs; the
     pairs that kill are ranked, and a covering set P (cover twice + most frequent, <= PMAX) is exported.
  3. Tables for P on every usable class (valid pairs = filtered shell, basic L1, |(w,w')| in {0,4}), int8.
  4. skel_w.exe: Knuth estimate, then the exhaustive run (sharded); survivors (M1 + budgets over P passed) go to the full
     vector level (C1's domains + exact-k value filter + all budgets + sbdfs_w.exe).
Usage: python wfe2.py <tag> <pair> [threads] [knuth_walks] [max_est_leaves]"""
import json, os, sys, time, struct, subprocess, collections, random
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wfe
import rank17_skeleton as sk

SKEL = os.path.join(HERE, 'skel_w.exe')
WORK = wfe.WORK
PMAX = 40


def kmask(n):
    return sum(1 << k for k in wfe.K_ALLOWED[n])


class Tables:
    """per usable class: valid unordered pairs, their two a-masks, witness values (all witnesses)"""

    def __init__(self, case):
        self.case = case
        W = case.W.astype(np.int64)
        sh = 2 * np.arange(case.Kv.shape[1])
        keep = np.ones(len(W), dtype=bool)
        for s0 in range(0, len(W), 100000):
            k = (((W[s0:s0 + 100000] @ case.Kv) % 4) << sh).sum(axis=1)
            keep[s0:s0 + 100000] = ~np.isin(k, case.rootkeys4)
        W = W[keep]
        cls = case.cls[keep]
        order = np.argsort(cls, kind='stable')
        self.W, self.cls = W[order], cls[order]
        b = np.flatnonzero(np.diff(self.cls)) + 1
        st = np.concatenate([[0], b]); en = np.concatenate([b, [len(self.cls)]])
        self.rng = {int(self.cls[s]): (s, e) for s, e in zip(st, en)}
        self.G2 = (case.GK % 2).astype(np.int64)
        self.pw = 1 << np.arange(15)
        self.cache = {}
        nw = len(case.wn)
        self.pairs = [(t, t2) for t in range(nw) for t2 in range(t, nw)]
        self.PT = np.array(self.pairs, dtype=np.int64).reshape(-1, 2)
        self.target = np.array([int(4 * (case.wx[t] @ case.D @ case.wx[t2])) // 7 for t, t2 in self.pairs])
        assert all(int(4 * (case.wx[t] @ case.D @ case.wx[t2])) % 7 == 0 for t, t2 in self.pairs)

    def raw(self, C):
        """(masks per pair row (2 per pair), values of w, values of w') for class C"""
        if C in self.cache:
            return self.cache[C]
        s, e = self.rng.get(C, (0, 0))
        Wc = self.W[s:e]
        ip = Wc @ self.case.GK.astype(np.int64) @ Wc.T
        n = len(Wc)
        I, J = np.nonzero(np.isin(np.abs(ip), (0, 4)) & (np.arange(n)[:, None] < np.arange(n)[None, :]))
        vals = Wc @ self.case.wx.T
        m1 = ((((Wc[I] + Wc[J]) // 2) @ self.G2) % 2) @ self.pw
        m2 = ((((Wc[I] - Wc[J]) // 2) @ self.G2) % 2) @ self.pw
        r = (m1, m2, vals[I], vals[J])
        if len(self.cache) > 4000:
            self.cache.pop(next(iter(self.cache)))
        self.cache[C] = r
        return r

    def table(self, C, pidx):
        """(unique masks, min (nm, |pidx|), max) over the valid pairs of class C for the witness pairs pidx"""
        m1, m2, A, B = self.raw(C)
        if len(m1) == 0:
            return np.zeros(0, dtype=np.int64), None, None
        PT = self.PT[pidx]
        con = A[:, PT[:, 0]] * A[:, PT[:, 1]] + B[:, PT[:, 0]] * B[:, PT[:, 1]]
        M = np.concatenate([m1, m2]); CON = np.concatenate([con, con])
        um, inv = np.unique(M, return_inverse=True)
        mn = np.full((len(um), len(pidx)), 10 ** 6, dtype=np.int64)
        mx = np.full((len(um), len(pidx)), -10 ** 6, dtype=np.int64)
        np.minimum.at(mn, inv, CON)
        np.maximum.at(mx, inv, CON)
        return um, mn, mx


def leaf_eval(T, s, pidx):
    """M1-aware budget check of a full skeleton s (6 classes, search order) over witness pairs pidx:
    returns ('M1', cell) or ('BUDGET', killing pair indices) or ('PASS', [])"""
    g = [0] + [int(x) for x in s]
    G = np.array(s, dtype=np.int64)
    lo = np.zeros(len(pidx), dtype=np.int64); hi = np.zeros(len(pidx), dtype=np.int64)
    for c in range(7):
        for d in range(c + 1, 7):
            um, mn, mx = T.table(g[c] ^ g[d], pidx)
            if len(um) == 0:
                return 'M1', (c, d)
            par = np.bitwise_count(um[:, None] & G[None, :]) & 1
            eps = np.ones(6, dtype=np.int64) if c == 0 else (np.arange(1, 7) == c).astype(np.int64)
            ok = (par == eps).all(axis=1)
            if not ok.any():
                return 'M1', (c, d)
            lo += mn[ok].min(axis=0); hi += mx[ok].max(axis=0)
    tg = T.target[pidx]
    bad = np.flatnonzero((lo > tg) | (hi < tg))
    return ('BUDGET', [pidx[b] for b in bad]) if len(bad) else ('PASS', [])


def sample_skeletons(case, nsamp, seed=0, limit=20000):
    out = []
    rs = np.random.RandomState(seed)
    keys = sorted(case.usable)
    tries = 0
    while len(out) < nsamp and tries < 20 * nsamp:
        tries += 1
        perm = rs.permutation(len(keys))
        S, ex, nn = skeletons_random(case, [keys[i] for i in perm], limit)
        if S:
            out.append(sorted(S))
    return out


def skeletons_random(case, order, limit):
    """first skeleton found with the given candidate order (no H-reduction); None if none within limit"""
    usable = case.usable
    U = np.array(order, dtype=np.int64)
    n = len(U)
    umask = np.zeros(1 << 15, dtype=bool); umask[U] = True
    bits = ((U[:, None] >> np.arange(15)) & 1).astype(np.int64)
    Gm2 = (case.GK % 2).astype(np.int64)
    Fm = (((bits @ Gm2) % 2) << np.arange(15)).sum(axis=1)
    FK = [tuple(int(bin(int(u) & xm).count('1') & 1) for xm in case.xmask) for u in U]
    FEAS = {m: {(o, r) for o in range(7) for r in range(7) if any(o <= q <= o + r for q in wfe.o_allowed(m))}
            for m in set(case.wn)}
    FU = [tuple(int(((bits[i] @ d) % 2)) for d in case.units) for i in range(n)]
    nodes = [0]
    Ul = U.tolist()
    odd = lambda i, j: bin(int(Fm[i]) & Ul[j]).count('1') & 1

    def rec(ch, span, o, ou):
        nodes[0] += 1
        if nodes[0] > limit:
            return None
        depth = len(ch)
        # candidates: any index not chosen, odd with all chosen, sums usable
        for i in range(n):
            x = Ul[i]
            if x in span or i in ch:
                continue
            if any(not odd(j, i) or not umask[Ul[j] ^ x] for j in ch):
                continue
            o2 = tuple(a + b for a, b in zip(o, FK[i]))
            if any((o2[t], 5 - depth) not in FEAS[case.wn[t]] for t in range(len(o2))):
                continue
            ou2 = tuple(a + b for a, b in zip(ou, FU[i]))
            if depth == 5:
                if ou2 and min(ou2) < 1:
                    continue
                G = [Ul[j] for j in ch] + [x]
                if m1_ok(case, G):
                    return G
                continue
            r = rec(ch + [i], span | {s_ ^ x for s_ in span}, o2, ou2)
            if r:
                return r
            if nodes[0] > limit:
                return None
        return None
    r = rec([], {0}, tuple([0] * len(case.wn)), tuple([0] * len(case.units)))
    return r, False, nodes[0]


def m1_ok(case, G):
    Ga = np.array(G, dtype=np.int64)
    masks = case.masks
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


def choose_pairs(T, samples, pmax=PMAX):
    allp = list(range(len(T.pairs)))
    kills = []
    stats = collections.Counter()
    for s in samples:
        v, info = leaf_eval(T, s, allp)
        stats[v] += 1
        if v == 'BUDGET':
            kills.append(set(info))
    freq = collections.Counter(p for k in kills for p in k)
    chosen = []
    for rnd in range(2):                                   # cover every killed sample twice when possible
        unc = [k for k in kills if len(k & set(chosen)) <= rnd]
        while unc:
            cnt = collections.Counter(p for k in unc for p in k if p not in chosen)
            if not cnt:
                break
            p, _ = cnt.most_common(1)[0]
            chosen.append(p)
            unc = [k for k in unc if len(k & set(chosen)) <= rnd]
    for p, _ in freq.most_common():
        if len(chosen) >= pmax:
            break
        if p not in chosen:
            chosen.append(p)
    return chosen[:max(pmax, 1)], stats, freq


def export(case, T, pidx, path):
    usable = case.usable
    orbits = sk.root_orbits(case.GK, list(usable.keys()))
    U = sorted(usable, key=lambda x: (orbits[x], x))
    G2 = T.G2
    pw = T.pw
    with open(path, 'wb') as f:
        f.write(struct.pack('<4i', len(U), len(case.wn), len(case.units), len(pidx)))
        f.write(np.array(T.target[pidx], dtype='<i4').tobytes())
        f.write(np.array([kmask(n) for n in case.wn], dtype=np.int8).tobytes())
        for x in U:
            bits = np.array([(x >> i) & 1 for i in range(15)], dtype=np.int64)
            fm = int((((bits @ G2) % 2) * pw).sum())
            wp = sum((bin(x & xm).count('1') & 1) << t for t, xm in enumerate(case.xmask))
            up = sum(int((bits @ d) % 2) << u for u, d in enumerate(case.units))
            um, mn, mx = T.table(x, pidx)
            f.write(struct.pack('<iibQIi', x, fm, int(orbits[x] == x), wp, up, len(um)))
            if len(um):
                assert np.abs(mn).max() <= 127 and np.abs(mx).max() <= 127
                f.write(np.array(um, dtype='<i2').tobytes())
                f.write(np.ascontiguousarray(mn, dtype=np.int8).tobytes())
                f.write(np.ascontiguousarray(mx, dtype=np.int8).tobytes())
    return U


def run(tag, pid, threads=2, walks=2000, max_est=5e9, nsamp=300, log=sys.stderr):
    t0 = time.time()
    case = wfe.Wcase(tag, pid)
    out = {'tag': tag, 'pair': pid, **case.stats}
    print('[%s:%d] init %.0fs %s' % (tag, pid, time.time() - t0, case.stats), file=log, flush=True)
    T = Tables(case)
    samples = sample_skeletons(case, nsamp)
    pidx, st, freq = choose_pairs(T, samples)
    out['samples'] = dict(st)
    out['pairs'] = [T.pairs[p] for p in pidx]
    print('[%s:%d] samples %s, %d pairs chosen (%.0fs)' % (tag, pid, dict(st), len(pidx), time.time() - t0), file=log,
          flush=True)
    os.makedirs(WORK, exist_ok=True)
    path = os.path.join(WORK, 'skw_%s_%d.bin' % (tag, pid))
    export(case, T, pidx, path)
    out['export_s'] = round(time.time() - t0, 1)
    est = subprocess.run([SKEL, path, '-', '1', '0', str(walks), '1'], capture_output=True, text=True).stdout.split()
    e = dict(zip(est[1::2], est[2::2]))
    out['knuth'] = e
    print('[%s:%d] knuth %s (%.0fs)' % (tag, pid, e, time.time() - t0), file=log, flush=True)
    if float(e.get('est_leaves', 0)) > max_est:
        out['verdict'] = 'TOO_BIG'
        return out
    procs = []
    for sh in range(threads):
        sp = path + '.surv%d' % sh
        procs.append((sp, subprocess.Popen([SKEL, path, sp, str(threads), str(sh)], stdout=subprocess.PIPE, text=True)))
    tot = collections.Counter()
    survivors = []
    for sp, p in procs:
        o = p.communicate()[0].split()
        d = dict(zip(o[1::2], o[2::2]))
        for k, v in d.items():
            tot[k] += int(v)
        survivors += [list(map(int, l.split())) for l in open(sp) if l.strip()]
        os.remove(sp)
    out['skel_w'] = dict(tot)
    print('[%s:%d] skel_w %s (%.0fs)' % (tag, pid, dict(tot), time.time() - t0), file=log, flush=True)
    kinds = collections.Counter()
    allp = list(range(len(T.pairs)))
    for s in survivors:
        v, info = leaf_eval(T, s, allp)
        if v != 'PASS':
            kinds['leaf ' + v] += 1
            continue
        r = case.vector_level(s)
        key = r['verdict'] + (' ' + r['why'].split()[0] if 'why' in r else '')
        kinds[key] += 1
        if r['verdict'] in ('FOUND', 'LIMIT'):
            out['stop'] = {'skeleton': s, 'result': {k: v for k, v in r.items() if k != 'depth'}}
            break
    out['survivor_verdicts'] = dict(kinds)
    out['verdict'] = out['stop']['result']['verdict'] if 'stop' in out else 'EXCLUDED'
    out['s'] = round(time.time() - t0, 1)
    os.remove(path)
    return out


if __name__ == '__main__':
    tag, pid = sys.argv[1], int(sys.argv[2])
    th = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    walks = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
    mest = float(sys.argv[5]) if len(sys.argv) > 5 else 5e9
    o = run(tag, pid, th, walks, mest)
    print(json.dumps(o), flush=True)
    with open(os.path.join(HERE, 'results2_%s.jsonl' % tag), 'a') as f:
        f.write(json.dumps(o) + '\n')
