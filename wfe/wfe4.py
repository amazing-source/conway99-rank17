r"""wfe4.py -- per-class driver for wfe_c.exe (all levels in C; structures, 28/09 night).
Export: filtered shell (per-vector witness filter, basic L1) restricted to usable classes, witness values and norms, the
usable classes in search order (H-orbit key), root classes mod 2 (stabilizer pruning), refined-L1 functionals per vector
(v = (r - w)/2, r a root of V with r = w mod 2V, f_v = (v, K basis) mod 2).
Run: wfe_c.exe est (sampling + Knuth), then wfe_c.exe run on `threads` shards.  Outputs per shard: .limit (skeletons
whose vector search hit the node cap), .found (full frames, to be verified).
Usage: python wfe4.py <tag> <pair>[,<pair>...] [threads] [walks] [max_est_leaves] [est_only 0/1] [nodecap]"""
import json, os, sys, time, struct, subprocess, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wfe
import wfe2
import rank17_skeleton as sk
from rank17_ypairs import short_vectors_chunked

EXE = os.path.join(HERE, 'wfe_c.exe')
WORK = wfe.WORK
PW15 = 1 << np.arange(15)


def l1_masks(case, W):
    """per vector (K coordinates): sorted list of refined-L1 functionals f_v on K/2K"""
    WV = W @ case.Kv
    pwV = 1 << np.arange(case.nB)
    keys = ((WV % 2) * pwV).sum(axis=1)
    KG = case.Kv @ case.Gv
    out = [None] * len(W)
    uk, inv = np.unique(keys, return_inverse=True)
    for j, key in enumerate(uk):
        idx = np.flatnonzero(inv == j)
        R = case.roots_mod2.get(int(key), [])
        if not len(R):
            for i in idx:
                out[i] = []
            continue
        R = np.array(R, dtype=np.int64)
        Vv = (R[None, :, :] - WV[idx][:, None, :])
        assert (Vv % 2 == 0).all()
        Vv //= 2
        F = (Vv @ KG.T) % 2
        M = (F * PW15).sum(axis=2)
        for a, i in enumerate(idx):
            out[i] = sorted(set(int(m) for m in M[a]))
    return out


def export(case, T, path):
    orbits = sk.root_orbits(case.GK, list(case.usable.keys()))
    U = sorted(case.usable, key=lambda x: (orbits[x], x))
    pos = {x: i for i, x in enumerate(U)}
    inU = np.array([int(c) in pos for c in T.cls])
    W = T.W[inU]
    ci = np.array([pos[int(c)] for c in T.cls[inU]], dtype=np.int64)
    order = np.argsort(ci, kind='stable')
    W, ci = W[order], ci[order]
    vals = W @ case.wx.T
    assert np.abs(W).max() <= 127 and np.abs(vals).max() <= 127
    nw = len(case.wn)
    G2 = (case.GK % 2).astype(np.int64)
    L1 = l1_masks(case, W)
    with open(path, 'wb') as f:
        f.write(struct.pack('<4i', len(W), nw, len(U), len(case.units)))
        f.write(np.array(case.GK, dtype='<i4').tobytes())
        f.write(np.array(T.target, dtype='<i4').tobytes())
        f.write(np.array([wfe2.kmask(n) for n in case.wn], dtype=np.int8).tobytes())
        f.write(np.array(case.wn, dtype=np.int8).tobytes())
        rec = np.zeros(len(W), dtype=[('x', 'i1', 15), ('v', 'i1', nw), ('c', '<i4')])
        rec['x'] = W; rec['v'] = vals; rec['c'] = ci
        f.write(rec.tobytes())
        for x in U:
            bits = np.array([(x >> i) & 1 for i in range(15)], dtype=np.int64)
            fm = int((((bits @ G2) % 2) * PW15).sum())
            wp = sum((bin(x & xm).count('1') & 1) << t for t, xm in enumerate(case.xmask))
            up = sum(int((bits @ d) % 2) << u for u, d in enumerate(case.units))
            f.write(struct.pack('<iibQI', x, fm, int(orbits[x] == x), wp, up))
        X, nrm = short_vectors_chunked(case.GK.tolist(), 2)
        R = X[nrm == 2].astype(np.int64)
        rows = sorted(set((int(((r & 1) * PW15).sum()), int((((r @ G2) % 2) * PW15).sum())) for r in R))
        f.write(struct.pack('<i', len(rows)))
        for rb, rf in rows:
            f.write(struct.pack('<ii', rb, rf))
        for lst in L1:
            f.write(struct.pack('<i', len(lst)))
            if lst:
                f.write(np.array(lst, dtype='<i2').tobytes())
    return U


def parse(out):
    d = {}
    for line in out.split('\n'):
        p = line.split()
        if not p:
            continue
        if p[0] == 'PAIRS':
            d['pairs'] = p[2:]
        elif p[0] in ('KNUTH', 'DONE'):
            d[p[0]] = dict((a, float(b)) for a, b in zip(p[1::2], p[2::2]))
    return d


def run(tag, pid, threads=2, walks=20000, max_est=5e7, est_only=False, nodecap=20000000, log=sys.stderr):
    t0 = time.time()
    case = wfe.Wcase(tag, pid)
    out = {'tag': tag, 'pair': pid, **case.stats}
    T = wfe2.Tables(case)
    os.makedirs(WORK, exist_ok=True)
    path = os.path.join(WORK, 'wfec_%s_%d.bin' % (tag, pid))
    export(case, T, path)
    out['prep_s'] = round(time.time() - t0, 1)
    r = subprocess.run([EXE, path, 'est', '1', '0', path, str(walks), '1'], capture_output=True, text=True)
    d = parse(r.stdout)
    out['est'] = d.get('KNUTH')
    out['sampling'] = r.stderr.strip().split('\n')[-1] if r.stderr else ''
    print('[%s:%d] %s | %s | est %s (%.0fs)' % (tag, pid, case.stats, out['sampling'], out['est'], time.time() - t0),
          file=log, flush=True)
    if est_only or out['est'] is None:
        out['verdict'] = 'EST_ONLY' if out['est'] else 'ERROR: ' + r.stdout[-300:] + r.stderr[-300:]
        os.remove(path)
        return out
    if out['est']['est_leaves'] > max_est:
        out['verdict'] = 'TOO_BIG'
        os.remove(path)
        return out
    procs = []
    for sh in range(threads):
        procs.append(subprocess.Popen([EXE, path, 'run', str(threads), str(sh), path, str(walks), '1', '40', str(nodecap)],
                                      stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True))
    tot = collections.Counter()
    limits, founds = [], []
    for sh, p in enumerate(procs):
        so = p.communicate()[0]
        dd = parse(so).get('DONE')
        if dd is None:
            out['verdict'] = 'ERROR shard %d: %s' % (sh, so[-300:])
            return out
        for k, v in dd.items():
            tot[k] += int(v)
        for suf, lst in (('limit', limits), ('found', founds)):
            fn = '%s.%s%d' % (path, suf, sh)
            if os.path.exists(fn):
                lst += [l.rstrip('\n') for l in open(fn) if l.strip()]
                os.remove(fn)
    out['run'] = dict(tot)
    out['limit_skeletons'] = limits[:50]
    out['found'] = founds[:200]
    if founds:
        out['verdict'] = 'FOUND'
    elif limits:
        out['verdict'] = 'LIMIT'
    else:
        out['verdict'] = 'EXCLUDED'
    out['s'] = round(time.time() - t0, 1)
    print('[%s:%d] %s %s (%.0fs)' % (tag, pid, out['verdict'], dict(tot), time.time() - t0), file=log, flush=True)
    os.remove(path)
    return out


if __name__ == '__main__':
    tag = sys.argv[1]
    pids = [int(x) for x in sys.argv[2].split(',')]
    th = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    walks = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
    mest = float(sys.argv[5]) if len(sys.argv) > 5 else 5e7
    eo = len(sys.argv) > 6 and sys.argv[6] == '1'
    cap = int(sys.argv[7]) if len(sys.argv) > 7 else 20000000
    for pid in pids:
        o = run(tag, pid, th, walks, mest, eo, cap)
        print(json.dumps(o), flush=True)
        with open(os.path.join(HERE, 'results4_%s.jsonl' % tag), 'a') as f:
            f.write(json.dumps(o) + '\n')
