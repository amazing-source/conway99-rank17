r"""wfe3.py -- per-class driver for skel_w2.exe (structures, 28/09 night).
  1. Wcase (wfe.py): witnesses, filtered shell, usable classes (+ M1 masks), H-orbits.
  2. Export the raw data (filtered shell minus basic-L1 vectors, values, classes in search order) for skel_w2.exe.
  3. skel_w2 est: sampling (killing witness pairs) + Knuth estimate with the chosen tables.
  4. skel_w2 run (sharded): survivors -> exact all-pair leaf check (wfe2.leaf_eval) -> vector level (wfe.vector_level).
Usage: python wfe3.py <tag> <pair>[,<pair>...] [threads] [walks] [max_est_leaves] [est_only 0/1]"""
import json, os, sys, time, struct, subprocess, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wfe
import wfe2
import rank17_skeleton as sk

EXE = os.path.join(HERE, os.environ.get('SKELEXE', 'skel_w3.exe'))
WORK = wfe.WORK


def export_raw(case, T, path):
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
    pw = 1 << np.arange(15)
    with open(path, 'wb') as f:
        f.write(struct.pack('<4i', len(W), nw, len(U), len(case.units)))
        f.write(np.array(case.GK, dtype='<i4').tobytes())
        f.write(np.array(T.target, dtype='<i4').tobytes())
        f.write(np.array([wfe2.kmask(n) for n in case.wn], dtype=np.int8).tobytes())
        rec = np.zeros(len(W), dtype=[('x', 'i1', 15), ('v', 'i1', nw), ('c', '<i4')])
        rec['x'] = W; rec['v'] = vals; rec['c'] = ci
        f.write(rec.tobytes())
        for x in U:
            bits = np.array([(x >> i) & 1 for i in range(15)], dtype=np.int64)
            fm = int((((bits @ G2) % 2) * pw).sum())
            wp = sum((bin(x & xm).count('1') & 1) << t for t, xm in enumerate(case.xmask))
            up = sum(int((bits @ d) % 2) << u for u, d in enumerate(case.units))
            f.write(struct.pack('<iibQI', x, fm, int(orbits[x] == x), wp, up))
        # root classes mod 2 for the stabilizer pruning (skel_w3)
        from rank17_ypairs import short_vectors_chunked
        X, nrm = short_vectors_chunked(case.GK.tolist(), 2)
        R = X[nrm == 2].astype(np.int64)
        rows = sorted(set((int(((r & 1) << np.arange(15)).sum()), int((((r @ G2) % 2) * pw).sum())) for r in R))
        f.write(struct.pack('<i', len(rows)))
        for rb, rf in rows:
            f.write(struct.pack('<ii', rb, rf))
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


def run(tag, pid, threads=2, walks=3000, max_est=2e9, est_only=False, log=sys.stderr):
    t0 = time.time()
    case = wfe.Wcase(tag, pid)
    out = {'tag': tag, 'pair': pid, **case.stats}
    T = wfe2.Tables(case)
    os.makedirs(WORK, exist_ok=True)
    path = os.path.join(WORK, 'skw2_%s_%d.bin' % (tag, pid))
    export_raw(case, T, path)
    out['prep_s'] = round(time.time() - t0, 1)
    r = subprocess.run([EXE, path, 'est', '1', '0', '-', str(walks), '1'], capture_output=True, text=True)
    d = parse(r.stdout)
    out['est'] = d.get('KNUTH')
    out['pairs'] = d.get('pairs')
    out['sampling'] = r.stderr.strip().split('\n')[-1] if r.stderr else ''
    print('[%s:%d] %s | %s | est %s (%.0fs)' % (tag, pid, case.stats, out['sampling'], out['est'], time.time() - t0),
          file=log, flush=True)
    if est_only or out['est'] is None:
        out['verdict'] = 'EST_ONLY' if out['est'] else 'ERROR'
        os.remove(path)
        return out
    if out['est']['est_leaves'] > max_est:
        out['verdict'] = 'TOO_BIG'
        os.remove(path)
        return out
    procs = []
    for sh in range(threads):
        sp = path + '.surv%d' % sh
        procs.append((sp, subprocess.Popen([EXE, path, 'run', str(threads), str(sh), sp, str(walks), '1'],
                                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)))
    tot = collections.Counter()
    survivors = []
    for sp, p in procs:
        dd = parse(p.communicate()[0]).get('DONE', {})
        for k, v in dd.items():
            tot[k] += int(v)
        if os.path.exists(sp):
            survivors += [list(map(int, l.split())) for l in open(sp) if l.strip()]
            os.remove(sp)
    out['skel'] = dict(tot)
    print('[%s:%d] skel %s (%.0fs)' % (tag, pid, dict(tot), time.time() - t0), file=log, flush=True)
    kinds = collections.Counter()
    allp = list(range(len(T.pairs)))
    for s in survivors:
        v, info = wfe2.leaf_eval(T, s, allp)
        if v != 'PASS':
            kinds['leaf ' + v] += 1
            continue
        rr = case.vector_level(s)
        key = rr['verdict'] + (' ' + rr['why'].split()[0] if 'why' in rr else '')
        kinds[key] += 1
        if rr['verdict'] in ('FOUND', 'LIMIT'):
            out['stop'] = {'skeleton': s, 'result': {k: v for k, v in rr.items() if k != 'depth'}}
            break
    out['survivor_verdicts'] = dict(kinds)
    out['verdict'] = out['stop']['result']['verdict'] if 'stop' in out else 'EXCLUDED'
    out['s'] = round(time.time() - t0, 1)
    os.remove(path)
    return out


if __name__ == '__main__':
    tag = sys.argv[1]
    pids = [int(x) for x in sys.argv[2].split(',')]
    th = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    walks = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
    mest = float(sys.argv[5]) if len(sys.argv) > 5 else 2e9
    eo = len(sys.argv) > 6 and sys.argv[6] == '1'
    for pid in pids:
        o = run(tag, pid, th, walks, mest, eo)
        print(json.dumps(o), flush=True)
        with open(os.path.join(HERE, 'results3_%s.jsonl' % tag), 'a') as f:
            f.write(json.dumps(o) + '\n')
