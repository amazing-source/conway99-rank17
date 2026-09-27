"""Stage A + B in one pass for a V = Z^17 y-pair (structures, 27/09 night).
Stage A as in stageA.py (rows by partition refinement, S_7 breaking N_12 = min, N_13 <= ... <= N_17), plus an early prune:
every cell class met so far must be USABLE (holds w != +-w' of norm 10 in K, not (3,1), with (w,w') in {0,+-4}).
Each leaf (skeleton) goes to stageB: oriented cell domains (p(a) = eps_c, refined L1), then the exact DFS (Gram, stars, frame).
--sample P processes each leaf with probability P (PASS 1, label PATTERN); without it the run is exhaustive (label COMPUTED once
certificates are attached; here verdicts are exact DFS results, no LP).
Usage: python run_pair.py PAIR [--sample P] [--seed S] [--nodes N] [--out FILE]"""
import json, os, random, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stageA as A
import stageB as B

GOOD = {2, 6, 10}
STATS = []


def main():
    argv = sys.argv[1:]
    pid = int(argv[0])
    samp = float(argv[argv.index('--sample') + 1]) if '--sample' in argv else None
    seed = int(argv[argv.index('--seed') + 1]) if '--seed' in argv else 0
    nlim = int(argv[argv.index('--nodes') + 1]) if '--nodes' in argv else 2 * 10 ** 7
    out = argv[argv.index('--out') + 1] if '--out' in argv else None
    rng = random.Random(seed)
    z = B.Z17(pid)
    y1, y2 = z.y1, z.y2
    cls = A.classes_of(y1, y2)
    keys = sorted(cls)
    O1 = sum(1 << i for i in range(17) if y1[i] % 2)
    O2 = sum(1 << i for i in range(17) if y2[i] % 2)
    pe_items = sorted(A.pe_constraints(y1, y2).items(), key=lambda t: t[1])
    ucache = {}
    kmax_c = [3] * 17
    for Bm, km in pe_items:
        if bin(Bm).count('1') == 1:
            i = Bm.bit_length() - 1
            kmax_c[i] = min(kmax_c[i], km)
    tight = [(Bm, km, bin(Bm).count('1') == 1) for Bm, km in pe_items if km < 3]

    def kfeasible(rows):
        # every PE item (B, kmax): the column sum over B has weight W = #rows with odd |B & row|; final k = min(W, 7-W) must be
        # <= kmax (and >= 1 for a single coordinate: units).  Prune when no final weight in reach works.
        rem = 7 - len(rows)
        for Bm, km, single in tight:
            w = sum(1 for m in rows if bin(m & Bm).count('1') & 1)
            lo = 1 if single else 0
            if not any(lo <= min(W, 7 - W) <= km for W in range(w, w + rem + 1)):
                return False
        return True

    def usable(S):
        r = ucache.get(S)
        if r is None:
            W = z.class_vectors(S)
            if len(W) < 2:
                r = False
            else:
                G = W @ W.T
                r = bool(((G == 0) | (G == 4) | (G == -4)).any())
            ucache[S] = r
        return r

    st = {'leaves': 0, 'processed': 0, 'empty': 0, 'root': 0, 'deeper': 0, 'found': 0, 'limit': 0, 'pe_fail': 0}
    STATS.append(st)
    surv = []
    t0 = time.time()

    def choose_row(parts, need_min, prev, wmin, forced):
        np_ = len(parts)
        sel = [0] * np_
        fset = set(forced)

        def rec(a, wt):
            if wt > 10:
                return
            if a == np_:
                if wt not in GOOD or wt < wmin:
                    return
                m = 0
                for P, s in zip(parts, sel):
                    for x in P[:s]:
                        m |= 1 << x
                if bin(m & O1).count('1') % 2 or bin(m & O2).count('1') % 2:
                    return
                if not usable(m):
                    return
                for q in prev:
                    d = wt + bin(q).count('1') - 2 * bin(m & q).count('1')
                    if d not in GOOD or d < need_min or not usable(m ^ q):
                        return
                newparts = []
                for P, s in zip(parts, sel):
                    if s: newparts.append(P[:s])
                    if s < len(P): newparts.append(P[s:])
                yield m, newparts
                return
            rngs = (len(parts[a]),) if a in fset else range(len(parts[a]) + 1)
            for s in rngs:
                sel[a] = s
                yield from rec(a + 1, wt + s)
            sel[a] = 0

        yield from rec(0, 0)

    def leaf(rows):
        cols = []
        for i in range(17):
            c = 0
            for r, m in enumerate(rows):
                if m >> i & 1:
                    c |= 1 << r
            if c == 0:
                return
            cols.append(c)
        for Bm, kmax in pe_items:
            s = 0
            for i in range(17):
                if Bm >> i & 1:
                    s ^= cols[i]
            if A.kval(s) > kmax:
                st['pe_fail'] += 1
                return
        st['leaves'] += 1
        if samp is not None and rng.random() > samp:
            return
        st['processed'] += 1
        sk = B.Skeleton(z, list(rows)).build()
        if min(sk.sizes()) == 0:
            st['empty'] += 1
        else:
            s2 = B.solve_c(sk, nlim)
            v = s2['verdict']
            if v == 'EXCLUDED':
                st['root' if s2['nodes'] == 1 else 'deeper'] += 1
                if s2['nodes'] > 10 ** 5:
                    st['heavy'] = st.get('heavy', 0) + 1
                    print('  HEAVY-EXCLUDED %d nodes %s' % (s2['nodes'], cols), flush=True)
            else:
                st['found' if v == 'FOUND' else 'limit'] += 1
                surv.append({'cols': cols, 'verdict': v, 'nodes': s2['nodes'], 'sol': s2.get('sol')})
                print('  SURVIVOR', v, cols, flush=True)
        if time.time() - st.get('_last', t0) > 120:
            st['_last'] = time.time()
            print('  %s  %.0fs' % ({k: v for k, v in st.items() if k[0] != '_'}, time.time() - t0), flush=True)

    def dfs(rows, parts):
        c = len(rows)
        if c == 7:
            leaf(rows)
            return
        w2 = bin(rows[1]).count('1') if c >= 2 else 2
        wmin = w2 if c <= 2 else max(w2, bin(rows[-1]).count('1'))
        forced = ()
        if c == 6:
            cov = 0
            for m in rows:
                cov |= m
            forced = tuple(a for a, P in enumerate(parts) if not (cov >> P[0] & 1))
        for m, newparts in choose_row(parts, w2 if c >= 2 else 2, rows[1:], wmin if c >= 2 else 2, forced):
            if tight and not kfeasible(rows + [m]):
                continue
            dfs(rows + [m], newparts)

    dfs([0], [cls[k] for k in keys])
    res = {'pair': pid, 'mode': 'sample %.4f seed %d' % (samp, seed) if samp is not None else 'exhaustive',
           'stats': st, 'survivors': len(surv), 'usable_classes': sum(ucache.values()), 'classes_seen': len(ucache),
           'seconds': round(time.time() - t0, 1)}
    print(json.dumps(res), flush=True)
    if out:
        json.dump({'result': res, 'survivors': surv}, open(out, 'w'))


if __name__ == '__main__':
    main()
