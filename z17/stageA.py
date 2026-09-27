"""Stage A for the V = Z^17 subcase (structures, 27/09 night): enumerate the binary layer rho(E_7) in F_2^17.

Rows: pi_1 = 0, pi_2..pi_7 in F_2^17 (bitmasks), with w_cd = pi_c + pi_d mod 2 for the 21 cells.
Necessary conditions (RANK17_HANDOFF sections 5-6, Z17_SEARCH_DESIGN.md):
  (a) wt(pi_c + pi_d) in {2,6,10} for all c < d (norm-10 odd supports, q = 1);
  (b) columns non-constant (units: LIFT + min U >= 2); parity energy for every v in Z^17 with Q(v) = 7|v|^2 - Y(v) <= 5,
      where Q(v) >= k(p(v))(7-k)/2 and p(v) = sum of the columns over the odd support of v;
  (c) every pi_c has even intersection with the odd support of y1 and of y2.
Symmetry: rows by partition refinement on coordinate classes +-(y1[i], y2[i]) (each output is one orbit of labelled row tuples
under the within-class permutations), then S_7 (block relabelling) and the class maps of the y-pair symmetries, removed by a canonical
form on the column multisets.  Symmetry breaking used before the canonical form (both valid, since S_7 acts on the labels):
N_12 = min over the 21 cells, and N_13 <= N_14 <= ... <= N_17.
Usage: python stageA.py PAIR_ID [--limit N] [--out FILE]"""
import json, os, sys, time, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, r'C:/home/researcher/Desktop/involution_exploratoire/routeC_check')
from lattice import short_vectors

GOOD = {2, 6, 10}


def load_pair(pid):
    P = json.load(open(os.path.join(os.path.dirname(HERE), 'z17_y_pair_representatives.json')))['pairs']
    p = [q for q in P if q['id'] == pid][0]
    return p['y1'], p['y2']


def classes_of(y1, y2):
    key = {}
    for i in range(17):
        a, b = y1[i], y2[i]
        k = max((a, b), (-a, -b))
        key.setdefault(k, []).append(i)
    return key


def pe_constraints(y1, y2):
    """All (odd-support mask B, kmax) with kmax < 3 from parity energy: Q(v) = 7|v|^2 - (v,y1)^2 - (v,y2)^2 <= 5.
    Q is PSD with kernel span(y1,y2) (primitive); enumerate Z^17 / (Z y1 + Z y2) via an adapted basis."""
    n = 17
    piv = None
    for i in range(n):
        for j in range(n):
            if i != j and abs(y1[i] * y2[j] - y1[j] * y2[i]) == 1:
                piv = (i, j); break
        if piv: break
    assert piv, 'no unimodular 2x2 minor'
    idx = [k for k in range(n) if k not in piv]
    Q = [[7 * (i == j) - y1[i] * y1[j] - y2[i] * y2[j] for j in idx] for i in idx]
    out = {}
    for c in short_vectors(Q, 5):
        c = list(c)
        q = sum(Q[a][b] * c[a] * c[b] for a in range(15) for b in range(15))
        if q <= 0:
            continue
        kmax = 0 if q < 3 else (1 if q < 5 else 2)
        B = 0
        for a, k in enumerate(idx):
            if c[a] % 2:
                B |= 1 << k
        if B == 0:
            continue
        out[B] = min(out.get(B, 3), kmax)
    return out


def kval(col):
    w = bin(col).count('1')
    return min(w, 7 - w)


def run(pid, limit=None, verbose=True, sink=None):
    y1, y2 = load_pair(pid)
    cls = classes_of(y1, y2)
    keys = sorted(cls)
    O1 = sum(1 << i for i in range(17) if y1[i] % 2)
    O2 = sum(1 << i for i in range(17) if y2[i] % 2)
    pe = pe_constraints(y1, y2)
    pe_items = sorted(pe.items(), key=lambda t: t[1])
    t0 = time.time()
    raw = []
    stats = {'leaves': 0, 'pe_fail': 0, 'cover_fail': 0}

    def choose_row(parts, need_min, need_max, prev, wmin, forced=()):
        """Yield (mask, newparts) for the next row: n_P elements from each part (the first n_P, canonical).
        forced: indices of parts that must be taken entirely (never-covered parts at the last row)."""
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
                for q in prev:
                    d = wt + bin(q).count('1') - 2 * bin(m & q).count('1')
                    if d not in GOOD or d < need_min:
                        return
                newparts = []
                for P, s in zip(parts, sel):
                    if s: newparts.append(P[:s])
                    if s < len(P): newparts.append(P[s:])
                yield m, newparts
                return
            rng = (len(parts[a]),) if a in fset else range(len(parts[a]) + 1)
            for s in rng:
                sel[a] = s
                yield from rec(a + 1, wt + s)
            sel[a] = 0

        yield from rec(0, 0)

    def finish(rows):
        cols = []
        for i in range(17):
            c = 0
            for r, m in enumerate(rows):          # rows[0] = pi_1 = 0 -> bit 0 stays 0
                if m >> i & 1:
                    c |= 1 << r
            if c == 0:
                stats['cover_fail'] += 1
                return None
            cols.append(c)
        for B, kmax in pe_items:
            s = 0
            for i in range(17):
                if B >> i & 1:
                    s ^= cols[i]
            if kval(s) > kmax:
                stats['pe_fail'] += 1
                return None
        return cols

    def dfs(rows, parts):
        if limit and len(raw) >= limit:
            return
        c = len(rows)                             # next row index (0-based); rows[0] = pi_1 = 0
        if c == 7:
            stats['leaves'] += 1
            cols = finish(rows)
            if cols is not None:
                if sink is None:
                    raw.append(cols)
                else:
                    sink(cols)
                    stats['kept'] = stats.get('kept', 0) + 1
                    if stats['kept'] % 200000 == 0:
                        print('  kept %d (%.0fs)' % (stats['kept'], time.time() - t0), flush=True)
            return
        w2 = bin(rows[1]).count('1') if c >= 2 else 2
        wmin = w2 if c == 1 else max(w2, bin(rows[-1]).count('1') if c >= 3 else w2)
        forced = ()
        if c == 6:                                # last row: parts not yet covered by any row must be taken entirely
            cov = 0
            for m in rows:
                cov |= m
            forced = tuple(a for a, P in enumerate(parts) if not (cov >> P[0] & 1))
        for m, newparts in choose_row(parts, w2 if c >= 2 else 2, 10, rows[1:], wmin, forced):
            dfs(rows + [m], newparts)

    init = [cls[k] for k in keys]
    dfs([0], init)
    t1 = time.time()
    if verbose:
        print('pair %d: classes %s, PE constraints %d (kmax 0/1/2: %s), raw leaves %d, kept %d, %s, %.1fs'
              % (pid, {k: len(cls[k]) for k in keys}, len(pe), [sum(1 for v in pe.values() if v == t) for t in (0, 1, 2)],
                 stats['leaves'], len(raw), stats, t1 - t0), flush=True)
    return y1, y2, cls, keys, raw


# ---------- canonical form under S_7 and the y-pair class maps ----------
def perm_tables():
    perms = list(itertools.permutations(range(7)))
    T = np.zeros((len(perms), 128), dtype=np.int16)
    for a, p in enumerate(perms):
        for x in range(128):
            y = 0
            for b in range(7):
                if x >> b & 1:
                    y |= 1 << p[b]
            if y & 1:
                y ^= 127
            T[a, x] = y
    return T


def class_maps(y1, y2, cls):
    """Class permutations induced by the allowed y-pair symmetries (sign of y1, sign of y2, swap) that preserve class sizes."""
    keys = sorted(cls)
    maps = []
    for sw in (False, True):
        for s1 in (1, -1):
            for s2 in (1, -1):
                def f(k):
                    a, b = k
                    a, b = (b, a) if sw else (a, b)
                    a, b = s1 * a, s2 * b
                    return max((a, b), (-a, -b))
                img = [f(k) for k in keys]
                if sorted(img) == keys and all(len(cls[k]) == len(cls[f(k)]) for k in keys):
                    m = tuple(keys.index(f(k)) for k in keys)
                    if m not in maps:
                        maps.append(m)
    return maps


def canon(cols, cls, keys, T, maps):
    """cols: list of 17 columns; returns a hashable canonical form (min over S_7 x class maps)."""
    segs = [cls[k] for k in keys]
    arr = np.array(cols, dtype=np.int64)
    mapped = T[:, arr]                             # (5040, 17)
    best = None
    for m in maps:
        # class a of the image receives the columns of class m^-1(a); order the segments by image class
        order = [None] * len(keys)
        for a, b in enumerate(m):
            order[b] = a
        blocks = [np.sort(mapped[:, segs[a]], axis=1) for a in order]
        M = np.concatenate(blocks, axis=1)
        # lexicographic min row
        idx = np.lexsort(M.T[::-1])
        row = tuple(int(x) for x in M[idx[0]])
        if best is None or row < best:
            best = row
    return best


if __name__ == '__main__' and '--raw' in sys.argv:
    # stream raw configs (17 columns as 2-hex-digit fields) to a file; canonicalise later with canon_file.py
    argv = sys.argv[1:]
    pid = int(argv[0])
    out = argv[argv.index('--raw') + 1]
    f = open(out, 'w')
    y1, y2, cls, keys, raw = run(pid, None, True, sink=lambda cols: f.write(''.join('%02x' % c for c in cols) + chr(10)))
    f.close()
    sys.exit(0)

if __name__ == '__main__':
    argv = sys.argv[1:]
    pid = int(argv[0])
    limit = int(argv[argv.index('--limit') + 1]) if '--limit' in argv else None
    out = argv[argv.index('--out') + 1] if '--out' in argv else None
    y1, y2, cls, keys, raw = run(pid, limit)
    T = perm_tables()
    maps = class_maps(y1, y2, cls)
    t0 = time.time()
    seen = {}
    for cols in raw:
        cf = canon(cols, cls, keys, T, maps)
        if cf not in seen:
            seen[cf] = cols
    print('pair %d: %d raw -> %d orbits (class maps %d), canon %.1fs' % (pid, len(raw), len(seen), len(maps), time.time() - t0), flush=True)
    if out:
        json.dump({'pair': pid, 'keys': [list(k) for k in keys], 'classes': {str(k): cls[k] for k in keys},
                   'configs': [list(cf) for cf in seen]}, open(out, 'w'))
