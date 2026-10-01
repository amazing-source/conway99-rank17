import numpy as np, sys, pickle, struct, subprocess, os, time
import evenpipe as ep, skel15
cells = [(c, d) for c in range(7) for d in range(c + 1, 7)]
def firstpos(X):
    nz = X != 0; f = X[np.arange(len(X)), nz.argmax(axis=1)]; return f > 0
def gens_for(L, Ib):
    """reflection vectors preserving every class set: norm-1 vectors of U, and roots whose class lies in I (up to sign)"""
    import skel15 as sk
    I = sk.span(list(Ib))
    g = [v for v in L.S1 if firstpos(v[None, :])[0]]
    rc = ep.cls(L.S2)
    g += [v for v, c in zip(L.S2, rc) if int(c) in I and firstpos(v[None, :])[0]]
    return np.array(g, dtype=np.int8).reshape(-1, 15)
def budgets_for(L, Ib):
    """orthogonal components of L = {u in U : u mod 2 in I^perp}, found from its vectors of norm <= 4 (decomposable ones
    removed).  Returns a list of (dim, B) with B an integer matrix whose rows span the component (coefficient coords).
    Any subspace gives a valid budget (sum over the frame of |P_F beta|^2 = 14 dim F); the components are just good ones."""
    from fractions import Fraction as Fr
    import skel15 as sk
    if not hasattr(L, 'S3'):
        L.S3 = ep.shell(L.G, 3); L.S4 = ep.shell(L.G, 4)
    I = list(Ib)
    par = lambda a, b: np.bitwise_count(a & b) & 1
    T = []
    for S in (L.S1, L.S2, L.S3, L.S4):
        if len(S) == 0: continue
        c = ep.cls(S); ok = np.ones(len(S), dtype=bool)
        for i in I: ok &= par(c, int(L.FM[i])) == 0
        X = S[ok]; X = X[firstpos(X)]
        T.append(X)
    T = np.vstack(T)
    GT = T @ L.G @ T.T
    nrm = np.diag(GT).copy()
    # v decomposable if some x (x != +-v) has (x, v) = |x|^2 or -(x, v) = |x|^2 with v -+ x in T (norm |v|^2 - |x|^2 > 0)
    key = {tuple(v): i for i, v in enumerate(T)}
    dec = np.zeros(len(T), dtype=bool)
    for i in range(len(T)):
        for sgn in (1, -1):
            js = np.flatnonzero((sgn * GT[i] == nrm) & (nrm < nrm[i]))
            for j in js:
                w = T[i] - sgn * T[j]
                if tuple(w) in key or tuple(-w) in key: dec[i] = True; break
            if dec[i]: break
    keep = np.flatnonzero(~dec)
    parent = list(range(len(keep)))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    sub = GT[np.ix_(keep, keep)]
    for a in range(len(keep)):
        for b in np.flatnonzero(sub[a, a + 1:] != 0) + a + 1:
            ra, rb = find(a), find(int(b))
            if ra != rb: parent[rb] = ra
    comps = {}
    for a in range(len(keep)): comps.setdefault(find(a), []).append(keep[a])
    out = []
    for idx in comps.values():
        V = T[idx]
        # row-reduce over Q to a basis
        rows = []
        M = [[Fr(int(x)) for x in v] for v in V]
        basis = []; piv = []
        for v, raw in zip(M, V):
            w = v[:]
            for bvec, p in zip(basis, piv):
                if w[p] != 0:
                    f = w[p] / bvec[p]; w = [a - f * b_ for a, b_ in zip(w, bvec)]
            nzp = next((i for i, a in enumerate(w) if a != 0), None)
            if nzp is not None: basis.append(w); piv.append(nzp); rows.append(raw)
        out.append(np.array(rows, dtype=np.int64))
    return out

def root_components(L, Ib=None):
    """spans of the irreducible components of the root system of L (Ib given) or of U (Ib None), and the norm-1 lines"""
    par = lambda a, b: np.bitwise_count(a & b) & 1
    out = []
    R = L.S2
    if Ib is not None:
        rc = ep.cls(R); ok = np.ones(len(R), dtype=bool)
        for i in Ib: ok &= par(rc, int(L.FM[i])) == 0
        R = R[ok]
    else:
        for e in (L.S1[firstpos(L.S1)] if len(L.S1) else []): out.append(e[None, :].copy())
    R = R[firstpos(R)] if len(R) else R
    if len(R) == 0: return out
    GR = R @ L.G @ R.T
    parent = list(range(len(R)))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a in range(len(R)):
        for b in np.flatnonzero(GR[a, a + 1:] != 0) + a + 1:
            ra, rb = find(a), find(int(b))
            if ra != rb: parent[rb] = ra
    comps = {}
    for a in range(len(R)): comps.setdefault(find(a), []).append(a)
    import sympy
    for idx in comps.values():
        V = sympy.Matrix(R[idx].tolist())
        piv = V.T.rref()[1]
        out.append(R[idx][list(piv)])
    return out

def budget_values(L, B, X):
    """|P_F x|^2 for the rows x of X, F = row span of B; exact, returns integers or None if not integral"""
    from fractions import Fraction as Fr
    import sympy
    M = sympy.Matrix((B @ L.G @ B.T).tolist()); Mi = M.inv()
    b = (X @ L.G @ B.T)
    den = sympy.ilcm(1, 1, *[sympy.fraction(e)[1] for e in Mi])
    Mn = np.array((Mi * den).tolist(), dtype=object)
    vals = ((b.astype(object) @ Mn) * b.astype(object)).sum(axis=1)
    return np.array([int(v) for v in vals], dtype=np.int64), int(den)        # den * |P_F x|^2

def group_for(L, info, Ib, omit):
    """label-permuting symmetries: reflections in the roots of L (class in I^perp) whose point pair avoids `omit`
    (they swap two lines), reflections in norm-1 vectors, and -1.  Returns (list of (M, old-line transposition or None),
    components of the root graph on the 8 points)."""
    par = lambda a, b: np.bitwise_count(a & b) & 1
    I = np.array(info['I'])
    canon = lambda x: int(np.min(I ^ int(x)))
    pairof = {canon(v): k for k, v in info['pair'].items()}
    rc = ep.cls(L.S2); inC = np.ones(len(L.S2), dtype=bool)
    for i in Ib: inC &= par(rc, int(L.FM[i])) == 0
    R = L.S2[inC]; R = R[firstpos(R)]
    pts = [i for i in range(8) if i != omit]
    comp = list(range(8))
    def f(a):
        while comp[a] != a: a = comp[a]
        return a
    out = []
    Id = np.eye(15, dtype=np.int64)
    for r in R:
        pr_ = pairof.get(canon(int(ep.cls(r[None, :])[0])))
        if pr_ is None: continue
        a, b = f(pr_[0]), f(pr_[1])
        if a != b: comp[max(a, b)] = min(a, b)
        if omit in pr_: continue
        M = Id - np.outer(L.G @ r, r)            # norm 2: v -> v - (v, r) r
        out.append((M, (pts.index(pr_[0]), pts.index(pr_[1]))))
    for e in L.S1[firstpos(L.S1)] if len(L.S1) else []:
        out.append((Id - 2 * np.outer(L.G @ e, e), None))
    out.append((-Id, None))
    return out, [f(a) for a in range(8)]

def write_instance(L, scls, cc, path, capm, reorder=True, gens=None, buds=None, grp=None):
    size = {p: int(sum((scls == v).sum() for v in cc[p])) for p in cells}
    if reorder:
        w = [sum(size[p] for p in cells if c in p) for c in range(7)]
        perm = sorted(range(7), key=lambda c: (w[c], c))          # new line i = old line perm[i]
    else:
        perm = list(range(7))
    sizes = []; doms = []
    with open(path, 'wb') as f:
        f.write(struct.pack('<i', 15)); f.write(L.G.astype('<i4').tobytes()); f.write(struct.pack("<ii", capm, 1))
        for (c, d) in cells:
            old = tuple(sorted((perm[c], perm[d])))
            X = np.vstack([L.S5[scls == v] for v in cc[old]]); P = X[firstpos(X)]; X = np.vstack([P, -P])
            sizes.append(len(X)); doms.append(X)
            f.write(struct.pack('<i', len(X))); f.write(X.astype(np.int8).tobytes())
        f.write(struct.pack('<i', 0 if gens is None else len(gens)))
        if gens is not None: f.write(gens.astype(np.int8).tobytes())
        blist = []
        for B in (buds or []):
            if len(B) >= 15: continue
            vd = [budget_values(L, B, X) for X in doms]
            den = vd[0][1]
            if 5 * den > 30000: continue
            blist.append((14 * len(B) * den, [v for v, _ in vd]))
        f.write(struct.pack('<i', len(blist)))
        for tgt, vals in blist:
            f.write(struct.pack('<i', tgt))
            for v in vals: f.write(v.astype('<i2').tobytes())
        inv = {old: new for new, old in enumerate(perm)}
        f.write(struct.pack('<i', 0 if grp is None else len(grp)))
        for M, tr in (grp or []):
            assert np.abs(M).max() <= 127
            f.write(M.astype(np.int8).tobytes())
            lp = list(range(7))
            if tr is not None:
                a, b = inv[tr[0]], inv[tr[1]]; lp[a], lp[b] = b, a
            f.write(bytes(lp))
    return sizes
if __name__ == '__main__':
    nm, gi, capm = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    L = ep.Lat(nm)
    good = pickle.load(open('skel_%s.pkl' % nm, 'rb'))
    Ib, omit, _ = good[gi]
    info = skel15.analyse(L, list(Ib)); cc = skel15.cell_classes(L, info, omit)
    scls = ep.cls(L.S5)
    path = 'one_%s_%d.bin' % (nm, gi)
    G_ = gens_for(L, Ib) if os.environ.get('NOGEN') is None else None
    buds = (budgets_for(L, Ib) + root_components(L, Ib) + root_components(L, None)) if os.environ.get('NOBUD') is None else None
    print('components (dims)', None if buds is None else sorted(len(b) for b in buds))
    grp, comp = group_for(L, info, Ib, omit) if os.environ.get('NOGRP') is None else (None, None)
    print('group generators', None if grp is None else len(grp), 'point components', comp)
    sizes = write_instance(L, scls, cc, path, capm, gens=G_, buds=buds, grp=grp)
    print('generators', 0 if G_ is None else len(G_))
    print('I', Ib, 'omit', omit, 'sizes', sizes)
    t0 = time.time()
    out = subprocess.run([os.environ.get('EXE','./evenframe3.exe'), path] + sys.argv[4:], capture_output=True, text=True).stdout
    print(out.split('\n')[0], round(time.time() - t0, 1), 's')
