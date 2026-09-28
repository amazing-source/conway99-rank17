"""n = 6, k = 3 rung (structures, 28/09): a = +-1 exactly on the 24 orbits of cells (i,j), i in {1,2,3}, j in {4..7}, 0 elsewhere,
with N a = 4a and R^T a = 0 (Sol's J3 pattern; a = W kappa, 7|kappa|^2 = 6). Exact model: orbit_sat.build(7) (lambda/mu exact) +
the 42 eigen-equations as cardinality constraints (encoding as in C1's n5_class.py: N_UV = b_UV - a_UV).
Sign patterns are enumerated from R^T a = 0 and reduced modulo the symmetry group of the setting, computed from vertex
relabellings: block permutations preserving {1,2,3} (S3 x S4), inner flips (j,s) -> (j, 1-s) per block, and a -> -a.
Usage: python n6_k3.py classes            -> list class representatives
       python n6_k3.py solve IDX [secs]    -> exact SAT test of class IDX (CaDiCaL), prints SAT / UNSAT / UNKNOWN"""
import sys, itertools, time, json
sys.path.insert(0, 'C:/home/researcher/Desktop/conway99-involution/code/z2')

m = 7
outer = [(i, j, e) for i in range(1, m + 1) for j in range(i + 1, m + 1) for e in (0, 1)]
N = len(outer)
idx = {o: t for t, o in enumerate(outer)}
P = {1, 2, 3}; Qs = {4, 5, 6, 7}
supp = [U for U, (i, j, e) in enumerate(outer) if (i in P) != (j in P)]
assert len(supp) == 24


def r(U):
    i, j, e = outer[U]
    v = [0] * (m + 1); v[i] = 1; v[j] = 1 if e == 0 else -1
    return v


def act(perm, flip, U):
    """image of orbit U under the vertex relabelling (block perm, inner flips): returns (U', sign) with
    rep(U) -> t^k rep(U'), sign = +1 if it lands on the representative, -1 if on its t-image"""
    i, j, e = outer[U]
    a = (perm[i], 0 ^ flip[perm[i]]); b = (perm[j], e ^ flip[perm[j]])
    (p, s), (q, t) = sorted([a, b])
    # the unordered pair {(p,s),(q,t)} with p<q; representative has first coordinate state 0
    if s == 0:
        return idx[(p, q, t)], 1
    return idx[(p, q, 1 - t)], -1


def classes():
    # enumerate: per block j in 4..7, the six signs of cells (1,j),(2,j),(3,j); constraint at j and at 1,2,3
    per = []
    for sg in itertools.product((1, -1), repeat=6):
        if sum(sg[2 * t] - sg[2 * t + 1] for t in range(3)) == 0:
            per.append(sg)
    vecs = []
    for combo in itertools.product(per, repeat=4):
        a = [0] * N
        for jb, sg in zip(range(4, 8), combo):
            for t, i in enumerate((1, 2, 3)):
                a[idx[(i, jb, 0)]] = sg[2 * t]; a[idx[(i, jb, 1)]] = sg[2 * t + 1]
        ok = True
        for c in range(1, 8):
            if sum(a[U] * r(U)[c] for U in range(N)) != 0:
                ok = False; break
        if ok:
            vecs.append(tuple(a))
    print('valid sign vectors:', len(vecs), flush=True)
    group = []
    for pp in itertools.permutations((1, 2, 3)):
        for qq in itertools.permutations((4, 5, 6, 7)):
            perm = [0] + list(pp) + list(qq)
            for fl in itertools.product((0, 1), repeat=7):
                flip = [0] + list(fl)
                images = [act(perm, flip, U) for U in range(N)]
                group.append(images)
    print('group size (without a -> -a):', len(group), flush=True)
    seen = set(); reps = []
    for v in vecs:
        if v in seen:
            continue
        reps.append(v)
        for images in group:
            w = [0] * N
            for U in range(N):
                U2, sg = images[U]
                w[U2] = sg * v[U]
            seen.add(tuple(w)); seen.add(tuple(-x for x in w))
    return reps


if __name__ == '__main__':
    if sys.argv[1] == 'classes':
        reps = classes()
        print('classes:', len(reps))
        json.dump(reps, open('n6_k3_classes.json', 'w'))
    else:
        from orbit_sat import build
        from pysat.card import CardEnc, EncType
        from pysat.solvers import Solver
        reps = json.load(open('n6_k3_classes.json'))
        k = int(sys.argv[2]); a = reps[k]
        t0 = time.time()
        cnf, pool, outer_, rep, trep, A, B = build(7, 'none')
        assert outer_ == outer
        for U in range(N):
            pos = [B[U, V] for V in range(N) if V != U and a[V] == 1] + [A[U, V] for V in range(N) if V != U and a[V] == -1]
            neg = [A[U, V] for V in range(N) if V != U and a[V] == 1] + [B[U, V] for V in range(N) if V != U and a[V] == -1]
            cnf.extend(CardEnc.equals(lits=pos + [-l for l in neg], bound=4 * a[U] + len(neg), vpool=pool,
                                      encoding=EncType.seqcounter).clauses)
        print('class', k, 'vars', pool.top, 'clauses', len(cnf.clauses), 'build %.1fs' % (time.time() - t0), flush=True)
        with Solver(name='cadical195', bootstrap_with=cnf.clauses) as s:
            t1 = time.time(); res = s.solve()
            print('RESULT n6k3 class', k, {True: 'SAT', False: 'UNSAT'}[res], 'solve %.1fs' % (time.time() - t1), flush=True)
