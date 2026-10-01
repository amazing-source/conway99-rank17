"""CHECK + construction attempt (small, exact): the rootless class L = Lambda(Hamming) in the model F2^4 = F2 + F8.
1. Verifies the hand claims: 7 forms sigma_x with zero sum and non-degenerate pairwise sums, 6 ovoids per form, the
   intersection pattern of ovoids.
2. Looks for a frame invariant under the element g of order 7 (multiplication by zeta on F8)."""
import itertools, numpy as np, collections, sys
# F8 = F2[z]/(z^3+z+1), elements as 3-bit ints
def mul(a, b):
    r = 0
    for i in range(3):
        if (b >> i) & 1: r ^= a << i
    for i in (4, 3):
        if (r >> i) & 1: r ^= 0b1011 << (i - 3)
    return r
def pw(a, n):
    r = 1
    for _ in range(n % 7 if a else n): r = mul(r, a)
    return r
def tr(a): return (a ^ mul(a, a) ^ mul(mul(a, a), mul(a, a))) & 1 if True else 0
def trace(a):
    s = a ^ pw(a, 2) ^ pw(a, 4)
    assert s in (0, 1); return s
Z = 2                                                   # zeta
V = [(t, u) for t in (0, 1) for u in range(8) if (t, u) != (0, 0)]
idx = {v: i for i, v in enumerate(V)}
def sigma(x):
    a = pw(Z, (-x) % 7); b = pw(Z, (-3 * x) % 7)
    def f(p, q):
        (t, u), (t2, u2) = p, q
        c = trace(mul(a, (u2 if t else 0) ^ (u if t2 else 0)))
        d = trace(mul(b, mul(u, pw(u2, 2)) ^ mul(pw(u, 2), u2)))
        return c ^ d
    return np.array([[f(p, q) for q in V] for p in V], dtype=np.int64)
S = [sigma(x) for x in range(7)]
add = lambda p, q: (p[0] ^ q[0], p[1] ^ q[1])
def nondeg(W):       # no point orthogonal to all points
    return not any((W[i] == 0).all() for i in range(15))
def bilinear(W):
    return all(W[idx[add(V[i], V[j])]][k] == (W[i][k] ^ W[j][k]) for i in range(15) for j in range(15) if i != j for k in range(15))
assert all(bilinear(s) for s in S) and all(nondeg(s) for s in S)
assert (sum(S) % 2 == 0).all()
assert all(nondeg((S[x] + S[y]) % 2) for x in range(7) for y in range(x + 1, 7))
assert not any(nondeg((S[x] + S[y] + S[z]) % 2) for x, y, z in itertools.combinations(range(7), 3))
print('CHECK forms: 7 forms, zero sum, pairwise sums non-degenerate, triple sums degenerate: OK')
def ovoids(W):
    out = []
    for c in itertools.combinations(range(15), 5):
        if all(W[i][j] for i, j in itertools.combinations(c, 2)): out.append(frozenset(c))
    return out
cells = [(x, y) for x in range(7) for y in range(x + 1, 7)]
OV = {c: ovoids((S[c[0]] + S[c[1]]) % 2) for c in cells}
assert all(len(v) == 6 for v in OV.values())
assert all(len(a & b) == 1 for v in OV.values() for a, b in itertools.combinations(v, 2))
pat = collections.Counter()
for c1, c2 in itertools.combinations(cells, 2):
    s = len(set(c1) & set(c2))
    for a in OV[c1]:
        for b in OV[c2]: pat[(s, len(a & b))] += 1
print('CHECK ovoids: 6 per cell form, same form meet in 1; (overlap, |meet|) counts:', dict(pat))
# pentagon claim: for O and an adjacent form, 1 ovoid disjoint, 5 meeting in 2 forming a 5-cycle on O
okp = True
for c1 in cells[:3]:
    for O in OV[c1]:
        for c2 in cells:
            if len(set(c1) & set(c2)) != 1: continue
            m = [len(O & b) for b in OV[c2]]
            if sorted(m) != [0, 2, 2, 2, 2, 2]: okp = False
            deg = collections.Counter(p for b in OV[c2] for p in O & b)
            if any(deg[p] != 2 for p in O): okp = False
print('CHECK pentagons:', okp)
# ---- Z7 action
def g(v): return (v[0], mul(Z, v[1]))
perm = [idx[g(v)] for v in V]                       # point i -> perm[i]
def rot(vec, k):
    out = vec
    for _ in range(k):
        n = np.zeros(15, dtype=np.int64); n[perm] = out; out = n
    return out
def cell_rot(c, k): return tuple(sorted(((c[0] + k) % 7, (c[1] + k) % 7)))
# g maps sigma_x to sigma_{x+1}? check: ovoid images land in the rotated cell
for c in cells:
    for O in OV[c]:
        img = frozenset(perm[i] for i in O)
        assert img in OV[cell_rot(c, 1)], 'rotation does not match the cell labelling'
print('CHECK symmetry: g maps the ovoids of cell {x,y} to those of {x+1,y+1}: OK')
base = [(0, 1), (0, 2), (0, 3)]
cand = []                                           # (base index, vector)
for bi, c in enumerate(base):
    for O in OV[c]:
        O = sorted(O)
        for signs in itertools.product((1, -1), repeat=5):
            v = np.zeros(15, dtype=np.int64); v[O] = signs
            cand.append((bi, v))
allowed = {0: (1, -1, -3), 1: (0, -2, -4), 2: (-1, -3, -5)}
def ok_pair(i, j):
    bi, v = cand[i]; bj, w = cand[j]
    for k in range(7):
        c2 = cell_rot(base[bj], k)
        if i == j and k == 0: continue
        s = len(set(base[bi]) & set(c2))
        if s == 2 and (v == rot(w, k)).all(): return False        # same vector twice
        if int(v @ rot(w, k)) not in allowed[s]: return False
    return True
n = len(cand)
ROT = np.array([[rot(v, k) for k in range(7)] for _, v in cand])        # (n, 7, 15)
self_ok = [ok_pair(i, i) for i in range(n)]
print('candidates', n, 'self-compatible', sum(self_ok))
good = [i for i in range(n) if self_ok[i]]
comp = {}
for a in good:
    for b in good:
        if a < b: comp[(a, b)] = ok_pair(a, b) and ok_pair(b, a)
byb = {bi: [i for i in good if cand[i][0] == bi] for bi in range(3)}
def cmp_(a, b): return comp[(a, b)] if a < b else comp[(b, a)]
sols = 0; tested = 0
for a1, a2 in itertools.combinations(byb[0], 2):
    if not cmp_(a1, a2): continue
    for b1, b2 in itertools.combinations(byb[1], 2):
        if not (cmp_(b1, b2) and cmp_(a1, b1) and cmp_(a1, b2) and cmp_(a2, b1) and cmp_(a2, b2)): continue
        for c1, c2 in itertools.combinations(byb[2], 2):
            six = [a1, a2, b1, b2, c1, c2]
            if not all(cmp_(p, q) for p, q in itertools.combinations(six, 2)): continue
            tested += 1
            F = np.vstack([ROT[i] for i in six])                # 42 x 15
            if (F.T @ F == 14 * np.eye(15, dtype=np.int64)).all():
                sols += 1
                if sols <= 2: print('FRAME FOUND', six); np.save('z7_frame_%d.npy' % sols, F)
print('pairwise-compatible sextuples', tested, 'tight frames', sols)
