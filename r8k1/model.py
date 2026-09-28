"""Vertex model of Y' = {x} u N(x) u star(c) u W2 for R8 k=1 (49 vertices), bitmask adjacency, partial lambda/mu check."""
C = 6
X = 0
def inner(t, s): return 1 + 2*t + (0 if s == 1 else 1)
def starv(al, j, be): return 15 + (0 if al == 1 else 12) + 2*j + (0 if be == 1 else 1)
def vW(k): return 39 + k
def sW(k): return 44 + k
NV = 49
STARV = list(range(15, 39)); WV = list(range(39, 49))
def sig(v):
    if v == 0: return 0
    if v < 15: return v + (1 if (v - 1) % 2 == 0 else -1)
    if v < 39:
        w = v - 15; al = w // 12; r = w % 12; return 15 + (1 - al)*12 + (r ^ 1)
    return v + 5 if v < 44 else v - 5
def star_root(v):
    """(alpha, j, beta) of a star vertex"""
    w = v - 15; al = 1 if w < 12 else -1; w %= 12; return al, w // 2, (1 if w % 2 == 0 else -1)
def root_vec(v, rhos):
    r = [0]*7
    if 15 <= v < 39:
        al, j, be = star_root(v); r[C] = al; r[j] = be
    elif 39 <= v < 44: r = list(rhos[v - 39])
    elif 44 <= v < 49: r = [-x for x in rhos[v - 44]]
    return r
class G:
    def __init__(self, rhos):
        self.rhos = rhos
        self.adj = [0]*NV
        for t in range(7):
            self.add(X, inner(t, 1)); self.add(X, inner(t, -1)); self.add(inner(t, 1), inner(t, -1))
        for v in STARV + WV:
            r = root_vec(v, rhos)
            for t in range(7):
                if r[t] != 0: self.add(v, inner(t, r[t]))
    def add(self, a, b):
        self.adj[a] |= 1 << b; self.adj[b] |= 1 << a
    def rem(self, a, b):
        self.adj[a] &= ~(1 << b); self.adj[b] &= ~(1 << a)
    def has(self, a, b): return (self.adj[a] >> b) & 1
    def pair_ok(self, a, b):
        c = bin(self.adj[a] & self.adj[b]).count('1')
        return c <= (1 if self.has(a, b) else 2)
    def full_ok(self):
        for a in range(NV):
            if bin(self.adj[a]).count('1') > 14: return False
            for b in range(a+1, NV):
                if not self.pair_ok(a, b): return False
        return True
    def add_checked(self, a, b):
        """add a-b and sigma image; return True if all affected pairs ok (edge kept), else revert and return False"""
        sa, sb = sig(a), sig(b)
        if a == sb: return False            # p ~ sigma p forbidden
        if self.has(a, b): return True
        self.add(a, b); self.add(sa, sb)
        if self.local_ok((a, b, sa, sb)): return True
        self.rem(a, b); self.rem(sa, sb); return False
    def local_ok(self, touched):
        for u in set(touched):
            if bin(self.adj[u]).count('1') > 14: return False
        # pairs whose common-neighbour count may have changed: (u, w) for u in touched, w any; plus pairs (w1, w2) with
        # both adjacent to a touched vertex -- they gained a common neighbour only if the new edge made u adjacent to both,
        # i.e. pairs (other endpoint, w) for w ~ u.  Checking all pairs (u, w) and (nbr, nbr') of touched vertices:
        cand = set(touched)
        for u in touched:
            m = self.adj[u]
            while m:
                b = m & -m; cand.add(b.bit_length() - 1); m ^= b
        for u in touched:
            for w in range(NV):
                if w != u and not self.pair_ok(u, w): return False
        # new common neighbour u for pairs (a, w) where a is new neighbour of u: covered above since a in touched
        return True
