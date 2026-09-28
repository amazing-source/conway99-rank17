"""Probe: number of star-W incidence configurations for the T1 cyclic skeleton (count only, with pruning)."""
import itertools, sys, time
from model import G, vW, sW, starv, sig, star_root, C, STARV
cyc = tuple(tuple((1 if t == k else (-1 if t == (k+1) % 5 else 0)) for t in range(7)) for k in range(5))
def hstar(v, a):
    al, j, be = star_root(v); return al * a[(j, al*be)]
def run(eps):
    a = {(j, s): eps[j] for j in range(6) for s in (1, -1)}
    g = G(cyc)
    PAIRS = list(itertools.combinations(range(5), 2))
    for k in range(5):
        g.add_checked(vW(k), sW((k+1) % 5)); g.add_checked(vW((k+1) % 5), sW(k))
    assert g.full_ok()
    Op = [v for v in STARV if star_root(v)[0] == 1]; Om = [v for v in STARV if star_root(v)[0] == -1]
    sigma_k = -4   # delta=2
    cands = []
    for k in range(5):
        cl = []
        for p in itertools.combinations(Op, 2):
            for m in itertools.combinations(Om, 2):
                if sum(hstar(v, a) for v in p + m) == sigma_k: cl.append(p + m)
        cands.append(cl)
    count = [0, 0]
    wsum = {v: 0 for v in STARV}
    def rec(k):
        if k == 5:
            count[0] += 1
            # star roots: W2-sum must be in -4h + {-2,0,2}
            for v in STARV:
                hv = hstar(v, a)
                if wsum[v] not in (-4*hv - 2, -4*hv, -4*hv + 2): return
            count[1] += 1; return
        for q in cands[k]:
            added = []; ok = True
            for v in q:
                if g.add_checked(vW(k), v): added.append(v)
                else: ok = False; break
            if ok:
                for v in q: wsum[v] += 2; wsum[sig(v)] -= 2
                # prune: a star root with h=+1 can't have positive wsum ... (T1: v's only on -1 roots so fine)
                rec(k+1)
                for v in q: wsum[v] -= 2; wsum[sig(v)] += 2
            for v in added: g.rem(vW(k), v); g.rem(sig(vW(k)), sig(v))
    t0 = time.time(); rec(0)
    return count, time.time() - t0
for eps in ([1,1,1,-1,-1,-1], [1,-1,1,-1,1,-1], [1,1,-1,-1,1,-1]):
    print(eps, run(eps))
