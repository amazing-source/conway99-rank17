"""Step C: W-W adjacency states per skeleton.  state[k,l] in {0 none, 1 aligned (v_k~sv_l), 2 anti (v_k~v_l), 3 D (both)}"""
import itertools, sys
from skel import shapes, signed_skeletons
from model import G, vW, sW
PAIRS = list(itertools.combinations(range(5), 2))
def ww_states(rhos):
    g = G(rhos)
    if not g.full_ok(): return []
    out = []
    st = {}
    def rec(i):
        if i == len(PAIRS):
            # delta_k in 2..6
            for k in range(5):
                al = sum(1 for p in PAIRS if k in p and st[p] == 1)
                an = sum(1 for p in PAIRS if k in p and st[p] == 2)
                if not (2 <= al - an <= 6): return
            out.append(dict(st)); return
        k, l = PAIRS[i]
        for s in (0, 1, 2, 3):
            edges = []
            if s in (1, 3): edges += [(vW(k), sW(l)), (vW(l), sW(k))]
            if s in (2, 3): edges += [(vW(k), vW(l))]
            added = []; ok = True
            for a, b in edges:
                if g.has(a, b): continue
                if g.add_checked(a, b): added.append((a, b))
                else: ok = False; break
            if ok:
                # prune: remaining pairs can still give delta_k >= 2?
                st[(k, l)] = s
                feas = True
                for m in range(5):
                    al = sum(1 for p in PAIRS[:i+1] if m in p and st[p] == 1)
                    an = sum(1 for p in PAIRS[:i+1] if m in p and st[p] == 2)
                    rem = sum(1 for p in PAIRS[i+1:] if m in p)
                    if al + rem - an < 2: feas = False
                if feas: rec(i+1)
                del st[(k, l)]
            from model import sig
            for a, b in added: g.rem(a, b); g.rem(sig(a), sig(b))
    rec(0)
    return out
if __name__ == '__main__':
    tot = 0; per = []
    sk = [s for sh in shapes() for s in signed_skeletons(sh)]
    for rhos in sk:
        w = ww_states(rhos); tot += len(w); per.append(len(w))
    print('skeletons', len(sk), 'W-W states total', tot, 'max per skeleton', max(per), 'skeletons with >0:', sum(1 for p in per if p))
