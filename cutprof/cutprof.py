"""Cut profiling (structures, 28/09 late).  For an open class (V, y-pair) with lattice K (rank 15, det 49):
 - witnesses: kappa in K* with n = 7|kappa|^2 <= 8, as z with z^T H z = n (H = 7 G^-1); (kappa, v) = z . x for v = x in K-basis;
 - rho: E7 -> K/2K is determined by the classes g_c = rho(eps_c + eps_7), c = 1..6.  Necessary: each of the 21 classes
   rho(eps_c + eps_d) contains a norm-10 vector of K (the frame vectors w_cd), and b(g_c, g_d) = 1 (mod 2).
 - the cut of kappa under rho: p_c = z . g_c (mod 2), p_7 = 0;  k = min(wt p, 7 - wt p).
 A class is closed if EVERY candidate rho has a witness in a closed cell."""
import sys, json, itertools
import numpy as np
sys.path.insert(0, 'C:/home/researcher/Desktop/conway99-rank17')
import kstar_scan as S
D = 'C:/home/researcher/Desktop/c1_recherche/data/rank17/'
CLOSED = lambda n, k: n <= 4 or (n == 5 and k <= 1) or (n == 6 and k <= 1) or (n == 7 and k == 0) or (n == 8 and k == 0)
def load_class(key):
    V, i = key.split(':'); i = int(i)
    d = json.load(open(D + 'K_' + V + '.json', encoding='utf-8'))
    for rec in d['K']:
        if rec['id'] == i: return rec['K_gram']
def profile(G, verbose=False):
    n = len(G)
    H = S.inv_times7(G)
    wit = [(v, z) for v, z in S.short_vectors(H, 8)]
    vec10 = [x for v, x in S.short_vectors(G, 10) if v == 10]
    Gm = np.array(G, dtype=np.int64)
    cls = {}
    for x in vec10:
        key = tuple(int(t) % 2 for t in x)
        cls.setdefault(key, []).append(x)
    C = list(cls.keys()); idx = {c: i for i, c in enumerate(C)}
    Ca = np.array(C, dtype=np.int64)
    B = (Ca @ Gm @ Ca.T) % 2
    m = len(C)
    adj = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if B[i, j] == 1:
                s = tuple((a + b) % 2 for a, b in zip(C[i], C[j]))
                if s in idx: adj[i].add(j); adj[j].add(i)
    # 6-cliques (as sorted tuples)
    cliques = []
    def ext(cl, cand):
        if len(cl) == 6: cliques.append(tuple(cl)); return
        for v in sorted(cand):
            if cl and v < cl[-1]: continue
            ext(cl + [v], cand & adj[v])
    ext([], set(range(m)))
    res = []
    for cl in cliques:
        g = Ca[list(cl)]                    # 6 x 15
        # independence over F2
        prof = set()
        for val, z in wit:
            p = (g @ np.array(z, dtype=np.int64)) % 2
            w = int(p.sum()); k = min(w, 7 - w)
            prof.add((val, k))
        closed = any(CLOSED(n_, k_) for n_, k_ in prof)
        res.append((cl, prof, closed))
    return wit, len(vec10), m, cliques, res
if __name__ == '__main__':
    key = sys.argv[1]
    G = load_class(key)
    wit, n10, m, cliques, res = profile(G)
    from collections import Counter
    print(key, 'witnesses n<=8:', Counter(v for v, z in wit), '| norm-10 vectors', n10, 'classes', m, '| 6-cliques', len(cliques))
    print('  cliques NOT closed by a closed cell:', sum(1 for r in res if not r[2]))
    opencells = Counter()
    for cl, prof, closed in res:
        if not closed: opencells[tuple(sorted(prof))] += 1
    for k, v in opencells.most_common(8): print('   ', v, k)
