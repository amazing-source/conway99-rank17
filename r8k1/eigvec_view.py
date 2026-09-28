import numpy as np
exec(open('analyze_leaves.py').read().split("lines = open")[0])
lines = open('leafdump.txt').read().strip().split('\n')
for ln in lines[:6]:
    rho, adj, h = load(ln); M, N, R = MU(rho, adj)
    w, V = np.linalg.eigh(M); v = V[:, 0]; v = v / np.abs(v).max()
    # positive-rep orientation: multiply star coordinate u by a_u (label value = h of rep e_c + s e_j = vertex 15+u)
    a = [h[u] for u in range(12)]   # h list is for vertices 15..38; vertex 15+u is the rep of label u
    Npos = N * np.outer(a + [1]*5, a + [1]*5)
    vp = v * np.array(a + [1]*5)
    att = [tuple(k for k in range(5) if Npos[u, 12+k] != 0) for u in range(12)]
    print('min eig %.3f' % w[0], ' W part:', ' '.join('%.2f' % x for x in vp[12:]))
    for u in range(12):
        j, s = u // 2, (1 if u % 2 == 0 else -1)
        stars = [(w_, int(Npos[u, w_])) for w_ in range(12) if Npos[u, w_] != 0]
        print('   label (line %d, s%+d) a=%+d att=%s  x=%.2f  star-nbrs(label,sign)=%s' % (j, s, a[u], att[u], vp[u], stars))
