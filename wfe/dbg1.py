import sys, json, numpy as np, collections
sys.path.insert(0, '.')
import wfe
case = wfe.Wcase('E8', 5)
S, ex, nn = case.skeletons(max_count=40, limit=10**6)
print('witness norms', case.wn)
D = case.D
G7 = np.array([[int(case.wx[a] @ D @ case.wx[b]) // 7 for b in range(len(case.wn))] for a in range(len(case.wn))])
print('7(k_a,k_b):'); print(G7)
for s in S[:6]:
    doms = case.domains(s)
    g = [0] + list(s)
    p = np.array([[bin(gc & xm).count('1') & 1 for gc in g] for xm in case.xmask])
    o = p.sum(1); k = np.minimum(o, 7 - o)
    print('skel', s, 'k per witness', k.tolist())
    print('  C1 domain sizes', [len(P) for P, _ in doms])
    r = case.vector_level(s)
    print('  ', {kk: v for kk, v in r.items() if kk != 'depth'})
