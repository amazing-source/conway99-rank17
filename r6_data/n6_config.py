# Data check (exact short vectors, no solver): for each open R6 class, the n=6 dual vectors kappa (7|kappa|^2 = 6)
# and their pairwise 7*(kappa,kappa'). Reports how often linked pairs (|7(k,k')| = 3, so k -+ k' also has n = 6) occur.
import sys, json
from collections import Counter
sys.path.insert(0, 'C:/home/researcher/Desktop/conway99-rank17')
import kstar_scan as S
D = 'C:/home/researcher/Desktop/c1_recherche/data/rank17/'
open_ = set(json.load(open(D + 'open_after_sol.json')))
lad = {'%s:%d' % (r['V'], r['id']): r for r in json.load(open('C:/home/researcher/Desktop/conway99-rank17/out_kstar_ladder.json'))}
files = ['K_Z17', 'K_E8', 'K_D12p', 'K_E72p', 'K_D82p', 'K_A15p', 'K_A11E6p']
stat = Counter(); linked = 0; tot = 0; nvec = Counter()
for f in files:
    d = json.load(open(D + f + '.json', encoding='utf-8'))
    tag = f[2:]
    for rec in d['K']:
        key = '%s:%d' % (tag, rec['id'])
        if key not in open_ or not lad[key]['counts'].get('6', 0):
            continue
        G = rec['K_gram']; H = S.inv_times7(G)
        sv = [z for v, z in S.short_vectors(H, 6) if v == 6]   # z with z^T H z = 6, i.e. 7|kappa|^2 = 6 (kappa = G^-1 z)
        tot += 1; nvec[len(sv)] += 1
        ips = set()
        for i in range(len(sv)):
            for j in range(i + 1, len(sv)):
                a = sum(sv[i][p] * H[p][q] * sv[j][q] for p in range(15) for q in range(15))  # = 7*(k_i,k_j)
                ips.add(abs(a))
        if 3 in ips: linked += 1
        stat[tuple(sorted(ips))] += 1
print('open R6 classes checked:', tot, ' with a linked pair (|7(k,k\')|=3):', linked)
print('number of n=6 vectors (up to sign) per class:', sorted(nvec.items()))
print('most common inner-product sets |7(k,k\')|:', stat.most_common(8))
