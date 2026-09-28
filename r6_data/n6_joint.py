# Data check (exact short vectors, no solver): in each open R6 class, for every n=6 dual vector k, list its partners k'
# with 7|k'|^2 = m <= 8 and 7(k,k') = t; a partner with m + 6 - 2|t| <= 8 gives a third short vector k -+ k'.
import sys, json
from collections import Counter
sys.path.insert(0, 'C:/home/researcher/Desktop/conway99-rank17')
import kstar_scan as S
D = 'C:/home/researcher/Desktop/c1_recherche/data/rank17/'
open_ = set(json.load(open(D + 'open_after_sol.json')))
lad = {'%s:%d' % (r['V'], r['id']): r for r in json.load(open('C:/home/researcher/Desktop/conway99-rank17/out_kstar_ladder.json'))}
files = ['K_Z17', 'K_E8', 'K_D12p', 'K_E72p', 'K_D82p', 'K_A15p', 'K_A11E6p']
pat = Counter(); tri = 0; tot = 0
for f in files:
    d = json.load(open(D + f + '.json', encoding='utf-8'))
    tag = f[2:]
    for rec in d['K']:
        key = '%s:%d' % (tag, rec['id'])
        if key not in open_ or not lad[key]['counts'].get('6', 0):
            continue
        H = S.inv_times7(rec['K_gram'])
        sv = S.short_vectors(H, 8)
        six = [z for v, z in sv if v == 6]
        ip = lambda x, y: sum(x[p] * H[p][q] * y[q] for p in range(15) for q in range(15))
        tot += 1
        best = None; has_tri = False
        for z in six:
            parts = []
            for v, w in sv:
                if w == z: continue
                t = abs(ip(z, w))
                parts.append((v, t))
                if v + 6 - 2 * t <= 8 and not (v == 6 and t == 6): has_tri = True
            parts = tuple(sorted(Counter(parts).items()))
            best = parts if best is None else min(best, parts)
        tri += has_tri
        pat[best] += 1
print('open R6 classes:', tot, '  with a short triple k, k\', k-+k\' (all n<=8):', tri)
for p, c in pat.most_common(6):
    print(c, 'classes; partners (n, |7(k,k\')|) -> count:', p)
