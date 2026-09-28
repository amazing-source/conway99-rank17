# Data check (exact, no solver): distribution of short-triple types (n(k), n(k'), n(k-k')) with n(k)=6, all <= 8,
# over the open R6 classes; and how many open classes (all 754) contain SOME triple with all three norms in {5,...,8}.
import sys, json
from collections import Counter
sys.path.insert(0, 'C:/home/researcher/Desktop/conway99-rank17')
import kstar_scan as S
D = 'C:/home/researcher/Desktop/c1_recherche/data/rank17/'
open_ = set(json.load(open(D + 'open_after_sol.json')))
files = ['K_Z17', 'K_E8', 'K_D12p', 'K_E72p', 'K_D82p', 'K_A15p', 'K_A11E6p']
typ = Counter(); cls_with = Counter(); any_tri = 0; n_open = 0
for f in files:
    d = json.load(open(D + f + '.json', encoding='utf-8'))
    tag = f[2:]
    for rec in d['K']:
        key = '%s:%d' % (tag, rec['id'])
        if key not in open_: continue
        n_open += 1
        H = S.inv_times7(rec['K_gram'])
        sv = S.short_vectors(H, 8)
        ip = lambda x, y: sum(x[p] * H[p][q] * y[q] for p in range(15) for q in range(15))
        found = set()
        for i in range(len(sv)):
            for j in range(i + 1, len(sv)):
                (a, x), (b, y) = sv[i], sv[j]
                t = ip(x, y)
                for s in (1, -1):
                    c = a + b - 2 * s * t
                    if 5 <= c <= 8:
                        found.add(tuple(sorted((a, b, c))))
        if found: any_tri += 1
        for tp in found: cls_with[tp] += 1
print('open classes:', n_open, ' containing a short triple with all norms in 5..8:', any_tri)
for tp, c in sorted(cls_with.items(), key=lambda z: -z[1])[:10]:
    print('  triple type (n,n\',n\'\') =', tp, ':', c, 'classes')
