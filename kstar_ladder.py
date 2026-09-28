"""For the pairs not killed by min K* >= 5/7: the exact minimum n = 7|kappa|^2 (up to bound 14) and the count of vectors
attaining it. Tells how far the min-K* bound must be pushed to kill everything."""
import json, sys
from collections import Counter
import kstar_scan as S
D = 'C:/home/researcher/Desktop/c1_recherche/data/rank17/'
files = ['K_Z17', 'K_E8', 'K_D12p', 'K_E72p', 'K_D82p', 'K_A15p', 'K_A11E6p', 'K_E82', 'K_D16p']
prev = {(r['V'], r['id']): r for r in json.load(open('out_kstar_all.json'))}
out = []
for f in files:
    d = json.load(open(D + f + '.json', encoding='utf-8'))
    tag = f[2:]
    for rec in d['K']:
        if prev[(tag, rec['id'])]['excluded']:
            continue
        H = S.inv_times7(rec['K_gram'])
        sv = S.short_vectors(H, 14)
        vals = sorted(v for v, _ in sv)
        mn = vals[0] if vals else None
        out.append({'V': tag, 'id': rec['id'], 'min7': mn, 'n_at_min': sum(1 for v in vals if v == mn),
                    'counts': dict(Counter(vals))})
    print(tag, 'done', flush=True)
json.dump(out, open('out_kstar_ladder.json', 'w'))
c = Counter(r['min7'] for r in out)
print('min 7|kappa|^2 over the', len(out), 'surviving pairs:', dict(sorted(c.items(), key=lambda x: (x[0] is None, x[0] or 0))))
