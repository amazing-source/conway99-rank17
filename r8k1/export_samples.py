"""export a few T1 samples (from the scratchpad pickles) into r8k1's --check format for cross-validation"""
import pickle, glob, sys
sys.path.insert(0, 'C:/home/researcher/AppData/Local/Temp/claude/C--Users-RESEARCHER-Desktop/c88c221a-efb9-4373-af95-1e500b2e81aa/scratchpad')
from t1c import build_graph, STAR, V, SV, C
import numpy as np
S = pickle.load(open(glob.glob('C:/home/researcher/AppData/Local/Temp/claude/C--Users-RESEARCHER-Desktop/c88c221a-efb9-4373-af95-1e500b2e81aa/scratchpad/t1_samples_0_20.pkl')[0], 'rb'))
def cidx(v):
    if v == 'x': return 0
    if isinstance(v, tuple) and v[0] == 'i': return 1 + 2*v[1] + (0 if v[2] == 1 else 1)
    r = v
    if r[C] != 0:
        j = [t for t in range(6) if r[t] != 0][0]
        return 15 + (0 if r[C] == 1 else 12) + 2*j + (0 if r[j] == 1 else 1)
    if r in V: return 39 + V.index(r)
    return 44 + SV.index(r)
out = open('check_in.txt', 'w'); ref = []
for s in S[:12]:
    adj = build_graph(s['eps'], s['G'], s['Q'], s['P'])
    masks = [0]*49
    for a in adj:
        for b in adj[a]: masks[cidx(a)] |= 1 << cidx(b)
    rho = [list(v[:6]) for v in V]
    out.write(' '.join(str(x) for r in rho for x in r) + ' ' + ' '.join('%x' % m for m in masks) + '\n')
    ref.append(np.linalg.eigvalsh(s['M']).min())
out.close(); print('python M_U min eigs:', [round(x, 4) for x in ref])
