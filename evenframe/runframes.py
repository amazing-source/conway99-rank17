import numpy as np, sys, pickle, struct, subprocess, os, time, collections
import evenpipe as ep, skel15, one
nm = sys.argv[1]; capm = int(sys.argv[2]) if len(sys.argv) > 2 else 20
first = int(sys.argv[3]) if len(sys.argv) > 3 else 0; last = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
L = ep.Lat(nm)
good = pickle.load(open('skel_%s.pkl' % nm, 'rb'))
scls = ep.cls(L.S5)
res = collections.Counter(); t0 = time.time(); tot = 0; found = []; mx = (0, None)
cache = (None, None, None)
log = open('run_%s_%d.log' % (nm, first), 'w')
for gi, (Ib, omit, sizes) in enumerate(good[first:last]):
    if cache[0] != Ib: cache = (Ib, skel15.analyse(L, list(Ib)), one.gens_for(L, Ib), one.budgets_for(L, Ib) + one.root_components(L, Ib) + one.root_components(L, None))
    grp, comp = one.group_for(L, cache[1], Ib, omit)
    if comp[omit] != omit:
        res['SKIP equivalent omit'] += 1; continue
    cc = skel15.cell_classes(L, cache[1], omit)
    path = 'tmp_%s_%d.bin' % (nm, os.getpid())
    one.write_instance(L, scls, cc, path, capm, gens=cache[2], buds=cache[3], grp=grp)
    ts = time.time()
    out = subprocess.run(['./evenframe3.exe', path], capture_output=True, text=True).stdout
    head = out.split('\n')[0].split()
    if not head or head[0] not in ('EXCLUDED', 'FOUND', 'LIMIT'):
        res['ERROR'] += 1; print('ERROR', first + gi, out[:200], file=log, flush=True); continue
    res[head[0]] += 1; tot += int(head[1])
    print(first + gi, Ib, omit, head[0], head[1], round(time.time() - ts, 1), file=log, flush=True)
    if int(head[1]) > mx[0]: mx = (int(head[1]), first + gi)
    if head[0] != 'EXCLUDED':
        found.append((first + gi, Ib, omit, out[:4000]))
        print('!!', first + gi, Ib, omit, out.split('\n')[0][:120], flush=True)
    if gi % 250 == 0: print(gi, dict(res), 'nodes', tot, 'max', mx, round(time.time() - t0), 's', flush=True)
os.remove(path)
print(nm, 'DONE', len(good[first:last]), dict(res), 'nodes', tot, 'max', mx, round(time.time() - t0), 's', flush=True)
pickle.dump(found, open('found_%s_%d.pkl' % (nm, first), 'wb'))
