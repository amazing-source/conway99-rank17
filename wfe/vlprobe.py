import sys, os, subprocess, collections, time
sys.path.insert(0, '.')
import wfe, wfe2, wfe3
tag, pid, walks = sys.argv[1], int(sys.argv[2]), sys.argv[3]
case = wfe.Wcase(tag, pid); T = wfe2.Tables(case)
path = os.path.join(wfe.WORK, 'vlp_%s_%d.bin' % (tag, pid)); os.makedirs(wfe.WORK, exist_ok=True)
wfe3.export_raw(case, T, path)
env = dict(os.environ, STABMAX='5', SAMPFILE=path + '.samp')
r = subprocess.run([wfe3.EXE, path, 'est', '1', '0', '-', walks, '3'], capture_output=True, text=True, env=env)
print(r.stderr.strip().split('\n')[-1])
S = [list(map(int, l.split())) for l in open(path + '.samp') if l.strip()]
allp = list(range(len(T.pairs)))
res = collections.Counter(); t0 = time.time()
for s in S:
    v, info = wfe2.leaf_eval(T, s, allp)
    rr = case.vector_level(s)
    key = (v, rr['verdict'] + (' ' + rr['why'].split()[0] + ' ' + rr['why'].split()[-1] if 'why' in rr else ' nodes=%d' % rr['nodes']))
    res[key] += 1
for k, v in sorted(res.items(), key=str): print(v, k)
print('per leaf s', round((time.time() - t0) / max(1, len(S)), 3))
