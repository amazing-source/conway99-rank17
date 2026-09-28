"""cross-validation: wfe_c.exe check mode vs Python vector_level (wfe.py) on sampled skeletons"""
import sys, os, subprocess, collections, time
sys.path.insert(0, '.')
import wfe, wfe2, wfe4
tag, pid, walks = sys.argv[1], int(sys.argv[2]), sys.argv[3]
case = wfe.Wcase(tag, pid); T = wfe2.Tables(case)
path = os.path.join(wfe.WORK, 'cval_%s_%d.bin' % (tag, pid)); os.makedirs(wfe.WORK, exist_ok=True)
t0 = time.time(); wfe4.export(case, T, path); print('export', round(time.time() - t0, 1), 's', flush=True)
env = dict(os.environ, STABMAX='5', SAMPFILE=path + '.samp')
# sampled skeletons via skel_w3 (same walk logic)
import wfe3
wfe3.export_raw(case, T, path + '.w3')
r = subprocess.run([os.path.join(wfe.HERE if hasattr(wfe, 'HERE') else '.', 'skel_w3.exe'), path + '.w3', 'est', '1', '0', '-', walks, '5'], capture_output=True, text=True, env=env)
S = [l for l in open(path + '.samp') if l.strip()]
open(path + '.skels', 'w').write(''.join(S))
t1 = time.time()
r = subprocess.run([wfe4.EXE, path, 'check', '1', '0', path, '20000', '1', '40', '20000000', path + '.skels'], capture_output=True, text=True, env=env)
cl = [l for l in r.stdout.split('\n') if l.startswith('CHECK')]
print('C check', len(cl), 'skeletons in', round(time.time() - t1, 1), 's')
res = collections.Counter()
for line, s in zip(cl, S):
    sk_ = list(map(int, s.split()))
    rr = case.vector_level(sk_)
    py = rr['verdict'] + (' ' + rr['why'].split()[0] + ' ' + rr['why'].split()[-1] if 'why' in rr else ' search')
    res[(line.replace('CHECK ', '').split(' nodes')[0], py)] += 1
for k, v in sorted(res.items(), key=str): print(v, k)
