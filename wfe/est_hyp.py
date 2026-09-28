"""HYPOTHETICAL gain measurement only (not a proof step): Knuth estimates with restricted k-rules
(e.g. n=8 witnesses k=1 only, as if R8 k=2,3 were proven; n=6 k=3 only, as if J2 were closed)."""
import sys, os, json, subprocess, time
sys.path.insert(0, '.')
import wfe, wfe2, wfe4
cls_, mode = sys.argv[1], sys.argv[2]
tag, pid = cls_.split(':')
if mode in ('r8', 'both'):
    wfe.K_ALLOWED[8] = {1}
if mode in ('j2', 'both'):
    wfe.K_ALLOWED[6] = {3}
t0 = time.time()
case = wfe.Wcase(tag, int(pid)); T = wfe2.Tables(case)
path = os.path.join(wfe.WORK, 'hyp_%s_%s_%s.bin' % (tag, pid, mode)); os.makedirs(wfe.WORK, exist_ok=True)
wfe4.export(case, T, path)
r = subprocess.run([wfe4.EXE.replace('wfe_c.exe', 'wfe_c2.exe'), path, 'est', '1', '0', path, '20000', '1'], capture_output=True, text=True, env=dict(os.environ, STABMAX='5'))
d = wfe4.parse(r.stdout)
print(json.dumps({'cls': cls_, 'mode': mode, 'wit': case.stats['witnesses'], 'usable': len(case.usable), 'est': d.get('KNUTH'), 'samp': r.stderr.strip().split('\n')[-1][:160], 's': round(time.time() - t0)}), flush=True)
os.remove(path)
