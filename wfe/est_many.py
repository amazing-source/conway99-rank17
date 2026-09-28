import sys, os, json, subprocess, time
sys.path.insert(0, '.')
import wfe, wfe2, wfe4
jobs = [a.split(':') for a in sys.argv[1].split(',')]
for tag, pid in jobs:
    t0 = time.time()
    try:
        case = wfe.Wcase(tag, int(pid)); T = wfe2.Tables(case)
        path = os.path.join(wfe.WORK, 'est_%s_%s.bin' % (tag, pid)); os.makedirs(wfe.WORK, exist_ok=True)
        wfe4.export(case, T, path)
        r = subprocess.run([wfe4.EXE, path, 'est', '1', '0', path, '20000', '1'], capture_output=True, text=True, env=dict(os.environ, STABMAX='5'))
        d = wfe4.parse(r.stdout)
        print(json.dumps({'cls': '%s:%s' % (tag, pid), 'wit': case.stats['witnesses'], 'usable': len(case.usable), 'est': d.get('KNUTH'), 'samp': r.stderr.strip().split('\n')[-1][:160], 's': round(time.time() - t0)}), flush=True)
        os.remove(path)
    except Exception as e:
        print(json.dumps({'cls': '%s:%s' % (tag, pid), 'error': repr(e)[:300]}), flush=True)
