"""Run run_pair.py sequentially on several pairs (PASS 1 sampling), each with a wall-clock cap.
Usage: python batch_pairs.py SAMPLE CAPSEC PAIR [PAIR ...]"""
import subprocess, sys, time, os
HERE = os.path.dirname(os.path.abspath(__file__))
samp, cap = sys.argv[1], int(sys.argv[2])
for p in sys.argv[3:]:
    t0 = time.time()
    out = os.path.join(HERE, 'out', 'run_p%s_s.json' % p)
    log = open(os.path.join(HERE, 'out', 'run_p%s_s.txt' % p), 'w')
    try:
        subprocess.run([sys.executable, os.path.join(HERE, 'run_pair.py'), p, '--sample', samp, '--seed', '1', '--out', out],
                       stdout=log, stderr=subprocess.STDOUT, timeout=cap)
        status = 'done'
    except subprocess.TimeoutExpired:
        status = 'CAP'
    log.close()
    print('pair %s: %s in %.0fs' % (p, status, time.time() - t0), flush=True)
