"""Run stage B on one skeleton (raw line of 17 hex columns) with a large node cap; save any solution.
Usage: python one_skel.py PAIR HEXLINE NODECAP OUT"""
import json, sys, time
import numpy as np
import stageB as B

pid, line, cap, out = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), sys.argv[4]
z = B.Z17(pid)
cols, pis = B.parse_raw(line)
sk = B.Skeleton(z, pis).build()
t0 = time.time()
st = {}
v = B.solve(sk, cap, st)
res = {'pair': pid, 'line': line, 'verdict': v, 'nodes': st['nodes'], 'depth': st.get('depth'), 'seconds': round(time.time() - t0, 1), 'dom': sk.sizes()}
print(json.dumps(res), flush=True)
if 'sol' in st:
    res['sol'] = st['sol']
json.dump(res, open(out, 'w'))
