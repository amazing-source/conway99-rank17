"""Write the text input of stageA_c.exe for a Z^17 pair (same classes / PE items as run_pair.py)."""
import sys
import stageA as A
pid = int(sys.argv[1]); out = sys.argv[2]
y1, y2 = A.load_pair(pid)
cls = A.classes_of(y1, y2); keys = sorted(cls)
pe = sorted(A.pe_constraints(y1, y2).items(), key=lambda t: t[1])
items = [(B, k, int(bin(B).count('1') == 1)) for B, k in pe if k < 3]
with open(out, 'w') as f:
    f.write(' '.join(map(str, y1)) + '\n' + ' '.join(map(str, y2)) + '\n')
    f.write('%d\n' % len(items))
    for B, k, s in items:
        f.write('%d %d %d\n' % (B, k, s))
    f.write('%d\n' % len(keys))
    for k in keys:
        f.write('%d %s\n' % (len(cls[k]), ' '.join(map(str, cls[k]))))
