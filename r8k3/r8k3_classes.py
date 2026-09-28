"""R8 k=3 witness patterns (orbit level, orbit_sat conventions), exact enumeration + reduction by symmetry (C2 solo 28/09).
h = +-1 on the 24 cut orbits (cells (i,d), i in I={1,2,3}, d in J={4..7}), +-2 on exactly two non-cut orbits, 0 elsewhere,
block balance (INN) at the 14 inner vertices.  Symmetry: S3(I) x S4(J) x per-block flips x global sign.
Output: data/r8k3/classes.json (one representative per class, with orbit size)."""
import itertools, json, os, numpy as np
m = 7
outer = [(i, j, e) for i in range(1, m + 1) for j in range(i + 1, m + 1) for e in (0, 1)]
N = len(outer)
rep = [frozenset([(i, 0), (j, e)]) for (i, j, e) in outer]
trep = [frozenset([(i, 1), (j, 1 - e)]) for (i, j, e) in outer]
where = {}
for U in range(N): where[rep[U]] = (U, 1); where[trep[U]] = (U, -1)
inner = [(i, s) for i in range(1, m + 1) for s in (0, 1)]
I, J = (1, 2, 3), (4, 5, 6, 7)
cut = [U for U, (i, j, e) in enumerate(outer) if (i in I) != (j in I)]
zer = [U for U in range(N) if U not in cut]
INN = np.zeros((len(inner), N), dtype=np.int64)
for t, c in enumerate(inner):
    for U in range(N):
        if c in rep[U]: INN[t, U] += 1
        if c in trep[U]: INN[t, U] -= 1

def act(hv, perm, flips, sg):
    out = [0] * N
    for U in range(N):
        if not hv[U]: continue
        img = frozenset((perm[i], s ^ flips[perm[i]]) for (i, s) in rep[U])
        V, sgn = where[img]; out[V] = sg * sgn * hv[U]
    return tuple(out)

group = []
for pI in itertools.permutations(I):
    for pJ in itertools.permutations(J):
        perm = {**dict(zip(I, pI)), **dict(zip(J, pJ))}
        for fl in itertools.product((0, 1), repeat=m):
            flips = {i + 1: fl[i] for i in range(m)}
            for sg in (1, -1): group.append((perm, flips, sg))

import sys
def canon(p):
    return min(act(p, *g) for g in group)
placements = {}
for (A, B) in itertools.combinations(zer, 2):
    for sa in (2, -2):
        for sb in (2, -2):
            z = [0] * N; z[A] = sa; z[B] = sb; z = tuple(z)
            placements.setdefault(canon(z), []).append(z)
preps = sorted(placements)
if "--placements" in sys.argv:
    for t, z in enumerate(preps):
        nz = [(outer[U], z[U]) for U in range(N) if z[U]]
        print(t, nz, "orbit", len(placements[z]))
    sys.exit()
t = int(sys.argv[sys.argv.index("--place") + 1]); z0 = np.array(preps[t])
stab = [g for g in group if act(tuple(z0), *g) == tuple(z0)]
Tcols = np.array(cut); M = INN[:, Tcols]; rhs = -(INN @ z0); pats = []
for hi in range(16):
    idx = np.arange(hi << 20, (hi + 1) << 20, dtype=np.int64)
    S = ((idx[:, None] >> np.arange(24)) & 1) * 2 - 1
    for row in S[np.all(S @ M.T == rhs[None, :], axis=1)]:
        hh = z0.copy(); hh[Tcols] = row; pats.append(tuple(int(x) for x in hh))
pset = set(pats); seen = set(); reps_ = []
for q in pats:
    if q in seen: continue
    orb = set(act(q, *g) for g in stab); assert orb <= pset
    seen |= orb; reps_.append((q, len(orb)))
os.makedirs("data/r8k3", exist_ok=True)
json.dump({"placement": t, "placement_orbit": len(placements[preps[t]]), "stab": len(stab), "patterns": len(pats),
           "classes": [{"h": list(q), "orbit": s_} for q, s_ in reps_]}, open("data/r8k3/place_%d.json" % t, "w"))
print("placement %d: stab %d, patterns %d, classes %d" % (t, len(stab), len(pats), len(reps_)))
