"""Exact tournament/characteristic-anchor enumeration; no solver.

A bit 1 on ij (i<j) means an arrow i -> j.  After the characteristic
triangle-parity gauge, Phi^-_ij = +1 for i->j, -1 otherwise.
Type-II Phi^+ = +2 on the positive defect edge and -2 on the negative.
"""
from itertools import combinations, permutations, product
from collections import Counter
import json, time, argparse
N=7
EDGES=list(combinations(range(N),2))
INDEX={e:i for i,e in enumerate(EDGES)}
ALL=(1<<len(EDGES))-1

def enumerate_scores(target):
    """All labelled tournaments with the specified outdegree vector.
    Each unordered edge is assigned once; prune only impossible degrees.
    """
    deg=[0]*N
    rem=[N-1]*N
    out=[]
    def rec(k,mask):
        if k==len(EDGES):
            assert deg==list(target)
            out.append(mask); return
        i,j=EDGES[k]
        rem[i]-=1; rem[j]-=1
        for winner,bit in ((i,1),(j,0)):
            deg[winner]+=1
            if all(deg[v]<=target[v]<=deg[v]+rem[v] for v in (i,j)):
                rec(k+1,mask | (bit<<k))
            deg[winner]-=1
        rem[i]+=1; rem[j]+=1
    rec(0,0)
    return out

def transform(mask,p,reverse=False):
    result=0
    for e,(i,j) in enumerate(EDGES):
        a,b=p[i],p[j]
        bit=(mask>>e)&1
        if a>b:
            a,b=b,a; bit^=1
        if reverse: bit^=1
        result|=bit<<INDEX[a,b]
    return result

def stabilizer_permutations(kind):
    if kind=='type0':
        return [(p,False) for p in permutations(range(N))]+[(p,True) for p in permutations(range(N))]
    # In the intersecting case the shared vertex is 0, low 1, high 2.
    if kind=='intersecting':
        ans=[]
        for tail in permutations((3,4,5,6)):
            ans.append(((0,1,2)+tail,False))
            ans.append(((0,2,1)+tail,True))
        return ans
    if kind=='disjoint':
        ans=[]
        for rev in (False,True):
            lo=(2,3) if rev else (0,1)
            hi=(0,1) if rev else (2,3)
            for a,b,c in product(permutations(lo),permutations(hi),permutations((4,5,6))):
                ans.append((a+b+c,rev))
        return ans
    raise ValueError(kind)

def data(kind):
    if kind=='type0': return [3]*N,[],[]
    if kind=='intersecting': return [3,2,4,3,3,3,3],[0,1],[0,2]
    if kind=='disjoint': return [2,2,4,4,3,3,3],[0,1],[2,3]
    raise ValueError(kind)

def enumerate_orbits(kind):
    scores,plus,minus=data(kind)
    boards=set(enumerate_scores(scores)); total=len(boards)
    group=stabilizer_permutations(kind)
    rest=set(boards); reps=[]
    while rest:
        rep=min(rest)
        orbit={transform(rep,p,r) for p,r in group}
        assert orbit<=boards
        assert orbit<=rest # distinct group orbits cannot partially overlap
        rest.difference_update(orbit)
        reps.append({'mask':rep,'orbit_size':len(orbit),
                     'signs':[1 if rep>>k&1 else -1 for k in range(21)]})
    return {'kind':kind,'scores':scores,'positive_defect':plus,'negative_defect':minus,
            'labelled_count':total,'group_size':len(group),'orbit_count':len(reps),
            'orbit_size_histogram':dict(Counter(x['orbit_size'] for x in reps)),
            'representatives':reps}

def verify_anchor(record,kind):
    scores,plus,minus=data(kind)
    signs=record['signs']; balance=[0]*N
    norm=0
    for (i,j),s in zip(EDGES,signs):
        balance[i]+=s; balance[j]-=s; norm+=s*s
    for edge,value in ((plus,2),(minus,-2)):
        if edge:
            for i in edge: balance[i]+=value
            norm+=value*value
    assert balance==[0]*N
    assert norm=={'type0':21,'intersecting':29,'disjoint':29}[kind]

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',default='anchors.json')
    args=ap.parse_args(); t=time.perf_counter()
    result={'definition':'1 bit on ij means i -> j; quotient by labelled normal-form stabilizer plus overall sign'}
    for kind in ('type0','intersecting','disjoint'):
        result[kind]=enumerate_orbits(kind)
        for rec in result[kind]['representatives']: verify_anchor(rec,kind)
        print(kind,{k:v for k,v in result[kind].items() if k!='representatives'},flush=True)
    result['elapsed_seconds']=time.perf_counter()-t
    with open(args.output,'w') as f: json.dump(result,f,indent=2)
    print('seconds',result['elapsed_seconds'])
