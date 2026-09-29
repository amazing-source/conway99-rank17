"""Exact, standalone checks and 42-coordinate export. No graph existence claimed."""
from itertools import combinations, product
import json, pathlib
D=pathlib.Path(__file__).resolve().parent
E=list(combinations(range(7),2)); I={e:k for k,e in enumerate(E)}
R=[]
for i,j in E:
 for sign in (1,-1):
  r=[0]*7;r[i]=1;r[j]=sign;R.append(r)

def balance(phi): return [sum(phi[u]*R[u][c] for u in range(42)) for c in range(7)]
def triangles(phi):
 return [sum(phi[2*I[tuple(sorted(e))]] for e in ((i,j),(i,k),(j,k)))%2
         for i,j,k in combinations(range(7),3)]
def switch(phi,eps):
 out=[0]*42
 for k,(i,j) in enumerate(E):
  swap=eps[i]*eps[j]<0
  out[2*k+(1 if swap else 0)]=eps[i]*phi[2*k]
  out[2*k+(0 if swap else 1)]=eps[i]*phi[2*k+1]
 return out

def gauge(phi):
 bit=[0]+[phi[2*I[0,j]]%2 for j in range(1,7)]
 return switch(phi,[1-2*b for b in bit])

def main():
 anchors=json.load(open(D/'anchors.json')); independent=json.load(open(D/'independent_counts.json'))
 out={'source_commit':'d421b842dc4eafa501c704c8d4822353b36f5aa9',
      'indexing':'vertices 0..6; cells lexicographic; Same then Opp',
      'warning':'These are complete characteristic ANCHORS, not completed N matrices or frame candidates.',
      'R':R,'anchors':[]}
 transformed=0
 for kind in ('type0','intersecting','disjoint'):
  rec=anchors[kind]; verify=independent[kind]
  assert rec['labelled_count']==verify['labelled_count']
  assert rec['orbit_count']==verify['orbit_count']
  assert rec['orbit_size_histogram']==verify['orbit_size_histogram']
  assert [x['mask'] for x in rec['representatives']]==verify['representative_masks']
  for j,a in enumerate(rec['representatives']):
   phi=[0]*42
   for k,s in enumerate(a['signs']): phi[2*k+1]=s
   if kind!='type0':
    phi[2*I[tuple(rec['positive_defect'])]]=2
    phi[2*I[tuple(rec['negative_defect'])]]=-2
   norm=sum(x*x for x in phi)
   assert norm==(21 if kind=='type0' else 29)
   assert balance(phi)==[0]*7
   assert not any(triangles(phi))
   assert all((phi[2*k]-phi[2*k+1])%2 for k in range(21))
   for eps in product((-1,1),repeat=7):
    moved=switch(phi,eps)
    assert sum(x*x for x in moved)==norm
    assert balance(moved)==[0]*7
    assert not any(triangles(moved))
    g=gauge(moved)
    assert all(g[2*k]%2==0 and g[2*k+1]%2==1 for k in range(21))
    assert sum(g[2*k] for k in range(21))==0
    assert norm%8==5
    transformed+=1
   out['anchors'].append({'id':f'{kind}:{j:02d}','norm':norm,'Phi':phi,
                         'orbit_size_in_normal_form':a['orbit_size']})
 # Strictness control: balanced half-edges need not satisfy triangle parity.
 bad=[0]*42
 for k,(i,j) in enumerate(E):
  hi=1 if (j-i)%7 in (1,2,5) else -1
  hj=1 if (i-j)%7 in (1,2,5) else -1
  label=0 if hi==hj else 1
  bad[2*k+label]=hi
 assert sum(x*x for x in bad)==21
 assert balance(bad)==[0]*7
 assert all((bad[2*k]-bad[2*k+1])%2 for k in range(21))
 assert any(triangles(bad))
 # Cell-pattern enum by excess above 21, with even-coordinate sum-zero.
 shapes=[]
 for n3 in range(3):
  for n2 in range(5):
   for n4 in range(2):
    if n2+n4>21: continue
    m=21+8*n3+4*n2+16*n4
    if m>37: continue
    if not any(sum(s*v for s,v in zip(signs,[2]*n2+[4]*n4))==0
               for signs in product((-1,1),repeat=n2+n4)):
     continue
    # A 3-entry needs 12 in a row of 10, ignoring its own coordinate.
    vals=[3]*n3+[1]*(21-n3)+[2]*n2+[4]*n4+[0]*(21-n2-n4)
    row_ok=all(4*abs(v)<=sum(sorted((abs(w) for q,w in enumerate(vals) if q!=i),reverse=True)[:10])
               for i,v in enumerate(vals) if v)
    shapes.append({'m':m,'odd_threes':n3,'even_twos':n2,'even_fours':n4,'row_capacity_survives':row_ok})
 result={'anchors_verified':len(out['anchors']),'all_coordinate_switches_verified':transformed,
         'independent_full_scan_agrees':True,'labelled_tournaments_scanned':independent['boards_scanned'],
         'invalid_balanced_Type0_rejected':True,
         'invalid_anchor_first_triangle':{'vertices':[0,1,2], 'sum_Same':sum(bad[2*I[e]] for e in ((0,1),(0,2),(1,2)))},
         'small_magnitude_patterns':shapes,
         'scope':'No N-matrix completion was searched or certified.'}
 (D/'anchors_42.json').write_text(json.dumps(out,indent=2)+'\n')
 (D/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
