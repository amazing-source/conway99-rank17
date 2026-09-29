"""Optional exact experiment, not a completeness reduction of the full task.
Search Type-0 N invariant under the affine order-21 group of the Paley
regular tournament. Vertex transitivity is an EXTRA assumption.
No LP/SAT/ILP/CP solver. All tests integral.
"""
import itertools as it, numpy as np, json, time
from collections import defaultdict
from pathlib import Path
D=Path(__file__).resolve().parent
def main():
 t0=time.perf_counter()
 group=[(a,b) for a in (1,2,4) for b in range(7)]; idx={g:i for i,g in enumerate(group)}
 def mul(g,h):a,b=g;c,d=h;return(a*c%7,(a*d+b)%7)
 inv=[idx[(pow(a,-1,7),(-pow(a,-1,7)*b)%7)] for a,b in group]
 rel=np.array([[idx[mul(group[inv[i]],group[j])] for j in range(21)] for i in range(21)])
 pairs=[];seen={0}
 for i in range(1,21):
  if i not in seen:pairs.append((i,inv[i]));seen|={i,inv[i]}
 assert len(pairs)==10
 R=np.zeros((42,7),dtype=np.int64)
 for i,(a,b) in enumerate(group):
  R[i,b]=R[i,(a+b)%7]=1;R[i+21,b]=1;R[i+21,(a+b)%7]=-1
 power=np.array([32**k for k in range(7)],dtype=np.int64)
 cplus=R[:21]@power;cminus=R[21:]@power
 cplus_inv=cplus[inv]
 def patterns(npairs,mpairs,root):
  out=[];by=defaultdict(list)
  for pos in it.combinations(range(10),npairs):
   avail=[k for k in range(10) if k not in pos]
   for neg in it.combinations(avail,mpairs):
    a=np.zeros(21,dtype=np.int64)
    for k in pos:a[list(pairs[k])]=1
    for k in neg:a[list(pairs[k])]=-1
    by[int(a@root)].append(len(out));out.append(a)
  return out,by
 def mat(a,b,c):
  N=np.zeros((42,42),dtype=np.int64)
  N[:21,:21]=a[rel];N[21:,21:]=b[rel]
  N[:21,21:]=c[rel];N[21:,:21]=N[:21,21:].T
  return N
 stats=[];solutions=[]
 for mode,(samep,samen,oppp,oppn,ncross) in enumerate(((1,1,2,0,3),(2,2,3,1,1))):
  As,ai=patterns(samep,samen,cplus);Bs,bi=patterns(oppp,oppn,cminus)
  ct=nr=local=0;minerr=None
  choices=list(it.combinations(range(21),ncross))
  data=[(x,sum(int(cminus[i]) for i in x),sum(int(cplus_inv[i]) for i in x),sum(1<<i for i in x)) for x in choices]
  for pp,pc,pci,pm in data:
   for mm,mc,mci,m_mask in data:
    if pm&m_mask:continue
    ct+=1
    ac=ai.get(mc-pc);bc=bi.get(mci-pci)
    if not ac or not bc:continue
    c=np.zeros(21,dtype=np.int64);c[list(pp)]=1;c[list(mm)]=-1
    for ia in ac:
     for ib in bc:
      nr+=1
      N=mat(As[ia],Bs[ib],c)
      assert np.all(N@R==0)
      deg=np.abs(N[[0,21]])@np.abs(R)
      cap=4-2*np.abs(R[[0,21]])
      if np.any(deg>cap):continue
      local+=1
      residual=N[[0,21]]@N-N[[0,21]]-12*np.eye(42,dtype=np.int64)[[0,21]]+R[[0,21]]@R.T
      err=int(np.sum(residual*residual))
      minerr=err if minerr is None else min(minerr,err)
      if err==0:
       assert np.all(N@N-N-12*np.eye(42,dtype=np.int64)+R@R.T==0)
       solutions.append(N.tolist())
  stats.append({'mode':mode,'within_same_patterns':len(As),'within_opp_patterns':len(Bs),'cross_patterns':ct,
                'NR_feasible':nr,'row_caps_feasible':local,'smallest_squared_polynomial_residual':minerr})
  print(stats[-1],flush=True)
 res={'extra_assumption':'N invariant under simultaneous affine group x -> a*x+b, a in {1,2,4}, on the two Paley edge orbits',
      'not_a_reduction_of_all_Type0':True,'solutions':len(solutions),'statistics':stats,
      'seconds':time.perf_counter()-t0}
 (D/'paley_invariant_experiment.json').write_text(json.dumps(res,indent=2)+'\n')
 if solutions:(D/'paley_invariant_matrices.json').write_text(json.dumps({'R':R.tolist(),'matrices':solutions}))
 print(res,flush=True)
if __name__=='__main__':main()
