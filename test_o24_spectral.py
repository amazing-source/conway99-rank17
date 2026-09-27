"""Limited numerical O24 spectral-cut experiment; NOT a proof/certificate.

Usage: python test_o24_spectral.py [number_of_seeded_pairs]
Writes the numerical JSON log next to this script. Requires NumPy and SciPy.
The model constructs integer coordinates with squared norm a.a/8.
"""
import numpy as np, time, json, sys
from pathlib import Path
OUT_DIR = Path(__file__).resolve().parent
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

def golay():
    g=0xAE3
    words=[]
    for a in range(4096):
        w=0
        for i in range(12):
            if a>>i&1: w ^= g<<i
        w |= (w.bit_count()&1)<<23
        words.append(w)
    from collections import Counter
    assert Counter(w.bit_count() for w in words)=={0:1,8:759,12:2576,16:759,24:1}
    return words

def setup():
    words=golay()
    bits=np.array([[(w>>i)&1 for i in range(24)] for w in words],dtype=np.int64)
    us=(1-2*bits)[bits[:,0]==0]
    octads=bits[bits.sum(axis=1)==8]
    signs=np.array([[1]+[1-2*((w>>i)&1) for i in range(7)] for w in range(128) if w.bit_count()%2==0],dtype=np.int64)
    hs=np.zeros((759*64,24),dtype=np.int64)
    for j,octad in enumerate(octads):
        hs[j*64:(j+1)*64,np.flatnonzero(octad)]=2*signs
    y1=np.ones(24,dtype=np.int64); y1[:4]=-3
    return bits,us,hs,y1

def random_pair(bits,y1,rng):
    while True:
        T=rng.choice(np.arange(4,24),4,replace=False)
        y2=np.ones(24,dtype=np.int64); y2[T]=-3
        cand=(1-2*bits)*y2
        ids=np.flatnonzero(cand@y1==0)
        if len(ids):
            y2=cand[rng.choice(ids)]
            assert y2@y2==56 and y1@y2==0
            return y2

def data(us,hs,y1,y2):
    U=us[(us@y1==0)&(us@y2==0)]
    H=hs[(hs@y1==0)&(hs@y2==0)]
    iu=np.triu_indices(24)
    vecs=np.vstack([U,H])
    A=(vecs[:,iu[0]]*vecs[:,iu[1]]).T
    A=np.vstack([A, np.r_[np.ones(len(U),dtype=int),np.zeros(len(H),dtype=int)],np.r_[np.zeros(len(U),dtype=int),np.ones(len(H),dtype=int)]])
    R=56*np.eye(24,dtype=np.int64)-np.outer(y1,y1)-np.outer(y2,y2)
    b=np.r_[R[iu],42,7]
    return U,H,A,b

def run(A,b,U,H,y1,y2, maxround=25):
    t0=time.monotonic(); cuts=[];rhs=[];records=[]
    n=A.shape[1];nu=len(U)
    # Feasibility objective; any returned feasible vertex/point is acceptable.
    obj=np.zeros(n)
    for k in range(maxround+1):
        lp=linprog(obj,A_eq=csr_matrix(A),b_eq=b,A_ub=csr_matrix(np.array(cuts)) if cuts else None,b_ub=rhs or None,bounds=(0,1),method='highs-ipm')
        if lp.status!=0:
            return {'status':lp.status,'round':k,'seconds':round(time.monotonic()-t0,3),'records':records,'cuts':len(cuts)},cuts,rhs
        Hsum=(H.T*lp.x[nu:])@H
        bound=32*np.eye(24)-(4/7)*(np.outer(y1,y1)+np.outer(y2,y2))
        ev,Q=np.linalg.eigh(bound-Hsum)
        records.append({'round':k,'mineig':float(ev[0]),'negative':int((ev< -1e-6).sum())})
        if ev[0]>-1e-7:
            return {'status':'spectral feasible','round':k,'seconds':round(time.monotonic()-t0,3),'records':records,'cuts':len(cuts)},cuts,rhs
        for j in np.flatnonzero(ev< -1e-6):
            q=np.rint(Q[:,j]*1000).astype(np.int64)
            coeff=7*(H@q)**2
            right=int(224*(q@q)-4*((y1@q)**2+(y2@q)**2))
            if coeff@lp.x[nu:]>right+1e-5:
                cuts.append(np.r_[np.zeros(nu,dtype=np.int64),coeff]);rhs.append(right)
    return {'status':'round limit','round':maxround,'seconds':round(time.monotonic()-t0,3),'records':records,'cuts':len(cuts)},cuts,rhs

if __name__=='__main__':
    bits,us,hs,y1=setup();rng=np.random.default_rng(20260927)
    out=[]
    num=int(sys.argv[1]) if len(sys.argv)>1 else 5
    if num < 1:
        raise SystemExit('The number of pairs must be positive.')
    for i in range(num):
        y2=random_pair(bits,y1,rng)
        U,H,A,b=data(us,hs,y1,y2)
        res,cuts,rhs=run(A,b,U,H,y1,y2,maxround=10)
        entry={'sample':i,'y1':y1.tolist(),'y2':y2.tolist(),'nu':len(U),'nh':len(H),**res}
        print(json.dumps(entry),flush=True);out.append(entry)
        if res['status']==2 and res['round']>0:
            np.savez_compressed(OUT_DIR / ('o24_spectral_sample_%d.npz'%i),A=A,b=b,U=U,H=H,y1=y1,y2=y2,cuts=np.array(cuts),rhs=np.array(rhs))
    with (OUT_DIR / 'o24_spectral_results.json').open('w', encoding='utf-8') as stream:
        json.dump(out, stream, indent=2)
