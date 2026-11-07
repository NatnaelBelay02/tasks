"""Author-side completeness check of the gyrotropic held-out guides.

usage: python3 exhaust.py heldout_guides.json certified.json nmax Bmax

For every gyrotropic guide and test frequency and every family n = 0 .. nmax, the zeros of D_n(gamma) with
|Re gamma| <= 1.02 A and |Im gamma| <= Bmax are counted with the argument principle (96 strips per sign,
the double-precision evaluator of solution/modes.py) and located; those with alpha <= A are compared with
the certified modes (beta to 1e-4).  Zeros with |alpha - A| < 0.03 A are printed (the instruction
guarantees none within 0.02 A).  This uses neither program's search region nor its stopping rule in n.
"""
import json
import math
import os
import sys
import time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'solution'))
import modes as solve   # noqa: E402
hybrid = solve
gs=json.load(open(sys.argv[1]))['guides'];orc=json.load(open(sys.argv[2]));orc=orc.get('guides',orc)
nmax=int(sys.argv[3]);Bmax=float(sys.argv[4])
A=solve.A;A1=solve.A1
tot=0;bad=0;t0=time.time()
for g in gs:
    G=solve.Guide(g)
    if not G.gyro: continue
    for j,fs in enumerate(g['test_frequencies']):
        for n in range(nmax+1):
            H=hybrid.Hyb(G,n,fs)
            zs=hybrid._strips(H,0.0,Bmax,A1,96)
            fam='LSE0' if n==0 else 'HYB%d'%n
            got={'forward':sorted(z.imag for z in zs if 0<z.real<=A),'reverse':sorted(-z.imag for z in zs if -A<=z.real<0)}
            near=[z for z in zs if abs(abs(z.real)-A)<0.03*A]
            for dr in got:
                ref=sorted(float(m[0]) for m in orc[g['id']]['modes'][j][dr].get(fam,[]))
                ok=len(ref)==len(got[dr]) and all(abs(a-b)<=1e-4*max(abs(a),1) for a,b in zip(ref,got[dr]))
                tot+=len(ref)
                if not ok:
                    bad+=1;print('MISMATCH',g['id'],fs,dr,fam,'oracle',ref,'scan',got[dr],flush=True)
            if near: print('NEAR LIMIT',g['id'],fs,fam,near,flush=True)
        print(g['id'],fs,'families 0..%d scanned to |beta| %g'%(nmax,Bmax),'t',round(time.time()-t0),flush=True)
print('modes',tot,'mismatches',bad)
