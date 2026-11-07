"""cost of the whitened errors-in-variables likelihood (states profiled out) at given parameters."""
import json, sys, numpy as np, model as M
def cost(d, params, cal_d, frames):
    net=M.load_json(f"{d}/network.json"); tr=M.load_json(f"{d}/truth.json")
    Z=M.read_pmu(f"{d}/pmu.csv") if hasattr(M,"read_pmu") else None
    Y=M.ybus(net["candidates"], tr["status"], params, net["shunt"])
    rows=[b-1 for b in M.PMU_BUSES]; free=[i for i in range(M.NB) if i!=M.Z7]
    S=np.zeros((M.NB,13),complex); S[free,np.arange(13)]=1; S[M.Z7]=-Y[M.Z7,free]/Y[M.Z7,M.Z7]
    cal=np.ones(26,complex)
    for k,(c,b) in enumerate(M.CHANNELS):
        if f"{c}{b}" in cal_d: e=cal_d[f"{c}{b}"]; cal[k]=(1+e["ratio_error"])*np.exp(1j*e["phase_error"])
    A=np.concatenate([np.eye(13),(Y@S)[rows]],0)*cal[:,None]
    sm=np.array([net["noise"][c]["magnitude_rel_std"] for c,_ in M.CHANNELS]); sa=np.array([net["noise"][c]["angle_std_rad"] for c,_ in M.CHANNELS])
    tot=0.0
    for s0 in range(0,len(frames),1000):
        z=Z[frames[s0:s0+1000]]; n=len(z); u=np.exp(-1j*np.angle(z)); sc=np.abs(z)
        def white(a):  # a (n,26,...) complex
            r=a*u[...,None] if a.ndim==3 else a*u
            m=sc[...,None] if a.ndim==3 else sc
            return np.concatenate([r.real/(m*(sm[:,None] if a.ndim==3 else sm)), r.imag/(m*(sa[:,None] if a.ndim==3 else sa))],1)
        Ab=np.broadcast_to(A,(n,26,13)); Jx=white(np.concatenate([Ab,1j*Ab],2)); y=white(z)
        x=np.linalg.lstsq(Jx[0],y[0],rcond=None) if False else None
        Q,_=np.linalg.qr(Jx); res=y-(Q@(Q.transpose(0,2,1)@y[:,:,None]))[...,0]
        tot+=float(np.sum(res**2))
    return tot
if __name__=="__main__":
    d,out=sys.argv[1],sys.argv[2]
    tr=M.load_json(f"{d}/truth.json"); bad=np.zeros(9000,bool)
    for s in tr["corrupted"]: bad[s["first_frame"]-5:s["last_frame"]+6]=True
    fr=np.flatnonzero(~bad)
    O=json.load(open(out)); p={}
    for b in O["branches"]:
        if str(b["id"]) in tr["params"]: p[str(b["id"])]={k:b[k] for k in tr["params"][str(b["id"])]}
    ct=cost(d,tr["params"],tr["cal"],fr); ce=cost(d,p,O["channels"],fr)
    print(d, "frames",len(fr),"cost truth",ct,"cost estimate",ce,"diff",ct-ce, "per frame", ct/len(fr))
