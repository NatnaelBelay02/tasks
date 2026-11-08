"""Second implementation of the synthesis, written separately from solution/synthesize.py.

usage: ALT_DPS=120 python3 alt_synthesize.py spec.json network.json

Passband polynomial from the closed-form generalized-Chebyshev phase sum in w^2 (bisection),
transmission zeros by Newton's method on the stopband equal-minima conditions, spectral
factorization with mpmath.polyroots, and a depth-first ladder extraction with partial removals
at 0 and infinity.  ALT_DPS sets the working precision (default 120 digits).  The README
reports its results at several precisions.
"""
import sys, json, math, os
import mpmath as mp
DPS = int(os.environ.get('ALT_DPS', '120'))
mp.mp.dps = DPS

def main():
    spec = json.load(open(sys.argv[1], encoding='utf-8'))
    Rs = float(spec['source_resistance_ohm'])
    f1, f2 = [float(v) for v in spec['passband_hz']]
    RLdb = float(spec['return_loss_db']); n = int(spec['degree'])
    a = int(spec['transmission_zeros']['dc'])
    bands = spec['transmission_zeros']['stopbands']
    f0 = math.sqrt(f1*f2) if f1 > 0 else f2
    w0 = 2*math.pi*f0
    y1 = (mp.mpf(f1)/f0)**2; y2 = (mp.mpf(f2)/f0)**2
    eps = 1/mp.sqrt(mp.power(10, mp.mpf(RLdb)/10)-1)
    m = sum(int(b['zeros']) for b in bands)
    Ninf = n - a - 2*m
    r = n % 2; d = (n-r)//2
    def xof(y): return (2*y-y1-y2)/(y2-y1)
    x0 = xof(mp.mpf(0))
    def clamp(v): return max(mp.mpf(-1), min(mp.mpf(1), v))
    def theta(y, yk):
        x = clamp(xof(y))
        T = Ninf*mp.acos(x)
        if a: T += a*mp.acos(clamp((1-x*x0)/(x-x0)))
        for q in yk:
            xp = xof(q); T += 2*mp.acos(clamp((1-x*xp)/(x-xp)))
        return T
    def resp(yk):
        yi = []
        for j in range(d):
            tgt = (2*j+1)*mp.pi
            lo, hi = y1, y2  # theta decreasing in y
            for it in range(200):
                mid = (lo+hi)/2
                if theta(mid, yk) > tgt: lo = mid
                else: hi = mid
                if hi-lo < mp.mpf(10)**(-mp.mp.dps+3)*y2: break
            yi.append((lo+hi)/2)
        num = y2**(mp.mpf(r)/2); den = y2**(mp.mpf(a)/2)
        for q in yi: num *= abs(y2-q)
        for q in yk: den *= abs(y2-q)
        c = eps*den/num
        return yi, c
    def logK(y, yi, c, yk):
        v = mp.log(c) + mp.mpf(r-a)/2*mp.log(y)
        for q in yi: v += mp.log(abs(y-q))
        for q in yk: v -= mp.log(abs(y-q))
        return v
    def dlogK(y, yi, yk):
        v = mp.mpf(r-a)/(2*y)
        for q in yi: v += 1/(y-q)
        for q in yk: v -= 1/(y-q)
        return v
    def bis(g, lo, hi):
        for it in range(300):
            mid = mp.sqrt(lo*hi)
            if g(mid) < 0: lo = mid
            else: hi = mid
            if hi/lo - 1 < mp.mpf(10)**(-mp.mp.dps+5): break
        return mp.sqrt(lo*hi)
    binfo = []
    for b in bands:
        ye = (mp.mpf(float(b['edge_hz']))/f0)**2
        binfo.append((ye, int(b['zeros']), ye > y2))
    def zeros_from_u(u):
        yk = []; idx = 0
        for ye, mb, up in binfo:
            for j in range(mb):
                t = 1+mp.exp(u[idx]); idx += 1
                yk.append(ye*t if up else ye/t)
        return yk
    def eqs(u):
        yk = zeros_from_u(u)
        yi, c = resp(yk)
        g = lambda y: dlogK(y, yi, yk)
        out = []; idx = 0
        for ye, mb, up in binfo:
            zs = sorted(yk[idx:idx+mb]); idx += mb
            le = logK(ye, yi, c, yk)
            pts = []
            for i in range(mb-1):
                pts.append(bis(g, zs[i]*(1+mp.mpf(10)**-30), zs[i+1]*(1-mp.mpf(10)**-30)))
            if up:
                hi = zs[-1]*2
                while g(hi) < 0: hi *= 2
                pts.append(bis(g, zs[-1]*(1+mp.mpf(10)**-30), hi))
            else:
                lo = zs[0]/2
                while g(lo) > 0: lo /= 2
                pts.append(bis(g, lo, zs[0]*(1-mp.mpf(10)**-30)))
            for p in pts: out.append(logK(p, yi, c, yk)-le)
        return out
    yk = []
    if m:
        best = None
        for scale in [1, 0.3, 3, 0.1, 10]:
            u = []
            for ye, mb, up in binfo:
                for j in range(mb):
                    u.append(mp.log(mp.mpf(0.1)*scale*(4**j)))
            try:
                mp.mp.dps = 25
                u = [mp.mpf(v) for v in u]
                ok = False
                for prec in [25, DPS]:
                    mp.mp.dps = prec
                    u = [mp.mpf(v) for v in u]
                    for it in range(100):
                        F0 = eqs(u); nr = mp.norm(mp.matrix(F0))
                        if nr < mp.mpf(10)**(-prec+8): break
                        h = mp.mpf(10)**(-prec//2)
                        J = mp.matrix(m, m)
                        for k in range(m):
                            uu = list(u); uu[k] += h
                            Fk = eqs(uu)
                            for i in range(m): J[i, k] = (Fk[i]-F0[i])/h
                        du = mp.lu_solve(J, -mp.matrix(F0))
                        lam = mp.mpf(1)
                        while lam > 1e-4:
                            un = [u[i]+lam*du[i] for i in range(m)]
                            try:
                                if mp.norm(mp.matrix(eqs(un))) < nr: break
                            except Exception: pass
                            lam /= 2
                        u = un
                    ok = nr < mp.mpf(10)**(-prec+10)
                if ok:
                    best = u; break
            except Exception as e:
                mp.mp.dps = DPS
                continue
        mp.mp.dps = DPS
        yk = zeros_from_u(best)
    yi, c = resp(yk)
    def pmul(p, q):
        o = [mp.mpf(0)]*(len(p)+len(q)-1)
        for i, x in enumerate(p):
            for j, z in enumerate(q): o[i+j] += x*z
        return o
    def padd(p, q):
        o = [mp.mpf(0)]*max(len(p), len(q))
        for i, x in enumerate(p): o[i] += x
        for i, x in enumerate(q): o[i] += x
        return o
    def psc(p, k): return [k*x for x in p]
    def pev(p, s):
        v = 0
        for x in reversed(p): v = v*s+x
        return v
    F = [c] if r == 0 else [mp.mpf(0), c]
    for q in yi: F = pmul(F, [q, 0, 1])
    P = [mp.mpf(0)]*a+[mp.mpf(1)]
    for q in yk: P = pmul(P, [q, 0, 1])
    Qt = [((-1)**r)*c*c] ; Qt = [mp.mpf(0)]*r + Qt
    for q in yi: Qt = pmul(Qt, pmul([q, 1], [q, 1]))
    Pt = [mp.mpf(0)]*a + [mp.mpf((-1)**a)]
    for q in yk: Pt = pmul(Pt, pmul([q, 1], [q, 1]))
    Qt = padd(Qt, Pt)
    roots = mp.polyroots(list(reversed(Qt)), maxsteps=500, extraprec=400)
    E = [mp.mpc(c)]
    for t in roots: E = pmul(E, [mp.sqrt(t), 1])
    E = [mp.re(x) for x in E]
    def Kmag2(y):
        s = mp.mpc(0, mp.sqrt(y))
        return abs(pev(F, s)/pev(P, s))**2
    def trim(p):
        p = list(p); M = max(abs(x) for x in p)
        while len(p) > 1 and abs(p[-1]) < M*mp.mpf(10)**-28: p.pop()
        return p
    def divquad(p, w2):
        p = list(p); q = [mp.mpf(0)]*(len(p)-2)
        for i in range(len(p)-1, 1, -1):
            co = p[i]; q[i-2] = co; p[i] -= co; p[i-2] -= co*w2
        return q
    tol = mp.mpf(10)**-20
    def verify(arms, R):
        fs = []
        lo = 0.0 if f1 == 0 else f1/30
        import numpy as np
        for t in np.linspace(0, 1, 120): fs.append(lo+(f2-lo)*t)
        for t in np.geomspace(f2, 30*f2, 120): fs.append(t)
        if f1 > 0:
            for t in np.geomspace(f1/30, f1, 80): fs.append(t)
        for q in yk:
            for dd in [1e-3, -1e-3, 1e-6, -1e-6]: fs.append(float(mp.sqrt(q))*f0*(1+dd))
        for ff in fs:
            if ff <= 0: ff = f2*1e-6
            y = (mp.mpf(ff)/f0)**2; s = mp.mpc(0, mp.sqrt(y))
            A, B, C, D = 1, 0, 0, 1
            for pos, typ, L, Cc in arms:
                if typ == 'L': z = s*L
                elif typ == 'C': z = 1/(s*Cc)
                elif typ == 'series_LC': z = s*L+1/(s*Cc)
                else: z = 1/(1/(s*L)+s*Cc)
                if pos == 'series': A, B, C, D = A, A*z+B, C, C*z+D
                else: A, B, C, D = A+B/z, B, C+D/z, D
            G = 4/R/abs(A+B/R+C+D/R)**2
            K2 = Kmag2(y)
            if abs(10*mp.log10(G*(1+K2))) > 1e-7: return False
        return True
    sys.setrecursionlimit(10000)
    best = [None]           # (inductors, arms, R) of the best ladder so far
    def nind(arms): return sum(1 for e in arms if e[1] in ('L', 'series_LC', 'parallel_LC'))
    def dfs(N, D, isZ, ninf, ndc, fin, arms):
        N = trim(N); D = trim(D)
        if len(N) == 1 and len(D) == 1:
            if ninf or ndc or fin: return None
            R = N[0]/D[0] if isZ else D[0]/N[0]
            if R <= 0 or not (0.5 <= R <= 2): return None
            if (best[0] is None or nind(arms) < best[0][0]) and verify(arms, R):
                best[0] = (nind(arms), arms, R)
            return None
        if len(arms) > n+m: return None
        # branch and bound: each remaining finite zero still needs a resonator (one inductor)
        if best[0] is not None and nind(arms) + len(fin) >= best[0][0]: return None
        for kindZ in (isZ, not isZ):
            num, den = (N, D) if kindZ == isZ else (D, N)
            pinf = len(num) == len(den)+1
            Md = max(abs(x) for x in den)
            pdc = abs(den[0]) < Md*tol
            if pdc: den = [mp.mpf(0)]+den[1:]
            for w2 in sorted(set(fin)):
                w = mp.sqrt(w2); s = mp.mpc(0, w)
                val = pev(num, s)/pev(den, s)
                for mode in ('inf', 'dc'):
                    if mode == 'inf' and pinf:
                        k = mp.re(val/s); kf = num[-1]/den[-1]
                        if not (k > kf*1e-12 and k < kf*(1-1e-12)): continue
                        n2 = padd(num, psc([0]+den, -k)); el = (('series', 'L', k, None) if kindZ else ('shunt', 'C', None, k))
                    elif mode == 'dc' and pdc:
                        k = mp.re(val*s); d1 = den[1:]; kf = num[0]/d1[0]
                        if not (k > kf*1e-12 and k < kf*(1-1e-12)): continue
                        n2 = padd([0]+num, psc(den[1:], -k))[0:]
                        n2 = padd(psc(num, 1), psc(d1, -k)); n2 = n2
                        el = (('series', 'C', None, 1/k) if kindZ else ('shunt', 'L', 1/k, None))
                    else: continue
                    dd2 = den
                    qn = divquad(n2, w2)
                    kr = mp.re(pev(dd2, s)/(s*pev(qn, s)))
                    if kr <= 0: continue
                    newnum = divquad(padd(dd2, psc([0]+qn, -kr)), w2)
                    if kindZ: el2 = ('shunt', 'series_LC', 1/kr, kr/w2)
                    else: el2 = ('series', 'parallel_LC', kr/w2, 1/kr)
                    fl = list(fin); fl.remove(w2)
                    res = dfs(newnum, qn, not kindZ, ninf, ndc, tuple(fl), arms+[el, el2])
            if pinf and ninf > 0:
                k = num[-1]/den[-1]
                n2 = padd(num, psc([0]+den, -k))[:-1]
                el = ('series', 'L', k, None) if kindZ else ('shunt', 'C', None, k)
                if k > 0:
                    res = dfs(n2, den, kindZ, ninf-1, ndc, fin, arms+[el])
            if pdc and ndc > 0:
                d1 = den[1:]; k = num[0]/d1[0]
                n2 = padd(num, psc(d1, -k))[1:]
                el = ('series', 'C', None, 1/k) if kindZ else ('shunt', 'L', 1/k, None)
                if k > 0:
                    res = dfs(n2, d1, kindZ, ninf, ndc-1, fin, arms+[el])
        return None
    for sg in (1, -1):
        Fs = psc(F, sg)
        # E and F have the same leading coefficient, so one of E + F, E - F loses its top term
        N = trim(padd(E, Fs)); D = trim(padd(E, psc(Fs, -1)))
        dfs(N, D, True, Ninf, a, tuple(yk), [])
    if best[0] is None:
        sys.exit('synthesis failed')
    _, arms, R = best[0]
    out = []
    for pos, typ, L, Cc in arms:
        o = {'position': pos, 'type': typ}
        if L is not None: o['L'] = float(L*Rs/w0)
        if Cc is not None: o['C'] = float(Cc/(Rs*w0))
        out.append(o)
    with open(sys.argv[2], 'w', encoding='utf-8') as fh:
        json.dump({'load_resistance_ohm': float(R*Rs), 'arms': out}, fh, indent=1)

if __name__ == '__main__':
    main()
