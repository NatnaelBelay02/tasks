#!/usr/bin/env python3
"""Reference filter synthesis program.

usage: python3 synthesize.py spec.json network.json

1. Approximation: the characteristic function with an equiripple passband (Remez exchange
   algorithm in u = (f / f_2)^2 for given transmission zeros) and equiripple stopbands: the
   finite transmission zeros are moved until the attenuation minima of each stopband are
   equal (scipy's fsolve on a 30-digit evaluation, then Newton's method at the working precision).
2. Feldtkeller equation E E_* = F F_* + P P_*, solved for the Hurwitz polynomial E.
3. Darlington synthesis: the input impedance Z = (E + F) / (E - F) (or with -F) is
   reduced to a resistor by removing poles at infinity and at zero, fully or partially
   (zero shifting), and finite poles as resonators; a depth-first search over the order of
   the removals keeps every element positive.
All arithmetic is done in mpmath, in units normalised to R_s and 2 pi f_2, at a precision set
from the degree and the relative bandwidth (set_precision), at least 100 digits.
"""
import json
import math
import sys

import mpmath as mp
from scipy.optimize import fsolve

mp.mp.dps = 100
EPS = mp.mpf(10) ** -60          # relative size below which a coefficient is rounding noise
TOL = {'div': mp.mpf(10) ** -40, 'dc': mp.mpf(10) ** -45, 'jw': mp.mpf(10) ** -30, 'newton': mp.mpf(10) ** -110}
ADMISSIBLE_ONLY = False          # True: keep only ladders whose load lies in [R_s / 2, 2 R_s]


def set_precision(spec):
    """Working precision for this specification.

    The polynomials are held as coefficients in s.  For a band-pass filter of relative bandwidth
    b = (f2 - f1) / sqrt(f1 f2) their roots crowd into two clusters of width ~ b, and every step
    of the extraction then loses about n log10(1 / b) digits to cancellation.  The precision is
    set to 40 + 3 times that loss (at least 100 digits), and every tolerance is a fixed power of
    the precision, so the margin between rounding noise and the tolerances grows with the loss.
    """
    global EPS
    f1, f2 = [float(x) for x in spec['passband_hz']]
    loss = int(spec['degree']) * math.log10(math.sqrt(f1 * f2) / (f2 - f1)) if f1 > 0 else 0.0
    dps = max(100, 10 * math.ceil((40 + 3 * max(loss, 0.0)) / 10))
    mp.mp.dps = dps
    EPS = mp.mpf(10) ** -int(0.6 * dps)
    TOL.update(div=mp.mpf(10) ** -int(0.4 * dps), dc=mp.mpf(10) ** -int(0.45 * dps),
               jw=mp.mpf(10) ** -int(0.3 * dps), newton=mp.mpf(10) ** -int(1.1 * dps))
    return dps


# ------------------------------------------------------------ polynomials ---
# coefficient lists, lowest power first

def trim(p, rel=None):
    p = list(p)
    rel = EPS if rel is None else rel
    big = max([abs(c) for c in p] + [mp.mpf(0)])
    while len(p) > 1 and abs(p[-1]) <= rel * big:
        p.pop()
    return p


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pscale(a, c):
    return [c * x for x in a]


def pmul(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pshift(a, k=1):
    return [mp.mpf(0)] * k + list(a)


def peval(a, x):
    v = 0
    for c in reversed(a):
        v = v * x + c
    return v


def pdiv(a, b):
    """Quotient of a / b; the remainder (checked to be negligible) is dropped."""
    a = list(a)
    q = [mp.mpf(0)] * max(1, len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        q[i] = a[i + len(b) - 1] / b[-1]
        for j, y in enumerate(b):
            a[i + j] -= q[i] * y
    rem = max([abs(x) for x in a[:len(b) - 1]] + [mp.mpf(0)])
    scale = max(abs(x) for x in q + list(b))
    if rem > TOL['div'] * scale * max(1, max(abs(x) for x in b)):
        raise ArithmeticError('inexact division')
    return q


def from_roots(roots):
    p = [mp.mpf(1)]
    for r in roots:
        p = pmul(p, [-r, 1])
    return p


# ---------------------------------------------------------- approximation ---

class Approx:
    """|K(j W)|^2 = c^2 u^r prod (u - ur_i)^2 / (u^a prod (u - uz_k)^2),  u = W^2, W = f / f2."""

    def __init__(self, spec):
        f1, f2 = [mp.mpf(str(x)) for x in spec['passband_hz']]
        self.n = int(spec['degree'])
        self.a = int(spec['transmission_zeros']['dc'])
        self.bands = [((mp.mpf(str(b['edge_hz'])) / f2) ** 2, int(b['zeros']))
                      for b in spec['transmission_zeros']['stopbands'] if int(b['zeros']) > 0]
        self.uz = []
        self.r = self.n % 2
        self.d = (self.n - self.r) // 2
        self.u1 = (f1 / f2) ** 2
        self.al = (1 + self.u1) / 2
        self.be = (1 - self.u1) / 2
        rho2 = mp.power(10, -mp.mpf(str(spec['return_loss_db'])) / 10)
        self.eps2 = rho2 / (1 - rho2)

    def u(self, x):
        return self.al + self.be * x

    def W(self, x):
        u = self.u(x)
        w = mp.power(u, mp.mpf(self.r - self.a) / 2)
        for z in self.uz:
            w /= abs(u - z)
        return w

    def dlogW(self, x):
        u = self.u(x)
        g = mp.mpf(self.r - self.a) / (2 * u) if self.r != self.a else mp.mpf(0)
        for z in self.uz:
            g -= 1 / (u - z)
        return self.be * g

    def newton_roots(self, f, df):
        """Refine the roots of the previous solve by Newton's method; None unless they come out real,
        increasing, inside (-1, 1) and with f changing sign at each one (then they are all d roots)."""
        r0 = getattr(self, 'warm_roots', None)
        if r0 is None or len(r0) != len(f) - 1:
            return None
        tol = mp.mpf(10) ** -int(mp.mp.dps * 0.85)
        out = []
        for x in r0:
            for _ in range(60):
                dx = peval(f, x) / peval(df, x)
                x -= dx
                if abs(dx) < tol:
                    break
            else:
                return None
            out.append(x)
        if any(not -1 < x < 1 for x in out) or any(b <= a for a, b in zip(out, out[1:])):
            return None
        edges = [mp.mpf(-1)] + [(a + b) / 2 for a, b in zip(out, out[1:])] + [mp.mpf(1)]
        for k in range(len(out)):
            if peval(f, edges[k]) * peval(f, edges[k + 1]) >= 0:
                return None
        return out

    def remez(self):
        d = self.d
        tol_gap = mp.mpf(10) ** -int(mp.mp.dps * 0.5)
        tol_x = mp.mpf(10) ** -int(mp.mp.dps * 0.8)
        tol_level = mp.mpf(10) ** -int(mp.mp.dps * 0.6)
        tol_stall = mp.mpf(10) ** -int(mp.mp.dps * 0.25)
        prev_spread, stall = mp.inf, 0
        if d == 0:
            self.f = [mp.mpf(1)]
            self.ext = [mp.mpf(1)]
            return
        lowpass_odd = self.u1 == 0 and self.r - self.a > 0
        if getattr(self, 'warm', None) is not None:
            pts = list(self.warm)          # extremal points of the previous solve (zeros moved a little)
        elif lowpass_odd:
            pts = [-mp.cos(mp.pi * (i + mp.mpf(1) / 2) / (d + mp.mpf(1) / 2)) for i in range(d + 1)]
        else:
            pts = [-mp.cos(mp.pi * i / d) for i in range(d + 1)]
        for it in range(100):
            A = mp.matrix(d + 1, d + 1)
            rhs = mp.matrix(d + 1, 1)
            for i, x in enumerate(pts):
                w = self.W(x)
                for j in range(d):
                    A[i, j] = w * x ** j
                A[i, d] = -(-1) ** i
                rhs[i] = -w * x ** d
            sol = mp.lu_solve(A, rhs)
            f = [sol[j] for j in range(d)] + [mp.mpf(1)]
            df = [j * f[j] for j in range(1, len(f))]
            roots = self.newton_roots(f, df)
            if roots is None:
                roots = sorted(mp.re(z) for z in mp.polyroots(f[::-1], maxsteps=500, extraprec=400))
            self.warm_roots = roots
            g = lambda x: peval(df, x) / peval(f, x) + self.dlogW(x)
            ext = []
            edges = [mp.mpf(-1)] + roots + [mp.mpf(1)]
            for k in range(d + 1):
                lo, hi = edges[k], edges[k + 1]
                span = hi - lo
                if k == 0 and not lowpass_odd and g(lo) <= 0:
                    ext.append(lo)
                    continue
                if k == d and g(hi) >= 0:
                    ext.append(hi)
                    continue
                a_ = lo + span * tol_gap if k > 0 or lowpass_odd else lo
                b_ = hi - span * tol_gap if k < d else hi
                for _ in range(400):
                    m = (a_ + b_) / 2
                    if g(m) > 0:
                        a_ = m
                    else:
                        b_ = m
                    if b_ - a_ < tol_x:
                        break
                ext.append((a_ + b_) / 2)
            vals = [self.W(x) * abs(peval(f, x)) for x in ext]
            pts = ext
            spread = (max(vals) - min(vals)) / max(vals)
            if spread < tol_level:
                break
            # at low working precision (the double-precision zero search runs at 30 digits) a very wide
            # band reaches a rounding floor above tol_level: stop once the spread no longer shrinks
            if it >= 3 and spread < tol_stall and spread > prev_spread / 4:
                stall += 1
                if stall >= 2:
                    break
            else:
                stall = 0
            prev_spread = spread
        self.f = f
        self.ext = ext
        self.roots = roots
        self.warm = ext

    # ---- stopbands: log |W f| as a function of u, and its minima ----
    def logwf(self, u):
        v = mp.mpf(self.r - self.a) / 2 * mp.log(u)
        for q in self.ur:
            v += mp.log(abs(u - q))
        for z in self.uz:
            v -= mp.log(abs(u - z))
        return v

    def dlogwf(self, u):
        g = mp.mpf(self.r - self.a) / (2 * u)
        for q in self.ur:
            g += 1 / (u - q)
        for z in self.uz:
            g -= 1 / (u - z)
        return g

    def stationary(self, lo, hi):
        """The point where d/du log |W f| changes from negative to positive in (lo, hi)."""
        if hi is None:
            hi = 2 * lo
            while self.dlogwf(hi) <= 0:
                hi *= 2
        if lo is None:
            lo = hi / 2
            while self.dlogwf(lo) >= 0:
                lo /= 2
        a, b = mp.log(lo), mp.log(hi)
        for _ in range(int(mp.mp.dps * 4)):
            m = (a + b) / 2
            if self.dlogwf(mp.exp(m)) < 0:
                a = m
            else:
                b = m
        return mp.exp((a + b) / 2)

    def solve_for(self, uz):
        self.uz = list(uz)
        self.remez()
        self.ur = [self.u(z) for z in getattr(self, 'roots', [])]

    def zeros_from(self, p):
        """Transmission zeros of every stopband, ordered and outside its edge, from parameters p."""
        out, k = [], 0
        for ue, m in self.bands:
            t = mp.mpf(0)
            for j in range(m):
                t += mp.exp(p[k + j])
                out.append(ue * mp.exp(t) if ue > 1 else ue * mp.exp(-t))
            k += m
        return out

    def stopband_residual(self, p):
        """log |W f| at each interior stopband minimum minus its value at the band edge."""
        self.solve_for(self.zeros_from(p))
        res, k = [], 0
        for ue, m in self.bands:
            zs = sorted(self.uz[k:k + m])
            k += m
            segs = list(zip(zs, zs[1:] + [None])) if ue > 1 else list(zip([None] + zs[:-1], zs))
            edge = self.logwf(ue)
            for lo, hi in segs:
                lo2 = lo * (1 + mp.mpf(10) ** -12) if lo is not None else None
                hi2 = hi * (1 - mp.mpf(10) ** -12) if hi is not None else None
                res.append(self.logwf(self.stationary(lo2, hi2)) - edge)
        return res

    def place_zeros(self):
        if not self.bands:
            self.solve_for([])
            return
        p0 = []
        for ue, m in self.bands:
            s = abs(float(mp.log(ue)))
            p0 += [math.log(0.05 * s)] + [math.log(0.5 * s)] * (m - 1)
        with mp.workdps(30):
            p1 = fsolve(lambda p: [float(x) for x in self.stopband_residual([mp.mpf(float(v)) for v in p])], p0,
                        xtol=1e-13)
        p = mp.findroot(lambda *v: self.stopband_residual(list(v)), [mp.mpf(float(v)) for v in p1],
                        tol=TOL['newton'], maxsteps=30)
        p = [p] if not hasattr(p, 'rows') else [p[i] for i in range(p.rows)]
        self.solve_for(self.zeros_from(p))

    def build(self):
        """Normalised polynomials F(s), P(s), E(s) (s in units of 2 pi f2)."""
        self.place_zeros()
        ur = self.ur

        def k2(x):
            u = self.u(x)
            v = u ** self.r / u ** self.a
            for q in ur:
                v *= (u - q) ** 2
            for z in self.uz:
                v /= (u - z) ** 2
            return v
        c2 = self.eps2 / max(k2(x) for x in self.ext)
        self.c = mp.sqrt(c2)
        # F(s) = c s^r prod (s^2 + ur),  P(s) = s^a prod (s^2 + uz)
        F = [self.c]
        for q in ur:
            F = pmul(F, [q, 0, 1])
        F = pshift(F, self.r)
        P = [mp.mpf(1)]
        for z in self.uz:
            P = pmul(P, [z, 0, 1])
        P = pshift(P, self.a)
        # Q(w) = P(s)P(-s) + F(s)F(-s), w = s^2
        Pw = [mp.mpf(1)]
        for z in self.uz:
            Pw = pmul(Pw, pmul([z, 1], [z, 1]))
        Pw = pscale(pshift(Pw, self.a), (-1) ** self.a)
        Fw = [c2]
        for q in ur:
            Fw = pmul(Fw, pmul([q, 1], [q, 1]))
        Fw = pscale(pshift(Fw, self.r), (-1) ** self.r)
        Q = padd(Pw, Fw)
        wr = mp.polyroots(Q[::-1], maxsteps=2000, extraprec=1000)
        sr = [-mp.sqrt(w) for w in wr]
        E = [mp.mpf(1)]
        for z in sr:
            E = pmul_c(E, [-z, 1])
        E = [mp.re(x) * self.c for x in E]
        self.F, self.P, self.E = F, P, E
        return F, P, E


def pmul_c(a, b):
    out = [mp.mpc(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


# ------------------------------------------------------------- extraction ---

def lead(p):
    return trim(p)[-1]


def deg(p):
    return len(trim(p)) - 1


def is_zero_at0(p):
    p = trim(p)
    return abs(p[0]) <= TOL['dc'] * max(abs(c) for c in p)


def search(N, D, zeros, path, out, limit):
    """N/D is the impedance looking into the rest of the ladder."""
    if len(out) >= limit:
        return
    N, D = trim(N), trim(D)
    if deg(N) == 0 and deg(D) == 0:
        if not zeros:
            R = N[0] / D[0]
            if R > 0 and (not ADMISSIBLE_ONLY or admissible(R)):
                out.append((list(path), R))
        return
    for view in ('Z', 'Y'):
        num, den = (N, D) if view == 'Z' else (D, N)
        if deg(num) == deg(den) + 1:
            Lfull = lead(num) / lead(den)
            # full removal of the pole at infinity
            rest = trim(padd(num, pscale(pshift(den), -Lfull)))
            rest = rest[:len(num) - 1] if len(rest) >= len(num) else rest
            el = ('series', 'L', Lfull) if view == 'Z' else ('shunt', 'C', Lfull)
            nxt = (rest, den) if view == 'Z' else (den, rest)
            path.append(el)
            search(nxt[0], nxt[1], zeros, path, out, limit)
            path.pop()
            for uz in sorted(set(zeros)):
                w = mp.sqrt(uz)
                X = peval(num, 1j * w) / peval(den, 1j * w)
                if abs(mp.re(X)) > TOL['jw'] * abs(X):
                    continue
                Lp = mp.im(X) / w
                if not (0 < Lp < Lfull * (1 - mp.mpf(10) ** -20)):
                    continue
                try:
                    q1 = pdiv(trim(padd(num, pscale(pshift(den), -Lp))), [uz, 0, 1])
                    k2 = mp.re(peval(den, 1j * w) / (1j * w * peval(q1, 1j * w)))
                    if not k2 > 0:
                        continue
                    q2 = pdiv(trim(padd(den, pscale(pshift(q1), -k2))), [uz, 0, 1])
                except ArithmeticError:
                    continue
                if view == 'Z':
                    els = [('series', 'L', Lp), ('shunt', 'series_LC', 1 / k2, k2 / uz)]
                    nxt = (q1, q2)
                else:
                    els = [('shunt', 'C', Lp), ('series', 'parallel_LC', k2 / uz, 1 / k2)]
                    nxt = (q2, q1)
                rem = list(zeros)
                rem.remove(uz)
                path.extend(els)
                search(nxt[0], nxt[1], rem, path, out, limit)
                del path[-2:]
        if is_zero_at0(den) and not is_zero_at0(num):
            den1 = trim(den)[1:]
            K0 = num[0] / den1[0]
            rest = trim(padd(num, pscale(den1, -K0)))
            rest = rest[1:]
            el = ('series', 'C', 1 / K0) if view == 'Z' else ('shunt', 'L', 1 / K0)
            nxt = (rest, den1) if view == 'Z' else (den1, rest)
            path.append(el)
            search(nxt[0], nxt[1], zeros, path, out, limit)
            path.pop()
            for uz in sorted(set(zeros)):
                w = mp.sqrt(uz)
                X = peval(num, 1j * w) / peval(den, 1j * w)
                if abs(mp.re(X)) > TOL['jw'] * abs(X):
                    continue
                Kp = -w * mp.im(X)
                if not (0 < Kp < K0 * (1 - mp.mpf(10) ** -20)):
                    continue
                try:
                    q1 = pdiv(trim(padd(num, pscale(den1, -Kp))), [uz, 0, 1])
                    k2 = mp.re(peval(den1, 1j * w) / peval(q1, 1j * w))
                    if not k2 > 0:
                        continue
                    q2 = pshift(pdiv(trim(padd(den1, pscale(q1, -k2))), [uz, 0, 1]))
                except ArithmeticError:
                    continue
                if view == 'Z':
                    els = [('series', 'C', 1 / Kp), ('shunt', 'series_LC', 1 / k2, k2 / uz)]
                    nxt = (q1, q2)
                else:
                    els = [('shunt', 'L', 1 / Kp), ('series', 'parallel_LC', k2 / uz, 1 / k2)]
                    nxt = (q2, q1)
                rem = list(zeros)
                rem.remove(uz)
                path.extend(els)
                search(nxt[0], nxt[1], rem, path, out, limit)
                del path[-2:]


def admissible(R):
    return mp.mpf(1) / 2 * (1 - mp.mpf(10) ** -12) <= R <= 2 * (1 + mp.mpf(10) ** -12)


def synthesize(spec, limit=10 ** 6, sigmas=(1, -1)):
    set_precision(spec)
    ap = Approx(spec)
    F, P, E = ap.build()
    sols = []
    for sg in sigmas:
        N = padd(E, pscale(F, sg))
        D = padd(E, pscale(F, -sg))
        N, D = trim(N), trim(D)
        out = []
        search(N, D, list(ap.uz), [], out, limit)
        sols += [(sg, p, R) for p, R in out]
        if len(sols) >= limit:
            break
    return ap, sols


def denormalise(spec, path, R):
    Rs = float(spec['source_resistance_ohm'])
    w0 = 2 * math.pi * float(spec['passband_hz'][1])
    arms = []
    for el in path:
        pos, typ = el[0], el[1]
        if typ == 'L':
            arms.append({'position': pos, 'type': 'L', 'L': float(el[2]) * Rs / w0})
        elif typ == 'C':
            arms.append({'position': pos, 'type': 'C', 'C': float(el[2]) / (Rs * w0)})
        else:
            arms.append({'position': pos, 'type': typ, 'L': float(el[2]) * Rs / w0,
                         'C': float(el[3]) / (Rs * w0)})
    return {'load_resistance_ohm': float(R) * Rs, 'arms': arms}


def choose(spec, sols):
    """Fewest inductors among the ladders whose load lies in [R_s / 2, 2 R_s]; then the load closest to R_s."""
    best = None
    for sg, path, R in sols:
        if not admissible(R):
            continue
        key = (sum(1 for e in path if 'L' in e[1]), abs(mp.log(R)))
        if best is None or key < best[0]:
            best = (key, path, R)
    return best


def main():
    global ADMISSIBLE_ONLY
    spec = json.load(open(sys.argv[1]))
    # the first 40 ladders with an admissible load (from both signs of F) are enough to choose from
    ADMISSIBLE_ONLY = True
    ap, sols = synthesize(spec, limit=40)
    best = choose(spec, sols)
    if best is None:
        raise SystemExit('no ladder with an admissible load found')
    json.dump(denormalise(spec, best[1], best[2]), open(sys.argv[2], 'w'), indent=1)


if __name__ == '__main__':
    main()
