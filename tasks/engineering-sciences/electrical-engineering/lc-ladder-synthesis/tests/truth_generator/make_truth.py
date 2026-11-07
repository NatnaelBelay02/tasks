"""Reference characteristic functions and their certificates.

usage: python3 make_truth.py specs_dir truth.json [solution_dir]

For every specification the characteristic function (equiripple passband, equiripple
stopbands) is computed here independently of the reference solution:
  1. starting values in double precision: for trial transmission zeros a discrete Remez
     exchange on a 20 000-point grid in u = (f / f_2)^2 gives the passband polynomial, and
     scipy's fsolve moves the zeros until the attenuation minima of each stopband are equal;
  2. one Newton solve in mpmath (60 digits) of all the conditions together: the passband
     equioscillation W(xi_j) f(xi_j) = (-1)^(d-j) h at d + 1 extremal points with
     d/du log|W f| = 0 at the interior ones, and in every stopband log|W f| = lambda_b at the
     band edge and at the m_b interior minima, with d/du log|W f| = 0 there.  The unknowns
     are the reflection zeros, the transmission zeros, the extremal points, h and lambda_b;
  3. certificates: (a) |K| is evaluated on 6000 points of the passband and every local
     maximum is refined by golden-section search; the maximum over the passband must equal
     epsilon to 1e-30 and be attained, with alternating sign of K, at d + 1 points (by the
     alternation theorem the passband polynomial is then the unique minimax one for these
     zeros); (b) in every stopband |K| is evaluated on 4000 logarithmically spaced points
     between the edge and beyond the last zero, every local minimum is refined, and the
     minimum over the band must be attained at exactly m_b + 1 points (the edge and one
     between consecutive zeros, or between the last zero and dc or infinity), all equal to
     1e-30.
If a solution directory is given, its ladder search is used to enumerate the canonical
ladders of the reference function, for the statistics in the README; one ladder with an
admissible load is stored as a witness that every specification is realisable (the ablation
script checks it with the verifier's own network analysis).  The verifier uses only the
characteristic function.
"""
import json
import math
import os
import sys

import mpmath as mp
import numpy as np
from scipy.optimize import fsolve

DPS = 60


def setup(spec):
    f1, f2 = float(spec['passband_hz'][0]), float(spec['passband_hz'][1])
    n = int(spec['degree'])
    a = int(spec['transmission_zeros']['dc'])
    r = n % 2
    rho2 = mp.power(10, -mp.mpf(str(spec['return_loss_db'])) / 10)
    bands = []
    for b in spec['transmission_zeros']['stopbands']:
        ue = (mp.mpf(str(b['edge_hz'])) / mp.mpf(str(f2))) ** 2
        bands.append(dict(ue=ue, m=int(b['zeros']), up=ue > 1))
    return dict(n=n, a=a, r=r, d=(n - r) // 2, u1=(mp.mpf(str(f1)) / mp.mpf(str(f2))) ** 2,
                eps2=rho2 / (1 - rho2), bands=bands)


# ------------------------------------------------ double-precision start ---

def float_remez(P, uz):
    """Discrete weighted Remez on a grid for given transmission zeros: start values."""
    d, r, a = P['d'], P['r'], P['a']
    u1 = float(P['u1'])
    uz = np.array(uz, float)
    t = np.cos(np.linspace(np.pi, 0, 20001))
    u = u1 + (1 - u1) * (t + 1) / 2
    if r > a and u1 == 0:
        u = u[1:]
    W = u ** ((r - a) / 2)
    for z in uz:
        W = W / np.abs(u - z)
    x = 2 * (u - u1) / (1 - u1) - 1
    V = np.polynomial.chebyshev.chebvander(x, d)
    idx = np.linspace(0, len(u) - 1, d + 1).astype(int)
    for _ in range(200):
        A = np.column_stack([W[idx, None] * V[idx, :d], -(-1.0) ** np.arange(d + 1)])
        sol = np.linalg.solve(A, -W[idx] * V[idx, d])
        coef = np.append(sol[:d], 1.0)
        E = W * (V @ coef)
        s = np.sign(E)
        runs = np.split(np.arange(len(E)), np.nonzero(s[1:] != s[:-1])[0] + 1)
        new = np.array([run[np.argmax(np.abs(E[run]))] for run in runs if len(run)])
        if len(new) != d + 1 or np.array_equal(new, idx):
            break
        idx = new
    roots = np.polynomial.chebyshev.chebroots(coef).real
    return np.sort(u1 + (1 - u1) * (roots + 1) / 2), u[idx]


def mp_remez(P, uz, npts=4001, dps=40):
    """The same discrete Remez exchange at `dps` digits, for wide band-pass filters (f2 / f1 > 10).
    There several reflection zeros lie within ~1e-7 of the lower band edge in x, where f is of order
    (1e-7)^3 or smaller and lost in the double-precision rounding of its Chebyshev coefficients; the
    grid adds points spaced geometrically in u towards u1, which a grid in x cannot resolve."""
    d, r, a = P['d'], P['r'], P['a']
    with mp.workdps(dps):
        u1 = mp.mpf(P['u1'])
        t = [-mp.cos(mp.pi * k / (npts - 1)) for k in range(npts)]
        grid = [u1 + (1 - u1) * (x + 1) / 2 for x in t]
        grid += [u1 * mp.power(1 / u1, (x + 1) / 2) for x in t]
        u = sorted(set(grid))
        uzm = [mp.mpf(float(z)) for z in uz]
        W = []
        for x in u:
            w = mp.power(x, mp.mpf(r - a) / 2)
            for z in uzm:
                w /= abs(x - z)
            W.append(w)
        X = [2 * (x - u1) / (1 - u1) - 1 for x in u]
        V = []
        for x in X:
            row = [mp.mpf(1), x]
            while len(row) < d + 1:
                row.append(2 * x * row[-1] - row[-2])
            V.append(row[:d + 1])
        idx = [int(round(k * (len(u) - 1) / d)) for k in range(d + 1)]
        coef = None
        for _ in range(200):
            A = mp.matrix(d + 1, d + 1)
            rhs = mp.matrix(d + 1, 1)
            for i, k in enumerate(idx):
                for j in range(d):
                    A[i, j] = W[k] * V[k][j]
                A[i, d] = -(-1) ** i
                rhs[i] = -W[k] * V[k][d]
            sol = mp.lu_solve(A, rhs)
            coef = [sol[j] for j in range(d)] + [mp.mpf(1)]
            E = [W[k] * sum(c * v for c, v in zip(coef, V[k])) for k in range(len(u))]
            new, run = [], [0]
            for k in range(1, len(E)):
                if (E[k] > 0) != (E[run[0]] > 0):
                    new.append(max(run, key=lambda i: abs(E[i])))
                    run = [k]
                else:
                    run.append(k)
            new.append(max(run, key=lambda i: abs(E[i])))
            if len(new) != d + 1 or new == idx:
                break
            idx = new
        # Chebyshev to power series in x, then the roots back to u
        Tk = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]]
        while len(Tk) < d + 1:
            p1, p0 = Tk[-1], Tk[-2]
            nxt = [mp.mpf(0)] + [2 * c for c in p1]
            for i, c in enumerate(p0):
                nxt[i] -= c
            Tk.append(nxt)
        pw = [mp.mpf(0)] * (d + 1)
        for c, tk in zip(coef, Tk):
            for i, v in enumerate(tk):
                pw[i] += c * v
        xr = sorted(mp.re(z) for z in mp.polyroots(pw[::-1], maxsteps=400, extraprec=2 * dps))
        roots = [float(u1 + (1 - u1) * (x + 1) / 2) for x in xr]
        return np.array(roots), np.array([float(u[k]) for k in idx])


def glog(P, q, z, u):
    """log |W f| at u (double precision)."""
    return ((P['r'] - P['a']) / 2 * math.log(u) + sum(math.log(abs(u - x)) for x in q)
            - sum(math.log(abs(u - x)) for x in z))


def golden_min(f, lo, hi):
    """Minimum of f on [lo, hi], searched in log u.  The start grid is clustered geometrically
    towards both ends, where a minimum next to a transmission zero of a narrow-band filter lies."""
    t = np.geomspace(1e-10, 0.5, 400)
    xs = np.unique(np.log(np.concatenate([lo + (hi - lo) * t, hi - (hi - lo) * t])))[1:-1]
    k = int(np.argmin([f(math.exp(x)) for x in xs]))
    a, b = xs[max(k - 1, 0)], xs[min(k + 1, len(xs) - 1)]
    g = (math.sqrt(5) - 1) / 2
    for _ in range(100):
        m1, m2 = b - g * (b - a), a + g * (b - a)
        if f(math.exp(m1)) < f(math.exp(m2)):
            b = m2
        else:
            a = m1
    return math.exp((a + b) / 2)


def band_zeros(P, v):
    """Transmission zeros of every stopband from unconstrained parameters (ordered, outside the edge)."""
    out, k = [], 0
    for b in P['bands']:
        ue, m = float(b['ue']), b['m']
        steps = np.cumsum(np.exp(v[k:k + m]))
        k += m
        out.append(sorted(ue * np.exp(steps if b['up'] else -steps)))
    return out


def interior_minima(P, q, z, b, zb):
    """One attenuation minimum between consecutive zeros, and between the outermost zero and dc or infinity."""
    f = lambda u: glog(P, q, z, u)
    if b['up']:
        segs = list(zip(zb, zb[1:] + [zb[-1] * 1e4]))
    else:
        segs = list(zip([zb[0] * 1e-4] + zb[:-1], zb))
    return [golden_min(f, lo * (1 + 1e-9), hi * (1 - 1e-9)) for lo, hi in segs]


def float_stage(P):
    v0 = []
    for b in P['bands']:
        s = abs(math.log(float(b['ue'])))
        v0 += list(np.log([0.05 * s] + [0.5 * s] * (b['m'] - 1)))

    remez = mp_remez if 0 < float(P['u1']) < 1e-2 else float_remez     # f2 / f1 > 10: extended precision

    def residual(v):
        zs = band_zeros(P, v)
        z = [x for zb in zs for x in zb]
        q, _ = remez(P, z)
        out = []
        for b, zb in zip(P['bands'], zs):
            g0 = glog(P, q, z, float(b['ue']))
            out += [glog(P, q, z, y) - g0 for y in interior_minima(P, q, z, b, zb)]
        return out
    v = fsolve(residual, v0) if v0 else []
    zs = band_zeros(P, v)
    z = [x for zb in zs for x in zb]
    q, ext = remez(P, z)
    ys = [interior_minima(P, q, z, b, zb) for b, zb in zip(P['bands'], zs)]
    return zs, q, ext, ys


# ------------------------------------------------------ 60-digit Newton ---

def S(P, rho, z, u):
    """Signed W f: u^((r - a)/2) prod (u - rho_i) / prod |u - z_k|."""
    v = mp.power(u, mp.mpf(P['r'] - P['a']) / 2)
    for x in z:
        v /= abs(u - x)
    for x in rho:
        v *= (u - x)
    return v


def dlog(P, rho, z, u):
    g = mp.mpf(P['r'] - P['a']) / (2 * u) if P['r'] != P['a'] else mp.mpf(0)
    for x in z:
        g -= 1 / (u - x)
    for x in rho:
        g += 1 / (u - x)
    return g


def root_dlog(P, rho, z, lo, hi):
    """The zero of d/du log|W f| in (lo, hi), by bisection (the derivative changes sign there)."""
    glo = dlog(P, rho, z, lo)
    for _ in range(4 * DPS):
        mid = (lo + hi) / 2
        if (dlog(P, rho, z, mid) > 0) == (glo > 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def polish(P, q0, z0, zs, ext0, inner):
    """Move the grid estimates of the stationary points onto the zeros of d/du log|W f| for the
    starting zeros, so that Newton starts close; needed for narrow bands, where the double-precision
    grid and bracketing are too coarse."""
    ext = list(ext0)
    nodes = [P['u1']] + sorted(q0) + [mp.mpf(1)]
    for j in inner:
        ext[j] = root_dlog(P, q0, z0, nodes[j] * (1 + mp.mpf(10) ** -30), nodes[j + 1] * (1 - mp.mpf(10) ** -30)) \
            if j > 0 or P['u1'] > 0 else ext[j]
    ys = []
    for b, zb in zip(P['bands'], zs):
        zb = sorted(mp.mpf(float(x)) for x in zb)
        d = mp.mpf(10) ** -30
        segs = list(zip(zb, zb[1:])) if b['up'] else list(zip(zb, zb[1:]))
        pts = [root_dlog(P, q0, z0, lo * (1 + d), hi * (1 - d)) for lo, hi in segs]
        if b['up']:
            lo = zb[-1] * (1 + d)
            hi = zb[-1] * 2
            while (dlog(P, q0, z0, hi) > 0) == (dlog(P, q0, z0, lo) > 0):
                hi *= 2
            pts.append(root_dlog(P, q0, z0, lo, hi))
        else:
            hi = zb[0] * (1 - d)
            lo = zb[0] / 2
            while (dlog(P, q0, z0, lo) > 0) == (dlog(P, q0, z0, hi) > 0):
                lo /= 2
            pts.insert(0, root_dlog(P, q0, z0, lo, hi))
        ys.append(pts)
    return ext, ys


def newton(P, zs, q0, ext0, ys0):
    mp.mp.dps = DPS
    d, r, a, u1 = P['d'], P['r'], P['a'], P['u1']
    ext0 = [mp.mpf(float(e)) for e in ext0]
    lo_end = abs(ext0[0] - u1) < 1e-9 and not (r > a and u1 == 0)
    inner = list(range(1 if lo_end else 0, d))
    nb = [b['m'] for b in P['bands']]

    def unpack(v):
        v = list(v)
        rho, k = v[:d], d
        xs = list(ext0)
        for j in inner:
            xs[j] = v[k]
            k += 1
        if lo_end:
            xs[0] = mp.mpf(u1)
        xs[d] = mp.mpf(1)
        h = v[k]
        k += 1
        zb, yb, lb = [], [], []
        for m in nb:
            zb.append(v[k:k + m])
            yb.append(v[k + m:k + 2 * m])
            lb.append(v[k + 2 * m])
            k += 2 * m + 1
        return rho, xs, h, zb, yb, lb

    def eqs(*v):
        rho, xs, h, zb, yb, lb = unpack(v)
        z = [x for zz in zb for x in zz]
        out = [S(P, rho, z, xs[j]) - (-1) ** (d - j) * h for j in range(d + 1)]
        out += [dlog(P, rho, z, xs[j]) for j in inner]
        for b, yy, lam in zip(P['bands'], yb, lb):
            out.append(mp.log(abs(S(P, rho, z, b['ue']))) - lam)
            for y in yy:
                out += [mp.log(abs(S(P, rho, z, y))) - lam, dlog(P, rho, z, y)]
        return out
    q0 = [mp.mpf(float(x)) for x in q0]
    z0 = [mp.mpf(float(x)) for zb in zs for x in zb]
    ext0, ys0 = polish(P, q0, z0, zs, ext0, inner)
    v0 = q0 + [ext0[j] for j in inner] + [abs(S(P, q0, z0, mp.mpf(1)))]
    for b, zb, yy in zip(P['bands'], zs, ys0):
        v0 += [mp.mpf(float(x)) for x in zb] + [mp.mpf(float(y)) for y in yy]
        v0.append(mp.log(abs(S(P, q0, z0, b['ue']))))
    # plain Newton with the Jacobian from mp.jacobian; mp.findroot's own step control stalls on
    # the wide band-pass filters, whose reflection zeros crowd the lower band edge
    v, best = mp.matrix(v0), None
    floor = mp.mpf(10) ** (-DPS + 12)
    for _ in range(40):
        F = mp.matrix(eqs(*v))
        nF = mp.norm(F)
        if best is not None and nF < floor and nF > best[0] / 10:
            break
        if best is None or nF < best[0]:
            best = (nF, v)
        v = v + mp.lu_solve(mp.jacobian(lambda *x: eqs(*x), list(v)), -F)
    assert best[0] < floor, ('Newton did not converge', mp.nstr(best[0], 5))
    sol = best[1]
    return unpack([sol[i] for i in range(len(v0))])


# --------------------------------------------------------- certificates ---

def golden(f, lo, hi, maximum):
    g = (mp.sqrt(5) - 1) / 2
    for _ in range(250):
        m1, m2 = hi - g * (hi - lo), lo + g * (hi - lo)
        if (f(m1) > f(m2)) == maximum:
            hi = m2
        else:
            lo = m1
    return (lo + hi) / 2


def certify_passband(P, rho, z):
    """Max of |W f| over the passband: its value, the number of points attaining it, alternation."""
    u1 = P['u1']
    N = 6000
    us = [u1 + (1 - u1) * (1 - mp.cos(mp.pi * k / N)) / 2 for k in range(N + 1)]
    if 0 < u1 < 1e-2:
        # wide band-pass: extrema crowd the lower edge closer than the grid above resolves
        us = sorted(set(us + [u1 * mp.power(1 / u1, mp.mpf(k) / N) for k in range(N + 1)]))
        N = len(us) - 1
    if P['r'] > P['a'] and u1 == 0:
        us[0] = (us[0] + us[1]) / 2
    f = lambda u: abs(S(P, rho, z, u))
    vals = [f(u) for u in us]
    peaks = []
    for k in range(N + 1):
        left = vals[k - 1] if k > 0 else -1
        right = vals[k + 1] if k < N else -1
        if vals[k] >= left and vals[k] >= right:
            peaks.append(us[k] if k in (0, N) else golden(f, us[k - 1], us[k + 1], True))
    pv = [f(u) for u in peaks]
    top = max(pv)
    at = [u for u, v in zip(peaks, pv) if v > top * (1 - mp.mpf(10) ** -30)]
    signs = [mp.sign(S(P, rho, z, u)) for u in at]
    return top, len(at), all(signs[i] != signs[i + 1] for i in range(len(signs) - 1))


def certify_stopband(P, rho, z, b, zb):
    """Min of |W f| over a stopband and the number of points attaining it."""
    lo, hi = (b['ue'], max(zb) * 100) if b['up'] else (min(zb) / 100, b['ue'])
    N = 4000
    us = [lo * (hi / lo) ** (mp.mpf(k) / N) for k in range(N + 1)]
    # for narrow bands the zeros sit within a fraction of a percent of the edge: add 4000 points
    # log-spaced between the edge and the outermost zero, and 300 on either side of every zero,
    # geometrically closer towards it
    far = max(zb) if b['up'] else min(zb)
    us += [b['ue'] * (far / b['ue']) ** (mp.mpf(k) / N) for k in range(1, N)]
    for x in zb:
        for t in np.geomspace(1e-9, 0.5, 300):
            us += [x * (1 + mp.mpf(t)), x * (1 - mp.mpf(t))]
    us = sorted(set(u for u in us if lo <= u <= hi and u not in zb))
    N = len(us) - 1
    f = lambda u: abs(S(P, rho, z, u))
    vals = [f(u) for u in us]
    dips = [b['ue']]
    for k in range(1, N):
        if vals[k] <= vals[k - 1] and vals[k] <= vals[k + 1]:
            dips.append(golden(f, us[k - 1], us[k + 1], False))
    dv = [f(u) for u in dips]
    low = min(dv)
    # the band minimum must not be undercut anywhere on the grid either
    assert min(vals) >= low * (1 - mp.mpf(10) ** -30)
    return low, sum(1 for v in dv if v < low * (1 + mp.mpf(10) ** -30))


def main():
    specs_dir, out = sys.argv[1], sys.argv[2]
    soldir = sys.argv[3] if len(sys.argv) > 3 else None
    if soldir:
        sys.path.insert(0, soldir)
        import synthesize as SY
    truth = {}
    for name in sorted(os.listdir(specs_dir)):
        if not name.endswith('.json'):
            continue
        spec = json.load(open(os.path.join(specs_dir, name)))
        f2 = mp.mpf(str(spec['passband_hz'][1]))
        P = setup(spec)
        zs, q0, ext0, ys0 = float_stage(P)
        rho, xs, h, zb, yb, lb = newton(P, zs, q0, ext0, ys0)
        z = [x for zz in zb for x in zz]
        top, nat, alt = certify_passband(P, rho, z)
        assert nat == P['d'] + 1 and alt, (spec['id'], nat, alt)
        assert all(P['u1'] < q < 1 for q in rho)
        c2 = P['eps2'] / top ** 2
        levels = []
        for b, zz in zip(P['bands'], zb):
            assert all(x > b['ue'] for x in zz) if b['up'] else all(0 < x < b['ue'] for x in zz)
            low, count = certify_stopband(P, rho, z, b, zz)
            assert count == b['m'] + 1, (spec['id'], count)
            levels.append(10 * mp.log10(1 + c2 * low ** 2))
        rec = {'degree': P['n'], 'r': P['r'], 'dc': P['a'],
               'reflection_zeros_u': [mp.nstr(q, 45) for q in sorted(rho)],
               'transmission_zeros_u': [mp.nstr(x, 45) for x in sorted(z)],
               'transmission_zeros_hz': [float(f2 * mp.sqrt(x)) for x in sorted(z)],
               'stopband_min_attenuation_db': [round(float(x), 6) for x in levels],
               'c2': mp.nstr(c2, 45), 'f2_hz': spec['passband_hz'][1]}
        line = '%-10s d=%d alternation points %d  zeros %s  A_min %s dB' % (
            spec['id'], P['d'], nat, ['%.6g' % x for x in rec['transmission_zeros_hz']],
            ['%.2f' % x for x in rec['stopband_min_attenuation_db']])
        if soldir:
            # every canonical ladder up to degree 12; above that the first 40 with an admissible load
            full = P['n'] <= 12
            SY.ADMISSIBLE_ONLY = not full
            ap, sols = SY.synthesize(spec) if full else SY.synthesize(spec, limit=40)
            dz = max([abs(mp.mpf(q) - p) / p for q, p in zip(sorted(rho) + sorted(z), sorted(ap.ur) + sorted(ap.uz))]
                     + [mp.mpf(0)])
            best = SY.choose(spec, sols)
            adm = [R for _, p, R in sols if 0.5 <= R <= 2]
            rec['witness'] = SY.denormalise(spec, best[1], best[2])
            rec['ladders'] = {'canonical': len(sols), 'admissible_load': len(adm)} if full else \
                {'canonical': None, 'admissible_load': None}
            line += '  |dzero| vs solution %.1e  ladders %s, admissible load %s' % (float(dz), rec['ladders']['canonical'], rec['ladders']['admissible_load'])
        mp.mp.dps = DPS
        print(line, flush=True)
        truth[spec['id']] = rec
    json.dump(truth, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main()
