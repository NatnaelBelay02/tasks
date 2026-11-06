"""Reference solution: modes of lossy layer-loaded rectangular guides.

usage: python3 solve.py [guides.json [modes.json]]

For each guide, test frequency and family kind (LSE, LSM):
  1. the zeros q of D(q, f) are counted with the argument principle (Delves and Lyness) on
     rectangles covering the region where some family has attenuation <= 200 dB/m (2 % margin):
     Re q >= -A^2 + (Im q)^2 / (4 A^2), Re q <= Q.  arg D is followed along each edge with
     adaptive steps; the transfer-matrix state is renormalised by positive real factors after
     every layer, which leaves arg D unchanged and avoids overflow;
  2. a rectangle with more than one zero is split; a rectangle with one zero gives its zero by
     Newton's method on D'/D (the log-derivative is carried with the state and is also
     unaffected by the renormalisation), started from the first moment of the contour integral;
  3. each zero is polished in mpmath (50 digits) and mapped to gamma = sqrt(k_n^2 - q) with
     Re gamma > 0 for every family n; modes with attenuation <= 200 dB/m are reported.
The LSE region ends at the Rayleigh bound Re q <= k0^2 max Re eps.  For LSM the region is
extended until D_LSM shows the large-q behaviour of evanescent layers on a whole strip
(no zero in a strip four times longer than the last zero found, and beyond every interface
surface-wave estimate k0^2 eps_i eps_j / (eps_i + eps_j)).
"""
import cmath
import copy
import json
import math
import sys

import mpmath as mp
import numpy as np

import hybrid

C = 299792458.0
DB = 20 / math.log(10)
A_DB = 200.0
A = A_DB / DB                      # Np/m
A1 = A * 1.02


class Guide:
    def __init__(self, rec):
        self.id = rec['id']
        self.a = float(rec['a'])
        self.b = float(rec['b'])
        self.layers = []
        self.mlayers = []
        for L in rec['layers']:
            m = L['material']
            d = L['thickness']
            if m['type'] == 'dielectric':
                self.layers.append((float(d), 'd', float(m['eps']), float(m['loss_tangent'])))
                self.mlayers.append((mp.mpf(str(d)), 'd', mp.mpf(str(m['eps'])), mp.mpf(str(m['loss_tangent']))))
            elif m['type'] == 'plasma':
                self.layers.append((float(d), 'p', float(m['f_p']), float(m['collision_frequency'])))
                self.mlayers.append((mp.mpf(str(d)), 'p', mp.mpf(str(m['f_p'])), mp.mpf(str(m['collision_frequency']))))
            elif m['type'] == 'magnetized_plasma':
                sg = 1 if m['bias_direction'] == '+y' else -1
                vals = [m['f_p'], m['collision_frequency'], m['cyclotron_frequency']]
                self.layers.append((float(d), 'm') + tuple(float(v) for v in vals) + (sg,))
                self.mlayers.append((mp.mpf(str(d)), 'm') + tuple(mp.mpf(str(v)) for v in vals) + (sg,))
            else:
                sg = 1 if m['bias_direction'] == '+y' else -1
                vals = [m['eps'], m['loss_tangent'], m['saturation_magnetization_gauss'], m['bias_field_oe'], m['linewidth_oe']]
                self.layers.append((float(d), 'f') + tuple(float(v) for v in vals) + (sg,))
                self.mlayers.append((mp.mpf(str(d)), 'f') + tuple(mp.mpf(str(v)) for v in vals) + (sg,))
        self.gyro = any(L[1] in ('f', 'm') for L in self.layers)
        self.b_str = rec['b']
        self.tests = rec['test_frequencies']

    @staticmethod
    def eps_of(L, f, lib=cmath):
        d, kind, p1, p2 = L[:4]
        if kind in ('d', 'f'):
            return p1 * (1 - 1j * p2)
        pi = mp.pi if lib is mp else math.pi
        return 1 - p1 * p1 / (f * (f - 1j * p2 / (2 * pi)))

    # ---- double precision: value up to a positive factor, and D'/D ----
    def logd(self, kind, q, f):
        """(D / s, D' / D) for a positive real s."""
        k02 = (2 * math.pi * f / C) ** 2
        if kind == 'LSE':
            v, p, dv, dp = 0j, 1 + 0j, 0j, 0j
        else:
            v, p, dv, dp = 1 + 0j, 0j, 0j, 0j
        ep = None
        for L in self.layers:
            d = L[0]
            e = self.eps_of(L, f)
            if kind == 'LSM' and ep is not None:
                r = e / ep
                p, dp = p * r, dp * r
            z = k02 * e - q
            kap = cmath.sqrt(-z)
            if kap.real < 0:
                kap = -kap
            if abs(kap * d) < 1e-3:
                w = -z * d * d
                Cc = 1 + w / 2 + w * w / 24
                S = d * (1 + w / 6 + w * w / 120)
                Cz = -d * d * (0.5 + w / 12)
                Sz = -d ** 3 * (1 / 6 + w / 60)
                scale = 1.0
            else:
                # scaled by exp(-Re kappa d) (a positive real factor)
                e1 = cmath.exp(kap * d - kap.real * d)
                e2 = cmath.exp(-kap * d - kap.real * d)
                Cc = (e1 + e2) / 2
                S = (e1 - e2) / (2 * kap)
                Cz = -d / 2 * S
                Sz = (d * Cc - S) / (2 * z)
                scale = 1.0
            # d/dq of z is -1
            nv = Cc * v + S * p
            np_ = -z * S * v + Cc * p
            ndv = Cc * dv + S * dp - Cz * v - Sz * p
            ndp = -z * S * dv + Cc * dp + (S + z * Sz) * v - Cz * p
            m = max(abs(nv), abs(np_), 1e-300)
            v, p, dv, dp = nv / m, np_ / m, ndv / m, ndp / m
            ep = e
        if kind == 'LSE':
            return v, (dv / v if v != 0 else complex('inf'))
        return p / ep, (dp / p if p != 0 else complex('inf'))

    # ---- mpmath polish ----
    def dfun_mp(self, kind, q, f):
        k02 = (2 * mp.pi * f / C) ** 2
        v, p = (mp.mpc(0), mp.mpc(1)) if kind == 'LSE' else (mp.mpc(1), mp.mpc(0))
        ep = None
        for L in self.mlayers:
            d = L[0]
            e = self.eps_of(L, f, lib=mp)
            if kind == 'LSM' and ep is not None:
                p = p * e / ep
            z = k02 * e - q
            s = mp.sqrt(z)
            Cc = mp.cos(s * d)
            S = mp.sin(s * d) / s if s != 0 else d
            v, p = Cc * v + S * p, -z * S * v + Cc * p
            ep = e
        return v if kind == 'LSE' else p / ep

    def polish(self, kind, q, f):
        with mp.workdps(50):
            fm = mp.mpf(str(f))
            r = mp.findroot(lambda x: self.dfun_mp(kind, x, fm), mp.mpc(q.real, q.imag), tol=mp.mpf(10) ** -40)
            return complex(r)


# ------------------------------------------------------------ argument principle ---
def winding(G, kind, f, corners, maxpts=200000):
    """Winding number of D around the rectangle (counter-clockwise corners), adaptive."""
    total = 0.0
    npts = 0
    for a, b in zip(corners, corners[1:] + corners[:1]):
        v0, _ = G.logd(kind, a, f)
        if v0 == 0:
            return None
        # walk the edge with step control on the phase change
        t, vt = 0.0, v0
        h = 1 / 64
        while t < 1:
            h = min(h, 1 - t)
            z1 = a + (b - a) * (t + h)
            v1, _ = G.logd(kind, z1, f)
            npts += 1
            if v1 == 0 or npts > maxpts:
                return None
            dphi = cmath.phase(v1 / vt)
            if abs(dphi) > 0.4 and h > 1e-12:
                h /= 2
                continue
            total += dphi
            t += h
            vt = v1
            if abs(dphi) < 0.1:
                h *= 1.5
    return int(round(total / (2 * math.pi))), total


def rect_zero(G, kind, f, r0, r1, i0, i1):
    """The single zero in a rectangle: first moment of the contour integral, then Newton on D'/D."""
    c = complex((r0 + r1) / 2, (i0 + i1) / 2)
    # contour moment sum z * (D'/D) dz / (2 pi i), trapezoid on 4 x 400 points
    corners = [complex(r0, i0), complex(r1, i0), complex(r1, i1), complex(r0, i1)]
    num = 0j
    N = 400
    for a, b in zip(corners, corners[1:] + corners[:1]):
        for k in range(N):
            z = a + (b - a) * (k + 0.5) / N
            _, ld = G.logd(kind, z, f)
            num += z * ld * (b - a) / N
    z = num / (2j * math.pi)
    if not (r0 <= z.real <= r1 and i0 <= z.imag <= i1):
        z = c
    for _ in range(60):
        _, ld = G.logd(kind, z, f)
        if ld == 0 or not cmath.isfinite(ld):
            break
        step = 1 / ld
        z -= step
        if abs(step) < 1e-14 * max(abs(z), 1):
            break
    return z


def zeros_in_strip(G, kind, f, r0, r1, depth=0):
    """All zeros with r0 <= Re q <= r1 inside the attenuation parabola."""
    h = 2 * A1 * math.sqrt(max(r1, 0) + A1 * A1)
    return zeros_in_rect(G, kind, f, r0, r1, -h, h, depth)


def zeros_in_rect(G, kind, f, r0, r1, i0, i1, depth=0):
    # drop rectangles entirely outside the parabola
    if r1 < -A1 * A1:
        return []
    m = 0.0 if i0 <= 0 <= i1 else min(abs(i0), abs(i1))
    if m * m / (4 * A1 * A1) - A1 * A1 > r1:
        return []
    corners = [complex(r0, i0), complex(r1, i0), complex(r1, i1), complex(r0, i1)]
    w = winding(G, kind, f, corners)
    if w is None:
        # a zero sits on the contour: nudge the rectangle
        dr, di = (r1 - r0) * 1e-3, (i1 - i0) * 1e-3
        return zeros_in_rect(G, kind, f, r0 - dr, r1 + dr, i0 - di, i1 + di, depth + 1)
    n = w[0]
    if n == 0:
        return []
    if n == 1 and depth > 0:
        return [rect_zero(G, kind, f, r0, r1, i0, i1)]
    if depth > 60:
        raise RuntimeError('cannot separate zeros in %s %s' % (G.id, kind))
    if (r1 - r0) >= (i1 - i0):
        mid = (r0 + r1) / 2
        return zeros_in_rect(G, kind, f, r0, mid, i0, i1, depth + 1) + zeros_in_rect(G, kind, f, mid, r1, i0, i1, depth + 1)
    mid = (i0 + i1) / 2
    return zeros_in_rect(G, kind, f, r0, r1, i0, mid, depth + 1) + zeros_in_rect(G, kind, f, r0, r1, mid, i1, depth + 1)


def all_zeros(G, kind, f):
    k02 = (2 * math.pi * f / C) ** 2
    eps = [Guide.eps_of(L, f) for L in G.layers]
    if kind == 'LSE':
        qmax = k02 * max(e.real for e in eps) * 1.0001 + 1
    else:
        est = [k02 * max(e.real for e in eps)]
        for x, y in zip(eps, eps[1:]):
            if (x + y) != 0:
                est.append(abs(k02 * x * y / (x + y)))
        qmax = 4 * max(est)
    out = []
    r = -A1 * A1 * 1.0001
    last_zero = 0.0
    while True:
        width = max(k02 / 4, abs(r) / 16, 1.0)
        hi = r + width
        zs = zeros_in_strip(G, kind, f, r, hi)
        out += zs
        if zs:
            last_zero = max(last_zero, max(z.real for z in zs))
        r = hi
        if r >= qmax and (kind == 'LSE' or r >= 4 * max(last_zero, k02)):
            break
    return out


def gyro_guide(G, rec):
    """Gyrotropic guide: LSE0 (n = 0, TE) and the hybrid families HYB<n>, forward and reverse differ."""
    res = []
    for fs in rec['test_frequencies']:
        entry = {'forward': {}, 'reverse': {}}
        n = 0
        while True:
            H, fwd, rev = hybrid.family_modes(G, n, fs, A, A1)
            name = 'LSE0' if n == 0 else 'HYB%d' % n
            for dr, lst, sign in (('forward', fwd, 1), ('reverse', rev, -1)):
                if lst:
                    entry[dr][name] = [[g.imag, g.real, H.power(sign * g)] for g in lst]
            # past the last family that can have alpha <= A: k_y^2 beyond k0^2 max|eps mu| + A^2
            if n >= 1 and not fwd and not rev and H.ky ** 2 > H.k02 * hybrid_reach(H) + A * A:
                break
            n += 1
        res.append(entry)
        print(G.id, fs, {d: {k: len(v) for k, v in e.items()} for d, e in entry.items()}, flush=True)
    return {'modes': res}


def hybrid_reach(H):
    m = 1.0
    for e1, e3, e2, mu, ka in H.ten:
        mue = (mu * mu - ka * ka) / mu
        m = max(m, max(abs(e1), abs(e3)) * max(abs(mu), abs(mue), 1))
    return m


def solve_guide(rec):
    G = Guide(rec)
    if G.gyro:
        return gyro_guide(G, rec)
    res = []
    for fs in rec['test_frequencies']:
        f = float(fs)
        entry = {}
        for kind in ('LSE', 'LSM'):
            zs = [G.polish(kind, z, fs) for z in all_zeros(G, kind, f)]
            uniq = []
            for z in zs:
                if all(abs(z - w) > 1e-9 * max(abs(z), 1) for w in uniq):
                    uniq.append(z)
            n = 0 if kind == 'LSE' else 1
            qmax = max([q.real for q in uniq] + [0.0])
            # Re gamma^2 = k_n^2 - Re q; once k_n^2 - A^2 > max Re q every alpha exceeds A
            while (n * math.pi / G.b) ** 2 - A * A <= qmax:
                kn2 = (n * math.pi / G.b) ** 2
                lst = []
                for q in uniq:
                    g = cmath.sqrt(kn2 - q)
                    if g.real < 0:
                        g = -g
                    if g.real <= A:
                        lst.append(g)
                if lst:
                    H = hybrid.Hyb(G, n, fs)
                    entry['%s%d' % (kind, n)] = [[g.imag, g.real, H.power(g)] for g in lst]
                n += 1
        # reciprocal guide: the modes towards -z have the same gamma and the same power split
        res.append({'forward': entry, 'reverse': copy.deepcopy(entry)})
        print(G.id, fs, {k: len(v) for k, v in entry.items()}, flush=True)
    return {'modes': res}


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else '/root/data/guides.json'
    dst = sys.argv[2] if len(sys.argv) > 2 else '/root/results/modes.json'
    guides = json.load(open(src))['guides']
    out = {}
    for rec in guides:
        out[rec['id']] = solve_guide(rec)
    with open(dst, 'w') as fh:
        json.dump({'guides': out}, fh)


if __name__ == '__main__':
    main()
