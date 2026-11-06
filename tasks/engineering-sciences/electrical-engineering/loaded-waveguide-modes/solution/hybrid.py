"""Gyrotropic layers (y-biased ferrite, y-magnetised plasma): hybrid families and modal power.

Fields of family n vary across the height as
    E_x, E_z, H_y ~ sin(k_y y),   E_y, H_x, H_z ~ cos(k_y y),   k_y = n pi / b,
and along the guide as exp(-gamma z).  With the permeability [[mu, 0, -j ka], [0, 1, 0], [j ka, 0, mu]]
and the permittivity [[e1, 0, -j e2], [0, e3, 0], [j e2, 0, e1]] of a layer, the x components follow
algebraically from the x components of the curl equations,
    e_x = (gamma h_y - k_y h_z + j u e2 e_z) / (u e1),   h_x = (-(k_y e_z + gamma e_y) / w + j ka h_z) / mu,
(w = j omega mu0, u = j omega eps0), and the tangential state psi = (e_y, e_z, h_y, h_z) obeys
psi' = M(gamma) psi.  psi is continuous at every interface; e_y = e_z = 0 on both side walls.  So a
family-n mode is a zero of D(gamma) = det of the e-rows of T(a) [0 0; 0 0; 1 0; 0 1], T(a) the product
of the layer transfer matrices exp(M d).  In an isotropic layer M^2 = -K^2 I (K^2 = k0^2 eps + gamma^2 -
k_y^2), so exp(M d) = cos(K d) I + sin(K d) / K M.  The family n = 0 is the TE problem (e_y, h_z).

D depends on gamma itself, not on gamma^2: zeros with Re gamma > 0 are forward modes, zeros with
Re gamma < 0 reverse modes (reported as -gamma).
"""
import cmath
import math

import mpmath as mp
import numpy as np
import scipy.linalg as sl

C = 299792458.0
MU0 = 4e-7 * math.pi            # only the products w u = -omega^2 / c^2 and w / u enter the zeros
GYRO = 2.8e6


def tensors(L, f, lib=cmath):
    """(e1, e3, e2, mu, ka) of a layer at frequency f; see the module docstring for the tensor forms."""
    kind = L[1]
    pi = mp.pi if lib is mp else math.pi
    j = mp.mpc(0, 1) if lib is mp else 1j
    if kind in ('d', 'f'):
        e = L[2] * (1 - j * L[3])
        if kind == 'd':
            return e, e, 0, 1, 0
        _, _, _, _, Ms, H0, dH, sg = L
        g = mp.mpf('2.8e6') if lib is mp else GYRO
        f0 = g * H0 + j * g * dH / 2
        fm = g * Ms
        den = f0 * f0 - f * f
        return e, e, 0, 1 + f0 * fm / den, sg * f * fm / den
    if kind == 'p':
        e = 1 - L[2] * L[2] / (f * (f - j * L[3] / (2 * pi)))
        return e, e, 0, 1, 0
    # 'm': magnetised plasma (d, 'm', f_p, nu, f_c, s)
    _, _, fp, nu, fc, sg = L
    fn = f - j * nu / (2 * pi)
    den = f * (fn * fn - fc * fc)
    e1 = 1 - fp * fp * fn / den
    e3 = 1 - fp * fp / (f * fn)
    gg = fp * fp * fc / den
    # [eps] = [[e1, 0, j s g], [0, e3, 0], [-j s g, 0, e1]]  ->  e2 = -s g in the [[e1,0,-j e2],...] form
    return e1, e3, -sg * gg, 1, 0


def gyrotropic(L):
    return L[1] in ('f', 'm')


def m_matrix(g, f, ky, ten, lib=np):
    e1, e3, e2, mu, ka = ten
    if lib is mp:
        two_pi = 2 * mp.pi
        w = mp.mpc(0, 1) * two_pi * f * mp.mpf(4) * mp.pi * mp.mpf(10) ** -7
        u = mp.mpc(0, 1) * two_pi * f / (mp.mpf(4) * mp.pi * mp.mpf(10) ** -7 * mp.mpf(C) ** 2)
        J = mp.mpc(0, 1)
        Mx = mp.matrix(4, 4)
    else:
        w = 2j * math.pi * f * MU0
        u = 2j * math.pi * f / (MU0 * C * C)
        J = 1j
        Mx = np.zeros((4, 4), complex)
    ex = [0, J * u * e2 / (u * e1), g / (u * e1), -ky / (u * e1)]
    hx = [-g / w / mu, -ky / w / mu, 0, J * ka / mu]
    for c in range(4):
        Mx[1, c] = -g * ex[c]
        Mx[0, c] = ky * ex[c] - w * J * ka * hx[c]
        Mx[3, c] = -g * hx[c]
        Mx[2, c] = -ky * hx[c] + u * J * e2 * ex[c]
    Mx[1, 2] += w
    Mx[0, 3] -= w * mu
    Mx[3, 0] -= u * e3
    Mx[2, 1] += u * e1
    return Mx, ex, hx


class Hyb:
    """Double-precision and mpmath evaluation of D for one guide, family n and frequency."""

    def __init__(self, G, n, fs):
        self.G = G
        self.n = n
        self.f = float(fs)
        self.fm = mp.mpf(str(fs))
        self.ky = n * math.pi / G.b
        self.k02 = (2 * math.pi * self.f / C) ** 2
        self.ten = [tensors(L, self.f) for L in G.layers]

    # ---- value up to a positive real factor (for the argument principle) ----
    def value(self, g):
        if self.n == 0:
            return self._te(g)
        S = np.zeros((4, 2), complex)
        S[2, 0] = S[3, 1] = 1
        for L, ten in zip(self.G.layers, self.ten):
            d = L[0]
            M, _, _ = m_matrix(g, self.f, self.ky, ten)
            if gyrotropic(L):
                E = sl.expm(M * d)
            else:
                K2 = self.k02 * ten[0] + g * g - self.ky ** 2
                kap = cmath.sqrt(-K2)
                if kap.real < 0:
                    kap = -kap
                if abs(kap * d) < 1e-8:
                    c, s = 1.0, d
                else:
                    e1 = cmath.exp(kap * d - kap.real * d)
                    e2 = cmath.exp(-kap * d - kap.real * d)
                    c, s = (e1 + e2) / 2, (e1 - e2) / (2 * kap)
                E = c * np.eye(4) + s * M
            S = E @ S
            S /= max(np.abs(S).max(), 1e-300)
        return S[0, 0] * S[1, 1] - S[0, 1] * S[1, 0]

    def _te(self, g):
        """n = 0: state (e_y, h_z) of the TE problem."""
        s = np.array([0, 1], complex)
        for L, ten in zip(self.G.layers, self.ten):
            d = L[0]
            M, _, _ = m_matrix(g, self.f, 0.0, ten)
            T = M[np.ix_([0, 3], [0, 3])]
            s = sl.expm(T * d) @ s
            s /= max(np.abs(s).max(), 1e-300)
        return s[0]

    # ---- mpmath ----
    def layer_mp(self, ML, g, d, ky):
        """(exp(M d), M, e_x row, h_x row) of the mpmath layer ML at gamma = g."""
        ten = tensors(ML, self.fm, lib=mp)
        M, ex, hx = m_matrix(g, self.fm, ky, ten, lib=mp)
        if gyrotropic(ML):
            return mp.expm(M * d), M, ex, hx
        k02 = (2 * mp.pi * self.fm / C) ** 2
        K2 = k02 * ten[0] + g * g - ky * ky
        K = mp.sqrt(K2)
        c = mp.cos(K * d)
        s = mp.sin(K * d) / K if K != 0 else d
        return c * mp.eye(4) + s * M, M, ex, hx

    def value_mp(self, g):
        ky = self.n * mp.pi / mp.mpf(str(self.G.b_str))
        if self.n == 0:
            s = mp.matrix([[0], [0], [0], [1]])
        else:
            s = mp.matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
        for ML in self.G.mlayers:
            E, _, _, _ = self.layer_mp(ML, g, ML[0], ky)
            s = E * s
        if self.n == 0:
            return s[0, 0]
        return s[0, 0] * s[1, 1] - s[0, 1] * s[1, 0]

    def polish(self, g):
        with mp.workdps(45):
            r = mp.findroot(self.value_mp, mp.mpc(g.real, g.imag), tol=mp.mpf(10) ** -32)
            return complex(r)

    # ---- modal power per layer ----
    def power(self, g, nodes=24):
        """Fractions P_i / P of the time-averaged power flow through the cross-section in each layer."""
        with mp.workdps(40):
            gm = mp.mpc(g.real, g.imag)
            ky = self.n * mp.pi / mp.mpf(str(self.G.b_str))
            if self.n == 0:
                s = mp.matrix([[0], [0], [0], [1]])
            else:
                S = mp.matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
                for ML in self.G.mlayers:
                    E, _, _, _ = self.layer_mp(ML, gm, ML[0], ky)
                    S = E * S
                # null vector of the 2x2 e-block
                a, b, c, d = S[0, 0], S[0, 1], S[1, 0], S[1, 1]
                v = mp.matrix([[-b], [a]]) if abs(a) + abs(b) >= abs(c) + abs(d) else mp.matrix([[-d], [c]])
                s = mp.matrix([[0], [0], [v[0]], [v[1]]])
            xs, ws = gauss(nodes)
            P = []
            for ML in self.G.mlayers:
                dl = ML[0]
                tot = mp.mpf(0)
                for x, wq in zip(xs, ws):
                    t = dl * (x + 1) / 2
                    E, M, ex, hx = self.layer_mp(ML, gm, t, ky)
                    p = E * s
                    exv = sum(ex[k] * p[k] for k in range(4))
                    hxv = sum(hx[k] * p[k] for k in range(4))
                    tot += wq * dl / 2 * mp.re(exv * mp.conj(p[2]) - p[0] * mp.conj(hxv))
                P.append(tot)
                E, _, _, _ = self.layer_mp(ML, gm, dl, ky)
                s = E * s
            tot = sum(P)
            return [float(x / tot) for x in P]


_GAUSS = {}


def gauss(n):
    if n not in _GAUSS:
        x, w = np.polynomial.legendre.leggauss(n)
        _GAUSS[n] = ([mp.mpf(float(v)) for v in x], [mp.mpf(float(v)) for v in w])
        # refine nodes in mpmath (Newton on P_n) for 40-digit quadrature
        xs = []
        for x0 in _GAUSS[n][0]:
            xr = mp.findroot(lambda t: mp.legendre(n, t), x0)
            xs.append(xr)
        ws = [2 / ((1 - x * x) * mp.diff(lambda t: mp.legendre(n, t), x) ** 2) for x in xs]
        _GAUSS[n] = (xs, ws)
    return _GAUSS[n]


# ------------------------------------------------------------ argument principle in gamma ---
def winding(H, corners):
    total = 0.0
    for a, b in zip(corners, corners[1:] + corners[:1]):
        t, vt, h = 0.0, H.value(a), 1 / 64
        while t < 1:
            h = min(h, 1 - t)
            v1 = H.value(a + (b - a) * (t + h))
            dphi = cmath.phase(v1 / vt)
            if abs(dphi) > 0.3 and h > 1e-13:
                h /= 2
                continue
            total += dphi
            t += h
            vt = v1
            if abs(dphi) < 0.08:
                h *= 1.5
    return int(round(total / (2 * math.pi)))


def zeros_rect(H, r0, r1, i0, i1, depth=0):
    n = winding(H, [complex(r0, i0), complex(r1, i0), complex(r1, i1), complex(r0, i1)])
    if n == 0:
        return []
    size = max(r1 - r0, i1 - i0)
    if n == 1 and size < 1e-4 * max(abs(complex(r0, i0)), abs(complex(r1, i1)), 1.0):
        return [H.polish(complex((r0 + r1) / 2, (i0 + i1) / 2))]
    if depth > 120:
        raise RuntimeError('cannot isolate zeros in %s' % H.G.id)
    if (r1 - r0) >= (i1 - i0):
        m = r0 + 0.4999917 * (r1 - r0)          # off-centre split: never lands on a symmetry line
        return zeros_rect(H, r0, m, i0, i1, depth + 1) + zeros_rect(H, m, r1, i0, i1, depth + 1)
    m = i0 + 0.5000123 * (i1 - i0)
    return zeros_rect(H, r0, r1, i0, m, depth + 1) + zeros_rect(H, r0, r1, m, i1, depth + 1)


def bmax(H):
    k0 = math.sqrt(H.k02)
    m = 1.0
    for e1, e3, e2, mu, ka in H.ten:
        mue = (mu * mu - ka * ka) / mu
        m = max(m, max(abs(e1), abs(e3), abs(e2) + abs(e1)) * max(abs(mu), abs(mue), abs(ka) + abs(mu), 1))
    return 4 * k0 * math.sqrt(m) + 2 * H.ky + 200


def family_modes(G, n, fs, A, A1):
    """(forward, reverse) gamma lists (alpha <= A) of family n."""
    H = Hyb(G, n, fs)
    B = bmax(H)
    zs = []
    nstrip = 48
    for k in range(nstrip):
        b0 = -B + 2 * B * k / nstrip + 1e-3
        b1 = -B + 2 * B * (k + 1) / nstrip + 1e-3
        for r0, r1 in ((-A1, -1e-6), (1e-6, A1)):
            zs += zeros_rect(H, r0, r1, b0, b1)
    fwd = [z for z in zs if 0 < z.real <= A]
    rev = [-z for z in zs if -A <= z.real < 0]
    return H, fwd, rev
