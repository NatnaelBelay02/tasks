"""Reference solution: modes of lossy layer-loaded rectangular guides.

usage: python3 modes.py guides.json modes.json

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

Gyrotropic guides (y-biased ferrite, y-magnetised plasma): hybrid families and modal power.

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

Search region (family_modes): |Im gamma| <= B0 (a material scale) in strips, then shells [B, 2B], doubling
B until a shell is empty and B is past 20 / d_min of the thinnest ferrite or plasma layer: slow
quasi-static waves (magnetostatic surface waves of a film, surface plasmons and magnetoplasmons) have
|gamma| d of order one, so they can lie far beyond any material scale.  Families n = 1, 2, ... are
searched until one has no zero past k_y^2 > k0^2 max|eps mu| + A^2.

Power split (Hyb.power): the state is carried from the wall x = 0 across the guide in mpmath.  An error
delta in gamma grows like delta exp(|gamma| x), so the power split uses the polished gamma itself
(EXACT, 45 + |gamma| a / 2.3 digits), not its double-precision value, at 40 + growth / 2.3 digits, with
Gauss-Legendre quadrature on pieces of each layer over which the field grows by at most e^0.5.
"""
import cmath
import copy
import json
import math
import mpmath as mp
import numpy as np
import scipy.linalg as sl
import sys


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
    def _layer(self, g, d, gyro, ten, te):
        """exp(M d) up to a positive real factor: closed form in an isotropic layer, eigendecomposition
        in a gyrotropic one (scipy expm if the eigenvectors are ill-conditioned)."""
        M, _, _ = m_matrix(g, self.f, self.ky, ten)
        if te:
            M = M[np.ix_([0, 3], [0, 3])]
        if gyro:
            lam, V = np.linalg.eig(M * d)
            if np.linalg.cond(V) < 1e8:
                return (V * np.exp(lam - lam.real.max())) @ np.linalg.inv(V)
            m = max(1, int(math.ceil(np.abs(lam).max() / 2)))
            E1 = sl.expm(M * (d / m))
            E1 /= max(np.abs(E1).max(), 1e-300)
            return np.linalg.matrix_power(E1, m)
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
        return c * np.eye(M.shape[0]) + s * M

    def value(self, g):
        """D(gamma) up to a positive real factor.  The two state columns are re-orthonormalised after
        every layer (Gram-Schmidt: a unimodular column operation and positive scalings), which keeps
        them independent for |gamma| >> k0 and leaves arg D unchanged."""
        if self.n == 0:
            s = np.array([0, 1], complex)
            for L, ten in zip(self.G.layers, self.ten):
                s = self._layer(g, L[0], gyrotropic(L), ten, True) @ s
                s /= max(np.abs(s).max(), 1e-300)
            return s[0]
        S = np.zeros((4, 2), complex)
        S[2, 0] = S[3, 1] = 1
        for L, ten in zip(self.G.layers, self.ten):
            S = self._layer(g, L[0], gyrotropic(L), ten, False) @ S
            S[:, 0] /= np.linalg.norm(S[:, 0])
            S[:, 1] -= np.vdot(S[:, 0], S[:, 1]) * S[:, 0]
            S[:, 1] /= np.linalg.norm(S[:, 1])
        return S[0, 0] * S[1, 1] - S[0, 1] * S[1, 0]

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
        """Secant iteration in mpmath; converged when the step is below 1e-30 |gamma| (D itself can be
        of any size, so its value is not used as the stopping test)."""
        dps = 45 + int(abs(g) * self.G.a / 2.3)
        with mp.workdps(dps):
            x0 = mp.mpc(g.real, g.imag)
            r = mp.findroot(self.value_mp, x0, tol=mp.mpf(10) ** (-dps - 15), verify=False, maxsteps=80)
            r2 = mp.findroot(self.value_mp, r, tol=mp.mpf(10) ** (-dps - 15), verify=False, maxsteps=10)
            if abs(r2 - r) > mp.mpf(10) ** (20 - dps) * abs(r) or abs(r - x0) > 1e-3 * abs(x0):
                raise RuntimeError('polishing did not converge near %r in %s' % (g, self.G.id))
            EXACT[complex(r2)] = r2
            return complex(r2)

    # ---- modal power per layer ----
    def power(self, g, nodes=24):
        """Fractions P_i / P of the time-averaged power flow through the cross-section in each layer.
        Each layer is cut into pieces over which the field grows by at most e^0.5 (the integrand by e),
        with 24-point Gauss-Legendre quadrature on each piece; the working precision grows with the total
        growth of the field across the guide."""
        rates = []
        for ten in self.ten:
            M, _, _ = m_matrix(complex(g), self.f, self.ky, ten)
            rates.append(max(abs(np.linalg.eigvals(M).real).max(), 1.0))
        growth = sum(r * L[0] for r, L in zip(rates, self.G.layers))
        with mp.workdps(40 + int(growth / 2.3)):
            gm = EXACT.get(complex(g), mp.mpc(g.real, g.imag))
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
            for ML, rate in zip(self.G.mlayers, rates):
                dl = ML[0]
                m = max(1, int(math.ceil(rate * float(dl) / 0.5)))
                h = dl / m
                Eh, _, _, _ = self.layer_mp(ML, gm, h, ky)
                nodes_E = [self.layer_mp(ML, gm, h * (x + 1) / 2, ky) for x in xs]
                tot = mp.mpf(0)
                for _ in range(m):
                    for (E, M, ex, hx), wq in zip(nodes_E, ws):
                        p = E * s
                        exv = sum(ex[k] * p[k] for k in range(4))
                        hxv = sum(hx[k] * p[k] for k in range(4))
                        tot += wq * h / 2 * mp.re(exv * mp.conj(p[2]) - p[0] * mp.conj(hxv))
                    s = Eh * s
                P.append(tot)
            tot = sum(P)
            return [float(x / tot) for x in P]


_GAUSS = {}

# High-precision value of every polished propagation constant, keyed by its double-precision value.  The
# power split carries the field across the guide, where an error delta in gamma grows like
# delta exp(|gamma| x); it therefore uses this value, polished to 40 + |gamma| a / 2.3 digits.
EXACT = {}


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
def gwinding(H, corners):
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
    n = gwinding(H, [complex(r0, i0), complex(r1, i0), complex(r1, i1), complex(r0, i1)])
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
    """Material scale of |Im gamma|: 4 k0 sqrt(max |eps| |mu|) + 2 k_y + 200 rad/m."""
    k0 = math.sqrt(H.k02)
    m = 1.0
    for e1, e3, e2, mu, ka in H.ten:
        mue = (mu * mu - ka * ka) / mu
        m = max(m, max(abs(e1), abs(e3), abs(e2) + abs(e1)) * max(abs(mu), abs(mue), abs(ka) + abs(mu), 1))
    return 4 * k0 * math.sqrt(m) + 2 * H.ky + 200


def quasistatic_scale(G):
    """|Im gamma| beyond which every gyrotropic or plasma layer is at least 20 decay lengths thick
    (|gamma| d >= 20).  Slow quasi-static waves (magnetostatic surface waves of a ferrite film,
    surface plasmons and magnetoplasmons) have |gamma| d of order one; beyond this scale the interfaces
    are decoupled and a zero could only come from a single-interface resonance, which the empty last
    shell of the search excludes."""
    ds = [L[0] for L in G.layers if L[1] in ('f', 'm', 'p')]
    return 20.0 / min(ds) if ds else 0.0


def _strips(H, lo, hi, A1, nstrip):
    zs = []
    for sgn in (1, -1):
        for k in range(nstrip):
            b0 = lo + (hi - lo) * k / nstrip
            b1 = lo + (hi - lo) * (k + 1) / nstrip
            if sgn < 0:
                b0, b1 = -b1, -b0
            b0 += 1e-3
            b1 += 1e-3
            for r0, r1 in ((-A1, -1e-6), (1e-6, A1)):
                zs += zeros_rect(H, r0, r1, b0, b1)
    return zs


def family_modes(G, n, fs, A, A1):
    """(forward, reverse) gamma lists (alpha <= A) of family n.  |Im gamma| <= B0 (material scale) is
    searched in 24 strips per sign; then shells [B, 2B] in 16 strips each, doubling B until a shell has
    no zero and B is past the quasi-static scale of the guide."""
    H = Hyb(G, n, fs)
    B = bmax(H)
    zs = _strips(H, 0.0, B, A1, 24)
    Bq = quasistatic_scale(G)
    while True:
        new = _strips(H, B, 2 * B, A1, 16)
        zs += new
        B *= 2
        if not new and B >= Bq:
            break
    fwd = [z for z in zs if 0 < z.real <= A]
    rev = [-z for z in zs if -A <= z.real < 0]
    return H, fwd, rev


# ======================================================================= isotropic guides, main ===
DB = 20 / math.log(10)
A_DB = 200.0
A = A_DB / DB                      # Np/m
A1 = A * 1.02


EXACTQ = {}     # high-precision transverse eigenvalues q, keyed by their double values


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
        dps = 50 + int(abs(cmath.sqrt(q)) * self.a / 2.3)
        with mp.workdps(dps):
            fm = mp.mpf(str(f))
            r = mp.findroot(lambda x: self.dfun_mp(kind, x, fm), mp.mpc(q.real, q.imag),
                            tol=mp.mpf(10) ** (-dps - 15), verify=False, maxsteps=80)
            EXACTQ[complex(r)] = r
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
            H, fwd, rev = family_modes(G, n, fs, A, A1)
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
                        qx = EXACTQ.get(q)
                        if qx is not None:
                            with mp.workdps(mp.mp.dps + 60):
                                gx = mp.sqrt(mp.mpf(n) ** 2 * mp.pi ** 2 / mp.mpf(str(G.b_str)) ** 2 - qx)
                                if mp.re(gx) < 0:
                                    gx = -gx
                                EXACT[g] = gx
                if lst:
                    H = Hyb(G, n, fs)
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
