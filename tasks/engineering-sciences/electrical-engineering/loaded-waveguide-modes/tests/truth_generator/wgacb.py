"""Ball (arb/acb) evaluation of the transverse dispersion functions of a lossy layered guide.

D_LSE(q, f) = u(a) for u'' + (k0^2 eps(x, f) - q) u = 0, u(0) = 0, u'(0) = 1, u and u' continuous;
D_LSM(q, f) = w(a) for v'' + (k0^2 eps(x, f) - q) v = 0, w = v'/eps, v(0) = 1, w(0) = 0,
              v and w continuous;
k0 = 2 pi f / c, eps complex (dielectric eps'(1 - j tan d), collisional Drude plasma
1 - f_p^2 / (f (f - j nu / 2 pi))), time dependence exp(+j 2 pi f t).

Everything is carried as a first-order dual number (value, d/dq, d/df) with complex-ball
(acb) coefficients, so one evaluation on a ball gives enclosures of D, dD/dq and dD/df.
"""
from fractions import Fraction

import flint
from flint import acb, arb

flint.ctx.prec = 256
C_LIGHT = 299792458
flint.ctx.prec = 1100                 # constants at high precision (see make_truth: slow-wave refinement)
TWO_PI_C = 2 * arb.pi() / C_LIGHT
flint.ctx.prec = 256


def A(x):
    """Exact arb from a decimal string / Fraction / int."""
    if isinstance(x, arb):
        return x
    q = Fraction(str(x))
    return arb(q.numerator) / q.denominator


class Du:
    """value, d/dq, d/df (acb)."""
    __slots__ = ('v', 'q', 'f')

    def __init__(self, v, q=0, f=0):
        self.v = v if isinstance(v, acb) else acb(v)
        self.q = q if isinstance(q, acb) else acb(q)
        self.f = f if isinstance(f, acb) else acb(f)

    def __add__(s, o):
        if not isinstance(o, Du):
            return Du(s.v + o, s.q, s.f)
        return Du(s.v + o.v, s.q + o.q, s.f + o.f)

    __radd__ = __add__

    def __neg__(s):
        return Du(-s.v, -s.q, -s.f)

    def __sub__(s, o):
        return s + (-o)

    def __rsub__(s, o):
        return (-s) + o

    def __mul__(s, o):
        if not isinstance(o, Du):
            return Du(s.v * o, s.q * o, s.f * o)
        return Du(s.v * o.v, s.q * o.v + s.v * o.q, s.f * o.v + s.v * o.f)

    __rmul__ = __mul__

    def inv(s):
        r = 1 / s.v
        r2 = -r * r
        return Du(r, s.q * r2, s.f * r2)

    def __truediv__(s, o):
        if not isinstance(o, Du):
            return Du(s.v / o, s.q / o, s.f / o)
        return s * o.inv()


def _fact_table(n):
    t = [arb(1)]
    for k in range(1, n + 1):
        t.append(t[-1] * k)
    return t


_KS = 60
_FACT = _fact_table(2 * _KS + 4)


def _series(w, d):
    """C = sum w^k/(2k)!, S = d sum w^k/(2k+1)!, and their w-derivatives, |w| <= 16."""
    R = w.abs_upper()
    C = acb(0); S = acb(0); Cw = acb(0); Sw = acb(0)
    p = acb(1)
    for k in range(_KS):
        C += p / _FACT[2 * k]
        S += p / _FACT[2 * k + 1]
        if k + 1 < _KS:
            Cw += (k + 1) * p / _FACT[2 * k + 2]
            Sw += (k + 1) * p / _FACT[2 * k + 3]
        p *= w
    # remainder: sum_{k >= K} R^k / (2k)! <= 2 R^K / (2K)!  for R <= 16 (ratio < 1/2), same for the others
    tail = 2 * (R ** _KS) * _KS / _FACT[2 * _KS - 2]
    err = acb(arb(0, tail), arb(0, tail))
    return C + err, (S + err) * d, Cw + err, (Sw + err) * d


def layer_funcs(z, d):
    """C = cosh(k d), S = sinh(k d)/k with k^2 = -z (entire in z), and dC/dz, dS/dz.
    Returns None if no branch is safe for this ball (caller subdivides)."""
    w = -z * d * d
    if w.abs_upper() <= 16:
        C, S, Cw, Sw = _series(w, d)
        return C, S, -d * d * Cw, -d * d * Sw
    re, im = z.real, z.imag
    if (im > 0) or (im < 0) or (re < 0):
        # the principal sqrt(-z) is analytic on the ball: -z avoids (-inf, 0]
        k = (-z).sqrt()
        e1 = (k * d).exp()
        e2 = (-k * d).exp()
        C = (e1 + e2) / 2
        S = (e1 - e2) / (2 * k)
        return C, S, -d / 2 * S, (d * C - S) / (2 * z)
    if re > 0:
        # the principal sqrt(z) is analytic on the ball
        s = z.sqrt()
        C = (s * d).cos()
        S = (s * d).sin() / s
        return C, S, -d / 2 * S, (d * C - S) / (2 * z)
    return None


class Guide:
    def __init__(self, rec):
        self.id = rec['id']
        self.a = A(rec['a'])
        self.b = A(rec['b'])
        self.layers = []
        for L in rec['layers']:
            m = L['material']
            d = A(L['thickness'])
            if m['type'] == 'dielectric':
                self.layers.append((d, 'd', A(m['eps']), A(m['loss_tangent'])))
            elif m['type'] == 'plasma':
                self.layers.append((d, 'p', A(m['f_p']), A(m['collision_frequency'])))
            elif m['type'] == 'magnetized_plasma':
                # (d, 'm', f_p, nu, f_c, s): the TE (n = 0) problem sees only e3, the Drude permittivity
                sgn = 1 if m['bias_direction'] == '+y' else -1
                self.layers.append((d, 'm', A(m['f_p']), A(m['collision_frequency']), A(m['cyclotron_frequency']), sgn))
            else:
                sgn = 1 if m['bias_direction'] == '+y' else -1
                self.layers.append((d, 'f', A(m['eps']), A(m['loss_tangent']), A(m['saturation_magnetization_gauss']),
                                    A(m['bias_field_oe']), A(m['linewidth_oe']), sgn))
        self.ferrite = any(L[1] in ('f', 'm') for L in self.layers)      # gyrotropic layer present
        self.band = [A(x) for x in rec['band']] if 'band' in rec else None
        self.tests = [A(x) for x in rec['test_frequencies']]

    def _eps_point(self, L, f):
        d, kind, p1, p2 = L[:4]
        if kind in ('d', 'f'):
            return acb(p1, -p1 * p2), acb(0)
        c0 = acb(0, p2 / (2 * arb.pi()))             # j nu / 2 pi
        den = f * (f - c0)
        return 1 - p1 * p1 / den, p1 * p1 * (2 * f - c0) / (den * den)

    def eps(self, L, f):
        """eps and d eps / df at f (dual in f).  On a wide frequency ball the value is enclosed by
        the mean-value form eps(mid) + eps'(ball) (f - mid), which avoids the dependency blow-up."""
        e, de = self._eps_point(L, f)
        if L[1] == 'p' and not f.is_exact():
            m = acb(f.real.mid())
            em, _ = self._eps_point(L, m)
            e2 = em + de * (f - m)
            if e2.rad() < e.rad():
                e = e2
        return Du(e, 0, de)

    def k0eps(self, L, f):
        """k0^2 eps as a dual in f."""
        e = self.eps(L, f)
        k2 = TWO_PI_C * TWO_PI_C
        val = k2 * f * f * e.v
        der = k2 * (2 * f * e.v + f * f * e.f)
        if not f.is_exact():
            m = acb(f.real.mid())
            em = self.eps(L, m)
            v2 = k2 * m * m * em.v + der * (f - m)
            if v2.rad() < val.rad():
                val = v2
        return Du(val, 0, der)

    def D(self, kind, q, f, scaled=True, mask=None):
        """D_kind at complex ball q, real ball f, as Du (value, d/dq, d/df); None if a layer ball
        needs subdivision.  With scaled=True, every layer with Re(kappa d) > 1 on the whole ball
        (kappa = sqrt(q - k0^2 eps)) is divided by cosh(kappa d): the step becomes
        [[1, T/kappa], [kappa T, 1]], T = tanh(kappa d).  That factor is analytic and nonzero on
        the ball, so the zeros and their multiplicities are those of D, and the balls stay tight.
        mask: the set of scaled layers to use (taken from an enclosing ball); self.mask records
        the set used."""
        used = []
        f = acb(f) if not isinstance(f, acb) else f
        if kind == 'LSE':
            v, p = Du(0), Du(1)
        else:
            v, p = Du(1), Du(0)
        eprev = None
        for L in self.layers:
            d = L[0]
            e = self.eps(L, f)
            if kind == 'LSM' and eprev is not None:
                p = p * (e / eprev)
            ke = self.k0eps(L, f)
            z = Du(ke.v - q, acb(-1), ke.f)
            done = False
            want = mask[len(used)] if mask is not None else None
            if want or (want is None and scaled and ((z.v.real < 0) or (z.v.imag > 0) or (z.v.imag < 0))):
                kv = (-z.v).sqrt()
                if want or (kv.real * d) > 1:
                    dk = -1 / (2 * kv)                       # d kappa / dz
                    kap = Du(kv, dk * z.q, dk * z.f)
                    tv = (kv * d).tanh()
                    dt = (1 - tv * tv) * d                   # d T / d kappa
                    T = Du(tv, dt * kap.q, dt * kap.f)
                    v, p = v + (T / kap) * p, kap * T * v + p
                    done = True
            used.append(done)
            if not done:
                r = layer_funcs(z.v, d)
                if r is None:
                    return None
                C, S, Cz, Sz = r
                Cd = Du(C, Cz * z.q, Cz * z.f)
                Sd = Du(S, Sz * z.q, Sz * z.f)
                v, p = Cd * v + Sd * p, (-1) * z * Sd * v + Cd * p
            eprev = e
        self.mask = used
        if kind == 'LSE':
            return v
        return p / eprev


GYRO = A('2.8e6')        # gyromagnetic ratio, Hz per oersted


def polder(L, f):
    """(mu, kappa) of a y-biased ferrite layer at the exact frequency f (acb); identity for others.
    mu = 1 + f0 fm / (f0^2 - f^2), kappa = s f fm / (f0^2 - f^2), f0 = g H0 + j g dH / 2, fm = g 4 pi Ms,
    s = +1 for bias along +y; [mu] = [[mu, 0, -j kappa], [0, 1, 0], [j kappa, 0, mu]]."""
    if L[1] != 'f':
        return acb(1), acb(0)
    _, _, e, td, Ms, H0, dH, sgn = L
    f0 = acb(GYRO * H0, GYRO * dH / 2)
    fm = GYRO * Ms
    den = f0 * f0 - f * f
    return 1 + f0 * fm / den, sgn * f * fm / den


def DTE(G, g, f, mask=None):
    """y-independent TE modes (E = E_y(x) y, fields exp(-g z)) of a guide with ferrite layers:
    E'' + (k0^2 eps mu_e + g^2) E = 0 in each layer, mu_e = (mu^2 - kappa^2)/mu; E and
    P = (E' - j tau g E)/mu_e (tau = kappa/mu, P proportional to H_z) continuous; E(0) = 0, P(0) = 1;
    D = E(a).  Returns Du (value, d/dg, 0); scaled layers as in Guide.D."""
    f = acb(f) if not isinstance(f, acb) else f
    k2 = TWO_PI_C * TWO_PI_C * f * f
    E, P = Du(0), Du(1)
    gd = Du(g, acb(1), 0)
    used = []
    for L in G.layers:
        d = L[0]
        e = G._eps_point(L, f)[0]
        mu, ka = polder(L, f)
        mue = (mu * mu - ka * ka) / mu
        tau = ka / mu
        jt = acb(0, 1) * tau
        Ep = P * mue + gd * jt * E
        z = gd * gd + k2 * e * mue                      # K^2; the layer functions take z = K^2
        z = Du(z.v, z.q, acb(0))
        done = False
        want = mask[len(used)] if mask is not None else None
        if want or (want is None and ((z.v.real < 0) or (z.v.imag > 0) or (z.v.imag < 0))):
            kv = (-z.v).sqrt()
            if want or (kv.real * d) > 1:
                dk = -1 / (2 * kv)
                kap = Du(kv, dk * z.q, 0)
                tv = (kv * d).tanh()
                T = Du(tv, (1 - tv * tv) * d * kap.q, 0)
                E, Ep = E + (T / kap) * Ep, kap * T * E + Ep
                done = True
        used.append(done)
        if not done:
            r = layer_funcs(z.v, d)
            if r is None:
                return None
            C, S, Cz, Sz = r
            Cd = Du(C, Cz * z.q, 0)
            Sd = Du(S, Sz * z.q, 0)
            E, Ep = Cd * E + Sd * Ep, (-1) * z * Sd * E + Cd * Ep
        P = (Ep - gd * jt * E) / mue
    G.mask = used
    return E


def acb_box(r0, r1, i0, i1):
    return acb(arb((r0 + r1) / 2, (r1 - r0) / 2), arb((i0 + i1) / 2, (i1 - i0) / 2))
