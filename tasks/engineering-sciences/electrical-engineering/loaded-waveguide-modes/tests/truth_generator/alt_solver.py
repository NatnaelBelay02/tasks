"""Second, separately written solver, used to score a correct submission that does not come from
the oracle or from the certified pipeline.  Shipped as a program: run as /root/results/modes.py it
passes the verifier.

usage: python3 alt_solver.py guides.json out.json [guide_id ...]

It never evaluates the transfer-matrix dispersion function on a contour and uses no interval
arithmetic:
  * multi-domain Chebyshev collocation of v'' + k0^2 eps(x) v = q v (complex eps), with the wall
    and interface conditions as algebraic rows, gives a generalised eigenproblem A x = q B x;
    the discretisation is solved twice (N and 1.5 N points per layer) and an eigenvalue is
    kept if the two agree to 1e-4 (relative);
  * each kept eigenvalue is polished by secant iteration on a determinant built from 2 x 2 matrix
    exponentials in mpmath (50 + |gamma| a / 2.3 digits); convergence is judged by the step, since
    |D| grows like exp(|gamma| a);
  * gamma = sqrt(k_n^2 - q) with Re gamma > 0 for every family n; modes with attenuation
    <= 200 dB/m are reported.
Guides with gyrotropic layers (ferrite, magnetised plasma): family LSE0 from a quadratic pencil of
the TE problem; the hybrid families HYB<n> from a six-field collocation of the curl equations, which
is linear in gamma (A x = gamma B x); every eigenvalue in the region is polished on a determinant of
4 x 4 matrix exponentials and kept only if the polishing converges near it.  The grids are sized for a
target |beta| (2 s sqrt((beta_t + k_y) d) + 16 points per layer, which resolves the exp(-kappa x) boundary layers
on Chebyshev nodes): beta_t starts at the quasi-static scale 3 / d of the thinnest ferrite or plasma
layer (slow magnetostatic and electrostatic surface waves have |beta| d of order one) and at twice the
largest |beta| of the previous family, and is doubled while a zero lies beyond beta_t / 2.
Power split (every family, every guide): at the polished gamma, kept to high precision, the tangential
state is carried through the guide in mpmath (60 + |gamma| a / 2.3 digits: an error delta in gamma
grows like delta exp(|gamma| x) across the guide) by exp(M t) of the 4 x 4 matrix obtained from the
six curl equations, e_x and h_x follow from the two algebraic ones, and the Poynting flux is integrated
by Clenshaw-Curtis quadrature on pieces of each layer short enough for the field's growth rate.  (The
null vector of the collocation matrix itself gives the fractions only to about 1e-5.)
"""
import json
import math
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as sla

C = 299792458.0
A = 10 * math.log(10)               # 200 dB/m in Np/m


def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = c[-1] = 2
    c = c * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    D = np.outer(c, 1 / c) / (X - X.T + np.eye(N + 1))
    return D - np.diag(D.sum(axis=1)), x


def eps_of(m, f, lib=np):
    if m['type'] == 'dielectric':
        e, t = (mp.mpf(str(m['eps'])), mp.mpf(str(m['loss_tangent']))) if lib is mp else (float(m['eps']), float(m['loss_tangent']))
        return e * (1 - 1j * t)
    if lib is mp:
        fp, nu, f = mp.mpf(str(m['f_p'])), mp.mpf(str(m['collision_frequency'])), mp.mpf(f)
        return 1 - fp ** 2 / (f * (f - 1j * nu / (2 * mp.pi)))
    fp, nu = float(m['f_p']), float(m['collision_frequency'])
    return 1 - fp ** 2 / (f * (f - 1j * nu / (2 * math.pi)))


def pencil(kind, rec, f, scale):
    lay = rec['layers']
    eps = [eps_of(L['material'], f) for L in lay]
    k0sq = (2 * math.pi * f / C) ** 2
    bl = []
    for L in lay:
        d = float(L['thickness'])
        N = int(min(140, max(20, scale * d * math.sqrt(k0sq * max(abs(e) for e in eps)) + 20)))
        D, _ = cheb(N)
        D1 = D * 2 / d
        bl.append((N, D1, D1 @ D1))
    off = np.concatenate([[0], np.cumsum([N + 1 for N, _, _ in bl])])
    n = off[-1]
    Am = np.zeros((n, n), complex)
    Bm = np.zeros((n, n), complex)
    for i, (N, D1, D2) in enumerate(bl):
        o = off[i]
        for k in range(1, N):
            Am[o + k, o:o + N + 1] = D2[k]
            Am[o + k, o + k] += k0sq * eps[i]
            Bm[o + k, o + k] = 1
    last = len(bl) - 1
    N0, D10, _ = bl[0]
    NL, D1L, _ = bl[last]
    if kind == 'LSE':
        Am[off[0] + N0, off[0] + N0] = 1
        Am[off[last], off[last]] = 1
    else:
        Am[off[0] + N0, off[0]:off[1]] = D10[N0]
        Am[off[last], off[last]:off[last + 1]] = D1L[0]
    for i in range(last):
        Ni, D1i, _ = bl[i]
        Nj, D1j, _ = bl[i + 1]
        r1, r2 = off[i], off[i + 1] + Nj
        Am[r1, off[i]] = 1
        Am[r1, off[i + 1] + Nj] = -1
        wi = 1.0 if kind == 'LSE' else 1 / eps[i]
        wj = 1.0 if kind == 'LSE' else 1 / eps[i + 1]
        Am[r2, off[i]:off[i + 1]] = wi * D1i[0]
        Am[r2, off[i + 1]:off[i + 2]] -= wj * D1j[Nj]
    return Am, Bm


def det(kind, rec, q, f):
    """D_kind by 2 x 2 matrix exponentials (mpmath)."""
    k2 = (2 * mp.pi / C) ** 2 * mp.mpf(f) ** 2
    st = mp.matrix([[0], [1]]) if kind == 'LSE' else mp.matrix([[1], [0]])
    ep = None
    for L in rec['layers']:
        e = eps_of(L['material'], f, lib=mp)
        if kind == 'LSM' and ep is not None:
            st[1] = st[1] * e / ep
        M = mp.matrix([[0, 1], [q - k2 * e, 0]])
        st = mp.expm(M * mp.mpf(str(L['thickness']))) * st
        ep = e
    return st[0] if kind == 'LSE' else st[1] / ep


def eigen_q(kind, rec, f, b):
    k0sq = (2 * math.pi * f / C) ** 2
    out = []
    w1 = sla.eig(*pencil(kind, rec, f, 3.0), right=False)
    w2 = sla.eig(*pencil(kind, rec, f, 4.5), right=False)
    w1 = w1[np.isfinite(w1)]
    w2 = w2[np.isfinite(w2)]
    for q in w2:
        # only zeros that can have alpha <= A in some family: Re sqrt(-q) <= 1.05 A
        g = np.sqrt(complex(-q))
        if abs(g.real) > 1.05 * A or abs(q) > 5e4 * k0sq:
            continue
        if np.min(np.abs(w1 - q)) > 1e-4 * max(abs(q), k0sq):
            continue
        out.append(complex(q))
    return out


def polish(kind, rec, f, q):
    return _root(lambda x: det(kind, rec, x, f), q, abs(np.sqrt(complex(q))) * float(rec['a']))


# ------------------------------------------------------------------ ferrite guides
def polder(m, f, lib=np):
    if m['type'] != 'ferrite':
        return 1.0, 0.0
    if lib is mp:
        g = mp.mpf('2.8e6')
        H0, dH, Ms = (mp.mpf(str(m[k])) for k in ('bias_field_oe', 'linewidth_oe', 'saturation_magnetization_gauss'))
        f = mp.mpf(f)
    else:
        g = 2.8e6
        H0, dH, Ms = (float(m[k]) for k in ('bias_field_oe', 'linewidth_oe', 'saturation_magnetization_gauss'))
    s = 1 if m['bias_direction'] == '+y' else -1
    f0 = g * H0 + 1j * g * dH / 2
    fm = g * Ms
    den = f0 * f0 - f * f
    return 1 + f0 * fm / den, s * f * fm / den


def eps_any(m, f, lib=np):
    if m['type'] == 'ferrite':
        mm = dict(m, type='dielectric')
        return eps_of(mm, f, lib)
    return eps_of(m, f, lib)


def te_pencil(rec, f, scale, bt=0.0):
    """(A0 + g A1 + g^2 A2) x = 0 for E'' + (k0^2 eps mu_e + g^2) E = 0, E continuous,
    (E' - j tau g E)/mu_e continuous, E = 0 at the walls; linearised as a 2n x 2n pencil."""
    lay = rec['layers']
    k0sq = (2 * math.pi * f / C) ** 2
    prm = []
    for L in lay:
        mu, ka = polder(L['material'], f)
        mue = (mu * mu - ka * ka) / mu
        prm.append((eps_any(L['material'], f) * mue, mue, ka / mu))
    bl = []
    for L, (em, _, _) in zip(lay, prm):
        d = float(L['thickness'])
        N = int(min(260, max(20, scale * d * math.sqrt(k0sq * abs(em)) + 20, 3 * scale * math.sqrt(bt * d) + 20)))
        D, _ = cheb(N)
        D1 = D * 2 / d
        bl.append((N, D1, D1 @ D1))
    off = np.concatenate([[0], np.cumsum([N + 1 for N, _, _ in bl])])
    n = off[-1]
    A0 = np.zeros((n, n), complex)
    A1 = np.zeros((n, n), complex)
    A2 = np.zeros((n, n), complex)
    for i, (N, D1, D2) in enumerate(bl):
        o = off[i]
        for k in range(1, N):
            A0[o + k, o:o + N + 1] = D2[k]
            A0[o + k, o + k] += k0sq * prm[i][0]
            A2[o + k, o + k] = 1
    last = len(bl) - 1
    A0[off[0] + bl[0][0], off[0] + bl[0][0]] = 1
    A0[off[last], off[last]] = 1
    for i in range(last):
        Ni, D1i, _ = bl[i]
        Nj, D1j, _ = bl[i + 1]
        r1, r2 = off[i], off[i + 1] + Nj
        A0[r1, off[i]] = 1
        A0[r1, off[i + 1] + Nj] = -1
        (_, mL, tL), (_, mR, tR) = prm[i], prm[i + 1]
        A0[r2, off[i]:off[i + 1]] = D1i[0] / mL
        A0[r2, off[i + 1]:off[i + 2]] -= D1j[Nj] / mR
        A1[r2, off[i]] += -1j * tL / mL
        A1[r2, off[i + 1] + Nj] -= -1j * tR / mR
    Z = np.zeros((n, n), complex)
    I = np.eye(n)
    return np.block([[Z, I], [A0, A1]]), np.block([[I, Z], [Z, -A2]])


def te_det(rec, g, f):
    """E(a) by 2 x 2 matrix exponentials of (E, E'), with the jump of E' at each interface."""
    k2 = (2 * mp.pi / C) ** 2 * mp.mpf(f) ** 2
    st = mp.matrix([[0], [1]])
    prev = None
    for L in rec['layers']:
        m = L['material']
        mu, ka = polder(m, f, lib=mp)
        mue = (mu * mu - ka * ka) / mu
        tau = ka / mu
        if prev is not None:
            mL, tL = prev
            P = (st[1] - 1j * tL * g * st[0]) / mL
            st[1] = mue * P + 1j * tau * g * st[0]
        M = mp.matrix([[0, 1], [-(k2 * eps_any(m, f, lib=mp) * mue + g * g), 0]])
        st = mp.expm(M * mp.mpf(str(L['thickness']))) * st
        prev = (mue, tau)
    return st[0]


SCALE = 2.5       # collocation resolution factor (n >= 1)
EXACT = {}      # high-precision zeros (gamma or q), keyed by their double-precision values


def _root(fun, g, a=0.0):
    """Secant iteration (mpmath, 50 digits) from g.  |D| grows like exp(|gamma| a), so convergence is
    judged by the step, not by |D|: the result is accepted if a second run from it moves it by less
    than 1e-20 |r|.  Returns None if the iteration fails."""
    dps = 50 + int(a / 2.3)
    with mp.workdps(dps):
        try:
            x0 = mp.mpc(g.real, g.imag)
            r = mp.findroot(fun, x0, tol=mp.mpf(10) ** (-dps - 10), verify=False, maxsteps=80)
            r2 = mp.findroot(fun, r, tol=mp.mpf(10) ** (-dps - 10), verify=False, maxsteps=20)
        except (ValueError, ZeroDivisionError):
            return None
        if abs(r2 - r) > mp.mpf(10) ** (20 - dps) * abs(r):
            return None
        EXACT[complex(r2)] = r2
        return complex(r2)


def qs_target(rec):
    """Quasi-static |beta| scale 3 / d of the thinnest ferrite or plasma layer (0 without one): slow
    magnetostatic and electrostatic surface waves have |beta| d of order one."""
    ds = [float(L['thickness']) for L in rec['layers'] if L['material']['type'] != 'dielectric']
    return 3.0 / min(ds) if ds else 0.0


def te_solve(rec, f):
    k0 = 2 * math.pi * f / C
    bt = qs_target(rec)
    while True:
        zs = _te_zeros(rec, f, bt)
        top = max([abs(z.imag) for z in zs] + [0.0])
        if top <= bt / 2:
            break
        bt = 2 * max(bt, top)
    fwd = [[z.imag, z.real, power_split(rec, f, 0, z)] for z in zs if 0 < z.real <= A]
    rev = [[-z.imag, -z.real, power_split(rec, f, 0, z)] for z in zs if 0 < -z.real <= A]
    return fwd, rev


def _te_zeros(rec, f, bt):
    k0 = 2 * math.pi * f / C
    w2 = sla.eig(*te_pencil(rec, f, 4.5, bt), right=False)
    w2 = w2[np.isfinite(w2)]
    zs = []
    for g in w2:
        if abs(g.real) > 1.05 * A or abs(g) > 50 * k0 + 1.5 * bt:
            continue
        r = _root(lambda x: te_det(rec, x, f), g, abs(g) * float(rec['a']))
        if r is None or abs(r - g) > 1e-3 * max(abs(g), 1):
            continue
        if all(abs(r - z) > 1e-8 * max(abs(r), 1) for z in zs):
            zs.append(r)
    return zs


# ------------------------------------------------------------------ six-field collocation
def tensor6(m, f, lib=np):
    """(e1, e3, e2, mu, ka): [eps] = [[e1, 0, -j e2], [0, e3, 0], [j e2, 0, e1]],
    [mu] = [[mu, 0, -j ka], [0, 1, 0], [j ka, 0, mu]] of a material at frequency f."""
    J = mp.mpc(0, 1) if lib is mp else 1j
    num = (lambda v: mp.mpf(str(v))) if lib is mp else float
    f = num(f)
    t = m['type']
    if t == 'magnetized_plasma':
        fp, nu, fc = num(m['f_p']), num(m['collision_frequency']), num(m['cyclotron_frequency'])
        s = 1 if m['bias_direction'] == '+y' else -1
        pi = mp.pi if lib is mp else math.pi
        w = f - J * nu / (2 * pi)                       # f - j nu / (2 pi)
        e1 = 1 - fp ** 2 * w / (f * (w ** 2 - fc ** 2))
        e3 = 1 - fp ** 2 / (f * w)
        gx = fp ** 2 * fc / (f * (w ** 2 - fc ** 2))
        # eps_xz = +j s gx  ->  -j e2 = j s gx
        return e1, e3, -s * gx, 1, 0
    e = eps_any(m, f, lib)
    mu, ka = polder(m, f, lib)
    return e, e, 0, mu, ka


def cc_weights(N):
    """Clenshaw-Curtis weights on x_k = cos(pi k / N), k = 0..N, for the interval [-1, 1]."""
    th = np.pi * np.arange(N + 1) / N
    w = np.zeros(N + 1)
    v = np.ones(N - 1)
    if N % 2 == 0:
        w[0] = w[N] = 1 / (N * N - 1)
        for k in range(1, N // 2):
            v -= 2 * np.cos(2 * k * th[1:N]) / (4 * k * k - 1)
        v -= np.cos(N * th[1:N]) / (N * N - 1)
    else:
        w[0] = w[N] = 1 / (N * N)
        for k in range(1, (N - 1) // 2 + 1):
            v -= 2 * np.cos(2 * k * th[1:N]) / (4 * k * k - 1)
    w[1:N] = 2 * v / N
    return w


class Six:
    """Fields (e_x, e_y, e_z, h_x, h_y, h_z) of family n (normalised h = eta0 H) on Chebyshev nodes of
    every layer; the six curl equations are A0 x + gamma A1 x = 0.  Row bookkeeping: the e_y' and
    e_z' equations are dropped at the right end of every layer and at the left end of layer 0, the
    h_y' and h_z' equations at the left end of layers 1.., and those rows carry the wall conditions
    e_y = e_z = 0 and the continuity of e_y, e_z, h_y, h_z."""

    def __init__(self, rec, f, n, scale, bt=0.0):
        self.lay = rec['layers']
        k0 = 2 * math.pi * f / C
        ky = n * math.pi / float(rec['b'])
        self.ten = [tensor6(L['material'], f) for L in self.lay]
        self.bl = []
        for L, (e1, e3, e2, mu, ka) in zip(self.lay, self.ten):
            d = float(L['thickness'])
            ref = k0 * math.sqrt(max(abs(e1), abs(e3), abs(e1) + abs(e2)) * max(abs(mu), abs(mu) + abs(ka), 1.0))
            N = int(min(200, max(16, scale * d * ref + 16, 2 * scale * math.sqrt((bt + ky) * d) + 16)))
            D, _ = cheb(N)
            self.bl.append((N, D * 2 / d, d))
        self.off = np.concatenate([[0], np.cumsum([6 * (N + 1) for N, _, _ in self.bl])])
        nn = self.off[-1]
        A0 = np.zeros((nn, nn), complex)
        A1 = np.zeros((nn, nn), complex)
        jk = 1j * k0
        last = len(self.bl) - 1

        def ix(i, c, k):                      # unknown: layer i, field c (ex ey ez hx hy hz), node k
            return self.off[i] + c * (self.bl[i][0] + 1) + k

        for i, ((N, D1, d), (e1, e3, e2, mu, ka)) in enumerate(zip(self.bl, self.ten)):
            EX, EY, EZ, HX, HY, HZ = range(6)
            for k in range(N + 1):
                # (1) ky ez + g ey + jk (mu hx - j ka hz) = 0      [x component of curl E]
                r = ix(i, 0, k)
                A0[r, ix(i, EZ, k)] += ky
                A1[r, ix(i, EY, k)] += 1
                A0[r, ix(i, HX, k)] += jk * mu
                A0[r, ix(i, HZ, k)] += k0 * ka
                # (4) g hy - ky hz - jk (e1 ex - j e2 ez) = 0      [x component of curl h]
                r = ix(i, 3, k)
                A1[r, ix(i, HY, k)] += 1
                A0[r, ix(i, HZ, k)] -= ky
                A0[r, ix(i, EX, k)] -= jk * e1
                A0[r, ix(i, EZ, k)] -= k0 * e2
                drop_e = (k == 0) or (i == 0 and k == N)
                drop_h = (i > 0 and k == N)
                if not drop_e:
                    # (2) -g ex - ez' + jk hy = 0
                    r = ix(i, 2, k)
                    A1[r, ix(i, EX, k)] -= 1
                    A0[r, ix(i, EZ, 0):ix(i, EZ, 0) + N + 1] -= D1[k]
                    A0[r, ix(i, HY, k)] += jk
                    # (3) ey' - ky ex + jk (j ka hx + mu hz) = 0
                    r = ix(i, 1, k)
                    A0[r, ix(i, EY, 0):ix(i, EY, 0) + N + 1] += D1[k]
                    A0[r, ix(i, EX, k)] -= ky
                    A0[r, ix(i, HX, k)] -= k0 * ka
                    A0[r, ix(i, HZ, k)] += jk * mu
                if not drop_h:
                    # (5) -g hx - hz' - jk e3 ey = 0
                    r = ix(i, 5, k)
                    A1[r, ix(i, HX, k)] -= 1
                    A0[r, ix(i, HZ, 0):ix(i, HZ, 0) + N + 1] -= D1[k]
                    A0[r, ix(i, EY, k)] -= jk * e3
                    # (6) hy' + ky hx - jk (j e2 ex + e1 ez) = 0
                    r = ix(i, 4, k)
                    A0[r, ix(i, HY, 0):ix(i, HY, 0) + N + 1] += D1[k]
                    A0[r, ix(i, HX, k)] += ky
                    A0[r, ix(i, EX, k)] += k0 * e2
                    A0[r, ix(i, EZ, k)] -= jk * e1
            # boundary rows
            if i == 0:
                A0[ix(0, 1, N), ix(0, EY, N)] = 1          # e_y(0) = 0
                A0[ix(0, 2, N), ix(0, EZ, N)] = 1          # e_z(0) = 0
            if i == last:
                A0[ix(i, 1, 0), ix(i, EY, 0)] = 1          # e_y(a) = 0
                A0[ix(i, 2, 0), ix(i, EZ, 0)] = 1          # e_z(a) = 0
            else:
                Nj = self.bl[i + 1][0]
                for row, c in ((1, EY), (2, EZ)):
                    A0[ix(i, row, 0), ix(i, c, 0)] = 1
                    A0[ix(i, row, 0), ix(i + 1, c, Nj)] = -1
                for row, c in ((4, HY), (5, HZ)):
                    A0[ix(i + 1, row, Nj), ix(i, c, 0)] = 1
                    A0[ix(i + 1, row, Nj), ix(i + 1, c, Nj)] = -1
        self.A0, self.A1, self.ix = A0, A1, ix

    def gammas(self):
        w = sla.eig(self.A0, -self.A1, right=False)
        return w[np.isfinite(w)]

    def power(self, g):
        """Per-layer fractions of the time-averaged power flow Re int (e_x h_y* - e_y h_x*) dx."""
        Mt = self.A0 + g * self.A1
        lu = sla.lu_factor(Mt)
        rng = np.random.default_rng(7)
        x = rng.standard_normal(Mt.shape[0]) + 1j * rng.standard_normal(Mt.shape[0])
        for _ in range(3):
            x = sla.lu_solve(lu, x)
            x /= np.abs(x).max()
        P = []
        for i, (N, _, d) in enumerate(self.bl):
            fld = [x[self.ix(i, c, 0):self.ix(i, c, 0) + N + 1] for c in range(6)]
            ex, ey, _, hx, hy, _ = fld
            P.append(float(np.real(cc_weights(N) @ (ex * np.conj(hy) - ey * np.conj(hx)))) * d / 2)
        tot = sum(P)
        return [p / tot for p in P]


def state_matrix(L, g, f, n, b):
    """(M, e_x row, h_x row) in mpmath for the state (e_y, e_z, h_y, h_z) of family n in layer L."""
    k0 = 2 * mp.pi * mp.mpf(f) / C
    ky = n * mp.pi / mp.mpf(str(b))
    jk = mp.mpc(0, 1) * k0
    e1, e3, e2, mu, ka = tensor6(L['material'], f, lib=mp)
    # e_x and h_x from the two algebraic equations, as rows acting on (e_y, e_z, h_y, h_z)
    ex = [0, -k0 * e2 / (jk * e1), g / (jk * e1), -ky / (jk * e1)]
    hx = [-g / (jk * mu), -ky / (jk * mu), 0, -k0 * ka / (jk * mu)]
    M = mp.matrix(4, 4)
    for c in range(4):
        M[0, c] = ky * ex[c] + k0 * ka * hx[c]
        M[1, c] = -g * ex[c]
        M[2, c] = -ky * hx[c] - k0 * e2 * ex[c]
        M[3, c] = -g * hx[c]
    M[0, 3] -= jk * mu
    M[1, 2] += jk
    M[2, 1] += jk * e1
    M[3, 0] -= jk * e3
    return M, ex, hx


def hyb_det(rec, g, f, n):
    """det of the (e_y, e_z) rows of prod exp(M d) [0; I] for the state (e_y, e_z, h_y, h_z)."""
    S = mp.matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
    for L in rec['layers']:
        M, _, _ = state_matrix(L, g, f, n, rec['b'])
        S = mp.expm(M * mp.mpf(str(L['thickness']))) * S
    return S[0, 0] * S[1, 1] - S[0, 1] * S[1, 0]


def power_split(rec, f, n, g, nodes=24):
    """Per-layer fractions of Re int (e_x h_y* - e_y h_x*) dx for the mode gamma = g of family n."""
    with mp.workdps(60 + int(abs(g) * float(rec['a']) / 2.3)):
        gm = EXACT.get(complex(g), mp.mpc(g.real, g.imag))
        lay = [(L, state_matrix(L, gm, f, n, rec['b']), mp.mpf(str(L['thickness']))) for L in rec['layers']]
        if n == 0:
            s = mp.matrix([[0], [0], [0], [1]])
        else:
            S = mp.matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
            for _, (M, _, _), d in lay:
                S = mp.expm(M * d) * S
            a, b2, c, d2 = S[0, 0], S[0, 1], S[1, 0], S[1, 1]
            v = (-b2, a) if abs(a) + abs(b2) >= abs(c) + abs(d2) else (-d2, c)
            s = mp.matrix([[0], [0], [v[0]], [v[1]]])
        xk = [mp.cos(mp.pi * k / nodes) for k in range(nodes + 1)]
        wk = [mp.mpf(x) for x in cc_weights(nodes)]
        P = []
        for L, (M, ex, hx), d in lay:
            Mf = np.array([[complex(M[r, c]) for c in range(4)] for r in range(4)])
            rate = max(abs(np.linalg.eigvals(Mf)))
            pieces = max(1, int(math.ceil(float(rate * d) / 2.0)))
            h = d / pieces
            tot = mp.mpf(0)
            sp = s
            Enode = [mp.expm(M * (h * (x + 1) / 2)) for x in xk]     # the same on every piece
            Eh = mp.expm(M * h)
            for _ in range(pieces):
                for En, w in zip(Enode, wk):
                    p = En * sp
                    exv = sum(ex[k] * p[k] for k in range(4))
                    hxv = sum(hx[k] * p[k] for k in range(4))
                    tot += w * h / 2 * mp.re(exv * mp.conj(p[2]) - p[0] * mp.conj(hxv))
                sp = Eh * sp
            P.append(tot)
            s = mp.expm(M * d) * s
        T = sum(P)
        return [float(x / T) for x in P]


def hyb_family(rec, f, n, bt0=0.0):
    """Zeros of family n: the collocation is resolved for |beta| up to bt (at least the quasi-static target
    and bt0, twice the largest |beta| of the previous family), and bt is doubled while a zero lies beyond
    bt / 2."""
    bt = max(qs_target(rec), bt0)
    while True:
        zs = _hyb_zeros(rec, f, n, bt)
        top = max([abs(z.imag) for z in zs] + [0.0])
        if top <= bt / 2:
            return zs
        bt = 2 * max(bt, top)


def _hyb_zeros(rec, f, n, bt):
    k0 = 2 * math.pi * f / C
    ky = n * math.pi / float(rec['b'])
    w2 = Six(rec, f, n, SCALE, bt).gammas()
    zs = []
    for g in w2:
        if abs(g.real) > 1.05 * A or abs(g.real) < 1e-9 or abs(g) > 50 * k0 + 4 * ky + 1.5 * bt:
            continue
        r = _root(lambda x: hyb_det(rec, x, f, n), g, abs(g) * float(rec['a']))
        if r is None or abs(r - g) > 1e-3 * max(abs(g), 1):
            continue
        if all(abs(r - z) > 1e-8 * max(abs(r), 1) for z in zs):
            zs.append(r)
    return zs


def gyro_solve(rec, f):
    k0 = 2 * math.pi * f / C
    b = float(rec['b'])
    fwd, rev = te_solve(rec, f)
    out_f = {'LSE0': fwd} if fwd else {}
    out_r = {'LSE0': rev} if rev else {}
    m = 1.0
    for L in rec['layers']:
        e1, e3, e2, mu, ka = tensor6(L['material'], f)
        m = max(m, max(abs(e1), abs(e3), abs(e1) + abs(e2)) * max(abs(mu), abs(mu) + abs(ka), 1.0))
    n, empty, bt_prev = 1, 0, 0.0
    while True:
        ky = n * math.pi / b
        zs = hyb_family(rec, f, n, bt_prev)
        bt_prev = 2 * max([abs(z.imag) for z in zs] + [bt_prev / 2])
        if zs:
            fw = sorted([[z.imag, z.real, power_split(rec, f, n, z)] for z in zs if 0 < z.real <= A])
            rv = sorted([[-z.imag, -z.real, power_split(rec, f, n, z)] for z in zs if 0 < -z.real <= A])
            if fw:
                out_f['HYB%d' % n] = fw
            if rv:
                out_r['HYB%d' % n] = rv
        empty = 0 if zs else empty + 1
        if empty >= 2 and ky * ky > k0 * k0 * m + A * A:
            break
        n += 1
    return {'forward': out_f, 'reverse': out_r}


def solve(rec):
    b = float(rec['b'])
    if any(L['material']['type'] in ('ferrite', 'magnetized_plasma') for L in rec['layers']):
        return {'modes': [gyro_solve(rec, float(f)) for f in rec['test_frequencies']]}
    res = []
    for f in rec['test_frequencies']:
        f = float(f)
        entry = {}
        for kind in ('LSE', 'LSM'):
            zs = []
            for q in eigen_q(kind, rec, f, b):
                p = polish(kind, rec, f, q)
                if p is not None and all(abs(p - z) > 1e-8 * max(abs(p), 1) for z in zs):
                    zs.append(p)
            n = 0 if kind == 'LSE' else 1
            qmax = max([z.real for z in zs] + [0.0])
            while (n * math.pi / b) ** 2 - A * A <= qmax:
                lst = []
                for q in zs:
                    g = np.sqrt(complex((n * math.pi / b) ** 2 - q))
                    if g.real < 0:
                        g = -g
                    if g.real <= A:
                        lst.append(g)
                        qx = EXACT.get(q)
                        if qx is not None:
                            with mp.workdps(60 + int(abs(g) * float(rec['a']) / 2.3)):
                                gx = mp.sqrt((n * mp.pi / mp.mpf(str(rec['b']))) ** 2 - qx)
                                EXACT[g] = -gx if mp.re(gx) < 0 else gx
                if lst:
                    entry['%s%d' % (kind, n)] = [[g.imag, g.real, power_split(rec, f, n, g)] for g in lst]
                n += 1
        res.append({'forward': entry, 'reverse': {k: [[p[0], p[1], list(p[2])] for p in v] for k, v in entry.items()}})
    return {'modes': res}


def main():
    guides = json.load(open(sys.argv[1]))['guides']
    only = sys.argv[3:]
    out = {}
    for rec in guides:
        if only and rec['id'] not in only:
            continue
        out[rec['id']] = solve(rec)
        print(rec['id'], 'done', flush=True)
    with open(sys.argv[2], 'w') as fh:
        json.dump({'guides': out}, fh)


if __name__ == '__main__':
    main()
