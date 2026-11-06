"""Second, separately written solver, used to score a correct submission that does not come from
the oracle or from the certified pipeline.

usage: python3 alt_solver.py guides.json out.json [guide_id ...]

It never evaluates the transfer-matrix dispersion function on a contour and uses no interval
arithmetic:
  * multi-domain Chebyshev collocation of v'' + k0^2 eps(x) v = q v (complex eps), with the wall
    and interface conditions as algebraic rows, gives a generalised eigenproblem A x = q B x;
    the discretisation is solved twice (N and 1.5 N points per layer) and an eigenvalue is
    kept if the two agree to 1e-4 (relative); the polishing below gives the accuracy;
  * each kept eigenvalue is polished by complex Newton iteration on a determinant built from
    2 x 2 matrix exponentials in mpmath (30 digits);
  * gamma = sqrt(k_n^2 - q) with Re gamma > 0 for every family n; modes with attenuation
    <= 200 dB/m are reported.
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
    with mp.workdps(30):
        try:
            r = mp.findroot(lambda x: det(kind, rec, x, f), mp.mpc(q.real, q.imag), tol=mp.mpf(10) ** -24)
        except (ValueError, ZeroDivisionError):
            return None
        return complex(r)


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


def te_pencil(rec, f, scale):
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
        N = int(min(140, max(20, scale * d * math.sqrt(k0sq * abs(em)) + 20)))
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


def te_solve(rec, f):
    k0 = 2 * math.pi * f / C
    w1 = sla.eig(*te_pencil(rec, f, 3.0), right=False)
    w2 = sla.eig(*te_pencil(rec, f, 4.5), right=False)
    w1 = w1[np.isfinite(w1)]
    w2 = w2[np.isfinite(w2)]
    zs = []
    for g in w2:
        if abs(g.real) > 1.05 * A or abs(g) > 50 * k0:
            continue
        if np.min(np.abs(w1 - g)) > 1e-4 * max(abs(g), k0):
            continue
        with mp.workdps(30):
            try:
                r = complex(mp.findroot(lambda x: te_det(rec, x, f), mp.mpc(g.real, g.imag), tol=mp.mpf(10) ** -24))
            except (ValueError, ZeroDivisionError):
                continue
        if all(abs(r - z) > 1e-8 * max(abs(r), 1) for z in zs):
            zs.append(r)
    fwd = [[z.imag, z.real] for z in zs if 0 < z.real <= A]
    rev = [[-z.imag, -z.real] for z in zs if 0 < -z.real <= A]
    return {'forward': {'LSE0': fwd} if fwd else {}, 'reverse': {'LSE0': rev} if rev else {}}


def solve(rec):
    b = float(rec['b'])
    if any(L['material']['type'] == 'ferrite' for L in rec['layers']):
        return {'modes': [te_solve(rec, float(f)) for f in rec['test_frequencies']]}
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
                        lst.append([g.imag, g.real])
                if lst:
                    entry['%s%d' % (kind, n)] = lst
                n += 1
        res.append({'forward': entry, 'reverse': {k: [list(p) for p in v] for k, v in entry.items()}})
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
