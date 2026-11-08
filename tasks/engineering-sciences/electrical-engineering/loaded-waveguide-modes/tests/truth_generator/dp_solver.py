"""Third solver: everything in double precision (numpy / scipy), to show that the tolerances need
accurate numerics, not extended precision.

usage: python3 dp_solver.py guides.json out.json [guide_id ...]

The candidate zeros come from the Chebyshev collocation of alt_solver.py (its search, grids and
family loop are reused unchanged).  The two stages that alt_solver.py carries out in mpmath are
replaced by double-precision ones:
  * polishing: secant iteration on a determinant that is built from both walls.  The tangential
    state of the wall solutions is carried from each wall towards a matching interface, in steps
    with |kappa| h <= 2, and re-orthonormalised (QR) after every step, so nothing overflows and the
    span of each set stays accurate in the direction in which it grows.  The determinant of the two
    sets at the matching interface vanishes at a mode.  The matching interface is the one at which
    the two sets are closest to dependent at the collocation estimate (for a surface wave, a face
    of the film it is bound to);
  * power split: the null vector of the two sets at the matching interface fixes the mode; its
    coefficients are carried back towards each wall through the stored QR factors (a division by
    the growth, which is stable), and the Poynting flux Re int (e_x h_y* - e_y h_x*) dx is integrated
    by Clenshaw-Curtis quadrature on every step.  The field is never carried across the film in the
    direction in which it decays, which is what makes a one-sided shot need extra digits.
"""
import cmath
import json
import math
import os
import sys

import numpy as np
import scipy.linalg as sla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alt_solver as S  # noqa: E402

C = S.C
A = S.A
NODES = 12          # Clenshaw-Curtis nodes per step
STEP = 2.0          # |kappa| h per step


# ------------------------------------------------------------------ layer matrices
def state4(L, g, f, n, b):
    """(M, ex row, hx row) of the state (e_y, e_z, h_y, h_z) of family n (as alt_solver.state_matrix)."""
    k0 = 2 * math.pi * f / C
    ky = n * math.pi / b
    jk = 1j * k0
    e1, e3, e2, mu, ka = S.tensor6(L['material'], f)
    ex = np.array([0, -k0 * e2 / (jk * e1), g / (jk * e1), -ky / (jk * e1)], complex)
    hx = np.array([-g / (jk * mu), -ky / (jk * mu), 0, -k0 * ka / (jk * mu)], complex)
    M = np.zeros((4, 4), complex)
    M[0] = ky * ex + k0 * ka * hx
    M[1] = -g * ex
    M[2] = -ky * hx - k0 * e2 * ex
    M[3] = -g * hx
    M[0, 3] -= jk * mu
    M[1, 2] += jk
    M[2, 1] += jk * e1
    M[3, 0] -= jk * e3
    return M, ex, hx


def system(rec, f, n, g):
    """Per layer (M, ex, hx, d) restricted to the state of family n: for n = 0 the TE subsystem
    (e_y, h_z) (e_z = h_y = 0 for a field uniform in y between the top and bottom walls)."""
    b = float(rec['b'])
    idx = [0, 3] if n == 0 else [0, 1, 2, 3]
    out = []
    for L in rec['layers']:
        M, ex, hx = state4(L, g, f, n, b)
        out.append((M[np.ix_(idx, idx)], ex[idx], hx[idx], float(L['thickness'])))
    wall = np.array([[0.0], [1.0]], complex) if n == 0 else np.array([[0, 0], [0, 0], [1, 0], [0, 1]], complex)
    return out, wall, idx


def q_system(kind, rec, f, q):
    """Isotropic guides, family-independent: state (psi, psi' / w), w = 1 (LSE) or eps (LSM)."""
    k2 = (2 * math.pi * f / C) ** 2
    out = []
    for L in rec['layers']:
        e = S.eps_of(L['material'], f)
        w = 1.0 if kind == 'LSE' else e
        M = np.array([[0, w], [(q - k2 * e) / w, 0]], complex)
        out.append((M, None, None, float(L['thickness'])))
    wall = np.array([[0.0], [1.0]], complex) if kind == 'LSE' else np.array([[1.0], [0.0]], complex)
    return out, wall


def steps(lay):
    """[(layer index, M, h, E = expm(M h))] in the order of x, each step with |kappa| h <= STEP."""
    out = []
    for i, (M, _, _, d) in enumerate(lay):
        rate = max(abs(np.linalg.eigvals(M)).max(), 1e-12)
        k = max(1, int(math.ceil(rate * d / STEP)))
        h = d / k
        E = sla.expm(M * h)
        out += [(i, M, h, E)] * k
    return out


def carry(Q, mats):
    """Carry the basis Q through the step matrices; returns the bases at every step start, the R
    factors, the final orthonormal basis and the sum of log det R (for an analytic normalisation)."""
    bases, Rs, slog = [], [], 0j
    for E in mats:
        bases.append(Q)
        Q, R = np.linalg.qr(E @ Q)
        Rs.append(R)
        slog += cmath.log(np.linalg.det(R))
    return bases, Rs, Q, slog


def two_sided(lay, wall, m):
    """Wall bases carried to interface m from the left (forward) and from the right (backward)."""
    st = steps(lay)
    cut = sum(1 for s in st if s[0] < m)
    left, right = st[:cut], st[cut:]
    bl, Rl, QL, sl = carry(wall, [E for _, _, _, E in left])
    br, Rr, QR, sr = carry(wall, [sla.expm(-M * h) for _, M, h, _ in reversed(right)])
    return (left, bl, Rl, QL, sl), (right, br, Rr, QR, sr)


def match_value(lay, wall, m, ref=0j):
    (_, _, _, QL, sl), (_, _, _, QR, sr) = two_sided(lay, wall, m)
    return np.linalg.det(np.hstack([QL, QR])) * cmath.exp(sl + sr - ref), sl + sr


def best_interface(lay, wall):
    """The interface at which the two wall sets are closest to dependent."""
    best = None
    for m in range(len(lay) + 1):
        (_, _, _, QL, _), (_, _, _, QR, _) = two_sided(lay, wall, m)
        sv = np.linalg.svd(np.hstack([QL, QR]), compute_uv=False)
        r = sv[-1] / sv[0]
        if best is None or r < best[0]:
            best = (r, m)
    return best[1]


def secant(fun, x0):
    """Secant iteration in double precision; None if it does not settle."""
    x1 = x0 * (1 + 1e-7) + 1e-9
    f0, f1 = fun(x0), fun(x1)
    for _ in range(60):
        if f1 == f0:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if not np.isfinite(x2):
            return None
        x0, f0, x1 = x1, f1, x2
        f1 = fun(x1)
        if abs(x1 - x0) <= 4e-16 * abs(x1):
            return x1
    return x1 if abs(x1 - x0) <= 1e-12 * abs(x1) else None


def polish_gamma(rec, f, n, g):
    lay, wall, _ = system(rec, f, n, g)
    m = best_interface(lay, wall)
    _, ref = match_value(lay, wall, m)

    def fun(x):
        return match_value(system(rec, f, n, x)[0], wall, m, ref)[0]
    return secant(fun, complex(g))


def polish_q(kind, rec, f, q):
    lay, wall = q_system(kind, rec, f, q)
    m = best_interface(lay, wall)
    _, ref = match_value(lay, wall, m)

    def fun(x):
        return match_value(q_system(kind, rec, f, x)[0], wall, m, ref)[0]
    return secant(fun, complex(q))


# ------------------------------------------------------------------ power split
def cc(N):
    return np.cos(np.pi * np.arange(N + 1) / N), S.cc_weights(N)


def power_split(rec, f, n, g, nodes=NODES):
    """Per-layer fractions of Re int (e_x h_y* - e_y h_x*) dx for the mode gamma = g of family n."""
    g = complex(g)
    lay, wall, idx = system(rec, f, n, g)
    m = best_interface(lay, wall)
    (left, bl, Rl, QL, _), (right, br, Rr, QR, _) = two_sided(lay, wall, m)
    k = QL.shape[1]
    _, _, vh = np.linalg.svd(np.hstack([QL, -QR]))
    v = np.conj(vh[-1])
    cL, cR = v[:k], v[k:]
    xk, wk = cc(nodes)
    P = np.zeros(len(lay))
    pos = [idx.index(c) if c in idx else None for c in range(4)]

    def flux(i, s):
        """e_x h_y* - e_y h_x* of the state s (restricted) in layer i."""
        _, ex, hx, _ = lay[i]
        full = [s[p] if p is not None else 0 for p in pos]
        exv, hxv = ex @ s, hx @ s
        return (exv * np.conj(full[2]) - full[0] * np.conj(hxv)).real

    # left part: coefficients at every step start by back-substitution through R
    c = cL
    for j in range(len(left) - 1, -1, -1):
        c = sla.solve_triangular(Rl[j], c)
        i, M, h, _ = left[j]
        s0 = bl[j] @ c
        for x, w in zip(xk, wk):
            t = h * (x + 1) / 2
            P[i] += w * h / 2 * flux(i, sla.expm(M * t) @ s0)
    # right part: carried backward from x = a; step j covers the j-th step from the right
    c = cR
    rr = list(reversed(right))
    for j in range(len(rr) - 1, -1, -1):
        c = sla.solve_triangular(Rr[j], c)
        i, M, h, _ = rr[j]
        s0 = br[j] @ c                       # state at the right end of this step
        for x, w in zip(xk, wk):
            t = h * (x + 1) / 2
            P[i] += w * h / 2 * flux(i, sla.expm(-M * t) @ s0)
    T = P.sum()
    return [float(p / T) for p in P]


# ------------------------------------------------------------------ the alt_solver search, re-polished
def _te_zeros(rec, f, bt):
    k0 = 2 * math.pi * f / C
    w2 = sla.eig(*S.te_pencil(rec, f, 4.5, bt), right=False)
    w2 = w2[np.isfinite(w2)]
    zs = []
    for g in w2:
        if abs(g.real) > 1.05 * A or abs(g) > 50 * k0 + 1.5 * bt:
            continue
        r = polish_gamma(rec, f, 0, g)
        if r is None or abs(r - g) > 1e-3 * max(abs(g), 1):
            continue
        if all(abs(r - z) > 1e-8 * max(abs(r), 1) for z in zs):
            zs.append(r)
    return zs


def _hyb_zeros(rec, f, n, bt):
    k0 = 2 * math.pi * f / C
    ky = n * math.pi / float(rec['b'])
    w2 = S.Six(rec, f, n, S.SCALE, bt).gammas()
    zs = []
    for g in w2:
        if abs(g.real) > 1.05 * A or abs(g.real) < 1e-9 or abs(g) > 50 * k0 + 4 * ky + 1.5 * bt:
            continue
        r = polish_gamma(rec, f, n, g)
        if r is None or abs(r - g) > 1e-3 * max(abs(g), 1):
            continue
        if all(abs(r - z) > 1e-8 * max(abs(r), 1) for z in zs):
            zs.append(r)
    return zs


def polish(kind, rec, f, q):
    return polish_q(kind, rec, f, q)


S._te_zeros = _te_zeros
S._hyb_zeros = _hyb_zeros
S.polish = polish
S.power_split = power_split


def main():
    guides = json.load(open(sys.argv[1]))['guides']
    only = sys.argv[3:]
    out = {}
    for rec in guides:
        if only and rec['id'] not in only:
            continue
        out[rec['id']] = S.solve(rec)
        print(rec['id'], 'done', flush=True)
    with open(sys.argv[2], 'w') as fh:
        json.dump({'guides': out}, fh)


if __name__ == '__main__':
    main()
