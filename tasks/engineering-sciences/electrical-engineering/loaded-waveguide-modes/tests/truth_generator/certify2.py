"""Certified enumeration for the lossy guides (ball arithmetic, python-flint).

modes(G, kind, f, A1):  every zero q of D_kind(., f) with Re sqrt(-q) <= A1 (the parabolic region
    Re q >= -A1^2 + (Im q)^2 / (4 A1^2)), as certified disjoint boxes, each holding exactly one
    simple zero (complex Krawczyk test).  Since Re sqrt(k_n^2 - q) >= Re sqrt(-q), this contains
    every zero that has attenuation <= A1 in any family.
crossings(G, kind, n, A):  every (beta, f), f in the band, with D_kind(k_n^2 - (A + j beta)^2, f) = 0,
    i.e. every frequency at which a mode of family n has attenuation exactly A (2-D real Krawczyk).
Upper bounds on Re q: LSE by the Rayleigh quotient (Re q <= k0^2 max Re eps); LSM by tail_free(),
    a proof that D_LSM has no zero with Re q >= Q in the region (impedance-ratio recursion).
"""
import math

from flint import acb, arb

from wgacb import A as ARB, TWO_PI_C, Du, Guide


def box(r0, r1, i0, i1):
    return acb(arb((r0 + r1) / 2, (r1 - r0) / 2), arb((i0 + i1) / 2, (i1 - i0) / 2))


def lohi(x):
    return float(x.lower()), float(x.upper())


# ------------------------------------------------------------------ LSM tail
def _t_final(G, eps, T, ratio):
    """Impedance-ratio recursion for D_LSM: t = v'/(kappa v), t = 0 at x = 0; D = 0 iff t(a) = 0.
    Returns False if 0 may lie in t(a) (or a denominator may vanish)."""
    t = acb(0)
    for i in range(len(G.layers)):
        if i:
            t = t * (eps[i] / eps[i - 1]) * ratio[i - 1]
        den = 1 + t * T[i]
        if den.contains(0):
            return False
        t = (t + T[i]) / den
    return not t.contains(0)


def tail_free(G, f, Q, A1, nimag=8, ratio_step=1.05, span=1e3):
    """True if D_LSM has no zero with Re q >= Q and Re sqrt(-q) <= A1 (so |Im q| <= 2 A1 sqrt(Re q + A1^2)),
    for every f in the ball f.  [Q, span Q] is covered by narrow boxes on which tanh(kappa d) is
    enclosed directly; beyond span Q, |tanh(kappa d) - 1| <= 2e/(1 - e) with e = exp(-2 d Re kappa)
    and Re kappa >= sqrt(Re(q - k0^2 eps))."""
    f = acb(f)
    k2 = TWO_PI_C * TWO_PI_C
    ke = [k2 * f * f * G.eps(L, f).v for L in G.layers]
    eps = [G.eps(L, f).v for L in G.layers]
    lo = arb(Q)
    top = arb(Q) * span
    while lo < top:
        hi = lo * ratio_step
        him = 2 * A1 * (hi + A1 * A1).sqrt()
        for k in range(nimag):
            i0 = -him + 2 * him * k / nimag
            i1 = -him + 2 * him * (k + 1) / nimag
            qb = acb(arb((lo + hi) / 2, (hi - lo) / 2), arb((i0 + i1) / 2, ((i1 - i0) / 2).upper()))
            kap = [(qb - x).sqrt() for x in ke]
            if not all(kk.real > 0 for kk in kap):
                return False
            T = [(kk * L[0]).tanh() for kk, L in zip(kap, G.layers)]
            # kappa_i / kappa_{i+1} = sqrt(1 + (ke_{i+1} - ke_i) / (q - ke_{i+1})), a tight enclosure
            ratio = [(1 + (ke[i + 1] - ke[i]) / (qb - ke[i + 1])).sqrt() for i in range(len(kap) - 1)]
            if not _t_final(G, eps, T, ratio):
                return False
        lo = hi
    # unbounded rest, Re q >= lo
    r = (1 / lo).upper()
    sig = acb(arb(0, r), arb(0, r))
    W = [(1 - x * sig).sqrt() for x in ke]
    T = []
    for x, L in zip(ke, G.layers):
        m = lo - x.real.upper()
        if not m > 0:
            return False
        e = (-2 * L[0] * m.sqrt()).exp()
        rr = 2 * e / (1 - e)
        T.append(acb(arb(1, rr.upper()), arb(0, rr.upper())))
    ratio = [W[i] / W[i + 1] for i in range(len(W) - 1)]
    return _t_final(G, eps, T, ratio)


def lsm_bound(G, f, A1):
    """Smallest Q = 2^k Q0 for which tail_free succeeds."""
    f = acb(f)
    k2 = TWO_PI_C * TWO_PI_C
    m = max(float(abs(k2 * f * f * G.eps(L, f).v).upper()) for L in G.layers)
    dmin = min(float(L[0]) for L in G.layers)
    Q = 2 * m
    for _ in range(40):
        if tail_free(G, f, Q, A1):
            return Q
        Q *= 1.5
    raise RuntimeError('no LSM tail bound for %s' % G.id)


def lse_bound(G, f):
    f = acb(f)
    k2 = TWO_PI_C * TWO_PI_C
    return max(float((k2 * f * f * G.eps(L, f).v).real.upper()) for L in G.layers) * 1.0001 + 1


# ------------------------------------------------------------------ modes at one frequency
def _outside(r0, r1, i0, i1, A1):
    """Box entirely outside Re q >= -A1^2 + (Im q)^2/(4 A1^2)?"""
    if r1 < -A1 * A1:
        return True
    m = 0.0 if i0 <= 0 <= i1 else min(abs(i0), abs(i1))
    return m * m / (4 * A1 * A1) - A1 * A1 > r1


def modes(G, kind, f, A1, Q, maxbox=2_000_000):
    """Certified zeros of D_kind(., f) in the region Re sqrt(-q) <= A1, Re q <= Q."""
    f = ARB(f) if not isinstance(f, arb) else f
    fa = acb(f)
    A1 = float(A1)
    r0, r1 = -A1 * A1 * 1.0001, Q
    im = 2 * A1 * math.sqrt(Q + A1 * A1) * 1.0001
    # initial strips, so the first boxes are not absurdly elongated
    stack = []
    nstrip = max(1, int((r1 - r0) / max(im, 1.0)))
    nstrip = min(nstrip, 4000)
    edges = [r0 + (r1 - r0) * k / nstrip for k in range(nstrip + 1)]
    for a, b in zip(edges, edges[1:]):
        h = 2 * A1 * math.sqrt(max(b, 0) + A1 * A1) * 1.0001
        stack.append((a, b, -h, h))
    found, nbox = [], 0
    while stack:
        r0, r1, i0, i1 = stack.pop()
        nbox += 1
        if nbox > maxbox:
            raise RuntimeError('too many boxes %s %s %s' % (G.id, kind, f))
        if _outside(r0, r1, i0, i1, A1):
            continue
        B = box(r0, r1, i0, i1)
        res = _test(G, kind, B, fa)
        if res == 'empty':
            continue
        if isinstance(res, acb):
            found.append(res)
            continue
        # bisect the longer side, measured in units of the local imaginary extent
        if (r1 - r0) >= (i1 - i0):
            m = (r0 + r1) / 2
            stack += [(r0, m, i0, i1), (m, r1, i0, i1)]
        else:
            m = (i0 + i1) / 2
            stack += [(r0, r1, i0, m), (r0, r1, m, i1)]
    out = []
    for z in found:
        z = refine(G, kind, z, fa)
        if not any(z.overlaps(w) for w in out):
            out.append(z)
    return out, nbox


def _test(G, kind, B, fa):
    c = acb(B.real.mid(), B.imag.mid())
    DB = G.D(kind, B, fa)
    if DB is None:
        return None
    Dc = G.D(kind, c, fa, mask=G.mask)
    if not DB.v.contains(0):
        return 'empty'
    mv = Dc.v + DB.q * (B - c)
    if not mv.contains(0):
        return 'empty'
    if DB.q.contains(0):
        return None
    Y = 1 / acb(Dc.q.real.mid(), Dc.q.imag.mid())
    K = c - Y * Dc.v + (1 - Y * DB.q) * (B - c)
    if B.contains_interior(K):
        return B
    return None


def refine(G, kind, B, fa, rounds=60):
    """Newton from the box midpoint, then a Krawczyk check on a tight box around the result."""
    z = acb(B.real.mid(), B.imag.mid())
    for _ in range(rounds):
        d = G.D(kind, z, fa)
        step = d.v / d.q
        z = acb((z - step).real.mid(), (z - step).imag.mid())
        if float(step.abs_upper()) < 1e-70 * max(1.0, float(z.abs_upper())):
            break
    r = arb(1e-55) * max(1.0, float(z.abs_upper()))
    T = acb(arb(z.real.mid(), r), arb(z.imag.mid(), r))
    DT = G.D(kind, T, fa)
    Dz = G.D(kind, z, fa, mask=G.mask)
    Y = 1 / acb(Dz.q.real.mid(), Dz.q.imag.mid())
    K = z - Y * Dz.v + (1 - Y * DT.q) * (T - z)
    assert T.contains_interior(K), 'refinement not certified'
    assert B.contains(K), 'refinement left the certified box'
    return K


# ------------------------------------------------------------------ attenuation crossings
def crossings(G, kind, n, A, pieces, maxbox=3_000_000):
    """Certified (beta, f) boxes with D(k_n^2 - (A + j beta)^2, f) = 0, f in the band,
    searched piece by piece with the q bound of each piece (piece_bounds)."""
    out, nbox = [], 0
    for a, b, Q in pieces:
        X, nb = _crossings_piece(G, kind, n, A, Q, float(a.mid()), float(b.mid()), maxbox)
        nbox += nb
        for x in X:
            if not any(x[0].overlaps(y[0]) and x[1].overlaps(y[1]) for y in out):
                out.append(x)
    return out, nbox


def _crossings_piece(G, kind, n, A, Q, f0, f1, maxbox):
    kn = arb.pi() * n / G.b
    kn2 = kn * kn
    Aa = ARB(A) if not isinstance(A, arb) else A
    bmax2 = Q - float(kn2.mid()) + float(Aa.mid()) ** 2
    if bmax2 <= 0:
        return [], 0
    bmax = math.sqrt(bmax2) * 1.0001
    stack = []
    nb = 64
    for k in range(nb):
        stack.append((-bmax + 2 * bmax * k / nb, -bmax + 2 * bmax * (k + 1) / nb, f0, f1))
    found, nbox = [], 0
    while stack:
        b0, b1, g0, g1 = stack.pop()
        nbox += 1
        if nbox > maxbox:
            raise RuntimeError('too many crossing boxes %s %s%d' % (G.id, kind, n))
        res = _test2(G, kind, kn2, Aa, arb((b0 + b1) / 2, (b1 - b0) / 2), arb((g0 + g1) / 2, (g1 - g0) / 2))
        if res == 'empty':
            continue
        if res is not None:
            found.append(res)
            continue
        if (b1 - b0) / bmax >= (g1 - g0) / (f1 - f0):
            m = (b0 + b1) / 2
            stack += [(b0, m, g0, g1), (m, b1, g0, g1)]
        else:
            m = (g0 + g1) / 2
            stack += [(b0, b1, g0, m), (b0, b1, m, g1)]
    out = []
    for X in found:
        X = refine2(G, kind, kn2, Aa, X)
        if not any(X[0].overlaps(Y[0]) and X[1].overlaps(Y[1]) for Y in out):
            out.append(X)
    return out, nbox


def _F(G, kind, kn2, Aa, beta, f, mask=None):
    """D at q = kn2 - (A + j beta)^2 as (value, d/dbeta, d/df)."""
    g = acb(Aa, beta)
    q = kn2 - g * g
    d = G.D(kind, q, acb(f), mask=mask)
    if d is None:
        return None
    _F.mask = G.mask
    dqdb = -2 * g * acb(0, 1)
    return d.v, d.q * dqdb, d.f


def _test2(G, kind, kn2, Aa, bB, fB):
    bc, fc = arb(bB.mid()), arb(fB.mid())
    FB = _F(G, kind, kn2, Aa, bB, fB)
    if FB is None:
        return None
    Fc = _F(G, kind, kn2, Aa, bc, fc, mask=G.mask)
    if not FB[0].contains(0):
        return 'empty'
    mv = Fc[0] + FB[1] * (bB - bc) + FB[2] * (fB - fc)
    if not mv.contains(0):
        return 'empty'
    # 2-D real Krawczyk on (Re D, Im D)
    J = [[FB[1].real, FB[2].real], [FB[1].imag, FB[2].imag]]
    Jc = [[float(Fc[1].real.mid()), float(Fc[2].real.mid())], [float(Fc[1].imag.mid()), float(Fc[2].imag.mid())]]
    det = Jc[0][0] * Jc[1][1] - Jc[0][1] * Jc[1][0]
    if det == 0:
        return None
    Y = [[arb(Jc[1][1] / det), arb(-Jc[0][1] / det)], [arb(-Jc[1][0] / det), arb(Jc[0][0] / det)]]
    F0 = [Fc[0].real, Fc[0].imag]
    dx = [bB - bc, fB - fc]
    K = []
    for i in range(2):
        s = [bc, fc][i] - (Y[i][0] * F0[0] + Y[i][1] * F0[1])
        for j in range(2):
            m = (1 if i == j else 0) - (Y[i][0] * J[0][j] + Y[i][1] * J[1][j])
            s += m * dx[j]
        K.append(s)
    if bB.contains_interior(K[0]) and fB.contains_interior(K[1]):
        return (bB, fB)
    return None


def refine2(G, kind, kn2, Aa, X, rounds=60):
    b, f = arb(X[0].mid()), arb(X[1].mid())
    for _ in range(rounds):
        v, db, df = _F(G, kind, kn2, Aa, b, f)
        J = [[db.real, df.real], [db.imag, df.imag]]
        det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
        sb = (J[1][1] * v.real - J[0][1] * v.imag) / det
        sf = (-J[1][0] * v.real + J[0][0] * v.imag) / det
        b, f = arb((b - sb).mid()), arb((f - sf).mid())
        if float(abs(sf).upper()) < 1e-60 * float(abs(f).upper()):
            break
    rb = arb(1e-50) * max(1.0, abs(float(b.mid())))
    rf = arb(1e-50) * float(f.mid())
    res = _test2(G, kind, kn2, Aa, arb(b.mid(), rb.upper()), arb(f.mid(), rf.upper()))
    assert isinstance(res, tuple), 'crossing refinement not certified'
    assert X[0].contains(res[0]) and X[1].contains(res[1]), 'refinement left the certified box'
    return res


def piece_bounds(G, kind, A1, pieces=32, depth=6):
    """[(f0, f1, Q)] covering the band: Re q <= Q for every zero with Re sqrt(-q) <= A1 and f in [f0, f1].
    A piece on which the LSM tail proof fails is halved (up to `depth` times)."""
    f0, f1 = G.band
    todo = [(f0 + (f1 - f0) * k / pieces, f0 + (f1 - f0) * (k + 1) / pieces, 0) for k in range(pieces)]
    out = []
    while todo:
        a, b, dep = todo.pop()
        fb = arb((a + b) / 2, ((b - a) / 2).upper())
        if kind == 'LSE':
            out.append((a, b, lse_bound(G, fb)))
            continue
        try:
            out.append((a, b, lsm_bound(G, fb, A1)))
        except RuntimeError:
            if dep >= depth:
                raise
            m = (a + b) / 2
            todo += [(a, m, dep + 1), (m, b, dep + 1)]
    return sorted(out, key=lambda x: float(x[0].mid()))


# ------------------------------------------------------------------ ferrite guides: TE modes in the gamma plane
from wgacb import DTE, polder  # noqa: E402


def _te_params(G, f):
    f = acb(f)
    k2 = TWO_PI_C * TWO_PI_C * f * f
    out = []
    for L in G.layers:
        e = G._eps_point(L, f)[0]
        mu, ka = polder(L, f)
        mue = (mu * mu - ka * ka) / mu
        out.append((k2 * e * mue, mue, ka / mu))
    return out


def _sigma_final(G, prm, gb, T, kap, ratio):
    """sigma = kappa E / E' through the layers (sigma = 0 at x = 0); D = 0 iff sigma(a) = 0.
    Interface L -> R: rho_R = a rho_L + b, a = mu_eR kappa_L / (mu_eL kappa_R),
    b = j (tau_R - tau_L mu_eR / mu_eL) gamma / kappa_R, sigma = 1/rho."""
    s = acb(0)
    for i in range(len(G.layers)):
        if i:
            (_, meL, tL), (_, meR, tR) = prm[i - 1], prm[i]
            a = meR / meL * ratio[i - 1]
            b = acb(0, 1) * (tR - tL * meR / meL) * gb[i]
            den = a + b * s
            if den.contains(0):
                return False
            s = s / den
        den = 1 + s * T[i]
        if den.contains(0):
            return False
        s = (s + T[i]) / den
    return not s.contains(0)


def tail_free_te(G, f, B, A1, nalpha=8, ratio_step=1.03, span=1e3):
    """True if D_TE has no zero with |Im gamma| >= B and |Re gamma| <= A1 at frequency f."""
    prm = _te_params(G, f)
    for sign in (1, -1):
        lo = arb(B)
        top = arb(B) * span
        while lo < top:
            hi = lo * ratio_step
            for k in range(nalpha):
                a0 = -A1 + 2 * A1 * k / nalpha
                a1 = -A1 + 2 * A1 * (k + 1) / nalpha
                g = acb(arb((a0 + a1) / 2, (a1 - a0) / 2), sign * arb((lo + hi) / 2, ((hi - lo) / 2).upper()))
                # kappa = sign (gamma / j) W with W = sqrt(1 + ke / gamma^2) near 1 (tight enclosures)
                ig2 = 1 / (g * g)
                W = [(1 + ke * ig2).sqrt() for ke, _, _ in prm]
                gj = g * acb(0, -1) * sign
                kap = [gj * w for w in W]
                if not all(kk.real > 0 for kk in kap):
                    return False
                T = [(kk * L[0]).tanh() for kk, L in zip(kap, G.layers)]
                ratio = [W[i] / W[i + 1] for i in range(len(W) - 1)]
                gb = [acb(0, sign) / w for w in W]
                if not _sigma_final(G, prm, gb, T, kap, ratio):
                    return False
            lo = hi
        # unbounded rest: |beta| >= lo; kappa = (sign gamma / j) W, W = sqrt(1 + ke / gamma^2)
        r = (1 / lo).upper()
        sg = acb(arb(0, r), arb(0, r))                 # 1 / gamma
        W = [(1 + ke * sg * sg).sqrt() for ke, _, _ in prm]
        T = []
        for (ke, _, _), L in zip(prm, G.layers):
            m = lo * lo - arb(A1) * arb(A1) - abs(ke).upper()
            if not m > 0:
                return False
            e = (-2 * L[0] * m.sqrt()).exp()
            rr = 2 * e / (1 - e)
            T.append(acb(arb(1, rr.upper()), arb(0, rr.upper())))
        ratio = [W[i] / W[i + 1] for i in range(len(W) - 1)]
        gb = [acb(0, sign) / w for w in W]              # gamma / kappa = sign j / W
        if not _sigma_final(G, prm, gb, T, None, ratio):
            return False
    return True


def te_bound(G, f, A1):
    prm = _te_params(G, f)
    m = max(float(abs(ke).upper()) for ke, _, _ in prm)
    B = 2 * math.sqrt(m)
    for _ in range(40):
        if tail_free_te(G, f, B, A1):
            return B
        B *= 1.3
    raise RuntimeError('no TE tail bound for %s' % G.id)


def _test_te(G, B, fa):
    c = acb(B.real.mid(), B.imag.mid())
    DB = DTE(G, B, fa)
    if DB is None:
        return None
    Dc = DTE(G, c, fa, mask=G.mask)
    if not DB.v.contains(0):
        return 'empty'
    if not (Dc.v + DB.q * (B - c)).contains(0):
        return 'empty'
    if DB.q.contains(0):
        return None
    Y = 1 / acb(Dc.q.real.mid(), Dc.q.imag.mid())
    K = c - Y * Dc.v + (1 - Y * DB.q) * (B - c)
    if B.contains_interior(K):
        return B
    return None


def refine_te(G, B, fa, rounds=80):
    z = acb(B.real.mid(), B.imag.mid())
    for _ in range(rounds):
        d = DTE(G, z, fa)
        step = d.v / d.q
        z = acb((z - step).real.mid(), (z - step).imag.mid())
        if float(step.abs_upper()) < 1e-70 * max(1.0, float(z.abs_upper())):
            break
    r = arb(1e-55) * max(1.0, float(z.abs_upper()))
    T = acb(arb(z.real.mid(), r), arb(z.imag.mid(), r))
    DT = DTE(G, T, fa)
    Dz = DTE(G, z, fa, mask=G.mask)
    Y = 1 / acb(Dz.q.real.mid(), Dz.q.imag.mid())
    K = z - Y * Dz.v + (1 - Y * DT.q) * (T - z)
    assert T.contains_interior(K), 'TE refinement not certified'
    assert B.contains(K), 'TE refinement left the certified box'
    return K


def modes_te(G, f, A1, Bmax, maxbox=3_000_000):
    """Certified zeros gamma of D_TE(., f) with |Re gamma| <= A1, |Im gamma| <= Bmax."""
    fa = acb(ARB(f) if not isinstance(f, arb) else f)
    A1 = float(A1)
    stack = []
    nb = max(4, int(2 * Bmax / max(4 * A1, 1.0)))
    nb = min(nb, 4000)
    for k in range(nb):
        b0 = -Bmax + 2 * Bmax * k / nb
        b1 = -Bmax + 2 * Bmax * (k + 1) / nb
        stack += [(-A1, 0.0, b0, b1), (0.0, A1, b0, b1)]
    found, nbox = [], 0
    while stack:
        r0, r1, i0, i1 = stack.pop()
        nbox += 1
        if nbox > maxbox:
            raise RuntimeError('too many TE boxes %s %s' % (G.id, f))
        Bx = box(r0, r1, i0, i1)
        res = _test_te(G, Bx, fa)
        if res == 'empty':
            continue
        if isinstance(res, acb):
            found.append(res)
            continue
        if (r1 - r0) >= (i1 - i0):
            m = (r0 + r1) / 2
            stack += [(r0, m, i0, i1), (m, r1, i0, i1)]
        else:
            m = (i0 + i1) / 2
            stack += [(r0, r1, i0, m), (r0, r1, m, i1)]
    out = []
    for z in found:
        z = refine_te(G, z, fa)
        if not any(z.overlaps(w) for w in out):
            out.append(z)
    return out, nbox
