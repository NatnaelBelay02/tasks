"""Ball (acb) evaluation and certified enumeration of the modes of family n >= 1 of a guide with
gyrotropic layers (y-biased ferrite, y-magnetised plasma), and of the modal power split.

Fields of family n:  E_x, E_z, H_y ~ sin(k_y y);  E_y, H_x, H_z ~ cos(k_y y);  k_y = n pi / b;
all ~ exp(-gamma z).  With normalised magnetic field h = eta0 H (so that j omega mu0 and j omega eps0
both become j k0), and the layer tensors
    [mu]  = [[mu, 0, -j ka], [0, 1, 0], [j ka, 0, mu]],
    [eps] = [[e1, 0, -j e2], [0, e3, 0], [j e2, 0, e1]],
the x components of the curl equations give
    e_x = (gamma h_y - k_y h_z + j k0 e2 e_z) / (j k0 e1),
    h_x = (-(k_y e_z + gamma e_y) / (j k0) + j ka h_z) / mu,
and the tangential state psi = (e_y, e_z, h_y, h_z) obeys psi' = M(gamma) psi with
    e_y' = k_y e_x - j k0 (j ka h_x + mu h_z),     e_z' = -gamma e_x + j k0 h_y,
    h_y' = -k_y h_x + j k0 (j e2 e_x + e1 e_z),    h_z' = -gamma h_x - j k0 e3 e_y.
psi is continuous at interfaces and e_y = e_z = 0 on both side walls, so the modes of family n are
the zeros of D(gamma) = det(e-rows of T(a) [0; I]), T(a) the product of the layer transfer
matrices exp(M d).  D is entire in gamma (M is polynomial in gamma).
  * isotropic layer: M^2 = -K^2 I, K^2 = k0^2 eps + gamma^2 - k_y^2, so exp(M d) = C I + S M with
    C = cos(K d), S = sin(K d)/K (wgacb.layer_funcs, entire in K^2).  If the layer is certainly
    evanescent on the ball (Re kappa d > 1, kappa^2 = -K^2) the factor is divided by cosh(kappa d):
    I + (tanh(kappa d)/kappa) M.  cosh(kappa d) is analytic and nonzero there, so the zeros, and
    their simplicity, are unchanged.
  * gyrotropic layer: exp of the 8x8 block [[M d, M' d], [0, M d]] (acb_mat.exp, rigorous) gives
    exp(M d) and its gamma-derivative in one evaluation.
  * balancing: M is far from normal for |gamma| >> k0 (entries ~ gamma^2 / k0 next to entries ~ k0),
    and a norm-based enclosure of exp(M d) then overestimates by up to exp(|M| d - Re kappa d).  The
    state is therefore carried in balanced coordinates psi = T psi~, T = diag(t) with exact powers of
    two chosen (Parlett-Reinsch, at the ball centre) for the first gyrotropic layer:
    every layer matrix becomes T^-1 exp(M d) T = exp(T^-1 M T d), and
    D = t_1 t_2 / (t_3 t_4) det(e-rows of prod T^-1 E T [e_3 e_4]) is the same function.
"""
import scipy.linalg as _sl
import numpy as _np
import math

from flint import acb, acb_mat, arb

from wgacb import A as ARB, TWO_PI_C, layer_funcs
import certify2 as CZ

GYRO = ARB('2.8e6')
J = acb(0, 1)


def tensors(L, f):
    """(e1, e3, e2, mu, ka) of a wgacb layer tuple at the exact frequency f (acb)."""
    kind = L[1]
    if kind == 'd':
        e = acb(L[2], -L[2] * L[3])
        return e, e, acb(0), acb(1), acb(0)
    if kind == 'f':
        _, _, e, td, Ms, H0, dH, sgn = L
        ee = acb(e, -e * td)
        f0 = acb(GYRO * H0, GYRO * dH / 2)
        fm = GYRO * Ms
        den = f0 * f0 - f * f
        return ee, ee, acb(0), 1 + f0 * fm / den, sgn * f * fm / den
    if kind == 'p':
        _, _, fp, nu = L
        ee = 1 - fp * fp / (f * (f - acb(0, nu / (2 * arb.pi()))))
        return ee, ee, acb(0), acb(1), acb(0)
    # magnetised plasma: (d, 'm', f_p, nu, f_c, s);  [eps] = [[e1, 0, j s g], [0, e3, 0], [-j s g, 0, e1]]
    _, _, fp, nu, fc, sgn = L
    fn = f - acb(0, nu / (2 * arb.pi()))
    den = f * (fn * fn - fc * fc)
    e1 = 1 - fp * fp * fn / den
    e3 = 1 - fp * fp / (f * fn)
    g = fp * fp * fc / den
    return e1, e3, -sgn * g, acb(1), acb(0)


def gyrotropic(L):
    return L[1] in ('f', 'm')


def m_rows(g, k0, ky, ten):
    """M(gamma) and dM/dgamma as nested lists, plus the e_x and h_x coefficient rows."""
    e1, e3, e2, mu, ka = ten
    jk = J * k0
    ex = [acb(0), J * e2 / e1, g / (jk * e1), -ky / (jk * e1)]
    hx = [-g / (jk * mu), -ky / (jk * mu), acb(0), J * ka / mu]
    dex = [acb(0), acb(0), 1 / (jk * e1), acb(0)]
    dhx = [-1 / (jk * mu), acb(0), acb(0), acb(0)]
    M = [[acb(0)] * 4 for _ in range(4)]
    dM = [[acb(0)] * 4 for _ in range(4)]
    for c in range(4):
        M[0][c] = ky * ex[c] - jk * J * ka * hx[c]
        M[1][c] = -g * ex[c]
        M[2][c] = -ky * hx[c] + jk * J * e2 * ex[c]
        M[3][c] = -g * hx[c]
        dM[0][c] = ky * dex[c] - jk * J * ka * dhx[c]
        dM[1][c] = -ex[c] - g * dex[c]
        dM[2][c] = -ky * dhx[c] + jk * J * e2 * dex[c]
        dM[3][c] = -hx[c] - g * dhx[c]
    M[0][3] -= jk * mu
    M[1][2] += jk
    M[2][1] += jk * e1
    M[3][0] -= jk * e3
    return M, dM, ex, hx


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PIDX = {p: k for k, p in enumerate(PAIRS)}


def additive_compound(M):
    """M^[2]: (u ^ v)' = M u ^ v + u ^ M v, in the basis PAIRS; (M^[2] y)_ij = sum_k M_ik y_kj + M_jk y_ik."""
    A = [[acb(0)] * 6 for _ in range(6)]
    for p, (i, j) in enumerate(PAIRS):
        for k in range(4):
            if k != j:                                  # M_ik y_kj
                q, sg = (PIDX[(k, j)], 1) if k < j else (PIDX[(j, k)], -1)
                A[p][q] += sg * M[i][k]
            if k != i:                                  # M_jk y_ik
                q, sg = (PIDX[(i, k)], 1) if i < k else (PIDX[(k, i)], -1)
                A[p][q] += sg * M[j][k]
    return A


def wedge2(E):
    """Second compound (2 x 2 minors) of a 4 x 4 matrix: (E u ^ E v)_p = sum_q W_pq (u ^ v)_q."""
    return [[E[i][k] * E[j][l] - E[i][l] * E[j][k] for (k, l) in PAIRS] for (i, j) in PAIRS]


def dwedge2(E, dE):
    return [[dE[i][k] * E[j][l] + E[i][k] * dE[j][l] - dE[i][l] * E[j][k] - E[i][l] * dE[j][k]
             for (k, l) in PAIRS] for (i, j) in PAIRS]


class HybGuide:
    def __init__(self, G, n):
        self.G = G
        self.n = n
        self.ky = arb.pi() * n / G.b

    def params(self, f):
        return [tensors(L, f) for L in self.G.layers]

    @staticmethod
    def _gyro_exp(M, dM, d, V=None):
        """exp(M d) and its gamma-derivative from the 8x8 block exponential.  With V (a point matrix,
        approximate eigenvectors of M at the ball centre) the exponential is taken of V^-1 M V,
        which is nearly diagonal: exp(M d) = V exp(V^-1 M V d) V^-1 holds exactly for any invertible
        V, and the ball enclosure no longer pays for the non-normality of M."""
        Mm, dMm = acb_mat(M), acb_mat(dM)
        if V is not None:
            Vb = acb_mat(V)
            Vi = Vb.inv()
            Mm, dMm = Vi * Mm * Vb, Vi * dMm * Vb
        X = acb_mat(8, 8)
        for r in range(4):
            for c in range(4):
                X[r, c] = Mm[r, c] * d
                X[r + 4, c + 4] = Mm[r, c] * d
                X[r, c + 4] = dMm[r, c] * d
        Y = X.exp()
        E = acb_mat(4, 4)
        dE = acb_mat(4, 4)
        for r in range(4):
            for c in range(4):
                E[r, c] = Y[r, c]
                dE[r, c] = Y[r, c + 4]
        if V is not None:
            E, dE = Vb * E * Vi, Vb * dE * Vi
        return E, dE

    @staticmethod
    def _eigvecs(M):
        Mf = _np.array([[complex(float(x.real.mid()), float(x.imag.mid())) for x in row] for row in M])
        if not _np.all(_np.isfinite(Mf)):
            return None
        w, V = _np.linalg.eig(Mf)
        if not _np.all(_np.isfinite(V)) or _np.linalg.cond(V) > 1e8:
            return None
        return [[acb(float(V[r, c].real), float(V[r, c].imag)) for c in range(4)] for r in range(4)]

    def balance(self, g, k0, ten):
        """Exact power-of-two scaling vector t for the first gyrotropic layer at the centre of g."""
        c = acb(g.real.mid(), g.imag.mid())
        for L, t in zip(self.G.layers, ten):
            if gyrotropic(L):
                M, _, _, _ = m_rows(c, k0, self.ky, t)
                Mf = _np.array([[complex(float(x.real.mid()), float(x.imag.mid())) for x in row] for row in M])
                _, (sc, _) = _sl.matrix_balance(Mf, permute=False, separate=True)
                return [arb(float(x)) for x in sc]
        return [arb(1)] * 4

    def layer(self, L, ten, g, k0, k02, d, scale, tv=None):
        """(E, E') = (exp(M d), d/dgamma exp(M d)) as acb_mat, possibly divided by cosh(kappa d);
        with tv given, in balanced coordinates (T^-1 E T, T^-1 E' T)."""
        M, dM, _, _ = m_rows(g, k0, self.ky, ten)
        if tv is not None:
            M = [[M[r][c] * tv[c] / tv[r] for c in range(4)] for r in range(4)]
            dM = [[dM[r][c] * tv[c] / tv[r] for c in range(4)] for r in range(4)]
        Mm, dMm = acb_mat(M), acb_mat(dM)
        if gyrotropic(L):
            c = acb(g.real.mid(), g.imag.mid())
            Mc, _, _, _ = m_rows(c, k0, self.ky, ten)
            if tv is not None:
                Mc = [[Mc[r][q] * tv[q] / tv[r] for q in range(4)] for r in range(4)]
            V = self._eigvecs(Mc)
            E, dE = self._gyro_exp(M, dM, d, V)
            if g.real.rad() > 0 or g.imag.rad() > 0:
                # mean-value form E(g) in E(c) + E'(ball) (g - c): the direct ball enclosure of exp(M d)
                # loses up to two digits per layer once |gamma| >> k0; the centre value is tight.
                Ec, _ = self._gyro_exp(Mc, Mc, d, V)
                E = Ec + dE * (g - c)
            return E, dE, False
        z = k02 * ten[0] + g * g - self.ky * self.ky     # K^2
        dz = 2 * g
        I4 = acb_mat([[1 if r == c else 0 for c in range(4)] for r in range(4)])
        if scale is None:
            scale = (z.real < 0) or (z.imag > 0) or (z.imag < 0)
            if scale:
                kv = (-z).sqrt()
                scale = bool(kv.real * d > 1)
        if scale:
            kv = (-z).sqrt()
            T = (kv * d).tanh()
            R = T / kv                                      # tanh(kappa d) / kappa
            dk = -g / kv                                    # kappa' = -gamma / kappa
            dR = ((1 - T * T) * d * kv - T) / (kv * kv) * dk
            return I4 + Mm * R, dMm * R + Mm * dR, True
        r = layer_funcs(z, d)
        if r is None:
            return None
        Cc, S, Cz, Sz = r
        return I4 * Cc + Mm * S, I4 * (Cz * dz) + Mm * (Sz * dz) + dMm * S, False

    # ---- second compound (exterior square): D = (psi_1 ^ psi_2)_{e_y e_z} without cancellation ----
    def layer2m(self, L, ten, g, k0, k02, d, tv, want):
        """(C, C', scaled): wedge^2 of the balanced layer matrix and its gamma-derivative (6x6 acb_mat).
        want (None = choose) fixes whether an isotropic layer is divided by cosh^2(kappa d)."""
        M, dM, _, _ = m_rows(g, k0, self.ky, ten)
        M = [[M[r][c] * tv[c] / tv[r] for c in range(4)] for r in range(4)]
        dM = [[dM[r][c] * tv[c] / tv[r] for c in range(4)] for r in range(4)]
        A2, dA2 = acb_mat(additive_compound(M)), acb_mat(additive_compound(dM))
        if gyrotropic(L):
            c = acb(g.real.mid(), g.imag.mid())
            Mc, _, _, _ = m_rows(c, k0, self.ky, ten)
            Mc = [[Mc[r][q] * tv[q] / tv[r] for q in range(4)] for r in range(4)]
            Ac = additive_compound(Mc)
            V = self._eigvecs6(Ac)
            C, dC = self._block_exp(A2, dA2, d, V)
            if g.real.rad() > 0 or g.imag.rad() > 0:
                Am = acb_mat(Ac)
                if V is not None:
                    Vb = acb_mat(V)
                    Cc = Vb * (Vb.inv() * Am * Vb * d).exp() * Vb.inv()
                else:
                    Cc = (Am * d).exp()
                C = Cc + dC * (g - c)
            return C, dC, False
        # isotropic: exp(M d) = Cs I + Sn M  ->  wedge^2 = Cs^2 I + Cs Sn M^[2] + Sn^2 wedge^2 M
        W2 = wedge2(M)
        dW2 = dwedge2(M, dM)
        z = k02 * ten[0] + g * g - self.ky * self.ky     # K^2
        dz = 2 * g
        kv = None
        if want is None:
            if (z.real < 0) or (z.imag > 0) or (z.imag < 0):
                kv = (-z).sqrt()
                if not bool(kv.real * d > 1):
                    kv = None
        elif want:
            kv = (-z).sqrt()            # tanh(kappa d) / kappa and cosh(kappa d) are even in kappa
        if kv is not None:
            # divided by cosh^2(kappa d) (analytic, nonzero): Cs = 1, Sn = tanh(kappa d) / kappa
            T = (kv * d).tanh()
            Sn = T / kv
            dk = -g / kv
            dSn = ((1 - T * T) * d * kv - T) / (kv * kv) * dk
            Cs, dCs = acb(1), acb(0)
        else:
            r = layer_funcs(z, d)
            if r is None:
                return None
            Cs, Sn, Cz, Sz = r
            dCs, dSn = Cz * dz, Sz * dz
        I6 = acb_mat([[1 if a == b else 0 for b in range(6)] for a in range(6)])
        W2m, dW2m = acb_mat(W2), acb_mat(dW2)
        C = I6 * (Cs * Cs) + A2 * (Cs * Sn) + W2m * (Sn * Sn)
        dC = (I6 * (2 * Cs * dCs) + A2 * (dCs * Sn + Cs * dSn) + dA2 * (Cs * Sn)
              + W2m * (2 * Sn * dSn) + dW2m * (Sn * Sn))
        return C, dC, kv is not None

    @staticmethod
    def _block_exp(A, dA, d, V):
        n = 6
        if V is not None:
            Vb = acb_mat(V)
            Vi = Vb.inv()
            A, dA = Vi * A * Vb, Vi * dA * Vb
        X = acb_mat(2 * n, 2 * n)
        for r in range(n):
            for c in range(n):
                X[r, c] = A[r, c] * d
                X[r + n, c + n] = A[r, c] * d
                X[r, c + n] = dA[r, c] * d
        Y = X.exp()
        E = acb_mat(n, n)
        dE = acb_mat(n, n)
        for r in range(n):
            for c in range(n):
                E[r, c] = Y[r, c]
                dE[r, c] = Y[r, c + n]
        if V is not None:
            E, dE = Vb * E * Vi, Vb * dE * Vi
        return E, dE

    @staticmethod
    def _eigvecs6(A):
        Af = _np.array([[complex(float(x.real.mid()), float(x.imag.mid())) for x in row] for row in A])
        if not _np.all(_np.isfinite(Af)):
            return None
        w, V = _np.linalg.eig(Af)
        if not _np.all(_np.isfinite(V)) or _np.linalg.cond(V) > 1e8:
            return None
        return [[acb(float(V[r, c].real), float(V[r, c].imag)) for c in range(6)] for r in range(6)]

    def D(self, g, f, mask=None):
        """(D, dD/dgamma) at the ball g.  D = t_1 t_2 y_{(e_y, e_z)}, y the exterior product of the two
        states started from h_y and h_z, carried through the layers by the second compounds of the
        balanced layer matrices (no cancellation between nearly parallel growing states).  Isotropic
        layers certainly evanescent on the ball are divided by cosh^2(kappa d); mask (from the box
        evaluation) makes the centre evaluation use the same choice, i.e. the same function."""
        f = acb(f)
        k0 = TWO_PI_C * f
        k02 = k0 * k0
        ten = self.params(f)
        tv = self.balance(g, k0, ten)
        y = acb_mat(6, 1)
        y[5, 0] = 1 / (tv[2] * tv[3])
        dy = acb_mat(6, 1)
        used = []
        for i, (L, t) in enumerate(zip(self.G.layers, ten)):
            want = None if mask is None else mask[i]
            r = self.layer2m(L, t, g, k0, k02, L[0], tv, want)
            if r is None:
                return None
            C, dC, sc = r
            used.append(sc)
            y, dy = C * y, dC * y + C * dy
        self.mask = used
        w = tv[0] * tv[1]
        return y[0, 0] * w, dy[0, 0] * w

# ------------------------------------------------------------- certified enumeration ---
def _test(H, B, fa):
    c = acb(B.real.mid(), B.imag.mid())
    r = H.D(B, fa)
    if r is None:
        return None
    vB, dB = r
    if not vB.contains(0):
        return 'empty'
    mask = H.mask
    vc, dc = H.D(c, fa, mask=mask)
    if not (vc + dB * (B - c)).contains(0):
        return 'empty'
    if dB.contains(0):
        return None
    Y = 1 / acb(dc.real.mid(), dc.imag.mid())
    K = c - Y * vc + (1 - Y * dB) * (B - c)
    if B.contains_interior(K):
        return B
    return None


def refine(H, B, fa, rounds=80):
    z = acb(B.real.mid(), B.imag.mid())
    for _ in range(rounds):
        v, dv = H.D(z, fa)
        step = v / dv
        z = acb((z - step).real.mid(), (z - step).imag.mid())
        if float(step.abs_upper()) < 1e-70 * max(1.0, float(z.abs_upper())):
            break
    rr = arb(1e-55) * max(1.0, float(z.abs_upper()))
    T = acb(arb(z.real.mid(), rr), arb(z.imag.mid(), rr))
    vT, dT = H.D(T, fa)
    vz, dz = H.D(z, fa, mask=H.mask)
    Y = 1 / acb(dz.real.mid(), dz.imag.mid())
    K = z - Y * vz + (1 - Y * dT) * (T - z)
    assert T.contains_interior(K), 'hybrid refinement not certified'
    assert B.contains(K), 'hybrid refinement left the certified box'
    return K


def modes(H, f, A1, Bmax, maxbox=4_000_000):
    """Certified zeros gamma of D_n(., f) with |Re gamma| <= A1, |Im gamma| <= Bmax (disjoint boxes)."""
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
    import flint
    prec0 = flint.ctx.prec
    flint.ctx.prec = 128          # enumeration; the refinement below runs at the full precision
    try:
        found, nbox = _enumerate(H, fa, stack, maxbox, f)
    finally:
        flint.ctx.prec = prec0
    out = []
    for z in found:
        z = refine(H, z, fa)
        if not any(z.overlaps(w) for w in out):
            out.append(z)
    return out, nbox


def _enumerate(H, fa, stack, maxbox, f):
    found, nbox = [], 0
    while stack:
        r0, r1, i0, i1 = stack.pop()
        nbox += 1
        if nbox > maxbox:
            raise RuntimeError('too many hybrid boxes %s %s n=%d' % (H.G.id, f, H.n))
        Bx = CZ.box(r0, r1, i0, i1)
        res = _test(H, Bx, fa)
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
    return found, nbox


def reach(G, f):
    """k0^2 max |eps| max(|mu|, |mu_e|, 1) over the layers (float), and k0."""
    fa = acb(ARB(f))
    k0 = float((TWO_PI_C * ARB(f)).mid())
    m = 1.0
    for L in G.layers:
        e1, e3, e2, mu, ka = tensors(L, fa)
        mue = (mu * mu - ka * ka) / mu
        ee = max(float(e1.abs_upper()), float(e3.abs_upper()), float((e1.abs_upper() + e2.abs_upper())))
        mm = max(float(mu.abs_upper()), float(mue.abs_upper()), float(mu.abs_upper() + ka.abs_upper()), 1.0)
        m = max(m, ee * mm)
    return k0, m
