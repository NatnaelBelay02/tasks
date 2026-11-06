"""Certified per-layer power split of a mode (ball arithmetic, acb.integral).

For a mode gamma of family n the tangential state psi = (e_y, e_z, h_y, h_z) (normalised h = eta0 H,
hybacb) starts at x = 0 from (0, 0, v) with v spanning the null space of the e-block of T(a) [0; I]
(n >= 1), or from (0, 0, 0, 1) for n = 0 (the TE problem: h_y stays 0).  Inside layer i,
psi(x_i + t) = exp(M t) psi(x_i), and the time-averaged power through the part of the cross-section
occupied by the layer is proportional to
    P_i = Re int_0^{d_i} (e_x conj(h_y) - e_y conj(h_x)) dt
(the factor 1/2, eta0 and the y integral, b/2 for n >= 1 and b for n = 0, are common to all layers).
The integrand is evaluated as e_x(t) c(h_y)(t) - e_y(t) c(h_x)(t) with c(g)(t) = conj(g(conj t)),
which equals the conjugate on the real path and is entire in t, so acb.integral gives a rigorous
enclosure.  Returned: the fractions P_i / sum_k P_k as arb balls, and sum |P_k| / |sum P_k|.
"""
from flint import acb, acb_mat, arb

from wgacb import TWO_PI_C, layer_funcs
import hybacb as H


def _state_exp(L, ten, g, k0, ky, t):
    """exp(M t) for a (possibly complex-ball) local coordinate t, and the e_x, h_x rows."""
    M, _, ex, hx = H.m_rows(g, k0, ky, ten)
    Mm = acb_mat(M)
    if H.gyrotropic(L):
        return (Mm * t).exp(), ex, hx
    z = k0 * k0 * ten[0] + g * g - ky * ky
    r = layer_funcs(z, t)
    if r is None:
        return None, ex, hx
    Cc, S, _, _ = r
    I4 = acb_mat([[1 if a == b else 0 for b in range(4)] for a in range(4)])
    return I4 * Cc + Mm * S, ex, hx


def _row(rowc, p):
    return sum((rowc[k] * p[k, 0] for k in range(4)), acb(0))


def fractions(G, n, gamma, f, rel_tol=arb(2) ** -90):
    fa = acb(f)
    k0 = TWO_PI_C * fa
    ky = arb.pi() * n / G.b
    ten = [H.tensors(L, fa) for L in G.layers]
    if n == 0:
        s = acb_mat([[0], [0], [0], [1]])
    else:
        S = acb_mat([[0, 0], [0, 0], [1, 0], [0, 1]])
        for L, t in zip(G.layers, ten):
            E, _, _ = _state_exp(L, t, gamma, k0, ky, acb(L[0]))
            S = E * S
        a, b, c, d = S[0, 0], S[0, 1], S[1, 0], S[1, 1]
        if (a.abs_lower() + b.abs_lower()) >= (c.abs_lower() + d.abs_lower()):
            v0, v1 = -b, a
        else:
            v0, v1 = -d, c
        s = acb_mat([[0], [0], [v0], [v1]])
    P = []
    for L, t in zip(G.layers, ten):
        d = acb(L[0])
        sc = s

        def integrand(x, analytic, L=L, t=t, sc=sc):
            E, ex, hx = _state_exp(L, t, gamma, k0, ky, x)
            Ec, exc, hxc = _state_exp(L, t, gamma, k0, ky, x.conjugate())
            if E is None or Ec is None:
                return acb('nan')
            p = E * sc
            pc = Ec * sc
            exv = _row(ex, p)
            hxv_c = _row(hxc, pc).conjugate()
            hy_c = pc[2, 0].conjugate()
            return exv * hy_c - p[0, 0] * hxv_c

        Ii = acb.integral(integrand, 0, d, rel_tol=rel_tol)
        P.append(Ii.real)
        E, _, _ = _state_exp(L, t, gamma, k0, ky, d)
        s = E * s
    tot = sum(P, arb(0))
    assert not tot.contains(0), (G.id, n, 'total power not certainly nonzero')
    return [p / tot for p in P], sum((abs(p) for p in P), arb(0)) / abs(tot)
