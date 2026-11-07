"""Independent author-side solver for pmu-grid-identification (shares no code with solution/identify.py).

usage: python3 independent.py <input_dir> <output.json>

Model: per frame the 13 free bus voltages are unknown states; the measurement of channel k is
h_k = c_k z_k(theta, x_t), and the residual is the exact polar log ratio
    r_k = [ln|m_k / h_k| / sigma_m,  arg(m_k / h_k) / sigma_a],
which is the maximum-likelihood residual of the multiplicative noise model.  States and parameters are
estimated jointly by Levenberg-Marquardt; the per-frame state blocks are eliminated by a Schur complement.
  * start: equation-error linear least squares of measured currents on measured voltages;
  * superset: every non-known corridor in series-admittance coordinates (g, b, charging multiplier);
  * malfunctions: exact repeats (stale), then a CUSUM change-point search on each PMU's per-frame
    chi-square block, the strongest non-overlapping run first, refitted up to four rounds;
  * topology: Wald-ordered backward elimination, a corridor is dropped when the likelihood ratio
    statistic of keeping it is below 30;
  * final fit in physical coordinates on all unmasked data.
"""
import json
import sys

import numpy as np

IN, OUT = sys.argv[1], sys.argv[2]
NET = json.load(open(f"{IN}/network.json"))
CANDS = NET["candidates"]
NB = len(NET["buses"])
ZB = NET["zero_injection_bus"] - 1
REF = set(NET["reference_pmus"])
PMU = [b for b in NET["buses"] if b - 1 != ZB]
CH = [("V", b) for b in PMU] + [("I", b) for b in PMU]
NC = len(CH)
FREE = [i for i in range(NB) if i != ZB]
SM = np.array([NET["noise"][c]["magnitude_rel_std"] for c, _ in CH])
SA = np.array([NET["noise"][c]["angle_std_rad"] for c, _ in CH])
CALK = [k for k, (c, b) in enumerate(CH) if b not in REF]

data = np.genfromtxt(f"{IN}/pmu.csv", delimiter=",", skip_header=1)
cols = open(f"{IN}/pmu.csv").readline().strip().split(",")
T = len(data)
MEAS = np.stack([data[:, cols.index(f"{c}{b}_mag")] * np.exp(1j * data[:, cols.index(f"{c}{b}_ang")]) for c, b in CH], 1)
LOGM = np.log(MEAS)


# ------------------------------------------------------------------ parameters
class Par:
    """Parameter vector: branch block (per active corridor) then 2 per calibrated channel."""

    def __init__(self, active, phys):
        self.active, self.phys = list(active), phys
        self.idx, n = {}, 0
        for c in CANDS:
            if c["known"] or c["id"] not in self.active:
                continue
            if c["model"] in ("line", "series"):
                keys = ["r", "x"] if phys else (["g", "b", "k"] if c["model"] == "line" else ["g", "b"])
            elif c["model"] == "transformer":
                keys = ["x", "tap"]
            else:
                keys = ["x"]
            for k in keys:
                self.idx[(c["id"], k)] = n
                n += 1
        self.nbr = n
        for k in CALK:
            self.idx[("cal", k, "rho")] = n
            self.idx[("cal", k, "phi")] = n + 1
            n += 2
        self.n = n

    def nominal(self, open_zero=True):
        th = np.zeros(self.n)
        for c in CANDS:
            if c["known"] or c["id"] not in self.active:
                continue
            p = c["nominal"]
            if c["model"] in ("line", "series"):
                if self.phys:
                    th[self.idx[(c["id"], "r")]], th[self.idx[(c["id"], "x")]] = p["r"], p["x"]
                else:
                    y = 1 / complex(p["r"], p["x"])
                    on = 0.0 if (open_zero and c["database_status"] == "open") else 1.0
                    th[self.idx[(c["id"], "g")]], th[self.idx[(c["id"], "b")]] = on * y.real, on * y.imag
                    if c["model"] == "line":
                        th[self.idx[(c["id"], "k")]] = on
            elif c["model"] == "transformer":
                th[self.idx[(c["id"], "x")]], th[self.idx[(c["id"], "tap")]] = p["x"], p["tap"]
            else:
                th[self.idx[(c["id"], "x")]] = p["x"]
        return th

    def branch_terms(self, th):
        """List of (i, j, [(dYii, dYij, dYji, dYjj) per parameter index]) and Y itself."""
        Y = np.zeros((NB, NB), complex)
        terms = []
        for c in CANDS:
            if not (c["known"] or c["id"] in self.active):
                continue
            i, j, p, m = c["from_bus"] - 1, c["to_bus"] - 1, c["nominal"], c["model"]
            bc = p.get("b", 0.0) if m == "line" else 0.0
            d = []
            if c["known"]:
                ys, t, kap = 1 / complex(p.get("r", 0.0), p["x"]), p.get("tap", 1.0), 1.0
            elif m in ("line", "series") and self.phys:
                r, x = th[self.idx[(c["id"], "r")]], th[self.idx[(c["id"], "x")]]
                ys, t, kap = 1 / complex(r, x), 1.0, 1.0
                dz = -ys * ys
                d += [(self.idx[(c["id"], "r")], dz, 0, 0), (self.idx[(c["id"], "x")], 1j * dz, 0, 0)]
            elif m in ("line", "series"):
                ys = complex(th[self.idx[(c["id"], "g")]], th[self.idx[(c["id"], "b")]])
                t = 1.0
                kap = th[self.idx[(c["id"], "k")]] if m == "line" else 1.0
                d += [(self.idx[(c["id"], "g")], 1.0, 0, 0), (self.idx[(c["id"], "b")], 1j, 0, 0)]
                if m == "line":
                    d += [(self.idx[(c["id"], "k")], 0, 0.5j * bc, 0)]
            elif m == "transformer":
                x, t = th[self.idx[(c["id"], "x")]], th[self.idx[(c["id"], "tap")]]
                ys, kap = 1 / complex(0, x), 1.0
                d += [(self.idx[(c["id"], "x")], 1j / (x * x), 0, 0), (self.idx[(c["id"], "tap")], 0, 0, 1.0)]
            else:
                x = th[self.idx[(c["id"], "x")]]
                ys, t, kap = 1 / complex(0, x), 1.0, 1.0
                d += [(self.idx[(c["id"], "x")], 1j / (x * x), 0, 0)]
            sh = 0.5j * bc * kap
            Y[i, i] += (ys + sh) / t ** 2
            Y[j, j] += ys + sh
            Y[i, j] -= ys / t
            Y[j, i] -= ys / t
            for pi, dys, dsh, dt in d:
                # derivative of the 2x2 branch stamp
                a = (dys + dsh) / t ** 2 - 2 * (ys + sh) / t ** 3 * dt
                b = -dys / t + ys / t ** 2 * dt
                e = dys + dsh
                terms.append((pi, i, j, a, b, e))
        Y[NET["shunt"]["bus"] - 1, NET["shunt"]["bus"] - 1] += 1j * NET["shunt"]["nominal_b"]
        return Y, terms

    def cal(self, th):
        c = np.ones(NC, complex)
        for k in CALK:
            c[k] = (1 + th[self.idx[("cal", k, "rho")]]) * np.exp(1j * th[self.idx[("cal", k, "phi")]])
        return c


def elim(Y):
    S = np.zeros((NB, len(FREE)), complex)
    S[FREE, np.arange(len(FREE))] = 1
    S[ZB] = -Y[ZB, FREE] / Y[ZB, ZB]
    return S


def model_channels(par, th, X):
    """Predicted channel phasors (n, 26) for states X (n, 13) and the map matrices."""
    Y, terms = par.branch_terms(th)
    S = elim(Y)
    Vf = X @ S.T                                      # (n, 14)
    I = Vf @ Y.T
    cal = par.cal(th)
    Z = np.concatenate([Vf[:, [b - 1 for b in PMU]], I[:, [b - 1 for b in PMU]]], 1) * cal
    Amap = np.concatenate([S[[b - 1 for b in PMU]], (Y @ S)[[b - 1 for b in PMU]]], 0) * cal[:, None]   # (26, 13)
    return Z, Amap, Vf, terms, cal


def residual(par, th, X, W):
    Z, Amap, Vf, terms, cal = model_channels(par, th, X)
    if RECT:
        d = 1.0 - Z / MEAS[FR]                 # first-order form of log(m / h): linear in the states
    else:
        d = LOGM[FR] - np.log(Z)
        d = d.real + 1j * np.angle(np.exp(1j * d.imag))
    r = np.concatenate([d.real / SM, d.imag / SA], 1) * W
    return r, Z, Amap, Vf, terms, cal


def jac(par, th, Z, Amap, Vf, terms, cal, W, states_only=False):
    """Jacobians of the whitened residual: Jx (n, 52, 26), Jt (n, 52, P)."""
    n = len(Z)
    inv = 1.0 / (MEAS[FR] if RECT else Z)                 # d(-log h) = -dh / h ; d(1 - h/m) = -dh / m
    G = -Amap[None] * inv[:, :, None]                     # (n, 26, 13) complex, w.r.t. complex states
    Gx = np.concatenate([G, 1j * G], 2)                   # real and imaginary parts of the states
    Jx = np.concatenate([Gx.real / SM[None, :, None], Gx.imag / SA[None, :, None]], 1) * W[:, :, None]
    if states_only:
        return Jx, None
    Jt = np.zeros((n, NC, par.n), complex)
    rows = {b - 1: 13 + k for k, b in enumerate(PMU)}
    for pi, i, j, a, b, e in terms:
        # dI_i = a V_i + b V_j ; dI_j = b V_i + e V_j
        if i in rows:
            Jt[:, rows[i], pi] += -(a * Vf[:, i] + b * Vf[:, j]) * cal[rows[i]] * inv[:, rows[i]]
        if j in rows:
            Jt[:, rows[j], pi] += -(b * Vf[:, i] + e * Vf[:, j]) * cal[rows[j]] * inv[:, rows[j]]
    for k in CALK:
        # h ~ (1 + rho) e^{j phi}: d(-log h) = -1 / (1 + rho), -j ; d(1 - h/m) = -(h/m) / (1 + rho), -j h/m
        f = Z[:, k] * inv[:, k] if RECT else 1.0
        Jt[:, k, par.idx[("cal", k, "rho")]] = -f / abs(cal[k])
        Jt[:, k, par.idx[("cal", k, "phi")]] = -1j * f
    Jt = np.concatenate([Jt.real / SM[None, :, None], Jt.imag / SA[None, :, None]], 1) * W[:, :, None]
    return Jx, Jt


def cross_term(par, th, Z, Amap, Vf, terms, cal, W, r):
    """Second-order state x branch-parameter block sum_i r_i d2 r_i / dx dp (n, 26, P): the part of the
    Hessian that Gauss-Newton drops; with many per-frame nuisance states it is what makes plain
    Gauss-Newton converge only linearly (the same term separates the Golub-Pereyra and Kaufman
    variable-projection Jacobians)."""
    n = len(Z)
    Y, _ = par.branch_terms(th)
    S = elim(Y)
    q = (r[:, :NC] / SM + 1j * r[:, NC:] / SA) * W[:, :NC]          # complex weighted residual (n, 26)
    C = np.zeros((n, 26, par.n))
    rows = {b - 1: 13 + k for k, b in enumerate(PMU)}
    for pi, i, j, a, b, e in terms:
        for bus, (ci, cj) in ((i, (a, b)), (j, (b, e))):
            if bus not in rows:
                continue
            k = rows[bus]
            dA = (ci * S[i] + cj * S[j]) * cal[k]                     # d A[k, :] / dp  (13,)
            N = (ci * Vf[:, i] + cj * Vf[:, j]) * cal[k]              # d h_k / dp      (n,)
            if RECT:
                d2 = dA[None, :] / MEAS[FR][:, k:k + 1] + 0 * N[:, None]
            else:
                d2 = dA[None, :] / Z[:, k:k + 1] - Amap[k][None, :] * (N / Z[:, k] ** 2)[:, None]   # d2 l / dx dp
            w = -np.conj(q[:, k])[:, None]
            C[:, :13, pi] += np.real(w * d2)
            C[:, 13:, pi] += np.real(w * 1j * d2)
    return C


FR = np.arange(T)
VERB = False
RECT = True
NEWTON = True


def lm(par, th, X, W, iters=30, tol=1e-9, fixed=False):
    """Joint Levenberg-Marquardt over parameters and per-frame states (Schur complement); with fixed=True
    only the states move (Gauss-Newton per frame)."""
    lam = 1e-3
    rel = 0.0
    r, *mc = residual(par, th, X, W)
    cost = float(np.sum(r * r))
    for it in range(iters):
        if not fixed:
            # re-solve the states for the current parameters (two Gauss-Newton steps per frame)
            for _ in range(2):
                Jx, _ = jac(par, th, *mc, W, states_only=True)
                R = r.reshape(len(r), -1)
                dx = -np.linalg.solve(np.einsum("nij,nik->njk", Jx, Jx) + 1e-12 * np.eye(26),
                                      np.einsum("nij,ni->nj", Jx, R)[:, :, None])[..., 0]
                X = X + dx[:, :13] + 1j * dx[:, 13:]
                r, *mc = residual(par, th, X, W)
            cost = float(np.sum(r * r))
        Jx, Jt = jac(par, th, *mc, W)
        R = r.reshape(len(r), -1)
        Axx = np.einsum("nij,nik->njk", Jx, Jx)
        Axt = np.einsum("nij,nik->njk", Jx, Jt)
        if not fixed and NEWTON:
            Axt = Axt + cross_term(par, th, *mc, W, R)
        gx = np.einsum("nij,ni->nj", Jx, R)
        Att = np.einsum("nij,nik->jk", Jt, Jt)
        gt = np.einsum("nij,ni->j", Jt, R)
        while True:
            D = Axx + lam * (np.einsum("nii->ni", Axx)[:, :, None] * np.eye(26) + 1e-12 * np.eye(26))
            Di = np.linalg.inv(D)
            if fixed:
                dt = np.zeros(par.n)
            else:
                DA = Di @ Axt                                                    # (n, 26, P)
                H = Att + lam * np.diag(np.diag(Att)) - np.tensordot(Axt, DA, axes=([0, 1], [0, 1]))
                g = gt - np.tensordot(DA, gx, axes=([0, 1], [0, 1]))
                dt = -np.linalg.solve(H, g)
            dx = -np.einsum("njk,nk->nj", Di, gx + Axt @ dt)
            th2 = th + dt
            X2 = X + dx[:, :13] + 1j * dx[:, 13:]
            r2, *mc2 = residual(par, th2, X2, W)
            c2 = float(np.sum(r2 * r2))
            if c2 < cost:
                rel = (cost - c2) / cost
                th, X, r, mc, cost = th2, X2, r2, mc2, c2
                lam = max(lam / 4, 1e-10)
                break
            lam *= 8
            if lam > 1e8:
                rel = 0.0
                break
        if VERB:
            print(f"    lm it {it} cost/frame {cost / len(r):.3f} lam {lam:.1e}", flush=True)
        if rel < tol or lam > 1e8:
            break
    if fixed:
        return th, X, cost, r, None
    DA = np.linalg.inv(Axx + 1e-12 * np.eye(26)) @ Axt
    Hu = Att - np.tensordot(Axt, DA, axes=([0, 1], [0, 1]))
    return th, X, cost, r, np.linalg.pinv(Hu)


def exact_states(par, th, W):
    """With the rectangular residual the states enter linearly: per-frame weighted least squares."""
    Z, Amap, *_ = model_channels(par, th, np.zeros((1, 13), complex))
    m = MEAS[FR]
    # r = [Re(1 - h/m) / sm, Im(1 - h/m) / sa] W, h = Amap x  ->  real linear system in (Re x, Im x)
    G = Amap[None] / m[:, :, None]                       # h / m = G x
    Gr = np.concatenate([G, 1j * G], 2)
    A = np.concatenate([Gr.real / SM[None, :, None], Gr.imag / SA[None, :, None]], 1) * W[:, :, None]
    b = np.concatenate([np.ones_like(m.real) / SM, np.zeros_like(m.real)], 1) * W
    Q, R = np.linalg.qr(A)
    x = np.linalg.solve(R, Q.transpose(0, 2, 1) @ b[:, :, None])[..., 0]
    return x[:, :13] + 1j * x[:, 13:]


def vp(par, th, W, iters=40, tol=1e-10):
    """Variable projection: the states are re-solved exactly for every parameter value, and the step
    uses the Golub-Pereyra Gram matrix  H = Jt'(I - Px)Jt + C' Axx^-1 C  with C = sum r d2r/dx dtheta."""
    lam = 1e-3
    X = exact_states(par, th, W)
    r, *mc = residual(par, th, X, W)
    cost = float(np.sum(r * r))
    for it in range(iters):
        Jx, Jt = jac(par, th, *mc, W)
        R = r.reshape(len(r), -1)
        Axx = np.einsum("nij,nik->njk", Jx, Jx)
        Axt = np.einsum("nij,nik->njk", Jx, Jt)
        C = cross_term(par, th, *mc, W, R)
        Ai = np.linalg.inv(Axx + 1e-12 * np.eye(26))
        H = (np.einsum("nij,nik->jk", Jt, Jt) - np.tensordot(Axt, Ai @ Axt, axes=([0, 1], [0, 1]))
             + np.tensordot(C, Ai @ C, axes=([0, 1], [0, 1])))
        g = np.einsum("nij,ni->j", Jt, R) - np.tensordot(Ai @ Axt, np.einsum("nij,ni->nj", Jx, R), axes=([0, 1], [0, 1]))
        rel = 0.0
        while True:
            dt = -np.linalg.solve(H + lam * np.diag(np.diag(H)), g)
            th2 = th + dt
            X2 = exact_states(par, th2, W)
            r2, *mc2 = residual(par, th2, X2, W)
            c2 = float(np.sum(r2 * r2))
            if c2 < cost:
                rel = (cost - c2) / cost
                th, X, r, mc, cost = th2, X2, r2, mc2, c2
                lam = max(lam / 5, 1e-12)
                break
            lam *= 10
            if lam > 1e10:
                break
        if VERB:
            print(f"    vp it {it} cost/frame {cost / len(r):.4f} lam {lam:.1e}", flush=True)
        if lam > 1e10 or (rel < tol and lam < 1e-3):
            break
    Hu = np.einsum("nij,nik->jk", Jt, Jt) - np.tensordot(Axt, Ai @ Axt, axes=([0, 1], [0, 1]))
    return th, X, cost, r, np.linalg.pinv(Hu)


def start_states(frames):
    return MEAS[frames][:, [k for k, (c, b) in enumerate(CH) if c == "V"]].copy()


def linear_states(par, th, frames, W):
    """Per-frame weighted linear least squares of the measured phasors on the 13 free voltages
    (isotropic weights 1 / (|m| sigma)); the start of the log-residual Gauss-Newton."""
    Z, Amap, *_ = model_channels(par, th, start_states(frames[:1]))
    m = MEAS[frames]
    w = W[:, :NC] / (np.abs(m) * np.hypot(SM, SA))
    A = w[:, :, None] * Amap[None]
    rhs = (w * m)[:, :, None]
    AH = np.conj(A).transpose(0, 2, 1)
    return np.linalg.solve(AH @ A + 1e-14 * np.eye(13), AH @ rhs)[..., 0]


def equation_error(par, frames, W):
    """Linear least squares of measured currents on measured voltages for the corridor admittances."""
    th = par.nominal()
    X = start_states(frames)
    Y0, terms = Par(par.active, False).branch_terms(th)
    Z, Amap, Vf, terms, cal = model_channels(par, th, X)
    rows = []
    rhs = []
    lin = [pi for pi, *_ in terms if pi < par.nbr and (CANDS_BY_PI(par, pi))]
    # residual of currents with the linear parameters at zero
    th0 = th.copy()
    for pi in set(lin):
        th0[pi] = 0.0
    Z0, *_ = model_channels(par, th0, X)
    meas_I = MEAS[frames][:, 13:]
    y = (meas_I - Z0[:, 13:]).reshape(-1)
    cols = []
    for pi in sorted(set(lin)):
        e = th0.copy(); e[pi] = 1.0
        Z1, *_ = model_channels(par, e, X)
        cols.append((Z1[:, 13:] - Z0[:, 13:]).reshape(-1))
    A = np.stack(cols, 1)
    wt = 1.0 / np.abs(meas_I).reshape(-1).clip(1e-3)
    Ar = np.concatenate([(A * wt[:, None]).real, (A * wt[:, None]).imag])
    yr = np.concatenate([(y * wt).real, (y * wt).imag])
    sol = np.linalg.lstsq(Ar, yr, rcond=None)[0]
    for pi, v in zip(sorted(set(lin)), sol):
        th[pi] = v
    return th


def CANDS_BY_PI(par, pi):
    for (key, v) in par.idx.items():
        if v == pi and key[0] != "cal":
            return key[1] in ("g", "b", "k")
    return False


# ------------------------------------------------------------------ malfunctions
def stale_runs():
    out = []
    for k, (c, b) in enumerate(CH):
        same = np.flatnonzero(np.all(data[1:, [cols.index(f"{c}{b}_mag"), cols.index(f"{c}{b}_ang")]] ==
                                     data[:-1, [cols.index(f"{c}{b}_mag"), cols.index(f"{c}{b}_ang")]], axis=1)) + 1
        if len(same):
            for run in np.split(same, np.flatnonzero(np.diff(same) > 1) + 1):
                out.append((b, k, int(run[0]), int(run[-1])))
    return out


def pmu_chi2(r):
    """Per-frame chi-square of each PMU's four residual components: (T, 13)."""
    q = r * r
    out = np.zeros((len(r), len(PMU)))
    for j, b in enumerate(PMU):
        ks = [k for k, (c, bb) in enumerate(CH) if bb == b]
        out[:, j] = q[:, ks].sum(1) + q[:, [NC + k for k in ks]].sum(1)
    return out


def changepoints(x, a, e, span=120):
    """Refine the two ends of a run by the CUSUM change-point estimate on the excess statistic."""
    lo, hi = max(0, a - span), min(len(x) - 1, e + span)
    seg = x[lo:hi + 1]
    best = (-np.inf, a, e)
    c = np.r_[0, np.cumsum(seg)]
    n = len(seg)
    tot = c[-1]
    for i in range(0, max(1, a - lo + 40), 1):
        for j in range(max(i + 5, e - lo - 40), n):
            inside = c[j + 1] - c[i]
            L = j - i + 1
            out_ = tot - inside
            if n - L <= 0:
                continue
            g = inside * inside / L + out_ * out_ / (n - L) - tot * tot / n
            if g > best[0]:
                best = (g, lo + i, lo + j)
    return best[1], best[2]


def detect(r, mask_frames):
    q = pmu_chi2(r)
    med = np.median(q[~mask_frames], 0)
    z = q / med
    excess = z - 1.0
    found = []
    for j, b in enumerate(PMU):
        x = np.where(mask_frames[:, j], 0.0, excess[:, j])
        # Page CUSUM with reference 1.5 (in units of the median statistic); a run starts where the sum leaves
        # zero and ends where it peaks; it alarms if the peak exceeds 60
        s_, start, peak, runs = 0.0, None, (0.0, None), []
        for t in range(T + 1):
            s_ = max(0.0, s_ + (x[t] if t < T else -np.inf) - 1.5)
            if s_ > 0 and start is None:
                start, peak = t, (0.0, None)
            if s_ > peak[0]:
                peak = (s_, t)
            if s_ == 0 and start is not None:
                if peak[0] > 60:
                    runs.append([start, peak[1]])
                start = None
        for a, e in runs:
            if e - a + 1 < 15:
                continue
            a2, e2 = changepoints(x, a, e)
            # extend each end while the local mean of the normalised statistic stays significant
            # (ramped and out-and-back timing errors fade in gently)
            zz = np.where(mask_frames[:, j], 1.0, z[:, j])
            while a2 > 8 and zz[a2 - 9:a2].mean() > 2.0:
                a2 -= 1
            while e2 < T - 9 and zz[e2 + 1:e2 + 10].mean() > 2.0:
                e2 += 1
            if e2 - a2 + 1 > 1300:                  # longer than any malfunction: model error, not a fault
                continue
            gain = float(np.sum(x[a2:e2 + 1]))
            found.append((gain, b, a2, e2))
    return sorted(found, reverse=True)


# ------------------------------------------------------------------ pipeline
def weights(segs, stale):
    W = np.ones((T, 2 * NC))
    for b, k, a, e in stale:
        W[a:e + 1, [k, NC + k]] = 0.0
    for b, a, e in segs:
        ks = [k for k, (c, bb) in enumerate(CH) if bb == b]
        W[max(0, a - 5):e + 6, ks + [NC + k for k in ks]] = 0.0
    return W


def fit(active, phys, th, X, W, frames, iters=40, fixed=False):
    global FR
    FR = frames
    par = Par(active, phys)
    if RECT:
        if fixed:
            X = exact_states(par, th, W[frames])
            r, *_ = residual(par, th, X, W[frames])
            return par, th, X, float(np.sum(r * r)), r, None
        return (par,) + vp(par, th, W[frames], iters=iters)
    if X is None or len(X) != len(frames):
        X = linear_states(par, th, frames, W[frames])
    return (par,) + lm(par, th, X, W[frames], iters=iters, fixed=fixed)


def to_phys(par_a, th_a, active):
    par = Par(active, True)
    th = par.nominal()
    for key, i in par.idx.items():
        if key[0] == "cal":
            th[i] = th_a[par_a.idx[key]]
        elif key[1] in ("x", "tap") and (key[0], "g") not in par_a.idx:
            th[i] = th_a[par_a.idx[key]]
        else:
            y = complex(th_a[par_a.idx[(key[0], "g")]], th_a[par_a.idx[(key[0], "b")]])
            zz = 1 / y
            th[i] = zz.real if key[1] == "r" else zz.imag
    return par, th


def sub(par_a, th_a, active):
    par = Par(active, False)
    th = np.array([th_a[par_a.idx[k]] for k in par.idx])
    return par, th


def adm_values(par, th):
    return {k: th[i] for k, i in par.idx.items()}


def score_stat(active, vals, cid, X, W, frames):
    """Rao score statistic for adding corridor cid (admittance and charging at zero) to the fit."""
    global FR
    FR = frames
    act2 = sorted(active + [cid])
    par = Par(act2, False)
    th = par.nominal(open_zero=False)
    for k, i in par.idx.items():
        if k in vals:
            th[i] = vals[k]
        elif k[0] == cid:
            th[i] = 0.0
    r, *mc = residual(par, th, X, W[frames])
    Jx, Jt = jac(par, th, *mc, W[frames])
    R = r.reshape(len(r), -1)
    Axx = np.einsum("nij,nik->njk", Jx, Jx)
    Axt = np.einsum("nij,nik->njk", Jx, Jt)
    DA = np.linalg.inv(Axx + 1e-12 * np.eye(26)) @ Axt
    H = np.einsum("nij,nik->jk", Jt, Jt) - np.tensordot(Axt, DA, axes=([0, 1], [0, 1]))
    g = np.einsum("nij,ni->j", Jt, R) - np.tensordot(DA, np.einsum("nij,ni->nj", Jx, R), axes=([0, 1], [0, 1]))
    return float(g @ np.linalg.solve(H, g))


def implausible(par, th, cid):
    c = CANDS[cid - 1]
    y = complex(th[par.idx[(cid, "g")]], th[par.idx[(cid, "b")]])
    if abs(y) < 1e-9:
        return True
    z = 1 / y
    nom = complex(c["nominal"]["r"], c["nominal"]["x"])
    return z.real <= 0 or z.imag <= 0 or abs(z - nom) > 0.5 * abs(nom)


def malfunction_rounds(active, th, W, segs, rounds=4):
    for rnd in range(rounds):
        par, th, X, cost, r, cov = fit(active, False, th, None, W, SUBF, iters=40)
        _, _, _, _, rf, _ = fit(active, False, th, None, W, ALLF, iters=4, fixed=True)
        rf = np.where(W > 0, rf.reshape(T, -1), 0.0)
        mask_frames = stale_mask.copy()
        for b, a, e in segs:
            mask_frames[max(0, a - 5):e + 6, PMU.index(b)] = True
        new = list(segs)
        for gain, b, a, e in detect(rf, mask_frames):
            if gain < 200:
                continue
            if any(a <= y and e >= x for _, x, y in new) or any(a <= f1 and e >= f0 for _, _, f0, f1 in stale):
                continue
            new.append((b, a, e))
            break                                            # one interval per round: masks reshape the rest
        if sorted(new) == sorted(segs):
            break
        segs = sorted(new)
        W = weights(segs, stale)
    return segs, W, th


def topology_search(active, th, W):
    """Forward-backward stepwise search: score test for database-open corridors, Wald / plausibility
    order for the in-service ones, each step confirmed by the likelihood ratio against a cut of 30."""
    par, th, X, cost, r, cov = fit(active, False, th, None, W, SUBF, iters=60)
    for step in range(8):
        vals = adm_values(par, th)
        moved = False
        # additions
        outside = [c["id"] for c in CANDS if not c["known"] and c["id"] not in active and c["model"] in ("line", "series")]
        scores = sorted(((score_stat(active, vals, cid, X, W, SUBF), cid) for cid in outside), reverse=True)
        for sc, cid in scores[:2]:
            if sc < 30:
                break
            act2 = sorted(active + [cid])
            best = None
            for frac in (1.0, 0.0):
                p2 = Par(act2, False)
                t2 = p2.nominal(open_zero=False)
                for k, i in p2.idx.items():
                    if k in vals:
                        t2[i] = vals[k]
                    elif k[0] == cid:
                        t2[i] *= frac if frac > 0 else 1e-3
                res = fit(act2, False, t2, None, W, SUBF, iters=60)
                if best is None or res[3] < best[3]:
                    best = res
            if cost - best[3] > 30:
                active, (par, th, X, cost, r, cov) = act2, best
                moved = True
                break
        if moved:
            continue
        # removals: implausible corridors first, then the weakest by Wald statistic
        cands = []
        for cid in active:
            if CANDS[cid - 1]["model"] not in ("line", "series"):
                continue
            ig, ib = par.idx[(cid, "g")], par.idx[(cid, "b")]
            y = th[[ig, ib]]
            z2 = float(y @ np.linalg.solve(cov[np.ix_([ig, ib], [ig, ib])], y))
            cands.append((0 if implausible(par, th, cid) else 1, z2, cid))
        cands.sort()
        for flag, z2, cid in cands[:3]:
            if flag == 1 and z2 > 900:
                break
            act2 = [a for a in active if a != cid]
            best = None
            for start in ("current", "nameplate"):
                p2 = Par(act2, False)
                t2 = p2.nominal(open_zero=False)
                for k, i in p2.idx.items():
                    if k in vals and (start == "current" or k[0] == "cal"):
                        t2[i] = vals[k]
                res = fit(act2, False, t2, None, W, SUBF, iters=60)
                if best is None or res[3] < best[3]:
                    best = res
            if best[3] - cost < 30:
                active, (par, th, X, cost, r, cov) = act2, best
                moved = True
                break
        if not moved:
            break
    return active, par, th


stale = stale_runs()
stale_mask = np.zeros((T, len(PMU)), bool)
for b, k, a, e in stale:
    stale_mask[a:e + 1, PMU.index(b)] = True
ALLF = np.arange(T)
SUBF = np.arange(0, T, 4)
segs = []
W = weights(segs, stale)
active = [c["id"] for c in CANDS if not c["known"] and c["database_status"] == "closed"]
par0 = Par(active, False)
FR = SUBF
th = equation_error(par0, SUBF, W[SUBF])
for outer in range(3):
    segs, W, th = malfunction_rounds(active, th, W, segs)
    old = list(active)
    active, par, th = topology_search(active, th, W)
    if active == old and outer > 0:
        break
segs, W, th = malfunction_rounds(active, th, W, segs)
parp, thp = to_phys(Par(active, False), th, active)
parp, thp, Xp, cp, rp, covp = fit(active, True, thp, None, W, ALLF, iters=40)
out = {"branches": [], "channels": {}, "corrupted_segments": []}
for c in CANDS:
    ins = c["known"] or c["id"] in active
    e = {"id": c["id"], "in_service": bool(ins)}
    if ins and not c["known"]:
        for k in {"line": ["r", "x"], "series": ["r", "x"], "transformer": ["x", "tap"]}[c["model"]]:
            e[k] = float(thp[parp.idx[(c["id"], k)]])
    out["branches"].append(e)
for k in CALK:
    c, b = CH[k]
    out["channels"][f"{c}{b}"] = {"ratio_error": float(thp[parp.idx[("cal", k, "rho")]]),
                                  "phase_error": float(np.angle(np.exp(1j * thp[parp.idx[("cal", k, "phi")]])))}
st = {}
for b, k, a, e in stale:
    st.setdefault(b, []).append((a, e))
for b, a, e in segs:
    st.setdefault(b, []).append((a, e))
for b, lst in st.items():
    lst.sort()
    merged = []
    for a, e in lst:
        if merged and a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([a, e])
    for a, e in merged:
        out["corrupted_segments"].append({"bus": int(b), "first_frame": int(a), "last_frame": int(e)})
json.dump(out, open(OUT, "w"), indent=1)
print("topology", active, "segments", out["corrupted_segments"], flush=True)
