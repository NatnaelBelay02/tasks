"""Reference solution: PMU network identification (topology, branch parameters, instrument
calibration, malfunction intervals) by robust errors-in-variables maximum likelihood.

Usage: python3 identify.py <data_dir> <output.json>
<data_dir> holds pmu.csv and network.json.

Method: for fixed parameters theta the true phasors of every frame enter linearly, so they are projected
out per frame (variable projection) and the projected, noise-whitened cost is minimised over theta by
Levenberg-Marquardt with the exact Golub-Pereyra Jacobian.  Topology: superset model converged from
several starts (lowest cost kept), then backward elimination on a penalised likelihood.  Malfunctions:
robust (Cauchy) reweighting per frame and PMU, rolling-median residual runs, endpoints refined by fitting
step / ramp / out-and-back shapes to per-frame offset estimates.
"""
import json
import sys
import time

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

DATA, OUT = sys.argv[1], sys.argv[2]
net = json.load(open(f"{DATA}/network.json"))
NB = len(net["buses"])
Z7 = net["zero_injection_bus"] - 1
REF = set(net["reference_pmus"])
FREE = [i for i in range(NB) if i != Z7]
NOI = net["noise"]
CANDS = net["candidates"]
SHUNT = net["shunt"]

# ---------------------------------------------------------------- data
raw = np.genfromtxt(f"{DATA}/pmu.csv", delimiter=",", names=True)
T = len(raw)
CH = [("V", b) for b in net["buses"] if b - 1 != Z7] + [("I", b) for b in net["buses"] if b - 1 != Z7]
Z = np.stack([raw[f"{c}{b}_mag"] * np.exp(1j * raw[f"{c}{b}_ang"]) for c, b in CH], 1)      # (T, 26)
MAG = np.stack([raw[f"{c}{b}_mag"] for c, b in CH], 1); ANG = np.stack([raw[f"{c}{b}_ang"] for c, b in CH], 1)
NCH = len(CH)
CAL = [k for k, (c, b) in enumerate(CH) if b not in REF]          # channels with unknown calibration


def whitener(z):
    """Per-channel 2x2 inverse square-root covariance of polar noise at phasor z: (T, NCH, 2, 2)."""
    out = np.zeros(z.shape + (2, 2))
    for k, (c, _) in enumerate(CH):
        sm, sa = NOI[c]["magnitude_rel_std"], NOI[c]["angle_std_rad"]
        r, a = np.abs(z[:, k]), np.angle(z[:, k])
        ca, sa_ = np.cos(a), np.sin(a)
        out[:, k, 0, 0] = ca / (sm * r); out[:, k, 0, 1] = sa_ / (sm * r)
        out[:, k, 1, 0] = -sa_ / (sa * r); out[:, k, 1, 1] = ca / (sa * r)
    return out


# ---------------------------------------------------------------- parameterisation
class Model:
    """theta = [branch params..., (ratio, phase) for each calibrated channel]."""

    def __init__(self, active, coords):
        self.active = active            # list of candidate ids that are modelled as in service
        self.coords = coords            # "adm" (g, b series admittance) or "phys" (r, x / x, tap)
        self.names, th = [], []
        for c in CANDS:
            if c["known"] or c["id"] not in active:
                continue
            nom = c["nominal"]
            if coords == "adm" and c["model"] in ("line", "series"):
                y = 1 / complex(nom["r"], nom["x"])
                self.names += [(c["id"], "g"), (c["id"], "bs")]; th += [y.real, y.imag]
                if c["model"] == "line":                 # charging multiplier: lets an open line drop out
                    self.names += [(c["id"], "kappa")]; th += [1.0]
            elif c["model"] in ("line", "series"):
                self.names += [(c["id"], "r"), (c["id"], "x")]; th += [nom["r"], nom["x"]]
            elif c["model"] == "transformer":
                self.names += [(c["id"], "x"), (c["id"], "tap")]; th += [nom["x"], nom["tap"]]
            else:
                self.names += [(c["id"], "x")]; th += [nom["x"]]
        self.nb = len(th)
        for k in CAL:
            self.names += [("cal", k, "ratio"), ("cal", k, "phase")]; th += [0.0, 0.0]
        self.theta0 = np.array(th)

    def ybus(self, th):
        Y = np.zeros((NB, NB), complex)
        val = dict(zip(self.names, th))
        for c in CANDS:
            if not (c["known"] or c["id"] in self.active):
                continue
            nom, m = c["nominal"], c["model"]
            g = lambda k: val[(c["id"], k)]
            if c["known"]:
                p = nom
                ys = 1 / complex(p.get("r", 0.0), p["x"]); tap = p.get("tap", 1.0)
            elif m in ("line", "series"):
                ys = complex(g("g"), g("bs")) if self.coords == "adm" else 1 / complex(g("r"), g("x")); tap = 1.0
            elif m == "transformer":
                ys = 1 / complex(0, g("x")); tap = g("tap")
            else:
                ys = 1 / complex(0, g("x")); tap = 1.0
            bc = nom.get("b", 0.0) if m == "line" else 0.0
            if (c["id"], "kappa") in val:
                bc *= g("kappa")
            i, j = c["from_bus"] - 1, c["to_bus"] - 1
            ytt = ys + 0.5j * bc
            Y[i, i] += ytt / tap ** 2; Y[i, j] -= ys / tap; Y[j, i] -= ys / tap; Y[j, j] += ytt
        Y[SHUNT["bus"] - 1, SHUNT["bus"] - 1] += 1j * SHUNT["nominal_b"]
        return Y

    def A(self, th):
        """Real (2*NCH x 2*13) map from free bus voltages [Re, Im interleaved] to channel phasors."""
        Y = self.ybus(th)
        S = np.zeros((NB, len(FREE)), complex)
        S[FREE, np.arange(len(FREE))] = 1
        S[Z7] = -Y[Z7, FREE] / Y[Z7, Z7]
        YS = Y @ S
        cal = np.ones(NCH, complex)
        for idx, k in enumerate(CAL):
            e, ph = th[self.nb + 2 * idx], th[self.nb + 2 * idx + 1]
            cal[k] = (1 + e) * np.exp(1j * ph)
        Ac = np.stack([(S if c == "V" else YS)[b - 1] * cal[k] for k, (c, b) in enumerate(CH)])   # (NCH, 13)
        Ar = np.zeros((2 * NCH, 2 * len(FREE)))
        Ar[0::2, 0::2] = Ac.real; Ar[0::2, 1::2] = -Ac.imag
        Ar[1::2, 0::2] = Ac.imag; Ar[1::2, 1::2] = Ac.real
        return Ar


def whiten(Wh, v):
    """Apply per-channel 2x2 whitening Wh (n, NCH, 2, 2) to v (..., n, 2*NCH) with trailing channel pairs."""
    v2 = v.reshape(v.shape[:-1] + (NCH, 2))
    out = Wh[..., 0] * v2[..., :1] + Wh[..., 1] * v2[..., 1:]
    return out.reshape(v.shape)


def stale_segments():
    segs = []
    for k, (c, b) in enumerate(CH):
        same = (MAG[1:, k] == MAG[:-1, k]) & (ANG[1:, k] == ANG[:-1, k])
        idx = np.flatnonzero(same) + 1
        if len(idx):
            for run in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
                segs.append((b, c, int(run[0]), int(run[-1])))
    return segs


def theta_from(model, adm):
    """Initial theta for `model` from admittance-coordinate values {(id, 'g'|'bs'|...): v} (nominal elsewhere)."""
    th = model.theta0.copy()
    for i, nm in enumerate(model.names):
        if nm[0] == "cal":
            if nm in adm:
                th[i] = adm[nm]
            continue
        cid, key = nm
        if model.coords == "adm" or key in ("x", "tap") and (cid, "g") not in adm:
            if nm in adm:
                th[i] = adm[nm]
            continue
        if (cid, "g") in adm:
            zz = 1 / complex(adm[(cid, "g")], adm[(cid, "bs")])
            th[i] = zz.real if key == "r" else zz.imag
    return th


def ystat(model, theta, cov, cid):
    ig = model.names.index((cid, "g")); ib = model.names.index((cid, "bs"))
    y = np.array([theta[ig], theta[ib]]); C = cov[np.ix_([ig, ib], [ig, ib])]
    return float(np.hypot(*y)), float(np.sqrt(abs(y @ np.linalg.solve(C, y))))


def values(model, theta):
    return dict(zip(model.names, theta))


def PROJ(model, th):
    """Physical bounds: series admittance magnitude <= 1.5 x nameplate, charging multiplier in [0, 1.5]."""
    th = th.copy()
    if model.coords != "adm":
        return th
    for cid in model.active:
        c = CANDS[cid - 1]
        if c["known"] or c["model"] not in ("line", "series"):
            continue
        ig, ib = model.names.index((cid, "g")), model.names.index((cid, "bs"))
        ynom = abs(1 / complex(c["nominal"]["r"], c["nominal"]["x"]))
        mag = np.hypot(th[ig], th[ib])
        if mag > 1.5 * ynom:
            th[ig] *= 1.5 * ynom / mag; th[ib] *= 1.5 * ynom / mag
        if (cid, "kappa") in model.names:
            ik = model.names.index((cid, "kappa")); th[ik] = min(max(th[ik], 0.0), 1.5)
    return th


T0 = time.time()
PMUS = sorted({b for _, b in CH})
ROWS = {b: [k for k, (_, bb) in enumerate(CH) if bb == b] for b in PMUS}
ALL = np.arange(T)


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def rolling_median(x, w):
    pad = np.pad(x, w // 2, mode="edge")
    return np.median(sliding_window_view(pad, w), axis=1)[:len(x)]


def runs_of(flag, min_len=1):
    idx = np.flatnonzero(flag)
    if not len(idx):
        return []
    out = []
    for r in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
        if len(r) >= min_len:
            out.append([int(r[0]), int(r[-1])])
    return out


def merge(runs, gap):
    runs = sorted(runs)
    out = []
    for a, e in runs:
        if out and a - out[-1][1] <= gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([a, e])
    return out


# ---------------------------------------------------------------- exact variable-projection fit
def whiten_T(Wh, v):
    """Apply the transpose of the per-channel 2x2 whitening to v (..., n, 2*NCH)."""
    v2 = v.reshape(v.shape[:-1] + (NCH, 2))
    out = Wh[..., 0, :] * v2[..., :1] + Wh[..., 1, :] * v2[..., 1:]
    return out.reshape(v.shape)


def fit(model, frames, iters=40, verbose=True, weights=None, theta=None, chunk=2000):
    """Levenberg-Marquardt on the variable-projection cost with the exact (Golub-Pereyra) Jacobian:
    J = -(P_perp dA u + A^{+T} dA^T r) per frame, so far-from-optimum steps remain descent directions."""
    th = PROJ(model, (model.theta0 if theta is None else theta).copy())
    zr = np.stack([Z.real, Z.imag], -1).reshape(T, -1)[frames]
    Wh = whitener(Z[frames])
    if weights is not None:
        Wh = Wh * np.sqrt(weights)[:, :, None, None]
    y = whiten(Wh, zr)
    lam = 1e-3
    P = len(th)
    n = len(frames)

    def state(th):
        A = model.A(th)
        At = whiten(Wh[:, None], np.broadcast_to(A.T, (n,) + A.T.shape)).transpose(0, 2, 1)
        Q, R = np.linalg.qr(At)
        u = np.linalg.solve(R, (Q.transpose(0, 2, 1) @ y[:, :, None]))[..., 0]
        res = y - (At @ u[:, :, None])[..., 0]
        return A, Q, R, u, res

    A, Q, R, u, res = state(th)
    cost = np.sum(res ** 2)
    N = None
    for it in range(iters):
        dA = np.empty((P,) + A.shape)
        for p in range(P):
            hh = 1e-7 * max(1.0, abs(th[p]))
            tp, tm = th.copy(), th.copy(); tp[p] += hh; tm[p] -= hh
            dA[p] = (model.A(tp) - model.A(tm)) / (2 * hh)
        N = np.zeros((P, P)); grad = np.zeros(P)
        for s0 in range(0, n, chunk):
            sl = slice(s0, s0 + chunk)
            G = u[sl] @ dA.transpose(0, 2, 1)
            G = whiten(Wh[sl], G).transpose(1, 2, 0)                  # (m, 52, P) whitened dA_p u
            Qc, Rc = Q[sl], R[sl]
            J = -(G - Qc @ (Qc.transpose(0, 2, 1) @ G))
            sT = whiten_T(Wh[sl], res[sl])                             # W^T r   (m, 52)
            x = np.einsum("pij,mi->mjp", dA, sT)                       # dA_p^T W^T r   (m, 26, P)
            yv = np.linalg.solve(Rc.transpose(0, 2, 1), x)             # R^{-T} x
            J -= Qc @ yv                                               # - A^{+T} dA^T r
            Jf = J.reshape(-1, P)
            N += Jf.T @ Jf; grad += Jf.T @ res[sl].reshape(-1)
        c1 = cost
        while True:
            step = np.linalg.solve(N + lam * np.diag(np.diag(N)), -grad)
            th2 = PROJ(model, th + step)
            A2, Q2, R2, u2, res2 = state(th2)
            c2 = np.sum(res2 ** 2)
            if c2 < cost:
                th, A, Q, R, u, res, cost = th2, A2, Q2, R2, u2, res2, c2
                lam = max(lam / 5, 1e-12); break
            lam *= 10
            if lam > 1e10:
                break
        rel = (c1 - cost) / c1
        if verbose:
            print(f"   it {it:2d} cost {cost:.1f} rel decrease {rel:.2e} lam {lam:.0e}", flush=True)
        if lam > 1e10 or (rel < 1e-10 and lam < 1e-3):
            break
    if N is None:
        return th, None, res.reshape(n, NCH, 2), u
    return th, np.linalg.inv(N), res.reshape(n, NCH, 2), u


# ---------------------------------------------------------------- robust fitting
MASKW = 1e-6
PEN = 30.0
DEC = 4
FR = np.arange(0, T, DEC)          # search fits use every DEC-th frame; the final fit uses all frames


def raw_resid(model, theta, u):
    A = model.A(theta)
    zr = np.stack([Z.real, Z.imag], -1).reshape(T, -1)
    return whiten(whitener(Z), zr - u @ A.T).reshape(T, NCH, 2)


def states(model, theta, w):
    Wh = whitener(Z) * np.sqrt(w)[:, :, None, None]
    y = whiten(Wh, np.stack([Z.real, Z.imag], -1).reshape(T, -1))
    A = model.A(theta)
    At = whiten(Wh[:, None], np.broadcast_to(A.T, (T,) + A.T.shape)).transpose(0, 2, 1)
    Q, R = np.linalg.qr(At)
    return np.linalg.solve(R, (Q.transpose(0, 2, 1) @ y[:, :, None]))[..., 0]


def cauchy_fit(model, theta, w0, rounds=6, iters=10, c=3.0):
    """IRLS with Cauchy weights per (frame, PMU), scaled by the median statistic, from the first iteration."""
    u = states(model, theta, w0)
    for k in range(rounds):
        raw = raw_resid(model, theta, u)
        st = np.stack([np.mean(raw[:, ROWS[b], :] ** 2, axis=(1, 2)) for b in PMUS], 1)
        med = np.median(st)
        wb = 1.0 / (1.0 + st / (c * c * med))
        w = w0.copy()
        for j, b in enumerate(PMUS):
            w[:, ROWS[b]] *= wb[:, [j]]
        theta, cov, res, _ = fit(model, FR, iters=iters, weights=w[FR], theta=theta, verbose=False)
        u = states(model, theta, w)
    log(f"   robust fit: median PMU statistic {med:.2f}")
    return theta, cov, w, u


def clean_weights(segs):
    w = w_stale.copy()
    for b, a, e in segs:
        w[max(0, a - 3):e + 4, ROWS[b]] = MASKW
    return w


def converge(act, vals, w, iters):
    m = Model(act, "adm")
    th, cov, res, u = fit(m, FR, iters=iters, weights=w[FR], theta=theta_from(m, vals), verbose=False)
    return m, th, cov, float(np.sum(res ** 2)), u


class CalOnly:
    """Branch parameters frozen at given values; only the channel calibrations are free."""
    def __init__(self, inner, vals):
        self.inner, self.coords, self.active = inner, "cal", inner.active
        full = np.array([vals[n] for n in inner.names])
        self.fixed = full[:inner.nb]
        self.names = inner.names[inner.nb:]; self.nb = 0
        self.theta0 = full[inner.nb:].copy()

    def A(self, th):
        return self.inner.A(np.r_[self.fixed, th])


def staged_start(w):
    """Calibrations first (database topology at nameplate), then the database topology, then the superset."""
    db = [c["id"] for c in CANDS if not c["known"] and c["database_status"] == "closed"]
    mdb = Model(db, "adm"); vdb = values(mdb, mdb.theta0)
    mc = CalOnly(mdb, vdb)
    thc, _, _, _ = fit(mc, FR, iters=40, weights=w[FR], theta=mc.theta0, verbose=False)
    vdb.update(dict(zip(mc.names, thc)))
    m, th, cov, c, u = converge(db, vdb, w, 40)
    act, m0, th0 = superset_start()
    sv = values(m0, th0); sv.update({k: v for k, v in values(m, th).items() if k in sv})
    for _ in range(2):
        m, th, cov, c, u = converge(act, sv, w, 40); sv = values(m, th)
    log(f"   staged start: superset cost/frame {c / len(FR):.3f}")
    return sv, c


class AdmOnly:
    """Channel calibrations frozen; only the branch parameters are free."""
    def __init__(self, inner, th):
        self.inner, self.coords, self.active, self.nb = inner, "adm", inner.active, inner.nb
        self.fixed = th[inner.nb:].copy(); self.names = inner.names[:inner.nb]; self.theta0 = th[:inner.nb].copy()

    def A(self, th):
        return self.inner.A(np.r_[th, self.fixed])


def admittance_first_start(w):
    """Superset admittances fitted with all calibrations frozen at zero (database-open corridors start at 0)."""
    act, m0, th0 = superset_start()
    ma = AdmOnly(m0, th0)
    tha, _, _, _ = fit(ma, FR, iters=60, weights=w[FR], theta=ma.theta0, verbose=False)
    return values(m0, np.r_[tha, ma.fixed])


def best_superset(w, extra=()):
    """Converge the superset from several starts and keep the lowest cost (the cost surface has
    compensating local minima: calibrations and neighbouring corridors absorb a wrong initial topology)."""
    act, m0, th0 = superset_start()
    starts = list(extra) + [("nameplate", values(m0, th0)), ("staged", staged_start(w)[0]),
                            ("admittance-first", admittance_first_start(w))]
    best = None
    for name, vals in starts:
        r = converge(act, vals, w, 60)
        log(f"   superset fit ({name} start) cost/frame {r[3] / len(FR):.4f}, at bound {at_bounds(r[0], r[1])}")
        if best is None or r[3] < best[3]:
            best = r
    return act, best


def at_bounds(m, th):
    out = []
    for cid in m.active:
        c = CANDS[cid - 1]
        if c["known"] or c["model"] not in ("line", "series"):
            continue
        ynom = abs(1 / complex(c["nominal"]["r"], c["nominal"]["x"]))
        if ystat(m, th, np.eye(len(th)), cid)[0] > 1.499 * ynom:
            out.append(cid)
    return out


def superset_start():
    act = [c["id"] for c in CANDS if not c["known"]]
    m = Model(act, "adm")
    th = m.theta0.copy()
    for i, nm in enumerate(m.names):
        if nm[0] != "cal" and CANDS[nm[0] - 1]["database_status"] == "open" and nm[1] in ("g", "bs", "kappa"):
            th[i] = 0.0
    return act, m, th


def select_topology(w, start_vals):
    """Superset fit on clean channels (warm start from the robust superset estimate), then backward
    elimination on the penalised likelihood."""
    active, (m, th, cov, c0, u) = best_superset(w, [("warm", start_vals)])
    weak = [cid for cid in active if CANDS[cid - 1]["model"] in ("line", "series") and ystat(m, th, cov, cid)[1] < 4]
    if weak:
        act2 = [a for a in active if a not in weak]
        m2, th2, cov2, c2, u2 = converge(act2, values(m, th), w, 40)
        d = (c2 + PEN * len(act2)) - (c0 + PEN * len(active))
        log(f"   batch removal {weak}: penalised change {d:.1f}")
        if d < 0:
            active, m, th, cov, c0, u = act2, m2, th2, cov2, c2, u2
    while True:
        zs = sorted((ystat(m, th, cov, cid)[1], cid) for cid in active
                    if CANDS[cid - 1]["model"] in ("line", "series"))
        best = None
        for z, cid in zs[:3]:
            if z > 25:
                break
            act2 = [a for a in active if a != cid]
            m2, th2, cov2, c2, u2 = converge(act2, values(m, th), w, 40)
            d = (c2 + PEN * len(act2)) - (c0 + PEN * len(active))
            log(f"   remove {cid:2d} (z {z:.1f}): penalised change {d:.1f}")
            if d < 0 and (best is None or d < best[0]):
                best = (d, act2, m2, th2, cov2, c2, u2)
        if best is None:
            break
        _, active, m, th, cov, c0, u = best
        log(f"   -> topology {active}")
    # Physical plausibility.  The instruction bounds every unknown r and x to 10 % of nameplate; a corridor
    # fitted with a non-positive resistance or reactance, or more than 50 % from nameplate, is a
    # compensating artefact of a false minimum (an open corridor kept alive by its neighbours).  Refit
    # without it from fresh starts and keep the result if the penalised likelihood improves.
    while True:
        best = None
        for cid in active:
            c = CANDS[cid - 1]
            if c["model"] not in ("line", "series"):
                continue
            y = complex(th[m.names.index((cid, "g"))], th[m.names.index((cid, "bs"))])
            nom = complex(c["nominal"]["r"], c["nominal"]["x"])
            zz = 1 / y if abs(y) > 1e-9 else complex(np.inf, np.inf)
            if zz.real > 0 and zz.imag > 0 and abs(zz - nom) <= 0.5 * abs(nom):
                continue
            act2 = [a for a in active if a != cid]
            fresh = values(Model(act2, "adm"), Model(act2, "adm").theta0)
            fresh.update({k: v for k, v in values(m, th).items() if k[0] == "cal"})
            for name, vals in (("current", values(m, th)), ("nameplate", fresh)):
                r = converge(act2, vals, w, 60)
                d = (r[3] + PEN * len(act2)) - (c0 + PEN * len(active))
                log(f"   implausible {cid} (z = {zz.real:.4f} + {zz.imag:.4f}j): removal from {name} start, penalised change {d:.1f}")
                if d < 0 and (best is None or d < best[0]):
                    best = (d, act2) + r
        if best is None:
            break
        _, active, m, th, cov, c0, u = best
        log(f"   -> topology {active}")
    return active, m, th, cov, u


# ---------------------------------------------------------------- malfunction detection
def detect(raw, stale_mask):
    """Candidate runs per PMU from rolling-median residual statistics (PMU level and channel level)."""
    cand = {}
    for b in PMUS:
        ks = ROWS[b]
        okk = ~stale_mask[:, ks].any(1)
        st = np.mean(raw[:, ks, :] ** 2, axis=(1, 2))
        st = st / np.median(st[okk])
        flag = rolling_median(np.where(okk, st, 1.0), 31) > 3.0
        for k in ks:
            sc = np.mean(raw[:, k, :] ** 2, axis=1)
            sc = sc / np.median(sc[~stale_mask[:, k]])
            flag |= rolling_median(np.where(~stale_mask[:, k], sc, 1.0), 31) > 4.0
        r = [x for x in merge(runs_of(flag), 40) if x[1] - x[0] + 1 >= 15]
        if r:
            cand[b] = r
    return cand


def offset_series(model, theta, u, w, b, lo, hi):
    """Per-frame WLS estimates of a common phase offset of PMU b and of a ratio offset of each of its channels."""
    fr = np.arange(lo, hi + 1)
    A = model.A(theta)
    Wh = whitener(Z[fr])
    ww = w[fr].copy()
    ww[:, ROWS[b]] = (~stale_mask[fr][:, ROWS[b]]).astype(float)
    Wh = Wh * np.sqrt(ww)[:, :, None, None]
    h = (u[fr] @ A.T).reshape(len(fr), NCH, 2)
    hc = h[..., 0] + 1j * h[..., 1]
    cols = []
    dz = np.zeros((len(fr), NCH), complex)
    dz[:, ROWS[b]] = 1j * hc[:, ROWS[b]]
    cols.append(dz)
    for k in ROWS[b]:
        dz = np.zeros((len(fr), NCH), complex); dz[:, k] = hc[:, k]; cols.append(dz)
    X = np.concatenate([np.broadcast_to(A, (len(fr),) + A.shape)] +
                       [np.stack([c.real, c.imag], -1).reshape(len(fr), -1, 1) for c in cols], axis=2)
    Xw = np.stack([whiten(Wh, X[:, :, j]) for j in range(X.shape[2])], 2)
    yw = whiten(Wh, np.stack([Z[fr].real, Z[fr].imag], -1).reshape(len(fr), -1))
    Q, R = np.linalg.qr(Xw)
    d = np.abs(np.diagonal(R, axis1=1, axis2=2))
    R = R + np.eye(R.shape[1]) * (d < 1e-9 * d.max())[:, :, None]       # guard rank-deficient frames
    beta = np.linalg.solve(R, (Q.transpose(0, 2, 1) @ yw[:, :, None]))[..., 0]
    Rinv = np.linalg.inv(R)
    var = np.einsum("tij,tij->tj", Rinv, Rinv)
    m = A.shape[1]
    return fr, beta[:, m:], var[:, m:]


def shapes(n, a, e):
    L = e - a + 1
    i = np.arange(1, L + 1)
    g_step = np.ones(L)
    g_ramp = i / L
    g_tri = 1 - np.abs(2 * i / (L + 1) - 1)
    return {"step": g_step, "ramp": g_ramp, "tri": g_tri}


def best_shape(fr, x, var, a0, e0, kinds):
    """Grid search of (shape, first, last) maximising the weighted SSE reduction of a single-amplitude model."""
    wt = 1.0 / var
    best = (-1, None)

    def scan(a_range, e_range):
        nonlocal best
        for a in a_range:
            for e in e_range:
                if e - a < 4:
                    continue
                ia, ie = a - fr[0], e - fr[0]
                if ia < 0 or ie >= len(fr):
                    continue
                xs, ws = x[ia:ie + 1], wt[ia:ie + 1]
                for kname, g in shapes(len(fr), a, e).items():
                    if kname not in kinds:
                        continue
                    num = np.sum(ws * xs * g); den = np.sum(ws * g * g)
                    gain = num * num / den
                    if gain > best[0]:
                        best = (gain, (kname, a, e, num / den))
    scan(range(a0 - 150, a0 + 61, 3), range(e0 - 60, e0 + 151, 3))
    _, (k, a, e, s) = best
    scan(range(a - 4, a + 5), range(e - 4, e + 5))
    return best


def refine(model, theta, u, w, cand):
    """Endpoint refinement by fitting malfunction shapes to the per-frame offset estimates."""
    segs = []
    for b, runs in cand.items():
        for a0, e0 in runs:
            lo, hi = max(0, a0 - 250), min(T - 1, e0 + 250)
            fr, beta, var = offset_series(model, theta, u, w, b, lo, hi)
            cands = [best_shape(fr, beta[:, 0], var[:, 0], a0, e0, ("step", "ramp", "tri"))]
            for j in range(1, beta.shape[1]):
                cands.append(best_shape(fr, beta[:, j], var[:, j], a0, e0, ("step",)))
            gain, (kind, a, e, s) = max(cands, key=lambda c: c[0])
            log(f"   malfunction PMU {b}: run {a0}-{e0} -> {kind} {a}-{e} amplitude {s:.2e} (chi2 gain {gain:.0f})")
            segs.append((b, a, e, gain))
    return segs


# ---------------------------------------------------------------- pipeline
log(f"{T} frames, {NCH} channels")
stale = stale_segments()
stale_mask = np.zeros((T, NCH), bool)
for b, c, f0, f1 in stale:
    stale_mask[f0:f1 + 1, CH.index((c, b))] = True
log("stale runs:", stale)
w_stale = np.where(stale_mask, MASKW, 1.0)
active, mA, th = superset_start()
log("superset:", active)
_, (mA, th, _, _, _) = best_superset(w_stale)
th, cov, w, u = cauchy_fit(mA, th, w_stale)
segs = None
for rnd in range(3):
    raw = raw_resid(mA, th, u)
    cand = detect(raw, stale_mask)
    for b in list(cand):
        cand[b] = [r for r in cand[b] if not any(sb == b and r[0] <= f1 and r[1] >= f0 for sb, _, f0, f1 in stale)]
    log(f"round {rnd}: candidate runs", cand)
    found = refine(mA, th, u, w, cand)
    acc = [(b, f0, f1) for b, _, f0, f1 in stale]
    new_segs = []
    for b, a, e, g in sorted(found, key=lambda t: -t[3]):
        if g <= 200.0 or any(a <= y and e >= x for _, x, y in acc):
            continue
        acc.append((b, a, e)); new_segs.append((b, a, e))
    new_segs.sort()
    log("   accepted:", new_segs)
    if new_segs == segs:
        break
    segs = new_segs
    if rnd == 0:
        sup_vals = values(mA, th)          # robust superset estimate
    active, mA, th, cov, u = select_topology(clean_weights(segs), sup_vals)
    th, cov, w, u = cauchy_fit(mA, th, w_stale, rounds=3)
theta_adm = values(mA, th)

# final physical fit on clean channels
w_clean = clean_weights(segs)
m2 = Model(active, "phys")
th2, cov2, _, _ = fit(m2, ALL, iters=50, weights=w_clean, theta=theta_from(m2, theta_adm), verbose=False)
status = {c["id"]: (c["known"] or c["id"] in active) for c in CANDS}
allsegs = {}
for b, c, f0, f1 in stale:
    allsegs.setdefault(b, []).append([f0, f1])
for b, a, e in segs:
    allsegs.setdefault(b, []).append([a, e])
out = {"branches": [], "channels": {}, "corrupted_segments": []}
val = dict(zip(m2.names, th2))
for c in CANDS:
    entry = {"id": c["id"], "in_service": bool(status[c["id"]])}
    if status[c["id"]] and not c["known"]:
        for key in {"line": ["r", "x"], "series": ["r", "x"], "transformer": ["x", "tap"]}[c["model"]]:
            entry[key] = float(val[(c["id"], key)])
    out["branches"].append(entry)
for idx, k in enumerate(CAL):
    c, b = CH[k]
    out["channels"][f"{c}{b}"] = {"ratio_error": float(th2[m2.nb + 2 * idx]),
                                  "phase_error": float(np.angle(np.exp(1j * th2[m2.nb + 2 * idx + 1])))}
for b in sorted(allsegs):
    for a, e in merge(allsegs[b], 0):
        out["corrupted_segments"].append({"bus": int(b), "first_frame": int(a), "last_frame": int(e)})
# Report within the ranges the instruction states for every true value (r, x of lines and series branches
# within 10 % of nameplate, transformer x within 5 %, tap within 0.019; |ratio| <= 0.01, |phase| <= 0.0105).
# The truth lies inside, so clipping can only reduce an error.  It matters for weakly observed corridors,
# whose likelihood is curved on the scale of their uncertainty (a fitted r can then land far outside).
for entry in out["branches"]:
    c = CANDS[entry["id"] - 1]
    for key in ("r", "x", "tap"):
        if key in entry:
            nom = c["nominal"][key]
            half = 0.019 if key == "tap" else (0.05 if c["model"] == "transformer" else 0.10) * abs(nom)
            entry[key] = float(min(max(entry[key], nom - half), nom + half))
for ch in out["channels"].values():
    ch["ratio_error"] = float(min(max(ch["ratio_error"], -0.01), 0.01))
    ch["phase_error"] = float(min(max(ch["phase_error"], -0.0105), 0.0105))
json.dump(out, open(OUT, "w"), indent=1)
log("wrote", OUT, "topology", active, "segments", out["corrupted_segments"])
