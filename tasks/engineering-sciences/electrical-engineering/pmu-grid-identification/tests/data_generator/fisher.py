"""Author-side Cramer-Rao bound and tolerances (4 sigma) of a recording.

usage: python3 fisher.py <rec_dir> [out_tolerances.json]

The measurement of channel k in frame t is h_k = c_k z_k(theta, x_t), c_k = (1 + rho_k) exp(j phi_k), with
polar noise (relative magnitude sigma_m, angle sigma_a) around h_k.  The states x_t (real and imaginary
parts of the 13 free bus voltages, already rotated by the frequency drift) are nuisance parameters and are
removed per frame by a Schur complement:  F = sum_t J_th^T (I - P_x) J_th  over the frames free of
malfunctions, with whitened Jacobians evaluated at the true values.  Out-of-service candidates get a
counterfactual bound with the candidate added in service at nameplate (never used for grading)."""
import json
import sys

import numpy as np

import model as M


def unknowns(cands, status):
    names = []
    for c in cands:
        if c["known"] or not status[str(c["id"])]:
            continue
        keys = {"line": ["r", "x"], "series": ["r", "x"], "transformer": ["x", "tap"], "reactance": ["x"]}[c["model"]]
        names += [(c["id"], k) for k in keys]
    return names


def fisher(net, status, params, Vt, frames):
    cands = net["candidates"]
    names = unknowns(cands, status)
    cal = np.ones(26, complex)
    calnames = []
    for k, (c, b) in enumerate(M.CHANNELS):
        key = f"{c}{b}"
        if key in params["cal"]:
            e = params["cal"][key]
            cal[k] = (1 + e["ratio_error"]) * np.exp(1j * e["phase_error"])
            calnames.append((key, k))
    Y = M.ybus(cands, status, params["branch"], net["shunt"])
    free = [i for i in range(M.NB) if i != M.Z7]
    S = np.zeros((M.NB, 13), complex)
    S[free, np.arange(13)] = 1
    S[M.Z7] = -Y[M.Z7, free] / Y[M.Z7, M.Z7]
    rows = [b - 1 for b in M.PMU_BUSES]
    # dY/dp by central differences of the exact builder (Y is smooth in r, x, tap)
    dY = []
    for cid, key in names:
        base = dict(params["branch"][str(cid)])
        h = 1e-6 * abs(base[key])
        pp, pm = dict(params["branch"]), dict(params["branch"])
        pp[str(cid)] = dict(base, **{key: base[key] + h})
        pm[str(cid)] = dict(base, **{key: base[key] - h})
        dY.append((M.ybus(cands, status, pp, net["shunt"]) - M.ybus(cands, status, pm, net["shunt"])) / (2 * h))
    nP = len(names) + 2 * len(calnames)
    sm = np.array([net["noise"][c]["magnitude_rel_std"] for c, _ in M.CHANNELS])
    sa = np.array([net["noise"][c]["angle_std_rad"] for c, _ in M.CHANNELS])
    YS = (Y @ S)[rows]                                   # 13 x 13: free voltages -> PMU currents
    Ax = np.zeros((26, 13), complex)
    Ax[:13] = np.eye(13)
    Ax[13:] = YS
    Ax = Ax * cal[:, None]
    F = np.zeros((nP, nP))
    for s0 in range(0, len(frames), 1500):
        fr = frames[s0:s0 + 1500]
        V = Vt[fr]                                       # (n, 14) true rotated voltages
        n = len(fr)
        z = np.concatenate([V[:, rows], V @ Y.T[:, rows]], 1)          # true channel phasors (n, 26)
        m = z * cal
        # complex Jacobians
        Jt = np.zeros((n, 26, nP), complex)
        for p, dYp in enumerate(dY):
            Jt[:, 13:, p] = (V @ dYp.T)[:, rows] * cal[13:]
        for q, (key, k) in enumerate(calnames):
            Jt[:, k, len(names) + 2 * q] = z[:, k] * np.exp(1j * np.angle(cal[k]))     # d/d rho
            Jt[:, k, len(names) + 2 * q + 1] = 1j * m[:, k]                           # d/d phi
        # whitening at the measured phasor: rotate by -angle, scale by 1/(|m| sigma)
        u = np.exp(-1j * np.angle(m))[:, :, None]
        def white(Jc):
            r = Jc * u
            return np.concatenate([r.real / (np.abs(m)[:, :, None] * sm[None, :, None]),
                                   r.imag / (np.abs(m)[:, :, None] * sa[None, :, None])], 1)
        Wt = white(Jt)                                                        # (n, 52, nP)
        Jx = np.concatenate([np.broadcast_to(Ax, (n, 26, 13)), 1j * np.broadcast_to(Ax, (n, 26, 13))], 2)
        Wx = white(Jx)                                                        # (n, 52, 26)
        Q, _ = np.linalg.qr(Wx)
        Rt = Wt - Q @ (Q.transpose(0, 2, 1) @ Wt)
        F += np.einsum("tij,tik->jk", Rt, Rt)
    return names, calnames, F


def tolerances(net, truth, Vt):
    bad = np.zeros(len(Vt), bool)
    for s in truth["corrupted"]:
        bad[s["first_frame"]:s["last_frame"] + 1] = True
    frames = np.flatnonzero(~bad)
    params = {"branch": truth["params"], "cal": truth["cal"]}
    names, calnames, F = fisher(net, truth["status"], params, Vt, frames)
    C = np.linalg.inv(F)
    sd = np.sqrt(np.diag(C))
    ev = np.linalg.eigvalsh(F)
    out = {"branches": {}, "channels": {}}
    for (cid, key), v in zip(names, sd):
        out["branches"].setdefault(str(cid), {})[key] = 4 * float(v)
    for q, (key, k) in enumerate(calnames):
        out["channels"][key] = {"ratio_error": 4 * float(sd[len(names) + 2 * q]),
                                "phase_error": 4 * float(sd[len(names) + 2 * q + 1])}
    # counterfactual bounds for out-of-service candidates (nameplate, in service)
    for c in net["candidates"]:
        if c["known"] or truth["status"][str(c["id"])]:
            continue
        st = dict(truth["status"]); st[str(c["id"])] = True
        pb = dict(truth["params"]); pb[str(c["id"])] = {k: v for k, v in c["nominal"].items() if k in ("r", "x", "tap")}
        nm2, _, F2 = fisher(net, st, {"branch": pb, "cal": truth["cal"]}, Vt, frames)
        sd2 = np.sqrt(np.diag(np.linalg.inv(F2)))
        for (cid, key), v in zip(nm2, sd2):
            if cid == c["id"]:
                out["branches"].setdefault(str(cid), {})[key] = 4 * float(v)
    out["branches"] = dict(sorted(out["branches"].items(), key=lambda kv: int(kv[0])))
    return out, ev, len(names) + 2 * len(calnames)


def main():
    d = sys.argv[1]
    net = M.load_json(f"{d}/network.json")
    truth = M.load_json(f"{d}/truth.json")
    Vt = np.load(f"{d}/true_states.npz")["V"]
    out, ev, n = tolerances(net, truth, Vt)
    print(d, "unknowns", n, "Fisher eigenvalues min %.3g, condition %.3g" % (ev.min(), ev.max() / ev.min()))
    if len(sys.argv) > 2:
        json.dump(out, open(sys.argv[2], "w"), indent=1)


if __name__ == "__main__":
    main()
