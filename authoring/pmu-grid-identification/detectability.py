"""Author-side detectability of the planted malfunctions.

usage: python3 detectability.py <rec_dir> <name>

For every malfunction interval, the deviation it adds to the channels of its PMU is whitened with the
noise model and projected orthogonally to the per-frame state Jacobian (what the 13 complex states of
the frame cannot absorb, at the true parameters).  Reported: the largest per-frame chi-square and the
sum over the interval.  Stale data are reported by the exact-repeat length (they are exact copies)."""
import sys

import numpy as np

import model as M
import gen_recording as G


def main():
    d, name = sys.argv[1], sys.argv[2]
    net, tr = M.load_json(f"{d}/network.json"), M.load_json(f"{d}/truth.json")
    Vt = np.load(f"{d}/true_states.npz")["V"]
    Y = M.ybus(net["candidates"], tr["status"], tr["params"], net["shunt"])
    rows = [b - 1 for b in M.PMU_BUSES]
    free = [i for i in range(M.NB) if i != M.Z7]
    S = np.zeros((M.NB, 13), complex); S[free, np.arange(13)] = 1; S[M.Z7] = -Y[M.Z7, free] / Y[M.Z7, M.Z7]
    cal = np.ones(26, complex)
    for k, (c, b) in enumerate(M.CHANNELS):
        if f"{c}{b}" in tr["cal"]:
            e = tr["cal"][f"{c}{b}"]; cal[k] = (1 + e["ratio_error"]) * np.exp(1j * e["phase_error"])
    Ax = np.concatenate([np.eye(13), (Y @ S)[rows]], 0) * cal[:, None]
    sm = np.array([net["noise"][c]["magnitude_rel_std"] for c, _ in M.CHANNELS])
    sa = np.array([net["noise"][c]["angle_std_rad"] for c, _ in M.CHANNELS])
    intervals = tr["corrupted"]
    for iv, mf in zip(intervals, G.MALF[name]):
        f0, f1 = iv["first_frame"], iv["last_frame"]
        L = f1 - f0 + 1
        if mf["kind"] == "stale":
            print(f"  PMU {iv['bus']:2d} {f0}-{f1}: stale {mf['ch']}, {L} exact repeats")
            continue
        V = Vt[f0:f1 + 1]
        z = np.concatenate([V[:, rows], V @ Y.T[:, rows]], 1) * cal
        ks = [k for k, (c, b) in enumerate(M.CHANNELS) if b == mf["bus"]]
        dz = np.zeros_like(z)
        if mf["kind"] == "timing":
            i = np.arange(1, L + 1)
            g = {"step": np.ones(L), "ramp": i / L, "tri": 1 - np.abs(2 * i / (L + 1) - 1)}[mf["shape"]]
            dz[:, ks] = z[:, ks] * (np.exp(1j * mf["size"] * g)[:, None] - 1)
        else:
            k = [k for k in ks if M.CHANNELS[k][0] == mf["ch"]][0]
            dz[:, k] = z[:, k] * mf["size"]
        u = np.exp(-1j * np.angle(z))
        def white(a):
            r = a * u[..., None] if a.ndim == 3 else a * u
            return np.concatenate([r.real / (np.abs(z)[..., None] * sm[:, None] if a.ndim == 3 else np.abs(z) * sm),
                                   r.imag / (np.abs(z)[..., None] * sa[:, None] if a.ndim == 3 else np.abs(z) * sa)], 1)
        w = white(dz)                                                     # (L, 52)
        Jx = np.concatenate([np.broadcast_to(Ax, (L, 26, 13)), 1j * np.broadcast_to(Ax, (L, 26, 13))], 2)
        Wx = white(Jx)
        Q, _ = np.linalg.qr(Wx)
        res = w - (Q @ (Q.transpose(0, 2, 1) @ w[:, :, None]))[..., 0]
        chi = np.sum(res ** 2, 1)
        print(f"  PMU {iv['bus']:2d} {f0}-{f1}: {mf['kind']} {mf.get('shape', mf.get('ch'))} {mf['size']:.2g}: "
              f"unabsorbed chi2 per frame max {chi.max():.0f}, median {np.median(chi):.0f}, sum {chi.sum():.3g}")


if __name__ == "__main__":
    main()
