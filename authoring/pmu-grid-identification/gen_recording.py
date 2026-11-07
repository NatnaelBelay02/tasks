"""Author-side generator of the synchrophasor recordings (v3 load model).

usage: python3 gen_recording.py <name> <network.json> <truth.json> <out_dir>

Every frame (2 per second, 9000 frames = 75 min) is an exact AC power flow of the IEEE 14-bus system (MATPOWER
case14 bus loads, generator set points and voltage set points) with the recording's true topology and
branch parameters.  Then all phasors of the frame are rotated by the accumulated frequency deviation,
the instrument errors and the polar noise are applied, the malfunctions are injected and the angles are
reported as principal values.

Load and generation model (rates are those of real bulk systems; see README, "Operating conditions"):
  * system load level s(t): steady, a ramp up (morning pickup) or down (evening decline) of 0.3 %/min,
    or a rise and fall at 0.3 %/min (a peak passing);
  * every load bus: an Ornstein-Uhlenbeck (OU) factor with 2 % standard deviation and tau = 300 s on
    P and Q (minute-scale ambient variation of the aggregate substation demand), a second OU factor of
    1 % / 300 s on Q only (power-factor variation), and a small fast component of 0.2 % / 0.5 s (the sum
    of many small independent switching loads);
  * bus 3 serves a steel works: an AC electric arc furnace of 30 MW / 24 MVAr whose arc makes the
    active power fluctuate by 8 % and the reactive power by 15 % of rating from frame to frame
    (OU, tau = 0.2 s) while the arc is on, with a power-off period of a few minutes for back-charging
    (10 to 15 min) in some recordings;
  * the generator at bus 2 follows the 5-minute real-time dispatch (a new set point every 300 s, reached
    at 8 MW/min, kept within 20 to 100 MW) and takes 25 % of the load deviation through AGC (first
    order, tau = 60 s);
    the slack generator at bus 1 takes the rest; buses 2, 3, 6 and 8 hold their voltage set points;
  * system frequency deviation: OU with tau = 30 s and the recording's standard deviation.
"""
import json
import sys

import numpy as np

import model as M

DT = 0.5                                  # 2 frames per second: 9000 frames = 75 minutes
FRAMES = 9000
PD = np.array([0, 21.7, 94.2, 47.8, 7.6, 11.2, 0, 0, 29.5, 9.0, 3.5, 6.1, 13.5, 14.9]) / 100
QD = np.array([0, 12.7, 19.0, -3.9, 1.6, 7.5, 0, 0, 16.6, 5.8, 1.8, 1.6, 5.8, 5.0]) / 100
VSET = {1: 1.06, 2: 1.045, 3: 1.01, 6: 1.07, 8: 1.09}
PG2 = 0.40
EAF_BUS, EAF_P, EAF_Q = 3, 0.30, 0.24          # part of the bus-3 load (94.2 MW in case14)

CONFIG = {
    # profile: (kind, start level, change over the recording); eaf_off: furnace power-off [first s, last s]
    "example": dict(seed=1001, profile=("steady", 1.00, 0.0), freq_std=0.0, eaf_off=None),
    "rec_a": dict(seed=1101, profile=("ramp", 0.775, 0.225), freq_std=0.02, eaf_off=(2600, 3500)),
    "rec_b": dict(seed=1201, profile=("swing", 0.8875, 0.1125), freq_std=0.03, eaf_off=None),
    "rec_c": dict(seed=1301, profile=("ramp", 1.0, -0.225), freq_std=0.0, eaf_off=(750, 1650)),
    "rec_d": dict(seed=1401, profile=("steady", 1.00, 0.0), freq_std=0.02, eaf_off=(3050, 3950)),
}

# malfunction kinds of the planted intervals (bus, first, last as in the truth file)
MALF = {
    "example": [dict(bus=10, kind="stale", ch=["V"])],
    "rec_a": [dict(bus=12, kind="timing", shape="tri", size=0.03), dict(bus=9, kind="timing", shape="step", size=3e-3),
              dict(bus=4, kind="stale", ch=["I"])],
    "rec_b": [dict(bus=3, kind="ratio", ch="V", size=8e-4), dict(bus=12, kind="timing", shape="ramp", size=0.02),
              dict(bus=10, kind="stale", ch=["V"])],
    "rec_c": [dict(bus=13, kind="timing", shape="step", size=4.5e-3)],
    "rec_d": [dict(bus=2, kind="stale", ch=["I"]), dict(bus=14, kind="timing", shape="step", size=2.4e-3),
              dict(bus=11, kind="ratio", ch="V", size=1.6e-3), dict(bus=6, kind="timing", shape="tri", size=0.025)],
}
EXAMPLE_INTERVALS = [dict(bus=10, first_frame=2400, last_frame=2439)]


def ou(rng, n, sigma, tau, size=()):
    """Stationary OU samples at spacing DT: x_{k+1} = a x_k + sigma sqrt(1 - a^2) w_k."""
    a = np.exp(-DT / tau)
    x = np.empty((n,) + size)
    x[0] = sigma * rng.standard_normal(size)
    w = rng.standard_normal((n,) + size) * sigma * np.sqrt(1 - a * a)
    for k in range(1, n):
        x[k] = a * x[k - 1] + w[k]
    return x


def profile(kind, s0, ds, t):
    if kind == "steady":
        return np.full_like(t, s0)
    if kind == "ramp":
        return s0 + ds * t / t[-1]
    return s0 + ds * (1 - np.abs(2 * t / t[-1] - 1))           # rise and fall


def power_flow(Y, P, Q, V0):
    """Newton-Raphson (polar) for one frame; P, Q net injections (pu); returns complex bus voltages."""
    pv = [b - 1 for b in (2, 3, 6, 8)]
    pq = [i for i in range(M.NB) if i not in pv and i != 0]
    nonslack = pv + pq
    V = V0.copy()
    for _ in range(30):
        S = V * np.conj(Y @ V)
        mis = np.r_[P[nonslack] - S.real[nonslack], Q[pq] - S.imag[pq]]
        if np.max(np.abs(mis)) < 1e-13:
            return V
        # Jacobian by the standard complex derivative formulas
        Ibus = Y @ V
        dS_dVa = 1j * np.diag(V) @ np.conj(np.diag(Ibus) - Y @ np.diag(V))
        dS_dVm = np.diag(V) @ np.conj(Y @ np.diag(V / np.abs(V))) + np.diag(V / np.abs(V)) @ np.conj(np.diag(Ibus))
        J = np.block([[dS_dVa.real[np.ix_(nonslack, nonslack)], dS_dVm.real[np.ix_(nonslack, pq)]],
                      [dS_dVa.imag[np.ix_(pq, nonslack)], dS_dVm.imag[np.ix_(pq, pq)]]])
        dx = np.linalg.solve(J, mis)
        va, vm = np.angle(V), np.abs(V)
        va[nonslack] += dx[:len(nonslack)]
        vm[pq] += dx[len(nonslack):]
        V = vm * np.exp(1j * va)
    raise RuntimeError("power flow did not converge")


def main():
    name, net_p, truth_p, out = sys.argv[1:5]
    net, tr = M.load_json(net_p), M.load_json(truth_p)
    cfg = CONFIG[name]
    rng = np.random.default_rng(cfg["seed"])
    Y = M.ybus(net["candidates"], tr["status"], tr["params"], net["shunt"])
    t = np.arange(FRAMES) * DT
    s = profile(*cfg["profile"], t)
    loadb = PD > 0
    nl = int(loadb.sum())
    a = ou(rng, FRAMES, 0.02, 300.0, (nl,))
    aq = ou(rng, FRAMES, 0.01, 300.0, (nl,))
    fast = ou(rng, FRAMES, 0.002, 0.5, (nl,))
    eaf_p = 1 + ou(rng, FRAMES, 0.08, 0.2)
    eaf_q = 1 + ou(rng, FRAMES, 0.15, 0.2)
    on = np.ones(FRAMES)
    if cfg["eaf_off"]:
        t0, t1 = cfg["eaf_off"]
        on[(t >= t0) & (t <= t1)] = 0.0
    Pl = np.zeros((FRAMES, M.NB)); Ql = np.zeros((FRAMES, M.NB))
    base_p = PD.copy(); base_q = QD.copy()
    base_p[EAF_BUS - 1] -= EAF_P; base_q[EAF_BUS - 1] -= EAF_Q
    Pl[:, loadb] = base_p[loadb] * s[:, None] * (1 + a + fast)
    Ql[:, loadb] = base_q[loadb] * s[:, None] * (1 + a + aq + fast)
    Pl[:, EAF_BUS - 1] += EAF_P * eaf_p * on
    Ql[:, EAF_BUS - 1] += EAF_Q * eaf_q * on
    # generation at bus 2: 5-min dispatch set points reached at 8 MW/min, plus 25 % AGC share of the
    # load deviation; set points stay within 20 to 100 MW
    rate = 0.08 / 60.0
    sp = np.full(FRAMES, PG2)
    prev = PG2
    for tk in 300.0 * np.arange(1, int(t[-1] // 300) + 1):
        new = float(np.clip(prev + rng.normal(0, 0.04), 0.20, 1.00))
        dur = max(abs(new - prev) / rate, DT)
        sp = np.where(t >= tk, prev + (new - prev) * np.clip((t - tk) / dur, 0, 1), sp)
        prev = new
    dev = Pl.sum(1) - Pl[0].sum()
    agc = np.zeros(FRAMES)
    al = np.exp(-DT / 60.0)
    for k in range(1, FRAMES):
        agc[k] = al * agc[k - 1] + (1 - al) * 0.25 * dev[k]
    Pg2 = sp + agc
    # frequency deviation and its accumulated angle
    df = ou(rng, FRAMES, cfg["freq_std"], 30.0) if cfg["freq_std"] > 0 else np.zeros(FRAMES)
    rot = 2 * np.pi * np.cumsum(df) * DT
    # power flows
    V = np.array([VSET.get(b, 1.0) for b in range(1, M.NB + 1)], complex)
    Vt = np.zeros((FRAMES, M.NB), complex)
    for k in range(FRAMES):
        P = -Pl[k].copy(); Q = -Ql[k].copy()
        P[1] += Pg2[k]
        V = power_flow(Y, P, Q, V)
        Vt[k] = V
    It = Vt @ Y.T
    assert np.max(np.abs(It[:, M.Z7])) < 1e-10
    rotf = np.exp(1j * rot)[:, None]
    Vt, It = Vt * rotf, It * rotf
    # instruments
    z = np.stack([(Vt if c == "V" else It)[:, b - 1] for c, b in M.CHANNELS], 1)
    cal = np.ones(len(M.CHANNELS), complex)
    for k, (c, b) in enumerate(M.CHANNELS):
        if f"{c}{b}" in tr["cal"]:
            e = tr["cal"][f"{c}{b}"]
            cal[k] = (1 + e["ratio_error"]) * np.exp(1j * e["phase_error"])
    zt = z * cal
    sm = np.array([net["noise"][c]["magnitude_rel_std"] for c, _ in M.CHANNELS])
    sa = np.array([net["noise"][c]["angle_std_rad"] for c, _ in M.CHANNELS])
    zm = zt * (1 + sm * rng.standard_normal(zt.shape)) * np.exp(1j * sa * rng.standard_normal(zt.shape))
    # malfunctions
    intervals = tr["corrupted"] if name != "example" else EXAMPLE_INTERVALS
    for iv, mf in zip(intervals, MALF[name]):
        assert iv["bus"] == mf["bus"]
        f0, f1 = iv["first_frame"], iv["last_frame"]
        L = f1 - f0 + 1
        ks = [k for k, (c, b) in enumerate(M.CHANNELS) if b == mf["bus"]]
        if mf["kind"] == "stale":
            for k in ks:
                if M.CHANNELS[k][0] in mf["ch"]:
                    zm[f0:f1 + 1, k] = zm[f0 - 1, k]
        elif mf["kind"] == "timing":
            i = np.arange(1, L + 1)
            g = {"step": np.ones(L), "ramp": i / L, "tri": 1 - np.abs(2 * i / (L + 1) - 1)}[mf["shape"]]
            zm[f0:f1 + 1, ks] *= np.exp(1j * mf["size"] * g)[:, None]
        else:
            k = [k for k in ks if M.CHANNELS[k][0] == mf["ch"]][0]
            zm[f0:f1 + 1, k] *= 1 + mf["size"]
    # write
    cols = ["frame"]
    for b in M.PMU_BUSES:
        cols += [f"V{b}_mag", f"V{b}_ang", f"I{b}_mag", f"I{b}_ang"]
    idx = {ch: k for k, ch in enumerate(M.CHANNELS)}
    with open(f"{out}/pmu.csv", "w") as fh:
        fh.write(",".join(cols) + "\n")
        for f in range(FRAMES):
            row = [str(f)]
            for b in M.PMU_BUSES:
                for c in ("V", "I"):
                    v = zm[f, idx[(c, b)]]
                    row += [repr(float(abs(v))), repr(float(np.angle(v)))]
            fh.write(",".join(row) + "\n")
    np.savez_compressed(f"{out}/true_states.npz", V=Vt, I=It, Pl=Pl, Ql=Ql, Pg2=Pg2, df=df, s=s)
    print(name, "written; load level %.3f..%.3f, Pg2 %.3f..%.3f, rotation span %.1f rad" % (
        s.min(), s.max(), Pg2.min(), Pg2.max(), rot.max() - rot.min()))


if __name__ == "__main__":
    main()
