"""Author-side network and measurement model for pmu-grid-identification (shares no code with
solution/identify.py).  Branch admittances follow MATPOWER makeYbus with the off-nominal tap on the
from side; bus 7 is the zero-injection star point; per-unit on 100 MVA."""
import json
import numpy as np

NB = 14
Z7 = 6                                   # 0-based index of bus 7
REF = (1, 6)
PMU_BUSES = [b for b in range(1, NB + 1) if b != 7]
CHANNELS = [("V", b) for b in PMU_BUSES] + [("I", b) for b in PMU_BUSES]     # 26 channels


def load_json(p):
    with open(p) as fh:
        return json.load(fh)


def ybus(cands, status, params, shunt):
    """Ybus (14 x 14) for the in-service candidates; params[id] holds the true unknown values."""
    Y = np.zeros((NB, NB), complex)
    for c in cands:
        cid = c["id"]
        if not (c["known"] or status[str(cid)]):
            continue
        p = dict(c["nominal"])
        if not c["known"] and str(cid) in params:
            p.update(params[str(cid)])
        r = p.get("r", 0.0)
        ys = 1.0 / complex(r, p["x"])
        tap = p.get("tap", 1.0) if c["model"] == "transformer" else 1.0
        bc = c["nominal"].get("b", 0.0) if c["model"] == "line" else 0.0
        i, j = c["from_bus"] - 1, c["to_bus"] - 1
        Y[i, i] += (ys + 0.5j * bc) / tap ** 2
        Y[j, j] += ys + 0.5j * bc
        Y[i, j] -= ys / tap
        Y[j, i] -= ys / tap
    Y[shunt["bus"] - 1, shunt["bus"] - 1] += 1j * shunt["nominal_b"]
    return Y


def read_pmu(path):
    raw = np.genfromtxt(path, delimiter=",", names=True)
    Z = np.stack([raw[f"{c}{b}_mag"] * np.exp(1j * raw[f"{c}{b}_ang"]) for c, b in CHANNELS], 1)
    return Z
