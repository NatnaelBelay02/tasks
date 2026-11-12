"""Planted truth of a recording (the held-out heldout_truth.json; the example's truth is not shipped).

usage: python3 make_truth.py <recording> <network.json> <out.json>

Rules (the same for every recording, stated in the README): line and cable r, x drawn uniformly within
+-10 % of nameplate; transformer x within +-5 %, tap nameplate + k 0.625 % with k uniform in -3..3; ratio
errors uniform in +-1 %, phase errors uniform in +-10.5 mrad; the reference PMUs (buses 1 and 6) exact.
The draws use numpy's default_rng with the seed below, in candidate order and then channel order.
For the four held-out recordings (drawn in version 2 of the task) every unknown candidate gets its draws
and only those of the candidates in service are stored; for the example (drawn in version 3) only the
candidates in service are drawn.  Both rules are kept as they were used, so that the shipped files are
reproduced exactly (check with reproduce.py).
The topology (in-service status) is the database status except for the listed candidates, and the
malfunction intervals are those injected by gen_recording.py.  Written exactly as shipped (json.dumps,
indent=1, no trailing newline).
"""
import json
import sys

import numpy as np

import model as M

RECORDINGS = {
    # seed; in-service status that differs from the database; malfunction intervals (bus, first, last frame)
    "example": dict(seed=907, status={"13": False, "23": True}, corrupted=[(10, 2400, 2439)], draw_all=False),
    "rec_a": dict(seed=111, draw_all=True, status={"13": False, "23": True},
                  corrupted=[(12, 3100, 3499), (9, 5200, 5899), (4, 7000, 7059)]),
    "rec_b": dict(seed=121, draw_all=True, status={"18": False, "21": True},
                  corrupted=[(3, 800, 1399), (12, 4500, 4799), (10, 6000, 6039)]),
    "rec_c": dict(seed=131, draw_all=True, status={}, corrupted=[(13, 3000, 4199)]),
    "rec_d": dict(seed=141, draw_all=True, status={"22": True, "24": True},
                  corrupted=[(2, 1000, 1089), (14, 2500, 3299), (11, 4000, 4599), (6, 6500, 6899)]),
}


def truth(name, net):
    cfg = RECORDINGS[name]
    rng = np.random.default_rng(cfg["seed"])
    status = {str(c["id"]): (c["known"] or c["database_status"] == "closed") for c in net["candidates"]}
    status.update(cfg["status"])
    params = {}
    for c in net["candidates"]:
        in_service = status[str(c["id"])]
        if c["known"] or not (in_service or cfg["draw_all"]):
            continue
        n = c["nominal"]
        if c["model"] in ("line", "series"):
            v = {"r": n["r"] * (1 + rng.uniform(-0.1, 0.1)), "x": n["x"] * (1 + rng.uniform(-0.1, 0.1))}
        else:                                       # transformer (the reactances are all known)
            v = {"x": n["x"] * (1 + rng.uniform(-0.05, 0.05)), "tap": round(n["tap"] + 0.00625 * int(rng.integers(-3, 4)), 5)}
        if in_service:
            params[str(c["id"])] = v
    cal = {}
    for c, b in M.CHANNELS:
        if b not in M.REF:
            cal[f"{c}{b}"] = {"ratio_error": rng.uniform(-0.01, 0.01), "phase_error": rng.uniform(-0.0105, 0.0105)}
    return {"status": status, "params": params, "cal": cal,
            "corrupted": [{"bus": b, "first_frame": f, "last_frame": l} for b, f, l in cfg["corrupted"]]}


if __name__ == "__main__":
    with open(sys.argv[3], "w") as fh:
        fh.write(json.dumps(truth(sys.argv[1], M.load_json(sys.argv[2])), indent=1))
