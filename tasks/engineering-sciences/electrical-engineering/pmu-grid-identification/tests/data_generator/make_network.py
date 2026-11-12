"""Network description of a recording (network.json / heldout_network.json).

usage: python3 make_network.py <recording> <out.json>      recording: example, rec_a, rec_b, rec_c, rec_d

Candidates 1-20 are the branches of the IEEE 14-bus system (MATPOWER case14), with bus 7 as the explicit
star point of the three-winding transformer (4-7, 7-8 and 7-9; 4-7 and the two 7-x reactances are known);
candidates 21-24 are the extra corridors of the task (lines 1-4 and 3-5, cables 10-14 and 11-12), whose
nameplate values are typical of the neighbouring case14 branches.  Per recording only the operator's
database status of two candidates and the stated noise levels change.  The file is written exactly as
shipped (json.dumps with indent=1, no trailing newline).
"""
import json
import sys

CANDIDATES = [
    (1, 1, 2, "line", False, {"r": 0.01938, "x": 0.05917, "b": 0.0528}),
    (2, 1, 5, "line", False, {"r": 0.05403, "x": 0.22304, "b": 0.0492}),
    (3, 2, 3, "line", False, {"r": 0.04699, "x": 0.19797, "b": 0.0438}),
    (4, 2, 4, "line", False, {"r": 0.05811, "x": 0.17632, "b": 0.034}),
    (5, 2, 5, "line", False, {"r": 0.05695, "x": 0.17388, "b": 0.0346}),
    (6, 3, 4, "line", False, {"r": 0.06701, "x": 0.17103, "b": 0.0128}),
    (7, 4, 5, "series", False, {"r": 0.01335, "x": 0.04211}),
    (8, 4, 7, "transformer", True, {"x": 0.20912, "tap": 0.978}),
    (9, 4, 9, "transformer", False, {"x": 0.55618, "tap": 0.969}),
    (10, 5, 6, "transformer", False, {"x": 0.25202, "tap": 0.932}),
    (11, 6, 11, "series", False, {"r": 0.09498, "x": 0.1989}),
    (12, 6, 12, "series", False, {"r": 0.12291, "x": 0.25581}),
    (13, 6, 13, "series", False, {"r": 0.06615, "x": 0.13027}),
    (14, 7, 8, "reactance", True, {"x": 0.17615}),
    (15, 7, 9, "reactance", True, {"x": 0.11001}),
    (16, 9, 10, "series", False, {"r": 0.03181, "x": 0.0845}),
    (17, 9, 14, "series", False, {"r": 0.12711, "x": 0.27038}),
    (18, 10, 11, "series", False, {"r": 0.08205, "x": 0.19207}),
    (19, 12, 13, "series", False, {"r": 0.22092, "x": 0.19988}),
    (20, 13, 14, "series", False, {"r": 0.17093, "x": 0.34802}),
    (21, 1, 4, "line", False, {"r": 0.0482, "x": 0.1845, "b": 0.041}),
    (22, 3, 5, "line", False, {"r": 0.0612, "x": 0.1899, "b": 0.039}),
    (23, 10, 14, "series", False, {"r": 0.1025, "x": 0.2173}),
    (24, 11, 12, "series", False, {"r": 0.144, "x": 0.249}),
]
DATABASE_OPEN_DEFAULT = {21, 22, 23, 24}        # corridors 1-20 closed, the extra ones open

# per recording: candidates whose database status is flipped from the default, and the noise levels
RECORDINGS = {
    "example": dict(flip=set(), noise=((2e-4, 3e-4), (5e-4, 8e-4))),
    "rec_a": dict(flip=set(), noise=((1e-4, 5e-4), (2.5e-4, 1.25e-3))),
    "rec_b": dict(flip={13, 23}, noise=((1e-4, 8e-4), (3e-4, 2.4e-3))),
    "rec_c": dict(flip={13, 23}, noise=((1.5e-4, 7.5e-4), (4e-4, 2e-3))),
    "rec_d": dict(flip={13, 23}, noise=((2e-4, 4e-4), (5e-4, 1e-3))),
}


def network(name):
    cfg = RECORDINGS[name]
    (vm, va), (im, ia) = cfg["noise"]
    cands = []
    for cid, fb, tb, model, known, nominal in CANDIDATES:
        is_open = (cid in DATABASE_OPEN_DEFAULT) != (cid in cfg["flip"])
        cands.append({"id": cid, "from_bus": fb, "to_bus": tb, "model": model, "known": known,
                      "database_status": "open" if is_open else "closed", "nominal": dict(nominal)})
    return {"buses": list(range(1, 15)), "frame_rate": 2, "frames": 9000, "zero_injection_bus": 7,
            "reference_pmus": [1, 6],
            "noise": {"V": {"magnitude_rel_std": vm, "angle_std_rad": va},
                      "I": {"magnitude_rel_std": im, "angle_std_rad": ia}},
            "shunt": {"bus": 9, "nominal_b": 0.19},
            "candidates": cands}


if __name__ == "__main__":
    with open(sys.argv[2], "w") as fh:
        fh.write(json.dumps(network(sys.argv[1]), indent=1))
