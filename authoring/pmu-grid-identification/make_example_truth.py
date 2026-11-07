"""Planted truth of the example recording (the package ships its tolerances, not its truth).
Same rules as the held-out recordings: line and cable r, x within +-10 % of nameplate, transformer x
within +-5 %, taps nameplate +- up to three 0.625 % steps, ratio errors uniform in +-1 %, phase errors
uniform in +-10.5 mrad, reference PMUs 1 and 6 exact.  Topology: 6-13 out (database says closed),
10-14 in (database says open), the other database-open corridors 1-4, 3-5, 11-12 out."""
import json, sys
import numpy as np
import model as M
net = M.load_json(sys.argv[1])
rng = np.random.default_rng(907)
status = {str(c["id"]): (c["known"] or c["database_status"] == "closed") for c in net["candidates"]}
status["13"] = False; status["23"] = True
params = {}
for c in net["candidates"]:
    if c["known"] or not status[str(c["id"])]:
        continue
    n = c["nominal"]
    if c["model"] in ("line", "series"):
        params[str(c["id"])] = {"r": n["r"] * (1 + rng.uniform(-0.1, 0.1)), "x": n["x"] * (1 + rng.uniform(-0.1, 0.1))}
    elif c["model"] == "transformer":
        params[str(c["id"])] = {"x": n["x"] * (1 + rng.uniform(-0.05, 0.05)), "tap": round(n["tap"] + 0.00625 * int(rng.integers(-3, 4)), 5)}
cal = {}
for c, b in M.CHANNELS:
    if b not in M.REF:
        cal[f"{c}{b}"] = {"ratio_error": rng.uniform(-0.01, 0.01), "phase_error": rng.uniform(-0.0105, 0.0105)}
json.dump({"status": status, "params": params, "cal": cal,
           "corrupted": [{"bus": 10, "first_frame": 2400, "last_frame": 2439}]}, open(sys.argv[2], "w"), indent=1)
print({k: v for k, v in status.items() if not v})
