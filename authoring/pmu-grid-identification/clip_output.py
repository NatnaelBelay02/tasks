"""clip_output.py <network.json> <in.json> <out.json>: the output clipping of the programs, applied offline
to outputs written before it was added (the clipping is a pure post-processing step)."""
import json, sys
net, o = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
C = {c["id"]: c for c in net["candidates"]}
for e in o["branches"]:
    c = C[e["id"]]
    for k in ("r", "x", "tap"):
        if k in e:
            h = 0.019 if k == "tap" else (0.05 if c["model"] == "transformer" else 0.10) * abs(c["nominal"][k])
            e[k] = float(min(max(e[k], c["nominal"][k] - h), c["nominal"][k] + h))
for v in o["channels"].values():
    v["ratio_error"] = float(min(max(v["ratio_error"], -0.01), 0.01))
    v["phase_error"] = float(min(max(v["phase_error"], -0.0105), 0.0105))
json.dump(o, open(sys.argv[3], "w"), indent=1)
