"""Build probe programs that write crafted outputs from the planted truth (author-side only)."""
import json, os, sys
from decimal import Decimal
T = sys.argv[1]; out = sys.argv[2]
recs = {}
for r in sorted(os.listdir(f"{T}/tests/data")):
    d = f"{T}/tests/data/{r}"
    first = open(f"{d}/heldout_pmu.csv").read().split("\n")[1]
    recs[first] = dict(truth=json.load(open(f"{d}/heldout_truth.json")),
                       tol=json.load(open(f"{d}/heldout_tolerances.json")),
                       net=json.load(open(f"{d}/heldout_network.json")), name=r)
PROG = r'''
import json, sys
from decimal import Decimal, getcontext
getcontext().prec = 400
RECS = json.loads(__RECS__)
MODE = "__MODE__"
inp, outp = sys.argv[1], sys.argv[2]
first = open(inp + "/pmu.csv").read().split("\n")[1]
R = RECS[first]; tr, tol, net = R["truth"], R["tol"], R["net"]
special = {}
br = []
for c in net["candidates"]:
    e = {"id": c["id"], "in_service": bool(tr["status"][str(c["id"])])}
    if e["in_service"] and not c["known"]:
        for k, v in tr["params"][str(c["id"])].items():
            if k in ("r", "x", "tap"):
                e[k] = v
    br.append(e)
ch = {k: {"ratio_error": v["ratio_error"], "phase_error": v["phase_error"]} for k, v in tr["cal"].items()}
segs = [{"bus": s["bus"], "first_frame": s["first_frame"] - 30, "last_frame": s["last_frame"] + 30} for s in tr["corrupted"]]
# one branch value exactly at true + tolerance (exact decimal); on rec_a only in mode "over" it is pushed 1e-60 further
e = next(x for x in br if "x" in x and not net["candidates"][x["id"] - 1]["known"])
edge = Decimal(e["x"]) + Decimal(tol["branches"][str(e["id"])]["x"])
if MODE == "over" and R["name"] == "rec_a":
    edge += Decimal("1e-60")
e["x"] = "@@EDGE@@"
if MODE == "over" and R["name"] == "rec_b":
    segs[0]["first_frame"] -= 1          # 31 frames off
txt = json.dumps({"branches": br, "channels": ch, "corrupted_segments": segs}).replace('"@@EDGE@@"', str(edge))
open(outp, "w").write(txt)
'''
os.makedirs(out, exist_ok=True)
for mode in ("edge", "over"):
    open(f"{out}/probe_{mode}.py", "w").write(PROG.replace("__RECS__", repr(json.dumps(recs))).replace("__MODE__", mode))
print("ok", [v["name"] for v in recs.values()])
