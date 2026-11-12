"""score.py <rec_dir with truth.json, tolerances.json, network.json> <identification.json>"""
import json, sys
d, out = sys.argv[1], sys.argv[2]
tr, tol, net = (json.load(open(f"{d}/{n}")) for n in ("truth.json", "tolerances.json", "network.json"))
o = json.load(open(out))
cands = {c["id"]: c for c in net["candidates"]}
br = {b["id"]: b for b in o["branches"]}
topo = [cid for cid in cands if bool(br[cid]["in_service"]) != (cands[cid]["known"] or tr["status"][str(cid)])]
worst = (0, None)
for cid, ps in tr["params"].items():
    for k, v in ps.items():
        e = abs(br[int(cid)].get(k, float("nan")) - v) / (tol["branches"][cid][k] / 4)
        if not e <= worst[0]:
            worst = (e, f"branch {cid} {k}")
for ch, ps in tr["cal"].items():
    for k, v in ps.items():
        e = abs(o["channels"][ch][k] - v) / (tol["channels"][ch][k] / 4)
        if not e <= worst[0]:
            worst = (e, f"{ch} {k}")
segs = sorted((s["bus"], s["first_frame"], s["last_frame"]) for s in o["corrupted_segments"])
true = sorted((s["bus"], s["first_frame"], s["last_frame"]) for s in tr["corrupted"])
seg_ok = len(segs) == len(true) and all(a[0] == b[0] and abs(a[1] - b[1]) <= 30 and abs(a[2] - b[2]) <= 30 for a, b in zip(segs, true))
maxend = max([max(abs(a[1] - b[1]), abs(a[2] - b[2])) for a, b in zip(segs, true)] or [0]) if len(segs) == len(true) else None
ok = not topo and worst[0] <= 4 and seg_ok
print(f"{d.split('/')[-1]}: topology {'ok' if not topo else 'WRONG ' + str(topo)}, worst |err|/sigma {worst[0]:.2f} ({worst[1]}), "
      f"segments {'ok' if seg_ok else 'WRONG'} (max endpoint error {maxend}) -> {'PASS' if ok else 'FAIL'}")
if not seg_ok:
    print("   reported", segs, "true", true)
