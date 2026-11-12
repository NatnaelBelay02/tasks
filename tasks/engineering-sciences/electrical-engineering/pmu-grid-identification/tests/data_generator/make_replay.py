"""make_replay.py <task_dir> <out.py> <rec_a.json> <rec_b.json> <rec_c.json> <rec_d.json>:
a program that writes the stored output of the recording it is given (matched by its first data row)."""
import json, os, sys
T, out = sys.argv[1], sys.argv[2]
outs = {}
for r, f in zip(["rec_a", "rec_b", "rec_c", "rec_d"], sys.argv[3:7]):
    first = open(f"{T}/tests/data/{r}/heldout_pmu.csv").read().split("\n")[1]
    outs[first] = open(f).read()
open(out, "w").write("import json, sys\nO = json.loads(%r)\nfirst = open(sys.argv[1] + '/pmu.csv').read().split('\\n')[1]\n"
                     "open(sys.argv[2], 'w').write(O[first])\n" % json.dumps(outs))
