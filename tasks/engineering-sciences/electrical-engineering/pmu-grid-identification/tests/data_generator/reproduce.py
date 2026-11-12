"""Rebuild every data file of the task from the scripts and seeds in this directory, and compare
with the shipped files.

usage: python3 reproduce.py [work_dir]          (run from any directory; about 5 minutes on 5 cores)

Pipeline per recording (example, rec_a ... rec_d):
  make_network.py   network description            -> environment/data/network.json, tests/data/rec_*/heldout_network.json
  make_truth.py     planted truth (seeds 907, 111, 121, 131, 141)
                                                   -> tests/data/rec_*/heldout_truth.json (the example's is not shipped)
  gen_recording.py  the PMU recording (seeds 1001, 1101, 1201, 1301, 1401) and the true bus voltages
                                                   -> environment/data/pmu.csv, tests/data/rec_*/heldout_pmu.csv
  fisher.py         Cramer-Rao bounds, tolerances = 4 sigma
                                                   -> environment/data/tolerances.json, tests/data/rec_*/heldout_tolerances.json

Networks, truths and recordings must be byte-identical (SHA-256).  The tolerances are the diagonal of the
inverse of a Fisher matrix with a condition number up to 6e8; OpenBLAS picks its kernels by CPU, so their
last digits depend on the machine.  They are compared numerically and must agree to 1e-6 relative (the
shipped files and a rebuild agree to about 2e-9; the tolerance test compares errors with 4 sigma, so
this changes no verdict).  Exit status 0 only if everything agrees.
"""
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
TASK = os.path.dirname(os.path.dirname(HERE))
RECS = ["example", "rec_a", "rec_b", "rec_c", "rec_d"]


def shipped(rec, kind):
    if rec == "example":
        return None if kind == "truth" else os.path.join(TASK, "environment", "data", {
            "network": "network.json", "pmu": "pmu.csv", "tolerances": "tolerances.json"}[kind])
    return os.path.join(TASK, "tests", "data", rec, {"network": "heldout_network.json", "truth": "heldout_truth.json",
                                                    "pmu": "heldout_pmu.csv", "tolerances": "heldout_tolerances.json"}[kind])


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def run(args):
    subprocess.run([sys.executable] + args, cwd=HERE, check=True, stdout=subprocess.DEVNULL)


def build(rec, work):
    d = os.path.join(work, rec)
    os.makedirs(d, exist_ok=True)
    run(["make_network.py", rec, f"{d}/network.json"])
    run(["make_truth.py", rec, f"{d}/network.json", f"{d}/truth.json"])
    run(["gen_recording.py", rec, f"{d}/network.json", f"{d}/truth.json", d])
    run(["fisher.py", d, f"{d}/tolerances.json"])
    return rec


def leaves(x, pre=""):
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            out.update(leaves(v, f"{pre}/{k}"))
        return out
    return {pre: x}


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="pmu_reproduce_")
    with ThreadPoolExecutor(len(RECS)) as ex:
        list(ex.map(lambda r: build(r, work), RECS))
    ok = True
    for rec in RECS:
        d = os.path.join(work, rec)
        for kind, name in (("network", "network.json"), ("truth", "truth.json"), ("pmu", "pmu.csv")):
            ref = shipped(rec, kind)
            if ref is None:
                print(f"{rec:8s} {kind:10s} {sha(f'{d}/{name}')[:16]}  (not shipped)")
                continue
            same = sha(f"{d}/{name}") == sha(ref)
            ok &= same
            print(f"{rec:8s} {kind:10s} {sha(ref)[:16]}  {'identical' if same else 'DIFFERS'}")
        a = leaves(json.load(open(f"{d}/tolerances.json")))
        b = leaves(json.load(open(shipped(rec, "tolerances"))))
        same_keys = a.keys() == b.keys()
        rel = max(abs(a[k] / b[k] - 1) for k in b) if same_keys else float("inf")
        good = same_keys and rel <= 1e-6
        ok &= good
        print(f"{rec:8s} tolerances {len(b)} values, largest relative difference {rel:.1e}  {'agree' if good else 'DIFFER'}")
    print("all data files reproduced" if ok else "REPRODUCTION FAILED", f"(work directory {work})")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
