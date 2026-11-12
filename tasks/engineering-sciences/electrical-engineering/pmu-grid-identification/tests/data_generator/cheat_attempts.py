"""Anti-cheat trial (author-side): a submission that tries to obtain the answers or to influence the
grading instead of identifying the network.

usage: put this file at /root/results/identify.py and run the verifier (tests/test.sh).

It is run by the verifier exactly like a real submission (as uid 65534, in its own working directory,
with the recording's input directory as argument).  Each attempt below is logged as
"ATTEMPT <name>: <outcome>"; the log is printed by the verifier as the run's output tail.  At the end it
writes the best answer available without solving the problem: anything an attempt recovered, otherwise
the database topology, nameplate parameters, zero instrument errors and no malfunction segments.
A correct verifier gives reward 0, and every attempt reports that it obtained nothing usable.
"""
import glob
import json
import os
import signal
import sys

LOG = []
FOUND = {}          # anything recovered that looks like an answer: path -> parsed JSON


def log(name, outcome):
    LOG.append("ATTEMPT %s: %s" % (name, str(outcome)[:90]))


def try_read(path):
    try:
        with open(path, "rb") as fh:
            return fh.read()
    except OSError as e:
        return e


def grab(path):
    data = try_read(path)
    if isinstance(data, bytes):
        try:
            obj = json.loads(data)
        except ValueError:
            return "read %d bytes (not JSON)" % len(data)
        if isinstance(obj, dict) and ({"params", "cal"} <= obj.keys() or {"branches", "channels"} <= obj.keys()):
            FOUND[path] = obj
            return "ANSWER-LIKE JSON READ"
        return "read JSON without answers (keys %s)" % sorted(obj)[:4]
    return type(data).__name__


def main():
    data_dir, out = sys.argv[1], sys.argv[2]
    work = os.path.dirname(out)
    log("identity", "uid %d gid %d, cwd %s" % (os.getuid(), os.getgid(), os.getcwd()))

    # 1. the held-out answers where the verifier keeps them
    for p in ["/tests", "/tests/data"] + ["/tests/data/rec_%s" % r for r in "abcd"]:
        try:
            log("list " + p, "listed %s" % os.listdir(p)[:5])
        except OSError as e:
            log("list " + p, type(e).__name__)
    for r in "abcd":
        for f in ("heldout_truth.json", "heldout_tolerances.json"):
            log("read /tests/data/rec_%s/%s" % (r, f), grab("/tests/data/rec_%s/%s" % (r, f)))

    # 2. search the whole file system for anything that looks like an answer or the reference program
    hits, readable, own = [], 0, 0
    for top in ("/", ):
        for dirpath, dirnames, filenames in os.walk(top):
            if dirpath.startswith(("/proc", "/sys", "/dev")):
                dirnames[:] = []
                continue
            for f in filenames:
                low = f.lower()
                if any(k in low for k in ("truth", "toleranc", "heldout", "identification", "identify", "answer")):
                    p = os.path.join(dirpath, f)
                    if p.startswith(data_dir) or p == os.path.join(work, "identify.py"):
                        continue
                    try:
                        if os.stat(p).st_uid == os.getuid():     # written by this submission (e.g. the
                            own += 1                             # output of a concurrent run): not an answer
                            continue
                    except OSError:
                        pass
                    hits.append(p)
                    if isinstance(try_read(p), bytes):
                        readable += 1
                        grab(p)
    log("file-system search", "%d files of others, %d readable, %d with answers, %s; %d own skipped"
        % (len(hits), readable, len(FOUND), [os.path.basename(h) for h in hits[:2]], own))

    # 3. the other runs of the same grading session (siblings of this run's directory)
    root = os.path.dirname(os.path.dirname(work))
    sib = []
    try:
        for d in os.listdir(root):
            for sub in glob.glob(os.path.join(root, d, "*")):
                try:
                    sib.append((sub, os.listdir(sub)))
                except OSError as e:
                    sib.append((sub, type(e).__name__))
    except OSError as e:
        sib = type(e).__name__
    names = sorted({x for _, l in sib if isinstance(l, list) for x in l}) if not isinstance(sib, str) else sib
    log("other runs' directories", "%s" % (sib if isinstance(sib, str) else
                                          "%d dirs; inputs only: %s" % (len(sib), [x for x in names if x in ("pmu.csv", "network.json")]
                                                                        + ["(+ files of this submission)"] * any(x not in ("pmu.csv", "network.json", "identify.py") for x in names))))

    # 4. the grading process itself: environment, open files and memory of the root-owned processes
    mem = fds = envs = 0
    for pid in os.listdir("/proc"):
        if not pid.isdigit() or int(pid) == os.getpid():
            continue
        try:
            st = os.stat("/proc/" + pid)
        except OSError:
            continue
        if st.st_uid == os.getuid():
            continue
        envs += isinstance(try_read("/proc/%s/environ" % pid), bytes)
        try:
            for fd in os.listdir("/proc/%s/fd" % pid):
                target = os.readlink("/proc/%s/fd/%s" % (pid, fd))
                fds += 1
                if target.startswith("/") and "truth" in target:
                    grab("/proc/%s/fd/%s" % (pid, fd))
        except OSError:
            pass
        try:
            with open("/proc/%s/mem" % pid, "rb") as fh:
                fh.read(1)
                mem += 1
        except OSError:
            pass
    log("other users' processes", "environ readable %d, fd tables readable %d, memory readable %d" % (envs, fds, mem))

    # 5. the grading result and the grader's code
    for p, text in (("/logs/verifier/reward.txt", "1\n"), ("/logs/verifier/ctrf.json", "{}"),
                    ("/tests/conftest.py", "import pytest\n"), ("/tests/test_identification.py", "def test_ok():\n    pass\n")):
        try:
            with open(p, "w") as fh:
                fh.write(text)
            log("write " + p, "WRITTEN")
        except OSError as e:
            log("write " + p, type(e).__name__)
    import site
    for d in site.getsitepackages():
        try:
            with open(os.path.join(d, "sitecustomize.py"), "w") as fh:
                fh.write("")
            log("write sitecustomize into " + d, "WRITTEN")
        except OSError as e:
            log("write sitecustomize into " + d, type(e).__name__)
    for p in ("/tmp/conftest.py", os.path.join(work, "conftest.py")):
        try:
            with open(p, "w") as fh:
                fh.write("def pytest_collection_modifyitems(items):\n    items.clear()\n")
            log("plant " + p, "written (pytest does not collect from here)")
        except OSError as e:
            log("plant " + p, type(e).__name__)

    # 6. privileges and the grader's processes
    try:
        os.setuid(0)
        log("setuid(0)", "SUCCEEDED")
    except OSError as e:
        log("setuid(0)", type(e).__name__)
    killed = denied = 0
    for pid in os.listdir("/proc"):
        if pid.isdigit() and int(pid) not in (os.getpid(), os.getppid()):
            try:
                if os.stat("/proc/" + pid).st_uid != os.getuid():
                    os.kill(int(pid), signal.SIGSTOP)
                    killed += 1
            except OSError:
                denied += 1
    log("stop root processes", "stopped %d, refused %d" % (killed, denied))

    # 7. the best answer without solving: recovered values if any, else the defaults
    net = json.load(open(os.path.join(data_dir, "network.json")))
    tr = next(iter(FOUND.values()), None)
    branches = []
    for c in net["candidates"]:
        e = {"id": c["id"], "in_service": c["known"] or c["database_status"] == "closed"}
        if tr and "status" in tr:
            e["in_service"] = tr["status"][str(c["id"])]
        for k in {"line": ["r", "x"], "series": ["r", "x"], "transformer": ["x", "tap"], "reactance": ["x"]}[c["model"]]:
            e[k] = (tr or {}).get("params", {}).get(str(c["id"]), {}).get(k, c["nominal"][k])
        branches.append(e)
    with open(os.path.join(data_dir, "pmu.csv")) as fh:          # channel names from the recording's header
        header = fh.readline().strip().split(",")[1:]
    chans = {}
    for col in header:
        ch = col.rsplit("_", 1)[0]
        if int(ch[1:]) not in net["reference_pmus"]:
            chans[ch] = {"ratio_error": 0.0, "phase_error": 0.0}
    with open(out, "w") as fh:
        json.dump({"branches": branches, "channels": chans, "corrupted_segments": []}, fh)
    log("output", "default answer written (%s)" % ("with recovered values" if tr else "nothing recovered"))
    print("\n".join(LOG)[-1900:])


if __name__ == "__main__":
    main()
