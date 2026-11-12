"""Verifier for the PMU grid-identification task.

The submitted program /root/results/identify.py is run once on each of four held-out recordings, all four
concurrently, as an unprivileged user (uid/gid 65534) in its own writable working directory, with the
recording in a read-only input directory.  Each output is then graded against the planted truth of that
recording.  The tests map to the numbered requirements of instruction.md and every test must pass for
reward 1.  The reference values are the planted truth used to generate each recording; the tolerances are
4 x the Cramer-Rao standard deviation of each quantity (heldout_tolerances.json of each recording)."""
import decimal
import hashlib
import itertools
import json
import os
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
from decimal import Decimal

import pytest

PROGRAM = "/root/results/identify.py"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
RECORDINGS = sorted(os.listdir(DATA))
WALL_SEC = 5400                 # per run; all runs are concurrent
MAX_OUTPUT_BYTES = 1_000_000
SANDBOX_UID = SANDBOX_GID = 65534
SANDBOX_ENV = {"HOME": "/tmp", "PATH": "/usr/local/bin:/usr/bin:/bin", "LANG": "C.UTF-8",
               "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
               "NUMEXPR_NUM_THREADS": "1"}
PARAM_KEYS = {"line": ["r", "x"], "series": ["r", "x"], "transformer": ["x", "tap"], "reactance": ["x"]}
SEG_TOL = 30

# Exact decimal arithmetic: sums of two doubles need at most ~1100 significant digits; any inexact
# operation raises instead of rounding.
EXACT = decimal.Context(prec=4000, Emax=decimal.MAX_EMAX, Emin=decimal.MIN_EMIN, traps=[decimal.Inexact])


class Recording:
    def __init__(self, name):
        d = os.path.join(DATA, name)
        self.name = name
        self.truth = json.load(open(os.path.join(d, "heldout_truth.json")))
        self.tol = json.load(open(os.path.join(d, "heldout_tolerances.json")))
        self.net = json.load(open(os.path.join(d, "heldout_network.json")))
        self.cands = {c["id"]: c for c in self.net["candidates"]}
        self.channels = sorted(self.tol["channels"])
        self.last_frame = self.net["frames"] - 1


REC = {name: Recording(name) for name in RECORDINGS}


# ---------------------------------------------------------------- sandboxed execution
def _demote():
    os.setsid()
    os.setgroups([])
    os.setgid(SANDBOX_GID)
    os.setuid(SANDBOX_UID)


def _sandbox_pids():
    pids = []
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            with open(f"/proc/{pid}/status") as fh:
                uid = next(line for line in fh if line.startswith("Uid:")).split()[1]
        except (OSError, StopIteration):
            continue
        if int(uid) == SANDBOX_UID:
            pids.append(int(pid))
    return pids


def _kill_sandbox_user():
    """Kill every process owned by the sandbox uid (background children included).  A helper process
    that has become the sandbox user calls kill(-1, SIGKILL), which signals all of that user's processes
    in one call, so processes that keep forking cannot escape; repeated until none is left."""
    for _ in range(200):
        if not _sandbox_pids():
            return
        pid = os.fork()
        if pid == 0:
            try:
                os.setgroups([])
                os.setgid(SANDBOX_GID)
                os.setuid(SANDBOX_UID)
                os.kill(-1, signal.SIGKILL)
            finally:
                os._exit(0)
        os.waitpid(pid, 0)
        time.sleep(0.05)
    raise RuntimeError("processes of the sandbox user survived repeated kill(-1, SIGKILL)")


def _read_regular(path, limit):
    """Read a regular file without following a symlink; None if absent, not regular, or too large.
    O_NONBLOCK makes the open return at once for a FIFO or device (rejected by the type check), and at
    most limit + 1 bytes are read, so the reader cannot block or grow without bound."""
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC)
    except OSError:
        return None
    try:
        st = os.fstat(fd)
        if not stat.S_ISREG(st.st_mode) or st.st_size > limit:
            return None
        chunks, size = [], 0
        while size <= limit:
            b = os.read(fd, min(1 << 20, limit + 1 - size))
            if not b:
                break
            chunks.append(b)
            size += len(b)
        return None if size > limit else b"".join(chunks)
    finally:
        os.close(fd)


def _stage(program, name, root):
    """Fresh run directory: read-only input directory with the recording, writable working directory
    holding only a copy of the program."""
    base = tempfile.mkdtemp(prefix="run_", dir=root)
    os.chmod(base, 0o755)
    # a name the program cannot know in advance: derived from the private temporary directory
    data_dir = os.path.join(base, "input_" + hashlib.sha256((base + name).encode()).hexdigest()[:8])
    os.mkdir(data_dir)
    for f in ("pmu.csv", "network.json"):             # stored as heldout_<name>, staged under the contract name
        shutil.copyfile(os.path.join(DATA, name, "heldout_" + f), os.path.join(data_dir, f))
        os.chmod(os.path.join(data_dir, f), 0o444)
    os.chmod(data_dir, 0o555)
    work = os.path.join(base, "work")
    os.mkdir(work)
    os.chmod(work, 0o1777)
    with open(os.path.join(work, "identify.py"), "wb") as fh:
        fh.write(program)
    os.chmod(os.path.join(work, "identify.py"), 0o444)
    return data_dir, work


def _check_sandbox(data_dir, work):
    """The stated execution contract holds for the sandbox user: the working directory and /tmp are
    writable, the input directory is not, and /tests is unreadable."""
    probe = ("import os,sys\n"
             "open('probe.txt','w').write('x'); os.remove('probe.txt')\n"
             "p='/tmp/probe_%d' % os.getpid(); open(p,'w').write('x'); os.remove(p)\n"
             "d=sys.argv[1]\n"
             "assert open(os.path.join(d,'network.json')).read(1)\n"
             "try:\n    open(os.path.join(d,'w.txt'),'w'); sys.exit('input directory writable')\n"
             "except PermissionError:\n    pass\n"
             "try:\n    os.listdir('/tests'); sys.exit('/tests readable')\n"
             "except PermissionError:\n    pass\n")
    r = subprocess.run([sys.executable, "-c", probe, data_dir], cwd=work, env=SANDBOX_ENV, preexec_fn=_demote,
                       stdin=subprocess.DEVNULL, capture_output=True, timeout=60)
    assert r.returncode == 0, f"sandbox self-check failed: {r.stderr.decode(errors='replace')[-500:]}"


@pytest.fixture(scope="session")
def runs():
    """Run the submitted program once on every held-out recording (concurrently) and collect the outputs."""
    program = _read_regular(PROGRAM, 50_000_000)
    if program is None:
        return {name: {"error": f"{PROGRAM} is missing, not a regular file, or a symlink"} for name in RECORDINGS}
    root = tempfile.mkdtemp(prefix="pgi_")
    os.chmod(root, 0o755)
    procs = {}
    for name in RECORDINGS:                  # fixed (sorted) launch order: the verdict does not depend on chance
        data_dir, work = _stage(program, name, root)
        _check_sandbox(data_dir, work)
        out = os.path.join(work, "identification.json")
        log = tempfile.TemporaryFile()                  # unlinked file, not a pipe
        p = subprocess.Popen([sys.executable, "identify.py", data_dir, out], cwd=work, env=SANDBOX_ENV,
                             preexec_fn=_demote, stdin=subprocess.DEVNULL, stdout=log, stderr=log)
        procs[name] = (p, out, log, time.monotonic())
    results = {}
    for name, (p, out, log, t0) in procs.items():
        try:
            rc = p.wait(timeout=max(1.0, t0 + WALL_SEC - time.monotonic()))
            elapsed = time.monotonic() - t0
            results[name] = {"returncode": rc, "elapsed": elapsed}
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGKILL)
            p.wait()
            results[name] = {"error": f"program did not finish within {WALL_SEC} s"}
    _kill_sandbox_user()
    for name, (p, out, log, t0) in procs.items():
        log.seek(0, os.SEEK_END)
        log.seek(max(0, log.tell() - 2000))
        tail = log.read().decode(errors="replace")
        log.close()
        r = results[name]
        print(f"{name}: {r}; output tail:\n{tail}")
        if "error" in r:
            continue
        if r["returncode"] != 0:
            r["error"] = f"program exited with status {r['returncode']}; last output:\n{tail[-800:]}"
            continue
        text = _read_regular(out, MAX_OUTPUT_BYTES)
        if text is None:
            r["error"] = "identification.json missing, not a regular file, a symlink, or larger than 1 MB"
            continue
        r["text"] = text
    shutil.rmtree(root, ignore_errors=True)
    return results


# ---------------------------------------------------------------- strict JSON loading
class _BadJSON(ValueError):
    pass


_MAX_DOUBLE = Decimal(sys.float_info.max)      # (2 - 2**-52) * 2**1023, exact


def _check_literal(text):
    """Every number literal must not exceed the largest finite double in absolute value (exact comparison)."""
    if Decimal(text).copy_abs() > _MAX_DOUBLE:   # copy_abs is exact; abs() would round
        raise _BadJSON(f"number {text[:30]} exceeds the largest finite double in absolute value")


def _float(text):
    """Non-integer literals are kept as exact Decimals, so tolerance checks use the value as written."""
    _check_literal(text)
    return Decimal(text)


def _int(text):
    _check_literal(text)
    return int(text)


def _constant(token):
    raise _BadJSON(f"non-finite constant {token}")


def _unique(pairs):
    obj = {}
    for k, v in pairs:
        if k in obj:
            raise _BadJSON(f"object repeats key {k!r}")
        obj[k] = v
    return obj


def _load(runs, name):
    r = runs[name]
    if "error" in r:
        pytest.fail(f"{name}: {r['error']}")
    problem = None
    try:
        data = json.loads(r["text"].decode("utf-8"), parse_constant=_constant, parse_float=_float,
                          parse_int=_int, object_pairs_hook=_unique)
    except (ValueError, UnicodeDecodeError, RecursionError) as exc:
        problem = f"{type(exc).__name__}: {str(exc)[:200]}"
    if problem is not None:              # reported outside the except block: no chained traceback
        pytest.fail(f"{name}: identification.json is not strict UTF-8 JSON with in-range numbers and unique keys: {problem}")
    assert isinstance(data, dict), f"{name}: top level of identification.json must be an object"
    return data


def _is_int(v):
    """A JSON integer literal (no decimal point or exponent): decoded as int, never bool."""
    return isinstance(v, int) and not isinstance(v, bool)


def _is_num(v):
    """Any JSON number: int (integer literal) or Decimal (other literals)."""
    return _is_int(v) or isinstance(v, Decimal)


def _within(got, true, tol):
    """|exact(got) - true| <= tol, with true and tol the exact values of the stored doubles."""
    t, d = Decimal(true), Decimal(tol)
    lo, hi = EXACT.subtract(t, d), EXACT.add(t, d)
    return lo <= Decimal(got) <= hi


def _branches(rec, data):
    br = data.get("branches")
    assert isinstance(br, list), f"{rec.name}: branches must be a list"
    out = {}
    for j, b in enumerate(br):
        assert isinstance(b, dict), f"{rec.name}: branches[{j}] must be an object"
        assert _is_int(b.get("id")) and b["id"] in rec.cands, \
            f"{rec.name}: branches[{j}].id must be a candidate id written as a JSON integer"
        assert b["id"] not in out, f"{rec.name}: candidate {b['id']} appears more than once in branches"
        assert isinstance(b.get("in_service"), bool), f"{rec.name}: branches[{j}].in_service must be true or false"
        out[b["id"]] = b
    assert sorted(out) == sorted(rec.cands), \
        f"{rec.name}: branches must have exactly one entry per candidate id {sorted(rec.cands)}"
    return out


def _channels(rec, data):
    ch = data.get("channels")
    assert isinstance(ch, dict), f"{rec.name}: channels must be an object"
    assert sorted(ch) == rec.channels, f"{rec.name}: channels keys {sorted(ch)} must be exactly {rec.channels}"
    for k, v in ch.items():
        assert isinstance(v, dict), f"{rec.name}: channels[{k}] must be an object"
        for key in ("ratio_error", "phase_error"):
            assert _is_num(v.get(key)), f"{rec.name}: channels[{k}].{key} must be a JSON number"
    return ch


def _segments(rec, data):
    segs = data.get("corrupted_segments")
    assert isinstance(segs, list), f"{rec.name}: corrupted_segments must be a list"
    for j, s in enumerate(segs):
        assert isinstance(s, dict), f"{rec.name}: corrupted_segments[{j}] must be an object"
        for key in ("bus", "first_frame", "last_frame"):
            assert _is_int(s.get(key)), f"{rec.name}: corrupted_segments[{j}].{key} must be a JSON integer"
        assert 0 <= s["first_frame"] <= s["last_frame"] <= rec.last_frame, \
            f"{rec.name}: corrupted_segments[{j}] must satisfy 0 <= first_frame <= last_frame <= {rec.last_frame}"
    return segs


# ---------------------------------------------------------------- tests
@pytest.mark.parametrize("name", RECORDINGS)
def test_program_run(runs, name):
    """Requirement 1: on this held-out recording the program exits with status 0 within the wall-clock limit and
    leaves identification.json as a regular file (not a symlink) of at most 1 MB in its working directory."""
    r = runs[name]
    assert "error" not in r, f"{name}: {r.get('error')}"


@pytest.mark.parametrize("name", RECORDINGS)
def test_schema(runs, name):
    """Requirement 2: strict JSON (in-range number literals, no repeated keys); one branches entry per candidate
    with an integer id and boolean in_service; parameters of every in-service, non-known entry present as
    JSON numbers; exactly the 22 channels with numeric ratio/phase errors; segments of JSON integers with
    0 <= first_frame <= last_frame <= frames - 1."""
    rec, data = REC[name], _load(runs, name)
    br = _branches(rec, data)
    for cid, b in br.items():
        c = rec.cands[cid]
        if b["in_service"] and not c["known"]:
            for key in PARAM_KEYS[c["model"]]:
                assert _is_num(b.get(key)), f"{name}: branch {cid}: {key} must be a JSON number when in service"
    _channels(rec, data)
    _segments(rec, data)


@pytest.mark.parametrize("name", RECORDINGS)
def test_topology(runs, name):
    """Requirement 3: in_service equals the true status of every candidate."""
    rec = REC[name]
    br = _branches(rec, _load(runs, name))
    wrong = [cid for cid in sorted(rec.cands) if br[cid]["in_service"] != rec.truth["status"][str(cid)]]
    assert not wrong, f"{name}: wrong in_service for candidates " + ", ".join(
        f"{cid} ({rec.cands[cid]['from_bus']}-{rec.cands[cid]['to_bus']}: reported {br[cid]['in_service']})"
        for cid in wrong)


@pytest.mark.parametrize("name", RECORDINGS)
def test_branch_parameters(runs, name):
    """Requirement 4: every unknown parameter of every truly in-service, non-known candidate is within its
    tolerance (exact decimal absolute difference, inclusive)."""
    rec = REC[name]
    br = _branches(rec, _load(runs, name))
    bad = []
    for cid_s, params in sorted(rec.truth["params"].items(), key=lambda kv: int(kv[0])):
        cid = int(cid_s)
        b = br[cid]
        assert b["in_service"], f"{name}: branch {cid} is in service; its parameters cannot be graded when reported out of service"
        for key, true in params.items():
            got = b.get(key)
            assert _is_num(got), f"{name}: branch {cid}: {key} must be a JSON number"
            tol = rec.tol["branches"][cid_s][key]
            if not _within(got, true, tol):
                bad.append(f"branch {cid} {key}: off by {abs(float(got) - true):.3g} (tolerance {tol:.3g})")
    assert not bad, f"{name}: " + "; ".join(bad)


@pytest.mark.parametrize("name", RECORDINGS)
def test_calibration(runs, name):
    """Requirement 5: every reported ratio and phase error is within its tolerance (exact decimal absolute
    difference, inclusive, no angle wrapping)."""
    rec = REC[name]
    ch = _channels(rec, _load(runs, name))
    bad = []
    for cname in rec.channels:
        for key in ("ratio_error", "phase_error"):
            true, tol = rec.truth["cal"][cname][key], rec.tol["channels"][cname][key]
            if not _within(ch[cname][key], true, tol):
                bad.append(f"{cname} {key}: off by {abs(float(ch[cname][key]) - true):.3g} (tolerance {tol:.3g})")
    assert not bad, f"{name}: " + "; ".join(bad)


@pytest.mark.parametrize("name", RECORDINGS)
def test_corrupted_segments(runs, name):
    """Requirement 6: exactly one segment per true malfunction interval, paired one-to-one with the intervals on
    the same bus with both endpoints within 30 frames."""
    rec = REC[name]
    segs = _segments(rec, _load(runs, name))
    true = rec.truth["corrupted"]
    assert len(segs) == len(true), \
        f"{name}: {len(segs)} corrupted_segments reported; the recording has {len(true)} malfunction intervals"

    def match(s, t):
        return (s["bus"] == t["bus"] and abs(s["first_frame"] - t["first_frame"]) <= SEG_TOL
                and abs(s["last_frame"] - t["last_frame"]) <= SEG_TOL)
    ok = any(all(match(segs[i], t) for i, t in zip(perm, true)) for perm in itertools.permutations(range(len(segs))))
    unmatched = [t["bus"] for t in true if not any(match(s, t) for s in segs)]
    assert ok, (f"{name}: no one-to-one pairing of reported segments with the {len(true)} malfunction intervals; "
                f"intervals with no matching segment are on buses {sorted(unmatched)}")
