"""Verifier for lc-ladder-synthesis.

The agent's program /root/results/synthesize.py is run on each held-out specification in
/tests/specs (never shown to the agent) as an unprivileged user, without network, in an
empty working directory, and the network file it writes is read as data.  The network is analysed
here in 40-digit arithmetic and compared with the reference characteristic function in
/tests/truth.json (generated and certified by /tests/truth_generator/make_truth.py).
"""
import json
import math
import os
import resource
import shutil
import signal
import stat
import subprocess
import tempfile

import mpmath as mp
import pytest

PROGRAM = '/root/results/synthesize.py'
SPECS = '/tests/specs'
TRUTH = '/tests/truth.json'

TIME_LIMIT = 300            # seconds per specification
MAX_BYTES = 1 << 20
TOL_DB = 1e-5               # insertion-loss tolerance
NOBODY = 65534
POSITIONS = {'series', 'shunt'}
TYPES = {'L': ('L',), 'C': ('C',), 'series_LC': ('L', 'C'), 'parallel_LC': ('L', 'C')}

SPEC_IDS = sorted(f[:-5] for f in os.listdir(SPECS) if f.endswith('.json')) if os.path.isdir(SPECS) else []


# ------------------------------------------------------------- running ---

def _limits():
    """Applied in the child before exec, and inherited by every thread and process it starts:
    one CPU core, and a 3 GB address-space limit per process."""
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    resource.setrlimit(resource.RLIMIT_AS, (3 << 30, 3 << 30))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


def _kill_group(proc):
    """Kill the program and every process it started (they share its process group)."""
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


def run_program(spec_id, program=PROGRAM, time_limit=TIME_LIMIT, as_nobody=True):
    """Run the program on one specification; return the path of the network file it wrote.

    The program, the specification and an empty working directory are in separate folders;
    the program is called with absolute paths.  Its output goes to files, not pipes, and the
    whole process group is killed when the time limit expires and again after it exits, so a
    left-over child process cannot hold the verifier.
    """
    if not os.path.isfile(program):
        pytest.fail('%s: %s does not exist' % (spec_id, program))
    base = tempfile.mkdtemp(prefix='synth-')
    prog_dir, in_dir, work = [os.path.join(base, d) for d in ('prog', 'in', 'work')]
    for d in (prog_dir, in_dir, work):
        os.mkdir(d)
    prog = os.path.join(prog_dir, 'synthesize.py')
    spec = os.path.join(in_dir, spec_id + '.json')
    out = os.path.join(work, spec_id + '.network.json')
    shutil.copyfile(program, prog)
    shutil.copyfile(os.path.join(SPECS, spec_id + '.json'), spec)
    kw = {}
    if as_nobody and os.geteuid() == 0:
        for root, dirs, files in os.walk(base):
            for x in dirs + files:
                os.chown(os.path.join(root, x), NOBODY, NOBODY)
        os.chown(base, NOBODY, NOBODY)
        kw = dict(user=NOBODY, group=NOBODY, extra_groups=[])
    env = {'PATH': '/usr/local/bin:/usr/bin:/bin', 'HOME': work, 'TMPDIR': work,
           'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
           'PYTHONDONTWRITEBYTECODE': '1'}
    log = tempfile.TemporaryFile()
    proc = subprocess.Popen(['python3', prog, spec, out], cwd=work, env=env, stdin=subprocess.DEVNULL,
                            stdout=subprocess.DEVNULL, stderr=log, preexec_fn=_limits,
                            start_new_session=True, **kw)
    try:
        code = proc.wait(timeout=time_limit)
    except subprocess.TimeoutExpired:
        _kill_group(proc)
        proc.wait()
        pytest.fail('%s: the program did not finish within %d s' % (spec_id, time_limit))
    finally:
        _kill_group(proc)
    if code != 0:
        log.seek(0)
        pytest.fail('%s: the program exited with status %d: %s' % (spec_id, code,
                                                                   log.read()[-600:].decode(errors='replace')))
    return out


# ------------------------------------------------------------- reading ---

def _reject_constant(name):
    raise ValueError('non-finite JSON constant %s' % name)


def _unique_keys(pairs):
    keys = [k for k, _ in pairs]
    if len(set(keys)) != len(keys):
        raise ValueError('repeated key in a JSON object')
    return dict(pairs)


def read_network(path):
    st = os.lstat(path)
    if not stat.S_ISREG(st.st_mode):
        raise AssertionError('network.json is not a regular file')
    if st.st_size > MAX_BYTES:
        raise AssertionError('network.json is larger than 1 MiB')
    fd = os.open(path, os.O_RDONLY | os.O_NONBLOCK | getattr(os, 'O_NOFOLLOW', 0))
    try:
        data = os.read(fd, MAX_BYTES + 1)
    finally:
        os.close(fd)
    return json.loads(data.decode('utf-8'), parse_constant=_reject_constant, object_pairs_hook=_unique_keys)


def _value(x, what):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or not x > 0:
        raise AssertionError('%s must be a positive finite number, got %r' % (what, x))
    return float(x)


def check_contract(net):
    if not isinstance(net, dict) or set(net) != {'load_resistance_ohm', 'arms'}:
        raise AssertionError('network.json must be an object with exactly the keys load_resistance_ohm and arms')
    RL = _value(net['load_resistance_ohm'], 'load_resistance_ohm')
    arms = net['arms']
    if not isinstance(arms, list) or not arms:
        raise AssertionError('arms must be a non-empty list')
    out = []
    for i, arm in enumerate(arms):
        if not isinstance(arm, dict) or arm.get('position') not in POSITIONS or arm.get('type') not in TYPES:
            raise AssertionError('arm %d: position must be series or shunt and type one of %s' % (i, sorted(TYPES)))
        need = TYPES[arm['type']]
        if set(arm) != {'position', 'type'} | set(need):
            raise AssertionError('arm %d: a %s arm has exactly the keys position, type, %s' % (i, arm['type'], ', '.join(need)))
        out.append((arm['position'], arm['type'], {k: _value(arm[k], 'arm %d %s' % (i, k)) for k in need}))
    return RL, out


# ------------------------------------------------------------ analysis ---

def gain_network(arms, Rs, RL, f):
    """Transducer power gain 4 Rs RL / |A RL + B + C Rs RL + D Rs|^2 of the ladder."""
    s = 2j * mp.pi * mp.mpf(f)
    A, B, C, D = mp.mpc(1), mp.mpc(0), mp.mpc(0), mp.mpc(1)
    for pos, typ, v in arms:
        if typ == 'L':
            z = s * v['L']
        elif typ == 'C':
            z = 1 / (s * v['C'])
        elif typ == 'series_LC':
            z = s * v['L'] + 1 / (s * v['C'])
        else:
            z = 1 / (1 / (s * v['L']) + s * v['C'])
        if pos == 'series':
            B, D = A * z + B, C * z + D
        else:
            y = 1 / z
            A, C = A + B * y, C + D * y
    den = A * RL + B + C * Rs * RL + D * Rs
    return 4 * Rs * RL / abs(den) ** 2


def gain_reference(rec, f):
    u = (mp.mpf(f) / mp.mpf(rec['f2_hz'])) ** 2
    k2 = mp.mpf(rec['c2']) * u ** rec['r'] / u ** rec['dc']
    for q in rec['reflection_zeros_u']:
        k2 *= (u - mp.mpf(q)) ** 2
    for z in rec['transmission_zeros_u']:
        k2 /= (u - mp.mpf(z)) ** 2
    return 1 / (1 + k2)


def near_zero_window(spec):
    """Relative half-width w of the window around a finite transmission zero in which grid points
    are omitted: 0.5 %, scaled down in proportion for band-pass filters narrower than 1.6 %
    (w = 0.5 % b / 1.6 %, with b = (f2 - f1) / sqrt(f1 f2)), so that it stays a fixed fraction of
    the band and never covers the passband."""
    f1, f2 = float(spec['passband_hz'][0]), float(spec['passband_hz'][1])
    if f1 == 0:
        return 5e-3
    b = (f2 - f1) / math.sqrt(f1 * f2)
    return 5e-3 * min(1.0, b / 0.016)


def _frequencies(spec, rec):
    """1001 passband points (401 uniform, 600 spaced geometrically towards f1 and f2); 120
    log-spaced points from f2 up to 30 f2 and, for a band-pass filter, 120 from f1 down to f1 / 30;
    grid points within w of a reference zero omitted (w = 0.5 %, smaller for bands narrower than
    1.6 %, see near_zero_window); points 2 w and 6 w (1 % and 3 % when w = 0.5 %) on either side of
    each reference transmission zero; and the stopband edges, which are always checked, also when
    an edge lies within w of a zero."""
    w = near_zero_window(spec)
    f1, f2 = float(spec['passband_hz'][0]), float(spec['passband_hz'][1])
    fs = [f1 + (f2 - f1) * k / 400 for k in range(401)]
    # the ripples crowd towards the band edges (for the widest band-pass filters the lowest lies
    # within a few % of f1, 1e-5 of the band): 400 points spaced geometrically in f - f1 and 200
    # in f2 - f, from 1e-8 (f2 - f1)
    fs += [f1 + (f2 - f1) * 10 ** (-8 + 8 * k / 400) for k in range(400)]
    fs += [f2 - (f2 - f1) * 10 ** (-8 + 7 * k / 200) for k in range(200)]
    if f1 == 0:
        fs[0] = f2 * 1e-4
    else:
        fs += [f1 * 30 ** (-k / 120) for k in range(1, 121)]
    fs += [f2 * 30 ** (k / 120) for k in range(1, 121)]
    zs = rec['transmission_zeros_hz']
    fs = [f for f in fs if all(abs(f / fz - 1) >= w for fz in zs)]
    for fz in zs:
        fs += [fz * (1 + e) for e in (-6 * w, -2 * w, 2 * w, 6 * w)]
    # a point placed next to one zero must not fall inside the window of another
    fs = [f for f in fs if all(abs(f / fz - 1) >= w * (1 - 1e-9) for fz in zs)]
    # the stopband edges are required points (the instruction names them): never omitted
    return fs + [float(b['edge_hz']) for b in spec['transmission_zeros']['stopbands']]


def check_network(spec, rec, RL, arms):
    """Raise AssertionError with a reason if the ladder does not meet the specification."""
    Rs = float(spec['source_resistance_ohm'])
    n_el = sum(len(TYPES[t]) for _, t, _ in arms)
    canonical = int(spec['degree']) + sum(int(b['zeros']) for b in spec['transmission_zeros']['stopbands'])
    if n_el > canonical:
        raise AssertionError('%d reactive elements; a canonical ladder of this specification has %d'
                             % (n_el, canonical))
    if not Rs / 2 <= RL <= 2 * Rs:
        raise AssertionError('load resistance %.6g ohm is outside [R_s / 2, 2 R_s]' % RL)
    worst, where = 0.0, None
    with mp.workdps(40):
        for f in _frequencies(spec, rec):
            g = gain_network(arms, Rs, RL, f)
            g0 = gain_reference(rec, f)
            if not g > 0:
                raise AssertionError('zero transmission at %.9g Hz' % f)
            err = abs(float(10 * mp.log10(g / g0)))
            if err > worst:
                worst, where = err, (f, float(-10 * mp.log10(g0)), float(-10 * mp.log10(g)))
    if worst > TOL_DB:
        raise AssertionError('insertion loss off by %.3g dB at %.9g Hz (reference %.6f dB, network %.6f dB)'
                             % (worst, where[0], where[1], where[2]))
    return worst


# --------------------------------------------------------------- tests ---

def test_program_present():
    """The deliverable /root/results/synthesize.py exists."""
    assert os.path.isfile(PROGRAM), '%s does not exist' % PROGRAM


@pytest.mark.parametrize('spec_id', SPEC_IDS)
def test_filter(spec_id):
    """The program's ladder for a held-out specification has the reference insertion loss to 1e-5 dB,
    at most the canonical number n + m of elements, and a load in [R_s / 2, 2 R_s]."""
    spec = json.load(open(os.path.join(SPECS, spec_id + '.json')))
    rec = json.load(open(TRUTH))[spec_id]
    path = run_program(spec_id)
    try:
        net = read_network(path)
        RL, arms = check_contract(net)
        check_network(spec, rec, RL, arms)
    except (AssertionError, ValueError, OSError, UnicodeDecodeError) as e:
        pytest.fail('%s: %s' % (spec_id, e))
