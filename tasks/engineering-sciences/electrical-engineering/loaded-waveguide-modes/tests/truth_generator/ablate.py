"""Ablations: the real verifier functions applied to correct, perturbed and wrong submissions.

usage: python3 ablate.py oracle.json alt.json [tests_dir]

oracle.json and alt.json are two correct submissions (solution/solve.py and alt_solver.py).
Wrong submissions are built in two ways:
  * targeted edits of the oracle output (one physics shortcut each);
  * whole submissions from the oracle's solver under a wrong assumption: first-order loss
    perturbation of the nearly lossless guide, the Np/m - dB/m confusion in the threshold,
    the opposite time convention.
The verifier module is imported from tests_dir (default: the parent folder) with its file paths
redirected to temporary files.
"""
import cmath
import copy
import json
import math
import os
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(HERE)
sys.path.insert(0, TESTS)
sys.path.insert(0, os.path.join(os.path.dirname(TESTS), 'solution'))
import test_modes as T  # noqa: E402
import solve as SO       # noqa: E402

T.GUIDES = os.path.join(TESTS, 'data', 'guides.json')
T.TRUTH = os.path.join(TESTS, 'truth.json')
TRUTH = json.load(open(T.TRUTH))
GUIDES = json.load(open(T.GUIDES))['guides']
GBYID = {g['id']: g for g in GUIDES}
TEST_FUNCS = [T.test_mode_families, T.test_propagation_constants]
A = 10 * math.log(10)


def verdict(path):
    """True when every verifier test passes on the file at path."""
    T.OUTPUT = path
    try:
        out = T._read_output()
        T._check_contract(out, GUIDES)
        res = out['guides']
        for f in TEST_FUNCS:
            args = {'output': res, 'truth': TRUTH, 'guides': GUIDES}
            names = f.__code__.co_varnames[:f.__code__.co_argcount]
            f(**{k: args[k] for k in names})
        return True, ''
    except BaseException as e:  # pytest.fail raises Failed (BaseException)
        return False, str(e).splitlines()[0][:110] if str(e) else type(e).__name__


def each_mode(S):
    for gid, r in S.items():
        for j, e in enumerate(r['modes']):
            for dr in ('forward', 'reverse'):
                for fam, lst in e[dr].items():
                    for p in lst:
                        yield gid, j, fam, p


def kn2(gid, fam):
    n = int(fam[3:])
    return (n * math.pi / float(GBYID[gid]['b'])) ** 2


def drop_where(S, pred):
    for gid, r in S.items():
        for j, e in enumerate(r['modes']):
          for m in (e['forward'], e['reverse']):
            for fam in list(m):
                m[fam] = [p for p in m[fam] if not pred(gid, j, fam, p)]
                if not m[fam]:
                    del m[fam]


def surface(gid, j, fam, p):
    """Re q above k0^2 max Re eps of the dielectric layers: a surface wave of a plasma interface."""
    g = GBYID[gid]
    f = float(g['test_frequencies'][j])
    q = kn2(gid, fam) - complex(p[1], p[0]) ** 2
    emax = max(float(L['material']['eps']) for L in g['layers'] if L['material']['type'] == 'dielectric')
    return q.real > (2 * math.pi * f / 299792458.0) ** 2 * emax


def evanescent(gid, j, fam, p):
    """|beta| < alpha: an evanescent mode below its cutoff (loss small compared with the decay)."""
    return abs(p[0]) < p[1]


# ------------------------------------------------------------------ alternative physics
def scaled_guides(s):
    out = []
    for g in GUIDES:
        g = copy.deepcopy(g)
        for L in g['layers']:
            m = L['material']
            if m['type'] in ('dielectric', 'ferrite'):
                m['loss_tangent'] = float(m['loss_tangent']) * s
            if m['type'] == 'ferrite':
                m['linewidth_oe'] = float(m['linewidth_oe']) * s
            if m['type'] == 'plasma':
                m['collision_frequency'] = float(m['collision_frequency']) * s
        out.append(g)
    return out


def solve_all(guides, a_lim):
    SO.A, SO.A1 = a_lim, a_lim * 1.02
    try:
        return {g['id']: SO.solve_guide(g) for g in guides}
    finally:
        SO.A, SO.A1 = A, A * 1.02


def perturbative():
    """First-order loss perturbation: gamma(1) ~ gamma(s1) + (gamma(s2) - gamma(s1)) (1 - s1) / (s2 - s1)
    from two nearly lossless solves; then the 200 dB/m limit."""
    s1, s2 = 1e-3, 2e-3
    r1 = solve_all(scaled_guides(s1), 3 * A)
    r2 = solve_all(scaled_guides(s2), 3 * A)
    out = {}
    for gid in r1:
        modes = []
        for e1, e2 in zip(r1[gid]['modes'], r2[gid]['modes']):
            ent = {}
            for dr in ('forward', 'reverse'):
                m1, m2, e = e1[dr], e2[dr], {}
                for fam, l1 in m1.items():
                    l2 = m2.get(fam, [])
                    for b1, a1 in l1:
                        g1 = complex(a1, b1)
                        if not l2:
                            continue
                        g2 = min((complex(a, b) for b, a in l2), key=lambda z: abs(z - g1))
                        g = g1 + (g2 - g1) * (1 - s1) / (s2 - s1)
                        if 0 < g.real <= A:
                            e.setdefault(fam, []).append([g.imag, g.real])
                ent[dr] = e
            modes.append(ent)
        out[gid] = {'modes': modes}
    return out


FERR = ['G13', 'G14', 'G15', 'G16']


def no_linewidth():
    """Ferrite linewidth ignored (Delta H = 0): only the dielectric loss of the ferrite remains."""
    gs = []
    for g in GUIDES:
        g = copy.deepcopy(g)
        for L in g['layers']:
            if L['material']['type'] == 'ferrite':
                L['material']['linewidth_oe'] = 0
        gs.append(g)
    return solve_all(gs, A)


def reciprocal(S):
    for gid in FERR:
        for e in S[gid]['modes']:
            e['reverse'] = copy.deepcopy(e['forward'])


def swap_dirs(S):
    for gid in FERR:
        for e in S[gid]['modes']:
            e['forward'], e['reverse'] = e['reverse'], e['forward']


def plan(base):
    P = []
    P.append(('gamma x (1 + 5e-8) for every mode (inside 1e-7)', True,
              lambda S: [p.__setitem__(0, p[0] * (1 + 5e-8)) or p.__setitem__(1, p[1] * (1 + 5e-8)) for *_, p in each_mode(S)]))
    P.append(('one gamma x (1 + 3e-7) (G04, 14 GHz)', False,
              lambda S: (lambda p: p.__setitem__(0, p[0] * (1 + 3e-7)))(S['G04']['modes'][1]['forward']['LSE0'][0])))
    P.append(('beta > 0 for every mode (lossless habit: backward waves get the wrong sign)', False,
              lambda S: [p.__setitem__(0, abs(p[0])) for *_, p in each_mode(S)]))
    P.append(('beta of the backward waves only flipped in G11', False,
              lambda S: [p.__setitem__(0, abs(p[0])) for gid, j, fam, p in each_mode(S) if gid == 'G11']))
    P.append(('evanescent modes below cutoff left out (|beta| < alpha)', False, lambda S: drop_where(S, evanescent)))
    P.append(('evanescent modes left out in G06 only', False,
              lambda S: drop_where({'G06': S['G06']}, evanescent)))
    P.append(('surface waves above k0^2 max eps not searched', False, lambda S: drop_where(S, surface)))
    P.append(('G12: the 5.15 GHz surface wave left out', False, lambda S: [S['G12']['modes'][1][d].clear() for d in ('forward', 'reverse')]))
    P.append(('G04 at 18 GHz: one of the three LSM1 modes left out', False,
              lambda S: S['G04']['modes'][3]['forward']['LSM1'].pop()))
    P.append(('G07 at 8.3 GHz: an evanescent LSE0 mode above the limit reported', False,
              lambda S: S['G07']['modes'][2]['forward'].update(LSE0=[[4.0, 1.05 * A]])))
    P.append(('LSM2 of G08 at 4.4 GHz reported as LSM1 (family index slip)', False,
              lambda S: (lambda m: m.setdefault('LSM1', []).extend(m.pop('LSM2')))(S['G08']['modes'][0]['forward'])))
    P.append(('opposite time convention exp(-j w t) (beta of every mode flipped)', False,
              lambda S: [p.__setitem__(0, -p[0]) for *_, p in each_mode(S)]))
    P.append(('ferrite guides treated as reciprocal (reverse modes copied from forward)', False, reciprocal))
    P.append(('Polder tensor with the wrong sign of kappa (forward and reverse swapped)', False, swap_dirs))
    P.append(('G14 only: reverse copied from forward (twin-slab differential phase missed)', False,
              lambda S: [e.__setitem__('reverse', copy.deepcopy(e['forward'])) for e in S['G14']['modes']]))
    P.append(('G15 at 11 GHz: the isolated reverse direction given the forward mode', False,
              lambda S: S['G15']['modes'][2].__setitem__('reverse', copy.deepcopy(S['G15']['modes'][2]['forward']))))
    P.append(('ferrite guide: a y-dependent LSE1 family reported (G13, 12 GHz)', False,
              lambda S: S['G13']['modes'][3]['forward'].update(LSE1=[[100.0, 1.0]])))
    return P


FORMAT = [
    ('extra top-level key', False, lambda d: d.update(note=1)),
    ('missing guide', False, lambda d: d['guides'].pop('G16')),
    ('extra key in a guide', False, lambda d: d['guides']['G01'].update(extra=0)),
    ('LSM0 family name', False, lambda d: d['guides']['G01']['modes'][0]['forward'].update(LSM0=[[100.0, 1.0]])),
    ('empty family list', False, lambda d: d['guides']['G01']['modes'][0]['forward'].update(LSM1=[])),
    ('mode entry with three numbers', False, lambda d: d['guides']['G01']['modes'][0]['forward']['LSE0'][0].append(1.0)),
    ('reverse key missing', False, lambda d: d['guides']['G01']['modes'][0].pop('reverse')),
    ('entry given as a bare family map (no forward/reverse)', False,
     lambda d: d['guides']['G01']['modes'].__setitem__(0, d['guides']['G01']['modes'][0]['forward'])),
    ('modes list too short', False, lambda d: d['guides']['G01']['modes'].pop()),
    ('alpha < 0 (the -gamma member; alpha and beta both negated)', False,
     lambda d: (lambda p: (p.__setitem__(0, -p[0]), p.__setitem__(1, -p[1])))(d['guides']['G02']['modes'][0]['forward']['LSE0'][0])),
    ('family name with a trailing newline', False,
     lambda d: d['guides']['G01']['modes'][0]['forward'].update({'LSE0\n': d['guides']['G01']['modes'][0]['forward'].pop('LSE0')})),
    ('family name with a leading zero', False,
     lambda d: d['guides']['G01']['modes'][0]['forward'].update({'LSE00': d['guides']['G01']['modes'][0]['forward'].pop('LSE0')})),
    ('beta as a string', False, lambda d: d['guides']['G01']['modes'][0]['forward']['LSE0'][0].__setitem__(0, '133.9')),
    ('boolean as alpha', False, lambda d: d['guides']['G01']['modes'][0]['forward']['LSE0'][0].__setitem__(1, True)),
]


def main():
    base = json.load(open(sys.argv[1]))
    alt = json.load(open(sys.argv[2]))
    tmp = tempfile.mkdtemp()
    p = os.path.join(tmp, 'case.json')
    rows = []

    def run(name, expect, sub):
        json.dump(sub, open(p, 'w'))
        rows.append((name, expect) + verdict(p))

    run('oracle output (solution/solve.py)', True, base)
    run('second implementation (alt_solver.py, Chebyshev collocation)', True, alt)
    for name, expect, fn in plan(base):
        d = copy.deepcopy(base)
        fn(d['guides'])
        run(name, expect, d)
    t0 = time.time()
    run('first-order loss perturbation of the nearly lossless guide', False, {'guides': perturbative()})
    run('limit read as 200 Np/m instead of 200 dB/m', False, {'guides': solve_all(GUIDES, 200.0)})
    run('limit 220 dB/m (no mode lies between 200 and 220 dB/m at any test frequency)', True,
        {'guides': solve_all(GUIDES, 220 / (20 / math.log(10)))})
    run('limit 300 dB/m instead of 200 dB/m', False, {'guides': solve_all(GUIDES, 300 / (20 / math.log(10)))})
    run('ferrite linewidth ignored (Delta H = 0)', False, {'guides': no_linewidth()})
    print('(alternative-physics solves: %.0f s)' % (time.time() - t0), flush=True)
    for name, expect, fn in FORMAT:
        d = copy.deepcopy(base)
        fn(d)
        run(name, expect, d)
    d = copy.deepcopy(base)
    d['guides']['G01']['modes'][0]['forward'] = {('LSE%d' % k): [[1.0 + 1e-3 * i, 1.0] for i in range(1000)] for k in range(400)}
    hp = os.path.join(tmp, 'huge.json')
    json.dump(d, open(hp, 'w'))
    t0 = time.time()
    v = verdict(hp)
    rows.append(('400 families x 1000 modes (%.1f MB, %.2fs)' % (os.path.getsize(hp) / 1e6, time.time() - t0), False) + v)
    q = os.path.join(tmp, 'nan.json')
    open(q, 'w').write(json.dumps(base).replace('"LSE0": [[', '"LSE0": [[NaN, 1.0], [', 1))
    rows.append(('NaN in a list', False) + verdict(q))
    q = os.path.join(tmp, 'dup.json')
    open(q, 'w').write(json.dumps(base).replace('{"guides": {', '{"guides": {"G01": {}, ', 1))
    rows.append(('repeated guide key (JSON object with a duplicate key)', False) + verdict(q))
    q = os.path.join(tmp, 'link.json')
    json.dump(base, open(p, 'w'))
    os.symlink(p, q)
    rows.append(('symlink to a correct file', False) + verdict(q))
    q = os.path.join(tmp, 'fifo.json')
    os.mkfifo(q)
    rows.append(('FIFO (must not block)', False) + verdict(q))
    q = os.path.join(tmp, 'big.json')
    with open(q, 'w') as fh:
        fh.write(' ' * (16 * 1024 * 1024 + 1))
    rows.append(('file over 16 MiB', False) + verdict(q))
    rows.append(('no file (nop)', False) + verdict(os.path.join(tmp, 'absent.json')))
    q = os.path.join(tmp, 'empty.json')
    open(q, 'w').close()
    rows.append(('empty file', False) + verdict(q))
    ok = 0
    for name, expect, got, why in rows:
        good = expect == got
        ok += good
        print('%-4s %-74s expect=%-5s got=%-5s %s' % ('ok' if good else 'BAD', name, expect, got, why))
    print('%d/%d ablations behave as expected' % (ok, len(rows)))


if __name__ == '__main__':
    main()
