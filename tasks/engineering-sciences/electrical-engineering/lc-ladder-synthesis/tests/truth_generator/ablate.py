"""Ablations: the verifier's own checks applied to correct and to deliberately wrong ladders.

usage: python3 ablate.py solution_dir [tests_dir [spec_id ...]]

For every held-out specification the reference solution's search is used to enumerate the
canonical ladders; the verifier functions (imported from tests_dir, paths redirected) are
then applied to:
  - the reference program's output (run through the verifier's sandboxed runner);
  - the stored witness, and the witness with every element scaled by 1 + 1e-12 (must pass);
  - one inductor scaled by 1 + 1e-4, and the reference function with its reflection zeros
    moved by 1e-6 relative (must fail);
  - ladders synthesised exactly for other transmission zeros (must fail): the reference
    zeros moved by 1e-7 relative, and the zeros of the starting guess (spaced out from the
    stopband edge but not equiripple);
  - the first ladder a depth-first search finds, without any choice of the realisation;
  - a canonical ladder of the right response whose load is outside [R_s / 2, 2 R_s];
  - format errors.
It also reports, per specification, whether any ladder exists when zeros are shifted only
by partial removal of poles at infinity (the textbook elliptic-filter procedure).
"""
import copy
import json
import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SOL = sys.argv[1]
TESTS = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(HERE)
sys.path.insert(0, TESTS)
sys.path.insert(1, SOL)
import test_synthesis as T  # noqa: E402
import synthesize as SY     # noqa: E402
import mpmath as mp         # noqa: E402

T.SPECS = os.path.join(TESTS, 'specs')
T.TRUTH = os.path.join(TESTS, 'truth.json')
TRUTH = json.load(open(T.TRUTH))
IDS = sys.argv[3:] or sorted(TRUTH)


def verdict(spec, rec, net):
    try:
        RL, arms = T.check_contract(net)
        T.check_network(spec, rec, RL, arms)
        return True, ''
    except (AssertionError, ValueError) as e:
        return False, str(e)[:100]


def resynthesise(spec, uz, limit=4000):
    """Canonical ladder (fewest inductors, admissible load) of the passband-equiripple function
    with the transmission zeros uz (normalised u = (f / f_2)^2), using the reference solution."""
    SY.set_precision(spec)
    SY.ADMISSIBLE_ONLY = True
    limit = min(limit, 40)
    ap = SY.Approx(spec)
    ap.place_zeros = lambda: ap.solve_for(uz)
    F, P, E = ap.build()
    sols = []
    for sg in (1, -1):
        N = SY.trim(SY.padd(E, SY.pscale(F, sg)))
        D = SY.trim(SY.padd(E, SY.pscale(F, -sg)))
        out = []
        SY.search(N, D, list(ap.uz), [], out, limit)
        sols += [(sg, p, R) for p, R in out]
    SY.ADMISSIBLE_ONLY = False
    best = SY.choose(spec, sols)
    return SY.denormalise(spec, best[1], best[2]) if best else None


def uses_dc_shift(path):
    for i in range(len(path) - 1):
        if path[i + 1][1] in ('series_LC', 'parallel_LC') and path[i][:2] in (('series', 'C'), ('shunt', 'L')):
            return True
    return False


def main():
    rows = []
    tally = {}
    for sid in IDS:
        spec = json.load(open(os.path.join(T.SPECS, sid + '.json')))
        rec = TRUTH[sid]
        out = T.run_program(sid, program=os.path.join(SOL, 'synthesize.py'))
        net = T.read_network(out)
        rows.append((sid, 'reference program output', True) + verdict(spec, rec, net))
        w = rec['witness']
        rows.append((sid, 'stored witness', True) + verdict(spec, rec, w))
        x = copy.deepcopy(w)
        for arm in x['arms']:
            for k in ('L', 'C'):
                if k in arm:
                    arm[k] *= 1 + 1e-12
        x['load_resistance_ohm'] *= 1 + 1e-12
        rows.append((sid, 'witness, every value x (1 + 1e-12)', True) + verdict(spec, rec, x))
        x = copy.deepcopy(w)
        next(a for a in x['arms'] if 'L' in a)['L'] *= 1 + 1e-4
        rows.append((sid, 'one inductor x (1 + 1e-4)', False) + verdict(spec, rec, x))
        r2 = copy.deepcopy(rec)
        r2['reflection_zeros_u'] = [str(mp.mpf(q) * (1 + mp.mpf('1e-6'))) for q in rec['reflection_zeros_u']]
        rows.append((sid, 'reflection zeros moved by 1e-6 (approximation error)', False) + verdict(spec, r2, w))
        if rec['transmission_zeros_u']:
            moved = [mp.mpf(z) * (1 + mp.mpf('1e-7')) for z in rec['transmission_zeros_u']]
            rows.append((sid, 'ladder for the transmission zeros moved by 1e-7', False)
                        + verdict(spec, rec, resynthesise(spec, moved)))
            SY.set_precision(spec)
            ap0 = SY.Approx(spec)
            p0 = []
            for ue, m in ap0.bands:
                sc = abs(float(mp.log(ue)))
                p0 += [mp.log(0.05 * sc)] + [mp.log(0.5 * sc)] * (m - 1)
            rows.append((sid, 'ladder for the starting-guess zeros (not equiripple)', False)
                        + verdict(spec, rec, resynthesise(spec, ap0.zeros_from(p0))))
        # all canonical ladders up to degree 12; for degree 14 and 16 the first 3000 per sign
        full = spec['degree'] <= 12
        ap, sols = SY.synthesize(spec) if full else SY.synthesize(spec, limit=3000)
        first = SY.denormalise(spec, sols[0][1], sols[0][2])
        ok, why = verdict(spec, rec, first)
        tally.setdefault('first ladder found passes', []).append(ok)
        bad_load = [(p, R) for _, p, R in sols if not 0.5 <= R <= 2]
        if bad_load:
            p_, R_ = max(bad_load, key=lambda t: abs(mp.log(t[1])))
            ok2, why2 = verdict(spec, rec, SY.denormalise(spec, p_, R_))
            rows.append((sid, 'canonical ladder with load %.3g R_s (outside [1/2, 2])' % float(R_), False, ok2, why2))
        if full:
            inf_only = [p for _, p, _ in sols if not uses_dc_shift(p)]
            tally.setdefault('a ladder exists with zero shifting from infinity only (degree <= 12)', []).append(bool(inf_only))
        rows.append((sid, 'first ladder of the search (no choice of realisation): %s' % ('passes' if ok else 'fails'),
                     None, ok, why))
    # format errors on one specification
    sid = IDS[0]
    spec = json.load(open(os.path.join(T.SPECS, sid + '.json')))
    rec, w = TRUTH[sid], TRUTH[sid]['witness']
    bad = [
        ('extra top-level key', lambda d: d.update(note=1)),
        ('missing load resistance', lambda d: d.pop('load_resistance_ohm')),
        ('negative element', lambda d: d['arms'][0].update({k: -v for k, v in d['arms'][0].items() if k in 'LC'})),
        ('boolean element value', lambda d: d['arms'][0].update({k: True for k in d['arms'][0] if k in 'LC'})),
        ('unknown arm type', lambda d: d['arms'][0].update(type='transformer')),
        ('extra key in an arm', lambda d: d['arms'][0].update(R=1.0)),
        ('empty ladder', lambda d: d.update(arms=[])),
        ('element value as a string', lambda d: d['arms'][0].update({k: str(v) for k, v in d['arms'][0].items() if k in 'LC'})),
    ]
    for name, fn in bad:
        d = copy.deepcopy(w)
        fn(d)
        rows.append((sid, name, False) + verdict(spec, rec, d))
    tmp = tempfile.mkdtemp()
    p = os.path.join(tmp, 'n.json')
    open(p, 'w').write(json.dumps(w).replace('"load_resistance_ohm": ', '"load_resistance_ohm": NaN, "x": ', 1))
    try:
        T.read_network(p)
        rows.append((sid, 'NaN in the file', False, True, ''))
    except ValueError as e:
        rows.append((sid, 'NaN in the file', False, False, str(e)))
    q = os.path.join(tmp, 'l.json')
    json.dump(w, open(p, 'w'))
    os.symlink(p, q)
    try:
        T.read_network(q)
        rows.append((sid, 'symlink to a correct file', False, True, ''))
    except AssertionError as e:
        rows.append((sid, 'symlink to a correct file', False, False, str(e)))
    good = 0
    graded = 0
    for sid, name, expect, got, why in rows:
        if expect is None:
            print('info %-10s %-62s %s' % (sid, name, why[:60]))
            continue
        graded += 1
        good += expect == got
        print('%-4s %-10s %-62s expect=%-5s got=%-5s %s' % ('ok' if expect == got else 'BAD', sid, name, expect, got, why[:70]))
    for k, v in tally.items():
        print('%s: %d of %d specifications' % (k, sum(v), len(v)))
    print('%d/%d ablations behave as expected' % (good, graded))


if __name__ == '__main__':
    main()
