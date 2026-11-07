"""Ablations on the held-out guides: the real verifier functions applied to correct and to edited outputs.

usage: python3 ablate.py tests_dir oracle_dir

tests_dir holds test_modes.py and data/run_*/ (heldout_guides.json, heldout_truth.json); oracle_dir holds a
correct output oracle_<run>.json for every run.  Each case edits the correct outputs in one way and is judged
by the verifier's own parsing and value tests (the program-level cases, such as the oracle itself, the second
solver, narrower search boxes, the v4 reference program and no program, are run through the full verifier).
"""
import copy
import json
import math
import os
import sys

TESTS, ORC = sys.argv[1], sys.argv[2]
sys.path.insert(0, TESTS)
import test_modes as T  # noqa: E402

A = 10 * math.log(10)
OUT = {r: json.load(open(os.path.join(ORC, 'oracle_%s.json' % r))) for r in T.RUNS}
GBYID = {g['id']: g for r in T.RUNS for g in T.GUIDES[r]}
RUN_OF = {g['id']: r for r in T.RUNS for g in T.GUIDES[r]}


def verdict(outs):
    """(passes, first failure line) of the verifier's contract and value tests on {run: output dict}."""
    for r in T.RUNS:
        raw = json.dumps(outs[r]).encode() if not isinstance(outs[r], bytes) else outs[r]
        runs = {r: {'returncode': 0, 'elapsed': 1.0, 'raw': raw}}
        for f in (T.test_output_contract, T.test_mode_families, T.test_propagation_constants, T.test_power_split):
            try:
                f(runs, r)
            except BaseException as e:   # AssertionError, or pytest's Failed
                return False, '%s: %s' % (f.__name__, str(e).splitlines()[0][:150] if str(e) else type(e).__name__)
    return True, ''


def edited(fn):
    outs = copy.deepcopy(OUT)
    fn(outs)
    return outs


def each(outs, gids=None):
    for r in T.RUNS:
        for gid, rec in outs[r]['guides'].items():
            if gids and gid not in gids:
                continue
            for j, e in enumerate(rec['modes']):
                yield gid, j, e


def drop(outs, pred):
    for gid, j, e in each(outs):
        for dr in ('forward', 'reverse'):
            for fam in list(e[dr]):
                e[dr][fam] = [m for m in e[dr][fam] if not pred(gid, j, dr, fam, m)]
                if not e[dr][fam]:
                    del e[dr][fam]


def k0(gid, j):
    return 2 * math.pi * float(GBYID[gid]['test_frequencies'][j]) / 299792458.0


def nfam(fam):
    return int(fam[3:])


CASES = []


def case(expect):
    def deco(fn):
        CASES.append((fn.__doc__.strip(), expect, fn))
        return fn
    return deco


# ---------------------------------------------------------------- must pass
@case(True)
def c_identity(o):
    """the correct outputs unchanged"""


@case(True)
def c_scale_gamma(o):
    """every gamma scaled by 1 + 5e-8"""
    for gid, j, e in each(o):
        for dr in e:
            for fam in e[dr]:
                for m in e[dr][fam]:
                    m[0] *= 1 + 5e-8
                    m[1] *= 1 + 5e-8


@case(True)
def c_scale_power(o):
    """every power fraction scaled by 1 + 3e-7"""
    for gid, j, e in each(o):
        for dr in e:
            for fam in e[dr]:
                for m in e[dr][fam]:
                    m[2] = [x * (1 + 3e-7) for x in m[2]]


@case(True)
def c_shuffle(o):
    """modes of every family listed in reverse order, families in reverse order"""
    for gid, j, e in each(o):
        for dr in e:
            e[dr] = {f: list(reversed(e[dr][f])) for f in reversed(list(e[dr]))}


# ---------------------------------------------------------------- search bounds and families
@case(False)
def c_box_material(o):
    """modes with |beta| > 10 k0 sqrt(max eps') + 200 rad/m left out (a search box sized by the materials)"""
    for gid, j, e in each(o):
        g = GBYID[gid]
        mx = max([1.0] + [abs(float(L['material'].get('eps', 1))) for L in g['layers']])
        kk = k0(gid, j)
        for dr in e:
            for fam in list(e[dr]):
                e[dr][fam] = [m for m in e[dr][fam] if abs(m[0]) <= 10 * kk * math.sqrt(mx) + 200]
                if not e[dr][fam]:
                    del e[dr][fam]


@case(False)
def c_box_50k0(o):
    """modes with |beta| > 50 k0 left out (the magnetostatic surface waves of the YIG films)"""
    drop(o, lambda gid, j, dr, fam, m: abs(m[0]) > 50 * k0(gid, j))


@case(False)
def c_last_family(o):
    """the highest family of every spectrum of H1 left out (search in n stopped one family early)"""
    for gid, j, e in each(o, {'H1'}):
        for dr in e:
            hyb = [f for f in e[dr] if f.startswith('HYB')]
            if hyb:
                del e[dr][max(hyb, key=nfam)]


@case(False)
def c_families_upto_8(o):
    """families beyond n = 8 left out (a fixed cap on the family index)"""
    drop(o, lambda gid, j, dr, fam, m: nfam(fam) > 8)


@case(False)
def c_one_way(o):
    """H5 treated as reciprocal: the reverse direction given the one-way forward modes"""
    for gid, j, e in each(o, {'H5'}):
        e['reverse'] = copy.deepcopy(e['forward'])


@case(False)
def c_reciprocal_h2(o):
    """the ground-plane guide H2 treated as reciprocal: reverse copied from forward"""
    for gid, j, e in each(o, {'H2'}):
        e['reverse'] = copy.deepcopy(e['forward'])


@case(False)
def c_bias_sign(o):
    """the sign of the bias ignored in H2 and H5 (forward and reverse swapped)"""
    for gid, j, e in each(o, {'H2', 'H5'}):
        e['forward'], e['reverse'] = e['reverse'], e['forward']


@case(False)
def c_beta_positive(o):
    """beta > 0 for every mode (the backward magnetoplasmons of H5 at 3.5 GHz given positive beta)"""
    for gid, j, e in each(o):
        for dr in e:
            for fam in e[dr]:
                for m in e[dr][fam]:
                    m[0] = abs(m[0])


@case(False)
def c_slow_plasmon(o):
    """the slow surface plasmons of H4 (|beta| > 8 k0) left out"""
    drop(o, lambda gid, j, dr, fam, m: gid == 'H4' and abs(m[0]) > 8 * k0(gid, j))


@case(False)
def c_ceramic_lsm(o):
    """H3: the LSM families left out (a single TE search)"""
    drop(o, lambda gid, j, dr, fam, m: gid == 'H3' and fam.startswith('LSM'))


@case(False)
def c_near_limit(o):
    """the mode with the largest attenuation of every H1 spectrum left out"""
    for gid, j, e in each(o, {'H1'}):
        for dr in e:
            best = max(((fam, k) for fam in e[dr] for k in range(len(e[dr][fam]))), key=lambda t: e[dr][t[0]][t[1]][1])
            del e[dr][best[0]][best[1]]
            if not e[dr][best[0]]:
                del e[dr][best[0]]


# ---------------------------------------------------------------- values
@case(False)
def c_gamma_off(o):
    """one gamma (H1, highest family) off by 3e-7"""
    for gid, j, e in each(o, {'H1'}):
        fam = max((f for f in e['forward'] if f.startswith('HYB')), key=nfam)
        e['forward'][fam][0][0] *= 1 + 3e-7
        break


@case(False)
def c_power_unsigned(o):
    """|P_i| / sum |P_k| (unsigned power fractions)"""
    for gid, j, e in each(o):
        for dr in e:
            for fam in e[dr]:
                for m in e[dr][fam]:
                    s = sum(abs(x) for x in m[2])
                    m[2] = [abs(x) / s for x in m[2]]


@case(False)
def c_power_reverse_from_forward(o):
    """reverse modes of H2 given the power split of the forward mode of the same family"""
    for gid, j, e in each(o, {'H2'}):
        for fam in e['reverse']:
            if fam in e['forward']:
                for m in e['reverse'][fam]:
                    m[2] = list(e['forward'][fam][0][2])


@case(False)
def c_power_off(o):
    """one power fraction off by 3e-6 sum |p|"""
    for gid, j, e in each(o, {'H2'}):
        m = e['forward']['LSE0'][0]
        m[2][1] += 3e-6 * sum(abs(x) for x in m[2])
        break


# ---------------------------------------------------------------- format
@case(False)
def c_missing_guide(o):
    """a guide missing from its output file"""
    del o['run_c']['guides']['H4']


@case(False)
def c_nan(o):
    """a NaN in a power list"""
    for gid, j, e in each(o, {'H6'}):
        e['forward']['LSE0'][0][2][0] = float('nan')
        break


@case(False)
def c_hyb0(o):
    """an HYB0 family"""
    for gid, j, e in each(o, {'H6'}):
        e['forward']['HYB0'] = e['forward'].pop('LSE0')
        break


@case(False)
def c_short_power(o):
    """a power list one number short"""
    for gid, j, e in each(o, {'H1'}):
        e['forward']['LSE0'][0][2].pop()
        break


def main():
    good = 0
    for name, expect, fn in CASES:
        ok, why = verdict(edited(fn))
        good += ok == expect
        print('%-6s %-5s %s%s' % ('ok' if ok == expect else 'WRONG', 'pass' if ok else 'fail', name,
                                  '' if ok else '  [' + why + ']'))
    print('%d of %d cases behave as expected' % (good, len(CASES)))


if __name__ == '__main__':
    main()
