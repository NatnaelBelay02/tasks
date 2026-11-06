"""Certified reference values (truth.json) from guides.json, independent of solution/.

usage: python3 make_truth.py guides.json truth.json [diagnostics.json] [guide_id ...]

For every guide, test frequency and kind (LSE, LSM), certify2.modes() encloses every zero q of
D_kind(., f) with Re sqrt(-q) <= 1.02 A (A = 200 dB/m = 10 ln 10 Np/m) in disjoint boxes, each
holding exactly one simple zero (complex Krawczyk test), after excluding the rest of the region
with ball arithmetic.  The region ends at Re q = Q: for LSE the Rayleigh bound
Re q <= k0^2 max Re eps, for LSM a proof that D_LSM has no zero beyond Q (certify2.tail_free).
Each zero is refined to a 1e-55 box.  For every family n, gamma = sqrt(k_n^2 - q) with
Re gamma > 0 is enclosed; the mode is reported if alpha = Re gamma <= A.  The generator stops if
an alpha lies within 2 % of A, if Re gamma is not certainly positive, or if two modes of a
family are closer than 1e-4 |gamma| (the guarantees stated in the instruction).
"""
import copy
import json
import math
import sys
import time

from flint import acb, arb, ctx

import certify2 as CZ
from wgacb import Guide

ctx.prec = 256
A = 10 * arb(10).log()                # 200 dB/m in Np/m
A1 = float((A * arb('1.02')).mid())
MARGIN = arb('0.02')


def s25(x):
    return x.mid().str(25, radius=False)


def families(G, kind, zeros):
    n = 0 if kind == 'LSE' else 1
    qmax = max([float(z.real.upper()) for z in zeros] + [0.0])
    out = {}
    while True:
        kn = arb.pi() * n / G.b
        kn2 = kn * kn
        if float((kn2 - arb(A1) ** 2).lower()) > qmax:
            break
        lst = []
        for z in zeros:
            g = (kn2 - z).sqrt()
            if g.real < 0:
                g = -g
            assert g.real > 0, (G.id, kind, n, 'Re gamma not certainly positive', g)
            al = g.real
            assert abs(al - A) > MARGIN * A, (G.id, kind, n, 'alpha within 2 % of the limit', al)
            if al < A:
                lst.append(g)
        for i in range(len(lst)):
            for j in range(i):
                d = (lst[i] - lst[j]).abs_lower()
                assert d > arb('1e-4') * lst[i].abs_upper(), (G.id, kind, n, 'modes closer than 1e-4 |gamma|')
        if lst:
            out['%s%d' % (kind, n)] = lst
        n += 1
    return out


def main():
    src, dst = sys.argv[1], sys.argv[2]
    diag_path = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3].endswith('.json') else None
    only = [x for x in sys.argv[3:] if not x.endswith('.json')]
    guides = json.load(open(src), parse_float=str, parse_int=str)['guides']
    truth, diag = {}, {}
    for rec in guides:
        if only and rec['id'] not in only:
            continue
        G = Guide(rec)
        modes, dg = [], []
        for fs in rec['test_frequencies']:
            f = G.tests[len(modes)]
            k0 = 2 * math.pi * float(fs) / 299792458
            if G.ferrite:
                # y-independent TE modes; D_TE depends on gamma itself (nonreciprocal)
                t = time.time()
                B = CZ.te_bound(G, f, A1)
                zeros, nbox = CZ.modes_te(G, f, A1, B)
                fwd, rev = [], []
                for z in zeros:
                    assert z.real > 0 or z.real < 0, (G.id, fs, 'Re gamma not certainly nonzero')
                    g = z if z.real > 0 else -z
                    assert abs(g.real - A) > MARGIN * A, (G.id, fs, 'alpha within 2 % of the limit', g)
                    if g.real < A:
                        (fwd if z.real > 0 else rev).append(g)
                for lst in (fwd, rev):
                    for i in range(len(lst)):
                        for j in range(i):
                            assert (lst[i] - lst[j]).abs_lower() > arb('1e-4') * lst[i].abs_upper(), (G.id, fs, 'close modes')
                entry = {'forward': {'LSE0': [[s25(g.imag), s25(g.real)] for g in fwd]} if fwd else {},
                         'reverse': {'LSE0': [[s25(g.imag), s25(g.real)] for g in rev]} if rev else {}}
                dg.append({'f': fs, 'kind': 'TE (ferrite)', 'B_over_k0': B / k0, 'zeros': len(zeros), 'boxes': nbox,
                           'seconds': round(time.time() - t, 1)})
                modes.append(entry)
                print(rec['id'], fs, 'forward', len(fwd), 'reverse', len(rev), flush=True)
                continue
            entry = {}
            for kind in ('LSE', 'LSM'):
                t = time.time()
                Q = CZ.lse_bound(G, f) if kind == 'LSE' else CZ.lsm_bound(G, f, A1)
                zeros, nbox = CZ.modes(G, kind, f, A1, Q)
                fam = families(G, kind, zeros)
                for k, lst in fam.items():
                    entry[k] = [[s25(g.imag), s25(g.real)] for g in lst]
                dg.append({'f': fs, 'kind': kind, 'Q_over_k0sq': Q / k0 ** 2, 'zeros': len(zeros), 'boxes': nbox,
                           'seconds': round(time.time() - t, 1)})
            entry = dict(sorted(entry.items(), key=lambda kv: (kv[0][:3], int(kv[0][3:]))))
            # reciprocal guide: the mode set towards -z is the same
            modes.append({'forward': entry, 'reverse': copy.deepcopy(entry)})
            print(rec['id'], fs, {k: len(v) for k, v in entry.items()}, flush=True)
        truth[rec['id']] = {'modes': modes}
        diag[rec['id']] = dg
    with open(dst, 'w') as fh:
        json.dump(truth, fh, indent=1)
    if diag_path:
        with open(diag_path, 'w') as fh:
            json.dump(diag, fh, indent=1)


if __name__ == '__main__':
    main()
