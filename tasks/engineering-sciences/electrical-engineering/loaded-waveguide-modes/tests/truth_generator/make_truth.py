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
import hybacb as HY
import power as PW
from wgacb import Guide

ctx.prec = 256
A = 10 * arb(10).log()                # 200 dB/m in Np/m
A1 = float((A * arb('1.02')).mid())
MARGIN = arb('0.02')
COND = 30          # guarantee stated in the instruction: sum |P_i| <= 30 |P| for every reported mode


def with_power(G, n, f, lst, sign=1):
    """[beta, alpha, fractions] strings for the modes gamma in lst (fields of the zero sign * gamma)."""
    out, worst = [], 0.0
    for g in lst:
        fr, cond = PW.fractions(G, n, g * sign, f)
        assert cond < COND, (G.id, n, 'power split too ill-conditioned', cond)
        for x in fr:
            assert float(x.rad()) < 1e-9, (G.id, n, 'power fraction not resolved', x)
        worst = max(worst, float(cond.upper()))
        out.append([s25(g.imag), s25(g.real), [x.mid().str(17, radius=False) for x in fr]])
    return out, worst


def hybrid_families(G, f, fs, dg):
    """Certified forward / reverse modes of the hybrid families n = 1, 2, ... of a gyrotropic guide."""
    fwd_all, rev_all = {}, {}
    k0, m = HY.reach(G, fs)
    n, empty = 1, 0
    worst = 0.0
    while True:
        H = HY.HybGuide(G, n)
        ky = float(H.ky.mid())
        B = 4 * k0 * math.sqrt(m) + 2 * ky + 400
        t = time.time()
        zeros, nbox = HY.modes(H, f, A1, B)
        fwd, rev = [], []
        for z in zeros:
            assert z.real > 0 or z.real < 0, (G.id, fs, n, 'Re gamma not certainly nonzero')
            g = z if z.real > 0 else -z
            assert abs(g.real - A) > MARGIN * A, (G.id, fs, n, 'alpha within 2 % of the limit', g)
            if g.real < A:
                (fwd if z.real > 0 else rev).append(g)
        for lst in (fwd, rev):
            for i in range(len(lst)):
                for j in range(i):
                    assert (lst[i] - lst[j]).abs_lower() > arb('1e-4') * lst[i].abs_upper(), (G.id, fs, n, 'close modes')
        if fwd:
            fwd_all['HYB%d' % n], w1 = with_power(G, n, f, fwd, 1)
            worst = max(worst, w1)
        if rev:
            rev_all['HYB%d' % n], w2 = with_power(G, n, f, rev, -1)
            worst = max(worst, w2)
        dg.append({'f': fs, 'kind': 'HYB%d' % n, 'B_over_k0': B / k0, 'zeros': len(zeros), 'boxes': nbox,
                   'seconds': round(time.time() - t, 1)})
        print(G.id, fs, 'HYB%d' % n, 'forward', len(fwd), 'reverse', len(rev), 'zeros', len(zeros), 'boxes', nbox, flush=True)
        # stop after two consecutive families without a zero in the region, past k_y^2 > k0^2 m + A1^2
        if not zeros and ky * ky > k0 * k0 * m + A1 * A1:
            empty += 1
            if empty == 2:
                break
        else:
            empty = 0
        n += 1
    return fwd_all, rev_all, worst


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
                entry = {'forward': {}, 'reverse': {}}
                worst = 0.0
                if fwd:
                    entry['forward']['LSE0'], w = with_power(G, 0, f, fwd, 1)
                    worst = max(worst, w)
                if rev:
                    entry['reverse']['LSE0'], w = with_power(G, 0, f, rev, -1)
                    worst = max(worst, w)
                dg.append({'f': fs, 'kind': 'TE (n = 0)', 'B_over_k0': B / k0, 'zeros': len(zeros), 'boxes': nbox,
                           'seconds': round(time.time() - t, 1)})
                print(rec['id'], fs, 'LSE0 forward', len(fwd), 'reverse', len(rev), flush=True)
                hf, hr, w = hybrid_families(G, f, fs, dg)
                entry['forward'].update(hf)
                entry['reverse'].update(hr)
                dg.append({'f': fs, 'max_power_condition': max(worst, w)})
                modes.append(entry)
                continue
            entry = {}
            for kind in ('LSE', 'LSM'):
                t = time.time()
                Q = CZ.lse_bound(G, f) if kind == 'LSE' else CZ.lsm_bound(G, f, A1)
                zeros, nbox = CZ.modes(G, kind, f, A1, Q)
                fam = families(G, kind, zeros)
                for k, lst in fam.items():
                    entry[k], w = with_power(G, int(k[3:]), f, lst, 1)
                    dg.append({'f': fs, 'family': k, 'max_power_condition': w})
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
