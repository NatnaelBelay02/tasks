"""Reference insertion loss of the example specifications, for the solver to check against.

usage: python3 make_example_losses.py examples_dir truth_examples.json out.json
The values are -10 log10 G of the reference characteristic function (see make_truth.py) at
60 frequencies per example: 30 in the passband, 30 in the stopbands.
"""
import json
import os
import sys

import mpmath as mp

mp.mp.dps = 40


def gain(rec, f):
    u = (mp.mpf(f) / mp.mpf(rec['f2_hz'])) ** 2
    k2 = mp.mpf(rec['c2']) * u ** rec['r'] / u ** rec['dc']
    for q in rec['reflection_zeros_u']:
        k2 *= (u - mp.mpf(q)) ** 2
    for z in rec['transmission_zeros_u']:
        k2 /= (u - mp.mpf(z)) ** 2
    return 1 / (1 + k2)


def main():
    exdir, truth, out = sys.argv[1:4]
    truth = json.load(open(truth))
    res = {}
    for sid, rec in sorted(truth.items()):
        spec = json.load(open(os.path.join(exdir, sid + '.json')))
        f1, f2 = spec['passband_hz']
        fs = [f1 + (f2 - f1) * (k + 0.5) / 30 for k in range(30)]
        if f1 > 0:
            fs += [f1 * 10 ** (-(k + 1) / 15) for k in range(15)] + [f2 * 10 ** ((k + 1) / 15) for k in range(15)]
        else:
            fs += [f2 * 10 ** ((k + 1) / 20) for k in range(30)]
        res[sid] = [[float('%.12g' % f), float(mp.nstr(-10 * mp.log10(gain(rec, float('%.12g' % f))), 15))] for f in sorted(fs)]
    json.dump(res, open(out, 'w'), indent=0)


if __name__ == '__main__':
    main()
