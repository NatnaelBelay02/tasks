"""Fraction of graded quantities that the default answer (nameplate parameters, zero calibration) passes."""
import json, sys
def frac(net, tr, tol):
    nom = {str(c['id']): c['nominal'] for c in net['candidates']}
    n = ok = 0
    for cid, d in tr['params'].items():
        for k, v in d.items():
            n += 1; ok += abs(v - nom[cid][k]) <= tol['branches'][cid][k]
    for ch, d in tr['cal'].items():
        for k, v in d.items():
            n += 1; ok += abs(v) <= tol['channels'][ch][k]
    return ok, n
if __name__ == '__main__':
    net, tr, tol = (json.load(open(p)) for p in sys.argv[1:4])
    print('%d of %d quantities pass with the default answer' % frac(net, tr, tol))
