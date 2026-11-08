import sys, json, time, collections
sys.path.insert(0, sys.argv[1])
import synthesize as SY
def nL(path): return sum(1 for e in path if 'L' in e[1])
for p in sys.argv[2:]:
    spec = json.load(open(p)); t = time.time()
    SY.ADMISSIBLE_ONLY = True
    ap, sols = SY.synthesize(spec)
    cnt = collections.Counter(nL(pth) for _, pth, R in sols)
    first = {sg: nL(pth) for sg, pth, R in reversed(sols)}
    print(spec['id'], 'adm', len(sols), 'Lcounts', dict(sorted(cnt.items())), 'first(+F)', nL(sols[0][1]) if sols else None, '%.0fs' % (time.time() - t), flush=True)
