"""Second check of min_inductors: the ladders of alt_synthesize.py (a separately written branch and
bound over the canonical ladders) are graded by the verifier and their inductor counts compared with
the stored minimum.

usage: python3 check_min_inductors.py tests_dir alt_output_dir
(alt_output_dir holds <spec_id>.json written by: python3 alt_synthesize.py spec.json out.json)
"""
import json
import os
import sys

TESTS, OUT = sys.argv[1], sys.argv[2]
sys.path.insert(0, TESTS)
import test_synthesis as T  # noqa: E402

truth = json.load(open(os.path.join(TESTS, 'truth.json')))
same = 0
for sid in sorted(truth):
    spec = json.load(open(os.path.join(TESTS, 'specs', sid + '.json')))
    rec = truth[sid]
    net = T.read_network(os.path.join(OUT, sid + '.json'))
    RL, arms = T.check_contract(net)
    n_ind = sum(1 for _, t, _ in arms if 'L' in T.TYPES[t])
    try:
        T.check_network(spec, rec, RL, arms)
        verdict = 'passes'
    except AssertionError as e:
        verdict = 'FAILS: %s' % str(e)[:80]
    same += n_ind == rec['min_inductors']
    print('%-11s alt %2d inductors, stored minimum %2d, verifier %s' % (sid, n_ind, rec['min_inductors'], verdict))
print('%d of %d specifications: the separate search finds the same minimum' % (same, len(truth)))
