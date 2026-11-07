"""Author-side ablations: string patches of solution/identify.py, one ingredient removed each.
usage: python3 make_ablations.py <identify.py> <out_dir>"""
import os, sys
src = open(sys.argv[1]).read()
out = sys.argv[2]
P = {}
# database topology trusted: no topology search at all
P["database_topology"] = [(
    'def select_topology(w, start_vals):\n',
    'def select_topology(w, start_vals):\n'
    '    db = [c["id"] for c in CANDS if not c["known"] and c["database_status"] == "closed"]\n'
    '    m0 = Model(db, "adm")\n'
    '    m, th, cov, c0, u = converge(db, values(m0, m0.theta0), w, 60)\n'
    '    return db, m, th, cov, u\n')]
# exact repeats only: no other malfunction detection
P["repeats_only"] = [('def detect(raw, stale_mask):\n', 'def detect(raw, stale_mask):\n    return {}\n')]
# superset fit from the nameplate start only
P["nameplate_start_only"] = [(
    '    starts = list(extra) + [("nameplate", values(m0, th0)), ("staged", staged_start(w)[0]),\n'
    '                            ("admittance-first", admittance_first_start(w))]\n',
    '    starts = [("nameplate", values(m0, th0))]\n')]
# Kaufman approximation of the variable-projection Jacobian
P["kaufman"] = [('            J -= Qc @ yv                                               # - A^{+T} dA^T r\n',
                 '            pass                                                       # Kaufman: term dropped\n')]
# no plausibility check (the v2 program)
P["no_plausibility"] = [('    while True:\n        best = None\n        for cid in active:\n            c = CANDS[cid - 1]\n',
                         '    while False:\n        best = None\n        for cid in active:\n            c = CANDS[cid - 1]\n')]
# per-frame whitening, but at each channel's mean angle over the recording (time-invariant)
P["mean_angle_whitening"] = [('        r, a = np.abs(z[:, k]), np.angle(z[:, k])\n',
                              '        r, a = np.abs(z[:, k]), np.full(len(z), np.angle(np.sum(Z[:, k])))\n')]
# isotropic noise: one sigma per channel, the rms of the magnitude and angle sigmas
P["isotropic_noise"] = [('        sm, sa = NOI[c]["magnitude_rel_std"], NOI[c]["angle_std_rad"]\n',
                         '        sm, sa = NOI[c]["magnitude_rel_std"], NOI[c]["angle_std_rad"]\n'
                         '        sm = sa = np.sqrt((sm * sm + sa * sa) / 2)\n')]
for name, edits in P.items():
    s = src
    for a, b in edits:
        assert s.count(a) == 1, (name, a[:50])
        s = s.replace(a, b)
    os.makedirs(f"{out}/{name}", exist_ok=True)
    open(f"{out}/{name}/identify.py", "w").write(s)
    print("wrote", name)
