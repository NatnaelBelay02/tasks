import sys, numpy as np, itertools
sys.path.insert(0, sys.argv[1])
import model as M
d = sys.argv[2]
net = M.load_json(d + "/heldout_network.json"); tr = M.load_json(d + "/heldout_truth.json")
Z = M.read_pmu(d + "/heldout_pmu.csv")[::10]
def run(calmode, tapside, status):
    cal = np.ones(26, complex)
    for k, (c, b) in enumerate(M.CHANNELS):
        if f"{c}{b}" in tr["cal"]:
            e = tr["cal"][f"{c}{b}"]; cal[k] = (1 + e["ratio_error"]) * np.exp(1j * e["phase_error"])
    Zc = Z / cal if calmode else Z
    cands = net["candidates"]
    if tapside:
        cands = [dict(c, from_bus=c["to_bus"], to_bus=c["from_bus"]) if c["model"] == "transformer" else c for c in cands]
    Y = M.ybus(cands, status, tr["params"], net["shunt"])
    V = np.zeros((len(Z), 14), complex)
    for k, (c, b) in enumerate(M.CHANNELS[:13]):
        V[:, b - 1] = Zc[:, k]
    free = [i for i in range(14) if i != M.Z7]
    V[:, M.Z7] = -(V[:, free] @ Y[Z7 := M.Z7, free]) / Y[M.Z7, M.Z7]
    Ip = (V @ Y.T)[:, [b - 1 for b in M.PMU_BUSES]]
    return np.median(np.abs(Zc[:, 13:] - Ip) / np.abs(Ip).clip(1e-3), 0)
st = tr["status"]
print('base', run(1, 0, st).round(4))
print('nocal', run(0, 0, st).round(4))
print('tap to-side', run(1, 1, st).round(4))
allc = {k: True for k in st}
print('all closed', run(1, 0, allc).round(4))
dbst = {str(c["id"]): c["database_status"] == "closed" for c in net["candidates"]}
print('database status', run(1, 0, dbst).round(4))
# predicted residual std from noise propagation (polar noise, magnitudes ~ |z|)
cal = np.ones(26, complex)
for k, (c, b) in enumerate(M.CHANNELS):
    if f"{c}{b}" in tr["cal"]:
        e = tr["cal"][f"{c}{b}"]; cal[k] = (1 + e["ratio_error"]) * np.exp(1j * e["phase_error"])
Zc = Z / cal
Y = M.ybus(net["candidates"], st, tr["params"], net["shunt"])
free = [i for i in range(14) if i != M.Z7]
S = np.zeros((14, 13), complex); S[free, np.arange(13)] = 1; S[M.Z7] = -Y[M.Z7, free] / Y[M.Z7, M.Z7]
YS = (Y @ S)[[b - 1 for b in M.PMU_BUSES]]          # 13 x 13 map from free V to currents
sv = net["noise"]["V"]; si = net["noise"]["I"]
Vabs = np.abs(Zc[:, :13]); Iabs = np.abs(Zc[:, 13:])
# variance of complex noise of a channel ~ |z|^2 (sm^2 + sa^2)
varV = Vabs ** 2 * (sv["magnitude_rel_std"] ** 2 + sv["angle_std_rad"] ** 2)
varI = Iabs ** 2 * (si["magnitude_rel_std"] ** 2 + si["angle_std_rad"] ** 2)
pred = varI + varV @ (np.abs(YS) ** 2).T
V = Zc[:, :13] @ np.eye(13)
Ip = Zc[:, :13] @ YS.T
obs = np.abs(Zc[:, 13:] - Ip) ** 2
print('observed / predicted residual variance per current channel:', (obs.mean(0) / pred.mean(0)).round(3))
bad = np.zeros(9000, bool)
for s in tr["corrupted"]:
    bad[s["first_frame"]:s["last_frame"] + 1] = True
ok = ~bad[::10]
print('clean frames only:', (obs[ok].mean(0) / pred[ok].mean(0)).round(3))
