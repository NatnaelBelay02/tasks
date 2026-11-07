# pmu-grid-identification

Write a program that identifies a 14-bus transmission network and its instrumentation from a 75-minute noisy synchrophasor (PMU) recording (2 frames per second, 9000 frames). It must find which candidate branches are in service, every unknown branch impedance and transformer tap, the ratio and phase error of every non-reference instrument channel, and the PMU malfunction intervals, to Cramér–Rao-level accuracy. The program is graded on four held-out recordings it has never seen.

| | |
|---|---|
| Domain / category | engineering-sciences / electrical-engineering (power-system measurement and identification) |
| Expert time estimate | 24 hours |
| Agent budget | 18000 s |
| Agent container | 4 CPUs, 8 GB RAM; one example recording with its tolerances, no truth |
| Grading container | 4 CPUs, 8 GB RAM, no network; the program runs once per held-out recording, all four at once, 3600 s wall limit each |

## Difficulty

Operators feed PMU data into state estimation, protection and stability studies that assume the network model and the instrument transformers are right. In the field, line impedances drift, nameplates are wrong, taps are not telemetered, the topology database lags switching, and voltage and current transformers carry ratio and phase errors of up to about 1 % and 0.6°. Validating the model and calibrating the instruments from synchrophasor data is routine work for transmission engineers and measurement researchers. It is hard because all of these unknowns enter the same measurements, and a program has to get them right on recordings it cannot inspect:

- **Errors in variables.** Every voltage and current is noisy and no true state is given (9000 × 26 complex nuisance phasors per recording). Regressing measured currents on measured voltages is biased by hundreds of standard deviations. The full errors-in-variables likelihood has to be maximised.
- **Noise that changes every frame.** The instrument noise is Gaussian in magnitude and angle, with the angle noise 2 to 8 times larger. The system frequency drifts, so every phasor rotates by up to about a hundred radians over a recording. In rectangular coordinates the noise covariance therefore changes from frame to frame. Compressing the recording into one sufficient statistic with one covariance per channel, which is enough on a steady, nearly isotropic recording, is no longer efficient (calibration table below).
- **Calibration is entangled with the parameters.** A current-ratio error looks like a scaled admittance, and a voltage-phase error looks like a rotated angle difference. The 44 instrument unknowns are identifiable only jointly with the branch parameters, and only through the two reference-grade PMUs.
- **Topology is not given, and the likelihood has false minima.** The database is wrong for 0 to 2 candidates per recording, in either direction. Fitting every candidate at once has compensating local minima: on rec_a the database lists 6–13 as closed although it is open, and every start of the superset fit keeps it alive with a non-physical impedance (negative resistance), while 6–12 and 12–13 bend to compensate (146σ off on 6–12). Only some starting points reach the global optimum, and which ones differs between recordings; recognising the false minimum needs the physical prior that the instruction states (impedances within 10 % of nameplate).
- **Unflagged, mostly subtle malfunctions.** Each recording has 1 to 4 intervals of stale data, timing errors (constant, ramping or out-and-back phase shifts of a whole PMU, 0.0024 to 0.03 rad) and ratio steps of one channel (0.0008 to 0.0016). After the per-frame states absorb what they can, the constant shifts and ratio steps add only 15 to 58 units of χ² per frame, spread over a PMU's four residual components. The large timing errors (up to 3600 per frame) also distort the residuals of neighbouring PMUs. Each interval must still be attributed to the right PMU and located to within 30 frames at both ends. Left in, they bias the parameters and the topology tests. Finding them needs a model that is already right, so detection, topology and estimation are coupled.
- **Tolerances are four Cramér–Rao standard deviations.** A merely consistent estimator does not pass. It has to be essentially efficient on every one of the four recordings.

The data are synthetic, generated from the public IEEE 14-bus system under realistic operating conditions (details under Verification, "Operating conditions").

**References.** A. G. Phadke, J. S. Thorp, *Synchronized Phasor Measurements and Their Applications*, Springer, 2nd ed. 2017. A. Abur, A. G. Expósito, *Power System State Estimation: Theory and Implementation*, Marcel Dekker 2004. G. H. Golub, V. Pereyra, "Separable nonlinear least squares: the variable projection method and its applications," Inverse Problems 19(2), 2003. L. Kaufman, "A variable projection method for solving separable nonlinear least squares problems," BIT 15, 1975. Q. Zhang, V. Vittal, G. T. Heydt, N. Logic, S. Sturgill, "The integrated calibration of synchronized phasor measurement data in power transmission systems," IEEE Trans. Power Delivery 26(4), 2011. E. S. Page, "Continuous inspection schemes," Biometrika 41, 1954. R. D. Zimmerman, C. E. Murillo-Sánchez, R. J. Thomas, "MATPOWER," IEEE Trans. Power Systems 26(1), 2011. F. Milano, R. Zárate-Miñano, "A systematic method to model power systems as stochastic differential algebraic equations," IEEE Trans. Power Systems 28(4), 2013. E. Hirst, B. Kirby, "Separating and measuring the regulation and load-following ancillary services," Utilities Policy 8(2), 1999. G. C. Montanari, M. Loggini, A. Cavallini, L. Pitti, D. Zaninelli, "Arc-furnace model for the study of flicker compensation in electrical networks," IEEE Trans. Power Delivery 9(4), 1994.

## Reference solution

`solution/solve.sh` copies `solution/identify.py` to `/root/results/identify.py`. The program uses numpy only and builds its own network model from `network.json`. For fixed parameters θ, the true phasors of each frame enter linearly: the 13 complex free voltages, with the zero-injection bus 7 eliminated exactly, map through a 52 × 26 real matrix A(θ) to the 26 complex channels. Each frame's states are projected out by a batched QR (variable projection). The polar noise is whitened per frame and channel by its 2 × 2 inverse square-root covariance, evaluated at the measured phasor. The projected cost is minimised over θ by Levenberg–Marquardt with the exact Golub–Pereyra Jacobian; with the Kaufman approximation one start stalls and the program fails one recording (calibration table). Series admittances are kept within 1.5 times nameplate magnitude and charging multipliers within [0, 1.5] during the search.

The pipeline:

1. **Stale data:** exact repeats of the previous frame, per channel.
2. **Superset model:** every non-known candidate in series-admittance coordinates (g, b_s, plus a charging multiplier for lines). It is converged from several starts, and the lowest cost is kept: nameplate values with database-open corridors at zero; calibrations fitted first on the database topology; admittances fitted first with calibrations frozen at zero; and, later, the previous robust estimate. Search fits use every fourth frame.
3. **Robust fit:** Cauchy reweighting per frame and PMU, scaled by the median residual statistic.
4. **Malfunctions:** runs of the rolling-median residual statistic per PMU and per channel are refined by fitting step, ramp and out-and-back shapes to per-frame estimates of a PMU phase offset and of channel ratio offsets. An interval is accepted if its χ² gain exceeds 200 and it overlaps no stronger one. The intervals are masked.
5. **Topology:** on the masked data, the superset is re-converged from all starts, then backward elimination on the penalised likelihood (cost + 30 per corridor). Then a plausibility check: a corridor whose fitted impedance has a non-positive resistance or reactance, or lies more than 50 % from nameplate (the instruction bounds every unknown to 10 %), is a compensating artefact of a false minimum; it is removed and the rest refitted from the current and from fresh nameplate starts, and the removal is kept if the penalised likelihood improves. Steps 3–5 repeat, at most three times, until the intervals stop changing.
6. **Final fit:** a maximum-likelihood fit in physical coordinates (r, x; x, tap; ratio, phase) on all frames, with the intervals masked, gives the reported values.

Measured in the grading container with nothing else running: all four held-out runs finished within 1100 s in the final harbor oracle run (within 1421 s in an earlier run), against the 3600 s limit. Each run peaked below 0.95 GB of memory (cgroup `memory.peak`, five runs at 1 CPU each).

## Verification

**How the program is run.** The grading container holds Python 3.12 with numpy 1.26.4 and scipy 1.13.1, the four held-out recordings and their truth (in `/tests/data/rec_*`, stored as `heldout_pmu.csv`, `heldout_network.json`, `heldout_tolerances.json` and `heldout_truth.json`; `/tests` is mode 700), and the submitted program. The verifier opens `/root/results/identify.py` without following symlinks and without blocking (so a FIFO is rejected at once), requires a regular file and copies its content. For each recording, in shuffled order, it stages a fresh run directory: a read-only input directory with a random name holding only `pmu.csv` and `network.json` (copied from the `heldout_` files), and a working directory (mode 1777) holding only the copy of the program. It first checks, as the sandbox user, that the working directory and `/tmp` are writable, the input directory is not, and `/tests` is unreadable. It then starts all four runs at once as uid/gid 65534 in their own sessions, with the environment stated in the instruction and stdout/stderr going to unlinked files. Each run must exit with status 0 within 3600 s; a run past the limit has its process group killed. Afterwards every process of uid 65534 is killed: a helper that has become uid 65534 calls `kill(-1, SIGKILL)`, repeated until `/proc` shows none, so processes that keep forking cannot outlive the run. The output is then opened the same way as the program and must be a regular file, not a symlink, of at most 1 MB; at most 1 MB + 1 byte is ever read.

**Checks.** The tests map one-to-one to the numbered requirements of the instruction, each parametrised over the four recordings (24 tests), and all must pass:

1. **Run:** exit status 0 within the limit, output present as stated.
2. **Schema:** one entry per candidate with an integer `id` and a boolean `in_service`; parameters of in-service, non-known entries present as JSON numbers; exactly the 22 channels; segments of integers with 0 ≤ first ≤ last ≤ frames − 1. While parsing, every number literal is compared exactly against the largest finite double, and repeated keys are rejected. Non-integer literals are kept as exact decimals.
3. **Topology:** `in_service` equals the true status of all 24 candidates.
4. **Branch parameters:** every unknown r, x, x and tap of the truly in-service, non-known candidates is within its tolerance.
5. **Calibration:** all 44 ratio and phase errors are within their tolerances, with no angle wrapping.
6. **Malfunctions:** as many segments as true intervals, paired one-to-one on the same bus with both endpoints within 30 frames.

Requirements 4 and 5 compare the exact decimal value written with the stored double, exactly and inclusive at the bound.

**Where the data come from.** The author-side generator (`gen_recording.py`, with `model.py`, kept outside the package) starts from the IEEE 14-bus system (MATPOWER `case14`: bus loads, generator and voltage set points). Every frame is an exact AC power flow (Newton–Raphson to 1e-13 pu) of that frame's loads and generation. Then:
- all phasors of the frame are rotated by the accumulated frequency deviation;
- the instrument errors and the polar noise are applied;
- the malfunctions are injected;
- angles are reported as principal values.

The rebuilt model was checked against the v2 recordings before it was used: with the planted truth, the measured currents of rec_a agree with Y·V to within the noise (variance ratio 0.95–1.16 per channel on the clean frames).

| item | value (all recordings) |
|---|---|
| network | bus 7 is the explicit star point of the three-winding transformer (4–7, 7–8, 7–9 known); extra corridors 1–4 and 3–5 (lines), 10–14 and 11–12 (cables) |
| planted parameters | line and cable r, x within ±10 % of nameplate; transformer x within ±5 %; taps nameplate ± up to three 0.625 % steps; new draw per recording |
| instrument errors | ratio uniform in ±1 %, phase uniform in ±10.5 mrad; buses 1 and 6 exact; new draw per recording |
| recording | 9000 frames at 2 frames per second (75 min), a decimated archive of the PMU stream |
| operating conditions | see the next table |
| frequency deviation | Ornstein–Uhlenbeck with τ = 30 s; its integral rotates every phasor of a frame |

**Operating conditions** (v3; the v2 review found the earlier load model physically unfaithful). v2 drove the identification with loads that no bulk system shows: an independent 5 % per-bus Ornstein–Uhlenbeck component with τ = 1.5 s, a 3 % common component with τ = 20 s, and system load profiles from 0.5 to 1.1 and back within 15 minutes (4–8 %/min). v3 replaces them with the mechanisms that actually move a transmission system's operating point, each at a realistic rate, and lengthens the recording to 75 minutes, so that ordinary ramps supply the excitation:

| mechanism | model | resulting rate | basis |
|---|---|---|---|
| system load level | steady; ramp up (morning pickup) or down (evening decline); or rise and fall (a peak passing) | at most 0.3 %/min | bulk-system ramps are sub-percent per minute |
| aggregate demand of each load bus | OU on P and Q, σ = 2 %, τ = 300 s; OU on Q only, σ = 1 %, τ = 300 s (power-factor drift); fast component σ = 0.2 %, τ = 0.5 s (many small switching loads) | demand of all buses except the steel works: 1-min change 0.50–0.62 % rms, net ramp ≤ 0.33 %/min | minute-scale "regulation" variability of aggregated load [Hirst & Kirby]; OU load models [Milano & Zárate-Miñano] |
| steel works at bus 3 | 30 MW / 24 MVAr AC electric arc furnace inside the 94.2 MW bus load: P fluctuates 8 % and Q 15 % of rating from frame to frame while the arc burns (OU, τ = 0.2 s); a 15-min power-off for back-charging in rec_a, rec_c and rec_d | bus 3 between 56 and 105 MW | the random, frame-to-frame P/Q fluctuation of arc furnaces [Montanari et al.]; power-off periods between heats are routine |
| generator at bus 2 | 5-minute real-time dispatch: a new set point every 300 s (random step, σ = 4 MW, kept within 20–100 MW), reached at 8 MW/min; plus 25 % of the load deviation through AGC (first order, τ = 60 s); bus 1 is the slack | 7–97 MW over a recording | 5-minute security-constrained dispatch of real-time markets; AGC participation |
| voltage control | buses 1, 2, 3, 6, 8 hold their set points (1.06, 1.045, 1.01, 1.07, 1.09 pu) | – | case14 |
| frequency | OU, τ = 30 s, standard deviation 0, 20 or 30 mHz (largest deviation 0.1 Hz) | angles of a frame rotate by ≤ 0.62 rad/s, 58–99 rad over a recording | ordinary interconnection frequency wander |

The only fast, large variation is at the steel works, a named industrial load whose arc physically produces it. All other demand moves at sub-percent per minute.

| recording | system load level | arc furnace off | freq. deviation std | noise σ_m / σ_a, V; I | database errors | malfunctions (bus, frames, kind, size) |
|---|---|---|---|---|---|---|
| example (`/root/data`) | steady 1.00 | – | 0 | 2e-4 / 3e-4; 5e-4 / 8e-4 | 6–13 listed closed, open; 10–14 listed open, closed | 10: 2400–2439, V stale |
| rec_a | ramp 0.775 → 1.0 | 2600–3500 s | 20 mHz (99 rad span) | 1e-4 / 5e-4; 2.5e-4 / 1.25e-3 | same two as the example | 12: 3100–3499 out-and-back timing error, peak 0.03 rad; 9: 5200–5899 phase shift 3e-3 rad; 4: 7000–7059 I stale |
| rec_b | 0.8875 → 1.0 → 0.8875 | – | 30 mHz (92 rad) | 1e-4 / 8e-4; 3e-4 / 2.4e-3 | 1–4 listed open, closed; 10–11 listed closed, open | 3: 800–1399 V ratio +8e-4; 12: 4500–4799 timing ramp to 0.02 rad, then jump back; 10: 6000–6039 V stale |
| rec_c | ramp 1.0 → 0.775 | 750–1650 s | 0 | 1.5e-4 / 7.5e-4; 4e-4 / 2e-3 | none | 13: 3000–4199 phase shift 4.5e-3 rad |
| rec_d | steady 1.00 | 3050–3950 s | 20 mHz (58 rad) | 2e-4 / 4e-4; 5e-4 / 1e-3 | 3–5 and 11–12 listed open, closed | 2: 1000–1089 I stale; 14: 2500–3299 phase shift 2.4e-3 rad; 11: 4000–4599 V ratio +1.6e-3; 6: 6500–6899 out-and-back timing error, peak 0.025 rad |

The planted parameters, instrument errors, topology errors, noise levels and malfunctions of the four held-out recordings are those of v2 (`heldout_truth.json` is unchanged); only the operating conditions, the frame rate and therefore the recordings and the tolerances are new. The example's truth is a new draw under the same rules (it is not shipped).

The example shows the format only. It has no frequency drift, nearly isotropic noise, a steady load level and one stale interval. Every held-out property is stated in the instruction with its range, and the instruction's ranges were checked against the four recordings as staged. Malfunction sizes were set so that every interval is detectable. The part of the injected deviation that the per-frame states cannot absorb, in units of the noise (χ²), is at least 17 per frame (the smallest ratio step, rec_b) and at least 9.1 × 10³ summed over the interval (author-side `detectability.py`).

**Tolerances.** Each recording's tolerance file (`tolerances.json` for the example, `tests/data/rec_*/heldout_tolerances.json` for the held-out recordings) holds 4σ, where σ is the Cramér–Rao standard deviation. It is computed from the exact Fisher information of the measurement model at that recording's planted truth, over its frames free of malfunctions. The true states of every frame are treated as nuisance parameters and removed by a Schur complement. There are 78 unknowns (82 on rec_d); per recording, the Fisher matrix has smallest eigenvalue 1.7e3 to 5.4e3 and condition number 1.7e8 to 6.0e8, so every unknown is identifiable. On the held-out recordings the tolerances range from 2.0e-4 to 6.7e-2 pu for branch parameters, 1.4e-4 to 6.6e-2 for ratio errors and 2.9e-4 to 8.6e-2 rad for phase errors. Out-of-service candidates also get tolerances, so the example's file does not reveal the topology. These come from a counterfactual bound with the candidate in service at nameplate and are never used for grading. The bound is computed by an author-side pipeline (generator, model and Fisher information, `fisher.py`) that shares no code with the reference solution.

**Why the recording is 75 minutes long.** With realistic rates, a 15-minute window does not move the operating point enough: at 10 frames per second over 15 minutes the bound on the calibrations grew by a median factor of 4–6 against v2, and the "default answer" (nameplate parameters, zero calibration) fell within tolerance for 46–60 of the ~80 graded quantities. Over 75 minutes ordinary ramps, dispatch steps and the furnace cycle provide the excitation: branch tolerances are back to the v2 level (median ratio to v2 0.86–1.14), calibration tolerances are 2.1–3.3 times the v2 ones (median), and the default answer passes 31, 47, 45 and 55 of the 78–82 quantities on rec_a–rec_d (v2: 30, 40, 44, 48), so it still fails every recording on the remaining ones, besides topology and malfunctions. The frame count, and with it every frame-indexed quantity of the contract (malfunction lengths, the 30-frame endpoint tolerance, the run time), is unchanged.

**Monte Carlo calibration.** The reference program was re-run on recordings with the same planted truth and malfunctions and new noise, load and frequency realisations (two per held-out recording). Each was scored against its own tolerances:

| recording | noise seed | topology | malfunctions matched (max endpoint error) | worst \|error\| / σ | result |
|---|---|---|---|---|---|
| rec_a | task data | correct | 3/3 (0) | 2.2 | pass |
| rec_a | 2111 | correct | 3/3 (0) | 2.2 | pass |
| rec_a | 2112 | correct | 3/3 (0) | 2.6 | pass |
| rec_b | task data | correct | 3/3 (2) | 1.9 | pass |
| rec_b | 2211 | correct | 3/3 (2) | 2.0 | pass |
| rec_b | 2212 | correct | 3/3 (2) | 2.9 | pass |
| rec_c | task data | correct | 1/1 (0) | 2.3 | pass |
| rec_c | 2311 | correct | 1/1 (0) | 1.7 | pass |
| rec_c | 2312 | correct | 1/1 (0) | 2.5 | pass |
| rec_d | task data | correct | 4/4 (1) | 2.2 | pass |
| rec_d | 2411 | correct | 4/4 (0) | 1.7 | pass |
| rec_d | 2412 | correct | 4/4 (0) | 2.1 | pass |
| example | task data | correct | 1/1 (0) | 2.4 | pass |

All 13 runs pass, with the worst error at 2.9σ against the 4σ limit.

**Independent implementation.** A second program (author-side `independent_v2.py`, about 470 lines, numpy only) was written without reusing code from the reference program. It differs at every stage:

| stage | reference program | independent program |
|---|---|---|
| residual | rectangular, whitened per frame by the 2 × 2 polar covariance | exact polar log residual, log(z / h), magnitude and angle scaled by their σ |
| per-frame states | projected out (variable projection, exact Jacobian) | estimated jointly with the parameters; Levenberg–Marquardt with the state blocks removed by a Schur complement |
| start | multi-start superset (nameplate, staged, admittance-first) | one equation-error start |
| outliers | Cauchy reweighting per frame and PMU | χ² trimming of frame × PMU blocks |
| malfunctions | step, ramp and out-and-back shape fits to per-frame offsets | one-sided CUSUM (Page) on each PMU's 4-row χ² block, refitted up to four rounds |
| topology | penalised backward elimination (cost + 30 per corridor) | Wald-ordered backward elimination with a likelihood-ratio cut of 30 |

On all five recordings it finds the true topology and every malfunction interval (endpoints within 16 frames). Its worst errors are 2.4σ (rec_a), 2.2σ (rec_b), 2.9σ (rec_c), 2.3σ (rec_d) and 2.6σ (example). Its four held-out outputs, replayed through the verifier image, score reward 1 (24/24). Submitted itself as `/root/results/identify.py` to the final verifier image (`tests/test.sh` running all 24 tests of `tests/test_identification.py` on the four held-out recordings), it also scores reward 1 (24/24); all four runs finished within 1592 s. The two programs report the same topology on every recording. Their reported values differ by at most 0.91σ, where σ is the CRLB value of that quantity (rec_d, branch 11 x). Two different efficient estimators therefore land within a fraction of the tolerance of each other and of the planted truth.

**Probes of the verifier.** Small programs that write crafted outputs, built author-side from the planted truth, were run through the verifier image:

| probe | reward |
|---|---|
| planted truth; one value exactly at true + tolerance (exact decimal); all segment endpoints 30 frames off | 1 |
| that value plus 1e-60; one endpoint 31 frames off; segment on the wrong bus; extra or missing segment | 0 |
| `1e309`, or an integer above the largest double, in an unused field; `NaN` | 0 |
| a repeated key whose first copy is wrong and last copy correct (so a parser that keeps the last copy would read the planted truth): `"branches":null` before the real `"branches"`; repeated `in_service` in a branch entry; repeated `ratio_error` in a channel entry; repeated extra top-level key | 0 |
| `id` written `1.0`; frame written with an exponent; inverted segment | 0 |
| extra channel `V1`; boolean in place of a number; phase error off by 2π; an open branch reported in service | 0 |
| extra keys everywhere; parameters on out-of-service entries | 1 |
| output written as a symlink to the verifier's truth file | 0 (not a regular file) |
| reading `/tests`; writing `/root/results` | 0 (permission denied) |
| correct output, then exit status 1 | 0 |
| correct output written only after the wall limit (verifier copy with a 20 s limit) | 0 (killed) |
| wrong output, then a background child that writes the correct output 5 s later | 0 (child killed) |
| wrong output, then a chain of processes that fork a new process every 5 ms and write the correct output after 3 s | 0 (all killed) |
| output left as a FIFO | 0 (not a regular file; the verifier does not block) |
| output of 200 000 nested brackets | 0 (rejected as invalid JSON) |

**Harbor runs on the final package** (harbor 0.24.0, nothing else running on the host):

| agent | reward | detail |
|---|---|---|
| `oracle` | 1.0 | 24/24 tests; all four held-out runs finished within 1100 s; total job 19 min |
| `nop` | 0.0 | no `/root/results/identify.py`; all 24 tests fail; no exceptions |

**Calibration: what weaker methods score.** Each method below was run author-side on the four held-out recordings. Its four outputs were then graded by the verifier image through a replay program that writes the stored output for each recording. Cells give the worst |error|/σ over branch parameters and calibrations; "topology" means a wrong `in_service`, and "segments" means requirement 6 fails. Rows 1–2 rebuild two complete pipelines that solved an earlier version of this task with a steady recording. Both are given the true topology and the true malfunction intervals, padded by 10 frames, which favours them. Rows 3–8 are the reference program with one ingredient removed (author-side `make_ablations.py`, string patches of `identify.py`).

| method | rec_a | rec_b | rec_c | rec_d | reward |
|---|---|---|---|---|---|
| 1. rectangular errors-in-variables fit, one noise covariance per channel for the whole recording | 6.3 | 4.9 | 3.4 | 2.7 | 0 |
| 2. one 26 × 26 complex scatter matrix of the recording, one isotropic σ per channel, subspace fit | 3.3 | 6.2 | 4.6 | 2.1 | 0 |
| 3. per-frame whitening, but at each channel's mean angle (time-invariant) | 4.8, topology, segments | 94.6, topology, segments | 2.7 | 2.6 | 0 |
| 4. isotropic noise (one σ per channel, rms of magnitude and angle σ) | 3.7 | 6.1, segments | 5.0 | 2.4 | 0 |
| 5. Kaufman approximation of the variable-projection Jacobian | 33.5, topology | 1.9 | 2.3 | 2.2 | 0 |
| 6. superset fit from the nameplate start only | 2.2 | 5.6, topology | 2.3 | 2.2 | 0 |
| 7. database topology trusted | 83.5, topology, segments | 112.9, topology, segments | 2.3 | 31.4, topology, segments | 0 |
| 8. exact repeats only, no other malfunction detection | 75.4, topology, segments | 17.5, segments | 12.2, segments | 71.8, topology, segments | 0 |
| reference program | 2.2 | 1.9 | 2.3 | 2.2 | 1 |

Each ingredient of the reference program is needed on at least one held-out recording. With the Kaufman approximation, the staged start on rec_a stalls at a cost of 4606 per frame, against 76 with the exact Jacobian. The search then keeps an open corridor. The verifier's failed tests agree with this table for every row.

## Version history

- **v3 (this version).** Rework after review: the v2 load model was not physically faithful (fast per-bus and common load fluctuations and 4–8 %/min system ramps with no stated mechanism). The operating conditions are now the mechanisms of the "Operating conditions" table, each at a realistic rate. The recording is 75 minutes at 2 frames per second, so the frame count and the contract are unchanged. The recordings, the tolerances and the example's truth were regenerated; the planted truth of the held-out recordings, the topology errors, the malfunctions and the noise levels are those of v2. The verifier code is unchanged. The reference program gained the plausibility check of pipeline step 5. On the new rec_a, every start of the superset fit keeps the open corridor 6–13 alive with a negative resistance, and without the check the program fails. The independent program was rewritten for the new data.
- **v2.** 15-minute recording at 10 frames per second, with the load model the review found unfaithful.
