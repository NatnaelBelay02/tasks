# pmu-grid-identification

Write a program that identifies a 14-bus transmission network and its instrumentation from a 75-minute noisy synchrophasor (PMU) recording (2 frames per second, 9000 frames). It must find which candidate branches are in service, every unknown branch impedance and transformer tap, the ratio and phase error of every non-reference instrument channel, and the PMU malfunction intervals, to Cramér–Rao-level accuracy. The program is graded on four held-out recordings it has never seen.

| | |
|---|---|
| Domain / category | engineering-sciences / electrical-engineering (power-system measurement and identification) |
| Expert time estimate | 24 hours |
| Agent budget | 18000 s |
| Agent container | 4 CPUs, 8 GB RAM; one example recording with its tolerances, no truth |
| Grading container | 4 CPUs, 8 GB RAM, no network; the program runs once per held-out recording, all four at once, 5400 s wall limit each |

## Difficulty

Operators feed PMU data into state estimation, protection and stability studies that assume the network model and the instrument transformers are right. In the field, line impedances drift, nameplates are wrong, taps are not telemetered, the topology database lags switching, and voltage and current transformers carry ratio and phase errors of up to about 1 % and 0.6°. Validating the model and calibrating the instruments from synchrophasor data is routine work for transmission engineers and measurement researchers. It is hard because all of these unknowns enter the same measurements, and a program has to get them right on recordings it cannot inspect:

- **Errors in variables.** Every voltage and current is noisy and no true state is given (9000 × 26 complex nuisance phasors per recording). Regressing measured currents on measured voltages is biased by hundreds of standard deviations. The full errors-in-variables likelihood has to be maximised.
- **Noise that changes every frame.** The instrument noise is Gaussian in magnitude and angle, with the angle noise 2 to 8 times larger. The system frequency drifts, so every phasor rotates by up to about a hundred radians over a recording. In rectangular coordinates the noise covariance therefore changes from frame to frame. Whitening each channel with one fixed covariance, or with one isotropic σ, which is enough on a steady, nearly isotropic recording, is no longer efficient (calibration table below).
- **Calibration is entangled with the parameters.** A current-ratio error looks like a scaled admittance, and a voltage-phase error looks like a rotated angle difference. The 44 instrument unknowns are identifiable only jointly with the branch parameters, and only through the two reference-grade PMUs.
- **Topology is not given, and the likelihood has false minima.** The database is wrong for 0 to 2 candidates per recording, in either direction. Fitting every candidate at once has compensating local minima: on rec_a the database lists 6–13 as closed although it is open, and every start of the superset fit keeps it alive with a non-physical impedance (negative resistance), while 6–12 and 12–13 bend to compensate (146σ off on 6–12). Only some starting points reach the global optimum, and which ones differs between recordings; recognising the false minimum needs the physical prior that the instruction states (impedances within 10 % of nameplate).
- **Unflagged, mostly subtle malfunctions.** Each recording has 1 to 4 intervals of stale data, timing errors (constant, ramping or out-and-back phase shifts of a whole PMU, 0.0024 to 0.03 rad) and ratio steps of one channel (0.0008 to 0.0016). After the per-frame states absorb what they can, the constant shifts and ratio steps add only 15 to 58 units of χ² per frame, spread over a PMU's four residual components. The large timing errors (up to 3600 per frame) also distort the residuals of neighbouring PMUs. Each interval must still be attributed to the right PMU and located to within 30 frames at both ends. Left in, they bias the parameters and the topology tests. Finding them needs a model that is already right, so detection, topology and estimation are coupled.
- **Tolerances are four Cramér–Rao standard deviations.** A merely consistent estimator does not pass. It has to be essentially efficient on every one of the four recordings.

The data are synthetic, generated from the public IEEE 14-bus system under realistic operating conditions (details under Verification, "Operating conditions").

**References.** A. G. Phadke, J. S. Thorp, *Synchronized Phasor Measurements and Their Applications*, Springer, 2nd ed. 2017. A. Abur, A. G. Expósito, *Power System State Estimation: Theory and Implementation*, Marcel Dekker 2004. G. H. Golub, V. Pereyra, "Separable nonlinear least squares: the variable projection method and its applications," Inverse Problems 19(2), 2003. L. Kaufman, "A variable projection method for solving separable nonlinear least squares problems," BIT 15, 1975. Q. Zhang, V. Vittal, G. T. Heydt, N. Logic, S. Sturgill, "The integrated calibration of synchronized phasor measurement data in power transmission systems," IEEE Trans. Power Delivery 26(4), 2011. E. S. Page, "Continuous inspection schemes," Biometrika 41, 1954. R. D. Zimmerman, C. E. Murillo-Sánchez, R. J. Thomas, "MATPOWER," IEEE Trans. Power Systems 26(1), 2011. F. Milano, R. Zárate-Miñano, "A systematic method to model power systems as stochastic differential algebraic equations," IEEE Trans. Power Systems 28(4), 2013. E. Hirst, B. Kirby, "Separating and measuring the regulation and load-following ancillary services," Utilities Policy 8(2), 1999. G. C. Montanari, M. Loggini, A. Cavallini, L. Pitti, D. Zaninelli, "Arc-furnace model for the study of flicker compensation in electrical networks," IEEE Trans. Power Delivery 9(4), 1994.

## Reference solution

`solution/solve.sh` copies `solution/identify.py` to `/root/results/identify.py`. The program uses numpy only and builds its own network model from `network.json`. For fixed parameters θ, the true phasors of each frame enter linearly: the 13 complex free voltages, with the zero-injection bus 7 eliminated exactly, map through a 52 × 26 real matrix A(θ) to the 26 complex channels. Each frame's states are projected out by a batched QR (variable projection). The polar noise is whitened per frame and channel by its 2 × 2 inverse square-root covariance, evaluated at the measured phasor. The projected cost is minimised over θ by Levenberg–Marquardt with the exact Golub–Pereyra Jacobian; with the Kaufman approximation the superset fit stops short of the minimum on rec_b and the program fails that recording (calibration table). Series admittances are kept within 1.5 times nameplate magnitude and charging multipliers within [0, 1.5] during the search.

The pipeline:

1. **Stale data:** exact repeats of the previous frame, per channel.
2. **Superset model:** every non-known candidate in series-admittance coordinates (g, b_s, plus a charging multiplier for lines). It is converged from several starts, and the lowest cost is kept: nameplate values with database-open corridors at zero; calibrations fitted first on the database topology; admittances fitted first with calibrations frozen at zero; and, later, the previous robust estimate. Search fits use every fourth frame.
3. **Robust fit:** Cauchy reweighting per frame and PMU, scaled by the median residual statistic.
4. **Malfunctions:** runs of the rolling-median residual statistic per PMU and per channel are refined by fitting step, ramp and out-and-back shapes to per-frame estimates of a PMU phase offset and of channel ratio offsets. An interval is accepted if its χ² gain exceeds 200 and it overlaps no stronger one. The intervals are masked.
5. **Topology:** on the masked data, the superset is re-converged from all starts, then backward elimination on the penalised likelihood (cost + 30 per corridor). Then a plausibility check: a corridor whose fitted impedance has a non-positive resistance or reactance, or lies more than 50 % from nameplate (the instruction bounds every unknown to 10 %), is a compensating artefact of a false minimum; it is removed and the rest refitted from the current and from fresh nameplate starts, and the removal is kept if the penalised likelihood improves. Steps 3–5 repeat, at most three times, until the intervals stop changing.
6. **Final fit:** a maximum-likelihood fit in physical coordinates (r, x; x, tap; ratio, phase) on all frames, with the intervals masked.
7. **Report within the stated ranges:** each value is clipped to the range the instruction states for the true value (r and x within 10 % of nameplate, transformer x within 5 %, tap within 0.019, |ρ| ≤ 0.01, |φ| ≤ 0.0105). The truth lies inside, so clipping can only reduce an error. It matters for the weakly observed 9–14 line, whose Cramér–Rao σ on r is 9 % of r. On that scale the likelihood is visibly curved. On one Monte Carlo draw of rec_d the unclipped fit put r at 0.081 pu (−4.4σ, outside the stated range); along the straight line from the planted truth to that fit, the cost rises by 165 halfway and falls by 50 at the end. Clipped, the value is within tolerance.

Measured with the real verifier (`tests/test_identification.py`, run directly in the authoring container with nothing else running; see "Harbor runs" below): the four simultaneous runs finished after 1002 s (rec_d), 1146 s (rec_c), 1796 s (rec_b) and 2388 s (rec_a), and all 24 tests passed. Each run peaked at 0.87 GB of resident memory. rec_a is the slowest because of the plausibility check, about 300 s. v3 raised the wall limit from 3600 s to 5400 s, so the slowest run has a margin of 2.3×. That run predates the clipping of step 7. The final program was run through the verifier again in the same way: again 24/24, with all four runs finished within 2339 s and a peak of 0.80 GB each.

## Verification

**How the program is run.** The grading container holds Python 3.12 with numpy 1.26.4 and scipy 1.13.1, the four held-out recordings and their truth (in `/tests/data/rec_*`, stored as `heldout_pmu.csv`, `heldout_network.json`, `heldout_tolerances.json` and `heldout_truth.json`; `/tests` is mode 700), and the submitted program. The verifier opens `/root/results/identify.py` without following symlinks and without blocking (so a FIFO is rejected at once), requires a regular file and copies its content. For each recording, in shuffled order, it stages a fresh run directory: a read-only input directory with a random name holding only `pmu.csv` and `network.json` (copied from the `heldout_` files), and a working directory (mode 1777) holding only the copy of the program. It first checks, as the sandbox user, that the working directory and `/tmp` are writable, the input directory is not, and `/tests` is unreadable. It then starts all four runs at once as uid/gid 65534 in their own sessions, with the environment stated in the instruction and stdout/stderr going to unlinked files. Each run must exit with status 0 within 5400 s; a run past the limit has its process group killed. Afterwards every process of uid 65534 is killed: a helper that has become uid 65534 calls `kill(-1, SIGKILL)`, repeated until `/proc` shows none, so processes that keep forking cannot outlive the run. The output is then opened the same way as the program and must be a regular file, not a symlink, of at most 1 MB; at most 1 MB + 1 byte is ever read.

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

**Monte Carlo calibration.** The reference program was re-run on recordings with the same planted truth and malfunctions but new noise, load, dispatch, furnace and frequency realisations (one per held-out recording; generator seed override, author-side). Each was scored against its own tolerances, recomputed at its own states:

| recording | realisation | topology | malfunctions matched (max endpoint error) | worst \|error\| / σ | result |
|---|---|---|---|---|---|
| example | task data | correct | 1/1 (0) | 2.2 | pass |
| rec_a | task data | correct | 3/3 (1) | 2.3 | pass |
| rec_a | seed 11710 | correct | 3/3 (1) | 2.7 | pass |
| rec_b | task data | correct | 3/3 (4) | 1.8 | pass |
| rec_b | seed 11810 | correct | 3/3 (5) | 2.6 | pass |
| rec_c | task data | correct | 1/1 (0) | 3.5 | pass |
| rec_c | seed 2010 | correct | 1/1 (0) | 1.8 | pass |
| rec_d | task data | correct | 4/4 (0) | 2.0 | pass |
| rec_d | seed 12010 | correct | 4/4 (0) | 2.9 | pass |

All 9 runs pass, with the worst error at 3.5σ against the 4σ limit; the largest one (rec_c, I14 ratio error, 3.48σ) is reproduced to within 0.01σ by the independent program below, so it is a property of that noise realisation, not of the estimator. Before the clipping of pipeline step 7 was added, the seed-12010 run of rec_d failed on the 9–14 line (r at −4.43σ, x at +3.43σ). The likelihood check that explains it, with author-side `costcheck.py`:
- On all nine runs, the cost at the planted truth exceeds the cost at the fit by 44 to 85, a χ² of the order its 78–82 degrees of freedom imply. Generator, model and fit therefore agree, and no fit stopped short of its optimum.
- Small steps of every parameter change the cost exactly as the Fisher matrix predicts, within 0.03 %.
- Large steps along the weakest directions do not: there the likelihood is curved, and an occasional realisation lands far along the valley.

**Independent implementation.** A second program (author-side `independent.py`, about 720 lines, numpy only) was written without reusing code from the reference program. It optimises the same likelihood by different means at every stage:

| stage | reference program | independent program |
|---|---|---|
| residual | measured minus model, whitened per frame by the 2 × 2 inverse square-root polar covariance | relative residual 1 − h/z̃, real and imaginary parts divided by σ_m and σ_a |
| per-frame states | projected out by batched QR; Levenberg–Marquardt on the projected cost with the exact Golub–Pereyra Jacobian | re-solved exactly for every trial parameter; Levenberg–Marquardt with the Golub–Pereyra Gram matrix, including the second-order term C′A⁻¹C (verified against finite differences) |
| start | multi-start superset fit (nameplate, staged, admittance-first, warm) | one equation-error fit on the database topology |
| topology | superset, penalised backward elimination (cost + 30 per corridor), plausibility check | forward–backward stepwise search from the database topology: Rao score test (> 30) for adding a corridor, removals tried in implausibility and Wald order, each confirmed by a likelihood-ratio cut of 30 |
| outliers | Cauchy reweighting per frame and PMU | none beyond the flagged intervals |
| malfunctions | runs of a rolling-median statistic, then step, ramp and out-and-back shape fits to per-frame offsets | Page CUSUM on each PMU's mean-normalised χ², with change-point refinement of the ends; every round refits without the flagged data, re-detects all intervals and adds at most one new one |
| report | clipped to the stated ranges | clipped to the stated ranges |

On all five recordings it finds the true topology and every malfunction interval. Its worst endpoint error is 28 frames, on the start of rec_b's ramped timing error, which fades in from zero. Its worst errors are 2.5σ (example), 2.3σ (rec_a), 1.8σ (rec_b), 3.5σ (rec_c) and 2.5σ (rec_d) (3.0σ and 3.0σ on the example and rec_a before clipping). These per-recording figures come from a run before its normal-equation sums were accumulated in chunks of frames to cut memory; on the example, the chunked version reports the same values to within 0.0001σ. Submitted itself as `/root/results/identify.py` to the real verifier (as for the reference program above), it scores 24/24. All four simultaneous runs finished within 1843 s, each peaking at 0.57 GB. The two programs report the same topology and malfunction intervals (within 30 frames) on every recording. On rec_a, rec_b and rec_c their values differ by at most 0.25σ, where σ is the Cramér–Rao value of that quantity. On the example and on rec_d they differ by up to 2.0σ (example, branch 9 tap) and 1.5σ (rec_d, V10 ratio). There the costs of the two fits differ by only 3.3 and 3.1, against 85 and 50 between the planted truth and either fit: the difference lies along the flat, curved direction described above, not in an unconverged fit.

**Probes of the verifier.** Small programs that write crafted outputs, built author-side from the planted truth, were run through the verifier. Apart from the wall limit (3600 → 5400 s), the verifier is unchanged from v2, where the full table below was run in the verifier image. On the v3 data, the first row and the first two cases of the second row were re-run with the real verifier: the planted truth, with one branch value exactly at true + tolerance and every endpoint 30 frames off, scores 24/24. Moving that value 1e-60 further fails exactly `test_branch_parameters[rec_a]`, and moving one rec_b endpoint to 31 frames fails exactly `test_corrupted_segments[rec_b]`.

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

**Harbor runs.** v2 passed harbor 0.24.0 with `oracle` (reward 1, 24/24) and `nop` (reward 0, all 24 tests fail). The v3 results above come from the authoring container, where harbor and Docker are not available. There, `tests/test_identification.py` was run directly as root on the four held-out recordings, with the same sandbox: uid 65534, the stated environment, the read-only input directory, the self-check, and the process-group and uid kills. That container has Python 3.13 and numpy 2.5 rather than the grading image's Python 3.12 and numpy 1.26.4. Both programs use only basic numpy linear algebra, but the harbor `oracle` and `nop` runs should be repeated on the final package.

**Calibration: what weaker methods score.** Rows 1–5 are the reference program with one ingredient removed (author-side `make_ablations.py`, string patches of `identify.py`). Each was run author-side on the four held-out recordings, with the reference program's clipping applied to its output. Its four outputs were then graded by the real verifier through a replay program that writes the stored output for each recording. Cells give the worst |error|/σ over branch parameters and calibrations; "topology" means a wrong `in_service`, and "segments" means requirement 6 fails. Rows 6–7 were not run: they fail by construction, as stated in their cells.

| method | rec_a | rec_b | rec_c | rec_d | reward | verifier tests failed |
|---|---|---|---|---|---|---|
| 1. per-frame whitening, but at each channel's mean angle (one fixed covariance per channel) | 16.2, topology, segments | 262.2, topology, segments | 3.5 | 2.1 | 0 | topology, branch, calibration, segments on rec_a and rec_b |
| 2. isotropic noise (one σ per channel, rms of magnitude and angle σ) | 3.4 | 11.0, segments | 2.5 | 2.3 | 0 | branch, calibration, segments on rec_b |
| 3. Kaufman approximation of the variable-projection Jacobian | 2.3 | 5.5 | 3.5 | 2.0 | 0 | branch on rec_b |
| 4. superset fit from the nameplate start only | 2.3 | 4.7 | 3.5 | 2.0 | 0 | calibration on rec_b |
| 5. no plausibility check (pipeline step 5) | 36.8, topology | 1.8 | 3.5 | 2.0 | 0 | topology, branch, calibration on rec_a |
| 6. database topology trusted | topology | topology | not run | topology | 0 | by construction: the database is wrong for 2 candidates on rec_a, rec_b and rec_d |
| 7. exact repeats only, no other malfunction detection | segments | segments | segments | segments | 0 | by construction: every recording has a timing or ratio fault of 300–1200 frames, which repeats cannot reveal (stale intervals last 40–90 frames) |
| reference program | 2.3 | 1.8 | 3.5 | 2.0 | 1 | none |

Each ingredient of the reference program is needed on at least one held-out recording:
- **Noise model (rows 1–2).** A fixed or isotropic covariance mis-weights the angle noise as the phasors rotate. The PMU-level malfunction statistic then sees spurious excess χ², so intervals are missed or invented, and the parameters follow.
- **Exact Jacobian (row 3).** With the Kaufman approximation, every start of the superset fit on rec_b stops at a cost of at least 26.78 per search frame, against 25.87 with the exact Jacobian. The topology is still right, but the final fit stays in that basin and branch 21's resistance ends 5.5σ off.
- **Several starts (row 4).** On rec_b the reference program's best superset fit comes from the admittance-first start; the nameplate start alone ends at 27.72 per search frame and fails on the V12 ratio error (4.7σ).
- **Plausibility check (row 5).** See pipeline step 5 and the topology bullet under Difficulty.

## Version history

- **v3 (this version).** Rework after review: the v2 load model was not physically faithful (fast per-bus and common load fluctuations and 4–8 %/min system ramps with no stated mechanism). The operating conditions are now the mechanisms of the "Operating conditions" table, each at a realistic rate. The recording is 75 minutes at 2 frames per second, so the frame count and the contract are unchanged. The recordings, the tolerances and the example's truth were regenerated; the planted truth of the held-out recordings, the topology errors, the malfunctions and the noise levels are those of v2. The verifier code is unchanged. The reference program gained the plausibility check of pipeline step 5. On the new rec_a, every start of the superset fit keeps the open corridor 6–13 alive with a negative resistance, and without the check the program fails. The independent program was rewritten for the new data. The calibration table was recomputed on the v3 data. Two v2 rows were not recomputed and are dropped: two complete pipelines that had solved an earlier, steady-recording version of this task. Rows 1–2 of the new table cover the same noise-model weakness.
- **v2.** 15-minute recording at 10 frames per second, with the load model the review found unfaithful.
