Write a program, `/root/results/identify.py`, that identifies a 14-bus transmission network and its instrumentation from a 75-minute synchrophasor (PMU) recording. For a recording it must determine which candidate branches are in service, the electrical parameters of the in-service branches, the ratio and phase error of every non-reference instrument channel, and the time intervals in which a PMU malfunctioned. The program must also work on four other recordings of this format that are not provided here. `/root/data/` holds an example recording in the same format: `pmu.csv`, `network.json`, and `tolerances.json` (the accuracy tolerances for that recording; the four other recordings come without them).

## Network and measurement model

- All quantities are per unit on a 100 MVA base. Bus voltages and bus current injections are phasors in rectangular form $ V = |V| e^{j \theta} $, angles in radians.
- `network.json` lists 24 candidate branches (`candidates`), each with an `id`, `from_bus`, `to_bus`, a `model` and `nominal` (nameplate) parameters. Branch admittance follows MATPOWER conventions (`makeYbus`), with an off-nominal tap on the from side:
  - `line`: series impedance $ r + jx $ and total charging susceptance $ b $ split equally between the two ends; $ b $ is exactly the nominal value, $ r $ and $ x $ are unknown.
  - `series`: series impedance $ r + jx $, no charging; $ r $ and $ x $ unknown.
  - `transformer`: series reactance $ jx $ and a real off-nominal ratio `tap` on the from side; $ x $ and `tap` unknown.
  - `reactance`: series reactance $ jx $ only.
  - A candidate with `"known": true` is in service and its nominal parameters are exact.
  - The unknown $ r $ and $ x $ of lines and series branches lie within 10 % of nameplate, transformer $ x $ within 5 %, and `tap` within 0.019 of nameplate.
- `database_status` is the operator's topology database (`closed` or `open`). In a recording, the database is wrong for 0 to 2 candidates, in either direction; the true status of every candidate is to be determined from the data. An out-of-service branch carries no current and contributes nothing, including no charging.
- The bus listed in `shunt` has a shunt susceptance exactly equal to `nominal_b`. There are no other bus shunts.
- Bus `zero_injection_bus` (7) is the star point of a three-winding transformer. Its injection current is exactly zero and it has no PMU. Every other bus has a PMU measuring the bus voltage (channel `V<bus>`) and the bus current injection (channel `I<bus>`, positive into the network, so that $ I = Y_{\mathrm{bus}} V $).
- Each channel $ k $ reports $ \tilde z_k = (1 + \rho_k)\, e^{j \varphi_k}\, z_k\, (1 + n^{m}) e^{j n^{a}} $, where $ z_k $ is the true phasor, $ \rho_k $ the constant ratio error and $ \varphi_k $ the constant phase error (radians) of the channel's instrument transformer, uniquely fixed by the conventions $ 1 + \rho_k > 0 $ and $ -\pi < \varphi_k \le \pi $ (all true values satisfy $ |\rho_k| \le 0.01 $ and $ |\varphi_k| \le 0.0105 $), and $ n^{m} \sim \mathcal{N}(0, \sigma_m^2) $, $ n^{a} \sim \mathcal{N}(0, \sigma_a^2) $ are independent across channels and frames. $ \sigma_m $ (`magnitude_rel_std`) and $ \sigma_a $ (`angle_std_rad`) are given in `noise` for voltage and current channels; in the four other recordings $ \sigma_m $ is 1e-4 to 2e-4 for voltages and 2.5e-4 to 5e-4 for currents, and $ \sigma_a $ is 2 to 8 times $ \sigma_m $. The instruments of the PMUs in `reference_pmus` (buses 1 and 6) have $ \rho = \varphi = 0 $ exactly; all other channels have unknown $ \rho, \varphi $.
- Every frame is an AC power-flow state of the network. The operating point changes as in normal operation: the demand of each load bus varies slowly and randomly (minute-scale variation of the aggregate substation demand, with a small fast component) on top of a system load level that is steady, ramps up or down, or rises and falls, at most 0.3 % per minute; bus 3 also serves a steel works whose electric arc furnace makes its demand fluctuate from frame to frame and which can be switched off for part of a recording; the generator at bus 2 follows a 5-minute dispatch and automatic generation control, bus 1 is the slack, and buses 2, 3, 6 and 8 hold their voltage set points. Angles are measured against a nominal-frequency time reference, and the system frequency deviates from nominal, so the angles of all phasors of a frame can drift together by up to 0.65 rad per second and by up to about a hundred radians over a recording. Loads, generation and the true phasors of each frame are not given.
- Each recording contains 1 to 4 malfunction intervals. A malfunction affects one PMU (reference PMUs included) over consecutive frames, and its channels then do not follow the model above in one of these ways:
  - one or more of its channels repeat, on every frame of the interval, the exact values of the frame before it (stale data; 40 to 90 frames);
  - all of its channels carry an extra common phase shift (a timing error) that is constant, or grows from zero over the interval and jumps back after it, or grows and returns to zero, with largest magnitude between 0.002 and 0.03 rad (300 to 1200 frames);
  - one of its channels has its magnitude scaled by an extra constant factor $ 1 + \delta $ with $ |\delta| $ between 0.0008 and 0.002 (600 frames).

  No two intervals overlap in time, even on different PMUs, and the intervals are not flagged.

`pmu.csv` has one row per frame (2 frames per second, 9000 frames) with columns `frame` (0-based index) and, for every PMU bus $ b $, `V<b>_mag`, `V<b>_ang`, `I<b>_mag`, `I<b>_ang`. Reported angles are principal values in $ (-\pi, \pi] $. `network.json` also gives `frames`, `frame_rate`, `buses`, `zero_injection_bus`, `reference_pmus`, `noise` and `shunt`.

## Program contract

The program is run as `python3 identify.py <input_dir> <output_file>`:

- `<input_dir>` is an absolute path to a read-only directory holding only `pmu.csv` and `network.json`, in the formats of the example.
- The working directory is a fresh directory that holds only a copy of the program and that the program may write to, as well as `/tmp`, which is shared by the four simultaneous runs. `<output_file>` is the absolute path of `identification.json` in that directory, which the program must create.
- Each run has a wall-clock limit of 3600 seconds and must exit with status 0. Processes that a run leaves behind are killed before its output is read. The output must be a regular file (not a symlink) of at most 1 MB.
- The program runs as an unprivileged user (uid and gid 65534) with no network, an empty standard input, and exactly the environment variables `HOME=/tmp`, `PATH=/usr/local/bin:/usr/bin:/bin`, `LANG=C.UTF-8` and `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS` all set to 1.
- The program must run with Python 3.12, numpy 1.26.4, scipy 1.13.1 and the standard library; no packages can be installed. The four recordings are processed by four runs started at the same time, which share 4 CPUs and 8 GB of memory.

## Output

`identification.json` is one JSON object, encoded in UTF-8, with:

- `branches`: a list with exactly one entry per candidate `id`, 24 entries in any order. Each entry is an object with `id` (a JSON integer, written without a decimal point or exponent), `in_service` (`true` or `false`) and, when `in_service` is `true` and the candidate is not `known`, its unknown parameters as JSON numbers: `r` and `x` for `line` and `series`, `x` and `tap` for `transformer`, in per unit (`tap` dimensionless). Parameter keys on other entries may be present and have no accuracy requirement.
- `channels`: an object with exactly one key per channel with unknown errors, 22 keys in all: `"V<b>"` and `"I<b>"` for every PMU bus $ b $ except the reference buses 1 and 6 (`"V2"`, `"I2"`, `"V3"`, …). Each maps to an object with `ratio_error` ($ \rho_k $, dimensionless, with $ 1 + \rho_k > 0 $) and `phase_error` ($ \varphi_k $ in radians, the principal value with $ -\pi < \varphi_k \le \pi $), both JSON numbers. The instrument error must be reported in this representation, with the principal phase value and a positive ratio factor $ 1 + \rho_k $, not in an equivalent one such as $ \varphi_k + 2\pi $, or $ -(2 + \rho_k) $ together with $ \varphi_k + \pi $.
- `corrupted_segments`: a list, in any order, with exactly one entry per malfunction interval. Each entry is an object with `bus` (the bus of the malfunctioning PMU), `first_frame` and `last_frame` (the first and last affected frame, inclusive), all three JSON integers written without a decimal point or exponent, with $ 0 \le $ `first_frame` $ \le $ `last_frame` $ \le $ `frames` $ - 1 $.

The top-level object must contain these three keys, and the `channels` object must contain exactly the 22 channel keys. Other objects (the top level, branch entries, channel entries, segment entries) may carry extra keys; these have no content requirement, but they must still satisfy the number and key rules in the next paragraph.

Every number literal anywhere in the file, including in any extra field, must be a JSON number whose exact decimal value, before any rounding, has absolute value at most the largest finite IEEE-754 double, $ (2 - 2^{-52}) \cdot 2^{1023} $. `NaN`, `Infinity` and larger literals (such as `1e309`, or an integer one above that bound) are not allowed. No object may repeat a key.

## Requirements

The task is complete only if every requirement below holds on all four recordings that are not provided here.

1. The run exits with status 0 within the wall-clock limit and leaves `identification.json` as stated in the program contract.
2. The file and schema are exactly as above: all 24 candidates, all 22 channels and the segments present with the stated types and ranges, every number literal within the bound above and no object key repeated.
3. `in_service` matches the true status of every candidate.
4. For every truly in-service candidate that is not `known`, every unknown parameter differs from its true value by at most its tolerance.
5. Every `ratio_error` and `phase_error` differs from its true value by at most its tolerance.
6. The number of segments equals the number of malfunction intervals, and the segments and intervals can be paired one-to-one so that each pair has the same bus, and the `first_frame` and the `last_frame` of each pair each differ by at most 30 frames.

In requirements 4 and 5 the difference is the absolute difference, with no wrapping of angles modulo $ 2\pi $, between the exact decimal value of the reported number and the true value, and a difference equal to the tolerance satisfies the requirement. Each tolerance is four standard deviations of the Cramér–Rao bound for that quantity, computed at the true values of that recording over its frames free of malfunctions, as in the example's `tolerances.json`. On the four other recordings the tolerances range from 2.0e-4 to 6.7e-2 per unit for branch parameters, 1.4e-4 to 6.6e-2 for ratio errors and 2.9e-4 to 8.6e-2 rad for phase errors.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
