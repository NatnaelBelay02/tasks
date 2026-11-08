# lc-ladder-synthesis

The agent writes a filter synthesis program (the insertion-loss method [1–5]). It is graded by running it on 67 held-out specifications of application filters, 50–600 Ω, from 20 Hz to 2.5 GHz.

- **Low-pass, degree 5 to 16:** transmitter harmonic filters with a stopband from the 2nd harmonic (one of them all-pole, with no finite zeros); video and DAC reconstruction filters; anti-alias filters of audio converters and ADCs (one with a single finite zero); an ADSL splitter; two cable-TV diplexer arms.
- **Band-pass:** receiver preselectors and IF filters with stopbands on neighbouring services (including a medium-wave preselector with a single zero at dc, and six receivers, from 77.5 kHz to 2.2 GHz, with one zero at dc and two zeros below the band), a CB-band transmitter filter with no finite transmission zeros, and three receivers with a single finite zero: a 49 m broadcast preselector of degree 6, an ADS-B receiver with only a lower stopband and a GPS L1 preselector with only an upper one.
- **Narrow IF filters of degree 14 and 16:** ten IF filters at 455 kHz, 10.7, 21.4, 45 and 70 MHz with 0.11 % to 0.71 % bandwidth (SSB, CW, narrow-band FM and land-mobile channel filters). These LC ladders are the prototypes from which crystal and ceramic ladder filters are designed.
- **Wide band-pass filters, $ f_2/f_1 $ from 200 to 3333:** ten front-end and instrument filters of degree 10 to 16: a CISPR band-B and two 9 kHz–30 MHz EMC-receiver preselectors, LF–HF and VLF–MF receiver front ends, an eddy-current instrument filter, sonar and instrumentation input filters, and an audio band-limiting filter in front of a 48 kHz converter.
- **Full Cauer filters, with as many finite zeros as the degree allows (6 or 7):** three converter low-pass filters of degree 13, 15 and 16 (LP-DAC30M, LP-SDR2M4, LP-ADC250), and ten band-pass filters of degree 16: narrow-band FM and roofing IF filters at 10.7 and 21.4 MHz (0.2 %), TV band III, DAB and UHF TV preselectors, a 2.4 GHz WLAN front-end filter, a telephone voice-band filter, and HF and LF/MF receiver front ends ($ f_2/f_1 $ up to 30).

For each specification the program must produce a doubly terminated, transformerless LC ladder that meets four conditions:

- it realizes the **generalized elliptic** response to $ 10^{-5} $ dB: equiripple passband, and equiripple stopbands whose finite transmission zeros the program has to place;
- it is **canonical** ($ n + m $ elements);
- its **load** is within a factor 2 of the source resistance;
- it has the **fewest inductors** of all ladders that meet the other three (version 4.2.0).

The specification gives, for each stopband, its edge (where the application needs rejection: a harmonic, an image or alias band, the next service) and the number of finite zeros in it. The zero frequencies are not given. They are fixed by the condition that every attenuation minimum in a stopband, including the one at its edge, is equal: the Cauer (elliptic) property, here for band-pass filters with an odd or even number of zeros at dc and an unequal number of zeros on the two sides.

## Difficulty

1. **Approximation: placing the transmission zeros.**
   - For given zeros the equiripple passband is a weighted minimax problem in $ u = \omega^2 $, with weight $ u^{(r-a)/2}/\prod|u - u_k| $, solved by a Remez exchange [6, 7]. The zeros themselves are unknowns: each stopband with $ m_b $ zeros adds $ m_b $ conditions (its $ m_b + 1 $ attenuation minima are equal), and these are coupled to the passband through the reflection zeros. The program has to solve the whole system.
   - Closed forms do not cover the held-out specifications. Elliptic functions give the low-pass case with equal-ripple stopband from one edge; the low-pass to band-pass transformation gives only zeros geometrically symmetric about the centre frequency, with as many zeros at dc as at infinity. The held-out band-pass filters have 1 to 4 zeros at dc, unequal zero counts on the two sides, and stopband edges that are not geometric images of each other.
   - The zeros must be accurate. Re-synthesising the ladder with the reference zeros moved by $ 10^{-7} $ (relative) fails for all 65 specifications with finite zeros ($ 4.4 \times 10^{-5} $ to $ 2.2 \times 10^{-2} $ dB). A double-precision solve (a discrete Remez on a 20 000-point grid with scipy's fsolve on the zeros, the start of `make_truth.py`) misses the zeros of BP-915, BP-GSM900 and BP-DCF77 by $ 3 \times 10^{-4} $ to $ 5 \times 10^{-4} $ (relative) and fails them by 0.2 to 0.4 dB.
   - Moving the reflection zeros by $ 10^{-6} $ (relative) changes the insertion loss by $ 1.6 \times 10^{-5} $ to $ 9.2 \times 10^{-2} $ dB, above the tolerance. Zeros placed by a reasonable guess (spaced out from the edge, not equiripple) are off by 1.2 to 156 dB.
2. **Zero shifting at dc.**
   - The textbook elliptic-filter procedure creates each finite zero by partial removal of a pole at infinity [3, 4].
   - It finds **no ladder at all** for 23 of the 37 held-out specifications of degree up to 12: every band-pass filter there with a finite zero below the passband (all with finite zeros except BP-49M and BP-GPSL1, whose single zero lies above it).
   - These need partial removal of the pole at $ s = 0 $ (a series capacitor or shunt inductor) before the resonator.
3. **Choice of realization: an admissible load and the fewest inductors.**
   - A specification has many canonical ladders. They differ in the order of the extractions and in the kind of each element: a zero at dc is a series capacitor or a shunt inductor, a zero at infinity a series inductor or a shunt capacitor, and a finite zero is shifted by a partial removal at $ 0 $ or at $ \infty $, on the impedance or on the admittance. They also differ in the load they need.
   - The first ladder a depth-first extraction finds fails for 27 of the 29 band-pass specifications other than the wide ones: it needs a load outside $ [R_s/2, 2R_s] $, in some cases by many orders of magnitude. Only 4 to 528 of the 20 to 9072 ladders of such a specification of degree up to 12 have an admissible load.
   - Among the ladders with an admissible load the number of inductors varies: for the 54 specifications of version 4.1.0 there are 2 to 432 432 such ladders (BP-10K10M), and for 52 of them more than one inductor count, up to 4 apart. Requirement 4 asks for the fewest. This is the usual practical choice, because the inductors are the lossy, bulky and costly parts; handbooks tabulate minimum-inductor forms for this reason [5]. The first admissible ladder of a depth-first search has the fewest inductors for only 11 of the 54 (8 of the 15 low-pass, 3 of the 39 band-pass). Between 6.5 % and 100 % of a specification's admissible ladders are minimal.
   - So the program has to search all the realizations, not stop at the first. Enumerating them one by one is too slow: listing the 432 432 admissible ladders of BP-10K10M takes 261 s in the reference code, and the full Cauer filters (six or seven finite zeros, degree 16) have far larger trees. The search has to be pruned. A partial removal exists only where the pole exists, and only by less than its full residue. A branch can be abandoned once its inductors plus one per finite zero still to be extracted (each needs a resonator, which contains one inductor) reach the best ladder found so far (branch and bound).
4. **Numerical conditioning: the narrow IF filters.**
   - The synthesis works on polynomial coefficients in $ s $. For a band-pass filter of relative bandwidth $ b = (f_2 - f_1)/\sqrt{f_1 f_2} $ the $ 2n $ roots of $ E E_* $ and of $ E \pm F $ crowd into two clusters of width about $ b $ around $ \pm j\omega_0 $. The coefficients then fix those roots only to about $ n \log_{10}(1/b) $ digits fewer than they carry. The spectral factorization $ E E_* = F F_* + P P_* $ and every step of the continued-fraction extraction, which subtracts nearly equal polynomials, work at that loss. So does each zero test the extraction makes ("is the residual zero at $ s = 0 $, or at $ j\omega_k $?").
   - $ n \log_{10}(1/b) $ is at most 13 for the examples, at most 20 for the band-pass specifications of degree 6 to 12, and 30 to 46 for the ten IF filters of degree 14 and 16 (0.11 % to 0.71 % bandwidth).
   - The tolerance needs the ladder to reproduce the response over more than 150 dB of attenuation. So the program needs a working precision chosen for the filter in hand, well above this loss, or a formulation that avoids the loss (transformed variables [8]). A fixed precision that works on the examples does not carry over. The reference sets its precision from $ n $ and $ b $ (100 to 180 digits) and scales every tolerance with it.
   - The precision is not an author-chosen hurdle. The tolerance is fixed by the response ($ 10^{-5} $ dB, see the threshold rationale), and the digits needed follow from the conditioning of the extraction, which is poor for narrow-band filters [8]. The instruction states the scope (degree up to 16, bandwidth down to 0.1 %). Narrow LC ladders of this degree are what a crystal or ceramic ladder filter is designed from. Any method that controls the error passes, whether by more precision or by a scaled formulation (see "Second implementation and model runs").
5. **Numerical conditioning: the wide band-pass filters.**
   - In $ u = (f/f_2)^2 $ the passband is $ [u_1, 1] $ with $ u_1 = (f_1/f_2)^2 $, here $ 9 \times 10^{-8} $ to $ 2.5 \times 10^{-5} $. With 2 to 4 zeros at dc the weight $ u^{(r-a)/2} $ is very large near $ u_1 $, and the lowest reflection zeros crowd geometrically towards the edge: the lowest lies at 1.06 to 1.10 $ u_1 $, the next at 1.8 to 2.7 $ u_1 $, and for eight of the ten a third at 12 to 14 $ u_1 $. In the variable $ x \in [-1, 1] $ in which a Remez exchange on $ [u_1, 1] $ is usually written, the lowest zero lies $ 10^{-8} $ to $ 3 \times 10^{-6} $ from $ x = -1 $.
   - Near $ x = -1 $ the polynomial is a product of these small distances, far below the rounding of its coefficients in double precision. Given the exact transmission zeros, a double-precision discrete Remez on a 20 000-point grid (the start of `make_truth.py` for the narrower filters) gets the reflection zeros of the ten wide filters wrong by $ 1.7 \times 10^{-4} $ to $ 5.8 \times 10^{3} $ (relative), against $ 10^{-10} $ to $ 4 \times 10^{-9} $ for the other band-pass specifications. A Remez exchange for these filters has to be carried out in extended precision, or in a variable that resolves the edge (for example $ \log u $).
   - The instruction states the scope (passbands from audio frequencies, $ f_2/f_1 $ up to 3500). Filters like these are the preselectors of EMC test receivers (CISPR 16-1-1, bands A and B) and the input filters of broadband receivers and instruments.
6. **Held-out grading.** The program runs on specifications it has never seen, so it must be general. The 67 specifications cover the whole stated scope:
   - low-pass of every degree from 5 to 16, odd and even; band-pass of every even degree from 6 to 16;
   - 1, 2, 3 and 4 zeros at dc for band-pass;
   - 0 to 7 finite zeros in total (an empty list of stopbands for a low-pass and for a band-pass filter), in a lower stopband only, an upper stopband only, or both, with up to 7 in one stopband (full Cauer filters: LP-SDR2M4, LP-ADC250);
   - bandwidths from 0.11 % to $ f_2/f_1 = 3333 $, passbands from 20 Hz to 2.5 GHz.

## Reference solution

`solution/synthesize.py` works in mpmath, in units normalized to $ R_s $ and $ 2\pi f_2 $. Its precision is $ 40 + 3 n \log_{10}(1/b) $ digits rounded up to a multiple of 10, and at least 100: 100 digits up to degree 12, 140 to 180 for the narrow IF filters. The tolerances below are quoted at 100 digits; each is a fixed power of the precision (`set_precision`).

1. **Approximation.** Remez exchange in $ x \in [-1, 1] $, with $ u = \tfrac{1+u_1}{2} + \tfrac{1-u_1}{2}x $:
   - The reference points are the $ d+1 $ extrema of $ W f $. They are found between consecutive roots of $ f $ by bisection on $ \tfrac{d}{dx}\log|Wf| $, with the passband edges taken when they are extremal.
   - The iteration runs until the extrema agree to $ 10^{-60} $.
   - Each exchange starts from the extremal points of the previous one, and the roots of $ f $ are refined from the previous roots by Newton's method (`polyroots` only when that fails). In the 30-digit zero search below, the widest band-pass filters reach a rounding floor above that tolerance; there the exchange stops once the spread of the extrema no longer shrinks.
   - **Zero placement.** The zeros of a stopband are parametrised so that they stay ordered and outside the edge. For trial zeros the Remez exchange is run, and each interior attenuation minimum is found by bisection on $ \tfrac{d}{du}\log|Wf| $ in $ \log u $. The residuals are $ \log|Wf| $ at those minima minus $ \log|Wf| $ at the edge. `mpmath.findroot` solves them at the working precision. Its starting values come from scipy's fsolve on a double-precision model of the same conditions, in which the passband reflection zeros are taken from the closed form of the generalized Chebyshev function (its phase $ \sum_p w_p \arccos M_p(u) $ is $ (k + \tfrac12)\pi $ at the $ k $-th reflection zero, with $ M_p $ the Möbius map of $ [u_1, 1] $ onto $ [-1, 1] $ with its pole at the transmission zero $ p $). Only if that model gives no solution does fsolve run on the Remez step at 30 digits (version 4.1.0 always did that, which stalls with four or more zeros in one stopband).
   - $ c $ is then set so that the maximum of $ |K| $ is $ \varepsilon $.
2. **Feldtkeller equation.** The roots of $ Q(w) = P(s)P(-s) + F(s)F(-s) $, $ w = s^2 $, give $ E $ from the left half-plane roots $ s = -\sqrt{w} $, with leading coefficient $ c $.
3. **Extraction.**
   - The input impedance is $ Z = (E \pm F)/(E \mp F) $. A depth-first search reduces it to a resistor; at each step it tries, on $ Z $ or $ Y $:
     - full removal of a pole at $ \infty $ or at $ 0 $;
     - for each remaining finite zero $ \omega_k $, partial removal of the pole at $ \infty $ or at $ 0 $ that makes the remainder vanish at $ j\omega_k $, followed by full removal of the resulting pole pair as a resonator. The resonator is a shunt series-LC or a series parallel-LC.
   - Each step requires a positive element. Polynomial divisions must leave a remainder below $ 10^{-40} $ relative.
4. **Choice: branch and bound.** The search keeps only ladders with $ R_s/2 \le R_L \le 2R_s $ and runs over the whole tree, for both signs of $ F $. A branch is cut as soon as its inductors plus the number of finite zeros still to be extracted reach the inductor count of the best ladder found so far, so each ladder it keeps has fewer inductors than the one before, and the last one kept has the fewest. The moves that add a capacitor (series C for a pole at $ 0 $, shunt C for a pole at $ \infty $) are tried first, so that a good bound is found early. Among ladders with equally few inductors the program would take the load closest to $ R_s $.

Run time on one core (four specifications at a time on a 4-core machine): up to 8 s per specification of degree up to 12 other than the wide ones, 6 s to 8 s for the low-pass filters of degree 13 to 16 with up to three zeros, 12 s to 89 s for the narrow IF filters (140 to 180 digits), 8 s to 93 s for the wide band-pass filters, 12 s to 37 s for the full Cauer low-pass filters and 68 s to 104 s for the full Cauer band-pass filters. The branch and bound accounts for most of the time of the larger ones; in version 4.1.0, which stopped after 40 admissible ladders, the maximum was 51 s.

## Verification

`tests/test_synthesis.py` runs the agent's program on each held-out specification:

- the program, the specification and an empty working directory are in separate temporary folders, and the program is called with absolute paths as uid 65534;
- there is no network and a 300 s limit; the program is pinned to one CPU core (`sched_setaffinity`, inherited by every thread and child process), and every process has a 3 GB address-space limit (`RLIMIT_AS`);
- the program's output goes to files, not pipes. Its process group is killed when the time limit expires and again after it exits, so a left-over child process cannot stall the verifier;
- `test.sh` makes `/tests` unreadable to other users first, so the program cannot read the held-out specifications or the reference.

The network file is then read as data:

- regular file, no symlinks, at most 1 MiB;
- strict JSON (no NaN, Infinity or repeated keys);
- exact keys, positive finite values, known arm types.

Each held-out specification is one test, and every condition must hold:

- **Elements:** at most $ n + m $ elements, computed from the specification.
- **Load:** $ R_s/2 \le R_L \le 2R_s $.
- **Inductors:** at most `min_inductors`, the fewest of any canonical ladder of the specification with an admissible load (counted from the arm types: `L`, `series_LC` and `parallel_LC` contain one inductor each).
- **Insertion loss:** within $ 10^{-5} $ dB of the reference at 1121 to 1270 frequencies:
  - 1001 in the passband: 401 uniform, 400 spaced geometrically in $ f - f_1 $ and 200 in $ f_2 - f $, both from $ 10^{-8}(f_2 - f_1) $;
  - 120 above the passband, log-spaced from $ f_2 $ up to $ 30 f_2 $, and for band-pass filters 120 below it, log-spaced from $ f_1 $ down to $ f_1/30 $;
  - the stopband edges, always. An edge is checked even when it lies within $ w $ of a finite zero (below), as 55 of the 107 edges do; for BP-IF21M4N the edge is $ 2.5 \times 10^{-5} $ (relative) from a zero;
  - 4 next to each reference finite zero, at $ \pm 2w $ and $ \pm 6w $ (relative);
  - other grid points within $ w $ of a reference finite zero are omitted. $ w = 0.5\,\% $, scaled down in proportion for band-pass filters narrower than 1.6 % ($ w = 0.5\,\% \cdot b / 1.6\,\% $, 0.034 % to 0.22 % for the narrow IF filters). So $ \pm 2w, \pm 6w $ are ±1 %, ±3 % for all other specifications.

  The ladder is analysed by its chain matrix in 40-digit arithmetic. The transducer gain is $ 4R_sR_L/|AR_L + B + CR_sR_L + DR_s|^2 $.

Any ladder that meets the conditions passes, not only the reference one. The reference values the verifier uses are the certified characteristic function (below), with its reflection zeros, transmission zeros and $ c $, and the minimum number of inductors.

**Why the thresholds are what they are.**

- **$ 10^{-5} $ dB.** A ladder with exact element values reproduces the reference to about $ 10^{-11} $ dB in the reference runs. Scaling every element of the stored witness by $ 1 + 10^{-11} $ still passes for every specification (the stopband edges included), by $ 1 + 10^{-10} $ for 64 of the 67 and by $ 1 + 10^{-9} $ for 51. The ones that fail at $ 10^{-10} $ are three narrow IF filters (BP-IF21M4N, BP-NBFM10M7, BP-ROOF21M4), whose response is about a thousand times steeper and whose stopband edges lie close to a zero. $ 10^{-11} $ is still four orders of magnitude above the rounding of element values to double precision; a ladder written with double-precision values is off by at most $ 4.8 \times 10^{-11} $ dB (the stored witness of BP-IF21M4N). One inductor off by $ 10^{-4} $, or reflection zeros off by $ 10^{-6} $, fails (see the ablations).
- **1001 passband points.** The passband loss has at most $ d + 1 = 9 $ ripple maxima and $ d $ reflection zeros. The ripples are narrowest at the band edges, and for the wide band-pass filters the lowest ones span only a few % of $ f_1 $, $ 10^{-5} $ of the band, which a uniform grid misses. With the points clustered towards both edges, every interval between neighbouring reflection zeros (and between the edges and the outermost ones) holds at least 27 points (BP-LFMF), for all 67 specifications. Both responses are rational functions of $ \omega^2 $ of degree at most 16. An error therefore cannot fit between sample points that are this dense; it has to be smooth over the band.
- **120 log-spaced points on each side of the passband, out to a factor 30 from its edges.** That is 81 points per decade. They cover the transition band and the stopband up to that factor; the stopband edges are added as fixed points. Between 64 (LP-POTS, whose stopband starts at $ 6.5 f_2 $) and 132 grid points fall inside each stopband, against at most $ m_b + 1 = 8 $ attenuation minima to resolve. The stopband loss varies smoothly: it is high near each zero and dips once between neighbouring zeros. Beyond the outermost zero it rises monotonically towards its asymptote $ \propto \omega^{2b} $, whose coefficient also fixes the loss at all the other points, so no points are needed beyond the factor 30.
- **Points next to a finite zero.**
  - Near a zero $ f_z $, the relative error of $ G $ caused by a relative error $ \delta $ in the resonance frequency is about $ 2\delta/|f/f_z - 1| $.
  - Grid points too close to a zero would make the test depend on how the element values are rounded. With $ 10^{-9} $ scaling, a grid point $ 4 \times 10^{-4} $ from a zero would change by about $ 2 \times 10^{-5} $ dB, and closer points by more.
  - Grid points within $ w $ are therefore dropped and replaced by fixed points at $ \pm 2w $ and $ \pm 6w $. With $ w = 0.5\,\% $, a resonance misplaced by more than about $ 3 \times 10^{-8} $ (relative) is still detected, but double-precision values are not penalised.
  - The stopband edges are the exception: the instruction names them, so they are always checked. That is stricter next to a zero, as intended. The witness of BP-IF21M4N with one resonator capacitor scaled by $ 1 + 10^{-10} $ is off by $ 2.8 \times 10^{-6} $ dB on the grid but by $ 1.8 \times 10^{-5} $ dB at its edge (21.43 MHz), and fails. Double-precision rounding is still far below the tolerance there (see $ 10^{-5} $ dB above).
  - For a narrow band-pass filter a fixed 0.5 % would cover the whole passband (the zeros of BP-IF21M4N are 0.07 % from its edge), and ±1 % would land next to other zeros. So $ w $ shrinks in proportion to the bandwidth below 1.6 %, and points placed next to one zero that fall within $ w $ of another are dropped. Every specification keeps all 401 passband points.
- **Minimum number of inductors.** `min_inductors` is found by the reference search run as an exhaustive branch and bound (both signs of $ F $, every extraction order, partial and full removals at $ 0 $ and $ \infty $ on $ Z $ and $ Y $). The ladder that attains it is stored as the witness, and it passes the verifier. A second, separately written search (`alt_synthesize.py`, with its own bookkeeping of the zeros at dc and at infinity) finds the same minimum for every specification (`check_min_inductors.py`). The test is "at most", so a program that found a ladder with even fewer inductors would pass.
- **300 s per specification.** The reference program needs at most 104 s (BP-LFMF). The branch and bound shown under difficulty 3 brings a separately written depth-first search under 200 s on every specification, so the limit leaves room for a correct approach other than the reference's, while it still rejects a search that lists the realizations one by one.
- **1 MiB network file** (stated in the instruction). A canonical ladder here has at most $ n + m = 23 $ elements, so at most 23 arms. The reference output is about 2 KB.

**Where the specifications come from** (`tests/truth_generator/gen_specs.py`). Each specification is an application filter:

- the passband is the service band;
- the return loss is 14–26 dB (14 dB for the telephone splitter, 16 dB for the powerline coupler);
- the terminations are the system impedance (50, 75 or 600 Ω);
- each stopband starts where the application needs rejection:
  - the 2nd harmonic of the carrier (6.78, 13.56, 27.12, 144 MHz);
  - the first image of a 13.5 MHz video DAC and of a 2.5 GSPS DAC;
  - the first frequency aliased into the passband of a 96 kHz audio converter and of 50, 76.8 and 122.88 MSPS ADCs;
  - the first ADSL tone above a telephone band;
  - channel 2 above a cable-TV upstream band, and the downstream band above a DOCSIS 3.0 mid-split upstream band;
  - the GNSS L-bands above an 868 MHz transmitter;
  - neighbouring services of a preselector or IF filter: TV channel 4, the airband, the 2 m amateur band, TV band III, LTE bands 38 and 40, GSM-850, GSM-900 (uplink and downlink), the 700 MHz band, weather satellites at 137 MHz, the CB band, UMTS band 34 and the satellite S-band, the 41 m broadcast band, SSR interrogations at 1030 MHz (below the ADS-B band), the Globalstar and Iridium uplinks above GPS L1, and the adjacent channels of IF and time-signal receivers;
- the number of zeros per stopband is the one a designer would allot for that rejection;
- the wide band-pass filters cover the band of their receiver or instrument (CISPR 16-1-1 bands A and B, LF–HF, VLF–MF, the audio band), with stopbands below and above it (mains harmonics, infrasonic rumble, the alias band of a 48 kHz converter, or the bands next to a receiver's tuning range);
- the narrow IF filters (455 kHz SSB and CW, 10.7, 21.4, 45 and 70 MHz) have the bandwidth of their channel (500 Hz to 500 kHz, 0.11 % to 0.71 %), stopbands on the adjacent channels or the opposite sideband, and degree 14 or 16. Ladders like these are the LC prototypes from which crystal and ceramic ladder filters are designed;
- the full Cauer filters (version 4.2.0) put every transmission zero the degree allows into the stopbands:
  - the alias band of a 2.4 MSPS SDR receiver and of a 250 MSPS ADC (seven zeros) and the first image of a 30.72 MSPS DAC (six);
  - the adjacent channels of a 10.7 MHz narrow-band FM IF and of a 21.4 MHz roofing filter;
  - the FM band and the upper part of band III around a band III preselector;
  - LTE bands 40 and 38 around the 2.4 GHz WLAN band;
  - the bands below and above the DAB band III, UHF TV (with the aeronautical DME band above 960 MHz), HF (medium-wave broadcast below, VHF above) and LF/MF receiver bands;
  - mains hum and the band above 5 kHz around the telephone voice band.

  Five of these (BP-NBFM10M7, BP-ROOF21M4, BP-BANDIII, BP-HF3M30M, BP-LFMF) are frequency-scaled copies of randomly generated full Cauer specifications whose search was found to be hard, placed on the band of an application (frequency scaling leaves the normalized synthesis problem unchanged).

The degree is the one such a design needs, and $ n - a - 2m \ge 1 $. Band-pass degrees are even and have at least one zero at dc. The minimum stopband attenuations of the reference functions are 23 to 199 dB. A transmitter harmonic filter of high degree is specified as all-pole (LP-6M), and three band-pass receivers have a single finite zero.

6 example specifications, with reference insertion loss at 60 frequencies each, are given to the agent. The 67 graded specifications are not.

**Where the reference values come from** (`tests/truth_generator/make_truth.py`). This is separate code from the oracle:

1. **start:** for trial zeros, a discrete weighted Remez exchange in double precision on a 20 000-point grid in $ u $ (for the wide band-pass filters, $ u_1 < 10^{-2} $, at 40 digits on 4001 Chebyshev points plus 4001 spaced geometrically towards $ u_1 $, see difficulty 5); scipy's fsolve moves the zeros until the attenuation minima of each stopband are equal. The minima are searched on a grid clustered geometrically towards the zeros, which for the narrow IF filters lie within a fraction of a percent of the band edge. The stationary points are then moved onto the zeros of $ \tfrac{d}{du}\log|Wf| $ by bisection in mpmath;
2. **one joint Newton solve** in mpmath at 60 digits (plain Newton steps with the Jacobian from `mp.jacobian`). The equations are:
   - the passband equioscillation ($ Wf = \pm h $ at the $ d+1 $ extrema, zero log-derivative at the interior ones);
   - in each stopband, $ \log|Wf| = \lambda_b $ at the edge and at the $ m_b $ interior minima, with zero log-derivative there.

   The unknowns are the reflection zeros, the transmission zeros, the extremal points, $ h $ and $ \lambda_b $. This formulation differs from the oracle's nested solve;
3. **certificates:**
   - **passband:** $ |K| $ is evaluated on 6000 passband points (plus 6000 spaced geometrically towards $ u_1 $ for the wide band-pass filters) and every local maximum is refined by golden-section search. The maximum must equal $ \varepsilon $ to $ 10^{-30} $ and be attained, with alternating sign, at exactly $ d+1 $ points. By the alternation theorem [7] the passband polynomial is then the unique minimax one for these zeros.
   - **stopbands:** $ |K| $ is evaluated on 4000 log-spaced points of each stopband, from the edge to a factor 100 beyond the outermost zero, plus 4000 between the edge and the outermost zero and 300 on either side of every zero (geometrically closer towards it), and every local minimum is refined. The minimum over the band must be attained at exactly $ m_b + 1 $ points, equal to $ 10^{-30} $, and no grid point may lie below it. Every zero must lie inside its stopband.

The reflection and transmission zeros agree with the oracle's to $ 4 \times 10^{-61} $ (relative).

The truth file also stores, per specification, `min_inductors` (requirement 4) and a ladder with an admissible load that has that many inductors. Both come from the reference search run as an exhaustive branch and bound. The witness shows that the minimum is attained, and the ablations check it with the verifier's own analysis. The verifier uses `min_inductors` but not the witness.

The ladders with an admissible load, listed one by one for the 54 specifications of version 4.1.0 (`tests/truth_generator/count_inductors.py`, which uses the reference search without the bound). The branch and bound gives the same minimum for all 54:

| | ladders with admissible load | fewest inductors | spread of inductor counts | first one found is minimal |
|---|---|---|---|---|
| low-pass, degree 5–16 (15) | 2–1980 | 2–8 | 0–4 | 8 of 15 |
| band-pass from 1.3 % wide to $ f_2/f_1 = 3.2 $, degree 6–12 (18) | 4–528 | 3–6 | 0–3 | 1 of 18 |
| band-pass below 1 % (narrow IF and BP-ADSB), degree 8–16 (11) | 16–7236 | 4–8 | 1–3 | 2 of 11 |
| wide band-pass, $ f_2/f_1 \ge 200 $, degree 10–16 (10) | 1260–432 432 | 5–8 | 2–3 | 0 of 10 |

"First one found" is the first ladder with an admissible load in the reference's depth-first order (version 4.1.0 order, $ +F $ first).
For the 13 full Cauer filters only the minimum is computed (6 to 11 inductors); listing their ladders one by one is not practical.

**Ablations** (`tests/truth_generator/ablate.py`):

594/594 behave as expected, with the verifier of version 4.2.0. The script was run in three parts, so the format errors are checked on three specifications.
- For each of the 67 specifications (the two zero-placement ablations for the 65 with finite zeros):
  - these pass:
    - the reference program's output, through the sandboxed runner;
    - the stored witness (a minimum-inductor ladder);
    - the witness with every value scaled by $ 1 + 10^{-12} $;
  - these fail:
    - one inductor scaled by $ 1 + 10^{-4} $: $ 8.7 \times 10^{-4} $ to 8.3 dB;
    - the reference function with its reflection zeros moved by $ 10^{-6} $: $ 1.6 \times 10^{-5} $ to $ 9.2 \times 10^{-2} $ dB;
    - an exact ladder for the transmission zeros moved by $ 10^{-7} $: $ 4.4 \times 10^{-5} $ to $ 2.2 \times 10^{-2} $ dB;
    - an exact ladder for the starting-guess zeros, which are spaced out from the edge but not equiripple: 1.2 to 156 dB.
- For 61 specifications, a canonical ladder with exactly the right response and an admissible load, but more inductors than the minimum, fails (requirement 4). The one tested is the admissible ladder with the fewest inductors above the minimum, usually one more. The six without this ablation are LP-6M and BP-CB, whose admissible ladders all have the same number of inductors, and four of degree 14 to 16 whose ladders are only partly listed (the first 3000 per sign), with none above the minimum among them.
- For 38 band-pass specifications, a canonical ladder with exactly the right response but a load outside $ [R_s/2, 2R_s] $ fails. The loads tested range from $ 10^{-27} R_s $ to $ 2 \times 10^{26} R_s $.
- 10 format errors are rejected (on each of the three specifications): an extra key, a missing load, negative, boolean and string values, an unknown arm type, an extra arm key, an empty ladder, NaN, and a symlink.
- Shortcuts:
  - the first ladder of a depth-first extraction passes for 9 of the 67 (LP-2M, LP-6M, LP-27M12, LP-868, LP-ADC50, LP-DAC2G5, LP-DOCSIS85, LP-VIDEO, LP-ADC250) and fails for the other 58, on the load or on the number of inductors;
  - the first ladder with an admissible load has the fewest inductors for 12 of the 65 specifications whose listed ladders include an admissible one;
  - zero shifting only from infinity finds no ladder for 23 of the 37 specifications of degree up to 12 (the band-pass ones with a finite zero below the passband).
- The edge check: the stored witness of BP-IF21M4N with one resonator capacitor scaled by $ 1 + 10^{-10} $ is within $ 2.8 \times 10^{-6} $ dB on the grid but off by $ 1.8 \times 10^{-5} $ dB at its stopband edge, which lies inside the window around a zero. It passed the verifier of version 4.0.0, which dropped such edges, and fails that of 4.1.0.

**Second implementation and model runs.**

- A second, separately written program (`tests/truth_generator/alt_synthesize.py`, working precision set by `ALT_DPS`) solves the task by a different route:
  - the passband polynomial comes from the closed-form generalized-Chebyshev phase sum in $ \omega^2 $, $ \theta = N_\infty \arccos x + \sum_k \arccos\frac{1 - x x_k}{x - x_k} $, solved by bisection, instead of a Remez exchange;
  - Newton's method places the zeros;
  - a depth-first extraction with partial removals at $ 0 $ and $ \infty $ builds the ladder. It keeps its own count of the zeros at dc and at infinity still to be realized. Since version 4.2.0 it is an exhaustive branch and bound for the fewest inductors, like the reference, but separately written. Version 4.2.0 also fixes its handling of the sign $ -F $: it had truncated the wrong one of $ E \pm F $, so it only ever searched the ladders of $ +F $. That did not matter while any admissible ladder was enough, but the minimum-inductor ladder of LP-POTS, for example, comes from $ -F $.

  With the verifier of version 4.2.0 it passes all 67 graded specifications at 120 digits. It finds the same minimum number of inductors as the reference for every one (`check_min_inductors.py`), in up to 334 s per specification (three jobs at once on a 4-core machine; it is a check, not bound by the time limit). The precision study below was made with the earlier, first-ladder version on the specifications of version 4.1.0.

  At 120 digits (its default) it passes all 6 examples and all 54 graded specifications (with the verifier of version 4.1.0; largest error $ 8.1 \times 10^{-11} $ dB). On the 24 specifications of degree up to 12 it passes at 60 digits, fails BP-15M, BP-455ADJ and BP-DCF77 at 50, and also BP-10M, BP-2M and BP-UMTS at 40. On the nine narrow IF filters it fails four at 60 digits (BP-CW455, BP-IF10M7C, BP-IF21M4N, BP-IF70M5) and passes all nine at 80 and at 100. Each failure is "no ladder found": the zero tests of the extraction no longer resolve the cancellation. The tests therefore accept a correct alternative method, and reject the same method when it is carried out with too few digits.
- Calibration runs with the terminus-2 agent on version 2.0.4 (19 specifications), using the same verifier:
  - GPT-6 solved the approximation correctly but extended its extraction search by trial and error. The search took more than 300 s on 8 band-pass specifications and found no ladder for BP-CB. Reward 0.
  - Opus 5.5 passed with the 40-digit version of the program above. The six receivers with one zero at dc and two zeros below the band were added in version 2.1.0 because that precision does not realize them.
- Calibration runs on version 2.1.2 (24 graded specifications):
  - GPT-6 placed the zeros correctly, but its extraction search ran past 300 s on the five receivers with one zero at dc and two zeros below the band (BP-10M, BP-15M, BP-455ADJ, BP-DCF77, BP-UMTS). Reward 0.
  - Opus 5.5 passed (reward 1). Its program used the same method as above at a fixed 50 digits, with fixed tolerances ($ 10^{-25} $ for its zero tests and $ 10^{-30} $ for trimming coefficients).
- Version 3.0.0 adds the nine narrow IF filters of degree 14 and 16, which need 30 to 46 digits of cancellation (difficulty 4). That Opus program, rebuilt from its trajectory and run with the verifier of this version:
  - passes the 24 specifications of version 2.1.2;
  - fails all nine narrow IF filters at its 50 digits ("synthesis failed"), so its reward is 0;
  - with its precision raised to 60 digits it still fails six of them, at 70 digits one (BP-IF21M4N), and at 80 digits none.

  So the program is correct apart from its fixed precision: choosing the precision for the filter in hand, as the reference does, is what these specifications test.
- Calibration run on version 3.0.0 (34 tests):
  - GPT-6 failed 14 of the 34 tests; its extraction search found no ladder for BP-10M and others. Reward 0.
  - Opus 5.5 passed (reward 1). Its program sets the precision of the extraction from $ n $ and $ b $ ($ 50 + 2.5\,n\log_{10}(1/b) + 2n $ digits), but runs the Remez exchange in double precision, in Chebyshev polynomials of $ x $ on $ [u_1, 1] $, and places the zeros with scipy's fsolve on top of it.
- Version 4.0.0 adds the ten wide band-pass filters (difficulty 5). That Opus program, rebuilt from its trajectory and run with the verifier of this version:
  - passes the 33 earlier specifications;
  - fails all ten wide filters. On six it stops with an error (fsolve gets no usable Remez solution, or the extraction finds no ladder for the wrong function); on four it writes a ladder whose loss is off by $ 5.8 \times 10^{-3} $ dB (BP-CISPRB), 0.044 dB (BP-SONAR), 8.1 dB (BP-30K30M) and 89 dB (BP-10K10M).

  The reference and the second implementation both pass the ten (errors near $ 10^{-13} $ dB), because their exchange works in extended precision.
- Calibration run on version 4.0.0 (44 tests): neither model passed (reward 0 for both).
- Version 4.1.0 changes the grading in two ways, both after review of version 4.0.0:
  - the verifier always checks the stopband edges. Version 4.0.0 dropped an edge that fell within $ w $ of a zero (40 of the 73 edges, in 22 specifications), so a ladder could be off at an edge without failing (see the edge check under the ablations);
  - 11 specifications are added so that every part of the stated scope is graded: low-pass filters of degree 10 to 16 (LP-6M, LP-AUDIO96, LP-ADC50, LP-HF30, LP-DOCSIS85, LP-ADC122, LP-DAC2G5), among them one with no finite zeros and one with a single zero; band-pass filters with a single finite zero (BP-49M of degree 6, BP-ADSB with only a lower stopband, BP-GPSL1 with only an upper one); and the narrowest band of the scope, a 500 Hz CW filter at 455 kHz (BP-CW455N, 0.11 %).

  Both changes only add requirements; a program that passes version 4.1.0 also passes 4.0.0. The reference and the second implementation pass all 54.
- Calibration run on version 4.1.0 (55 tests):
  - GPT-6 failed 14 of the 55 tests (BP-10M, BP-15M, BP-2M, BP-455ADJ, BP-49M, BP-ADSB, BP-CW455, BP-DCF77, BP-GPSL1, BP-IF21M4, BP-IF45M, BP-IF70M2, BP-MW, BP-UMTS), each with "no canonical ladder found": its extraction search missed the ladders with an admissible load. Reward 0.
  - Opus 5.5 passed (reward 1, 55/55). Its program has these parts:
    - the zeros come from the closed-form phase sum in $ (f/f_0)^2 $, with fsolve and then `mp.findroot` at 50 digits on the stopband conditions;
    - the spectral factorization and the extraction run at a fixed 110 digits;
    - a depth-first extraction (full removals, a flip between $ Z $ and $ Y $, and zero shifting by a partial removal at $ \infty $ or $ 0 $ chosen by the sign of $ \operatorname{Im} Z(j\omega_k) $) returns the first ladder with an admissible load.

    Rebuilt from its trajectory, that program also passes 70 further specifications of the version 4.1.0 scope, generated at its edges (degree 12 to 16, 0.1 % to 0.3 % and $ f_2/f_1 $ up to 3500, two or three zeros, close stopband edges), apart from three that have no ladder with an admissible load at all. Its zeros agree with an independent `make_truth.py` solve on all of them. Within that scope the task no longer separated a complete implementation from the reference.
- Version 4.2.0 adds requirement 4 (the fewest inductors) and widens the scope to as many finite zeros as the degree allows (up to 7). It adds 13 full Cauer specifications, and the verifier timeout is raised to cover 300 s for every specification. This tests the choice of realization (difficulty 3) as an optimization over the whole extraction tree, not just a search for the first ladder that works:
  - the Opus program of the 4.1.0 run returns, for 39 of the 54 earlier specifications, a ladder with 1 to 4 more inductors than the minimum (it was not asked for the minimum then);
  - made exhaustive (keeping the admissible ladder with the fewest inductors, and returning the best one found after 280 s), that same depth-first search finds the minimum for the 54 earlier specifications and for the three full Cauer low-pass filters. For all ten full Cauer band-pass filters (6 or 7 finite zeros, degree 16) it has no admissible ladder after 280 s (seven of them) or not the minimal one (BP-DAB, BP-LFMF, BP-UHFTV). With the obvious bound added (cut a branch whose inductor count already reaches the best one) it runs past 300 s on all ten, and on four of the 54 earlier ones (BP-10K10M, BP-IF10M7N, BP-IF21M4, BP-IF70M2);
  - with both pruning steps of difficulty 3 (partial removals only of poles that exist, by less than their residue; branch and bound), the same program passes all 67 graded specifications, in at most 195 s per specification (on a 4-core machine with four jobs running). So the requirement is met by a complete, pruned search, which is what the reference does.

**Verifier runs (version 4.2.0):** the verifier image built from `tests/Dockerfile` (Python 3.11.18 with the pinned packages), run with `tests/test.sh` without network, on 1 CPU and with 4 GB:

- with the reference program in `/root/results`: reward 1, 68/68 tests passed, 1822 s;
- with no program: reward 0, 68/68 tests failed.

**Harbor runs (version 4.0.0):**

- `harbor run -a oracle`: reward 1.0, 44/44 tests, 13 min 2 s. The verifier, which runs the program on the 43 specifications without network and on one CPU core, takes 705 s.
- `harbor run -a nop`: reward 0, 44/44 tests failed.

## References

1. S. Darlington, "Synthesis of reactance 4-poles which produce prescribed insertion loss characteristics," *J. Math. Phys.* 18 (1939) 257–353.
2. W. Cauer, *Synthesis of Linear Communication Networks*, McGraw-Hill, 1958.
3. R. Saal and E. Ulbrich, "On the design of filters by synthesis," *IRE Trans. Circuit Theory* CT-5 (1958) 284–327.
4. G. C. Temes and J. W. LaPatra, *Introduction to Circuit Synthesis and Design*, McGraw-Hill, 1977.
5. A. I. Zverev, *Handbook of Filter Synthesis*, Wiley, 1967.
6. R. J. Cameron, "General coupling matrix synthesis methods for Chebyshev filtering functions," *IEEE Trans. Microwave Theory Tech.* 47 (1999) 433–442.
7. E. W. Cheney, *Introduction to Approximation Theory*, McGraw-Hill, 1966 (alternation theorem, Remez algorithm).
8. H. J. Orchard and G. C. Temes, "Filter design using transformed variables," *IEEE Trans. Circuit Theory* CT-15 (1968) 385–408.
