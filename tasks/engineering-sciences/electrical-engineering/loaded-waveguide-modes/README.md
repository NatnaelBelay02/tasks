# loaded-waveguide-modes

Mode analysis of 18 rectangular waveguides loaded, layer by layer across their width, with microwave laminates, argon discharge cells, magnetised YIG slabs and magnetised argon discharges. Every material is lossy and taken from a data sheet or a stated discharge.

For each guide, each of its four test frequencies and each direction of propagation, the agent must report every mode whose attenuation is at most 200 dB/m:
- its propagation constant $ \gamma = \alpha + j\beta $ ($ \alpha > 0 $, fields $ \propto e^{\mp\gamma z} $ for the forward and reverse direction);
- its **power split**: the fraction of the mode's time-averaged power flow that passes through each layer, $ p_i = P_i / P $ with $ P_i = \operatorname{Re}\iint_{\text{layer } i} (\mathbf E \times \mathbf H^*)\cdot\hat z\,dx\,dy $. The fractions add up to 1 and are negative in a layer where the power flows backwards.

This is the transverse-resonance analysis of slab-loaded, plasma-loaded, ferrite-loaded and magnetoplasma-loaded guides [1–3, 14, 16–18]. It is used to design dielectric-slab and ferrite phase shifters, resonance isolators, field-displacement devices, plasma switches, limiters and phase shifters [7, 8, 14].

The instruction states the physics, not the answer's equations:
- Maxwell's equations with perfectly conducting walls;
- the four material models, including the Polder tensor of a y-biased ferrite and the gyroelectric tensor of a y-magnetised plasma, with their sign conventions;
- the time convention, the branch of $ \gamma $, the two directions, and the definition of the power split.

The agent derives the transverse problems itself: the LSE/LSM split and the interface conditions of the isotropic guides; for the gyrotropic guides the nonreciprocal TE problem ($ n = 0 $) and the coupled four-component hybrid problem of every family $ n \ge 1 $; and the transverse fields ($ E_x $, $ H_x $ follow algebraically from the tangential state) needed for the Poynting flux. It then finds every mode and integrates its power flow layer by layer.

**Versions.**
- **v1 (lossless).** The review of v1.2 found that its hardest answers (near-degenerate pairs, very slow surface waves) existed only without loss, that the plasma layers matched no real material, and that there were no references.
- **v2.0 (lossy).** Every layer became a cited material at its real loss. The instruction still gave the two dispersion functions explicitly. GPT-6 solved it (see "Model runs") by complex root-finding on those functions, which shows that, with the equations given, the task reduced to numerics.
- **v3.0.** The dispersion functions are no longer given, and four nonreciprocal ferrite guides are added (y-independent TE modes only). Their forward and reverse modes differ, so the usual $ \pm\gamma $ pairing and any search in $ q = k_n^2 - \gamma^2 $ do not apply to them. Calibration: GPT-6 failed, Opus 5.5 (low) passed (see "Model runs").
- **v4.0 (this version).** Two changes, both aimed at what the passing v3.0 run did not have to do.
  - **Hybrid modes of the gyrotropic guides.** The ferrite guides now ask for every mode, not only the y-independent ones: for each $ n \ge 1 $ the family HYB$ _n $, which is neither LSE nor LSM. It needs a four-component (e_y, e_z, h_y, h_z) transfer problem whose matrix is quadratic in $ \gamma $ and contains both off-diagonal tensor terms. Two guides with a **magnetised argon discharge** (gyroelectric tensor, G17, G18) are added. Their HYB modes are strongly nonreciprocal, and G18 carries one-way surface magnetoplasmons with no reverse partner under the limit.
  - **Power split.** Every reported mode also carries its per-layer power fractions, to $ 10^{-6} \sum|p^*| $. A mode list can no longer be right by root-finding alone: the fields, including the algebraic $ E_x, H_x $ of the hybrid modes and the conjugate in $ \mathbf E \times \mathbf H^* $, must be reconstructed and integrated correctly.

## The guides

| id | guide | loading |
|---|---|---|
| G01 | WR-90 | dielectric-slab phase shifter: 3 mm RT/duroid 6002 against a side wall |
| G02 | WR-90 | centred H-plane slab: 2.54 mm RO3010 |
| G03 | WR-90 | dual-slab phase shifter: 1.27 mm RO3006 against both side walls |
| G04 | WR-90 | five 1.524 mm RO4350B sheets with 2.286 mm air gaps |
| G05 | WR-62 | 6.35 mm PTFE against a side wall |
| G06 | WR-137 | 2 mm RO3010 against both side walls |
| G07 | WR-90 | plasma switch: 4 mm argon cell between 1 mm quartz plates, centred |
| G08 | WR-137 | plasma phase shifter: 7 mm argon cell between 1.5 mm quartz plates, centred |
| G09 | WR-284 | S-band plasma limiter: 9 mm argon cell between 1.5 mm quartz plates, centred |
| G10 | WR-187 | 6 mm argon discharge behind a 1.5 mm quartz window, then 4 mm RO3006 |
| G11 | WR-112 | low-pressure cell: 4.5 mm argon between 1.5 mm quartz plates, centred |
| G12 | WR-90 | plasma surface-wave guide: 4 mm RO3010, then a 4 mm argon cell between 1 mm quartz plates |
| G13 | WR-90 | dielectric-loaded ferrite phase shifter: 2 mm YIG (800 Oe, +y) 3 mm from a side wall, backed by 2.54 mm RO3010 |
| G14 | WR-90 | twin-slab ferrite phase shifter: two 2 mm YIG slabs magnetised in opposite directions (400 Oe) |
| G15 | WR-90 | resonance isolator: 1 mm YIG near the quarter-width point, biased near ferrimagnetic resonance (3200 Oe) |
| G16 | WR-137 | C-band ferrite phase shifter: 3 mm YIG (300 Oe, +y) 5 mm from a side wall, backed by 3 mm RO3006 |
| G17 | WR-90 | magnetoplasma surface-wave guide: 4 mm RO3010, then a 4 mm magnetised argon cell (1500 G, +y) between 1 mm quartz plates |
| G18 | WR-137 | magnetoplasma phase shifter: 6 mm magnetised argon discharge (1000 G, −y) against a side wall, behind a 1.5 mm quartz window, then 3 mm RO3006 |

Waveguide inner dimensions are the EIA standard ones [9]. Thicknesses are stock laminate, plate and ferrite-slab thicknesses or the remainder of the guide width. Test frequencies are round values, chosen so that each guide's four frequencies sample different regimes of its spectrum. No number is tuned to a root.

**Dielectrics** (relative permittivity and loss tangent at 10 GHz, used at every test frequency):

| material | $ \varepsilon' $ | $ \tan\delta $ | source |
|---|---|---|---|
| RT/duroid 6002 | 2.94 | 0.0012 | Rogers data sheet [4] |
| RO4350B | 3.48 | 0.0037 | Rogers data sheet [4] |
| RO3006 | 6.15 | 0.0020 | Rogers data sheet [4] |
| RO3010 | 10.2 | 0.0022 | Rogers data sheet [4] |
| PTFE | 2.08 | 0.0004 | Pozar, table of dielectric constants [3] |
| fused quartz | 3.78 | 0.0001 | Pozar, table of dielectric constants [3] |
| air | 1 | 0 | |

**Ferrite.** Polycrystalline yttrium iron garnet, Trans-Tech G-113 [15], as listed in Pozar's table of ferrite properties [3]:
- $ 4\pi M_s = 1780 $ G, linewidth $ \Delta H = 45 $ Oe, $ \varepsilon' = 15.0 $, $ \tan\delta = 0.0002 $.
- The bias fields are internal fields: 300–800 Oe for the phase shifters, which operate above the ferrimagnetic resonance and the magnetostatic surface-wave band, and 3200 Oe for the resonance isolator ($ f_0 = 8.96 $ GHz).
- The Polder tensor with the loss term $ f_0 \to f_0 + j\gamma_g\Delta H/2 $ and $ \gamma_g = 2.8 $ MHz/Oe is the standard form [3, 14].
- In a guide with a ferrite layer the modes that vary with $ y $ are hybrid [14, 16]: with the bias along $ y $ the problem stays separable in $ y $ (fields $ \propto \sin, \cos(n\pi y/b) $), but each family $ n \ge 1 $ couples all six field components. A bias along $ x $ or $ z $ would couple different $ n $; the guides use the standard transverse ($ y $) bias of slab phase shifters and isolators.

**Argon discharges.** Each plasma layer is a planar discharge cell: argon between two fused-quartz plates (or behind one quartz window), with uniform electron density. This is the slab model used in analyses of plasma-loaded guides [7, 8].
- **Permittivity.** The collisional Drude form $ \varepsilon = 1 - \omega_p^2 / (\omega(\omega - j\nu)) $ for $ e^{+j\omega t} $ [5].
- **Plasma frequency.** It follows from the electron density, $ f_p = 8.98\,{\rm kHz}\,\sqrt{n_e / {\rm cm^{-3}}} $ [6].
- **Collision frequency.** $ \nu = n_g K_m $, with $ K_m \approx 10^{-13} $ m$ ^3 $/s for argon at a few eV [5] and $ n_g $ the gas density at 300 K.

| guide | $ n_e $ (m$ ^{-3} $) | $ f_p $ (GHz) | $ \nu $ (s$ ^{-1} $) | argon pressure | $ \nu/\omega $ at the test frequencies |
|---|---|---|---|---|---|
| G07 | $ 2.4 \times 10^{18} $ | 14 | $ 5 \times 10^8 $ | 21 Pa (0.16 Torr) | 0.007–0.011 |
| G08 | $ 1.0 \times 10^{18} $ | 9 | $ 2 \times 10^8 $ | 8 Pa (0.06 Torr) | 0.004–0.007 |
| G09 | $ 2.5 \times 10^{17} $ | 4.5 | $ 8 \times 10^8 $ | 33 Pa (0.25 Torr) | 0.036–0.049 |
| G10 | $ 6.1 \times 10^{17} $ | 7 | $ 3 \times 10^8 $ | 12 Pa (0.09 Torr) | 0.008–0.014 |
| G11 | $ 1.24 \times 10^{18} $ | 10 | $ 1.5 \times 10^8 $ | 6 Pa (0.05 Torr) | 0.003–0.004 |
| G12 | $ 2.1 \times 10^{18} $ | 13 | $ 4 \times 10^8 $ | 17 Pa (0.12 Torr) | 0.006–0.013 |

These densities and pressures are those of low-pressure glow, RF and microwave discharges [5].

**Magnetised argon discharges** (G17, G18). The same cells in a static field $ B_0 $ along $ \pm y $, as in magnetised-plasma phase shifters and in the analysis of surface magnetoplasmons [17, 18]. The cold-electron equation of motion $ m\dot{\mathbf v} = -e(\mathbf E + \mathbf v\times\mathbf B_0) - m\nu\mathbf v $ gives, for $ e^{+j\omega t} $ and $ \mathbf B_0 = s B_0\hat y $, the tensor of the instruction: $ \varepsilon_1 = 1 - f_p^2 f_\nu / (f(f_\nu^2 - f_c^2)) $, $ \varepsilon_3 = 1 - f_p^2/(f f_\nu) $ along the field, and the gyrotropic term $ \varepsilon_{xz} = +j s g $, $ g = f_p^2 f_c / (f(f_\nu^2 - f_c^2)) $, $ f_\nu = f - j\nu/2\pi $ [17]. The cyclotron frequency is $ f_c = 2.799\,{\rm MHz/G} \times B_0 $ [6].

| guide | $ n_e $ (m$ ^{-3} $) | $ f_p $ (GHz) | $ \nu $ (s$ ^{-1} $) | argon pressure | $ B_0 $ | $ f_c $ (GHz) |
|---|---|---|---|---|---|---|
| G17 | $ 1.5 \times 10^{18} $ | 11 | $ 3 \times 10^8 $ | 12 Pa | 1500 G, $ +y $ | 4.2 |
| G18 | $ 6.1 \times 10^{17} $ | 7 | $ 2 \times 10^8 $ | 8 Pa | 1000 G, $ -y $ | 2.8 |

The test frequencies lie above the electron cyclotron resonance and away from the upper-hybrid resonance $ \sqrt{f_p^2 + f_c^2} $ (11.8 and 7.5 GHz), so $ \varepsilon_1 $, $ \varepsilon_3 $ and $ g $ are all of order one and $ \varepsilon_1 $ changes sign across the band.

**Walls.** The walls are perfect conductors, as in the textbook analysis. A finite wall conductivity couples LSE and LSM at the top and bottom walls. Copper walls add about 0.1 dB/m in WR-90 at X band [3], outside the model and far below the 4 dB/m (2 %) margin that separates every mode from the 200 dB/m limit.

## Difficulty

- **The transverse problems must be derived.**
  - For the isotropic guides: the LSE/LSM decomposition with respect to the layering direction, and the interface conditions. LSE needs $ E_y, E_z $ continuous ($ u, u' $); LSM needs $ H_y, H_z $ continuous, so $ v $ and $ v'/\varepsilon $ are continuous.
  - The relation $ \gamma^2 = k_n^2 - q $ links the transverse eigenvalue to the propagation constant.
  - A common slip is to carry the LSE interface condition over to LSM. The surface waves and backward waves of the plasma guides exist only because of the $ 1/\varepsilon $ in that condition.
- **The hybrid families must be derived** (G13–G18, families HYB$ _n $, $ n \ge 1 $).
  - With $ E_x, E_z, H_y \propto \sin k_y y $ and $ E_y, H_x, H_z \propto \cos k_y y $, the $ x $ components of the two curl equations are algebraic: $ e_x = (\gamma h_y - k_y h_z - k_0 e_2 e_z)/(jk_0\varepsilon_1) $ and $ h_x = -(k_y e_z + \gamma e_y + k_0\kappa h_z)/(jk_0\mu) $ (normalised $ h = \eta_0 H $, $ [\varepsilon]_{xz} = -je_2 $).
  - The other four give $ \psi' = M(\gamma)\psi $ for $ \psi = (e_y, e_z, h_y, h_z) $, continuous at every interface, with $ e_y = e_z = 0 $ on both side walls. A mode is a zero of the $ 2\times2 $ determinant of the $ e $-rows of $ T(a)[0; I] $.
  - $ M $ is quadratic in $ \gamma $ and has terms odd in $ \gamma $ wherever $ \kappa \ne 0 $ or $ e_2 \ne 0 $. Dropping either off-diagonal term, or the coupling between the two polarisations, gives reciprocal LSE/LSM-like modes and fails.
  - In an isotropic layer $ M^2 = -K^2 I $; in a gyrotropic layer the four exponents are distinct and $ \exp(Md) $ has no simple closed form, and it is stiff: $ \|M d\| $ reaches 30–40 times its spectral radius for $ |\gamma| \gg k_0 $.
  - The family index has no a-priori end: the search must go on in $ n $ until no family can have $ \alpha \le A $.
- **The power split must be computed from the fields** (every mode, 1e-6).
  - It needs the full transverse field, including the algebraic $ e_x $, $ h_x $ of the hybrid modes, the null vector of the $ 2\times2 $ wall condition, and $ \operatorname{Re}\,(e_x h_y^* - e_y h_x^*) $ integrated layer by layer.
  - The textbook shortcuts fail: the TE expression $ -e_y h_x^* $ for every mode; dropping the conjugate; the fractions of the reverse mode taken at $ +\gamma_r $; $ |P_i|/\sum|P_k| $, which loses the sign of the backward power in a plasma layer; and power proportional to the layer width.
  - In the plasma guides the fractions are large and of both signs (backward power in the plasma, forward power in the dielectrics), e.g. {{EX_POWER}}.
- **Gyrotropic guides are nonreciprocal** (G13–G18, 48 of the 144 direction-spectra).
  - With the Polder tensor, the y-independent TE problem has $ E'' + (k_0^2\varepsilon\mu_e + \gamma^2)E = 0 $ in each layer, $ \mu_e = (\mu^2 - \kappa^2)/\mu $. The continuous field is $ H_z \propto (E' - j(\kappa/\mu)\gamma E)/\mu_e $.
  - The term linear in $ \gamma $ makes the dispersion relation depend on $ \gamma $ itself, not on $ \gamma^2 $: the forward and reverse modes differ.
  - G13 at 8 GHz has its slab mode at $ \beta = 452.5 $ rad/m forward (8.0 dB/m) against 344.1 rad/m reverse (19.0 dB/m). G16 at 5 GHz has 228.0 against 157.0 rad/m. G14 (twin slabs, opposite bias) differs by 1–16 rad/m.
  - The resonance isolator G15 at 11 GHz has a forward mode at 166 dB/m and no reverse mode under the limit. At 11.5 GHz the two directions have different mode sets.
  - Assuming reciprocity, or getting the sign of $ \kappa $ or of $ g $ relative to the propagation direction wrong (it swaps the two directions), fails these guides. So does a search in $ q = -\gamma^2 $, which cannot represent the asymmetry.
  - The hybrid families are more nonreciprocal than LSE$ _0 $: {{EX_NONRECIP}}.
  - **One-way surface magnetoplasmons** (G18). {{EX_ONEWAY}}
- **The modes are complex.** With lossy layers the dispersion functions have no zeros on the real axis, so scanning for sign changes finds nothing. The search is over a region of the complex plane: the parabola $ \operatorname{Re} q \ge k_n^2 - A^2 + (\operatorname{Im} q)^2 / 4A^2 $, or for ferrite guides the strip $ |\operatorname{Re}\gamma| \le A $.
- **Backward waves** (G07, G08, G09, G11). Below the surface-plasma resonance of a plasma column the LSM surface wave carries power and phase in opposite directions [7, 8]. With $ \alpha > 0 $ it has $ \beta < 0 $: G08 at 4.4 GHz has LSM1 at $ -539 $ rad/m (132 dB/m) and LSM2 at $ -415 $ rad/m (172 dB/m).
- **Evanescent modes inside the limit** (G03, G06, G08, G09, G10). Just below a cutoff a mode decays with $ \beta \approx 0 $ but can still be under 200 dB/m: G06 at 4.15 GHz has LSE0 at 168 dB/m.
- **Slow surface waves** (G12 at 5 and 5.15 GHz, $ \operatorname{Re} q \approx 41 $–$ 47\,k_0^2 $). A search bounded by the dielectric permittivities misses them.
- **Loss is not a perturbation here.** A lossless analysis followed by first-order loss perturbation fails every guide, and near cutoffs it changes the mode set.
- **Bookkeeping.** There are {{NMODES}} reported modes in 144 direction-spectra, with up to {{MAXFAM}} families and {{MAXMODES}} modes per spectrum, and empty spectra. LSM$ _0 $ and HYB$ _0 $ do not exist, but LSE$ _0 $ does; each mode carries one power fraction per layer ({{NFRAC}} numbers in all).

**The answers are those of the lossy guide, not of an ideal one.** Halving all losses (dielectric, plasma) changes 4 of the 48 isotropic spectra, and doubling them changes 7. Only the most attenuated backward and surface waves cross the limit, as their attenuation scales with the loss. The backward waves of G08 at 4.8 GHz (61 dB/m) and G11 (37–45 dB/m) persist.

A single wrong item in any of the 18 guides fails the task.

## Reference solution

`solution/solve.py` (dispersion functions and the isotropic search) and `solution/hybrid.py` (gyrotropic families and the power split) use numpy, scipy and mpmath.

**Isotropic guides (G01–G12).** Per test frequency and kind:
1. **Counting.** The zeros of $ D_{\rm LSE}(q) $ and $ D_{\rm LSM}(q) $ (transfer matrices, $ q = k_n^2 - \gamma^2 $) are counted with the argument principle [10] on rectangles covering the region $ \alpha \le 1.02 A $ of the $ q $ plane ($ A = 200 $ dB/m $ = 10\ln 10 $ Np/m).
   - $ \arg D $ is followed with steps limited to 0.4 rad.
   - The state is renormalised by positive real factors, which leave $ \arg D $ unchanged.
2. **Locating.** Each single zero is located by the contour moment and Newton's method on $ D'/D $, then polished in mpmath at 50 digits.
3. **Bounds.** The Rayleigh bound for LSE. For LSM, beyond every interface surface-wave estimate and four times past the last zero.
4. **Reporting.** Each zero gives $ \gamma = \sqrt{k_n^2 - q} $ ($ \operatorname{Re}\gamma > 0 $) for every family. These guides are reciprocal, so the reverse list equals the forward one.

**Gyrotropic guides (G13–G18).** For $ n = 0, 1, 2, \dots $:
1. **Dispersion function.** $ D_n(\gamma) $: for $ n = 0 $ the TE state $ (e_y, h_z) $, for $ n \ge 1 $ the $ 2 \times 2 $ determinant above. Isotropic layers use $ \cos(Kd) I + \sin(Kd)/K\,M $, gyrotropic layers `scipy.linalg.expm`; the state is renormalised by positive real factors.
2. **Counting.** The argument principle on rectangles covering $ 0 < \operatorname{Re}\gamma \le 1.02A $ (forward) and $ -1.02A \le \operatorname{Re}\gamma < 0 $ (reverse), in 48 strips of $ |\operatorname{Im}\gamma| \le 4 k_0 \sqrt{m} + 2k_y + 200 $, $ m = \max|\varepsilon|\max(|\mu|, |\mu_e|, 1) $.
3. **Locating.** Rectangles are bisected (off-centre) until each holds one zero and is $ 10^{-4} $ (relative) small, then the zero is polished in mpmath at 45 digits.
4. **Stopping in $ n $.** After the first $ n \ge 1 $ without a zero once $ k_y^2 > k_0^2 m + A^2 $.
5. **Reporting.** A reverse zero $ \gamma $ is reported as $ -\gamma $, with the power split of its own field.

**Power split (every guide).** At the polished $ \gamma $ the null vector of the $ 2 \times 2 $ wall condition gives $ \psi(0) $; $ \psi $ is carried through each layer by $ \exp(M t) $ in mpmath (40 digits), $ e_x $ and $ h_x $ follow algebraically, and $ \operatorname{Re}(e_x h_y^* - e_y h_x^*) $ is integrated by 24-point Gauss–Legendre quadrature per layer.

All 18 guides take about {{ORACLE_TIME}} on one core.

## Verification

`tests/test_modes.py` reads only `/root/results/modes.json`, and only as data.

**Reading the file safely.**
- It uses `lstat`, a non-blocking `O_NOFOLLOW` open and a read capped at 16 MiB.
- It rejects NaN, Infinity and repeated keys.
- It checks the whole output contract before any value: guide ids, the single key `modes`, one entry per test frequency with exactly `forward` and `reverse`, family names that fully match `LSE<n>` ($ n \ge 0 $), `LSM<n>` or `HYB<n>` ($ n \ge 1 $), at most 400 families and 1000 modes, and modes `[beta, alpha, power]` of finite numbers with $ \alpha > 0 $ and `power` a list of exactly one number per layer.

**The value tests**, for every guide, test frequency and direction:
- **`test_mode_families`.** The families present and their mode counts equal the reference.
- **`test_propagation_constants`.** The reported $ \gamma $ pair one-to-one with the reference within $ |\gamma - \gamma^*| \le 10^{-7}|\gamma^*| $.
- **`test_power_split`.** For each pair, every fraction satisfies $ |p_i - p_i^*| \le 10^{-6} \sum_k |p_k^*| $.

**Why the tolerance is what it is.**
- **Correct solvers pass easily.** Both correct solvers reproduce the certified values to about $ 10^{-12} $. Every $ \gamma $ scaled by $ 1 + 5\times10^{-8} $ passes; one off by $ 3 \times 10^{-7} $ fails.
- **The textbook approximation fails.** The tolerance is coarse enough for any double-precision solution of the stated model, and fine enough to reject first-order loss perturbation. That approximation is off by $ 1.3 \times 10^{-7} $ to 0.3 for the laminate guides and up to 1 for the plasma guides.
- **The nonreciprocity is far above the tolerance.** The smallest forward–reverse difference, in G14 at 12 GHz, is 0.4 rad/m out of 626, $ 6 \times 10^{-4} $ relative, which is 6000 tolerances.
- **Pairing is unambiguous.** Modes of a family and direction are at least $ 10^{-4}|\gamma| $ apart (certified, stated in the instruction).
- **Membership never depends on rounding.** No mode is within 2 % of the limit (certified, stated); limits of 200 and 220 dB/m give the same answer.
- **The power tolerance.** The certified fractions are enclosed to better than $ 10^{-9} $, and both correct solvers reproduce them to {{POWER_AGREE}}. Every fraction scaled by $ 1 + 3\times10^{-7} $ passes; one fraction off by $ 3\times10^{-6}\sum|p^*| $ fails. Scaling by $ \sum|p^*| $ keeps the test meaningful where the fractions are large and of opposite sign; the instruction guarantees $ \sum|P_i| \le 30|P| $ (largest certified value {{MAXCOND}}).

**Where the guides come from** (`tests/truth_generator/gen_guides.py`). The script writes `guides.json` exactly as shipped; the agent copy `environment/data/guides.json` is identical to `tests/data/guides.json`. Its docstring repeats the sources above.

**Where the reference values come from** (`make_truth.py`, `certify2.py`, `wgacb.py`, `hybacb.py`, `power.py`). The generator is separate from `solution/` and never uses the argument principle. It certifies everything with ball arithmetic (Arb, through python-flint, 256 bits) [11].

- **Evaluation.**
  - $ D $ and its derivatives are carried as dual numbers with complex-ball coefficients.
  - The layer functions $ \cosh(\kappa d) $ and $ \sinh(\kappa d)/\kappa $ are evaluated as entire functions: by Taylor series with a remainder bound near 0, and elsewhere with a square-root branch that avoids its cut on the ball.
  - Layers certainly evanescent on a box are divided by $ \cosh(\kappa d) $. That factor is analytic and nonzero, so the zeros are unchanged and the balls stay tight.
- **Isotropic guides.**
  - Branch-and-bound over the region $ \alpha \le 1.02A $ of the $ q $ plane: mean-value exclusion, complex Krawczyk acceptance [12], and refinement to $ 10^{-55} $ boxes.
  - Bounds: the Rayleigh bound for LSE. For LSM, a proof that $ D_{\rm LSM} \ne 0 $ for $ \operatorname{Re} q \ge Q $, by the impedance recursion $ t \mapsto (t + T)/(1 + tT) $, $ T = \tanh\kappa d $, on narrow boxes and an explicit bound for the unbounded rest.
- **Ferrite guides.**
  - The same branch-and-bound in the $ \gamma $ plane over $ |\operatorname{Re}\gamma| \le 1.02A $, $ |\operatorname{Im}\gamma| \le B $.
  - $ B $ ($ 14 $–$ 64\,k_0 $) is proved by the recursion for $ \sigma = \kappa E/E' $: $ \sigma \mapsto (\sigma + T)/(1 + \sigma T) $ across a layer, and $ \sigma \mapsto \sigma/(a + b\sigma) $ across an interface, with $ a = \mu_{eR}\kappa_L/(\mu_{eL}\kappa_R) $ and the nonreciprocal term $ b = j(\tau_R - \tau_L\mu_{eR}/\mu_{eL})\gamma/\kappa_R $.
  - It is evaluated with $ \kappa = \pm(\gamma/j)\sqrt{1 + k_0^2\varepsilon\mu_e/\gamma^2} $ for tight enclosures. $ D = 0 $ only if $ \sigma(a) = 0 $.
- **Hybrid families** (`hybacb.py`).
  - The same branch-and-bound in the $ \gamma $ plane over $ |\operatorname{Re}\gamma| \le 1.02A $, $ |\operatorname{Im}\gamma| \le 4k_0\sqrt m + 2k_y + 400 $, family by family until two consecutive families have no zero in the region past $ k_y^2 > k_0^2 m + A_1^2 $.
  - $ D_n $ is evaluated on the **second compound** (exterior square) of the layer matrices [19]: the 6-vector $ \psi_1 \wedge \psi_2 $ started from $ h_y \wedge h_z $ is carried through the layers and $ D_n $ is its $ e_y \wedge e_z $ component. The two growing states are never formed separately, so their near-parallelism costs nothing.
  - Isotropic layers: $ \wedge^2(C I + S M) = C^2 I + CS\,M^{[2]} + S^2 \wedge^2 M $, with the $ \cosh^2(\kappa d) $ scaling where the layer is certainly evanescent. Gyrotropic layers: $ \exp(M^{[2]}d) $ and its $ \gamma $-derivative by a $ 12 \times 12 $ block exponential taken in the (numerical, exactly inverted) eigenbasis of $ M^{[2]} $ at the box centre, and the mean-value form $ E(c) + E'(\text{box})(\gamma - c) $. All matrices are balanced by exact powers of two [20].
  - Without these steps a box at $ |\gamma| \approx 30k_0 $ in G13 at 8 GHz, $ n = 3 $, overestimates $ D' $ by a factor of 500, and the run did not finish in an hour; with them it takes {{BOXES_G13_3}} boxes.
- **Power split** (`power.py`). The null vector of the certified wall condition, then $ P_i = \operatorname{Re}\int (e_x\, \overline{h_y(\bar t)} - e_y\, \overline{h_x(\bar t)})\,dt $ by `acb.integral` (rigorous Gauss–Legendre with error bounds on the entire integrand) to relative accuracy $ 2^{-90} $.
- **Guarantees.** The generator stops unless, for every mode, $ \operatorname{Re}\gamma $ is certainly nonzero and $ |\alpha - A| > 0.02A $, modes of a family and direction are at least $ 10^{-4}|\gamma| $ apart, $ \sum|P_i|/|P| < 30 $, and every fraction is enclosed to $ 10^{-9} $. All 144 direction-spectra pass, {{TRUTH_STATS}}. `truth_diagnostics.json` records every run.

**Score of an independent correct submission.** `tests/truth_generator/alt_solver.py` shares no code or method with the oracle or the generator.
- **Isotropic guides.** Multi-domain Chebyshev collocation [13] of $ v'' + k_0^2\varepsilon v = qv $.
- **Gyrotropic guides, $ n = 0 $.** Collocation of $ E'' + (k_0^2\varepsilon\mu_e + \gamma^2)E = 0$ with the interface rows linear in $ \gamma $. This is a quadratic eigenproblem, solved by companion linearisation.
- **Gyrotropic guides, $ n \ge 1 $.** Collocation of all six curl equations in the six field components, which is a linear pencil $ A_0 x + \gamma A_1 x = 0 $ (no elimination of $ e_x, h_x $, no transfer matrix).
- **Polishing.** Eigenvalues are polished by Newton's method on determinants of $ 2\times2 $ (and, for the hybrids, $ 4 \times 4 $) matrix exponentials in mpmath.
- **Power split.** From the same six-field collocation at the polished $ \gamma $: three steps of inverse iteration give the fields on the Chebyshev nodes, and the Poynting flux is integrated by Clenshaw–Curtis quadrature.

Scored with the real verifier logic, it passes all 18 guides.

While it was being developed, it missed the even–odd LSM1 pair of G03 at 18 GHz. A strict reproducibility filter between two collocation grids removed it; the filter now only screens.

**Ablations** (`tests/truth_generator/ablate.py`): {{ABL}} behave as expected.
- **These pass:** the oracle; the second implementation; every $ \gamma $ scaled by $ 1 + 5\times10^{-8} $; every power fraction scaled by $ 1 + 3\times10^{-7} $; a 220 dB/m limit.
- **Gyrotropic shortcuts fail:**
  - the hybrid families left out (the v3.0 contract), or only HYB$ _1 $ searched;
  - the off-diagonal tensor terms dropped ($ \kappa = g = 0 $, solved again);
  - the magnetised plasma taken as an isotropic Drude plasma (solved again);
  - gyrotropic guides treated as reciprocal (reverse copied from forward), everywhere or only for the twin-slab G14;
  - the wrong sign of $ \kappa $ and $ g $ (forward and reverse swapped), or of $ g $ only (G17, G18);
  - the isolated reverse direction of G15 at 11 GHz given the forward mode;
  - the ferrite linewidth ignored;
  - the hybrid families named LSE$ _n $; a hybrid family reported for an isotropic guide.
- **Power-split shortcuts fail:**
  - the TE expression $ -e_y h_x^* $ used for every mode, or the conjugate dropped (both solved again);
  - reverse modes given the split of the field $ e^{-\gamma_r z} $, or the forward split of the same family;
  - $ |P_i| / \sum|P_k| $ (unsigned backward power);
  - fractions proportional to the layer widths;
  - one fraction off by $ 3\times10^{-6}\sum|p^*| $.
- **Other physics shortcuts fail:**
  - $ \beta > 0 $ for every mode, or for the backward waves of G11 only;
  - the opposite time convention;
  - evanescent modes below cutoff left out, everywhere or only in G06;
  - surface waves above $ k_0^2\max\varepsilon' $ not searched, or only G12's 5.15 GHz wave left out;
  - one of three LSM1 modes of G04 left out;
  - a mode above the limit reported;
  - LSM2 reported as LSM1;
  - one $ \gamma $ off by $ 3\times10^{-7} $;
  - first-order loss perturbation of the nearly lossless guide (all guides);
  - the limit read as 200 Np/m;
  - a 300 dB/m limit.
- **Format and safety cases fail cleanly:**
  - malformed structure: an extra key, a missing guide, an extra guide key, a missing `reverse` key, an entry without the forward/reverse level, an LSM0 or HYB0 family, an empty list, a mode without its power list, a power list one number short, a missing test-frequency entry;
  - bad values: $ \alpha < 0 $, a family name with a trailing newline or a leading zero, a string or a boolean as a number, NaN, a repeated key;
  - bad files: a symlink, a FIFO, a file over 16 MiB, no file, an empty file;
  - a large file with 400 families of 1000 modes.

**Model runs.**
- **v1.2 (lossless):** both models failed.
- **v2.0 (12 isotropic guides, dispersion functions given in the instruction), calibration run #239175:**
  - **GPT-6 (max) passed** (reward 1, 2 min). It wrote the two given dispersion functions, counted their zeros in a complex $ q $ rectangle with the argument principle, refined them with SciPy, and checked the counts on larger contours.
  - **Opus 5.5 (low) failed** (reward 0). Its contour sampling missed the LSE1 mode of G12 at 9 GHz. Its Newton step, on a function renormalised by $ \max|\cdot| $ and therefore not analytic, converged to a wrong value for LSM1 of G01 at 13 GHz ($ 0.358 + 103.23j $ instead of $ 0.411 + 103.33j $).
  - The GPT-6 run shows that v2.0, with the equations given, tested numerics rather than the physics. v3.0 therefore removes the given dispersion functions, so the transverse problems, including the LSM interface condition, must be derived. It also adds the ferrite guides, whose nonreciprocal TE problem cannot be written as a function of $ q $ at all.
- **v3.0 (ferrite guides added, dispersion functions no longer given):**
  - **GPT-6 failed** (reward 0). It missed one member of the close LSM1 pair of G03 at 16 GHz.
  - **Opus 5.5 (low) passed** (reward 1). It derived the LSE/LSM transfer equations and the nonreciprocal ferrite TE problem with the right sign of $ \kappa $, counted zeros with a tiled argument principle, checked the counts on boxes four times larger, and cross-checked with a finite-element solver.
  - What the passing run shows: every v3.0 answer was a zero of a scalar dispersion function that a careful agent can derive from textbook LSE/LSM and ferrite-slab theory, and the answer was only the list of zeros. v4.0 therefore adds the hybrid families (a coupled four-component problem with both off-diagonal tensor terms, no textbook closed form) and the magnetised-plasma guides, and asks for a field quantity, the power split, that is wrong unless the fields themselves are right.
- **Calibration runs on v4.0:** pending.

**Harbor runs** (v4.0): {{HARBOR}}

## References

1. R. F. Harrington, *Time-Harmonic Electromagnetic Fields*, McGraw-Hill, 1961, ch. 4 (partially filled waveguides, LSE and LSM modes).
2. R. E. Collin, *Field Theory of Guided Waves*, 2nd ed., IEEE Press, 1991, ch. 6 (inhomogeneously filled waveguides, transverse resonance).
3. D. M. Pozar, *Microwave Engineering*, 4th ed., Wiley, 2012:
   - partially loaded rectangular waveguide;
   - attenuation of rectangular waveguide;
   - permeability tensor of a magnetised ferrite with loss;
   - ferrite-loaded rectangular waveguide, phase shifters and resonance isolators;
   - tables of dielectric constants and loss tangents, and of ferrite properties.
4. Rogers Corporation, data sheets *RT/duroid 6002*, *RO4000 Series* (RO4350B), *RO3000 Series* (RO3006, RO3010): $ \varepsilon_r $ and dissipation factor at 10 GHz.
5. M. A. Lieberman and A. J. Lichtenberg, *Principles of Plasma Discharges and Materials Processing*, 2nd ed., Wiley, 2005: plasma dielectric constant with collisions, electron–neutral collision frequency, argon rate constants, densities and pressures of low-pressure discharges.
6. J. D. Huba, *NRL Plasma Formulary*, Naval Research Laboratory, 2019: $ f_{pe} = 8.98 \times 10^3 n_e^{1/2} $ Hz ($ n_e $ in cm$ ^{-3} $).
7. A. W. Trivelpiece and R. W. Gould, "Space charge waves in cylindrical plasma columns," *J. Appl. Phys.* 30 (1959) 1784–1793.
8. A. S. Omar and K. F. Schünemann, "Complex and backward-wave modes in inhomogeneously and anisotropically filled waveguides," *IEEE Trans. Microwave Theory Tech.* 35 (1987) 268–275.
9. EIA RS-261-B, *Rectangular waveguides (WR3 to WR2300)*.
10. L. M. Delves and J. N. Lyness, "A numerical method for locating the zeros of an analytic function," *Math. Comp.* 21 (1967) 543–560.
11. F. Johansson, "Arb: efficient arbitrary-precision midpoint-radius interval arithmetic," *IEEE Trans. Computers* 66 (2017) 1281–1292.
12. R. Krawczyk, "Newton-Algorithmen zur Bestimmung von Nullstellen mit Fehlerschranken," *Computing* 4 (1969) 187–201.
13. L. N. Trefethen, *Spectral Methods in MATLAB*, SIAM, 2000.
14. B. Lax and K. J. Button, *Microwave Ferrites and Ferrimagnetics*, McGraw-Hill, 1962: Polder tensor, ferrite-slab loaded waveguides, nonreciprocal phase shift and resonance isolation.
15. Trans-Tech (Skyworks), *Microwave Magnetic Materials* data sheet, G-113 yttrium iron garnet: $ 4\pi M_s $ 1780 G, $ \Delta H $ 45 Oe, $ \varepsilon' $ 15.0, $ \tan\delta $ 0.0002.
16. F. E. Gardiol, "Anisotropic slabs in rectangular waveguides," *IEEE Trans. Microwave Theory Tech.* 18 (1970) 461–467 (hybrid modes of transversely magnetised ferrite slabs).
17. T. H. Stix, *Waves in Plasmas*, AIP, 1992, ch. 1–2: cold magnetised-plasma dielectric tensor with collisions.
18. Z. Yu, G. Veronis, Z. Wang and S. Fan, "One-way electromagnetic waveguide formed at the interface between a plasmonic metal under a static magnetic field and a photonic crystal," *Phys. Rev. Lett.* 100 (2008) 023902.
19. B. S. Ng and W. H. Reid, "An initial value method for eigenvalue problems using compound matrices," *J. Comput. Phys.* 30 (1979) 125–136.
20. B. N. Parlett and C. Reinsch, "Balancing a matrix for calculation of eigenvalues and eigenvectors," *Numer. Math.* 13 (1969) 293–304.
