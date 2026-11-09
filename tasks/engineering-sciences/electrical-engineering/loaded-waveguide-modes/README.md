# loaded-waveguide-modes

A program that finds the modes of rectangular waveguides loaded, layer by layer across their width, with microwave laminates and ceramics, argon discharge cells, magnetised YIG slabs and films, and magnetised argon discharges. Every material is lossy and taken from a data sheet or a stated discharge. The program is graded on six held-out guides in four guide files; `/root/data/guides.json` gives 18 example guides without their modes.

For each guide, each of its test frequencies and each direction of propagation, the program must report every mode whose attenuation is at most 200 dB/m:
- its propagation constant $ \gamma = \alpha + j\beta $ ($ \alpha > 0 $, fields $ \propto e^{\mp\gamma z} $ for the forward and reverse direction);
- its **power split**: the fraction of the mode's time-averaged power flow that passes through each layer, $ p_i = P_i / P $ with $ P_i = \operatorname{Re}\iint_{\text{layer } i} (\mathbf E \times \mathbf H^*)\cdot\hat z\,dx\,dy $. The fractions add up to 1 and are negative in a layer where the power flows backwards.

This is the transverse-resonance analysis of slab-loaded, plasma-loaded, ferrite-loaded and magnetoplasma-loaded guides [1–3, 14, 16–18], and of shielded magnetostatic-wave lines [21–24]. It is used to design dielectric-slab and ferrite phase shifters, resonance isolators, field-displacement devices, plasma switches, limiters and phase shifters, and magnetostatic-wave delay lines and filters [7, 8, 14, 23].

The instruction states the physics, not the answer's equations:
- Maxwell's equations with perfectly conducting walls;
- the four material models, including the Polder tensor of a y-biased ferrite and the gyroelectric tensor of a y-magnetised plasma, with their sign conventions;
- the time convention, the branch of $ \gamma $, the two directions, and the definition of the power split;
- the ranges within which the held-out guides lie.

The agent derives the transverse problems itself: the LSE/LSM split and the interface conditions of the isotropic guides; for the gyrotropic guides the nonreciprocal TE problem ($ n = 0 $) and the coupled four-component hybrid problem of every family $ n \ge 1 $; and the transverse fields ($ E_x $, $ H_x $ follow algebraically from the tangential state) needed for the Poynting flux. It must also derive **where the modes can be**: how far in $ |\beta| $ and in the family index $ n $ a search has to go before no mode with $ \alpha \le 200 $ dB/m can remain. It then writes a program that finds every mode and integrates its power flow, on guides it has not seen.

**Versions.**
- **v1 (lossless).** The review of v1.2 found that its hardest answers (near-degenerate pairs, very slow surface waves) existed only without loss, that the plasma layers matched no real material, and that there were no references.
- **v2.0 (lossy).** Every layer became a cited material at its real loss. The instruction still gave the two dispersion functions explicitly. GPT-6 solved it (see "Model runs") by complex root-finding on those functions, which shows that, with the equations given, the task reduced to numerics.
- **v3.0.** The dispersion functions are no longer given, and four nonreciprocal ferrite guides are added (y-independent TE modes only). Calibration: GPT-6 failed, Opus 5.5 (low) passed.
- **v4.0.** The hybrid families HYB$ _n $ of the gyrotropic guides, two magnetised-discharge guides, and the per-layer power split of every mode. Calibration: GPT-6 (max) and Opus 5.5 (low) both passed (see "Model runs").
- **v5.0.** The v4.0 runs show how both models passed: each graded guide was in the agent's own input file, so a search with a box sized by the material constants could be checked guide by guide and patched by hand. GPT-6's first pass missed five LSM1 pairs; a contour audit on the visible guides found them, and they were added by hand. v5.0 grades a **program on held-out guides**, which must therefore be right without inspection, and the held-out guides lie in regimes that the example guides do not show and that a material-sized search misses:
  - **thin-film magnetostatic surface waves** (Damon–Eshbach waves of a 100 µm YIG film [21, 22]) at up to 178 $ k_0 $, in every family up to $ n = 12 $;
  - **one-way surface magnetoplasmons** of a magnetised low-pressure discharge, in 14 to 18 consecutive families, backward waves at one frequency;
  - **slow surface plasmons** of a low-pressure discharge behind thick quartz, at up to 12 $ k_0 $;
  - a **ceramic-loaded** guide with up to 7 families per spectrum, and a conventional C-band ferrite phase shifter.

  The model, the output format, the tolerances and the guarantees are those of v4.0. What changes is the deliverable (a program, run four times at once under a wall-clock limit) and the guides it is graded on. Calibration: GPT-6 (max) and Opus 5.5 (low) both failed (see "Model runs").
- **v5.1 (this version).** No change to the guides, the reference values, the tests or the tolerances. The program-size limit of the verifier (50 000 000 bytes) is stated in the instruction as a bullet of its own, and the tolerances are calibrated against a third, double-precision implementation (`tests/truth_generator/dp_solver.py`, see "Why the tolerances are what they are").

## The held-out guides

Four guide files (`tests/data/run_a` … `run_d`, as `heldout_guides.json`), six guides, 13 test frequencies, 26 direction-spectra and 208 modes (142 of them hybrid).

| file | id | guide | loading | test frequencies (GHz) |
|---|---|---|---|---|
| run_a | H1 | WR-42 cross-section (10.668 × 4.318 mm) | magnetostatic-wave delay line: 100 µm single-crystal YIG film (LPE, 300 Oe, $ +y $, $ \Delta H = 1 $ Oe) on a 0.5 mm GGG substrate, 3 mm from a side wall | 2.75, 2.85 |
| run_b | H2 | WR-34 cross-section (8.636 × 4.318 mm) | the same film with a ground plane: 0.25 mm from a side wall, bias $ -y $, on 0.5 mm GGG | 2.75, 2.85 |
| run_c | H3 | WR-284 | two 5 mm zirconium tin titanate slabs ($ \varepsilon' = 37 $, $ \tan\delta = 2\times10^{-4} $), one against a side wall, one 20 mm from it | 3.0, 3.6 |
| run_c | H4 | WR-137 | low-pressure discharge cell: 5 mm argon at 1 Pa ($ n_e = 7.9\times10^{17} $ m$ ^{-3} $, $ f_p = 8 $ GHz, $ \nu = 2.4\times10^7 $ s$ ^{-1} $) between 4 mm fused-quartz plates, centred | 3.5, 3.6 |
| run_d | H5 | WR-62 cross-section | magnetised low-pressure discharge: 6 mm argon (same plasma, $ B_0 = 500 $ G along $ -y $, $ f_c = 1.4 $ GHz) against a side wall behind a 4 mm fused-quartz window | 3.5, 3.8, 4.0 |
| run_d | H6 | WR-112 | C-band ferrite phase shifter: 2.54 mm RO3006 against a side wall, then 2.5 mm YIG (Trans-Tech G-113, 600 Oe, $ -y $) | 7.5, 9.0 |

**Materials.**
- **YIG film** (H1, H2): single-crystal yttrium iron garnet grown by liquid-phase epitaxy on gadolinium gallium garnet, the standard medium of magnetostatic-wave devices [22, 23]. $ 4\pi M_s = 1780 $ G, $ \varepsilon' = 15 $, $ \tan\delta = 2\times10^{-4} $ as for the polycrystal; linewidth 1 Oe (LPE films: 0.3–1 Oe [22, 23]). The bias 300 Oe is an internal field.
- **GGG substrate**: $ \varepsilon' = 12.2 $, $ \tan\delta = 3\times10^{-4} $ [24].
- **Titanate ceramic** (H3): Zr–Sn–Ti oxide microwave ceramic, $ \varepsilon' = 37 $, $ \tan\delta = 2\times10^{-4} $ ($ Q \approx 5000 $ at 3–4 GHz) [25].
- **Low-pressure discharges** (H4, H5): argon at 1 Pa, $ n_g = 2.4\times10^{20} $ m$ ^{-3} $ at 300 K, $ \nu = n_g K_m = 2.4\times10^7 $ s$ ^{-1} $; $ n_e = 7.9\times10^{17} $ m$ ^{-3} $ ($ f_p = 8 $ GHz), the density of an inductively coupled or ECR source at this pressure [5]. $ f_c = 2.799 $ MHz/G $ \times $ 500 G $ = 1.4 $ GHz.

**What is in them.**
- **H1, H2: thin-film magnetostatic surface waves.** Between $ \sqrt{f_0(f_0 + f_m)} = 2.21 $ GHz and $ f_0 + f_m/2 = 3.33 $ GHz ($ f_0 = 0.84 $, $ f_m = 4.98 $ GHz) a film magnetised in its plane carries the Damon–Eshbach surface wave [21], with $ \beta d $ of order one: $ \beta $ is 2 800–10 300 rad/m (49–178 $ k_0 $) for $ d = 100 $ µm. The enclosure is far below cutoff (its empty TE$ _{10} $ cutoff is 14 GHz), so these are the only modes under the limit. With $ k_y = n\pi/b $ the wave also propagates obliquely to the bias, so **every family** $ n = 0, 1, 2, \dots $ carries one forward and one reverse mode, with $ |\beta| $ and $ \alpha $ growing with $ n $, until the attenuation passes 200 dB/m: at 2.75 GHz families $ n = 0 $–12 (26 modes, 108–194 dB/m), at 2.85 GHz $ n = 0 $–9. In H1 the two directions are localised on the two film faces and differ by a few rad/m; in H2 the side wall 0.25 mm from one face makes them strongly nonreciprocal at low $ n $ (LSE0 at 2.75 GHz: 2 605.8 rad/m and 118.9 dB/m forward, 940.7 rad/m and 46.5 dB/m reverse).
- **H3: many families.** With $ \varepsilon' = 37 $ the slabs carry LSE$ _n $ and LSM$ _n $ modes up to $ n = 3 $: 7 families and 11 modes per direction at 3.6 GHz.
- **H4: slow surface plasmons.** Just below the plasma–quartz surface resonance $ f_p/\sqrt{1 + 3.78} = 3.66 $ GHz, the 4 mm quartz plates are thick on the scale $ 1/\beta $, and LSM$ _1 $–LSM$ _4 $ are surface plasmons at up to 921 rad/m (12.2 $ k_0 $) at 3.6 GHz, with 96–176 dB/m.
- **H5: one-way surface magnetoplasmons.** The magnetised discharge supports a surface magnetoplasmon only in the forward direction (the reverse has no mode under the limit), in **every family** from $ n = 1 $ to 16 (3.5 GHz), 18 (3.8 GHz) and 14 (4.0 GHz). At 3.5 GHz they are backward waves ($ \beta $ from $ -119 $ to $ -2804 $ rad/m); at 3.8 GHz $ \beta $ stays between 64 and 78 rad/m while $ k_y $ reaches 7 160 rad/m; at 4.0 GHz $ \beta $ reaches 1 787 rad/m (21 $ k_0 $).
- **H6:** an ordinary nonreciprocal phase shifter (LSE0, HYB1, HYB2), so the held-out set is not made only of unusual guides.

**Why these break material-sized searches.** A search box set from the material constants, $ |\beta| \lesssim c\,k_0\sqrt{\max|\varepsilon|\max|\mu|} $, misses the slow waves: for H1 at 2.75 GHz $ 2k_0\sqrt{\max|\varepsilon|\max|\mu|} + 200 = 1\,600 $ rad/m against modes up to 10 271 rad/m; for H4 at 3.6 GHz $ 2k_0\sqrt{\max|\varepsilon|} + 200 = 500 $ rad/m against 921 rad/m. A family loop that stops at the first empty family beyond $ k_y^2 > k_0^2 \max|\varepsilon\mu| $, or after a fixed number of families, has to be right about where the ladders end. The scale that bounds the search is the quasi-static one: slow waves of a layer of thickness $ d $ have $ |\beta| d = O(1) $, and the instruction states $ d \ge 0.1 $ mm.

**How the test frequencies were chosen.** Materials and geometries come from devices; the frequencies were then chosen inside the surface-wave bands so that every mode satisfies the stated guarantees, in particular that no attenuation lies within 2 % of 200 dB/m. In the thin-film guides consecutive families differ by 7–13 dB/m near the end of a ladder, so roughly half of the frequencies tried had a mode in the excluded band (H1: 2.6, 2.65, 2.7, 2.8, 2.9, 2.95, 3.0 GHz were rejected; H5: 3.3, 3.4 and 3.6 GHz).

## The example guides (`/root/data/guides.json`)

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

Waveguide inner dimensions are the EIA standard ones [9]. Thicknesses are stock laminate, plate and ferrite-slab thicknesses or the remainder of the guide width. Test frequencies are round values, chosen so that each guide's four frequencies sample different regimes of its spectrum. No number is tuned to a root. These guides were the graded guides of v4.0; in v5.0 they are examples, and their certified v4.0 modes are no longer shipped.

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

- **The program is graded on guides it has not seen** (new in v5.0). Every check an agent can make on `/root/data/guides.json` (argument-principle audits, larger boxes, a second method) has to be built into the program, because nothing can be inspected or patched for the held-out guides. The held-out guides are in the stated model and ranges, but in regimes the example guides do not show.
- **The search region must be derived, not sized from the materials.**
  - The magnetostatic surface waves of a 100 µm YIG film (H1, H2) have $ \beta d = O(1) $, i.e. $ \beta $ of 2 800–10 300 rad/m, 49–178 $ k_0 $, with $ |\mu_e| \approx 6 $ and $ \varepsilon' = 15 $. A box $ |\beta| \le c\,k_0\sqrt{\max|\varepsilon|\max|\mu|} $ with any $ c $ of the usual size misses all of them; so does any search that scales with $ k_0 $. The relevant scale is $ 1/d $ of the thinnest layer, which the instruction bounds ($ d \ge 0.1 $ mm).
  - The surface plasmons of H4 (12 $ k_0 $) lie past $ k_0\sqrt{\max|\varepsilon|} $ because the plasma–quartz pair is near its surface resonance; the bound comes from $ k_0^2\varepsilon_p\varepsilon_d/(\varepsilon_p + \varepsilon_d) $, not from either permittivity.
- **The family index has no small bound.** In H1 and H2 every family $ n = 0, \dots, 12 $ carries a mode, and in H5 every family $ n = 1, \dots, 18 $; their $ k_y $ (up to 7 160 rad/m in H5) exceeds $ k_0\sqrt{\max|\varepsilon\mu|} $ by a factor of 20 or more. A loop that stops on a fixed $ n $, or as soon as $ k_y $ passes the material wavenumber, ends early. The ladders end only where $ \alpha $ crosses 200 dB/m.
- **One-way and backward slow waves together.** The H5 magnetoplasmons exist only forward, in 14–18 families, and are backward waves ($ \beta < 0 $) at 3.5 GHz, forward waves at 4.0 GHz, and nearly standing ($ \beta \approx 70 $ rad/m with $ k_y $ up to 7 160 rad/m) at 3.8 GHz.
- **Stiffness at large $ |\gamma| $.** For the films $ |\gamma| a $ reaches 110 (H1), so $ e^{\pm\gamma x} $ spans $ e^{220} $ across the guide: a transfer-matrix evaluation must renormalise and keep its two state columns independent, a collocation must resolve boundary layers of width $ 1/\beta \approx 0.1 $ mm, and polishing must not judge convergence by $ |D| $.
- **A correct $ \gamma $ does not give a correct power split.** The fields of a slow wave are found by carrying the tangential state from one wall across the guide. An error $ \delta $ in $ \gamma $ (or in the null vector of the wall condition) excites the growing solution, which grows like $ \delta\, e^{|\gamma| x} $ while the mode itself decays away from the film; its power grows like $ \delta^2 e^{2|\gamma| x} $. With $ \gamma $ rounded to double precision ($ \delta \approx 10^{-16} $) and $ |\gamma| a \approx 110 $, this artefact is $ 10^{63} $ times the far-layer power of the mode. Both programs written for v5.0 first made exactly this error: for the H2 mode HYB9 at 2.85 GHz they reported $ p = (0, 0, 0, 1) $ (all power in the far air gap) instead of $ (1.94, -0.98, 0.04, 10^{-6}) $, with $ \gamma $ correct to $ 10^{-16} $, and the v4.0 reference program does the same on 20 modes of H1, 26 of H2 and 27 of H5. The power split needs $ \gamma $ to about $ 0.87 |\gamma| a $ more digits than the required tolerance, or an integration from both walls.
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
  - In the plasma guides the fractions are large and of both signs (backward power in the plasma, forward power in the dielectrics), e.g. G08 at 4.4 GHz, LSM1: $ p = (-0.177, -0.733, 2.821, -0.733, -0.177) $ (air, quartz, plasma, quartz, air); 73 of the 411 modes have a fraction below $ -10^{-3} $.
- **Gyrotropic guides are nonreciprocal** (G13–G18, 48 of the 144 direction-spectra).
  - With the Polder tensor, the y-independent TE problem has $ E'' + (k_0^2\varepsilon\mu_e + \gamma^2)E = 0 $ in each layer, $ \mu_e = (\mu^2 - \kappa^2)/\mu $. The continuous field is $ H_z \propto (E' - j(\kappa/\mu)\gamma E)/\mu_e $.
  - The term linear in $ \gamma $ makes the dispersion relation depend on $ \gamma $ itself, not on $ \gamma^2 $: the forward and reverse modes differ.
  - G13 at 8 GHz has its slab mode at $ \beta = 452.5 $ rad/m forward (8.0 dB/m) against 344.1 rad/m reverse (19.0 dB/m). G16 at 5 GHz has 228.0 against 157.0 rad/m. G14 (twin slabs, opposite bias) differs by 1–16 rad/m.
  - The resonance isolator G15 at 11 GHz has a forward mode at 166 dB/m and no reverse mode under the limit. At 11.5 GHz the two directions have different mode sets.
  - Assuming reciprocity, or getting the sign of $ \kappa $ or of $ g $ relative to the propagation direction wrong (it swaps the two directions), fails these guides. So does a search in $ q = -\gamma^2 $, which cannot represent the asymmetry.
  - The hybrid families are more nonreciprocal than LSE$ _0 $: G13 at 12 GHz has HYB1 at 474.7 and 712.7 rad/m forward against 456.1 and 642.7 rad/m reverse; G17 at 7.5 GHz has HYB1 at 192.5 rad/m (130.6 dB/m) forward against 138.6 rad/m (121.6 dB/m) reverse, with power splits of opposite character ($ +1.85 $ in the RO3010 and $ -1.46 $ in the plasma forward, $ -0.16 $ and $ -1.04 $ reverse).
  - **One-way surface magnetoplasmons** (G18). At 4.5 GHz, below the upper-hybrid resonance of the discharge, HYB1 (432.2 rad/m, 156.6 dB/m) and HYB2 (486.5 rad/m, 188.7 dB/m) exist only forward: no reverse mode of these families is under the limit. Both carry backward power in the plasma ($ p = -1.29 $ and $ -0.88 $). This is the lossy, bounded version of the one-way surface magnetoplasmon [18].
- **The modes are complex.** With lossy layers the dispersion functions have no zeros on the real axis, so scanning for sign changes finds nothing. The search is over a region of the complex plane: the parabola $ \operatorname{Re} q \ge k_n^2 - A^2 + (\operatorname{Im} q)^2 / 4A^2 $, or for ferrite guides the strip $ |\operatorname{Re}\gamma| \le A $.
- **Backward waves** (G07, G08, G09, G11). Below the surface-plasma resonance of a plasma column the LSM surface wave carries power and phase in opposite directions [7, 8]. With $ \alpha > 0 $ it has $ \beta < 0 $: G08 at 4.4 GHz has LSM1 at $ -539 $ rad/m (132 dB/m) and LSM2 at $ -415 $ rad/m (172 dB/m).
- **Evanescent modes inside the limit** (G03, G06, G08, G09, G10). Just below a cutoff a mode decays with $ \beta \approx 0 $ but can still be under 200 dB/m: G06 at 4.15 GHz has LSE0 at 168 dB/m.
- **Slow surface waves** (G12 at 5 and 5.15 GHz, $ \operatorname{Re} q \approx 41 $–$ 47\,k_0^2 $). A search bounded by the dielectric permittivities misses them.
- **Loss is not a perturbation here.** A lossless analysis followed by first-order loss perturbation fails every guide, and near cutoffs it changes the mode set.
- **Bookkeeping.** There are 411 reported modes (64 of them hybrid) in 144 direction-spectra, with up to 4 families and 8 modes per spectrum, and empty spectra. LSM$ _0 $ and HYB$ _0 $ do not exist, but LSE$ _0 $ does; each mode carries one power fraction per layer (1843 numbers in all).

**The answers are those of the lossy guide, not of an ideal one.** Halving all losses (dielectric, plasma) changes 4 of the 48 isotropic spectra, and doubling them changes 7. Only the most attenuated backward and surface waves cross the limit, as their attenuation scales with the loss. The backward waves of G08 at 4.8 GHz (61 dB/m) and G11 (37–45 dB/m) persist.

A single wrong item in any held-out guide fails the task. The bullets after the first six describe the v4.0 difficulties, which the held-out guides H3, H4 and H6 keep and which the program must still handle on unseen guides.

## Reference solution

`solution/solve.sh` copies `solution/modes.py` to `/root/results/modes.py`. It is one file (numpy, scipy, mpmath) in two parts: the isotropic guides in the $ q = k_n^2 - \gamma^2 $ plane, and the gyrotropic guides in the $ \gamma $ plane.

**Isotropic guides.** Per test frequency and kind (LSE, LSM):
1. **Counting.** The zeros of $ D_{\rm LSE}(q) $ and $ D_{\rm LSM}(q) $ (transfer matrices) are counted with the argument principle [10] on rectangles covering the region $ \alpha \le 1.02 A $ of the $ q $ plane ($ A = 200 $ dB/m $ = 10\ln 10 $ Np/m). $ \arg D $ is followed with steps limited to 0.4 rad; the state is renormalised by positive real factors, which leave $ \arg D $ unchanged.
2. **Locating.** Each single zero is located by the contour moment and Newton's method on $ D'/D $, then polished in mpmath at 50 digits.
3. **Bounds.** The Rayleigh bound for LSE. For LSM, beyond every interface surface-wave estimate $ k_0^2\varepsilon_i\varepsilon_j/(\varepsilon_i + \varepsilon_j) $ and four times past the last zero; near a surface resonance ($ \varepsilon_i + \varepsilon_j \to 0 $, H4) the estimate grows and the search follows it.
4. **Reporting.** Each zero gives $ \gamma = \sqrt{k_n^2 - q} $ ($ \operatorname{Re}\gamma > 0 $) for every family. These guides are reciprocal, so the reverse list equals the forward one.

**Gyrotropic guides.** For $ n = 0, 1, 2, \dots $:
1. **Dispersion function.** $ D_n(\gamma) $: for $ n = 0 $ the TE state $ (e_y, h_z) $, for $ n \ge 1 $ the $ 2 \times 2 $ determinant of the $ e $-rows of $ T(a)[0; I] $. Isotropic layers use $ \cos(Kd) I + \sin(Kd)/K\,M $ scaled by $ e^{-|\operatorname{Re}\kappa| d} $; gyrotropic layers the eigendecomposition of $ Md $ (`scipy.linalg.expm` with substeps if the eigenvectors are ill-conditioned), scaled by $ e^{-\max\operatorname{Re}\lambda} $. After every layer the two state columns are re-orthonormalised by Gram–Schmidt (a unimodular column operation and positive scalings), so $ \arg D $ is unchanged and the columns stay independent for $ |\gamma| \gg k_0 $, where $ |\gamma| a $ reaches several hundred.
2. **Counting.** The argument principle on rectangles covering $ 0 < \operatorname{Re}\gamma \le 1.02A $ (forward) and $ -1.02A \le \operatorname{Re}\gamma < 0 $ (reverse): first $ |\operatorname{Im}\gamma| \le B_0 = 4 k_0 \sqrt{m} + 2k_y + 200 $ ($ m = \max|\varepsilon|\max(|\mu|, |\mu_e|, 1) $) in 24 strips per sign, then **shells** $ B \le |\operatorname{Im}\gamma| \le 2B $ in 16 strips each, doubling $ B $ until a shell is empty **and** $ B $ is past the quasi-static scale $ 20/d_{\min} $ of the thinnest ferrite or plasma layer. Slow quasi-static waves have $ |\beta| d = O(1) $; beyond $ 20/d_{\min} $ every interface is decoupled from the next, and the last, empty shell excludes a single-interface resonance.
3. **Locating.** Rectangles are bisected (off-centre) until each holds one zero and is $ 10^{-4} $ (relative) small, then the zero is polished in mpmath at $ 45 + |\gamma| a/2.3 $ digits (secant iteration judged by its step, since $ |D| $ grows like $ e^{|\gamma| a} $); the polished value is kept for the power split.
4. **Stopping in $ n $.** After the first $ n \ge 1 $ without a zero once $ k_y^2 > k_0^2 m + A^2 $. The thin-film and magnetoplasmon ladders are contiguous in $ n $ and end where $ \alpha $ passes $ A $; the author-side scans below search 40–60 families further and to $ |\beta| = 4\times10^5 $ rad/m.
5. **Reporting.** A reverse zero $ \gamma $ is reported as $ -\gamma $, with the power split of its own field.

**Power split (every guide).** At the polished $ \gamma $ (its high-precision value, not the double reported) the null vector of the $ 2 \times 2 $ wall condition gives $ \psi(0) $; $ \psi $ is carried through each layer by $ \exp(M t) $ in mpmath at $ 40 + G/2.3 $ digits, $ G = \sum_i \max|\operatorname{Re}\lambda_i| d_i $ the growth of the field across the guide; $ e_x $ and $ h_x $ follow algebraically, and $ \operatorname{Re}(e_x h_y^* - e_y h_x^*) $ is integrated by 24-point Gauss–Legendre quadrature on pieces over which the field grows by at most $ e^{0.5} $. With $ \gamma $ rounded to double precision instead, the far-layer fractions of the slow waves are dominated by the excited growing solution (Difficulty).

**Run time.** Through the real verifier, with the four held-out files run at once on 4 cores, all four runs finish within 186 s (alone: 7 s for run_c, 2.5–3.3 min for the others). The 18 example guides take 7.5 min and reproduce the certified v4.0 values (worst $ |\Delta\gamma|/|\gamma| $ $ 6\times10^{-14} $, power $ 10^{-15} $). The output is deterministic.

## Verification

**How the program is run** (`tests/test_modes.py`). The grading container holds Python 3.11 with numpy 2.2.6, scipy 1.15.3, sympy 1.13.3 and mpmath 1.3.0 (the agent environment's versions), the four held-out guide files and their certified modes (`/tests/data/run_*/heldout_guides.json`, `heldout_truth.json`; `test.sh` makes `/tests` mode 700), and the submitted program.
- The verifier opens `/root/results/modes.py` without following symlinks and without blocking (a FIFO is rejected at once), requires a regular file of at most 50 000 000 bytes (stated in the instruction) and copies its content.
- For each guide file, in a fixed order (sorted by name, so that nothing in the verdict depends on chance), it stages a fresh run directory: a read-only input directory, with a name derived from a hash of its private temporary directory, holding only `guides.json`, and a working directory (mode 1777) holding only the copy of the program. It first checks, as the sandbox user, that the working directory and `/tmp` are writable, the input directory is not, and `/tests` is unreadable.
- It starts all four runs at once as uid/gid 65534 in their own sessions, with the environment stated in the instruction and stdout/stderr going to unlinked files. Each run must exit with status 0 within 3600 s; a run past the limit has its process group killed. Afterwards every process of uid 65534 is killed (a helper that has become uid 65534 calls `kill(-1, SIGKILL)` until `/proc` shows none).
- The output is opened the same way as the program and must be a regular file of at most 16 MiB. It is parsed as data, never executed: NaN, Infinity and repeated keys are rejected, and the whole output contract is checked before any value.

**The tests**, each parametrised over the four runs (20 tests; reward 1 only if all pass):
- **`test_run_completes`** (requirement 1): exit status 0 within the limit, output a regular file.
- **`test_output_contract`** (requirement 2): guide ids, the single key `modes`, one entry per test frequency with exactly `forward` and `reverse`, family names that fully match `LSE<n>` ($ n \ge 0 $), `LSM<n>` or `HYB<n>` ($ n \ge 1 $), at most 400 families and 1000 modes, modes `[beta, alpha, power]` of finite numbers with $ \alpha > 0 $ and one fraction per layer.
- **`test_mode_families`** (requirement 3): for every guide, test frequency and direction the families and their mode counts equal the reference.
- **`test_propagation_constants`** (requirement 4): the reported $ \gamma $ pair one-to-one with the reference within $ |\gamma - \gamma^*| \le 10^{-7}|\gamma^*| $.
- **`test_power_split`** (requirement 5): for each pair, every fraction satisfies $ |p_i - p_i^*| \le 10^{-6} \sum_k |p_k^*| $.

**Why the tolerances are what they are** (unchanged from v4.0, rechecked on the held-out guides, and calibrated in v5.1 against a double-precision implementation).
- **Correct solvers pass easily.** The two correct programs reproduce the certified values of the held-out guides to $ 10^{-14} $ in $ \gamma $ and $ 6\times10^{-16} $ in the fractions. Every $ \gamma $ scaled by $ 1 + 5\times10^{-8} $ passes; one off by $ 3 \times 10^{-7} $ fails.
- **Pairing is unambiguous.** Modes of a family and direction are at least $ 10^{-4}|\gamma| $ apart (certified, stated in the instruction).
- **Membership never depends on rounding.** No mode is within 2 % of the limit (certified, stated in the instruction).
- **The power tolerance.** The certified fractions are enclosed to better than $ 10^{-9} $. Every fraction scaled by $ 1 + 3\times10^{-7} $ passes; one fraction off by $ 3\times10^{-6}\sum|p^*| $ fails. The instruction guarantees $ \sum|P_i| \le 30|P| $ (largest certified value on the held-out guides 15.6, H4).
- **The wall-clock limit.** The reference program needs at most 186 s for a guide file and the second correct program 734 s, against 3600 s; the four runs share the 4 CPUs of the grading container.
- **The tolerances need accurate numerics, not extended precision** (v5.1). `tests/truth_generator/dp_solver.py` is a third implementation in which every step is carried out in double precision (numpy and scipy). It takes the candidate zeros from the collocation of `alt_solver.py` and replaces its two mpmath stages:
  - **polishing:** secant iteration on a determinant built from both walls. The two wall solutions are carried towards a matching interface in steps with $ |\kappa| h \le 2 $ and re-orthonormalised after every step, so the span of each set stays accurate in the direction in which it grows. The matching interface is the one at which the two sets are closest to dependent (for a surface wave, a face of the film it is bound to);
  - **power split:** the null vector at the matching interface fixes the mode; its coefficients are carried back towards each wall through the stored triangular factors, and the flux is integrated step by step. The field is never carried across the film in the direction in which it decays, so the $ \delta\, e^{|\gamma| x} $ error of a one-sided shot (see "Difficulty") does not arise.

  With the four held-out files run at once on 4 cores it passes all 20 tests (runs take 36 s, 895 s, 970 s and 1117 s). Its largest errors against the certified values, over all 208 modes:

  | run | guides | $ \max |\gamma - \gamma^*| / |\gamma^*| $ | $ \max |p_i - p_i^*| / \sum_k |p_k^*| $ |
  |---|---|---|---|
  | run_a | H1 | $ 2.5 \times 10^{-12} $ | $ 1.4 \times 10^{-9} $ |
  | run_b | H2 | $ 1.4 \times 10^{-12} $ | $ 4.7 \times 10^{-9} $ |
  | run_c | H3, H4 | $ 7.6 \times 10^{-13} $ | $ 5.9 \times 10^{-13} $ |
  | run_d | H5, H6 | $ 5.7 \times 10^{-13} $ | $ 1.7 \times 10^{-11} $ |
  | tolerance | | $ 10^{-7} $ | $ 10^{-6} $ |

  So a double-precision program clears the tolerances by a factor of 40 000 in $ \gamma $ and 200 in the power split.
- **A collocation eigenproblem alone is not accurate enough, and refining its grid does not help.** The eigenvalues and null vectors of the multi-domain Chebyshev collocation of `alt_solver.py`, without its polishing, for the forward HYB3 modes of H1 at 2.75 GHz (the per-layer cap of 200 nodes lifted; scale $ s $ as in `alt_solver.py`, whose default is 2.5):

  | $ s $ | nodes | $ \max |\gamma - \gamma^*| / |\gamma^*| $ | $ \max |p_i - p_i^*| / \sum_k |p_k^*| $ |
  |---|---|---|---|
  | 1.0 | 124 | $ 3.9 \times 10^{-9} $ | $ 1.0 \times 10^{-5} $ |
  | 1.5 | 155 | $ 7.0 \times 10^{-10} $ | $ 1.7 \times 10^{-6} $ |
  | 2.5 | 216 | $ 6.4 \times 10^{-9} $ | $ 1.6 \times 10^{-5} $ |
  | 5.0 | 369 | $ 1.5 \times 10^{-8} $ | $ 1.6 \times 10^{-5} $ |
  | 7.5 | 523 | $ 6.9 \times 10^{-8} $ | $ 9.7 \times 10^{-5} $ |

  Beyond about 150 nodes the errors grow: the rounding error of the Chebyshev differentiation matrices grows with the number of nodes, and the slow waves have boundary layers of width $ 1/\beta \approx 0.1 $ mm in a 10.7 mm guide. The collocation therefore has an accuracy floor of about $ 10^{-9} $ to $ 10^{-7} $ in $ \gamma $ and $ 10^{-6} $ to $ 10^{-4} $ in the power split, at or above the tolerances. It locates the modes (both `alt_solver.py` and `dp_solver.py` use it for that), but its values have to be polished on the exact dispersion function and the fields recomputed, which `dp_solver.py` does in double precision.
- **Calibration summary.**

  | quantity | worst valid result | tolerance | nearest defective result | margin |
  |---|---|---|---|---|
  | $ \gamma $ (relative) | $ 10^{-14} $ (reference, `alt_solver.py`); $ 2.5 \times 10^{-12} $ (`dp_solver.py`) | $ 10^{-7} $ | $ 7 \times 10^{-8} $ to $ 2.5 \times 10^{-7} $ (unpolished collocation, above; the v5.0 Opus run); $ 3 \times 10^{-7} $ (ablation) | $ 4 \times 10^{4} $ above the double-precision program |
  | power split (of $ \sum|p^*| $) | $ 6 \times 10^{-16} $ (reference, `alt_solver.py`); $ 4.7 \times 10^{-9} $ (`dp_solver.py`) | $ 10^{-6} $ | $ 1.7 \times 10^{-6} $ to $ 10^{-4} $ (collocation null vectors); $ 3 \times 10^{-6} $ (ablation); $ O(1) $ (one-sided shot from a double-precision $ \gamma $) | 200 above the double-precision program |

**Where the held-out guides come from** (`tests/truth_generator/gen_heldout.py`). The script writes the four guide files exactly as shipped; its docstring repeats the sources above. The example file `environment/data/guides.json` is written by `gen_guides.py` (v4.0, unchanged).

**Where the reference values come from** (`make_truth.py`, `certify2.py`, `wgacb.py`, `hybacb.py`, `power.py`). The generator is separate from `solution/` and never uses the argument principle. It certifies with ball arithmetic (Arb, through python-flint, 256 bits) [11].
- **Isotropic guides** (H3, H4): branch-and-bound over the region $ \alpha \le 1.02A $ of the $ q $ plane, with mean-value exclusion, complex Krawczyk acceptance [12] and refinement to $ 10^{-55} $ boxes; the Rayleigh bound for LSE and, for LSM, a proof that $ D_{\rm LSM} \ne 0 $ for $ \operatorname{Re} q \ge Q $ (impedance recursion $ t \mapsto (t + T)/(1 + tT) $, $ T = \tanh\kappa d $, on narrow boxes and an explicit bound for the unbounded rest). This is exhaustive: every mode of H3 and H4 is found and certified without reference to any solver.
- **Gyrotropic guides, $ n = 0 $** (H1, H2, H6): the same branch-and-bound in the $ \gamma $ plane over $ |\operatorname{Re}\gamma| \le 1.02A $, $ |\operatorname{Im}\gamma| \le B $, with $ B $ proved by the recursion for $ \sigma = \kappa E/E' $ across layers and interfaces (including the nonreciprocal interface term) to the unbounded rest. For the films $ B $ is 179–199 $ k_0 $ (10 700–11 400 rad/m), past their LSE0 surface waves at 2 600–3 300 rad/m. Exhaustive as well.
- **Hybrid families, $ n \ge 1 $.** From candidates: every mode reported by the reference program and by `alt_solver.py` (which must agree on every family, every count, and every $ \gamma $ to $ 10^{-9} $) is certified individually by a Krawczyk box (existence and uniqueness of a simple zero), refined to $ 10^{-55} $, and checked against the guarantees.
- **Completeness of the hybrid families** is checked author-side by `exhaust.py`, which is independent of both programs' stopping rules: for every gyrotropic held-out guide and frequency it counts zeros with the argument principle over $ |\operatorname{Re}\gamma| \le 1.02A $ and $ |\operatorname{Im}\gamma| \le 4\times10^5 $ rad/m (twice the quasi-static scale $ 20/d $ of the film, 40 times the largest mode) in families $ n = 0 $ to 60 (H1, H2) or 80 (H5, H6), and compares them with the certified list. It finds exactly the certified modes, 158 of 158, and no other zero with $ \alpha \le A $; the closest approach to the limit is 2.4 % (H5 at 3.8 GHz, HYB18, 195.1 dB/m).
- **Power split** (`power.py`): the null vector of the certified wall condition, then $ P_i = \operatorname{Re}\int (e_x\, \overline{h_y(\bar t)} - e_y\, \overline{h_x(\bar t)})\,dt $ by `acb.integral` to relative accuracy $ 2^{-90} $, on pieces of each layer over which the field grows by at most $ e $, each started from the state carried to its start.
- **Precision for the slow waves.** With $ \kappa a $ up to 116 ($ \kappa = \sqrt{|\gamma|^2 + k_y^2} $), $ D $ loses about $ 2.9\kappa a $ bits to cancellation and the enclosures of the fields grow like $ e^{2\kappa x} $. For the gyrotropic guides the generator therefore converts the inputs to balls at 1024 bits and works at 1024 bits; a candidate is first refined by Newton's method and certified on a box of $ 10^{-30} $ to $ 10^{-6} $ relative size (at $ 10^{-11} $ the enclosure of $ D' $ is already 0.2 of $ |D'| $ for the film modes), then refined to $ 10^{-(55 + 0.9\kappa a)} $; the power split runs at $ 256 + 4\kappa a $ bits. Without these steps the certificate fails for the film modes from $ n = 2 $ on, and the power of the highest families is not resolved to $ 10^{-9} $.
- **Guarantees.** The generator stops unless, for every mode, $ \operatorname{Re}\gamma $ is certainly nonzero and $ |\alpha - A| > 0.02A $, modes of a family and direction are at least $ 10^{-4}|\gamma| $ apart, $ \sum|P_i|/|P| < 30 $, and every fraction is enclosed to $ 10^{-9} $. All 26 held-out direction-spectra pass, in about 5 200 s of certified computation, most of it the power split of the 48 H5 magnetoplasmons; run_d was certified in four parallel parts (H5 at each frequency, and H6) whose results were merged, which is the same computation (`tests/data/run_*/truth_diagnostics.json`).

**Score of an independent correct submission.** `tests/truth_generator/alt_solver.py` shares no code or method with the reference program or the generator: multi-domain Chebyshev collocation [13] (isotropic: $ v'' + k_0^2\varepsilon v = qv $; $ n = 0 $: the TE problem as a quadratic eigenproblem; $ n \ge 1 $: all six curl equations as a linear pencil $ A_0 x + \gamma A_1 x = 0 $), eigenvalues kept if two grids agree, polished by secant iteration on determinants of matrix exponentials in mpmath, and the power split from the collocation fields by Clenshaw–Curtis quadrature. For v5.0 its grids are sized for a target $ |\beta| $ that starts at the quasi-static scale $ 3/d_{\min} $ and doubles while a zero lies beyond half of it (each layer gets $ 2s\sqrt{\beta_t d} + 16 $ points, which resolves the $ e^{-\beta x} $ boundary layers on Chebyshev nodes), and its polishing is judged by the step. Run as `/root/results/modes.py` through the real verifier it passes all 20 tests (all four runs at once finish within 734 s).

**Ablations** (`tests/truth_generator/ablate.py` for edited outputs; whole programs through the real verifier).

*Whole programs*, each run as `/root/results/modes.py` through the real verifier (20 tests):

| program | result | tests failed |
|---|---|---|
| reference program (`solution/modes.py`) | reward 1 | none (within 186 s) |
| second correct program (`alt_solver.py`) | reward 1 | none (within 734 s) |
| no program | reward 0 | all 20 |
| the v4.0 reference program (material-sized search, power split from double-precision $ \gamma $) | reward 0 | 7: families, $ \gamma $ and power of run_a and run_b (the LSE0 surface waves of H1 and H2 at 2.85 GHz are missed); power of run_d (27 H5 modes) |
| the reference with the search box of the passing v4.0 Opus run, $ |\beta| \le 2k_0\sqrt{\max|\varepsilon|\max|\mu|} + 200 $, in every family and for the isotropic $ q $ search | reward 0 | 12: every run (no film mode, the slow plasmons of H4, the H5 families from $ n = 4 $ on at 3.5 and 4.0 GHz) |
| the reference with the widest box of the passing v4.0 GPT-6 run, $ |\beta| \le 16k_0 $ | reward 0 | 9: run_a, run_b, run_d |

*Edited correct outputs* (`ablate.py`, judged by the verifier's own functions): 23 of 23 behave as expected.
- **These pass:** the correct outputs; every $ \gamma $ scaled by $ 1 + 5\times10^{-8} $; every fraction scaled by $ 1 + 3\times10^{-7} $; families and modes listed in another order.
- **Search region and families fail:** modes beyond a material box $ 10k_0\sqrt{\max\varepsilon'} + 200 $ dropped; modes beyond $ 50k_0 $ dropped; the last family of each H1 spectrum dropped; families beyond $ n = 8 $ dropped; the slow plasmons of H4 dropped; the LSM families of H3 dropped; the most attenuated mode of each H1 spectrum dropped.
- **Nonreciprocity fails:** H5 treated as reciprocal; H2 treated as reciprocal; the bias sign ignored in H2 and H5; $ \beta > 0 $ forced (the backward waves of H5 at 3.5 GHz).
- **Values fail:** one $ \gamma $ off by $ 3\times10^{-7} $; unsigned fractions $ |P_i|/\sum|P_k| $; reverse modes of H2 given the forward split; one fraction off by $ 3\times10^{-6}\sum|p^*| $.
- **Format fails:** a guide missing; NaN; an HYB0 family; a power list one number short.

The first versions of both v5.0 programs, which computed the power split from $ \gamma $ rounded to double precision, reproduced every $ \gamma $ to $ 10^{-15} $ but failed `test_power_split` on runs a, b and d (37, 35 and 27 modes off by up to 1.9 $ \sum|p^*| $).

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
- **v4.0 calibration (run #7064700): both models passed** (reward 1, 4 of 4 tests).
  - **GPT-6 (max)**, 20 min. It wrote a 4 × 4 transfer-matrix solver for every family, located zeros from seeds on a grid of $ |\gamma| \le 5k_0 $ ($ 7k_0 $ for gyrotropic guides), and then audited the counts with the argument principle on boxes $ |\operatorname{Re}(\gamma/k_0)| \le 10 $, later 16. The audit found five LSM1 pairs that its search had missed (G08 at 4.4 GHz, G09 at 2.75 and 2.9 GHz, G12 at 5 and 5.15 GHz); it found them again by local searches and added them to the file by hand.
  - **Opus 5.5 (low)**, 36 min. It counted zeros with the argument principle over $ |\operatorname{Im}\gamma| \le 2k_0\sqrt{m} + 200 $ ($ m = \max|\varepsilon|\max|\mu| $), stopped in $ n $ after four empty families, and reran with a finer contour to check the counts. It did not check the hybrid modes or the power split independently.
  - What the passing runs show: the graded guides were the agent's own input, so a solver whose search box comes from the material constants could be audited and repaired guide by guide. v5.0 therefore grades a program on held-out guides, whose slow waves lie far outside such boxes (H1 and H2 up to 178 $ k_0 $, H4 12 $ k_0 $, H5 in up to 18 consecutive families). Both v4.0 search strategies, reimplemented on the reference program, fail the held-out guides (ablations above).
- **v5.0 calibration: both models failed** (reward 0).
  - **GPT-6 (max)**, 7 min of its 5 h. Its root search bounds $ \beta $ by a material wavenumber and stops the family index on a similar bound. It reported no family at all for the forward modes of H1 and H2 at 2.75 GHz (13 expected), 2 of the 16 forward families of H5 at 3.5 GHz, and only LSM3 of LSM1–LSM3 for H4 at 3.5 GHz: the slow waves the held-out guides were chosen for. Its power split, propagated in double precision across the guide, was also off (H3 at 3.6 GHz, forward LSM3: $ 4.8 \times 10^{-6} $ instead of $ 6.3 \times 10^{-11} $ in layer 3).
  - **Opus 5.5 (low)**, 2 h 7 min; runs of 621–693 s. It solved a Chebyshev collocation of the Maxwell equations, which needs no search box, and found every family of every spectrum. It reported the collocation eigenvalues and the power split of the collocation eigenvectors without polishing either and without a convergence check. The first failing mode of each run: $ \gamma $ off by $ 2.5 \times $, $ 1.2 \times $ and $ 1.8 \times $ the tolerance (H1 HYB3, H2 HYB7, H5 HYB12), power fractions off by $ 2.6 \times 10^{-6} $, $ 6.8 \times 10^{-6} $ and $ 1.2 \times 10^{-6} $ of $ \sum|p^*| $ (H1 HYB2, H2 LSE0, H5 HYB7 at 3.5 GHz); run_c passed. The tests stop at the first failing mode of a run, so these are not necessarily its largest errors.
  - What the Opus run shows: these errors are the accuracy floor of a collocation eigenproblem on these guides, which refining the grid does not lower (table under "Why the tolerances are what they are"). Polishing the eigenvalues on the exact dispersion function and recomputing the fields from both walls brings the same method to $ 2.5 \times 10^{-12} $ in $ \gamma $ and $ 4.7 \times 10^{-9} $ in the power split in double precision (`dp_solver.py`). The failure is a numerical-method error (unpolished values of an ill-conditioned discretisation, not checked for convergence), not a tolerance at the limit of double precision.

**Verifier runs** (v5.1.2, with the fixed launch order): the verifier image built from `tests/Dockerfile` (Python 3.11 with numpy 2.2.6, scipy 1.15.3, sympy 1.13.3, mpmath 1.3.0), run with `tests/test.sh` without network, on 4 CPUs and with 8 GB:
- with the reference program as `/root/results/modes.py`: reward 1, 20/20 tests passed, 136 s;
- with no program: reward 0, 20/20 tests failed.

(Harbor itself could not build the images in the authoring sandbox, whose network proxy is not reachable from inside a Docker build; the image above was built from a copy of `tests/` with the proxy's certificate added.)

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
21. R. W. Damon and J. R. Eshbach, "Magnetostatic modes of a ferromagnet slab," *J. Phys. Chem. Solids* 19 (1961) 308–320.
22. D. D. Stancil and A. Prabhakar, *Spin Waves: Theory and Applications*, Springer, 2009: magnetostatic surface waves of in-plane magnetised films, oblique propagation, YIG films on GGG, linewidths.
23. W. S. Ishak, "Magnetostatic wave technology: a review," *Proc. IEEE* 76 (1988) 171–187: LPE YIG films, ground planes, delay lines and filters.
24. Gadolinium gallium garnet substrate data sheets (supplier data): $ \varepsilon_r \approx 12.2 $, low microwave loss.
25. Zirconium tin titanate dielectric-resonator ceramics (e.g. Trans-Tech D-8300 series): $ \varepsilon_r \approx 36\text{–}37 $, $ Q \cdot f \ge 2\times10^{4} $ GHz.
26. J. J. Brion, R. F. Wallis, A. Hartstein and E. Burstein, "Theory of surface magnetoplasmons in semiconductors," *Phys. Rev. Lett.* 28 (1972) 1455–1458.
