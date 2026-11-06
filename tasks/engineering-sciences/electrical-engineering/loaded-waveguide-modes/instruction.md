Find the guided modes of each lossy, layer-loaded rectangular waveguide in `/root/data/guides.json` and write them to `/root/results/modes.json`. For each guide, each of its test frequencies and each direction of propagation, report every mode whose attenuation is at most 200 dB/m, with its phase constant, its attenuation constant and how its power flow divides among the layers.

## Input

`/root/data/guides.json` is a UTF-8 JSON object with the single key `"guides"`, a list of guide records. Each record has these keys:

- `id`: string, unique.
- `description`: what the guide is (informational).
- `a`, `b`: the inner width and height of the guide in metres.
- `layers`: a list of layers that fill the guide from $ x = 0 $ to $ x = a $ in order. Each layer has a `thickness` in metres (the thicknesses add up to `a` exactly) and a `material`. `name` says which material or discharge the layer represents; it is informational. A material is one of:
  - `{"type": "dielectric", "name": ..., "eps": ..., "loss_tangent": ...}`: relative permittivity $ \varepsilon = \varepsilon' (1 - j \tan\delta) $ with $ \varepsilon' = $ `eps` and $ \tan\delta = $ `loss_tangent`, the same at every frequency; relative permeability 1.
  - `{"type": "plasma", "name": ..., "f_p": ..., "collision_frequency": ...}`: a collisional (Drude) electron plasma, $ \varepsilon(f) = 1 - \dfrac{f_p^2}{f\,(f - j\nu / 2\pi)} $ with plasma frequency $ f_p = $ `f_p` (Hz) and electron collision frequency $ \nu = $ `collision_frequency` (s$ ^{-1} $); relative permeability 1.
  - `{"type": "ferrite", "name": ..., "eps": ..., "loss_tangent": ..., "saturation_magnetization_gauss": ..., "linewidth_oe": ..., "bias_field_oe": ..., "bias_direction": ...}`: a saturated ferrite magnetised along $ +y $ or $ -y $ (`bias_direction` is `"+y"` or `"-y"`). Its permittivity is $ \varepsilon'(1 - j\tan\delta) $ as for a dielectric. Its relative permeability is the Polder tensor, in the coordinates $ (x, y, z) $:

    $$
    [\mu] = \begin{pmatrix} \mu & 0 & -j\kappa \\ 0 & 1 & 0 \\ j\kappa & 0 & \mu \end{pmatrix}, \qquad \mu = 1 + \frac{f_0 f_m}{f_0^2 - f^2}, \qquad \kappa = s\,\frac{f f_m}{f_0^2 - f^2},
    $$

    with $ s = +1 $ for `"+y"` and $ s = -1 $ for `"-y"`, $ f_0 = \gamma_g H_0 + j \gamma_g \Delta H / 2 $, $ f_m = \gamma_g\, 4\pi M_s $, and $ \gamma_g = 2.8 $ MHz/Oe. Here $ H_0 = $ `bias_field_oe` is the internal bias field (Oe), $ 4\pi M_s = $ `saturation_magnetization_gauss` (G), and $ \Delta H = $ `linewidth_oe` the resonance linewidth (Oe).
  - `{"type": "magnetized_plasma", "name": ..., "f_p": ..., "collision_frequency": ..., "cyclotron_frequency": ..., "bias_direction": ...}`: a collisional electron plasma in a static magnetic field along $ +y $ or $ -y $ (`bias_direction` is `"+y"` or `"-y"`). Its relative permeability is 1 and its relative permittivity is the tensor

    $$
    [\varepsilon] = \begin{pmatrix} \varepsilon_1 & 0 & j s g \\ 0 & \varepsilon_3 & 0 \\ -j s g & 0 & \varepsilon_1 \end{pmatrix}, \qquad \varepsilon_1 = 1 - \frac{f_p^2\, f_\nu}{f\,(f_\nu^2 - f_c^2)}, \qquad \varepsilon_3 = 1 - \frac{f_p^2}{f\, f_\nu}, \qquad g = \frac{f_p^2\, f_c}{f\,(f_\nu^2 - f_c^2)},
    $$

    with $ f_\nu = f - j\nu / 2\pi $, $ s = +1 $ for `"+y"` and $ s = -1 $ for `"-y"`, plasma frequency $ f_p = $ `f_p` (Hz), collision frequency $ \nu = $ `collision_frequency` (s$ ^{-1} $) and electron cyclotron frequency $ f_c = $ `cyclotron_frequency` (Hz).
- `test_frequencies`: a list of frequencies in hertz.

Every number in the file is an exact decimal, and the guide is defined by these exact values.

## Model

The guide occupies $ 0 < x < a $, $ 0 < y < b $ and is uniform in $ z $. All four walls are perfect electric conductors, and the speed of light is $ c = 299\,792\,458 $ m/s. Fields vary in time as $ e^{+j 2\pi f t} $, so lossy materials have $ \operatorname{Im} \varepsilon < 0 $. Inside each layer they obey Maxwell's equations without sources, $ \nabla \times \mathbf E = -j\omega\mu_0 [\mu]\, \mathbf H $ and $ \nabla \times \mathbf H = j\omega\varepsilon_0 [\varepsilon]\, \mathbf E $ ($ [\mu] = 1 $ and $ [\varepsilon] = \varepsilon $ where the material above is isotropic) ($ \omega = 2\pi f $). The tangential components of $ \mathbf E $ and $ \mathbf H $ are continuous at every interface between layers.

A **mode** at frequency $ f $ is a nonzero source-free field of the guide in which every component varies along the guide as $ e^{-\gamma z} $, with a complex propagation constant $ \gamma = \alpha + j\beta $.

- **Forward** modes decay towards $ +z $: $ \alpha > 0 $.
- **Reverse** modes decay towards $ -z $. A reverse mode is reported by its constant $ \gamma_r = -\gamma $, so that its fields vary as $ e^{+\gamma_r z} $ with $ \gamma_r = \alpha + j\beta $ and $ \alpha > 0 $.

For every reported mode, $ \alpha $ (Np/m) is its **attenuation constant** and $ \beta $ (rad/m) its **phase constant**, with its sign: for a backward wave, whose phase travels opposite to its power, $ \beta < 0 $. The **attenuation** of the mode in dB/m is $ 20 \log_{10}(e)\, \alpha $.

**Families.**
- **Guides without a gyrotropic layer** (see below). Every mode is either an LSE mode ($ E_x = 0 $) or an LSM mode ($ H_x = 0 $), and its fields vary across the height as $ \cos(n\pi y / b) $ or $ \sin(n\pi y / b) $. The families are $ {\rm LSE}_n $ ($ n = 0, 1, 2, \dots $) and $ {\rm LSM}_n $ ($ n = 1, 2, \dots $).
- **Guides with a gyrotropic layer** (a `ferrite` or a `magnetized_plasma` layer). Every mode has fields varying across the height as $ E_x, E_z, H_y \propto \sin(n\pi y / b) $ and $ E_y, H_x, H_z \propto \cos(n\pi y / b) $ for one integer $ n \ge 0 $. The modes with $ n = 0 $ have $ \mathbf E = E_y(x)\, \hat{\mathbf y}\, e^{-\gamma z} $ and form the family $ {\rm LSE}_0 $; for each $ n \ge 1 $ the modes form the family $ {\rm HYB}_n $ (they are in general neither LSE nor LSM). The forward and reverse modes of such a guide are in general different.

Report a mode exactly when its attenuation is at most 200 dB/m, that is $ \alpha \le 200 / (20 \log_{10} e) $ Np/m ($ \approx 23.03 $ Np/m). This includes evanescent modes below their cutoff whose attenuation is still under the limit.

In every guide, at every test frequency and in each direction:
- each mode with attenuation at most 200 dB/m is non-degenerate within its family (a simple zero of the family's dispersion relation), so each mode is reported once;
- no mode has an attenuation within 2 % of 200 dB/m;
- no two modes of the same family have $ \gamma $ values closer than $ 10^{-4} |\gamma| $;
- the power split of every mode (below) satisfies $ \sum_i |P_i| \le 30\, |P| $.

**Power split.** For a reported mode let $ P_i = \operatorname{Re} \iint (\mathbf E \times \mathbf H^*) \cdot \hat{\mathbf z}\; dx\, dy $ over the part of the cross-section occupied by layer $ i $ (at any fixed $ z $), and $ P = \sum_i P_i $ the total. The **power fractions** of the mode are $ p_i = P_i / P $, one per layer in the order of `layers`. They do not depend on the amplitude of the field, they add up to 1, and a fraction is negative where the power in that layer flows against the total. For a reverse mode they are the fractions of its own field (the one varying as $ e^{+\gamma_r z} $).

## Output format

`/root/results/modes.json` must be a regular file (not a symbolic link) of at most 16 MiB containing one UTF-8 JSON object whose only key is `"guides"`. No key may appear twice in any object.

- `"guides"` maps every `id` of the input (all of them, no others) to an object with exactly one key, `modes`.
- `modes` is a list with one entry per test frequency, in the order of `test_frequencies`.
- Each entry is an object with exactly the keys `forward` and `reverse`.
- Each of these maps a family name (`"LSE0"`, `"LSE1"`, ..., `"LSM1"`, `"LSM2"`, ..., `"HYB1"`, `"HYB2"`, ...) to a list of modes, in any order. Each mode is a list `[beta, alpha, power]` with `beta` in rad/m, `alpha` in Np/m ($ \alpha > 0 $) and `power` the list of the mode's power fractions $ p_i $, one number per layer of the guide.
- A family appears exactly when it has at least one mode to report at that frequency in that direction. A direction with no mode to report is the empty object.

Each `forward` or `reverse` object has at most 400 families, and every list at most 1000 modes. Every number must be finite.

## Output requirements

Every requirement must hold for the submission to succeed. The deliverable must be complete in `/root/results/modes.json`; it must not depend on executing any code you write or reading any other file. The file must follow the format above, and for every guide, every test frequency and both directions:

- the families present are exactly the true ones;
- in each family the reported modes pair one-to-one with the true modes, and each pair satisfies $ |\gamma - \gamma^*| \le 10^{-7} |\gamma^*| $, where $ \gamma = \alpha + j\beta $ is the reported value and $ \gamma^* $ the true one;
- the power fractions of each paired mode satisfy $ |p_i - p_i^*| \le 10^{-6} \sum_k |p_k^*| $ for every layer $ i $, where $ p^* $ are the true fractions.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
