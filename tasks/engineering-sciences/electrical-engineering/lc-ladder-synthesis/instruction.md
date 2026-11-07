Write a filter synthesis program, `/root/results/synthesize.py`, that designs a doubly terminated LC ladder filter for a given specification. The program must work for any specification of the kind described below, not only for the examples supplied here.

## Usage

```
python3 /root/results/synthesize.py SPEC_PATH NETWORK_PATH
```

`SPEC_PATH` is the path of a specification file, for example `/root/data/examples/LP-40M.json`. `NETWORK_PATH` is the path of the file the program must write the ladder to, for example `/tmp/LP-40M.network.json`.

The program must work under the following conditions:

- it is run from an empty working directory, as an unprivileged user, with no network access;
- it handles one specification per run, within 300 s; it and every thread or process it starts run on a single CPU core, and each process has a 3 GB address-space (virtual memory) limit;
- it uses only Python 3.11, the standard library, `numpy` 2.2.6, `scipy` 1.15.3, `sympy` 1.13.3 and `mpmath` 1.3.0;
- it reads no file other than `SPEC_PATH`.

## Specification

The specification is a UTF-8 JSON object with these keys:

- `id`, `description`: strings.
- `source_resistance_ohm`: $ R_s $.
- `passband_hz`: $ [f_1, f_2] $. A low-pass filter has $ f_1 = 0 $; a band-pass filter has $ 0 < f_1 < f_2 $.
- `return_loss_db`: the passband return loss $ RL $ (dB).
- `degree`: $ n $. It is even for a band-pass filter.
- `transmission_zeros`: an object with two keys:
  - `dc`: the number $ a $ of transmission zeros at $ f = 0 $ (0 for a low-pass filter, at least 1 for a band-pass filter);
  - `stopbands`: a list of at most two stopbands, each an object with the keys `edge_hz` and `zeros`. A stopband with `edge_hz` $ = f_e < f_1 $ is the band $ 0 < f \le f_e $; a stopband with $ f_e > f_2 $ is the band $ f \ge f_e $. There is at most one stopband on each side of the passband. `zeros` is the number $ m_b \ge 1 $ of finite transmission zeros that lie in that stopband. Their frequencies are not given; they are fixed by the stopband conditions under "The required response".

  With $ m $ the total number of finite transmission zeros, the remaining $ n - a - 2m \ge 1 $ transmission zeros are at infinity.

Examples are in `/root/data/examples/*.json`. `/root/data/examples/reference_insertion_loss.json` gives the reference insertion loss (dB) of each example at sample frequencies, as `[f_hz, loss_db]` pairs.

The program must also handle other specifications of the same kind: low-pass and band-pass filters of degree 5 to 16, band-pass filters with 1 to 4 transmission zeros at dc, up to 3 finite transmission zeros in total, passbands from audio frequencies to GHz, and band-pass bandwidths from 0.1 % of the centre frequency to more than three decades ($ f_2 / f_1 $ up to 3500).

## The required response

Let $ s = j\omega $, with $ \omega = 2\pi f $. The transducer power gain is $ G = P_{\text{load}} / P_{\text{available}} $, with

$$
G(\omega) = \frac{1}{1 + |K(j\omega)|^2}, \qquad K(s) = \frac{F(s)}{P(s)}, \qquad P(s) = s^a \prod_{k=1}^{m} (s^2 + \omega_k^2), \quad \omega_k = 2\pi f_k .
$$

$ F $ is the generalized Chebyshev (equiripple) polynomial of degree $ n $:

$$
F(s) = c\, s^r \prod_{i=1}^{d} (s^2 + \omega_i^2), \qquad r = n \bmod 2, \quad d = (n - r)/2 .
$$

- All $ \omega_i $ lie in the passband.
- $ c > 0 $ and the $ \omega_i $ are such that the maximum of $ |K(j\omega)| $ over the passband equals $ \varepsilon = (10^{RL/10} - 1)^{-1/2} $.
- That maximum is attained at $ d + 1 $ passband frequencies.

The finite transmission zeros $ f_k $ are placed so that the attenuation is also equiripple in every stopband:

- each $ f_k $ lies inside its stopband;
- in a stopband with $ m_b $ zeros, the minimum of $ |K(j\omega)| $ over the stopband is attained at exactly $ m_b + 1 $ frequencies: at the band edge $ f_e $, once between each pair of consecutive zeros of that stopband, and once between the zero farthest from the passband and $ f = 0 $ (lower stopband) or $ f = \infty $ (upper stopband).

## The ladder

The ladder is a cascade of arms between the source, a voltage source with series resistance $ R_s $, and the load resistance $ R_L $. Each arm is either `series` or `shunt` and is one of:

- `L`: an inductor;
- `C`: a capacitor;
- `series_LC`: an inductor in series with a capacitor;
- `parallel_LC`: an inductor in parallel with a capacitor.

There are no transformers, coupled coils or resistors other than $ R_s $ and $ R_L $. The ladder must satisfy all of the following:

1. Its insertion loss $ -10 \log_{10} G $ agrees with that of the required response to within $ 10^{-5} $ dB at every frequency from $ f_2 $ down to $ f_1/30 $ (band-pass) or to dc (low-pass) and from $ f_2 $ up to $ 30 f_2 $, including the stopband edges and the frequencies next to the finite transmission zeros.
2. It has at most $ n + m $ inductors and capacitors (a canonical ladder).
3. $ R_s/2 \le R_L \le 2R_s $.

## Output

`NETWORK_PATH` must be a regular file (not a symbolic link) of at most 1 MiB, holding a UTF-8 JSON object with exactly two keys:

- `load_resistance_ohm`: the number $ R_L $.
- `arms`: a non-empty list of arms in order from the source to the load. Each arm is an object with exactly these keys:
  - `position`: `"series"` or `"shunt"`;
  - `type`: one of the four arm types above;
  - `L` (henry) for the arm types that contain an inductor, and `C` (farad) for those that contain a capacitor.

The values of `load_resistance_ohm`, `L` and `C` must be positive finite numbers.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
