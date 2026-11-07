"""Held-out guide files: tests/data/run_a ... run_d/heldout_guides.json.

usage: python3 gen_heldout.py tests/data

Six guides in four files (one or two guides each), from devices of the kinds the example guides do not show:
  H1  magnetostatic-wave delay line: a 100 um single-crystal YIG film grown by liquid-phase epitaxy on a 0.5 mm
      gadolinium gallium garnet (GGG) substrate, in a shielding enclosure of WR-42 cross-section far below the
      cutoff of its empty-guide modes; internal bias 300 Oe along +y.
  H2  the same film with a ground plane 0.25 mm from one face, WR-34 cross-section, bias along -y.
  H3  WR-284 loaded with two 5 mm zirconium tin titanate ceramic slabs.
  H4  WR-137 with a low-pressure argon discharge cell (1 Pa) between 4 mm fused-quartz plates.
  H5  a magnetised low-pressure argon discharge (500 G along -y) behind a 4 mm fused-quartz window, WR-62
      cross-section.
  H6  a WR-112 C-band ferrite phase shifter (Trans-Tech G-113 YIG, 600 Oe along -y, on RO3006).
Materials:
  LPE YIG film: 4 pi Ms 1780 G, eps' 15, tan d 2e-4 (as the polycrystal), linewidth 1 Oe (LPE films 0.3-1 Oe;
    Stancil and Prabhakar, Spin Waves, 2009; Ishak, Proc. IEEE 76 (1988) 171).
  GGG: eps' 12.2, tan d 3e-4.
  Zr-Sn-Ti oxide ceramic: eps' 37, tan d 2e-4 (Q ~ 5000 at 3-4 GHz).
  RO3006 6.15 / 0.002 (Rogers data sheet); fused quartz 3.78 / 1e-4 (Pozar); air 1 / 0.
  Argon at 1 Pa and 300 K: n_g = 2.4e20 m^-3, nu = n_g K_m = 2.4e7 s^-1 (K_m ~ 1e-13 m^3/s; Lieberman and
    Lichtenberg); n_e = 7.9e17 m^-3 gives f_p = 8.98 kHz sqrt(n_e / cm^-3) = 8 GHz (NRL Plasma Formulary), the
    density of an inductively coupled or ECR source at this pressure; f_c = 2.799 MHz/G x 500 G = 1.4 GHz.
  Trans-Tech G-113 YIG: 4 pi Ms 1780 G, linewidth 45 Oe, eps' 15, tan d 2e-4.
Waveguide cross-sections: EIA RS-261-B (WR-34, WR-42, WR-62, WR-112, WR-137, WR-284).
The materials and geometries come from the devices; the test frequencies were chosen inside the surface-wave
bands such that every mode satisfies the guarantees stated in the instruction (no attenuation within 2 % of
200 dB/m in particular; see the README for the frequencies that were rejected).
"""
import json
import os
import sys


def diel(name, eps, tand):
    return lambda d: {'thickness': d, 'material': {'type': 'dielectric', 'name': name, 'eps': eps, 'loss_tangent': tand}}


AIR = diel('air', 1, 0)
QUARTZ = diel('fused quartz', 3.78, 0.0001)
GGG = diel('gadolinium gallium garnet (GGG)', 12.2, 0.0003)
CERAMIC = diel('zirconium tin titanate ceramic', 37, 0.0002)
RO3006 = diel('RO3006', 6.15, 0.002)
LPE = 'single-crystal YIG film (liquid-phase epitaxy)'


def ferrite(d, name, H0, dH, sg):
    return {'thickness': d, 'material': {'type': 'ferrite', 'name': name, 'eps': 15, 'loss_tangent': 0.0002,
                                         'saturation_magnetization_gauss': 1780, 'linewidth_oe': dH,
                                         'bias_field_oe': H0, 'bias_direction': sg}}


def plasma(d, name, fp, nu):
    return {'thickness': d, 'material': {'type': 'plasma', 'name': name, 'f_p': fp, 'collision_frequency': nu}}


def mplasma(d, name, fp, nu, fc, sg):
    return {'thickness': d, 'material': {'type': 'magnetized_plasma', 'name': name, 'f_p': fp,
                                         'collision_frequency': nu, 'cyclotron_frequency': fc, 'bias_direction': sg}}


def guide(gid, desc, a, b, layers, freqs):
    assert abs(sum(L['thickness'] for L in layers) - a) < 1e-12, gid
    return {'id': gid, 'description': desc, 'a': a, 'b': b, 'layers': layers, 'test_frequencies': freqs}


RUNS = {
    'run_a': [guide('H1', 'magnetostatic-wave delay line: 100 um YIG film (300 Oe, +y) on a 0.5 mm GGG substrate, '
                    '3 mm from a side wall of a WR-42-size enclosure',
                    0.010668, 0.004318,
                    [AIR(0.003), ferrite(0.0001, LPE, 300, 1.0, '+y'), GGG(0.0005), AIR(0.007068)],
                    [2750000000, 2850000000])],
    'run_b': [guide('H2', 'magnetostatic-wave delay line with a ground plane: 100 um YIG film (300 Oe, -y) 0.25 mm '
                    'from a side wall, on a 0.5 mm GGG substrate, WR-34-size enclosure',
                    0.008636, 0.004318,
                    [AIR(0.00025), ferrite(0.0001, LPE, 300, 1.0, '-y'), GGG(0.0005), AIR(0.007786)],
                    [2750000000, 2850000000])],
    'run_c': [guide('H3', 'WR-284 loaded with two 5 mm titanate ceramic slabs, one against a side wall, one 20 mm from it',
                    0.072136, 0.034036,
                    [CERAMIC(0.005), AIR(0.02), CERAMIC(0.005), AIR(0.042136)],
                    [3000000000, 3600000000]),
              guide('H4', 'WR-137 low-pressure discharge cell: 5 mm argon at 1 Pa between 4 mm fused-quartz plates, centred',
                    0.034849, 0.015799,
                    [AIR(0.0109245), QUARTZ(0.004),
                     plasma(0.005, 'argon discharge, n_e = 7.9e17 m^-3, 1 Pa', 8000000000, 24000000),
                     QUARTZ(0.004), AIR(0.0109245)],
                    [3500000000, 3600000000])],
    'run_d': [guide('H5', 'magnetised low-pressure discharge: 6 mm argon (500 G, -y) against a side wall behind a 4 mm '
                    'fused-quartz window, WR-62-size guide',
                    0.015799, 0.007899,
                    [mplasma(0.006, 'magnetised argon discharge, n_e = 7.9e17 m^-3, 1 Pa, B0 = 500 G',
                             8000000000, 24000000, 1400000000, '-y'),
                     QUARTZ(0.004), AIR(0.005799)],
                    [3500000000, 3800000000, 4000000000]),
              guide('H6', 'WR-112 ferrite phase shifter: 2.54 mm RO3006 against a side wall, then 2.5 mm YIG (600 Oe, -y)',
                    0.028499, 0.012624,
                    [RO3006(0.00254), ferrite(0.0025, 'YIG (Trans-Tech G-113)', 600, 45, '-y'), AIR(0.023459)],
                    [7500000000, 9000000000])],
}


def main():
    root = sys.argv[1]
    for run, guides in RUNS.items():
        os.makedirs(os.path.join(root, run), exist_ok=True)
        with open(os.path.join(root, run, 'heldout_guides.json'), 'w') as fh:
            json.dump({'guides': guides}, fh, indent=1)


if __name__ == '__main__':
    main()
