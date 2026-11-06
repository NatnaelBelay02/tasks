"""Input data: the 18 guides of guides.json.

usage: python3 gen_guides.py guides.json

Materials (relative permittivity and loss tangent at 10 GHz):
  RT/duroid 6002 2.94 / 0.0012, RO4350B 3.48 / 0.0037, RO3006 6.15 / 0.0020, RO3010 10.2 / 0.0022
  (Rogers Corporation data sheets); PTFE 2.08 / 0.0004 and fused quartz 3.78 / 0.0001 (Pozar,
  Microwave Engineering, table of dielectric constants and loss tangents); air 1 / 0.
Gas discharges: argon at 6 to 33 Pa (0.05 to 0.25 Torr) with electron densities of 2.5e17 to
  2.4e18 m^-3, in a planar cell between two fused-quartz plates.  f_p = 8.98 kHz sqrt(n_e / cm^-3)
  (NRL Plasma Formulary); nu = n_g K_m with K_m ~ 1e-13 m^3/s for argon at T_e of a few eV
  (Lieberman and Lichtenberg).  The README lists n_e and the pressure of every discharge.
Waveguides: WR-62, WR-90, WR-112, WR-137, WR-187, WR-284 inner dimensions (EIA RS-261-B).
Every thickness is a stock laminate or plate thickness or the remainder of the guide width, and
every test frequency a round frequency; none is tuned to a root.
"""

import json
import sys
from decimal import Decimal as Dm

WR = {'WR-62': ('0.015799', '0.007899'), 'WR-90': ('0.02286', '0.01016'), 'WR-112': ('0.028499', '0.012624'),
      'WR-137': ('0.034849', '0.015799'), 'WR-187': ('0.047549', '0.022149'), 'WR-284': ('0.072136', '0.034036')}

MAT = {
    'air': ('air', '1', '0'),
    'quartz': ('fused quartz', '3.78', '0.0001'),
    'ptfe': ('PTFE', '2.08', '0.0004'),
    '6002': ('RT/duroid 6002', '2.94', '0.0012'),
    '4350': ('RO4350B', '3.48', '0.0037'),
    '3006': ('RO3006', '6.15', '0.0020'),
    '3010': ('RO3010', '10.2', '0.0022'),
}


def diel(t, key):
    name, e, tand = MAT[key]
    return {'thickness': t, 'material': {'type': 'dielectric', 'name': name, 'eps': e, 'loss_tangent': tand}}


def yig(t, h0, sign):
    """Polycrystalline YIG (Trans-Tech G-113): 4 pi Ms = 1780 G, linewidth 45 Oe, eps' 15.0, tan d 0.0002;
    h0 is the internal bias field (Oe), magnetised along +y or -y."""
    return {'thickness': t, 'material': {'type': 'ferrite', 'name': 'YIG (Trans-Tech G-113)', 'eps': '15.0',
            'loss_tangent': '0.0002', 'saturation_magnetization_gauss': '1780', 'linewidth_oe': '45',
            'bias_field_oe': h0, 'bias_direction': sign}}


def argon(t, fp, nu, ne, pa):
    return {'thickness': t, 'material': {'type': 'plasma',
            'name': 'argon discharge, n_e = %s m^-3, %s Pa' % (ne, pa), 'f_p': fp, 'collision_frequency': nu}}


def magargon(t, fp, nu, fc, sign, ne, pa, gauss):
    """Argon discharge in a static magnetic field B0 = gauss (G) along +y or -y: electron cyclotron
    frequency f_c = e B0 / (2 pi m_e) = 2.799 MHz/G x B0, rounded to the stated value."""
    return {'thickness': t, 'material': {'type': 'magnetized_plasma',
            'name': 'magnetised argon discharge, n_e = %s m^-3, %s Pa, B0 = %s G' % (ne, pa, gauss), 'f_p': fp,
            'collision_frequency': nu, 'cyclotron_frequency': fc, 'bias_direction': sign}}


def rest(wr, layers):
    """Thickness that fills the guide width exactly."""
    a = Dm(WR[wr][0])
    return str(a - sum(Dm(L['thickness']) for L in layers if L is not None))


def g(gid, wr, desc, layers, tests):
    a, b = WR[wr]
    assert sum(Dm(L['thickness']) for L in layers) == Dm(a), gid
    return {'id': gid, 'description': '%s: %s' % (wr, desc), 'a': a, 'b': b, 'layers': layers,
            'test_frequencies': tests}


def guides():
    G = []
    L = [diel('0.003', '6002')]
    L.append(diel(rest('WR-90', L), 'air'))
    G.append(g('G01', 'WR-90', 'dielectric-slab phase shifter, 3 mm RT/duroid 6002 slab against a side wall',
               L, ['9e9', '13e9', '16e9', '18e9']))
    G.append(g('G02', 'WR-90', 'centred H-plane slab, 2.54 mm RO3010 slab in the middle of the guide',
               [diel('0.01016', 'air'), diel('0.00254', '3010'), diel('0.01016', 'air')], ['8e9', '13e9', '15e9', '16e9']))
    G.append(g('G03', 'WR-90', 'dual-slab phase shifter, 1.27 mm RO3006 slabs against both side walls',
               [diel('0.00127', '3006'), diel('0.02032', 'air'), diel('0.00127', '3006')], ['12.75e9', '13.5e9', '16e9', '18e9']))
    L = [diel('0.001524', '4350'), diel('0.002286', 'air')] * 5
    L.append(diel(rest('WR-90', L), 'air'))
    G.append(g('G04', 'WR-90', 'stacked-laminate loading, five 1.524 mm RO4350B sheets with 2.286 mm air gaps',
               L, ['10e9', '14e9', '15e9', '18e9']))
    L = [diel('0.00635', 'ptfe')]
    L.append(diel(rest('WR-62', L), 'air'))
    G.append(g('G05', 'WR-62', 'partially filled guide, 6.35 mm PTFE slab against a side wall', L, ['12e9', '16.5e9', '20e9', '22e9']))
    L = [diel('0.002', '3010'), None, diel('0.002', '3010')]
    L[1] = diel(rest('WR-137', L), 'air')
    G.append(g('G06', 'WR-137', 'dual-slab loaded guide, 2 mm RO3010 slabs against both side walls', L, ['4.15e9', '6e9', '8.4e9', '10e9']))
    # gas-discharge loaded guides: planar discharge cell, quartz plate / argon / quartz plate
    G.append(g('G07', 'WR-90', 'plasma switch, 4 mm argon discharge cell between 1 mm quartz plates in the middle of the guide',
               [diel('0.00843', 'air'), diel('0.001', 'quartz'), argon('0.004', '14e9', '5e8', '2.4e18', '21'),
                diel('0.001', 'quartz'), diel('0.00843', 'air')], ['7.2e9', '7.6e9', '8.3e9', '11e9']))
    G.append(g('G08', 'WR-137', 'plasma phase shifter, 7 mm argon discharge cell between 1.5 mm quartz plates in the middle of the guide',
               [diel('0.0124245', 'air'), diel('0.0015', 'quartz'), argon('0.007', '9e9', '2e8', '1.0e18', '8'),
                diel('0.0015', 'quartz'), diel('0.0124245', 'air')], ['4.4e9', '4.8e9', '5.6e9', '8e9']))
    G.append(g('G09', 'WR-284', 'S-band plasma limiter, 9 mm argon discharge cell between 1.5 mm quartz plates in the middle of the guide',
               [diel('0.030068', 'air'), diel('0.0015', 'quartz'), argon('0.009', '4.5e9', '8e8', '2.5e17', '33'),
                diel('0.0015', 'quartz'), diel('0.030068', 'air')], ['2.6e9', '2.75e9', '2.9e9', '3.5e9']))
    L = [argon('0.006', '7e9', '3e8', '6.1e17', '12'), diel('0.0015', 'quartz'), diel('0.004', '3006')]
    L.append(diel(rest('WR-187', L), 'air'))
    G.append(g('G10', 'WR-187', 'plasma-tuned dielectric slab, 6 mm argon discharge against a side wall behind a 1.5 mm quartz window, then 4 mm RO3006',
               L, ['3.5e9', '4.75e9', '5.25e9', '6.2e9']))
    G.append(g('G11', 'WR-112', 'low-pressure plasma cell, 4.5 mm argon discharge between 1.5 mm quartz plates in the middle of the guide',
               [diel('0.0105', 'air'), diel('0.0015', 'quartz'), argon('0.0045', '10e9', '1.5e8', '1.24e18', '6'),
                diel('0.0015', 'quartz'), diel('0.010499', 'air')], ['5.5e9', '5.8e9', '7e9', '9e9']))
    L = [diel('0.004', '3010'), diel('0.001', 'quartz'), argon('0.004', '13e9', '4e8', '2.1e18', '17'), diel('0.001', 'quartz')]
    L.append(diel(rest('WR-90', L), 'air'))
    G.append(g('G12', 'WR-90', 'plasma surface-wave guide, 4 mm RO3010 slab against a side wall, then a 4 mm argon discharge cell between 1 mm quartz plates',
               L, ['5e9', '5.15e9', '9e9', '11.4e9']))
    # ferrite-loaded guides (y-biased YIG slabs): nonreciprocal
    L = [diel('0.003', 'air'), yig('0.002', '800', '+y'), diel('0.00254', '3010')]
    L.append(diel(rest('WR-90', L), 'air'))
    G.append(g('G13', 'WR-90', 'dielectric-loaded ferrite phase shifter, 2 mm YIG slab (800 Oe, +y) 3 mm from a side wall, backed by 2.54 mm RO3010',
               L, FT['G13']))
    L = [diel('0.0045', 'air'), yig('0.002', '400', '+y'), None, yig('0.002', '400', '-y'), diel('0.0045', 'air')]
    L[2] = diel(rest('WR-90', L), 'air')
    G.append(g('G14', 'WR-90', 'twin-slab ferrite phase shifter, two 2 mm YIG slabs magnetised in opposite directions (400 Oe, +y and -y)',
               L, FT['G14']))
    L = [diel('0.005', 'air'), yig('0.001', '3200', '+y')]
    L.append(diel(rest('WR-90', L), 'air'))
    G.append(g('G15', 'WR-90', 'resonance isolator, 1 mm YIG slab near the quarter-width point, biased near ferrimagnetic resonance (3200 Oe, +y)',
               L, FT['G15']))
    L = [diel('0.005', 'air'), yig('0.003', '300', '+y'), diel('0.003', '3006')]
    L.append(diel(rest('WR-137', L), 'air'))
    G.append(g('G16', 'WR-137', 'C-band ferrite phase shifter, 3 mm YIG slab (300 Oe, +y) 5 mm from a side wall, backed by 3 mm RO3006',
               L, FT['G16']))
    # magnetised discharges (static B0 along y): gyroelectric, nonreciprocal hybrid modes
    L = [diel('0.004', 'air'), diel('0.001', 'quartz'), magargon('0.004', '11e9', '3e8', '4.2e9', '+y', '1.5e18', '12', '1500'),
         diel('0.001', 'quartz')]
    L.append(diel(rest('WR-90', L), 'air'))
    G.append(g('G17', 'WR-90', 'magnetised-plasma phase shifter, 4 mm argon discharge between 1 mm quartz plates, 4 mm from a side wall, in a 1500 G field along +y',
               L, ['7.5e9', '8.5e9', '10.5e9', '12.5e9']))
    L = [magargon('0.006', '7e9', '2e8', '2.8e9', '-y', '6.1e17', '8', '1000'), diel('0.0015', 'quartz'), diel('0.003', '3006')]
    L.append(diel(rest('WR-137', L), 'air'))
    G.append(g('G18', 'WR-137', 'magnetised plasma against a side wall, 6 mm argon discharge behind a 1.5 mm quartz window, then 3 mm RO3006, in a 1000 G field along -y',
               L, ['4.5e9', '5.5e9', '6.5e9', '8.2e9']))
    return G


FT = {'G13': ['8e9', '9e9', '10.5e9', '12e9'], 'G14': ['8.5e9', '10e9', '11.5e9', '12e9'],
      'G15': ['9e9', '10.25e9', '11e9', '11.5e9'], 'G16': ['5e9', '6e9', '7e9', '8e9']}


def num(x):
    s = format(Dm(x), 'f')
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    return s


def dumps(gs):
    """JSON with every number written as the exact decimal it was given as."""
    out = []
    for gg in gs:
        lay = []
        for L in gg['layers']:
            m = L['material']
            if m['type'] == 'dielectric':
                mat = '{"type": "dielectric", "name": %s, "eps": %s, "loss_tangent": %s}' % (
                    json.dumps(m['name']), num(m['eps']), num(m['loss_tangent']))
            elif m['type'] == 'plasma':
                mat = '{"type": "plasma", "name": %s, "f_p": %s, "collision_frequency": %s}' % (
                    json.dumps(m['name']), num(m['f_p']), num(m['collision_frequency']))
            elif m['type'] == 'magnetized_plasma':
                mat = ('{"type": "magnetized_plasma", "name": %s, "f_p": %s, "collision_frequency": %s, '
                       '"cyclotron_frequency": %s, "bias_direction": "%s"}') % (
                    json.dumps(m['name']), num(m['f_p']), num(m['collision_frequency']), num(m['cyclotron_frequency']),
                    m['bias_direction'])
            else:
                mat = ('{"type": "ferrite", "name": %s, "eps": %s, "loss_tangent": %s, "saturation_magnetization_gauss": %s, '
                       '"linewidth_oe": %s, "bias_field_oe": %s, "bias_direction": "%s"}') % (
                    json.dumps(m['name']), num(m['eps']), num(m['loss_tangent']), num(m['saturation_magnetization_gauss']),
                    num(m['linewidth_oe']), num(m['bias_field_oe']), m['bias_direction'])
            lay.append('    {"thickness": %s, "material": %s}' % (num(L['thickness']), mat))
        out.append('  {"id": "%s", "description": %s, "a": %s, "b": %s,\n   "layers": [\n%s\n   ],\n   "test_frequencies": [%s]}'
                   % (gg['id'], json.dumps(gg['description']), num(gg['a']), num(gg['b']), ',\n'.join(lay),
                      ', '.join(num(x) for x in gg['test_frequencies'])))
    return '{"guides": [\n' + ',\n'.join(out) + '\n]}\n'


if __name__ == '__main__':
    with open(sys.argv[1], 'w') as fh:
        fh.write(dumps(guides()))
