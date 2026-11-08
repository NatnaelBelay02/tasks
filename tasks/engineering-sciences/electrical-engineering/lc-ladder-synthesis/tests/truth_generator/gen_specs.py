"""Write the filter specifications.

usage: python3 gen_specs.py examples_dir hidden_dir

Every specification is an application filter with its passband, required return loss,
terminations and degree, and with stopbands that begin at the frequencies the filter has to
reject: transmitter harmonics, the alias band of a converter, or a neighbouring service.
The number of finite transmission zeros in each stopband is part of the specification;
their frequencies are not (the program has to place them).  Nothing is fitted.
"""
import json
import os
import sys

M = 1e6
K = 1e3

# id, description, Rs, passband, return loss dB, degree, dc zeros, stopbands [(edge Hz, number of finite zeros)]
SPECS = [
    # ---- examples (given to the solver) ----
    ('ex', 'LP-40M', '7.0-7.3 MHz amateur transmitter: harmonic low-pass, stopband from the 2nd harmonic of 7.0 MHz',
     50, (0, 7.3 * M), 20, 5, 0, [(14.0 * M, 2)]),
    ('ex', 'LP-433', '433.92 MHz ISM transmitter: harmonic low-pass, stopband from the 2nd harmonic',
     50, (0, 450 * M), 20, 6, 0, [(867.84 * M, 2)]),
    ('ex', 'LP-AUDIO48', 'anti-alias filter of a 48 kHz audio converter: 20 kHz passband, stopband from 28 kHz (first frequency aliased into the passband)',
     600, (0, 20 * K), 26, 7, 0, [(28 * K, 2)]),
    ('ex', 'BP-455K', '455 kHz AM IF selectivity filter, 10 kHz bandwidth, no finite transmission zeros',
     50, (450 * K, 460 * K), 20, 8, 4, []),
    ('ex', 'BP-PLC-A', 'G3-PLC CENELEC-A coupling filter, 35.9-90.6 kHz, upper stopband from 148.5 kHz (upper limit of the CENELEC bands)',
     50, (35.9 * K, 90.6 * K), 16, 8, 3, [(148.5 * K, 1)]),
    ('ex', 'BP-AIR', 'VHF airband receiver preselector, 118-137 MHz, stopbands below 108 MHz (FM broadcast) and above 156.8 MHz (marine channel 16)',
     50, (118 * M, 137 * M), 20, 10, 3, [(108 * M, 1), (156.8 * M, 1)]),
    # ---- hidden (graded) ----
    ('hid', 'LP-6M78', '6.78 MHz wireless-power class-E amplifier: harmonic low-pass, stopband from the 2nd harmonic',
     50, (0, 7.0 * M), 20, 7, 0, [(13.56 * M, 2)]),
    ('hid', 'LP-13M56', '13.56 MHz plasma RF generator: harmonic low-pass, stopband from the 2nd harmonic',
     50, (0, 14.0 * M), 26, 9, 0, [(27.12 * M, 3)]),
    ('hid', 'LP-27M12', '27.12 MHz ISM generator: harmonic low-pass, stopband from the 2nd harmonic',
     50, (0, 28.0 * M), 20, 8, 0, [(54.24 * M, 2)]),
    ('hid', 'LP-CATV42', 'cable-TV diplexer, 5-42 MHz upstream low-pass arm (75 ohm), stopband from 54 MHz (lower edge of channel 2)',
     75, (0, 42 * M), 20, 9, 0, [(54 * M, 3)]),
    ('hid', 'LP-VIDEO', 'BT.601 luminance reconstruction filter (75 ohm), 5.75 MHz passband, stopband from 7.75 MHz (first image)',
     75, (0, 5.75 * M), 20, 8, 0, [(7.75 * M, 2)]),
    ('hid', 'LP-POTS', 'ADSL splitter low-pass (600 ohm), 4 kHz voice band, stopband from 25.875 kHz (first ADSL tone)',
     600, (0, 4 * K), 14, 5, 0, [(25.875 * K, 2)]),
    ('hid', 'LP-2M', '144-148 MHz amateur transmitter: harmonic low-pass, stopband from the 2nd harmonic of 144 MHz',
     50, (0, 148 * M), 20, 6, 0, [(288 * M, 2)]),
    ('hid', 'LP-868', '868 MHz LoRa transmitter: low-pass, stopband from 1215 MHz (start of the GNSS L-bands)',
     50, (0, 870 * M), 18, 8, 0, [(1215 * M, 3)]),
    ('hid', 'BP-10M7', '10.7 MHz FM IF filter, 300 kHz bandwidth, stopbands below 10.25 MHz and above 11.15 MHz',
     50, (10.55 * M, 10.85 * M), 20, 10, 3, [(10.25 * M, 1), (11.15 * M, 1)]),
    ('hid', 'BP-FM', 'FM broadcast-band preselector, 87.5-108 MHz, stopbands below 72 MHz (top of TV channel 4) and above 118 MHz (VHF airband)',
     50, (87.5 * M, 108 * M), 20, 10, 3, [(72 * M, 1), (118 * M, 1)]),
    ('hid', 'BP-MARINE', 'marine VHF receiver preselector, 156-163 MHz, stopbands below 148 MHz (2 m amateur band) and above 174 MHz (TV band III)',
     50, (156 * M, 163 * M), 20, 12, 4, [(148 * M, 1), (174 * M, 1)]),
    ('hid', 'BP-DCF77', 'DCF77 time-signal receiver, 76.5-78.5 kHz, stopbands below 75 kHz (two zeros) and above 81 kHz',
     50, (76.5 * K, 78.5 * K), 20, 12, 1, [(75 * K, 2), (81 * K, 1)]),
    ('hid', 'BP-ISM24', '2.4 GHz ISM-band filter, 2400-2483.5 MHz, stopbands below 2300 MHz (LTE band 40) and above 2570 MHz (LTE band 38)',
     50, (2400 * M, 2483.5 * M), 20, 8, 2, [(2300 * M, 1), (2570 * M, 1)]),
    ('hid', 'BP-915', '902-928 MHz receiver, stopbands below 894 MHz (GSM-850 downlink) and above 935 MHz (GSM-900 downlink)',
     50, (902 * M, 928 * M), 20, 10, 3, [(894 * M, 1), (935 * M, 1)]),
    ('hid', 'BP-20M', '14.0-14.35 MHz amateur receiver, stopbands below 13.6 MHz and above 14.8 MHz',
     50, (14.0 * M, 14.35 * M), 20, 10, 3, [(13.6 * M, 1), (14.8 * M, 1)]),
    ('hid', 'BP-GSM900', 'GSM-900 downlink receiver, 925-960 MHz, stopbands below 915 MHz (uplink) and above 980 MHz',
     50, (925 * M, 960 * M), 20, 12, 3, [(915 * M, 1), (980 * M, 2)]),
    ('hid', 'BP-LTE20', 'LTE band 20 downlink receiver, 791-821 MHz, stopbands below 758 MHz (700 MHz band) and above 832 MHz (band 20 uplink)',
     50, (791 * M, 821 * M), 20, 12, 3, [(758 * M, 2), (832 * M, 1)]),
    ('hid', 'BP-CB', '26.965-27.405 MHz CB-band transmitter filter, no finite transmission zeros',
     50, (26.965 * M, 27.405 * M), 20, 8, 4, []),
    ('hid', 'BP-MW', 'AM medium-wave receiver preselector, 530-1700 kHz, stopbands below 283.5 kHz (top of the long-wave band) and above 2.3 MHz (120 m tropical band)',
     50, (530 * K, 1700 * K), 20, 8, 1, [(283.5 * K, 1), (2.3 * M, 2)]),
    ('hid', 'BP-455ADJ', '455 kHz AM IF filter with adjacent-channel rejection, 450-460 kHz, stopbands below 440 kHz (two zeros) and above 470 kHz',
     50, (450 * K, 460 * K), 20, 12, 1, [(440 * K, 2), (470 * K, 1)]),
    ('hid', 'BP-2M', '144-148 MHz amateur receiver, stopbands below 137 MHz (weather satellites, two zeros) and above 156 MHz (marine VHF)',
     50, (144 * M, 148 * M), 20, 12, 1, [(137 * M, 2), (156 * M, 1)]),
    ('hid', 'BP-15M', '21.0-21.45 MHz amateur receiver, stopbands below 20.5 MHz (two zeros) and above 22.0 MHz',
     50, (21.0 * M, 21.45 * M), 20, 12, 1, [(20.5 * M, 2), (22.0 * M, 1)]),
    ('hid', 'BP-10M', '28.0-29.7 MHz amateur receiver, stopbands below 27.0 MHz (CB band, two zeros) and above 30.5 MHz',
     50, (28.0 * M, 29.7 * M), 20, 12, 1, [(27.0 * M, 2), (30.5 * M, 1)]),
    ('hid', 'BP-UMTS', 'UMTS band 1 downlink receiver, 2110-2170 MHz, stopbands below 2025 MHz (band 34, two zeros) and above 2200 MHz (satellite S-band)',
     50, (2110 * M, 2170 * M), 20, 12, 1, [(2025 * M, 2), (2200 * M, 1)]),
    # narrow IF filters of degree 14 and 16 (the LC prototypes of crystal and ceramic ladder filters)
    ('hid', 'BP-IF45M', '45 MHz first IF of a VHF receiver, 150 kHz bandwidth, stopbands below 44.8 MHz and above 45.2 MHz',
     50, (44.925 * M, 45.075 * M), 20, 14, 4, [(44.8 * M, 1), (45.2 * M, 1)]),
    ('hid', 'BP-IF70M5', '70 MHz IF, 500 kHz bandwidth, stopband above 71 MHz (two zeros)',
     50, (69.75 * M, 70.25 * M), 20, 14, 2, [(71 * M, 2)]),
    ('hid', 'BP-SSB455', '455 kHz SSB IF filter, 2.4 kHz bandwidth, stopbands below 452 kHz (carrier and opposite sideband, two zeros) and above 457.5 kHz',
     50, (453.6 * K, 456.0 * K), 20, 16, 2, [(452 * K, 2), (457.5 * K, 1)]),
    ('hid', 'BP-CW455', '455 kHz CW IF filter, 1.2 kHz bandwidth, stopbands below 453.6 kHz (two zeros) and above 456.4 kHz',
     50, (454.4 * K, 455.6 * K), 20, 16, 2, [(453.6 * K, 2), (456.4 * K, 1)]),
    ('hid', 'BP-IF10M7N', '10.7 MHz IF of a narrow-band FM data receiver, 40 kHz bandwidth, stopbands below 10.65 MHz and above 10.75 MHz (two zeros)',
     50, (10.68 * M, 10.72 * M), 20, 16, 4, [(10.65 * M, 1), (10.75 * M, 2)]),
    ('hid', 'BP-IF10M7C', '10.7 MHz IF of a 25 kHz-channel receiver, 25 kHz bandwidth, stopband above 10.73 MHz (upper adjacent channel, two zeros)',
     50, (10.6875 * M, 10.7125 * M), 20, 16, 4, [(10.73 * M, 2)]),
    ('hid', 'BP-IF21M4', '21.4 MHz first IF of a land-mobile receiver, 50 kHz bandwidth, stopbands below 21.35 MHz and above 21.45 MHz (two zeros)',
     50, (21.375 * M, 21.425 * M), 20, 16, 4, [(21.35 * M, 1), (21.45 * M, 2)]),
    ('hid', 'BP-IF70M2', '70 MHz satellite-modem IF, 200 kHz bandwidth, stopbands below 69.8 MHz and above 70.2 MHz (two zeros)',
     50, (69.9 * M, 70.1 * M), 20, 16, 4, [(69.8 * M, 1), (70.2 * M, 2)]),
    ('hid', 'BP-IF21M4N', '21.4 MHz first IF of a 25 kHz-channel receiver, 30 kHz bandwidth, stopband above 21.43 MHz (two zeros)',
     50, (21.385 * M, 21.415 * M), 20, 16, 4, [(21.43 * M, 2)]),
    # wide band-pass filters, f2 / f1 from 200 to 3300 (EMC, broadband receivers, instruments, audio)
    ('hid', 'BP-CISPRB', 'CISPR band B conducted-emission preselector, 150 kHz-30 MHz, stopbands below 90 kHz and above 45 MHz (two zeros)',
     50, (150 * K, 30 * M), 20, 12, 4, [(90 * K, 1), (45 * M, 2)]),
    ('hid', 'BP-EMC9K', 'EMC-receiver preselector for CISPR bands A and B, 9 kHz-30 MHz, stopbands below 5 kHz and above 45 MHz (two zeros)',
     50, (9 * K, 30 * M), 20, 12, 4, [(5 * K, 1), (45 * M, 2)]),
    ('hid', 'BP-EMC9K14', 'EMC-receiver preselector, 9 kHz-30 MHz, degree 14, stopbands below 6 kHz and above 40 MHz',
     50, (9 * K, 30 * M), 20, 14, 4, [(6 * K, 1), (40 * M, 1)]),
    ('hid', 'BP-30K30M', 'LF-HF communications receiver front end, 30 kHz-30 MHz, stopbands below 20 kHz and above 45 MHz',
     50, (30 * K, 30 * M), 20, 12, 3, [(20 * K, 1), (45 * M, 1)]),
    ('hid', 'BP-10K10M', 'LF-HF wideband receiver front end (75 ohm), 10 kHz-10 MHz, stopbands below 6 kHz and above 15 MHz (two zeros)',
     75, (10 * K, 10 * M), 20, 16, 4, [(6 * K, 1), (15 * M, 2)]),
    ('hid', 'BP-EDDY', 'eddy-current instrument filter (75 ohm), 10 kHz-10 MHz, stopbands below 5 kHz (two zeros) and above 15 MHz',
     75, (10 * K, 10 * M), 20, 12, 2, [(5 * K, 2), (15 * M, 1)]),
    ('hid', 'BP-3K3M', 'VLF-MF wideband receiver front end, 3 kHz-3 MHz, stopbands below 2 kHz and above 4.5 MHz (two zeros)',
     50, (3 * K, 3 * M), 20, 14, 4, [(2 * K, 1), (4.5 * M, 2)]),
    ('hid', 'BP-1K1M', 'wideband instrumentation input filter (600 ohm), 1 kHz-1 MHz, stopbands below 600 Hz (mains harmonics) and above 1.5 MHz (two zeros)',
     600, (1 * K, 1 * M), 26, 12, 4, [(600, 1), (1.5 * M, 2)]),
    ('hid', 'BP-SONAR', 'sonar receiver input filter (600 ohm), 1 kHz-1 MHz, stopbands below 500 Hz and above 1.4 MHz',
     600, (1 * K, 1 * M), 26, 10, 3, [(500, 1), (1.4 * M, 1)]),
    ('hid', 'BP-AUDIO', 'audio band-limiting filter at a 48 kHz converter input (600 ohm), 20 Hz-20 kHz, stopbands below 8 Hz (infrasonic rumble) and above 28 kHz (alias band, two zeros)',
     600, (20, 20 * K), 26, 12, 4, [(8, 1), (28 * K, 2)]),
    # low-pass filters of degree 10 to 16 (one without finite zeros, one with a single zero)
    ('hid', 'LP-6M', '50-54 MHz amateur transmitter: all-pole (Chebyshev) harmonic low-pass, no finite transmission zeros',
     50, (0, 54 * M), 20, 10, 0, []),
    ('hid', 'LP-AUDIO96', 'anti-alias filter of a 96 kHz audio converter (600 ohm): 40 kHz passband, stopband from 56 kHz (first frequency aliased into the passband)',
     600, (0, 40 * K), 26, 11, 0, [(56 * K, 2)]),
    ('hid', 'LP-ADC50', 'anti-alias low-pass of a 50 MSPS ADC: 20 MHz passband, stopband from 30 MHz (one zero)',
     50, (0, 20 * M), 20, 12, 0, [(30 * M, 1)]),
    ('hid', 'LP-HF30', 'HF direct-sampling receiver anti-alias low-pass: 30 MHz passband, stopband from 46.8 MHz (aliased into the passband at 76.8 MSPS)',
     50, (0, 30 * M), 20, 13, 0, [(46.8 * M, 3)]),
    ('hid', 'LP-DOCSIS85', 'DOCSIS 3.0 mid-split diplexer, 5-85 MHz upstream low-pass arm (75 ohm), stopband from 108 MHz (start of the downstream band)',
     75, (0, 85 * M), 20, 14, 0, [(108 * M, 3)]),
    ('hid', 'LP-ADC122', 'anti-alias low-pass of a 122.88 MSPS ADC: 50 MHz passband, stopband from 72.88 MHz (aliased into the passband)',
     50, (0, 50 * M), 20, 15, 0, [(72.88 * M, 3)]),
    ('hid', 'LP-DAC2G5', 'reconstruction filter of a 2.5 GSPS DAC: 1 GHz passband, stopband from 1.5 GHz (first image)',
     50, (0, 1000 * M), 20, 16, 0, [(1500 * M, 2)]),
    # band-pass filters with a single finite zero (degree 6, 8 and 10; one with only a lower stopband)
    ('hid', 'BP-49M', '49 m broadcast-band preselector, 5.9-6.2 MHz, stopband above 7.2 MHz (41 m band, one zero)',
     50, (5.9 * M, 6.2 * M), 20, 6, 3, [(7.2 * M, 1)]),
    ('hid', 'BP-ADSB', '1090 MHz ADS-B receiver, 1085-1095 MHz, stopband below 1030 MHz (SSR interrogations, one zero)',
     50, (1085 * M, 1095 * M), 20, 8, 2, [(1030 * M, 1)]),
    ('hid', 'BP-GPSL1', 'GPS L1 receiver preselector, 1565-1586 MHz, stopband above 1610 MHz (Globalstar and Iridium uplinks, one zero)',
     50, (1565 * M, 1586 * M), 20, 10, 3, [(1610 * M, 1)]),
    # the narrowest band: a 500 Hz CW filter at 455 kHz (0.11 %)
    ('hid', 'BP-CW455N', '455 kHz CW IF filter, 500 Hz bandwidth, stopbands below 454.25 kHz and above 455.75 kHz (two zeros)',
     50, (454.75 * K, 455.25 * K), 20, 14, 4, [(454.25 * K, 1), (455.75 * K, 2)]),
    # full Cauer filters: as many finite zeros as the degree allows (version 4.2.0)
    ('hid', 'LP-DAC30M', 'reconstruction filter of a 30.72 MSPS DAC: 10 MHz passband, stopband from 20.72 MHz (first image), six zeros',
     50, (0, 10 * M), 20, 13, 0, [(20.72 * M, 6)]),
    ('hid', 'LP-SDR2M4', 'anti-alias low-pass of a 2.4 MSPS SDR receiver: 1 MHz passband, stopband from 1.4 MHz (aliased into the passband), seven zeros',
     50, (0, 1 * M), 20, 15, 0, [(1.4 * M, 7)]),
    ('hid', 'LP-ADC250', 'anti-alias low-pass of a 250 MSPS ADC: 100 MHz passband, stopband from 150 MHz (aliased into the passband), seven zeros',
     50, (0, 100 * M), 20, 16, 0, [(150 * M, 7)]),
    ('hid', 'BP-NBFM10M7', '10.7 MHz narrow-band FM IF filter, 21.4 kHz bandwidth, stopbands below 10.67737 MHz and above 10.72207 MHz (adjacent channels; four and three zeros)',
     50, (10.68931 * M, 10.71071 * M), 14, 16, 1, [(10.67737 * M, 4), (10.72207 * M, 3)]),
    ('hid', 'BP-ROOF21M4', '21.4 MHz first-IF roofing filter, 42.8 kHz bandwidth, stopbands below 21.35546 MHz (one zero) and above 21.44532 MHz (six zeros)',
     50, (21.37861 * M, 21.42141 * M), 20, 16, 1, [(21.35546 * M, 1), (21.44532 * M, 6)]),
    ('hid', 'BP-BANDIII', 'TV band III preselector for the lower channels, 174-192.29 MHz, stopbands below 170.84 MHz (one zero) and above 194.62 MHz (six zeros)',
     50, (174 * M, 192.2917 * M), 20, 16, 1, [(170.8412 * M, 1), (194.6164 * M, 6)]),
    ('hid', 'BP-WLAN24', '2.4 GHz WLAN front-end filter, 2400-2483.5 MHz, stopbands below 2300 MHz (LTE band 40, two zeros) and above 2570 MHz (LTE band 38, five zeros)',
     50, (2400 * M, 2483.5 * M), 20, 16, 1, [(2300 * M, 2), (2570 * M, 5)]),
    ('hid', 'BP-DAB', 'DAB band III receiver preselector (75 ohm), 174-240 MHz, stopbands below 150 MHz (two zeros) and above 270 MHz (four zeros)',
     75, (174 * M, 240 * M), 20, 16, 3, [(150 * M, 2), (270 * M, 4)]),
    ('hid', 'BP-UHFTV', 'UHF TV tuner preselector (75 ohm), 470-862 MHz, stopbands below 400 MHz and above 960 MHz (aeronautical DME band; three zeros each)',
     75, (470 * M, 862 * M), 20, 16, 3, [(400 * M, 3), (960 * M, 3)]),
    ('hid', 'BP-VOICE', 'telephone voice-band filter (600 ohm), 300-3400 Hz, stopbands below 100 Hz (mains hum, two zeros) and above 5 kHz (five zeros)',
     600, (300, 3400), 20, 16, 1, [(100, 2), (5000, 5)]),
    ('hid', 'BP-HF1M6', '1.6-30 MHz HF receiver preselector, stopbands below 530 kHz (one zero) and above 54 MHz (VHF TV, six zeros)',
     50, (1.6 * M, 30 * M), 20, 16, 1, [(0.53 * M, 1), (54 * M, 6)]),
    ('hid', 'BP-HF3M30M', '3-30 MHz HF receiver front end, stopbands below 1.541635 MHz (medium-wave broadcast, one zero) and above 40.58231 MHz (six zeros)',
     50, (3 * M, 30 * M), 20, 16, 1, [(1.541635 * M, 1), (40.58231 * M, 6)]),
    ('hid', 'BP-LFMF', '100 kHz-3 MHz LF/MF receiver front end, stopbands below 40.8874 kHz (VLF transmitters, four zeros) and above 7.032252 MHz (two zeros)',
     50, (100 * K, 3 * M), 26, 16, 3, [(40.8874 * K, 4), (7.032252 * M, 2)]),
]


def num(x):
    x = float('%.12g' % x)
    return int(x) if x == int(x) else x


def main():
    ex, hid = sys.argv[1], sys.argv[2]
    for d in (ex, hid):
        os.makedirs(d, exist_ok=True)
    for grp, sid, desc, rs, pb, rl, n, a, sb in SPECS:
        spec = {'id': sid, 'description': desc, 'source_resistance_ohm': rs,
                'passband_hz': [num(pb[0]), num(pb[1])], 'return_loss_db': rl, 'degree': n,
                'transmission_zeros': {'dc': a, 'stopbands': [{'edge_hz': num(e), 'zeros': k} for e, k in sb]}}
        json.dump(spec, open(os.path.join(ex if grp == 'ex' else hid, sid + '.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
