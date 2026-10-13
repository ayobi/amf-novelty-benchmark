#!/usr/bin/env python3
"""How do identity thresholds behave on whole long reads? L. Tedersoo's question II.2-3 for the AMF paper (13 October
2026): "how would the entire 4.5-kb fragment work here?". Exploratory: it is not part of the Phase 0 protocol and
scores nothing that protocol decides.

Reads: the Glomeromycota reads of the EUKARYOME long-read set whose species field gives a named species (protocol
section 2 rules), with all three windows cut. Each read is rebuilt in rRNA orientation as SSU + ITS + LSU window
(cut_windows.py splits a read into exactly these three parts), so reads stored on the other strand are turned round.

Regions per read:
  its      the full ITS window
  kruger   from the SSUmAf1 site to the LSUmBr site (Kruger et al. 2009), as in the culture analysis of the paper
  read     the whole read

For every read and region: the distance to the nearest other read of its species, of another named species of
its genus, and of another genus (edlib infix edit distance as % of the shorter sequence, as in the paper). From
these: the within-species spread (99th percentile of the nearest-conspecific distance), and for each cut-off t the
share of reads split from their species (no other read of it within t) and joined to another species of the genus
(one within t), overall and outside the paper's four complexes.

These are environmental reads with EUKARYOME's classification as the labels, not cultures: there is no genome-level
view, and a label error looks like a species overlap.

  pip install edlib     # if the environment lacks it
  python3 scripts/longread_thresholds.py --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta \\
      --windows windows --out amf_longread --procs 32
"""
import argparse
import collections
import gzip
import multiprocessing as mp
import os
import re
import statistics
import sys

import edlib

TAX = re.compile(r'(?:^|[;|,])\s*([kpcofgs])__([^;|,]*)')
EPITHET = re.compile(r'^[a-z][a-z-]+$')
BINOMIAL = re.compile(r'^[A-Z][a-z]+[_ ][a-z][a-z-]+$')
PLACEHOLDER = re.compile(r'(?i)(^|[_ .])(sp|spp|cf|aff|nov|unidentified|uncultured|unclassified|incertae|'
                         r'sedis|environmental|metagenome)([_ .]|$)|\d|\.')
IUPAC = {'Y': 'CT', 'R': 'AG', 'D': 'AGT', 'N': 'ACGT', 'W': 'AT', 'S': 'CG', 'K': 'GT', 'M': 'AC', 'H': 'ACT',
         'B': 'CGT', 'V': 'ACG'}
EQ = [(k, c) for k, v in IUPAC.items() for c in v]
SSUMAF1 = 'TGGGTAATCTTTTGAAACTTYA'
LSUMBR_RC = 'TCTARCATRTRTGCGAGTGTT'          # reverse complement of the LSUmBr consensus AACACTCGCAYAYATGYTAGA
WINDOW_FILES = {'ssu': 'win.ssu.fasta', 'its': 'win.full.fasta', 'lsu': 'win.lsu.fasta'}
REGIONS = ['its', 'kruger', 'read']
CLASSES = ['same species', 'same genus', 'other genus']
# the paper's four complexes (paper/analysis/v4.complexes.tsv); EUKARYOME may spell vesiculiferum as vesiculifer
COMPLEXES = [{'Glomus_chinense', 'Glomus_macrocarpum', 'Glomus_rugosae'},
             {'Diversispora_aestuarii', 'Diversispora_varaderana'},
             {'Diversispora_arenaria', 'Diversispora_jakucsiae'},
             {'Rhizophagus_irregularis', 'Rhizoglomus_vesiculiferum', 'Rhizophagus_vesiculifer',
              'Rhizoglomus_vesiculifer', 'Rhizophagus_vesiculiferum', 'Rhizoglomus_irregulare',
              'Rhizoglomus_irregularis'}]
SWEEP = [round(0.05 * i, 2) for i in range(0, 241)]
GENUS_SYN = {'Rhizoglomus': 'Rhizophagus'}     # as in the paper: one genus for the classes
# EUKARYOME v2.1 spells one species two ways; the two labels are merged (reported in the output)
SPECIES_SYN = {'Rhizoglomus_irregularis': 'Rhizoglomus_irregulare'}


def read_id(header):
    return re.split(r'[\s;]', header, maxsplit=1)[0]


def fasta(path):
    head, seq = None, []
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as fh:
        for line in fh:
            if line.startswith('>'):
                if head is not None:
                    yield head, ''.join(seq)
                head, seq = line[1:].rstrip('\n'), []
            else:
                seq.append(line.strip())
    if head is not None:
        yield head, ''.join(seq)


def names(header):
    t = {}
    for r, v in TAX.findall(header):
        t.setdefault(r, v.strip())
    g, s = t.get('g', ''), t.get('s', '')
    named_g = bool(g) and g[0].isupper() and not PLACEHOLDER.search(g)
    if s and EPITHET.match(s) and named_g:
        s = f'{g}_{s}'
    named_s = bool(s) and bool(BINOMIAL.match(s)) and not PLACEHOLDER.search(s)
    return t.get('p', ''), (g if named_g else None), (s if named_s else None)


def site(primer, s):
    r = edlib.align(primer, s, mode='HW', task='locations', additionalEqualities=EQ)
    return r['locations'][0] if r['editDistance'] <= 3 else None


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode='HW', task='distance')['editDistance'] / len(q)


G = {}


def init(reads):
    G['reads'] = reads


def nearest(i):
    """Nearest read of each class, per region, for read i."""
    reads = G['reads']
    ri = reads[i]
    best = {}
    for j, rj in enumerate(reads):
        if j == i:
            continue
        rel = ('same species' if rj['species'] == ri['species'] else
               'same genus' if rj['genus'] == ri['genus'] else 'other genus')
        for m in REGIONS:
            a, b = ri[m], rj[m]
            if a is None or b is None:
                continue
            k = (m, rel)
            d = dist(a, b)
            if k not in best or d < best[k][0]:
                best[k] = (d, rj['species'])
    return i, best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--windows', required=True)
    ap.add_argument('--out', default='amf_longread')
    ap.add_argument('--procs', type=int, default=os.cpu_count())
    ap.add_argument('--phylum', default='Glomeromycota')
    ap.add_argument('--min-reads', type=int, default=2, help='reads a species needs to enter (default 2)')
    ap.add_argument('--kruger-to-end', action='store_true', help='test data only: end the Kruger region at the '
                    'read end when LSUmBr is absent')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    lab = {}
    for h, _ in fasta(a.fasta):
        p, g, s = names(h)
        if p == a.phylum and g and s:
            lab[read_id(h)] = (GENUS_SYN.get(g, g), SPECIES_SYN.get(s, s))
    win = {m: {} for m in WINDOW_FILES}
    for m, f in WINDOW_FILES.items():
        for h, s in fasta(os.path.join(a.windows, f)):
            r = read_id(h)
            if r in lab:
                win[m][r] = s.upper()
    n_sp = collections.Counter(s for g, s in lab.values())
    reads, miss = [], collections.Counter()
    for r, (g, s) in sorted(lab.items()):
        if n_sp[s] < a.min_reads:
            miss['species with too few reads'] += 1
            continue
        if not all(r in win[m] for m in WINDOW_FILES):
            miss['windows not cut'] += 1
            continue
        full = win['ssu'][r] + win['its'][r] + win['lsu'][r]
        a1 = site(SSUMAF1, full)
        b1 = site(LSUMBR_RC, full)
        kr = None
        if a1 and b1 and b1[0] > a1[0]:
            kr = full[a1[0]:b1[0]]
        elif a1 and a.kruger_to_end:
            kr = full[a1[0]:]
        else:
            miss['no Kruger region (primer site missing)'] += 1
        reads.append({'id': r, 'genus': g, 'species': s, 'its': win['its'][r], 'kruger': kr, 'read': full})
    sp_n = collections.Counter(r['species'] for r in reads)
    print(f'{len(reads)} reads of {len(sp_n)} named {a.phylum} species with >= {a.min_reads} reads; left out: '
          f'{dict(miss)}; species labels merged: {SPECIES_SYN}', flush=True)
    lim = []
    for m in REGIONS:
        L = [len(r[m]) for r in reads if r[m]]
        print(f'  {m}: {len(L)} reads, length median {statistics.median(L):.0f} ({min(L)}-{max(L)})', flush=True)
        lim.append({'region': m, 'reads': len(L), 'species': len({r['species'] for r in reads if r[m]}),
                    'median_length': statistics.median(L), 'min_length': min(L), 'max_length': max(L)})
    with open(os.path.join(a.out, 'longread.limits.tsv'), 'w') as fh:
        fh.write('\t'.join(lim[0]) + '\tleft_out\tlabels_merged\n')
        for x in lim:
            fh.write('\t'.join(str(v) for v in x.values()) + f'\t{dict(miss)}\t{SPECIES_SYN}\n')

    with mp.Pool(a.procs, initializer=init, initargs=(reads,)) as pool:
        res = dict(pool.imap_unordered(nearest, range(len(reads)), chunksize=4))
    cx = {s: i for i, grp in enumerate(COMPLEXES) for s in grp}
    rows = []
    for i, r in enumerate(reads):
        for (m, rel), (d, other) in res[i].items():
            rows.append({'read': r['id'], 'species': r['species'], 'region': m, 'class': rel,
                         'nearest_pct': round(d, 4), 'nearest_species': other})
    with open(os.path.join(a.out, 'longread.nearest.tsv'), 'w') as fh:
        fh.write('read\tspecies\tregion\tclass\tnearest_pct\tnearest_species\n')
        for x in rows:
            fh.write('\t'.join(str(x[k]) for k in ('read', 'species', 'region', 'class', 'nearest_pct',
                                                     'nearest_species')) + '\n')

    def pct(v, f):
        return 100 * sum(1 for x in v if f(x)) / len(v) if v else float('nan')

    summ, sweep = [], []
    for m in REGIONS:
        own = [res[i][(m, 'same species')][0] for i in range(len(reads)) if (m, 'same species') in res[i]]
        cong = [res[i][(m, 'same genus')][0] for i in range(len(reads)) if (m, 'same genus') in res[i]]
        # outside the complexes: the nearest congener that is not a complex partner of the read's species
        cong_x = []
        for i, r in enumerate(reads):
            if (m, 'same genus') not in res[i]:
                continue
            d, other = res[i][(m, 'same genus')]
            if r['species'] in cx and cx.get(other) == cx[r['species']]:
                # recompute the nearest congener outside the complex (rare: only complex species)
                best = None
                for rj in reads:
                    if rj['genus'] == r['genus'] and rj['species'] != r['species'] and cx.get(rj['species']) != \
                            cx[r['species']] and r[m] and rj[m]:
                        dd = dist(r[m], rj[m])
                        best = dd if best is None or dd < best else best
                if best is not None:
                    cong_x.append(best)
            else:
                cong_x.append(d)
        oth = [res[i][(m, 'other genus')][0] for i in range(len(reads)) if (m, 'other genus') in res[i]]
        spread = statistics.quantiles(own, n=100, method='exclusive')[98] if len(own) > 1 else float('nan')
        curve = []
        for t in SWEEP:
            row = {'region': m, 'threshold_pct': t, 'identity_pct': round(100 - t, 2),
                   'split_from_species': round(pct(own, lambda x: x > t), 2),
                   'joined_to_congener': round(pct(cong, lambda x: x <= t), 2),
                   'joined_to_congener_excl_complex': round(pct(cong_x, lambda x: x <= t), 2),
                   'joined_to_other_genus': round(pct(oth, lambda x: x <= t), 2)}
            sweep.append(row)
            curve.append(row)
        best = min(curve, key=lambda x: (max(x['split_from_species'], x['joined_to_congener']), x['threshold_pct']))
        best_x = min(curve, key=lambda x: (max(x['split_from_species'], x['joined_to_congener_excl_complex']),
                                           x['threshold_pct']))
        for name, t in [('98.5% identity', 1.5), ('97% identity', 3.0), ('within-species spread', spread),
                        ('best single cut-off', best['threshold_pct']),
                        ('best single cut-off, complexes aside', best_x['threshold_pct'])]:
            summ.append({'region': m, 'threshold': name, 'identity_pct': round(100 - t, 2),
                         'reads_with_conspecific': len(own), 'reads_with_congener': len(cong),
                         'split_from_species_pct': round(pct(own, lambda x: x > t), 1),
                         'joined_to_congener_pct': round(pct(cong, lambda x: x <= t), 1),
                         'joined_to_congener_excl_complex_pct': round(pct(cong_x, lambda x: x <= t), 1),
                         'joined_to_other_genus_pct': round(pct(oth, lambda x: x <= t), 1)})
    for name, data in (('longread.sweep.tsv', sweep), ('longread.summary.tsv', summ)):
        with open(os.path.join(a.out, name), 'w') as fh:
            fh.write('\t'.join(data[0]) + '\n')
            for x in data:
                fh.write('\t'.join(str(x[k]) for k in data[0]) + '\n')
    # which species pairs make the joins, at the database thresholds and at each region's best cut-off
    with open(os.path.join(a.out, 'longread.pairs.tsv'), 'w') as fh:
        fh.write('region\tcutoff\tidentity_pct\tspecies\tnearest_congener\treads_joined\n')
        for m in REGIONS:
            best = next(x for x in summ if x['region'] == m and x['threshold'] == 'best single cut-off')
            for name, t in (('97% identity', 3.0), ('98.5% identity', 1.5),
                            ('best single cut-off', round(100 - best['identity_pct'], 2))):
                pc = collections.Counter()
                for i, r in enumerate(reads):
                    if (m, 'same genus') in res[i] and res[i][(m, 'same genus')][0] <= t:
                        pc[tuple(sorted((r['species'], res[i][(m, 'same genus')][1])))] += 1
                for (s1, s2), n in pc.most_common():
                    fh.write(f'{m}\t{name}\t{100 - t:.2f}\t{s1}\t{s2}\t{n}\n')
    with open(os.path.join(a.out, 'longread.species.tsv'), 'w') as fh:
        fh.write('species\treads\n')
        for s, n in sp_n.most_common():
            fh.write(f'{s}\t{n}\n')
    for x in summ:
        print(f"{x['region']:7s} {x['threshold']:38s} {x['identity_pct']:6.2f}%  split {x['split_from_species_pct']:5.1f}"
              f"  congener {x['joined_to_congener_pct']:5.1f} (no complexes {x['joined_to_congener_excl_complex_pct']:5.1f})"
              f"  other genus {x['joined_to_other_genus_pct']:4.1f}")
    print(f'reads with a conspecific / a congener: {summ[0]["reads_with_conspecific"]} / {summ[0]["reads_with_congener"]}')


if __name__ == '__main__':
    main()
