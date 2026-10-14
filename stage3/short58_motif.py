#!/usr/bin/env python3
"""Exploratory: can a short sequence motif flag the short-5.8S paralog? (L. Tedersoo's "microbarcode", 14 Oct 2026)

In Rhizophagus and [Rhizoglomus] the short type differs from the main type in one hairpin near the 3' end of the
5.8S (main ...CACTCTCTGG(T/C)AACCCGGGGAGTATG..., short ...CACTCTCTT(T/C)GGAGAGTATG...). The motif tested here,
YTYGGRGARTATG (IUPAC), covers the short-type sequence of every short copy of those genera. It was written from
these copies and is tested on the same copies, so the counts show that such a flag can work, not its error rate on
new data.

For every retained copy of the culture reference: whether its 5.8S carries the motif, where, and whether the motif
occurs anywhere else in the whole copy. Short copies are those of copy_paralog_effects.py: the short-5.8S class
and every copy whose 5.8S is identical to one of it. The ITS2 primers fITS7 and ITS3 are located in the 5.8S
with at most two mismatches, to show which ITS2 amplicons carry the motif.

  python3 stage3/short58_motif.py --ref reference/culture_ref_v4 --out benchmarks/copy_checks/short58_motif

Writes motif.copies.tsv, motif.summary.tsv, primers.tsv and README.md.
"""
import argparse
import collections
import csv
import os
import re
import statistics

MOTIF = 'YTYGGRGARTATG'
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': '[AG]', 'Y': '[CT]', 'W': '[AT]', 'S': '[CG]',
         'K': '[GT]', 'M': '[AC]', 'N': '[ACGT]'}
PRIMERS = {'fITS7': 'GTGARTCATCGAATCTTTG', 'ITS3': 'GCATCGATGAAGAACGCAGC'}


def read_fasta(path):
    out, cur = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                cur = line[1:].split()[0]
                out[cur] = []
            elif cur:
                out[cur].append(line.upper())
    return {k: ''.join(v) for k, v in out.items()}


def iupac_ok(a, b):
    return b == 'N' or a == b or (b in IUPAC and len(IUPAC[b]) > 1 and a in IUPAC[b])


def find_primer(seq, primer, maxmm=2):
    """best (start, mismatches) of primer in seq, 1-based, or (None, None)"""
    best = (None, None)
    for i in range(len(seq) - len(primer) + 1):
        mm = sum(not iupac_ok(seq[i + j], primer[j]) for j in range(len(primer)))
        if mm <= maxmm and (best[1] is None or mm < best[1]):
            best = (i + 1, mm)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    ref = {r['accession']: r for r in csv.DictReader(open(os.path.join(a.ref, 'copies.tsv')), delimiter='\t')
           if r['status'] == 'retained'}
    s58 = read_fasta(os.path.join(a.ref, 'fasta', 's58.fasta'))
    full = read_fasta(os.path.join(a.ref, 'fasta', 'full.fasta'))
    cls58 = {s58[x] for x, r in ref.items() if r['div_kind'] == 'short_5.8S'}
    short = {x for x in ref if ref[x]['div_kind'] == 'short_5.8S' or s58.get(x) in cls58}
    pat = re.compile(''.join(IUPAC[c] for c in MOTIF))
    rows, summ = [], collections.defaultdict(lambda: [0, 0])
    for x in sorted(ref):
        r = ref[x]
        cls = 'short 5.8S' if x in short else ('same length' if r['div_kind'] == 'same_length' else 'main type')
        m = pat.search(s58.get(x, ''))
        n_full = len(pat.findall(full.get(x, '')))
        n_else = n_full - (1 if m else 0)
        genus = 'Rhizophagus' if r['genus'] in ('Rhizophagus', 'Rhizoglomus') else r['genus']
        rows.append((x, genus, r['organism'], r['culture_label'], cls, len(s58.get(x, '')), int(bool(m)),
                     m.start() + 1 if m else '', n_else))
        summ[(genus, cls)][0] += 1
        summ[(genus, cls)][1] += bool(m)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, 'motif.copies.tsv'), 'w') as fh:
        fh.write('accession\tgenus\torganism\tculture\tclass\ts58_length\tmotif_in_5.8S\tmotif_start\t'
                 'motif_elsewhere_in_copy\n')
        for row in rows:
            fh.write('\t'.join(map(str, row)) + '\n')
    with open(os.path.join(a.out, 'motif.summary.tsv'), 'w') as fh:
        fh.write('genus\tclass\tcopies\twith_motif\n')
        for (g, c), (n, k) in sorted(summ.items()):
            fh.write(f'{g}\t{c}\t{n}\t{k}\n')
    prow = []
    for cls in ('short 5.8S', 'main type'):
        xs = [row[0] for row in rows if (row[0] in short) == (cls == 'short 5.8S')
              and ref[row[0]]['genus'] in ('Rhizophagus', 'Rhizoglomus')]
        for name, p in PRIMERS.items():
            hits = [find_primer(s58[x], p) for x in xs]
            pos = [h[0] for h in hits if h[0]]
            mms = [h[1] for h in hits if h[0] is not None]
            prow.append((cls, name, len(xs), len(pos), statistics.median(pos) if pos else '',
                         min(pos) if pos else '', max(pos) if pos else '',
                         dict(collections.Counter(mms))))
        ms = [row[7] for row in rows if row[0] in xs and row[7] != '']
        prow.append((cls, 'motif', len(xs), len(ms), statistics.median(ms) if ms else '',
                     min(ms) if ms else '', max(ms) if ms else '', ''))
    with open(os.path.join(a.out, 'primers.tsv'), 'w') as fh:
        fh.write('copies\tsite\tn\tfound\tmedian_start\tmin_start\tmax_start\tmismatches\n')
        for row in prow:
            fh.write('\t'.join(map(str, row)) + '\n')
    for (g, c), (n, k) in sorted(summ.items()):
        if k or c == 'short 5.8S' or g == 'Rhizophagus':
            print(f'{g:16s} {c:12s} {k}/{n}')
    for row in prow:
        print(row)


if __name__ == '__main__':
    main()
