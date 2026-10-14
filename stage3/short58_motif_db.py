#!/usr/bin/env python3
"""Exploratory: how often does the short-5.8S motif occur in a public ITS database? (independent test of the
motif of short58_motif.py, which was written from the culture copies)

For every sequence of a FASTA with UNITE-style headers (taxonomy as k__...;g__...;s__..., species hypothesis as
SH<digits>.<digits>FU), the motif YTYGGRGARTATG is searched on the given strand. A hit counts as in the 5.8S when the
ITS2 primer site fITS7 (GTGARTCATCGAATCTTTG, up to two mismatches) starts 30-45 nt before it; in the culture
copies of Rhizophagus the offset is 37 nt. A sequence covers the site when fITS7 is found and at least 60 nt
follow its start. fITS7 is searched only in sequences of Rhizophagus, Rhizoglomus and Oehlia, and in sequences
carrying the motif.

  python3 stage3/short58_motif_db.py --fasta UNITE_public_19.02.2025.fasta.gz --out short58_motif_unite

Writes genus.tsv (per genus: sequences, covering the site (target genera only), motif in the 5.8S, motif
elsewhere), species.tsv and
sh.tsv (Rhizophagus and Rhizoglomus), hits.tsv (every sequence with the motif) and manifest.txt. No sequences are
written.
"""
import argparse
import collections
import datetime
import gzip
import hashlib
import os
import re

MOTIF = re.compile('[CT]T[CT]GG[AG]GA[AG]TATG')
FITS7 = 'GTGARTCATCGAATCTTTG'
TARGET = {'Rhizophagus', 'Rhizoglomus', 'Oehlia'}
OFFSET = (30, 45)
SH = re.compile(r'SH\d+\.\d+FU')
TAXG = re.compile(r'g__([^;|]*)')
TAXS = re.compile(r's__([^;|]*)')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fasta(path):
    op = gzip.open if path.endswith('.gz') else open
    head, buf = None, []
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if head is not None:
                    yield head, ''.join(buf).upper()
                head, buf = line[1:].rstrip('\n'), []
            else:
                buf.append(line.strip())
        if head is not None:
            yield head, ''.join(buf).upper()


def fits7_sites(seq, maxmm=2):
    """start positions (0-based) of fITS7 with <= maxmm mismatches (R = A/G)"""
    out, n = [], len(FITS7)
    for i in range(len(seq) - n + 1):
        mm = 0
        for j in range(n):
            p, c = FITS7[j], seq[i + j]
            if not (c == p or (p == 'R' and c in 'AG')):
                mm += 1
                if mm > maxmm:
                    break
        if mm <= maxmm:
            out.append(i)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    gen = collections.defaultdict(lambda: [0, 0, 0, 0])          # seqs, covering, motif in 5.8S, motif elsewhere
    spc = collections.defaultdict(lambda: [0, 0, 0])               # Rhizophagus species: seqs, covering, short
    shs = collections.defaultdict(lambda: [0, 0, 0, collections.Counter()])
    hits, n = [], 0
    for head, seq in fasta(a.fasta):
        n += 1
        g = (TAXG.search(head).group(1) if TAXG.search(head) else '') or 'NA'
        s = (TAXS.search(head).group(1) if TAXS.search(head) else '') or 'NA'
        sh = SH.search(head).group(0) if SH.search(head) else 'NA'
        acc = head.split('|')[0].split()[0]
        ms = [m.start() for m in MOTIF.finditer(seq)]
        gen[g][0] += 1
        if g not in TARGET and not ms:
            continue
        sites = fits7_sites(seq)
        cover = any(len(seq) - i >= 60 for i in sites)
        in58 = [m for m in ms if any(OFFSET[0] <= m - i <= OFFSET[1] for i in sites)]
        gen[g][1] += cover
        gen[g][2] += bool(in58)
        gen[g][3] += bool(ms) and not in58
        if g in ('Rhizophagus', 'Rhizoglomus'):
            spc[s][0] += 1
            spc[s][1] += cover
            spc[s][2] += bool(in58)
            shs[sh][0] += 1
            shs[sh][1] += cover
            shs[sh][2] += bool(in58)
            shs[sh][3][s] += 1
        if ms:
            hits.append((acc, g, s, sh, len(seq), ','.join(str(m + 1) for m in ms),
                         ','.join(str(i + 1) for i in sites), int(bool(in58))))
    with open(os.path.join(a.out, 'genus.tsv'), 'w') as fh:
        fh.write('genus\tsequences\tcovering_site\tmotif_in_5.8S\tmotif_elsewhere\n')
        for g, v in sorted(gen.items(), key=lambda x: (-x[1][2], -x[1][3], x[0])):
            if v[1] or v[2] or v[3]:
                cov = str(v[1]) if g in TARGET else ''          # site coverage is only checked in the target genera
                fh.write(f'{g}\t{v[0]}\t{cov}\t{v[2]}\t{v[3]}\n')
    with open(os.path.join(a.out, 'species.tsv'), 'w') as fh:
        fh.write('species\tsequences\tcovering_site\tshort_type\n')
        for s, v in sorted(spc.items(), key=lambda x: -x[1][0]):
            fh.write(f'{s}\t' + '\t'.join(map(str, v)) + '\n')
    with open(os.path.join(a.out, 'sh.tsv'), 'w') as fh:
        fh.write('sh\tsequences\tcovering_site\tshort_type\tshort_pct_of_covering\tspecies_labels\n')
        for sh, v in sorted(shs.items(), key=lambda x: -x[1][0]):
            pct = f'{100 * v[2] / v[1]:.1f}' if v[1] else ''
            labs = '; '.join(f'{k} ({c})' for k, c in v[3].most_common(4))
            fh.write(f'{sh}\t{v[0]}\t{v[1]}\t{v[2]}\t{pct}\t{labs}\n')
    with open(os.path.join(a.out, 'hits.tsv'), 'w') as fh:
        fh.write('accession\tgenus\tspecies\tsh\tlength\tmotif_starts\tfits7_starts\tin_5.8S\n')
        for h in hits:
            fh.write('\t'.join(map(str, h)) + '\n')
    with open(os.path.join(a.out, 'manifest.txt'), 'w') as fh:
        fh.write(f'date\t{datetime.datetime.now().isoformat(timespec="seconds")}\n')
        fh.write(f'script\tshort58_motif_db.py\t{sha256(os.path.abspath(__file__))}\n')
        fh.write(f'input\t{a.fasta}\t{sha256(a.fasta)}\nsequences\t{n}\n')
    rz = [gen[g] for g in ('Rhizophagus', 'Rhizoglomus')]
    print(f'{n} sequences; Rhizophagus + Rhizoglomus: {sum(v[0] for v in rz)} sequences, '
          f'{sum(v[1] for v in rz)} covering the site, {sum(v[2] for v in rz)} short type; '
          f'motif in the 5.8S of other genera: {sum(v[2] for g, v in gen.items() if g not in ("Rhizophagus", "Rhizoglomus"))}')
    print(f'see {a.out}/genus.tsv, species.tsv, sh.tsv')


if __name__ == '__main__':
    main()
