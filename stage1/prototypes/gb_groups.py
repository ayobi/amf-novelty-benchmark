#!/usr/bin/env python3
"""Group the Stefani et al. (2025) rDNA copies by library and by culture from GenBank source qualifiers.

  library  /isolate minus the per-copy suffix (_s<N> or _s<N>_<M>), i.e. one barcoded DNA extract
  culture  /strain, else /culture_collection, else /specimen_voucher, else the library label

Group labels are prefixed with the organism so identical labels in different species never merge.

Usage:
  python3 gb_groups.py stefani_rdna.gb -o stefani_groups.tsv
"""
import argparse
import collections
import csv
import re
import statistics

SUFFIX = re.compile(r'[_-]s\d+(?:_\d+)?$')
QUALS = ('isolate', 'strain', 'culture_collection', 'specimen_voucher', 'note')


def parse_gb(path):
    """Yield one dict per record: accession plus its source-feature qualifiers."""
    acc, quals, cur, in_src = None, None, None, False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if line.startswith('LOCUS'):
                acc, quals, cur, in_src = None, collections.defaultdict(list), None, False
            elif line.startswith('VERSION'):
                acc = line.split()[1]
            elif line.startswith('//'):
                if acc:
                    rec = {k: '; '.join(v.strip().strip('"') for v in vals) for k, vals in quals.items()}
                    rec['accession'] = acc
                    yield rec
                acc, in_src = None, False
            elif line.startswith('     source '):
                in_src = True
            elif in_src:
                if line.startswith(' ' * 21 + '/'):
                    cur, _, val = line[22:].partition('=')
                    quals[cur].append(val)
                elif line.startswith(' ' * 21) and cur:
                    # NCBI wraps at a space when it can and hard-breaks long tokens otherwise
                    sep = ' ' if ' ' in quals[cur][-1] else ''
                    quals[cur][-1] += sep + line.strip()
                else:
                    in_src = False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gb', help='GenBank flat file (efetch -format gb)')
    ap.add_argument('-o', '--out', required=True, help='output TSV')
    ap.add_argument('--min-copies', type=int, default=3,
                    help='smallest group counted as usable for divergence (default 3)')
    a = ap.parse_args()

    rows, shared = [], collections.Counter()
    for rec in parse_gb(a.gb):
        org = rec.get('organism', '')
        iso = rec.get('isolate', '')
        lib = SUFFIX.sub('', iso)
        if iso and lib == iso:
            shared[iso] += 1  # no copy suffix: the label names the culture and copies share it
        cul, src = next(((rec[k], k) for k in ('strain', 'culture_collection', 'specimen_voucher') if rec.get(k)),
                        (lib, 'isolate' if lib else ''))
        lib = lib or cul
        rows.append({
            'accession': rec['accession'],
            'organism': org,
            'genus': org.split()[0].strip('[]') if org else '',
            **{k: rec.get(k, '') for k in QUALS},
            'library': f'{org} | {lib}' if lib else '',
            'culture': f'{org} | {cul}' if cul else '',
            'culture_from': src,
        })
    if not rows:
        raise SystemExit(f'no records parsed from {a.gb}')

    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

    print(f'{len(rows)} records -> {a.out}')
    print('qualifier present in:', {k: sum(bool(r[k]) for r in rows) for k in QUALS})
    print(f'per-copy isolate labels: {len(rows) - sum(shared.values())} records; '
          f'shared culture labels: {len(shared)} labels over {sum(shared.values())} records '
          f'(max {max(shared.values(), default=0)} per label)')
    ungrouped = collections.Counter(r['organism'] for r in rows if not r['culture'])
    print(f'records with no grouping label: {sum(ungrouped.values())}', ungrouped.most_common(5))

    for level in ('library', 'culture'):
        size = collections.Counter(r[level] for r in rows if r[level])
        if not size:
            print(f'\n{level}: no groups')
            continue
        usable = [c for c in size.values() if c >= a.min_copies]
        print(f'\n{level}: {len(size)} groups, {len(usable)} with >= {a.min_copies} copies ({sum(usable)} copies); '
              f'group size median {statistics.median(size.values())}, max {max(size.values())}')

    print('\nRhizophagus irregularis (copies, culture, library):')
    ri = collections.Counter((r['culture'], r['culture_from'], r['library']) for r in rows
                             if r['organism'] == 'Rhizophagus irregularis')
    for (cul, src, lib), c in sorted(ri.items()):
        print(f'  {c:3d}  {cul.split(" | ")[-1]} ({src})  {lib.split(" | ")[-1]}')

    pooled = collections.Counter((r['organism'], r['culture_from'] or '-', r['culture'].split(' | ')[-1] or '-')
                                 for r in rows if not r['isolate'])
    print(f'\nrecords without /isolate (previously pooled by species): {sum(pooled.values())}')
    for (org, src, lab), c in pooled.most_common(30):
        print(f'  {c:3d}  {org:34s} {src:18s} {lab}')
    notes = sorted({r['note'] for r in rows if not r['isolate'] and r['note']})
    print(f'  distinct /note values among them: {len(notes)}', notes[:3])


if __name__ == '__main__':
    main()
