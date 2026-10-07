#!/usr/bin/env python3
"""Do protein-coding genes separate the rDNA species complexes?

Stefani et al. (2025) sequenced three protein-coding genes for 143 cultures of 39 species and deposited
them in GenBank: glomalin PV808915-PV809055, RPB1 PV809056-PV809196, H+-ATPase PV809197-PV809333. They
give indicative species thresholds of 1.0% (glomalin), 1.1% (RPB1) and 1.7% (H+-ATPase).

For every pair of named reference species in v4.pairs.tsv (species_groups.py) that both have records for
a gene: the smallest distance between a record of one and a record of the other, and the largest distance
between two records of each species. A pair counts as separated when that smallest distance exceeds the
gene's threshold, and as cleanly separated when it also exceeds both species' own spread. The complexes
(pairs inseparable on SSU, ITS and LSU) are compared with the ordinary congeneric pairs, which show what
the genes do where rDNA already separates the species.

Distances: edlib infix edit distance as a percentage of the shorter sequence, as everywhere in the
project. Stefani et al. used alignment-based distances, so their thresholds are a guide here, not an
exact cut-off.

Records are joined to species through the curated culture names. When a record's culture label is a
culture of the reference (--cultures), the record is used only if its GenBank organism and the reference
name of that culture agree; a disagreement is left out and listed, because one of the two deposits names
the culture wrongly. Records of cultures outside the reference keep their GenBank name. Names are matched
with Rhizoglomus read as Rhizophagus, and an epithet that differs only in its ending (vesiculiferum,
vesiculiferus) is matched and reported.

The GenBank records are downloaded once (NCBI E-utilities, batches of 100) and cached in --out.

Usage (from the paper directory; internet is needed for the first run only):
  python3 analysis/complex_pcgenes.py --pairs analysis/v4.pairs.tsv --complexes analysis/v4.complexes.tsv \\
      --out ../bench/pcgenes --email you@example.org
Writes {out}/pcgenes.pairs.tsv, {out}/pcgenes.records.tsv and {out}/pcgenes.summary.json; with --paper
also numbers_pcgenes.tex and tables/tableS2_pcgenes.tex in the paper directory.
"""
import argparse
from collections import defaultdict
import csv
import itertools
import json
from pathlib import Path
import re
import statistics
import sys
import time
import urllib.parse
import urllib.request

import edlib

GENES = {  # gene: (first accession, last accession, indicative threshold % from Stefani et al. 2025)
    'glomalin': ('PV808915', 'PV809055', 1.0),
    'RPB1': ('PV809056', 'PV809196', 1.1),
    'H+-ATPase': ('PV809197', 'PV809333', 1.7),
}
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
GENUS_SYN = {'Rhizoglomus': 'Rhizophagus', 'Claroideoglomus': 'Entrophospora'}


def acc_range(first, last):
    p1, n1 = re.match(r'([A-Z]+)(\d+)', first).groups()
    p2, n2 = re.match(r'([A-Z]+)(\d+)', last).groups()
    if p1 != p2:
        raise ValueError('accession prefixes differ')
    return [f'{p1}{i:0{len(n1)}d}' for i in range(int(n1), int(n2) + 1)]


def fetch(ids, email, api_key, batch=100):
    out = []
    for i in range(0, len(ids), batch):
        params = {'db': 'nuccore', 'id': ','.join(ids[i:i + batch]), 'rettype': 'gb', 'retmode': 'text',
                  'tool': 'complex_pcgenes', 'email': email}
        if api_key:
            params['api_key'] = api_key
        data = urllib.parse.urlencode(params).encode()
        for attempt in range(4):
            try:
                with urllib.request.urlopen(urllib.request.Request(EUTILS, data=data), timeout=120) as r:
                    out.append(r.read().decode())
                break
            except Exception as e:  # network hiccups: retry with a pause
                if attempt == 3:
                    raise SystemExit(f'download failed for batch {i // batch + 1}: {e}')
                time.sleep(3 * (attempt + 1))
        time.sleep(0.4 if not api_key else 0.12)
    return ''.join(out)


def parse_gb(text):
    recs = []
    for block in text.split('\n//'):
        if 'ACCESSION' not in block:
            continue
        acc = re.search(r'^ACCESSION\s+(\S+)', block, re.M)
        org = re.search(r'/organism="([^"]+)"', block)
        strain = (re.search(r'/strain="([^"]+)"', block) or re.search(r'/culture_collection="([^"]+)"', block)
                  or re.search(r'/isolate="([^"]+)"', block) or re.search(r'/specimen_voucher="([^"]+)"', block))
        origin = block.split('\nORIGIN', 1)
        seq = re.sub(r'[^acgtnACGTN]', '', re.sub(r'\d', '', origin[1])) if len(origin) == 2 else ''
        if acc and org and seq:
            recs.append({'accession': acc.group(1), 'organism': ' '.join(org.group(1).split()),
                         'culture': ' '.join(strain.group(1).split()) if strain else '', 'seq': seq.upper()})
    return recs


def norm(name):
    t = name.replace('[', '').replace(']', '').split()
    if len(t) < 2:
        return None
    return GENUS_SYN.get(t[0], t[0]), t[1].lower()


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode='HW', task='distance')['editDistance'] / len(q)


def tex_sp(name, full=False):
    """\\spn{G.\\ epithet}; [Rhizoglomus] -> [R.]; 'Genus sp.' -> \\spn{Genus} sp.; full=True keeps the genus."""
    t = name.split()
    if len(t) < 2:
        return f'\\spn{{{name}}}'
    g, ep = t[0], ' '.join(t[1:])
    if ep.rstrip('.') in ('sp', 'spp'):
        return f'\\spn{{{g.strip("[]")}}} sp.'
    if full:
        return f'\\spn{{{name}}}'
    if g.startswith('['):
        return f'\\spn{{[{g[1]}.]\\ {ep}}}'
    return f'\\spn{{{g[0]}.\\ {ep}}}'


def write_paper(paper, rows, summary, records):
    """numbers_pcgenes.tex and tables/tableS2_pcgenes.tex: every protein-coding number the paper cites."""
    lines = ['% GENERATED by analysis/complex_pcgenes.py from the Stefani et al. (2025) protein-coding records; do not edit.']
    add = lambda k, v: lines.append(f'\\newcommand{{\\{k}}}{{{v}}}')
    add('PcRecords', sum(len(v) for v in records.values()))
    add('PcConflicts', len(summary['conflicts']))
    cult = sorted({(c['culture'], c['reference'], c['genbank']) for c in summary['conflicts']})
    add('PcConflictDetail', '; '.join(f"{c.replace('_', chr(92) + '_')} ({tex_sp(ref, True)} in the reference, "
                                      f"{tex_sp(gb, True)} in GenBank)" for c, ref, gb in cult) or 'none')
    near = [c['nearest_reference_pct'] for c in summary['conflicts'] if c.get('nearest_reference_pct') is not None]
    add('PcConflictNearest', f'{min(near):.1f}--{max(near):.1f}\\%' if near else '--')
    base = summary['congeneric_baseline']
    add('PcBaseline', '; '.join(f"{v['separated']} of {v['pairs']} on {g.replace('+', '$^+$')}" for g, v in base.items()))
    shares = [100 * v['separated'] / v['pairs'] for v in base.values()]
    add('PcBaselinePct', f'{min(shares):.0f}--{max(shares):.0f}\\%' if round(min(shares)) != round(max(shares)) else f'{min(shares):.0f}\\%')
    cp = summary['complex_pairs']
    pairs = sorted({(r['species_a'], r['species_b']) for r in cp})
    add('PcComplexPairsN', len(pairs))
    add('PcComplexPairsAll', len(pairs) + len(summary['missing_pairs']))
    above = [r for r in cp if r['separated']]
    add('PcComplexAboveN', len(above))
    add('PcComplexAbove', '; '.join(f"{tex_sp(r['species_a'])} and {tex_sp(r['species_b'])} on {r['gene'].replace('+', '$^+$')} "
                                    f"({r['min_between']:.2f}\\% against {r['threshold']:.1f}\\%)" for r in above) or 'none')
    add('PcComplexClean', sum(r['cleanly_separated'] == 1 for r in cp))
    n_clean = sum(r['cleanly_separated'] == 1 for r in cp)
    add('PcComplexCleanText', 'none' if n_clean == 0 else str(n_clean))
    add('PcComplexComparisons', len(cp))
    add('PcMissing', '; '.join(' / '.join(tex_sp(x.strip()) for x in m.split(' / ')) for m in summary['missing_pairs']) or 'none')
    (paper / 'numbers_pcgenes.tex').write_text('\n'.join(lines) + '\n')
    body = []
    for sa, sb in pairs:
        cells = []
        for g in GENES:
            r = next((x for x in cp if x['species_a'] == sa and x['species_b'] == sb and x['gene'] == g), None)
            cells.append('--' if r is None else f"{r['min_between']:.2f} ({r['max_within'] if r['max_within'] != '' else '--'})")
        body.append(f"{tex_sp(sa)} / {tex_sp(sb)} & " + ' & '.join(cells) + '\\\\')
    thr = ' / '.join(f"{v[2]:.1f}\\%" for v in GENES.values())
    table = f"""\\begin{{table}}[htbp]\\centering
\\begin{{threeparttable}}
\\caption{{Protein-coding genes on the rDNA species complexes: smallest distance (\\%) between the two species, with the
largest distance within either species in brackets.}}
\\label{{tab:pcgenes}}
\\footnotesize
\\begin{{tabular}}{{p{{0.44\\linewidth}}p{{0.14\\linewidth}}p{{0.14\\linewidth}}p{{0.14\\linewidth}}}}
\\toprule
Pair & Glomalin & RPB1 & H$^+$-ATPase\\\\
\\midrule
{chr(10).join(body)}
\\bottomrule
\\end{{tabular}}
\\begin{{tablenotes}}\\footnotesize
\\item Sequences of Stefani et al.\\ (2025). Indicative species thresholds from that study: {thr} (glomalin / RPB1 /
H$^+$-ATPase). Distances are edlib infix edit distances including any intron, so they are comparable with each other but
only approximately with those thresholds. Records are joined to species through the reference culture names; records whose
GenBank name disagrees with the reference name of their culture are left out ({len(summary['conflicts'])}).
\\end{{tablenotes}}
\\end{{threeparttable}}
\\end{{table}}
"""
    (paper / 'tables').mkdir(exist_ok=True)
    (paper / 'tables' / 'tableS2_pcgenes.tex').write_text(table)
    print(f'wrote {paper}/numbers_pcgenes.tex and tables/tableS2_pcgenes.tex')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pairs', type=Path, required=True, help='v4.pairs.tsv from species_groups.py')
    ap.add_argument('--complexes', type=Path, required=True, help='v4.complexes.tsv from species_groups.py')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--email', required=True, help='contact e-mail for NCBI E-utilities')
    ap.add_argument('--api-key', help='NCBI API key (optional; faster)')
    ap.add_argument('--cultures', type=Path, default=Path(__file__).parent / 'xw' / 'crosswalk_v4final.cultures.tsv',
                    help='reference culture table with culture and organism columns')
    ap.add_argument('--paper', type=Path, help='paper directory: write numbers_pcgenes.tex and tables/tableS2_pcgenes.tex')
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    pairs = list(csv.DictReader(open(a.pairs), delimiter='\t'))
    ref_species = sorted({p['species_a'] for p in pairs} | {p['species_b'] for p in pairs})
    complex_of = {}
    for r in csv.DictReader(open(a.complexes), delimiter='\t'):
        for sp in r['species'].split(';'):
            complex_of[sp] = r['complex']
    by_norm = {norm(s): s for s in ref_species if norm(s)}

    def to_ref(org):
        n = norm(org)
        if not n:
            return None, ''
        if n in by_norm:
            return by_norm[n], ''
        for (g, ep), s in by_norm.items():  # same genus, epithet differing only in its ending
            if g == n[0] and len(ep) >= 6 and ep[:-2] == n[1][:len(ep) - 2] and abs(len(ep) - len(n[1])) <= 2:
                return s, f'{org} read as {s}'
        return None, ''

    label_org = {}
    if a.cultures and a.cultures.exists():
        for r in csv.DictReader(open(a.cultures), delimiter='\t'):
            if ' | ' in r['culture']:   # a mixed culture split into lineages keeps one label for both
                label_org.setdefault(r['culture'].split(' | ', 1)[1], set()).add(r['organism'])
    records, notes, unmapped, conflicts = {}, set(), defaultdict(int), []
    for gene, (first, last, _) in GENES.items():
        cache = a.out / f'{gene.replace("+", "plus")}.gb'
        if not cache.exists():
            print(f'downloading {gene} ({first}-{last}) ...', flush=True)
            cache.write_text(fetch(acc_range(first, last), a.email, a.api_key))
        recs = parse_gb(cache.read_text())
        for r in recs:
            r['gene'] = gene
            r['species'], note = to_ref(r['organism'])
            if note:
                notes.add(note)
            orgs = label_org.get(r['culture'], set())
            r['reference_organism'] = ' / '.join(sorted(orgs))
            if orgs:
                ref_sps = {o if o in ref_species else None for o in orgs}
                if r['species'] not in ref_sps:
                    conflicts.append({'gene': gene, 'accession': r['accession'], 'culture': r['culture'],
                                      'genbank': r['organism'], 'reference': r['reference_organism']})
                    r['species'] = None
                    continue
            if not r['species']:
                unmapped[(gene, r['organism'])] += 1
        records[gene] = recs
        print(f'{gene}: {len(recs)} records, {sum(1 for r in recs if r["species"])} mapped to reference species')
    for n in sorted(notes):
        print(f'  name match: {n}')
    if unmapped:
        print('  organisms not in the reference: ' + '; '.join(f'{o} ({n}, {g})' for (g, o), n in sorted(unmapped.items())))
    for c in conflicts:
        rec = next(r for r in records[c['gene']] if r['accession'] == c['accession'])
        refs = [r for r in records[c['gene']] if r['organism'] in c['reference'].split(' / ') and r['accession'] != c['accession']
                and r['culture'] != c['culture']]
        c['nearest_reference_pct'] = round(min(dist(rec['seq'], r['seq']) for r in refs), 2) if refs else None
    if conflicts:
        print(f'  left out, GenBank name differs from the reference name of the culture: {len(conflicts)} records')
        for c in conflicts:
            print(f"      {c['gene']} {c['accession']} {c['culture']}: GenBank {c['genbank']}, reference {c['reference']}"
                  f" (nearest record named as the reference: {c['nearest_reference_pct']}%)")

    with open(a.out / 'pcgenes.records.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['gene', 'accession', 'organism', 'culture', 'reference_culture_name', 'species_used', 'length'])
        for gene, recs in records.items():
            for r in recs:
                w.writerow([gene, r['accession'], r['organism'], r['culture'], r['reference_organism'],
                            r['species'] or '', len(r['seq'])])

    rows = []
    for gene, recs in records.items():
        thr = GENES[gene][2]
        by_sp = defaultdict(list)
        for r in recs:
            if r['species']:
                by_sp[r['species']].append(r)
        intra = {s: max((dist(x['seq'], y['seq']) for x, y in itertools.combinations(rs, 2)), default=None)
                 for s, rs in by_sp.items()}
        for p in pairs:
            sa, sb = p['species_a'], p['species_b']
            if sa not in by_sp or sb not in by_sp:
                continue
            inter = min(dist(x['seq'], y['seq']) for x in by_sp[sa] for y in by_sp[sb])
            spread = max([v for v in (intra[sa], intra[sb]) if v is not None], default=None)
            same_genus = norm(sa)[0] == norm(sb)[0]
            cx = complex_of.get(sa) if complex_of.get(sa) and complex_of.get(sa) == complex_of.get(sb) else ''
            rows.append({'gene': gene, 'species_a': sa, 'species_b': sb, 'complex': cx,
                         'rdna_inseparable': p['inseparable'], 'same_genus': int(same_genus),
                         'records_a': len(by_sp[sa]), 'records_b': len(by_sp[sb]),
                         'min_between': round(inter, 2), 'max_within': '' if spread is None else round(spread, 2),
                         'threshold': thr, 'separated': int(inter > thr),
                         'cleanly_separated': '' if spread is None else int(inter > thr and inter > spread)})
    cols = list(rows[0]) if rows else []
    with open(a.out / 'pcgenes.pairs.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

    print('\nrDNA species complexes (pairs inseparable on SSU, ITS and LSU)')
    print(f"  {'pair':58s}{'gene':11s}{'between':>9s}{'within':>8s}{'thr':>6s}  verdict")
    summary = {'complex_pairs': [], 'congeneric_baseline': {}, 'name_matches': sorted(notes), 'conflicts': conflicts,
               'unmapped': {f'{g}: {o}': n for (g, o), n in unmapped.items()}}
    for r in sorted((r for r in rows if r['rdna_inseparable'] == '1'), key=lambda r: (r['complex'], r['species_a'], r['species_b'], r['gene'])):
        verdict = ('separated beyond both spreads' if r['cleanly_separated'] == 1 else
                   'above threshold' if r['separated'] else 'not separated')
        print(f"  {r['species_a'] + ' / ' + r['species_b']:58s}{r['gene']:11s}{r['min_between']:>8.2f}%"
              f"{(str(r['max_within']) + '%') if r['max_within'] != '' else '-':>8s}{r['threshold']:>5.1f}%  {verdict}")
        summary['complex_pairs'].append({**r, 'verdict': verdict})
    missing = [p for p in pairs if p['inseparable'] == '1' and not any(
        r['species_a'] == p['species_a'] and r['species_b'] == p['species_b'] for r in rows)]
    for p in missing:
        print(f"  {p['species_a'] + ' / ' + p['species_b']:58s}no protein-coding records for one or both species")
    print('\nbaseline: congeneric pairs that rDNA already separates')
    for gene in GENES:
        base = [r for r in rows if r['gene'] == gene and r['same_genus'] and r['rdna_inseparable'] != '1']
        if base:
            sep = sum(r['separated'] for r in base)
            mins = sorted(r['min_between'] for r in base)
            print(f"  {gene:11s}{len(base):>4d} pairs: above threshold {sep} ({100 * sep / len(base):.0f}%), "
                  f"median distance {statistics.median(mins):.2f}%, lowest {mins[0]:.2f}%")
            summary['congeneric_baseline'][gene] = {'pairs': len(base), 'separated': sep,
                                                     'median_between': statistics.median(mins), 'lowest_between': mins[0]}
    summary['missing_pairs'] = [f"{p['species_a']} / {p['species_b']}" for p in missing]
    (a.out / 'pcgenes.summary.json').write_text(json.dumps(summary, indent=2, default=str) + '\n')
    if a.paper:
        write_paper(a.paper, rows, summary, records)
    print(f'\ntables: {a.out}/pcgenes.pairs.tsv, pcgenes.records.tsv, pcgenes.summary.json')


if __name__ == '__main__':
    sys.exit(main())
