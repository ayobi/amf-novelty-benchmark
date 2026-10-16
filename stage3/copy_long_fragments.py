#!/usr/bin/env python3
"""Thresholds on long fragments (exploratory, no decision rule): L. Tedersoo's follow-up to II.2-3 (13 October 2026),
"How would the entire Kruger fragment, taken together, work? And how would the entire 4.5-kb fragment work?"

Regions, on the main-type copies of culture reference v4:
  full_its  the full ITS, for comparison (as in copy_thresholds.py)
  kruger    the Kruger fragment as far as the copies reach: from the SSUmAf1 site (Kruger et al. 2009) to the end of
            the copy at LSUmBr. The copies stop at LSUmBr, so the last stretch of the SSUmAf-LSUmAr amplicon is missing.
  full      the whole copy (AML1 to LSUmBr, primers trimmed), the longest linked fragment the cultures give. It is not
            the 4.5-kb fragment, which only long reads such as EUKARYOME's cover.

Distances, classes and the within-species spread are those of copy_thresholds.py (edlib infix edit distance as % of
the shorter sequence; 99th percentile of the distance to the nearest copy of another culture of the same species).
Added here, a threshold sweep: for each distance t, the share of copies that t splits from their own genome (more
than t from every other copy of their culture) or from their species (more than t from every other copy of it), and
the share it joins to another species of the same genus (a congener within t), or to another genus. Congeners are
counted twice: all of them, and without the partners of the paper's four complexes (--complexes), which no marker
separates.

  PYTHONPATH=<benchmark repo>/stage3 python3 stage3/copy_long_fragments.py --ref reference/culture_ref_v4 \\
      --complexes paper/analysis/v4.complexes.tsv --out benchmarks/copy_checks/thresholds_long

Writes long.limits.tsv, long.nearest.tsv, long.quantiles.tsv, long.sweep.tsv, long.summary.tsv and
fig_thresholds_long.{pdf,png}.
"""
import argparse
import collections
import os
import statistics
import sys

import edlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arms import is_named, read_fasta, read_tsv, write_tsv  # noqa: E402
from copy_thresholds import GENUS_SYN, PAPER_SPREAD, dist, q, spread  # noqa: E402

REGIONS = [('full_its', 'Full ITS'), ('kruger', 'Krüger fragment (SSUmAf–LSUmBr)'),
           ('full', 'Whole sequence (AML1–LSUmBr)')]
CLASSES = ['same culture', 'same species', 'same complex', 'same genus', 'other genus']
SSUMAF1 = 'TGGGTAATCTTTTGAAACTTYA'            # Kruger et al. 2009; SSUmAf2 differs by 2 positions
IUPAC = {'Y': 'CT', 'R': 'AG', 'D': 'AGT', 'N': 'ACGT', 'W': 'AT', 'S': 'CG', 'K': 'GT', 'M': 'AC', 'H': 'ACT',
         'B': 'CGT', 'V': 'ACG'}
EQ = [(k, c) for k, v in IUPAC.items() for c in v]
SWEEP = [round(0.05 * i, 2) for i in range(0, 241)]          # 0-12%


def kruger(full, max_edits=3):
    """The copy from the SSUmAf1 site to its end; copies without the site (more than max_edits) are left out."""
    out, miss = {}, []
    for acc, s in full.items():
        r = edlib.align(SSUMAF1, s, mode='HW', task='locations', additionalEqualities=EQ)
        if r['editDistance'] <= max_edits:
            out[acc] = s[r['locations'][0][0]:]
        else:
            miss.append(acc)
    return out, miss


def load(ref):
    full = read_fasta(os.path.join(ref, 'fasta', 'full.fasta'))
    kr, miss = kruger(full)
    seq = {'full_its': read_fasta(os.path.join(ref, 'fasta', 'full_its.fasta')), 'kruger': kr, 'full': full}
    return seq, miss


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--complexes', required=True, help='species_groups.py complexes table (paper/analysis)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--replot', action='store_true',
                    help='redraw the figure from the tables already in --out (no distances recomputed)')
    a = ap.parse_args()
    if a.replot:
        def num(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return v
        SPREAD = {r['region']: float(r['within_species_spread_p99'])
                  for r in read_tsv(os.path.join(a.out, 'long.limits.tsv'))}
        sweep = [{k: (num(v) if k != 'region' else v) for k, v in r.items()}
                 for r in read_tsv(os.path.join(a.out, 'long.sweep.tsv'))]
        summ = [{k: (num(v) if k not in ('region', 'label', 'threshold') else v) for k, v in r.items()}
                for r in read_tsv(os.path.join(a.out, 'long.summary.tsv'))]
        figures(a, sweep, summ, SPREAD)
        return
    os.makedirs(a.out, exist_ok=True)
    cul = {r['culture']: r for r in read_tsv(os.path.join(a.ref, 'cultures.tsv'))}
    rows = [r for r in read_tsv(os.path.join(a.ref, 'copies.tsv'))
            if r['status'] == 'retained' and r['copy_class'] == 'main']
    seq, miss = load(a.ref)
    named = {c for c, r in cul.items() if is_named(r['organism']) and r.get('name_status', 'ok') == 'ok'}
    sp_of = {c: cul[c]['organism'] for c in named}
    cx = {}
    for r in read_tsv(a.complexes):
        for sp in r['species'].split(';'):
            cx[sp] = r['complex']
    caches = {m: {} for m, _ in REGIONS}
    SPREAD = {m: spread(rows, seq[m], sp_of, caches[m]) for m, _ in REGIONS}
    if abs(SPREAD['full_its'] - PAPER_SPREAD['full_its']) > 5e-4:
        raise SystemExit(f"full ITS spread {SPREAD['full_its']:.4f}, the paper has {PAPER_SPREAD['full_its']}")
    lens = {m: [len(s) for s in seq[m].values()] for m, _ in REGIONS}
    write_tsv(os.path.join(a.out, 'long.limits.tsv'),
              [{'region': m, 'copies_with_region': len(seq[m]), 'median_length': statistics.median(lens[m]),
                'min_length': min(lens[m]), 'max_length': max(lens[m]),
                'within_species_spread_p99': round(SPREAD[m], 4)} for m, _ in REGIONS],
              ['region', 'copies_with_region', 'median_length', 'min_length', 'max_length',
               'within_species_spread_p99'])

    def genus(r):
        return GENUS_SYN.get(r['genus'], r['genus'])

    def relation(r1, r2):
        if r1['culture'] == r2['culture']:
            return 'same culture'
        if r1['culture'] in named and r2['culture'] in named and r1['organism'] == r2['organism']:
            return 'same species'
        if (r1['culture'] in named and r2['culture'] in named and r1['organism'] in cx
                and cx[r1['organism']] == cx.get(r2['organism'])):
            return 'same complex'
        if genus(r1) == genus(r2):
            return 'same genus' if (r1['culture'] in named and r2['culture'] in named) else None
        return 'other genus'

    by_order = collections.defaultdict(list)
    for r in rows:
        by_order[r['order']].append(r)
    near = []
    for order, rs in by_order.items():
        for r1 in rs:
            best = {}
            for r2 in rs:
                if r2['accession'] == r1['accession']:
                    continue
                rel = relation(r1, r2)
                if rel is None:
                    continue
                for m, _ in REGIONS:
                    s1, s2 = seq[m].get(r1['accession']), seq[m].get(r2['accession'])
                    if s1 is None or s2 is None:
                        continue
                    d = dist(s1, s2, caches[m])
                    if (m, rel) not in best or d < best[(m, rel)]:
                        best[(m, rel)] = d
            for (m, rel), d in best.items():
                near.append({'accession': r1['accession'], 'culture': r1['culture'], 'organism': r1['organism'],
                             'order': order, 'region': m, 'class': rel, 'nearest_pct': round(d, 4)})
    write_tsv(os.path.join(a.out, 'long.nearest.tsv'), near,
              ['accession', 'culture', 'organism', 'order', 'region', 'class', 'nearest_pct'])

    nv = collections.defaultdict(list)                         # (region, class) -> nearest distances
    for x in near:
        nv[(x['region'], x['class'])].append(x['nearest_pct'])
    qrows = []
    for m, _ in REGIONS:
        for c in CLASSES:
            v = nv[(m, c)]
            if v:
                qrows.append({'region': m, 'class': c, 'copies': len(v), 'median': round(statistics.median(v), 3),
                              'p95': round(q(v, 95), 3), 'p99': round(q(v, 99), 3),
                              'within_spread': f"{100 * sum(x <= SPREAD[m] for x in v) / len(v):.1f}",
                              'within_1.5': f"{100 * sum(x <= 1.5 for x in v) / len(v):.1f}",
                              'within_3.0': f"{100 * sum(x <= 3.0 for x in v) / len(v):.1f}"})
    write_tsv(os.path.join(a.out, 'long.quantiles.tsv'), qrows,
              ['region', 'class', 'copies', 'median', 'p95', 'p99', 'within_spread', 'within_1.5', 'within_3.0'])

    # sweep: split from own genome (nearest culture-mate > t), joined to a congener / another genus (nearest <= t)
    def pct(v, f):
        return 100 * sum(1 for x in v if f(x)) / len(v)
    per = collections.defaultdict(dict)                        # (region, accession) -> class -> nearest
    for x in near:
        per[(x['region'], x['accession'])][x['class']] = x['nearest_pct']
        per[(x['region'], x['accession'])]['named'] = x['culture'] in named

    def vals(m, f):
        return [v for (mm, _), d in per.items() if mm == m for v in [f(d)] if v is not None]

    def mn(d, ks):
        v = [d[k] for k in ks if k in d]
        return min(v) if v else None
    sweep, summ = [], []
    for m, lab in REGIONS:
        mate = vals(m, lambda d: d.get('same culture'))
        own = vals(m, lambda d: mn(d, ['same culture', 'same species']) if d['named'] else None)
        cong = vals(m, lambda d: mn(d, ['same genus', 'same complex']))
        cong_x = vals(m, lambda d: d.get('same genus'))
        oth = vals(m, lambda d: d.get('other genus'))

        def at(t):
            return {'split': pct(mate, lambda x: x > t), 'split_sp': pct(own, lambda x: x > t),
                    'cong': pct(cong, lambda x: x <= t), 'cong_x': pct(cong_x, lambda x: x <= t),
                    'oth': pct(oth, lambda x: x <= t)}
        curve = []
        for t in SWEEP:
            v = at(t)
            row = {'region': m, 'threshold_pct': t, 'identity_pct': round(100 - t, 2),
                   'split_from_genome': round(v['split'], 2), 'split_from_species': round(v['split_sp'], 2),
                   'joined_to_congener': round(v['cong'], 2), 'joined_to_congener_excl_complex': round(v['cong_x'], 2),
                   'joined_to_other_genus': round(v['oth'], 2)}
            sweep.append(row)
            curve.append(row)
        best = min(curve, key=lambda r: (max(r['split_from_genome'], r['joined_to_congener']), r['threshold_pct']))
        best_x = min(curve, key=lambda r: (max(r['split_from_genome'], r['joined_to_congener_excl_complex']),
                                           r['threshold_pct']))
        for name, t in [('98.5% identity', 1.5), ('97% identity', 3.0), ('within-species spread', SPREAD[m]),
                        ('best single threshold', best['threshold_pct']),
                        ('best single threshold, complexes aside', best_x['threshold_pct'])]:
            v = at(t)
            summ.append({'region': m, 'label': lab, 'threshold': name, 'distance_pct': round(t, 3),
                         'identity_pct': round(100 - t, 2), 'copies_with_culture_mate': len(mate),
                         'copies_with_conspecific': len(own), 'copies_with_congener': len(cong),
                         'copies_with_congener_excl_complex': len(cong_x), 'copies_with_other_genus': len(oth),
                         'split_from_genome_pct': round(v['split'], 1), 'split_from_species_pct': round(v['split_sp'], 1),
                         'joined_to_congener_pct': round(v['cong'], 1),
                         'joined_to_congener_excl_complex_pct': round(v['cong_x'], 1),
                         'joined_to_other_genus_pct': round(v['oth'], 1)})
    write_tsv(os.path.join(a.out, 'long.sweep.tsv'), sweep, list(sweep[0]))
    write_tsv(os.path.join(a.out, 'long.summary.tsv'), summ, list(summ[0]))
    for r in summ:
        print(f"{r['region']:9s} {r['threshold']:38s} {r['identity_pct']:6.2f}%  split genome "
              f"{r['split_from_genome_pct']:5.1f} species {r['split_from_species_pct']:5.1f}  congener "
              f"{r['joined_to_congener_pct']:5.1f} (no complexes {r['joined_to_congener_excl_complex_pct']:5.1f})"
              f"  other genus {r['joined_to_other_genus_pct']:4.1f}")
    print('copies:', {r['region']: (r['copies_with_culture_mate'], r['copies_with_conspecific'],
                                    r['copies_with_congener'], r['copies_with_congener_excl_complex'],
                                    r['copies_with_other_genus']) for r in summ})
    if miss:
        print(f'{len(miss)} copies lack the SSUmAf1 site and are left out of the Kruger fragment')

    figures(a, sweep, summ, SPREAD)



def figures(a, sweep, summ, SPREAD):
    # figure: the two error curves per region
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'Inter', 'font.size': 8, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b',
                         'xtick.color': '#52514e', 'ytick.color': '#52514e', 'axes.grid': True,
                         'grid.color': '#e6e5e0', 'grid.linewidth': 0.6})
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), sharey=True)
    XMAX = {'full_its': 12, 'kruger': 6, 'full': 4}
    for ax, (m, lab) in zip(axes, REGIONS):
        cur = [r for r in sweep if r['region'] == m and r['threshold_pct'] <= XMAX[m]]
        xs = [r['threshold_pct'] for r in cur]
        ax.plot(xs, [r['split_from_genome'] for r in cur], color='#2a78d6', lw=1.8, zorder=3)
        ax.plot(xs, [r['joined_to_congener'] for r in cur], color='#eb6834', lw=1.8, zorder=3)
        b = next(r for r in summ if r['region'] == m and r['threshold'] == 'best single threshold')
        ax.axvline(SPREAD[m], color='#52514e', lw=0.9, ls='--', zorder=2)
        ax.plot([b['distance_pct']], [max(b['split_from_genome_pct'], b['joined_to_congener_pct'])], 'o', ms=5,
                color='#0b0b0b', mec='#fcfcfb', mew=0.8, zorder=4)
        ax.text(SPREAD[m], 101, f' limit {SPREAD[m]:.1f}%', fontsize=7, color='#52514e', va='bottom')
        ax.set_title(f"{lab}\nbest {b['identity_pct']:.2f}%: {b['split_from_genome_pct']:.0f}% split, "
                     f"{b['joined_to_congener_pct']:.0f}% joined", fontsize=7.8, loc='left', pad=14)
        ax.set_xlim(0, XMAX[m])
        ax.set_ylim(0, 100)
        ax.set_xlabel('distance threshold (%)')
    axes[0].set_ylabel('% of main-type sequences')
    axes[0].text(0.35, 62, 'split from\nown culture', color='#2a78d6', fontsize=7.5)
    axes[0].text(6.6, 12, 'joined to\nanother species\nof the genus', color='#eb6834', fontsize=7.5)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(a.out, f'fig_thresholds_long.{ext}'), dpi=220)


if __name__ == '__main__':
    main()
