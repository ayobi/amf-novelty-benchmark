#!/usr/bin/env python3
"""Threshold exploration (exploratory; no decision rule): how far apart are copies at each level of relatedness,
per region, and how do ITS, SSU and LSU divergence relate?

1. For every main-type copy, the distance to its nearest main-type copy in four classes, per region (the SSU
   NS31-AM1 window on which MaarjAM virtual taxa are defined, the whole SSU flank, full ITS, LSU LROR-FLR2): same
   culture; same species, other culture; other species of the same genus; other genus of the same order.
   Distances: edlib infix edit distance, % of the shorter sequence (as in the paper).
2. For every pair of cultures in one order, the distance between their main-type medoids on each region,
   labelled same species / same genus / other genus.

Reference lines: the within-species spread (99th percentile of the distance from a main-type copy to the nearest
main-type copy of another culture of the same species, unique sequences per culture, exclusive method), computed
exactly as species_groups.py computes it (the three paper markers are checked against the paper's limits), and
the database identity thresholds on the windows they apply to (VT 97% on NS31-AM1, SH 98.5% on full ITS, V18 97%
on LROR-FLR2). The VT threshold is not drawn on the whole SSU flank: VTs are defined on NS31-AM1, which varies
more than the whole flank.

  python3 stage3/copy_thresholds.py --ref reference/culture_ref_v4 --out benchmarks/copy_checks/thresholds

Writes thresholds.nearest.tsv, thresholds.pairs.tsv, thresholds.quantiles.tsv, thresholds.limits.tsv,
fig_thresholds_nearest.{pdf,png}, fig_thresholds_regions.{pdf,png}.
"""
import argparse
import collections
import itertools
import os
import statistics
import sys

import edlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arms import is_named, read_fasta, read_tsv, write_tsv  # noqa: E402

REGIONS = [('SSU_VT_NS31-AM1', 'SSU NS31–AM1 (VT fragment)'), ('ssu_flank', 'SSU region (whole)'),
           ('full_its', 'Full ITS'), ('LSU_LROR-FLR2', 'LSU LR0R–FLR2')]
CLASSES = ['same culture', 'same species', 'same genus', 'other genus']
COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']       # validated categorical slots 1-4 (light)
GENUS_SYN = {'Rhizoglomus': 'Rhizophagus'}
# the paper's within-species spread (species_groups.py on culture reference v4, 99th percentile); recomputed and checked
PAPER_SPREAD = {'full_its': 5.988, 'LSU_LROR-FLR2': 3.027, 'ssu_flank': 0.8065}
# database identity thresholds, on the windows they are defined on
DB = {'SSU_VT_NS31-AM1': (3.0, 'VT 97%'), 'full_its': (1.5, 'SH 98.5%'), 'LSU_LROR-FLR2': (3.0, 'V18 97%')}


def dist(a, b, cache):
    if a == b:
        return 0.0
    k = (a, b) if a < b else (b, a)
    v = cache.get(k)
    if v is None:
        q, t = (a, b) if len(a) <= len(b) else (b, a)
        v = cache[k] = 100 * edlib.align(q, t, mode='HW', task='distance')['editDistance'] / len(q)
    return v


def q(vals, p):
    return statistics.quantiles(vals, n=1000, method='exclusive')[int(p * 10) - 1] if len(vals) >= 2 else ''


def spread(rows, seq, sp_of, cache):
    """species_groups.py: unique (sequence, culture) items per species, nearest item of another culture."""
    by_sp = collections.defaultdict(set)
    for r in rows:
        if r['culture'] in sp_of and r['accession'] in seq:
            by_sp[sp_of[r['culture']]].add((seq[r['accession']], r['culture']))
    vals = []
    for items in by_sp.values():
        if len({c for _, c in items}) < 2:
            continue
        for s1, c1 in items:
            vals.append(min(dist(s1, s2, cache) for s2, c2 in items if c2 != c1))
    return statistics.quantiles(vals, n=100, method='exclusive')[98]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--replot', action='store_true',
                    help='redraw the figures from the tables already in --out (no distances recomputed)')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.replot:
        def num(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return v
        SPREAD = {r['region']: float(r['within_species_spread_p99'])
                  for r in read_tsv(os.path.join(a.out, 'thresholds.limits.tsv'))}
        near = [dict(r, nearest_pct=float(r['nearest_pct'])) for r in read_tsv(os.path.join(a.out, 'thresholds.nearest.tsv'))]
        pairs = [{k: (num(v) if k not in ('order', 'culture_a', 'culture_b', 'class') else v) for k, v in r.items()}
                 for r in read_tsv(os.path.join(a.out, 'thresholds.pairs.tsv'))]
        qrows = read_tsv(os.path.join(a.out, 'thresholds.quantiles.tsv'))
        figures(a, near, pairs, SPREAD, qrows)
        return
    cul = {r['culture']: r for r in read_tsv(os.path.join(a.ref, 'cultures.tsv'))}
    rows = [r for r in read_tsv(os.path.join(a.ref, 'copies.tsv'))
            if r['status'] == 'retained' and r['copy_class'] == 'main']
    seq = {m: read_fasta(os.path.join(a.ref, 'fasta', f'{m}.fasta')) for m, _ in REGIONS}
    named = {c for c, r in cul.items() if is_named(r['organism']) and r.get('name_status', 'ok') == 'ok'}
    sp_of = {c: cul[c]['organism'] for c in named}
    caches = {m: {} for m, _ in REGIONS}

    SPREAD = {m: spread(rows, seq[m], sp_of, caches[m]) for m, _ in REGIONS}
    for m, v in PAPER_SPREAD.items():
        if abs(SPREAD[m] - v) > 5e-4:
            raise SystemExit(f'within-species spread for {m} is {SPREAD[m]:.4f}, the paper has {v}')
    write_tsv(os.path.join(a.out, 'thresholds.limits.tsv'),
              [{'region': m, 'within_species_spread_p99': round(SPREAD[m], 4),
                'database_threshold': DB[m][1] if m in DB else ''} for m, _ in REGIONS],
              ['region', 'within_species_spread_p99', 'database_threshold'])

    def genus(r):
        g = r['genus']
        return GENUS_SYN.get(g, g)

    def relation(r1, r2):
        if r1['culture'] == r2['culture']:
            return 'same culture'
        if r1['culture'] in named and r2['culture'] in named and r1['organism'] == r2['organism']:
            return 'same species'
        if genus(r1) == genus(r2):
            return 'same genus' if (r1['culture'] in named and r2['culture'] in named) else None
        return 'other genus'

    by_order = collections.defaultdict(list)
    for r in rows:
        by_order[r['order']].append(r)

    # 1. nearest copy per class (a region is skipped for a copy that lacks that window)
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
                    k = (m, rel)
                    if k not in best or d < best[k]:
                        best[k] = d
            for (m, rel), d in best.items():
                near.append({'accession': r1['accession'], 'culture': r1['culture'], 'organism': r1['organism'],
                             'order': order, 'region': m, 'class': rel, 'nearest_pct': round(d, 3)})
    write_tsv(os.path.join(a.out, 'thresholds.nearest.tsv'), near,
              ['accession', 'culture', 'organism', 'order', 'region', 'class', 'nearest_pct'])
    qrows = []
    for m, _ in REGIONS:
        for c in CLASSES:
            v = [x['nearest_pct'] for x in near if x['region'] == m and x['class'] == c]
            if not v:
                continue
            qrows.append({'region': m, 'class': c, 'copies': len(v), 'median': round(statistics.median(v), 3),
                          'p95': round(q(v, 95), 3), 'p99': round(q(v, 99), 3),
                          'below_spread': f"{100 * sum(x <= SPREAD[m] for x in v) / len(v):.1f}",
                          'below_db': f"{100 * sum(x <= DB[m][0] for x in v) / len(v):.1f}" if m in DB else ''})
    write_tsv(os.path.join(a.out, 'thresholds.quantiles.tsv'), qrows,
              ['region', 'class', 'copies', 'median', 'p95', 'p99', 'below_spread', 'below_db'])

    # 2. culture medoids (by full ITS), pairs within an order
    med = {}
    for c in sorted({r['culture'] for r in rows}):
        rs = [r for r in rows if r['culture'] == c]
        med[c] = min(rs, key=lambda r: sum(dist(seq['full_its'][r['accession']], seq['full_its'][s['accession']],
                                                caches['full_its']) for s in rs))
    pairs = []
    for order, rs in by_order.items():
        cs = sorted({r['culture'] for r in rs})
        for c1, c2 in itertools.combinations(cs, 2):
            r1, r2 = med[c1], med[c2]
            rel = relation(r1, r2)
            if rel is None or rel == 'same culture':
                continue
            p = {'order': order, 'culture_a': c1, 'culture_b': c2, 'class': rel}
            for m, _ in REGIONS:
                s1, s2 = seq[m].get(r1['accession']), seq[m].get(r2['accession'])
                p[m] = round(dist(s1, s2, caches[m]), 3) if s1 and s2 else ''
            pairs.append(p)
    write_tsv(os.path.join(a.out, 'thresholds.pairs.tsv'), pairs,
              ['order', 'culture_a', 'culture_b', 'class'] + [m for m, _ in REGIONS])

    figures(a, near, pairs, SPREAD, qrows)



def figures(a, near, pairs, SPREAD, qrows):
    # figures
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e',
                         'ytick.color': '#52514e', 'axes.grid': True, 'grid.color': '#e6e5e0', 'grid.linewidth': 0.6})
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), sharey=True)
    XMAX = {'SSU_VT_NS31-AM1': 12, 'ssu_flank': 8, 'full_its': 40, 'LSU_LROR-FLR2': 30}
    YTOP = 1.16                                   # headroom above 1.0 for the label of the two leftmost curves
    panels = []
    for i, (ax, (m, lab)) in enumerate(zip(axes.flat, REGIONS)):
        curves = {}
        for c, col in zip(CLASSES, COLORS):
            v = sorted(x['nearest_pct'] for x in near if x['region'] == m and x['class'] == c)
            if not v:
                continue
            y = [(k + 1) / len(v) for k in range(len(v))]
            ax.step(v, y, where='post', color=col, lw=1.6, label=c, zorder=3)
            curves[c] = v
        lines = [SPREAD[m]]
        ax.vlines(SPREAD[m], 0, 1, color='#52514e', lw=0.9, ls='--', label='within-species limit (99th percentile)',
                  zorder=2)
        note = f'limit {SPREAD[m]:.1f}%'
        if m in DB:
            ax.vlines(DB[m][0], 0, 1, color='#52514e', lw=0.9, ls=':', label='database threshold', zorder=2)
            lines.append(DB[m][0])
            note += f' · {DB[m][1]}'
        ax.text(0.98, 0.03, note, transform=ax.transAxes, fontsize=6, color='#52514e', ha='right', va='bottom')
        ax.set_xlim(0, XMAX[m])
        ax.set_ylim(0, YTOP)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_title(f'({"abcd"[i]}) {lab}', fontsize=8.5, loc='left', pad=4)
        ax.set_xlabel('Distance to nearest sequence (%)')
        if i % 2 == 0:
            ax.set_ylabel('Share of main-type sequences')
        panels.append((ax, m, curves, lines))
    h, l = axes.flat[2].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=6.5, loc='lower center', ncol=6, bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()

    def xat(v, yy):
        """x at which the step CDF of sorted values v first reaches share yy."""
        k = min(len(v) - 1, max(0, int(-(-yy * len(v) // 1)) - 1))
        return v[k]

    def place(ax, text, curve, others, lines, pref):
        """Put a label just right of a curve at the y nearest pref where its box crosses no curve or line."""
        t = ax.text(0, 0, text, fontsize=6.5, color='#0b0b0b', va='top', ha='left')
        bb = t.get_window_extent(renderer=rend)
        t.remove()
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        w = bb.width * (x1 - x0) / ax.bbox.width
        hh = bb.height * (y1 - y0) / ax.bbox.height
        pad = 0.012 * (x1 - x0)
        for dy in sorted({round(d / 100, 2) for d in range(-90, 91)}, key=abs):
            yy = pref + dy
            if yy - hh < 0.12 or yy > 0.97:
                continue
            bx0 = xat(curve, yy) + pad
            bx1 = bx0 + w
            if bx1 > x1 * 0.98:
                continue
            band = (yy - hh, yy)
            hit = any(bx0 - pad / 2 <= xl <= bx1 + pad / 2 for xl in lines)
            for v in [curve] + others:
                lo, hi = xat(v, band[0]), xat(v, band[1])
                if v is curve:
                    hi = xat(v, min(1.0, band[1]))
                    hit |= hi > bx0 - pad / 2
                else:
                    hit |= not (hi < bx0 - pad / 2 or lo > bx1 + pad / 2)
            if not hit:
                ax.text(bx0, yy, text, fontsize=6.5, color='#0b0b0b', va='top', ha='left', zorder=4)
                return True
        return False

    for ax, m, curves, lines in panels:
        allc = [curves[c] for c in CLASSES if c in curves]
        for c, pref in (('same genus', 0.35), ('other genus', 0.6)):
            v = curves.get(c)
            if v and not place(ax, c, v, [u for u in allc if u is not v], lines, pref):
                print(f'no free place for the label {c!r} in {m}')
        # the two nearly identical leftmost curves share one label, in the headroom above their top ends
        ax.text(0.012 * ax.get_xlim()[1], 1.015, 'same culture,\nsame species', fontsize=6.5, color='#0b0b0b',
                va='bottom', ha='left')
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(a.out, f'fig_thresholds_nearest.{ext}'), dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for ax, (m, lab) in zip(axes, [('ssu_flank', 'SSU region (%)'), ('LSU_LROR-FLR2', 'LSU LR0R–FLR2 (%)')]):
        for c, col in zip(CLASSES[1:], COLORS[1:]):
            pp = [p for p in pairs if p['class'] == c and p[m] != '']
            ax.scatter([p['full_its'] for p in pp], [p[m] for p in pp], s=8, color=col, alpha=0.55, linewidths=0,
                       label=f'{c} ({len(pp)})')
        ax.axvline(SPREAD['full_its'], color='#52514e', lw=0.8, ls='--')
        ax.axhline(SPREAD[m], color='#52514e', lw=0.8, ls='--')
        tx, ty, ha = (0.98, 0.55, 'right') if m == 'ssu_flank' else (0.02, 0.97, 'left')
        ax.text(tx, ty, f'dashed: within-species limit\nITS {SPREAD["full_its"]:.1f}%, '
                        f'{lab.split(" (")[0]} {SPREAD[m]:.1f}%',
                transform=ax.transAxes, fontsize=6, color='#52514e', ha=ha, va='top')
        if m == 'ssu_flank':
            xx = [0, 40]
            ax.plot(xx, [x / 2 for x in xx], color='#52514e', lw=0.8, ls='--')
            ax.text(30, 15.5, 'SSU = ITS / 2', fontsize=6, color='#52514e', rotation=0)
        ax.set_xlabel('Full ITS (%)')
        ax.set_ylabel(lab)
        ax.set_xlim(0, 45)
    axes[0].legend(frameon=False, fontsize=6.5, loc='upper left', markerscale=1.5)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(a.out, f'fig_thresholds_regions.{ext}'), dpi=200)
    plt.close(fig)
    for m, _ in REGIONS:
        print(f'within-species spread {m}: {SPREAD[m]:.4f}')
    for r in qrows:
        print('\t'.join(str(r[k]) for k in r))
    print(f'pairs of culture medoids: {len(pairs)}')


if __name__ == '__main__':
    main()
