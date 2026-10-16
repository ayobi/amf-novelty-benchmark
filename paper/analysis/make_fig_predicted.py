#!/usr/bin/env python3
"""Figure: predicted against observed in the independent tests. Each point is one test sequence, placed by its
distance to the nearest culture-reference species in the full ITS (x) and in LROR-FLR2 (y); the shaded box is
the within-species spread of the locked predictions. Fill gives the prediction, fixed before the run: inside a
reference species on all three markers (ITS, LSU and SSU flank; filled) or not (open). Colour gives what the
rules called: a known species or species group (orange) or not (grey).

Row 1, LongGloDB: (a) ASVs of species the reference does not hold, each called with its species removed from the
LSU reference (condition B); (b) ASVs of species it holds, LSU reference intact (condition A). With --euk, row 2
draws the EUKARYOME test the same way: (c) reads of taxa v4 does not hold (B where planned, A for the taxa
reported apart), (d) reads of species it holds (A).

Counts are computed from the locked predictions and the score tables and must equal the macros in
numbers_lgdb.tex (and numbers_euk.tex, once written), so the figure cannot drift from the text.

  python3 analysis/make_fig_predicted.py [--euk ../benchmarks/eukaryome --euk-scores PATH/runs_euk/06_scores]
  (from paper/; writes figures/fig_predicted.pdf and .png)
"""
import argparse
import csv
import json
import math
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

INK, INK2, GRID, SURFACE, NONE = '#0b0b0b', '#52514e', '#e6e5e0', '#fcfcfb', '#b9b8b2'
BLUE, ORANGE = '#2a78d6', '#eb6834'
SPREAD = '#eaf1fb'                       # tint of BLUE for the within-species box
plt.rcParams.update({'font.family': 'Inter', 'font.size': 7.5, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2,
                     'mathtext.fontset': 'custom', 'mathtext.it': 'Inter:italic', 'mathtext.rm': 'Inter'})

NOVEL = ('novel species, genus in v4', 'novel lineage, genus not in v4')
CAP_X, CAP_Y = 40.0, 32.0                # where blank distances are drawn: over 30%, or region missing from the
                                         # files the predictions were computed on
XT = [0, 1, 3, 6, 10, 20]
YT = [0, 1, 3, 6, 10, 20]
EDGE = '>30 or\nnone'


def macros(path):
    return dict(re.findall(r'\\newcommand\{\\(\w+)\}\{([^{}]*)\}', open(path).read())) if os.path.exists(path) else {}


def read_tsv(p):
    with open(p, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def val(s, cap):
    return cap if s in ('', None) else float(s)


def is_known(r):
    return r['status'].startswith('known')


def called(r):
    return [s for s in r['placement'].split(';') if s]


def test_rows(pred_path, scores, novel_condition='B'):
    """One record per test sequence: distances, prediction, call (B for absent taxa where a round exists)."""
    pred = {r['asv']: r for r in read_tsv(pred_path)}
    a = {r['accession']: r for r in read_tsv(os.path.join(scores, 'A_primary_v2.copies.tsv'))}
    b = {r['accession']: r for r in read_tsv(os.path.join(scores, 'B_primary_v2.copies.tsv'))}
    out = []
    for x, p in pred.items():
        if p['role'] not in NOVEL + ('known species',):
            continue
        row = b.get(x) if (p['role'] in NOVEL and novel_condition == 'B' and x in b) else a[x]
        inside = set(filter(None, p['inside'].split('; ')))
        out.append({'id': x, 'role': p['role'], 'organism': p['organism'], 'code': p['code'],
                    'its': val(p['ITS_nearest_pct'], CAP_X), 'lsu': val(p['LSU_nearest_pct'], CAP_Y),
                    'inside': inside, 'own_inside': p['own_species_in_v4'] in inside,
                    'known': is_known(row), 'called': called(row), 'kept': is_known(row) and p['organism'] in called(row)})
    return out


def counts(rows):
    nov = [r for r in rows if r['role'] in NOVEL]
    kn = [r for r in rows if r['role'] == 'known species']
    c = {'novel': len(nov), 'novel_inside': sum(bool(r['inside']) for r in nov),
         'novel_inside_known': sum(bool(r['inside']) and r['known'] for r in nov),
         'novel_outside_known': sum(not r['inside'] and r['known'] for r in nov),
         'known': len(kn), 'own_inside': sum(r['own_inside'] for r in kn),
         'own_inside_kept': sum(r['own_inside'] and r['kept'] for r in kn),
         'other_only': sum(bool(r['inside']) and not r['own_inside'] for r in kn),
         'other_only_called_other': sum(bool(r['inside']) and not r['own_inside'] and r['known'] and not r['kept'] for r in kn),
         'none': sum(not r['inside'] for r in kn), 'none_kept': sum(not r['inside'] and r['kept'] for r in kn),
         'none_known': sum(not r['inside'] and r['known'] for r in kn),
         'called_outside_prediction': sum(r['known'] and not (set(r['called']) & r['inside']) for r in rows)}
    return c


def check_lgdb(c, M):
    """The figure's counts against numbers_lgdb.tex; any difference stops the build."""
    n = lambda k: int(M[k])  # noqa: E731
    want = {'novel': n('LgNovelN'), 'novel_inside': n('LgLockedN'), 'novel_inside_known': n('LgFK') - n('LgFKOutsideLocked'),
            'novel_outside_known': n('LgFKOutsideLocked'), 'known': n('LgKnownN'), 'own_inside': n('LgPredInsideOwn'),
            'own_inside_kept': n('LgPredInsideOwnKept'), 'other_only': n('LgWrongSpPredicted'),
            'none': n('LgPredOutsideOwn') - n('LgWrongSpPredicted'),
            'none_kept': n('LgPredOutsideOwn') - n('LgPredOutsideOwnNotKept')}
    bad = {k: (c[k], v) for k, v in want.items() if c[k] != v}
    if bad:
        sys.exit(f'figure counts differ from numbers_lgdb.tex (figure, macro): {bad}')


def sqrt_axis(ax, xcap=True, ycap=True):
    fwd = lambda v: [math.sqrt(max(0.0, t)) for t in v] if hasattr(v, '__iter__') else math.sqrt(max(0.0, v))  # noqa: E731
    import numpy as np
    f = (lambda v: np.sqrt(np.clip(v, 0, None)), lambda v: np.square(v))
    ax.set_xscale('function', functions=f)
    ax.set_yscale('function', functions=f)
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 40)
    ax.set_xticks(XT + [CAP_X], [str(t) for t in XT] + [EDGE], fontsize=6.5)
    ax.set_yticks(YT + [CAP_Y], [str(t) for t in YT] + [EDGE], fontsize=6.5)
    ax.grid(True, color=GRID, lw=0.5, zorder=0)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(GRID)
    ax.tick_params(length=0, pad=2)
    del fwd


def draw(ax, rows, T, size):
    ax.add_patch(Rectangle((0, 0), T['ITS'], T['LSU'], facecolor=SPREAD, edgecolor=BLUE, lw=0.7, ls=(0, (3, 2)), zorder=1))
    order = sorted(rows, key=lambda r: (r['known'], bool(r['inside'])))           # grey first, orange on top
    for r in order:
        col = ORANGE if r['known'] else NONE
        filled = bool(r['inside'])
        ax.scatter([r['its']], [r['lsu']], s=size, marker='o', facecolor=col if filled else 'white', edgecolor=col,
                   linewidths=0.9, zorder=4 if r['known'] else 3, clip_on=False)


def label(ax, text, xy, xytext, ha='left'):
    ax.annotate(text, xy=xy, xytext=xytext, textcoords='data', fontsize=6.4, color=INK, ha=ha, va='center',
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6, shrinkA=1, shrinkB=3), zorder=6)


def centroid(rs):
    return (sum(r['its'] for r in rs) / len(rs), sum(r['lsu'] for r in rs) / len(rs))


def ital(name):
    g, _, e = name.partition(' ')
    return f'$\\it{{{g[0]}.\\,{e}}}$' if e else name


def panel_title(fig, x, y, title, sub):
    fig.text(x, y, title, fontsize=8.5, weight='bold', color=INK, va='top')
    fig.text(x, y - 0.035, sub, fontsize=6.8, color=INK2, va='top')


def summary_box(ax, lines):
    ax.text(0.03, 0.97, '\n'.join(lines), transform=ax.transAxes, ha='left', va='top', fontsize=6.4, color=INK,
            linespacing=1.45, bbox=dict(boxstyle='round,pad=0.45', facecolor=SURFACE, edgecolor=GRID, lw=0.6), zorder=7)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--lgdb', default='../benchmarks/longglodb')
    ap.add_argument('--lgdb-scores', default='../benchmarks/m4_external_holdout/runs_lgdb/06_scores')
    ap.add_argument('--euk', help='benchmarks/eukaryome (locked predictions); with --euk-scores, draws row 2')
    ap.add_argument('--euk-scores', help='runs_euk/06_scores of the EUKARYOME run')
    ap.add_argument('--out', default='figures/fig_predicted')
    a = ap.parse_args()

    M = macros('numbers_lgdb.tex')
    lg = test_rows(os.path.join(a.lgdb, 'predictions.predictions.tsv'), a.lgdb_scores)
    T = json.load(open(os.path.join(a.lgdb, 'predictions.lock')))['T']
    c = counts(lg)
    check_lgdb(c, M)
    two_rows = bool(a.euk and a.euk_scores)
    fig = plt.figure(figsize=(7.2, 7.0 if two_rows else 3.7))
    top = 0.975 if not two_rows else 0.985
    h = 0.63 if not two_rows else 0.31
    y0 = 0.215 if not two_rows else 0.595

    # ---- (a) LongGloDB, species the reference does not hold
    ax = fig.add_axes([0.105, y0, 0.375, h])
    sqrt_axis(ax)
    draw(ax, [r for r in lg if r['role'] in NOVEL], T, 16)
    ks = [r for r in lg if r['role'] in NOVEL and r['inside']]
    if ks:
        x, y = centroid(ks)
        org = ks[0]['organism']
        to = sorted({s for r in ks for s in r['called']})
        label(ax, f"{ital(org)} {ks[0]['code']}\ncalled {', '.join(ital(s) for s in to)} ({len(ks)})", (x, y), (1.6, 0.18))
    ax.set_xlabel('ITS: distance to the nearest reference species (%)', fontsize=7)
    ax.set_ylabel('LSU: distance to the nearest\nreference species (%)', fontsize=7)
    summary_box(ax, [f"Inside a reference species: {c['novel_inside']}, called known {c['novel_inside_known']}",
                     f"Inside none: {c['novel'] - c['novel_inside']}, called known {c['novel_outside_known']}"])
    panel_title(fig, 0.02, top, '(a) LongGloDB: taxa the reference does not hold',
                f"{c['novel']} ASVs, each species removed from the LSU reference")

    # ---- (b) LongGloDB, species the reference holds
    ax = fig.add_axes([0.60, y0, 0.375, h])
    sqrt_axis(ax)
    kn = [r for r in lg if r['role'] == 'known species']
    draw(ax, kn, T, 16)
    zero = [r for r in kn if r['its'] == 0 and r['lsu'] == 0]
    if zero:
        label(ax, f'{len(zero)} ASVs at 0, 0', (0.02, 0.02), (11, 0.3))
    other = [r for r in kn if r['inside'] and not r['own_inside']]
    if other:
        x, y = centroid(other)
        to = sorted({s for r in other for s in r['called']})
        label(ax, f"{ital(other[0]['organism'])} {other[0]['code']}\ncalled {', '.join(ital(s) for s in to)} ({len(other)})",
              (x, y), (11, 1.25))
    ex = [r for r in kn if not r['inside'] and r['known']]
    for r in ex:
        why = 'no ITS distance in the prediction' if r['its'] == CAP_X else 'predicted inside none'
        label(ax, f"{ital(r['organism'])} {r['code']}: kept;\n{why}", (r['its'], r['lsu']), (9, 6.5))
    ax.set_xlabel('ITS: distance to the nearest reference species (%)', fontsize=7)
    summary_box(ax, [f"Inside own species: {c['own_inside']}, kept {c['own_inside_kept']}",
                     f"Inside another species only: {c['other_only']}, called it {c['other_only_called_other']}",
                     f"Inside none: {c['none']}, kept {c['none_kept']}"])
    panel_title(fig, 0.525, top, '(b) LongGloDB: species the reference holds', f"{c['known']} ASVs, LSU reference intact")

    if two_rows:
        eu = test_rows(os.path.join(a.euk, 'predictions.predictions.tsv'), a.euk_scores)
        Te = json.load(open(os.path.join(a.euk, 'predictions.lock')))['T']
        ce = counts(eu)
        y1 = 0.11
        ax = fig.add_axes([0.105, y1, 0.375, h])
        sqrt_axis(ax)
        draw(ax, [r for r in eu if r['role'] in NOVEL], Te, 7)
        ax.set_xlabel('ITS: distance to the nearest reference species (%)', fontsize=7)
        ax.set_ylabel('LSU: distance to the nearest\nreference species (%)', fontsize=7)
        summary_box(ax, [f"Inside a reference species: {ce['novel_inside']}, called known {ce['novel_inside_known']}",
                         f"Inside none: {ce['novel'] - ce['novel_inside']}, called known {ce['novel_outside_known']}"])
        panel_title(fig, 0.02, 0.495, '(c) EUKARYOME: taxa the reference does not hold',
                    f"{ce['novel']} reads; LSU reference without the species or genus where V18 holds it")
        ax = fig.add_axes([0.60, y1, 0.375, h])
        sqrt_axis(ax)
        draw(ax, [r for r in eu if r['role'] == 'known species'], Te, 12)
        ax.set_xlabel('ITS: distance to the nearest reference species (%)', fontsize=7)
        summary_box(ax, [f"Inside own species: {ce['own_inside']}, kept {ce['own_inside_kept']}",
                         f"Inside another species only: {ce['other_only']}, called it {ce['other_only_called_other']}",
                         f"Inside none: {ce['none']}, kept {ce['none_kept']}"])
        panel_title(fig, 0.525, 0.495, '(d) EUKARYOME: species the reference holds', f"{ce['known']} reads, LSU reference intact")
        print('EUKARYOME counts:', ce)

    # ---- shared legend
    handles = [Line2D([], [], ls='', marker='o', ms=4.6, mfc=INK2, mec=INK2, label='predicted inside a reference species'),
               Line2D([], [], ls='', marker='o', ms=4.6, mfc='white', mec=INK2, mew=0.9, label='predicted inside none'),
               Line2D([], [], ls='', marker='s', ms=5, mfc=ORANGE, mec=ORANGE, label='called a known species or group'),
               Line2D([], [], ls='', marker='s', ms=5, mfc=NONE, mec=NONE, label='not called known'),
               Rectangle((0, 0), 1, 1, facecolor=SPREAD, edgecolor=BLUE, lw=0.7, ls=(0, (3, 2)),
                         label='within-species spread (ITS and LSU)')]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=6.6, handletextpad=0.4,
               columnspacing=1.6, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(f'{a.out}.pdf')
    fig.savefig(f'{a.out}.png', dpi=300)
    print('LongGloDB counts:', c)
    print(f'written {a.out}.pdf and .png')


if __name__ == '__main__':
    main()
