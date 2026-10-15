#!/usr/bin/env python3
"""Figure: how the three naming systems relate. (a) Lumping and splitting by each system, over the named cultures
(each culture counted by its most common name, as in section 3.2). (b) The cultures of two genera followed from
their VT through their LSU name to their SH.

Counts in (a) are read from numbers_generated.tex (written by make_numbers.py), so the figure cannot drift from
the text; (b) is drawn from the crosswalk table those counts come from.

  python3 analysis/make_fig_nesting.py        (writes figures/fig_nesting.pdf and .png)
"""
import csv
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, PathPatch, Rectangle  # noqa: E402
from matplotlib.path import Path  # noqa: E402

INK, INK2, GRID, SURFACE, NONE = '#0b0b0b', '#52514e', '#e6e5e0', '#fcfcfb', '#b9b8b2'
BLUE, ORANGE = '#2a78d6', '#eb6834'
SPECIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']         # validated categorical slots 1-4 (light)
plt.rcParams.update({'font.family': 'Inter', 'font.size': 7.5, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2,
                     'mathtext.fontset': 'custom', 'mathtext.it': 'Inter:italic', 'mathtext.rm': 'Inter'})

M = dict(re.findall(r'\\newcommand\{\\(\w+)\}\{([^{}]*)\}', open('numbers_generated.tex').read()))
num = lambda k: int(M[k])  # noqa: E731

fig = plt.figure(figsize=(7.2, 5.3))

# ---------------------------------------------------------------- (a) lumping and splitting
ax = fig.add_axes([0.20, 0.715, 0.74, 0.17])
rows = [('VT (SSU)', 'VT', 'VTs'), ('LSU name', 'LSU', 'LSU names'), ('SH (ITS)', 'SH', 'SHs')]
for i, (label, arm, plural) in enumerate(rows):
    y = len(rows) - 1 - i
    lu, lu_n = num(f'{arm}UnitsMulti'), num(f'{arm}UnitsN')
    sp, sp_n = num(f'{arm}SpeciesMulti'), num(f'{arm}SpeciesN')
    lv, sv = 100 * lu / lu_n, 100 * sp / sp_n
    ax.barh(y, -lv, height=0.56, color=BLUE, lw=0, zorder=2)
    ax.barh(y, sv, height=0.56, color=ORANGE, lw=0, zorder=2)
    ax.text(-lv - 1.5, y, f'{lu} of {lu_n} {plural}', ha='right', va='center', fontsize=6.6, color=INK)
    if sv > 45:                                     # long bar: label inside it
        ax.text(sv - 1.5, y, f'{sp} of {sp_n} species', ha='right', va='center', fontsize=6.6, color='white',
                weight='medium')
    else:
        ax.text(sv + 1.5, y, f'{sp} of {sp_n} species', ha='left', va='center', fontsize=6.6, color=INK)
    ax.text(-75, y, label, ha='right', va='center', fontsize=7.5, color=INK, weight='medium')
ax.axvline(0, color=INK2, lw=0.8, zorder=3)
ax.set_xlim(-72, 72)
ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_yticks([])
ticks = [-60, -40, -20, 0, 20, 40, 60]
ax.set_xticks(ticks, [f'{abs(t)}%' for t in ticks], fontsize=6.5)
ax.grid(True, axis='x', color=GRID, lw=0.6, zorder=0)
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color(GRID)
ax.text(-36, len(rows) - 0.35, 'lumps: names that hold more than one species', ha='center', va='bottom', fontsize=7,
        color=BLUE, weight='medium')
ax.text(36, len(rows) - 0.35, 'splits: species spread over more than one name', ha='center', va='bottom',
        fontsize=7, color=ORANGE, weight='medium')
fig.text(0.02, 0.975, '(a) Each system errs in its own direction', fontsize=8.5, weight='bold', color=INK, va='top')
fig.text(0.02, 0.945, f'{M["XwNamedCultures"]} named cultures of {M["XwNamedSpecies"]} species, each counted by '
         'its most common name', fontsize=6.8, color=INK2, va='top')


# ---------------------------------------------------------------- (b) two genera, culture by culture
def short(sp):
    g, _, e = sp.partition(' ')
    return f'{g[0]}. {e}' if e else sp


cul = list(csv.DictReader(open('analysis/xw/crosswalk_v4final.cultures.tsv'), delimiter='\t'))
named = [r for r in cul if r['name_status'] == 'ok' and r['named'] == 'True']
COLS = [('VT', 'VT (SSU)'), ('LSU', 'LSU name'), ('SH', 'SH (ITS)')]


def node_name(r, arm):
    v = r[f'{arm}_modal']
    if not v:
        return 'no name'
    return v.split('.')[0] if arm == 'SH' else (short(v) if arm == 'LSU' else v)


def alluvial(ax, genus, title):
    rs = [r for r in named if r['genus'] == genus]
    species = sorted({r['organism'] for r in rs})
    color = {s: SPECIES[i] for i, s in enumerate(species)}
    rs.sort(key=lambda r: (species.index(r['organism']), node_name(r, 'LSU'), node_name(r, 'SH')))
    gap, x = 0.9, [0.0, 1.0, 2.0]
    w = 0.045
    # node order per column: by the mean position of their cultures in the sorted list (fewer crossings)
    pos = {}
    for c, (arm, _) in enumerate(COLS):
        names = {}
        for i, r in enumerate(rs):
            names.setdefault(node_name(r, arm), []).append(i)
        order = sorted(names, key=lambda n: (n == 'no name', sum(names[n]) / len(names[n])))
        y = 0.0
        for n in order:
            k = len(names[n])
            pos[(c, n)] = [y, y + k, y]                     # top, bottom, next free slot
            y += k + gap
        total = y - gap
        for n in order:                                     # draw nodes
            t, b, _ = pos[(c, n)]
            ax.add_patch(Rectangle((x[c] - w / 2, t), w, b - t, color=NONE if n == 'no name' else INK2, lw=0,
                                   zorder=3))
            ital = arm == 'LSU' and n != 'no name'
            if c == 0:
                ax.text(x[c] - w, (t + b) / 2, n, ha='right', va='center', fontsize=6.3, color=INK)
            elif c == 2:
                ax.text(x[c] + w, (t + b) / 2, n, ha='left', va='center', fontsize=5.8, color=INK2 if n == 'no name' else INK)
            else:
                ax.text(x[c] + w * 0.9, t + 0.05, n, ha='left', va='top', fontsize=6.0, color=INK,
                        style='italic' if ital else 'normal',
                        bbox=dict(boxstyle='round,pad=0.12', fc=SURFACE, ec='none', alpha=0.85), zorder=4)
        ax.text(x[c], -0.6, COLS[c][1], ha='center', va='bottom', fontsize=6.8, color=INK2, weight='medium')
    # ribbons: one per culture, VT -> LSU -> SH
    slot = {}
    for r in rs:
        ys = []
        for c, (arm, _) in enumerate(COLS):
            key = (c, node_name(r, arm))
            y0 = pos[key][2]
            pos[key][2] += 1
            ys.append(y0)
        for c in range(2):
            y0, y1 = ys[c], ys[c + 1]
            x0, x1 = x[c] + w / 2, x[c + 1] - w / 2
            xm = (x0 + x1) / 2
            verts = [(x0, y0), (xm, y0), (xm, y1), (x1, y1), (x1, y1 + 1), (xm, y1 + 1), (xm, y0 + 1), (x0, y0 + 1),
                     (x0, y0)]
            codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4,
                     Path.CURVE4, Path.CLOSEPOLY]
            ax.add_patch(PathPatch(Path(verts, codes), fc=color[r['organism']], ec=SURFACE, lw=0.35, alpha=0.85,
                                   zorder=2))
        slot[r['culture']] = ys
    height = max(v[1] for v in pos.values())
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(height + 0.5, -2.2)
    ax.axis('off')
    ax.set_title(title, fontsize=8, color=INK, loc='left', pad=4)
    return species, color, height


axg = fig.add_axes([0.07, 0.085, 0.39, 0.475])
axf = fig.add_axes([0.57, 0.085, 0.39, 0.475])
sg, cg, hg = alluvial(axg, 'Glomus', r'$\it{Glomus}$: one VT, four species')
sf, cf, hf = alluvial(axf, 'Funneliformis', r'$\it{Funneliformis}$: one species, six SHs')
for ax, sps, col in ((axg, sg, cg), (axf, sf, cf)):
    hs = [Patch(color=col[s], label=short(s)) for s in sps]
    leg = ax.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=4, frameon=False, fontsize=6.4,
                    handlelength=1.0, handleheight=0.8, columnspacing=1.2, labelcolor=INK)
    for t in leg.get_texts():
        t.set_style('italic')
fig.text(0.02, 0.648, '(b) Two genera, culture by culture: VT, then LSU name, then SH', fontsize=8.5, weight='bold',
         color=INK, va='top')
fig.text(0.02, 0.620, 'Each ribbon is one culture, coloured by its species; gray nodes: no name in that system',
         fontsize=6.8, color=INK2, va='top')
for ext in ('pdf', 'png'):
    fig.savefig(f'figures/fig_nesting.{ext}', dpi=220)
print('wrote figures/fig_nesting.pdf and .png')
# the counts the caption cites
gl = [r for r in named if r['genus'] == 'Glomus']
mo = [r for r in named if r['organism'] == 'Funneliformis mosseae']
print(f"Glomus: {len(gl)} cultures, {len({r['organism'] for r in gl})} species, VTs {sorted({r['VT_modal'] for r in gl})}, "
      f"no SH {sum(not r['SH_modal'] for r in gl)}")
print(f"F. mosseae: {len(mo)} cultures; LSU name F. coronatus {sum(r['LSU_modal'] == 'Funneliformis coronatus' for r in mo)}, "
      f"no LSU name {sum(not r['LSU_modal'] for r in mo)}; SHs {len({r['SH_modal'] for r in mo if r['SH_modal']})}")
vt65 = sorted({r['organism'] for r in named if r['VT_modal'] == 'VTX00065'})
print('VTX00065 holds', vt65)
