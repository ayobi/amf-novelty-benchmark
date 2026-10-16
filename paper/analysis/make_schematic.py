#!/usr/bin/env python3
"""Figure 1: the linked rRNA gene sequence with the three barcode windows, and how a sequence is named.

Coordinates are the medians over the 903 retained copies of culture reference v4 (copy 2,758 nt): SSU 1-1473,
ITS1 1474-1580, 5.8S 1581-1737, ITS2 1738-1939, LSU 1940-2758; NS31-AM1 232-739, LROR-FLR2 1981-2712.
Window lengths are the medians of the primer-trimmed windows (NS31-AM1 508 nt, full ITS 466 nt, LROR-FLR2 735 nt).

  python3 analysis/make_schematic.py      (writes figures/fig_schematic.pdf and .png)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

INK, INK2, SURF = '#0b0b0b', '#52514e', '#fcfcfb'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'      # validated categorical slots 1-3 (light)
GENE, GENE2 = '#dcdbd6', '#b4b3ad'
REG = [('SSU (18S)', 1, 1473, GENE), ('ITS1', 1474, 1580, 'white'), ('5.8S', 1581, 1737, GENE2),
       ('ITS2', 1738, 1939, 'white'), ('LSU (28S, 5′ part)', 1940, 2758, GENE)]
WIN = [('NS31–AM1', 232, 739, BLUE, 'MaarjAM', 'virtual taxon (VT)', 510),
       ('full ITS', 1474, 1939, ORANGE, 'UNITE', 'species hypothesis (SH)', 470),
       ('LR0R–FLR2', 1981, 2712, AQUA, 'LSU reference tree', 'species clade', 735)]

plt.rcParams.update({'font.size': 7, 'text.color': INK})
fig = plt.figure(figsize=(7.2, 4.87))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 7.2)
ax.set_ylim(-0.32, 4.55)
ax.axis('off')

# ---------------------------------------------------------------- (a) the copy and its windows
X0, X1, NT = 0.45, 6.95, 2758
sx = lambda nt: X0 + (X1 - X0) * (nt - 0.5) / NT  # noqa: E731
ax.text(0.1, 4.5, '(a) One rRNA gene sequence, three marker regions', fontsize=8.5, fontweight='bold', va='top')
yb, hb = 3.86, 0.24
for name, a, b, fc in REG:
    ax.add_patch(Rectangle((sx(a), yb), sx(b + 1) - sx(a), hb, fc=fc, ec=INK2, lw=0.7, zorder=2))
    if b - a > 600:
        ax.text((sx(a) + sx(b + 1)) / 2, yb + hb / 2, name, ha='center', va='center', fontsize=7, zorder=3)
    else:
        ax.text((sx(a) + sx(b + 1)) / 2, yb - 0.05, name, ha='center', va='top', fontsize=6.5)
# the deletion of the short-5.8S copies, straight above the 5.8S
xd = sx(1581 + 125)
ax.annotate('short-5.8S sequences: 3–6 nt deleted in a 5.8S hairpin', xy=(xd, yb + hb), xytext=(xd, 4.33),
            fontsize=6.5, color=INK, ha='center', va='center',
            arrowprops=dict(arrowstyle='-|>', color=INK2, lw=0.6, shrinkA=3, shrinkB=1, mutation_scale=6))
# scale bar
ax.plot([sx(1), sx(500)], [3.72, 3.72], color=INK, lw=0.8)
ax.text(sx(250), 3.67, '500 nt', ha='center', va='top', fontsize=6)
# windows
yw = 3.42
for w, a, b, col, db, unit, ln in WIN:
    x0, x1 = sx(a), sx(b + 1)
    ax.plot([x0, x0, x1, x1], [yw + 0.07, yw, yw, yw + 0.07], color=col, lw=1.6, solid_capstyle='butt')
    xm = (x0 + x1) / 2
    ax.text(xm, yw - 0.07, f'{w} (~{ln} nt)', ha='center', va='top', fontsize=7, fontweight='bold', color=INK)
    ax.text(xm, yw - 0.25, f'{db}\n{unit}', ha='center', va='top', fontsize=6.5, color=INK, linespacing=1.25)

# ---------------------------------------------------------------- (b) from three names to one call
ax.text(0.1, 2.62, '(b) From three names to one assignment', fontsize=8.5, fontweight='bold', va='top')


def box(x, y, w, h, text, ec=INK2, lw=0.7, fs=6.6, fc=SURF, bold_first=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.0,rounding_size=0.06', fc=fc, ec=ec, lw=lw,
                                zorder=2))
    if bold_first:
        head, body = text.split('\n', 1)
        ax.text(x + 0.09, y + h - 0.09, head, ha='left', va='top', fontsize=fs + 0.4, fontweight='bold', zorder=3)
        ax.text(x + 0.09, y + h - 0.29, body, ha='left', va='top', fontsize=fs, linespacing=1.3, zorder=3)
    else:
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs, linespacing=1.3, zorder=3)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=8, color=INK2, lw=0.8,
                                 shrinkA=0, shrinkB=0, zorder=1))


YM = 1.12                                       # the row the arrows run along
# names from the three databases
for (w, a, b, col, db, unit, ln), y in zip(WIN, [1.77, 1.07, 0.37]):
    nm = {'MaarjAM': 'VT name', 'UNITE': 'SH name', 'LSU reference tree': 'LSU species clade'}[db]
    sub = {'MaarjAM': 'from the SSU', 'UNITE': 'from the ITS', 'LSU reference tree': 'from the LSU'}[db]
    box(0.12, y, 1.25, 0.5, f'{nm}\n{sub}', ec=col, lw=1.4)
    arrow(1.37, y + 0.25, 1.62, YM)
# the culture bridge
box(1.62, 0.5, 1.62, 1.24, 'Link to cultures\nEach name is linked to the\nspecies and genera of the\nreference cultures that have\nit; the assignment keeps what\nall linked names share',
    bold_first=True)
arrow(3.24, YM, 3.42, YM)
# the checks
box(3.42, 0.37, 1.9, 1.5, 'Checks, in order\n1  SSU discordant? SSU at least half\n    as divergent as ITS\n2  Divergent sequence? short 5.8S\n3  Just below a VT threshold?\n4  New? a VT or LSU name that no\n    culture has; no species shared;\n    ITS farther from every culture than\n    the within-species limit',
    bold_first=True)
arrow(5.32, YM, 5.50, YM)
# outcomes
out = [('Known species', 'the names point to one species'), ('Species group', 'several species the\nnames cannot split'),
       ('New species, known genus', 'genus known, species not'), ('New lineage', 'no name is found in a culture')]
yo = 2.08
for head, sub in out:
    h = 0.44 if '\n' not in sub else 0.56
    yo -= h
    ax.add_patch(FancyBboxPatch((5.50, yo), 1.6, h - 0.06, boxstyle='round,pad=0.0,rounding_size=0.06', fc=SURF,
                                ec=INK2, lw=0.7, zorder=2))
    ax.text(5.58, yo + h - 0.15, head, ha='left', va='center', fontsize=6.8, fontweight='bold')
    ax.text(5.58, yo + h - 0.26, sub, ha='left', va='top', fontsize=6.3, color=INK2, linespacing=1.25)
ax.text(5.50, yo - 0.06, 'Flags: SSU discordance, divergent\nsequence, near a VT threshold,\nnames that disagree on genus',
        ha='left', va='top', fontsize=6.0, color=INK2, linespacing=1.25)

for ext in ('pdf', 'png'):
    fig.savefig(f'figures/fig_schematic.{ext}', dpi=220, facecolor='white')
print('figures/fig_schematic.pdf')
