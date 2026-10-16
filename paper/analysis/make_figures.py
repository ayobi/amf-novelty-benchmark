#!/usr/bin/env python3
"""Figure 2: species called known although no reference culture holds them ("false known"), by method.

Every number is read from numbers_generated.tex (written by make_numbers.py from the run outputs), so the figure
cannot drift from the text.

  python3 analysis/make_figures.py        (writes figures/fig_falseknown.pdf and .png)
"""
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

M = dict(re.findall(r'\\newcommand\{\\(\w+)\}\{([^{}]*)\}', open('numbers_generated.tex').read()))
num = lambda k: int(M[k])  # noqa: E731

INK, INK2 = '#0b0b0b', '#52514e'
BLUE, ORANGE = '#2a78d6', '#eb6834'          # validated categorical slots 1-2 (light)
plt.rcParams.update({'font.size': 7.5, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.grid': True, 'axes.grid.axis': 'y', 'grid.color': '#e6e5e0', 'grid.linewidth': 0.6,
                     'hatch.color': '#fcfcfb', 'hatch.linewidth': 0.8,
                     'mathtext.fontset': 'custom', 'mathtext.it': 'Inter:italic', 'mathtext.rm': 'Inter'})

fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0), gridspec_kw={'width_ratios': [1.0, 1.0, 1.25]})


def bars(ax, labels, vals, dens, cols, hatches, ymax):
    xs = range(len(vals))
    for x, v, d, c, h in zip(xs, vals, dens, cols, hatches):
        ax.bar(x, 100 * v / d, width=0.62, color=c, hatch=h, edgecolor='#fcfcfb', lw=0, zorder=2)
        ax.text(x, 100 * v / d + ymax * 0.015, f'{v}/{d}\n({100 * v / d:.1f}%)', ha='center', va='bottom',
                fontsize=6.5, color=INK)
    ax.set_xticks(list(xs), labels, fontsize=6.6)
    ax.set_ylim(0, ymax)
    ax.tick_params(axis='x', length=0)


# (a) species withheld from the cultures; the LSU reference intact
a = axes[0]
bars(a, ['nearest LSU\nreference', 'LSU species\nclade', 'clade +\nITS rule'],
     [num('FKBase'), num('FKClade'), num('FKITS')], [num('HeldN')] * 3, [BLUE] * 3, ['', '', '////'], 34)
a.set_ylabel('Assigned to a known species (% of sequences)')
a.set_xlabel('LSU reference intact', fontsize=7, color=INK2)
a.set_title('(a) Species withheld from\nthe cultures', fontsize=8, loc='left')

# (b) Glomus rugosae, also removed from the LSU reference
b = axes[1]
n = num('RugN')
bars(b, ['intact', 'removed', 'removed +\nITS rule'],
     [num('RugFull'), num('RugRemoved'), num('RugRemovedITS')], [n] * 3, [BLUE, ORANGE, ORANGE], ['', '', '////'], 115)
b.set_xlabel('LSU reference', fontsize=7, color=INK2)
b.set_title('(b) $\\it{Glomus\\ rugosae}$, also\nremoved from the LSU reference', fontsize=8, loc='left')

# (c) the 34 test species, pooled
c = axes[2]
h = num('MfourHeldN')
bars(c, ['intact', 'intact +\nITS rule', 'removed', 'removed +\nITS rule'],
     [num('MfourControlFKvOne'), num('MfourControlFKvTwo'), num('MfourFKvOne'), num('MfourFKvTwo')], [h] * 4,
     [BLUE, BLUE, ORANGE, ORANGE], ['', '////', '', '////'], 17)
c.set_xlabel('LSU reference', fontsize=7, color=INK2)
c.set_title(f"(c) {M['MfourFolds']} test species, each also\nremoved from the LSU reference", fontsize=8, loc='left')

fig.tight_layout(w_pad=1.2)
for ext in ('pdf', 'png'):
    fig.savefig(f'figures/fig_falseknown.{ext}', dpi=220)
print('figures/fig_falseknown.pdf')
