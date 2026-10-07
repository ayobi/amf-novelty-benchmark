#!/usr/bin/env python3
"""Figure 1: false-known calls by configuration. Numbers are read from the
benchmark summaries listed in numbers.tex (typed by hand in draft 0.1)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

cfg = ['Nearest-reference\nLSU', 'LSU species-clade\nplacement', 'Placement +\nITS rule (N6)']
false_known = [215, 71, 18]
n = 803
fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={'width_ratios': [1.25, 1]})
bars = a.bar(cfg, [100 * v / n for v in false_known], color=['#9a9a9a', '#6a8caf', '#2f5d8a'])
for r, v in zip(bars, false_known):
    a.text(r.get_x() + r.get_width() / 2, r.get_height() + 0.6, f'{v}/{n}\n({100 * v / n:.1f}%)', ha='center', va='bottom', fontsize=8)
a.set_ylabel('False known, % of held-out copies')
a.set_ylim(0, 34)
a.set_title('(a) Culture bridge held out; V18 intact', fontsize=9, loc='left')
labels = ['G. rugosae\nV18 intact', 'G. rugosae\nV18 removed', 'G. rugosae\nremoved + N6']
vals = [4 / 75, 75 / 75, 62 / 75]
bars = b.bar(labels, [100 * v for v in vals], color=['#6a8caf', '#c0504d', '#d98c8a'])
for r, num in zip(bars, ['4/75', '75/75', '62/75']):
    b.text(r.get_x() + r.get_width() / 2, r.get_height() + 1.5, num, ha='center', va='bottom', fontsize=8)
b.set_ylim(0, 112)
b.set_title('(b) Target also removed from V18', fontsize=9, loc='left')
for ax in (a, b):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='x', labelsize=7.5)
plt.tight_layout()
plt.savefig('figures/fig1_falseknown.pdf')
print('figures/fig1_falseknown.pdf')
