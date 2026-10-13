#!/usr/bin/env python3
"""Supplementary files and the tree figure from the copy trees (benchmarks/copy_checks/trees).

  supp    Supplementary Data S1: for every order and region, the alignment (FASTA) and the tree (Newick, branch
          supports SH-aLRT/UFBoot) with tips relabelled `accession|species|culture|class|copies`, plus the member
          table (which copies each tip stands for) and a README.
  figure  The Rhizophagus clade of the Glomerales linked tree, rooted on a Funneliformis copy: tips coloured by
          species, short-5.8S tips marked, clades made only of short-5.8S copies shaded, strongly supported
          branches (SH-aLRT >= 80 and UFBoot >= 95) dotted.
  summary The same clade with clades of main-type copies of one species, and clades of short-5.8S copies only,
          collapsed into triangles (the main-text figure).

  python3 stage3/copy_tree_outputs.py supp   --ref reference/culture_ref_v4 --trees benchmarks/copy_checks/trees \\
      --out benchmarks/copy_checks/supplementary_data_S1
  python3 stage3/copy_tree_outputs.py figure --ref reference/culture_ref_v4 --trees benchmarks/copy_checks/trees \\
      --out benchmarks/copy_checks/figures
"""
import argparse
import csv
import os
import re
import shutil
import sys

from Bio import Phylo

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arms import read_tsv  # noqa: E402

ORDERS = ['Archaeosporales', 'Paraglomerales', 'Diversisporales', 'Glomerales']
TREES = [('ssu_flank', 'SSU flank'), ('full_its', 'full ITS'), ('lsu_flank', 'LSU flank'),
         ('linked', 'linked copy (SSU flank + full ITS + LSU flank, one partition each)')]
TAXA = [('Rhizophagus irregularis', 'irregularis'), ('[Rhizoglomus] vesiculiferum', 'vesiculiferum'),
        ('Rhizophagus intraradices', 'intraradices'), ('Rhizophagus prolifer', 'prolifer'),
        ('Rhizophagus silesianus', 'silesianus'), ('Rhizophagus sp. (in: glomeromycetes)', 'sp.')]
COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']   # validated slots 1-6 (light)


def members(path):
    return {r['tip']: r['accessions'].split(',') for r in read_tsv(path)}


def klass(r):
    return {'short_5.8S': 'short5.8S', 'same_length': 'samelength'}.get(r['div_kind'], 'main')


def safe(s):
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s).strip('_')


def cmd_supp(a):
    ref = {r['accession']: r for r in read_tsv(os.path.join(a.ref, 'copies.tsv')) if r['status'] == 'retained'}
    os.makedirs(a.out, exist_ok=True)
    n = 0
    for order in ORDERS:
        d = os.path.join(a.out, order)
        os.makedirs(d, exist_ok=True)
        for t, _ in TREES:
            src = os.path.join(a.trees, order, t)
            mem = members(src + '.members.tsv')
            lab = {}
            for tip, accs in mem.items():
                rs = [ref[x] for x in accs]
                orgs = sorted({r['organism'] for r in rs})
                cults = sorted({r['culture_label'] for r in rs})
                classes = sorted({klass(r) for r in rs})
                lab[tip] = '|'.join([tip, safe(orgs[0] if len(orgs) == 1 else f'{len(orgs)}_taxa'),
                                     safe(cults[0] if len(cults) == 1 else f'{len(cults)}_cultures'),
                                     '+'.join(classes), f'n{len(accs)}'])
            txt = open(src + '.treefile').read()
            txt = re.sub(r'([(,])([^(),:;]+):', lambda m: m.group(1) + "'" + lab.get(m.group(2), m.group(2)) + "':", txt)
            with open(os.path.join(d, f'{t}.tree.nwk'), 'w') as fh:
                fh.write(txt)
            shutil.copy(src + '.aln.fasta', os.path.join(d, f'{t}.alignment.fasta'))
            with open(os.path.join(d, f'{t}.members.tsv'), 'w', newline='') as fh:
                w = csv.writer(fh, delimiter='\t', lineterminator='\n')
                w.writerow(['tip', 'label', 'accession', 'species', 'culture', 'copy_class'])
                for tip, accs in sorted(mem.items()):
                    for x in accs:
                        r = ref[x]
                        w.writerow([tip, lab[tip], x, r['organism'], r['culture'], klass(r)])
            if t == 'linked':
                shutil.copy(src + '.partitions.nex', os.path.join(d, 'linked.partitions.nex'))
            n += 1
    with open(os.path.join(a.out, 'README.md'), 'w') as fh:
        fh.write('# Supplementary Data S1: alignments and trees of the culture copies\n\n'
                 'One folder per order. For each region (SSU flank, full ITS, LSU flank) and the linked copy:\n\n'
                 '- `<region>.alignment.fasta`: MAFFT L-INS-i alignment of the unique sequences of the order.\n'
                 '- `<region>.tree.nwk`: IQ-TREE maximum-likelihood tree (ModelFinder; branch labels SH-aLRT/UFBoot '
                 'from 1,000 replicates each). Tip label: `accession|species|culture|copy class|n copies`, where the '
                 'accession is the first GenBank accession carrying that exact sequence.\n'
                 '- `<region>.members.tsv`: every copy each tip stands for (identical sequences share one tip).\n'
                 '- `linked.partitions.nex`: the region partitions of the linked alignment.\n\n'
                 'Copy classes: main type, short5.8S, samelength (the paper\'s divergent classes). Built by '
                 '`stage3/copy_trees.py` under `benchmarks/copy_checks/PROTOCOL_COPY_CHECKS.md`; relabelled by '
                 '`stage3/copy_tree_outputs.py`.\n')
    print(f'{n} alignments and trees written to {a.out}')


def cmd_figure(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches  # noqa: F401
    from matplotlib.lines import Line2D
    plt.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.it': 'Inter:italic', 'mathtext.rm': 'Inter'})
    ref = {r['accession']: r for r in read_tsv(os.path.join(a.ref, 'copies.tsv')) if r['status'] == 'retained'}
    src = os.path.join(a.trees, 'Glomerales', 'linked')
    mem = members(src + '.members.tsv')
    tree = Phylo.read(src + '.treefile', 'newick')
    info = lambda n: [ref[x] for x in mem[n]]  # noqa: E731
    out = next(tp for tp in tree.get_terminals() if all(r['genus'] == 'Funneliformis' for r in info(tp.name)))
    tree.root_with_outgroup(out)
    rh = [tp for tp in tree.get_terminals() if any(r['genus'] in ('Rhizophagus', 'Rhizoglomus') for r in info(tp.name))]
    clade = tree.common_ancestor(rh)
    tips = clade.get_terminals()
    isshort = lambda n: any(r['div_kind'] == 'short_5.8S' for r in info(n))  # noqa: E731
    # coordinates
    depth = {}

    def setx(cl, x):
        depth[cl] = x
        for ch in cl.clades:
            setx(ch, x + (ch.branch_length or 0))
    setx(clade, 0.0)
    ypos = {tp: i for i, tp in enumerate(tips)}

    def sety(cl):
        if cl.is_terminal():
            return ypos[cl]
        ys = [sety(ch) for ch in cl.clades]
        ypos[cl] = (min(ys) + max(ys)) / 2
        return ypos[cl]
    sety(clade)
    xmax = max(depth[tp] for tp in tips)
    fig, ax = plt.subplots(figsize=(6.3, 9.3))            # one A4 page at 100% (Figure S3)
    # shade clades made only of short-5.8S copies
    pure = []

    def walk(cl):
        if all(isshort(tp.name) for tp in cl.get_terminals()):
            pure.append(cl)
            return
        for ch in cl.clades:
            walk(ch)
    walk(clade)
    for cl in pure:
        ts = cl.get_terminals()
        if len(ts) < 2:
            continue
        y0, y1 = min(ypos[t] for t in ts) - 0.5, max(ypos[t] for t in ts) + 0.5
        lo, hi = -xmax * 0.02, xmax * 1.55
        ax.axhspan(y0, y1, xmin=(depth[cl] - xmax * 0.01 - lo) / (hi - lo), xmax=1, color='#ecebe6', zorder=0, lw=0)
    for cl in clade.find_clades():
        x, y = depth[cl], ypos[cl]
        for ch in cl.clades:
            ax.plot([x, depth[ch]], [ypos[ch], ypos[ch]], color='#52514e', lw=0.7, zorder=2)
        if cl.clades:
            ax.plot([x, x], [min(ypos[c] for c in cl.clades), max(ypos[c] for c in cl.clades)],
                    color='#52514e', lw=0.7, zorder=2)
            lab = cl.name or (str(cl.confidence) if cl.confidence is not None else '')
            m = re.match(r'([\d.]+)/([\d.]+)', lab or '')
            if m and float(m.group(1)) >= 80 and float(m.group(2)) >= 95 and cl is not clade:
                ax.plot(x, y, 'o', ms=2.2, color='#0b0b0b', zorder=3)
    col = {t: c for (t, _), c in zip(TAXA, COLORS)}
    short = {t: s for t, s in TAXA}
    labels = []
    for tp in tips:
        rs = info(tp.name)
        orgs = sorted({r['organism'] for r in rs})
        cults = sorted({r['culture_label'] for r in rs})
        x, y = depth[tp], ypos[tp]
        c = col.get(orgs[0], '#52514e') if len(orgs) == 1 else '#52514e'
        sh = isshort(tp.name)
        sl = any(r['div_kind'] == 'same_length' for r in rs)
        ax.plot(x, y, '^' if sh else ('s' if sl else 'o'), ms=3.6 if (sh or sl) else 3.0, color=c,
                mec='#fcfcfb', mew=0.4, zorder=4)
        # label segments: species epithets in italics ('sp.' upright), then culture, copy count and class
        segs = []
        for k, o in enumerate(orgs):
            e = short.get(o, o.split()[-1])
            if k:
                segs.append((' + ', False))
            segs.append((e, e != 'sp.'))
        cu = cults[0] if len(cults) == 1 else f'{cults[0]} +{len(cults) - 1}'
        rest = f'  {cu}' + (f'  ×{len(rs)}' if len(rs) > 1 else '') + ('  short 5.8S' if sh else ('  same length' if sl else ''))
        segs.append((rest, False))
        labels.append((x + xmax * 0.012, y, segs, 'bold' if (sh or sl) else 'normal'))
    ax.set_ylim(len(tips) - 0.5, -0.8)
    ax.set_xlim(-xmax * 0.02, xmax * 1.55)
    ax.axis('off')
    ax.plot([0, 0.01], [len(tips) + 0.2] * 2, color='#0b0b0b', lw=0.8, clip_on=False)
    ax.text(0.005, len(tips) + 1.0, '0.01 substitutions/site', fontsize=5.5, ha='center', va='top')
    fig.tight_layout()
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    from matplotlib.transforms import offset_copy
    for x, y, segs, weight in labels:
        dx = 0.0
        for txt, it in segs:
            t = ax.text(x, y, txt, fontsize=5.2, va='center', color='#0b0b0b', fontweight=weight,
                        fontstyle='italic' if it else 'normal',
                        transform=offset_copy(ax.transData, fig=fig, x=dx, y=0, units='points'))
            dx += t.get_window_extent(renderer=rend).width * 72 / fig.dpi
    def disp(t):
        t = t.replace(' (in: glomeromycetes)', '')
        return f'$\\mathit{{{t[:-4]}}}$ sp.' if t.endswith(' sp.') else t
    handles = [Line2D([], [], marker='o', ls='', color=c, ms=4, label=disp(t)) for (t, _), c in zip(TAXA, COLORS)]
    handles += [Line2D([], [], marker='^', ls='', color='#52514e', ms=4, label='short-5.8S copy (bold label)'),
                Line2D([], [], marker='s', ls='', color='#52514e', ms=3.5, label='same-length divergent copy'),
                Line2D([], [], marker='o', ls='', color='#0b0b0b', ms=2.5,
                       label='SH-aLRT ≥ 80 and UFBoot ≥ 95'),
                matplotlib.patches.Patch(color='#ecebe6', label='clade of short-5.8S copies only')]
    leg = ax.legend(handles=handles, frameon=False, fontsize=6, loc='upper left', bbox_to_anchor=(0.0, 1.0))
    for t in leg.get_texts()[:len(TAXA)]:               # species names in italics ('sp.' upright, via mathtext)
        if '$' not in t.get_text():
            t.set_fontstyle('italic')
    os.makedirs(a.out, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(a.out, f'fig_rhizophagus_linked_tree.{ext}'), dpi=220)
    print(f'{len(tips)} tips; short-only clades: {[len(c.get_terminals()) for c in pure]}')


def cmd_summary(a):
    """The same Rhizophagus clade with clades collapsed: main-type copies of one species, and clades made only of
    short-5.8S copies, each drawn as one triangle (for the main text; the full tree is the 'figure' output)."""
    import math
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon
    from matplotlib.transforms import offset_copy
    plt.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.it': 'Inter:italic', 'mathtext.rm': 'Inter'})
    ref = {r['accession']: r for r in read_tsv(os.path.join(a.ref, 'copies.tsv')) if r['status'] == 'retained'}
    src = os.path.join(a.trees, 'Glomerales', 'linked')
    mem = members(src + '.members.tsv')
    tree = Phylo.read(src + '.treefile', 'newick')
    info = lambda n: [ref[x] for x in mem[n]]  # noqa: E731
    out = next(tp for tp in tree.get_terminals() if all(r['genus'] == 'Funneliformis' for r in info(tp.name)))
    tree.root_with_outgroup(out)
    rh = [tp for tp in tree.get_terminals() if any(r['genus'] in ('Rhizophagus', 'Rhizoglomus') for r in info(tp.name))]
    clade = tree.common_ancestor(rh)
    isshort = lambda n: any(r['div_kind'] == 'short_5.8S' for r in info(n))  # noqa: E731

    def copies(cl):
        return [r for t in cl.get_terminals() for r in info(t.name)]

    def kind(cl):
        rs = copies(cl)
        if all(isshort(t.name) for t in cl.get_terminals()):
            return 'short'
        if all(r['copy_class'] == 'main' for r in rs) and len({r['organism'] for r in rs}) == 1:
            return 'main'
        return 'tip' if cl.is_terminal() else None

    depth = {}

    def setx(cl, x):
        depth[cl] = x
        for ch in cl.clades:
            setx(ch, x + (ch.branch_length or 0))
    setx(clade, 0.0)
    leaves = []

    def walk(cl):
        k = kind(cl)
        if k:
            leaves.append((cl, k))
            return
        for ch in cl.clades:
            walk(ch)
    walk(clade)
    ypos = {cl: i for i, (cl, _) in enumerate(leaves)}
    isleaf = {cl for cl, _ in leaves}

    def sety(cl):
        if cl in isleaf:
            return ypos[cl]
        ys = [sety(ch) for ch in cl.clades]
        ypos[cl] = (min(ys) + max(ys)) / 2
        return ypos[cl]
    sety(clade)
    xmax = max(depth[t] for t in clade.get_terminals())
    col = {t: c for (t, _), c in zip(TAXA, COLORS)}
    short = {t: s for t, s in TAXA}
    fig, ax = plt.subplots(figsize=(7.2, 0.27 * len(leaves) + 1.3))

    def strong(cl):
        lab = cl.name or (str(cl.confidence) if cl.confidence is not None else '')
        m = re.match(r'([\d.]+)/([\d.]+)', lab or '')
        return bool(m and float(m.group(1)) >= 80 and float(m.group(2)) >= 95)

    def draw(cl):
        x, y = depth[cl], ypos[cl]
        if cl in isleaf:
            return
        ys = [ypos[ch] for ch in cl.clades]
        ax.plot([x, x], [min(ys), max(ys)], color='#52514e', lw=0.8, zorder=2)
        for ch in cl.clades:
            ax.plot([x, depth[ch]], [ypos[ch], ypos[ch]], color='#52514e', lw=0.8, zorder=2)
            draw(ch)
        if strong(cl) and cl is not clade:
            ax.plot(x, y, 'o', ms=3.0, color='#0b0b0b', zorder=5)
    draw(clade)
    labels = []
    for cl, k in leaves:
        x, y = depth[cl], ypos[cl]
        rs = copies(cl)
        orgs = sorted({r['organism'] for r in rs})
        ncul = len({r['culture'] for r in rs})
        tips = cl.get_terminals()
        if k == 'short':
            fc, ec = '#b4b3ad', '#52514e'
        else:
            fc = col.get(orgs[0], '#b4b3ad') if len(orgs) == 1 else '#b4b3ad'
            ec = fc
        xe = max(max(depth[t] for t in tips), x + xmax * 0.012)
        if len(tips) > 1:
            h = min(0.82, 0.30 + 0.11 * math.log2(len(rs)))
            ax.add_patch(Polygon([(x, y), (xe, y - h / 2), (xe, y + h / 2)], closed=True, fc=fc, ec=ec, lw=0.6,
                                 alpha=0.9 if k == 'short' else 0.75, zorder=3))
            if strong(cl):
                ax.plot(x, y, 'o', ms=3.0, color='#0b0b0b', zorder=5)
        else:
            mk = '^' if k == 'short' else ('s' if any(r['div_kind'] == 'same_length' for r in rs) else 'o')
            ax.plot(x, y, mk, ms=4.2, color=fc if k != 'short' else '#52514e', mec='#fcfcfb', mew=0.4, zorder=4)
        segs = []
        for i, o in enumerate(orgs):
            e = short.get(o, o.split()[-1])
            if i:
                segs.append((' + ', False))
            segs.append((e, e != 'sp.'))
        nmain = sum(r['copy_class'] == 'main' for r in rs)
        if k == 'short':
            what = 'short 5.8S'
            extra = '\u2020' if nmain else ''
        elif any(r['div_kind'] == 'same_length' for r in rs):
            what, extra = 'same-length copy', ''
        else:
            what, extra = 'main type', ''
        segs.append((f'  {what}: {len(rs)} cop{"y" if len(rs) == 1 else "ies"}{extra}, '
                     f'{ncul} culture{"" if ncul == 1 else "s"}', False))
        labels.append((xe + xmax * 0.015, y, segs, 'bold' if k == 'short' else 'normal'))
    ax.set_ylim(len(leaves) - 0.4, -0.9)
    ax.set_xlim(-xmax * 0.03, xmax * 1.95)
    ax.axis('off')
    ax.plot([0, 0.01], [len(leaves) - 0.1] * 2, color='#0b0b0b', lw=0.8, clip_on=False)
    ax.text(0.005, len(leaves) + 0.25, '0.01 substitutions/site', fontsize=6.5, ha='center', va='top')
    fig.tight_layout()
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for x, y, segs, weight in labels:
        dx = 0.0
        for txt, it in segs:
            t = ax.text(x, y, txt, fontsize=7, va='center', color='#0b0b0b', fontweight=weight,
                        fontstyle='italic' if it else 'normal',
                        transform=offset_copy(ax.transData, fig=fig, x=dx, y=0, units='points'))
            dx += t.get_window_extent(renderer=rend).width * 72 / fig.dpi

    def disp(t):
        t = t.replace(' (in: glomeromycetes)', '')
        return f'$\\mathit{{{t[:-4]}}}$ sp.' if t.endswith(' sp.') else t
    handles = [Line2D([], [], marker='>', ls='', color=c, ms=5, alpha=0.75, label=disp(t))
               for (t, _), c in zip(TAXA, COLORS)]
    handles += [Line2D([], [], marker='>', ls='', color='#b4b3ad', mec='#52514e', ms=5,
                       label='clade of short-5.8S copies only (bold)'),
                Line2D([], [], marker='o', ls='', color='#0b0b0b', ms=3, label='SH-aLRT \u2265 80 and UFBoot \u2265 95')]
    leg = ax.legend(handles=handles, frameon=False, fontsize=6.5, loc='lower right', bbox_to_anchor=(1.0, 0.0))
    for t in leg.get_texts()[:len(TAXA)]:
        if '$' not in t.get_text():
            t.set_fontstyle('italic')
    os.makedirs(a.out, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(a.out, f'fig_rhizophagus_summary_tree.{ext}'), dpi=220)
    print(f'{len(leaves)} leaves ({sum(k == "short" for _, k in leaves)} short-only clades)')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('supp', 'figure', 'summary'):
        p = sub.add_parser(name)
        p.add_argument('--ref', required=True)
        p.add_argument('--trees', required=True)
        p.add_argument('--out', required=True)
    a = ap.parse_args()
    {'supp': cmd_supp, 'figure': cmd_figure, 'summary': cmd_summary}[a.cmd](a)


if __name__ == '__main__':
    main()
