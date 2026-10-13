#!/usr/bin/env python3
"""Numbers and tables for the copy checks (benchmarks/copy_checks).

Reads the outputs of PROTOCOL_COPY_CHECKS.md (trees/, artefacts/) and of ADDENDUM1_FUNCTION.md (function/), the
exploratory threshold outputs (thresholds/, thresholds_long/), the paralog effects (paralog_effects/) and the
Glomerales linked tree, and writes numbers_copychecks.tex, tables/table_monophyly.tex, tables/table_artefacts.tex
and tables/table_clusters.tex. The shared-deletion counts are a post hoc reading, not
part of the locked protocol, and the text labels them so.

  python3 analysis/copy_checks_numbers.py --checks ../benchmarks/copy_checks \\
      --ref ../reference/culture_ref_v4 --paper .
"""
import argparse
import collections
import csv
import json
import os
import re
import statistics


def rd(p):
    with open(p, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def lens(path):
    out, name, n = {}, None, 0
    for line in open(path):
        if line.startswith('>'):
            if name:
                out[name] = n
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
    if name:
        out[name] = n
    return out


def emph(name):
    name = name.replace('[', '{[}').replace(']', '{]}')
    if name.endswith(' sp.'):
        return '\\spn{' + name[:-4] + '} sp.'
    return '\\spn{' + name + '}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--checks', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--paper', required=True)
    a = ap.parse_args()
    m = []

    def add(k, v, c=None):
        if c:
            m.append(f'% {c}')
        m.append(f'\\newcommand{{\\{k}}}{{{v}}}')

    # ------------------------------------------------------------------ monophyly
    g = rd(os.path.join(a.checks, 'trees', 'monophyly.groups.tsv'))

    def cell(tree, cset, level):
        xs = [r for r in g if r['tree'] == tree and r['copy_set'] == cset and r['level'] == level]
        return sum(r['status'] == 'monophyletic' for r in xs), len(xs)

    trees = [('ssu_flank', 'SSU flank'), ('full_its', 'Full ITS'), ('lsu_flank', 'LSU flank'),
             ('linked', 'Linked copy')]
    rows = []
    for t, lab in trees:
        c = [cell(t, cs, lv) for lv in ('culture', 'species', 'complex') for cs in ('all', 'main')]
        rows.append(f'{lab} & ' + ' & '.join(f'{k}/{n}' for k, n in c) + '\\\\')
    tag = {'ssu_flank': 'SSU', 'full_its': 'ITS', 'lsu_flank': 'LSU', 'linked': 'Linked'}
    for t, _ in trees:
        for cs, cn in (('all', 'All'), ('main', 'Main')):
            add(f'MonoCult{tag[t]}{cn}', cell(t, cs, 'culture')[0])
            add(f'MonoSp{tag[t]}{cn}', cell(t, cs, 'species')[0])
    add('MonoCultN', cell('ssu_flank', 'all', 'culture')[1], 'cultures with >= 2 copies')
    pct = [100 * cell(t, cs, 'culture')[0] / cell(t, cs, 'culture')[1] for t, _ in trees for cs in ('all', 'main')]
    add('MonoCultPctMin', f'{min(pct):.0f}\\%')
    add('MonoCultPctMax', f'{max(pct):.0f}\\%')
    esc = [int(r['escapees_n']) / int(r['n_copies']) for r in g if r['tree'] == 'linked' and r['copy_set'] == 'main'
           and r['level'] == 'culture' and r['status'] == 'non-monophyletic']
    add('MonoCultEscMedian', f'{100 * statistics.median(esc):.0f}\\%',
        "median share of a non-monophyletic culture's copies outside its largest pure clade (linked, main type)")
    lm = {r['group']: r for r in g if r['tree'] == 'linked' and r['copy_set'] == 'main' and r['level'] == 'species'}
    la = {r['group']: r for r in g if r['tree'] == 'linked' and r['copy_set'] == 'all' and r['level'] == 'species'}
    fixed = sorted(k for k in la if la[k]['status'] != 'monophyletic' and lm[k]['status'] == 'monophyletic')
    add('MonoFixedN', len(fixed), 'species monophyletic on the linked tree only once divergent copies are pruned')
    add('MonoComplexLinkedMain', '%d/%d' % cell('linked', 'main', 'complex'))

    with open(os.path.join(a.paper, 'tables', 'table_monophyly.tex'), 'w') as fh:
        fh.write('''\\begin{table}[htbp]\\centering
\\begin{threeparttable}
\\caption{Cultures, named species and species complexes that form a single clade (are monophyletic) on trees of
the culture copies.}
\\label{tab:monophyly}
\\footnotesize\\setlength{\\tabcolsep}{4pt}
\\begin{tabular}{lrrrrrr}
\\toprule
 & \\multicolumn{2}{c}{Cultures} & \\multicolumn{2}{c}{Named species} & \\multicolumn{2}{c}{Complexes}\\\\
\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}
Tree & all & main & all & main & all & main\\\\
\\midrule
''' + '\n'.join(rows) + '''
\\bottomrule
\\end{tabular}
\\begin{tablenotes}\\footnotesize
\\item Monophyletic groups / groups tested; a group is tested when it has two or more copies. \\emph{all}: all copies;
\\emph{main}: the same trees with the divergent copies removed. One tree per order and region (maximum likelihood,
IQ-TREE); the linked copy joins the three regions. Alignments and trees are in Supplementary Data S1.
\\end{tablenotes}
\\end{threeparttable}
\\end{table}
''')

    # ------------------------------------------------------------------ artifact check
    c = rd(os.path.join(a.checks, 'artefacts', 'artefacts.copies.tsv'))
    dv = [r for r in c if r['copy_class'] == 'divergent']
    V = ['paralog', 'recurrent deleted copy', 'unresolved', 'probable artefact']
    arows = []
    for cls, lab in (('short_5.8S', 'Short 5.8S'), ('same_length', 'Same length')):
        xs = [r for r in dv if r['div_kind'] == cls]
        cnt = collections.Counter(r['verdict'] for r in xs)
        arows.append(f'{lab} ({len(xs)}) & ' + ' & '.join(str(cnt[v]) for v in V) + '\\\\')
        key = 'Short' if cls == 'short_5.8S' else 'Same'
        for v, vk in zip(V, ('Paralog', 'RecDel', 'Unres', 'Artefact')):
            add(f'Art{key}{vk}', cnt[v])
        add(f'Art{key}ArtefactTen', sum(r['verdict_del10nt'] == 'probable artefact' for r in xs))
    thr = {r['region']: float(r['L_R']) for r in rd(os.path.join(a.checks, 'artefacts', 'artefacts.thresholds.tsv'))}
    add('ArtLConserved', f"{thr['s58']:g}", 'long-deletion threshold, nt')
    add('ArtLSpacer', f"{thr['its1']:g}")
    add('ArtMainFlagged', sum(r['copy_class'] == 'main' and any(x in r['long_deletion_regions'].split(',')
                                                                   for x in ('ssu_flank', 's58', 'lsu_flank')) for r in c))

    # post hoc: does the 5.8S deletion itself recur under other labels?
    ssu = lens(os.path.join(a.ref, 'fasta', 'ssu_flank.fasta'))
    its1 = lens(os.path.join(a.ref, 'fasta', 'its1.fasta'))
    dels = collections.defaultdict(list)
    for r in c:
        for i in r['indels'].split(';'):
            if not i:
                continue
            k, reg, p, L = i.split(':')
            if k == 'del' and reg == 's58':
                cm = r['comparator']
                dels[r['accession']].append((int(L), int(p) - 1 - ssu[cm] - its1[cm]))
    by = {r['accession']: r for r in c}
    shared = {}
    for acc, ds in dels.items():
        labs = set()
        for L, o in ds:
            if L < thr['s58']:
                continue
            for acc2, ds2 in dels.items():
                if by[acc2]['culture_label'] != by[acc]['culture_label'] and any(L2 == L and abs(o2 - o) <= 2 for L2, o2 in ds2):
                    labs.add(by[acc2]['culture_label'])
        shared[acc] = labs
    sh = [r for r in dv if r['div_kind'] == 'short_5.8S']
    add('ShDelRecur', sum(bool(shared.get(r['accession'])) for r in sh),
        'post hoc: short-5.8S copies whose 5.8S deletion (same length, start within 2 nt) occurs under another label')
    add('ShDelRecurArtefact', sum(bool(shared.get(r['accession'])) for r in sh if r['verdict'] == 'probable artefact'))
    six = [r for r in c if any(L == 6 and 119 <= o <= 124 for L, o in dels.get(r['accession'], []))]
    add('ShSixCopies', sum(r['copy_class'] == 'divergent' for r in six))
    add('ShSixMain', sum(r['copy_class'] == 'main' for r in six))
    add('ShSixCultures', len({r['culture'] for r in six}))

    with open(os.path.join(a.paper, 'tables', 'table_artefacts.tex'), 'w') as fh:
        fh.write('''\\begin{table}[htbp]\\centering
\\begin{threeparttable}
\\caption{Divergent copies by whether the same copy type occurs in other cultures (recurring) and whether it has
long deletions.}
\\label{tab:artefacts}
\\footnotesize\\setlength{\\tabcolsep}{4pt}
\\begin{tabular}{lrrrr}
\\toprule
 & \\multicolumn{2}{c}{Recurring} & \\multicolumn{2}{c}{Single}\\\\
\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}
Class (copies) & intact & deleted & intact & deleted\\\\
\\midrule
''' + '\n'.join(arows) + f'''
\\bottomrule
\\end{{tabular}}
\\begin{{tablenotes}}\\footnotesize
\\item \\emph{{Recurring}}: a copy within 1.0\\% in full ITS occurs under another culture label. \\emph{{Deleted}}: compared
with the most similar main-type copy of the same culture, the copy lacks at least {thr['s58']:g}~nt in a row in the SSU,
5.8S or LSU, or {thr['its1']:g}~nt in ITS1 or ITS2 (longer than 99\\% of such gaps between main-type copies). The four
cells correspond to the verdicts fixed before the check: paralog, recurrent deleted copy, unresolved and probable
artifact.
\\end{{tablenotes}}
\\end{{threeparttable}}
\\end{{table}}
''')

    # ------------------------------------------------------------------ function check (Addendum 1)
    fs = {(r['region'], r['group']): r for r in rd(os.path.join(a.checks, 'function', 'function.summary.tsv'))}
    fj = json.load(open(os.path.join(a.checks, 'function', 'function.json')))
    fc = rd(os.path.join(a.checks, 'function', 'function.copies.tsv'))
    tag = {'ssu_flank': 'SSU', 's58': 'Five', 'lsu_flank': 'LSU'}
    for reg, t in tag.items():
        s = fs[(reg, 'short 5.8S')]
        b = fs[(reg, 'baseline')]
        add(f'Fn{t}Irr', s['irregular'], f'function check, {reg}, short-5.8S copies')
        add(f'Fn{t}Subs', f"{float(s['median_subs_per100']):.2f}")
        add(f'Fn{t}SubsBase', f"{float(b['median_subs_per100']):.2f}")
        add(f'Fn{t}Deg', s['degraded'])
        add(f'Fn{t}BrokenShort', f"{100 * float(s['median_broken_frac']):.1f}\\%")
        add(f'Fn{t}BrokenBase', f"{100 * float(b['median_broken_frac']):.1f}\\%")
        add(f'Fn{t}Diff', f"{100 * float(s['paired_median_diff']):+.1f}")
        p = float(s['paired_wilcoxon_p'])
        if p < 0.001:
            e = int(f'{p:.0e}'.split('e')[1])
            add(f'Fn{t}P', f'$<10^{{{e + 1}}}$')
        else:
            add(f'Fn{t}P', f'{p:.2f}')
    add('FnLSUDegPct', f"{100 * int(fs[('lsu_flank', 'short 5.8S')]['degraded']) / int(fs[('lsu_flank', 'short 5.8S')]['copies']):.0f}\\%")
    add('FnShortN', fs[('lsu_flank', 'short 5.8S')]['copies'])
    for reg, t in tag.items():
        add(f'FnSame{t}Deg', fs[(reg, 'same length')]['degraded'])
    add('FnVerdictShort', fj['verdicts']['short 5.8S'])
    add('FnVerdictSame', fj['verdicts']['same length'])
    d = rd(os.path.join(a.checks, 'function', 'function.deletions.tsv'))
    hp = [r for r in d if r['group'] == 'short 5.8S' and int(r['n_paired']) >= 4
          and set(r['deleted_consensus_positions'].split(',')) & {'121', '122', '123', '128', '129', '130'}]
    add('FnHairpinShort', len(hp), 'short copies whose 5.8S deletion removes paired positions of the 121-130 hairpin')
    add('FnHairpinMain', sum(r['group'] == 'flagged main' and int(r['n_paired']) >= 4 for r in d))
    sv = [r for r in d if r['group'] == 'short 5.8S' and int(r['n_paired']) == 0]
    add('FnUnpairedShort', len(sv), 'short copies whose 5.8S deletion removes only unpaired positions')

    # the shortened hairpin (post hoc, descriptive): stems and loop read from the 5.8S between the conserved
    # anchors CGCAAAUUGC (consensus 105-114) and AUGC (137-140)
    s58 = {}
    name = None
    for line in open(os.path.join(a.ref, 'fasta', 's58.fasta')):
        if line.startswith('>'):
            name = line[1:].split()[0]
            s58[name] = ''
        else:
            s58[name] += line.strip().upper().replace('T', 'U')
    wc = {'AU', 'UA', 'GC', 'CG'}
    loops = collections.Counter()
    allwc = 0
    for r in hp:
        x = s58[r['accession']]
        i = x.find('CGCAAAUUGC')
        j = x.find('AUGC', i + 22)
        st5, st3, loop = x[i + 10:i + 16], x[j - 6:j], x[i + 16:j - 6]
        loops[loop] += 1
        allwc += all(st5[k] + st3[5 - k] in wc for k in range(6))
    add('FnLoopUUCG', loops['UUCG'], 'loops of the shortened 5.8S hairpin in the short copies (post hoc)')
    add('FnLoopUUUG', loops['UUUG'])
    add('FnLoopCUCG', loops['CUCG'])
    add('FnLoopOther', sum(v for k, v in loops.items() if k not in ('UUCG', 'UUUG', 'CUCG')))
    add('FnStemAllWC', allwc, 'short copies whose shortened 6-bp stem is all Watson-Crick')

    def med(group, col):
        v = [float(r[col]) for r in fc if r['group'] == group and r[col] != '']
        return statistics.median(v)
    add('FnITSTwoBPShort', f"{med('short 5.8S', 'its2_bp_distance_per_nt'):.2f}")
    add('FnITSTwoBPBase', f"{med('baseline', 'its2_bp_distance_per_nt'):.2f}")
    add('FnITSTwoMFEShort', f"{med('short 5.8S', 'its2_mfe_per_nt'):.2f}")
    add('FnITSTwoMFEBase', f"{med('baseline', 'its2_mfe_per_nt'):.2f}")

    # ------------------------------------------------------------------ thresholds (exploratory)
    q = {(r['region'], r['class']): r for r in rd(os.path.join(a.checks, 'thresholds', 'thresholds.quantiles.tsv'))}
    lim = {r['region']: float(r['within_species_spread_p99'])
           for r in rd(os.path.join(a.checks, 'thresholds', 'thresholds.limits.tsv'))}
    for reg, t in (('ssu_flank', 'SSU'), ('full_its', 'ITS'), ('LSU_LROR-FLR2', 'LSU'), ('SSU_VT_NS31-AM1', 'VT')):
        add(f'Thr{t}Congener', f"{float(q[(reg, 'same genus')]['below_spread']):.1f}\\%",
            f'{reg}: share of main-type copies with a congener within the within-species spread')
        add(f'Thr{t}SameCulture', f"{float(q[(reg, 'same culture')]['median']):.2f}\\%")
        add(f'Thr{t}SameSpecies', f"{float(q[(reg, 'same species')]['median']):.2f}\\%")
        add(f'Thr{t}OtherGenus', f"{float(q[(reg, 'other genus')]['median']):.1f}\\%")
        add(f'Thr{t}OtherGenusSpread', f"{float(q[(reg, 'other genus')]['below_spread']):.1f}\\%")
    add('ThrVTSpread', f"{lim['SSU_VT_NS31-AM1']:.1f}\\%", 'within-species spread on the VT fragment (99th pct)')
    add('ThrVTCongenerDB', f"{float(q[('SSU_VT_NS31-AM1', 'same genus')]['below_db']):.1f}\\%",
        'VT fragment: share of main-type copies with a congener / another genus within 3%')
    add('ThrVTOtherGenusDB', f"{float(q[('SSU_VT_NS31-AM1', 'other genus')]['below_db']):.1f}\\%")
    nr = rd(os.path.join(a.checks, 'thresholds', 'thresholds.nearest.tsv'))
    gen = {r['accession']: r['genus'] for r in rd(os.path.join(a.ref, 'copies.tsv'))}
    vt3 = [r for r in nr if r['region'] == 'SSU_VT_NS31-AM1' and r['class'] == 'other genus'
           and float(r['nearest_pct']) <= 3.0]
    add('ThrVTOtherGenusN', len(vt3))
    add('ThrVTOtherGenusDen', sum(r['region'] == 'SSU_VT_NS31-AM1' and r['class'] == 'other genus' for r in nr))
    add('ThrVTOtherGenusNotGiga', sum(gen[r['accession']] not in ('Cetraspora', 'Gigaspora', 'Dentiscutata') for r in vt3),
        'of them, copies outside Cetraspora, Gigaspora and Dentiscutata')
    add('ThrITSSameCultureAboveSH', f"{100 - float(q[('full_its', 'same culture')]['below_db']):.1f}\\%",
        'full ITS: share of main-type copies more than 1.5% from every other main-type copy of their culture')
    add('ThrITSSameSpeciesAboveSH', f"{100 - float(q[('full_its', 'same species')]['below_db']):.1f}\\%")
    add('ThrITSCongenerDB', f"{float(q[('full_its', 'same genus')]['below_db']):.1f}\\%")
    pr = rd(os.path.join(a.checks, 'thresholds', 'thresholds.pairs.tsv'))
    above = [r for r in pr if float(r['ssu_flank']) > float(r['full_its']) / 2]
    add('ThrPairsN', f'{len(pr):,}'.replace(',', '{,}'))
    add('ThrPairsAbove', len(above), 'culture pairs with SSU distance above half the ITS distance')
    add('ThrPairsAboveMaxSSU', f"{max(float(r['ssu_flank']) for r in above):.2f}\\%" if above else '--')
    add('ThrPairsRuleA', sum(float(r['ssu_flank']) >= 1 and float(r['ssu_flank']) >= float(r['full_its']) / 2 for r in pr))
    rat = [float(r['ssu_flank']) / float(r['full_its']) for r in pr if float(r['full_its']) >= 2.0]
    add('ThrPairsRatio', f'{statistics.median(rat):.2f}', 'median SSU-flank / full-ITS distance, culture pairs >= 2% apart in ITS')

    # ------------------------------------------------------------------ Rhizophagus clade of the linked tree
    from Bio import Phylo
    ref = {r['accession']: r for r in rd(os.path.join(a.ref, 'copies.tsv')) if r['status'] == 'retained'}
    tdir = os.path.join(a.checks, 'trees', 'Glomerales')
    mem = {r['tip']: r['accessions'].split(',') for r in rd(os.path.join(tdir, 'linked.members.tsv'))}
    tree = Phylo.read(os.path.join(tdir, 'linked.treefile'), 'newick')
    info = lambda n: [ref[x] for x in mem[n]]  # noqa: E731
    out = next(tp for tp in tree.get_terminals() if all(r['genus'] == 'Funneliformis' for r in info(tp.name)))
    tree.root_with_outgroup(out)
    rh = [tp for tp in tree.get_terminals() if any(r['genus'] in ('Rhizophagus', 'Rhizoglomus') for r in info(tp.name))]
    clade = tree.common_ancestor(rh)
    isshort = lambda n: any(r['div_kind'] == 'short_5.8S' for r in info(n))  # noqa: E731
    pure = []

    def walk(cl):
        if all(isshort(tp.name) for tp in cl.get_terminals()):
            pure.append(cl)
            return
        for ch in cl.clades:
            walk(ch)
    walk(clade)
    parts = []
    for cl in sorted(pure, key=lambda c: -sum(len(mem[t.name]) for t in c.get_terminals())):
        rs = [r for t in cl.get_terminals() for r in info(t.name)]
        nshort = sum(r['div_kind'] == 'short_5.8S' for r in rs)
        nmain = sum(r['copy_class'] == 'main' for r in rs)
        orgs = sorted({r['organism'] for r in rs if r['div_kind'] == 'short_5.8S'})
        sup = cl.confidence if cl.confidence is not None else cl.name
        sup = sup if isinstance(sup, str) else ''
        mm = re.match(r'([\d.]+)/([\d.]+)', sup or '')
        who = ' with '.join(emph(o.replace(' (in: glomeromycetes)', '')) for o in orgs)
        extra = f' and {nmain} copy filed as main type' if nmain == 1 else (f' and {nmain} copies filed as main type' if nmain else '')
        if len(cl.get_terminals()) > 1 and mm:
            parts.append(f'{who} ({nshort} copies{extra}; {float(mm.group(1)):g}/{float(mm.group(2)):g})')
        else:
            parts.append(f'{who} ({nshort} cop{"y" if nshort == 1 else "ies"}{extra})')
    add('TreeShortClades', '; '.join(parts), 'clades made only of short-5.8S copies in the Rhizophagus clade (SH-aLRT/UFBoot)')
    add('TreeShortCladesN', len(pure))
    multi = [cl for cl in pure if len(cl.get_terminals()) > 1]
    nshort = lambda cl: sum(r['div_kind'] == 'short_5.8S' for t in cl.get_terminals() for r in info(t.name))  # noqa: E731
    ntaxa = lambda cl: len({r['organism'] for t in cl.get_terminals() for r in info(t.name)  # noqa: E731
                            if r['div_kind'] == 'short_5.8S'})
    add('TreeShortTotal', sum(r['div_kind'] == 'short_5.8S' for t in clade.get_terminals() for r in info(t.name)),
        'short copies in the Rhizophagus clade; of them, in short-only clades of two or more tips; in such clades of two taxa')
    add('TreeShortInMulti', sum(nshort(cl) for cl in multi))
    add('TreeShortMultiN', len(multi))
    add('TreeShortCross', sum(nshort(cl) for cl in multi if ntaxa(cl) > 1))
    add('TreeShortCrossN', sum(ntaxa(cl) > 1 for cl in multi))
    add('TreeRhizTips', len(clade.get_terminals()))
    add('TreeRhizCopies', sum(len(mem[t.name]) for t in clade.get_terminals()))

    # ------------------------------------------------------------------ long fragments (exploratory, 13 Oct 2026)
    ld = os.path.join(a.checks, 'thresholds_long')
    L = {r['region']: r for r in rd(os.path.join(ld, 'long.limits.tsv'))}
    S = {(r['region'], r['threshold']): r for r in rd(os.path.join(ld, 'long.summary.tsv'))}

    def nt(x):
        return f'{float(x):,.0f}'.replace(',', '{,}')
    add('LongKrMedian', nt(L['kruger']['median_length']), 'Kruger fragment (SSUmAf1 site to LSUmBr) and whole copy: lengths, nt')
    add('LongKrMin', nt(L['kruger']['min_length']))
    add('LongKrMax', nt(L['kruger']['max_length']))
    add('LongFullMedian', nt(L['full']['median_length']))
    add('LongKrSpread', f"{float(L['kruger']['within_species_spread_p99']):.1f}\\%", 'within-species spread (99th pct)')
    add('LongFullSpread', f"{float(L['full']['within_species_spread_p99']):.1f}\\%")
    for reg, tag in (('full_its', 'ITS'), ('kruger', 'Kr'), ('full', 'Full')):
        for thr, t2 in (('97% identity', 'NS'), ('98.5% identity', 'NE'), ('best single threshold', 'Best'),
                        ('best single threshold, complexes aside', 'BestX'), ('within-species spread', 'Spr')):
            r = S[(reg, thr)]
            add(f'Long{tag}{t2}Id', f"{float(r['identity_pct']):.2f}\\%",
                f'{reg} at {thr}: identity, split from culture, split from species, joined (all; outside complexes)')
            add(f'Long{tag}{t2}Split', f"{float(r['split_from_genome_pct']):.1f}\\%")
            add(f'Long{tag}{t2}SplitSp', f"{float(r['split_from_species_pct']):.1f}\\%")
            add(f'Long{tag}{t2}Join', f"{float(r['joined_to_congener_pct']):.1f}\\%")
            add(f'Long{tag}{t2}JoinX', f"{float(r['joined_to_congener_excl_complex_pct']):.1f}\\%")

    # ------------------------------------------------------------------ what the short paralog does (exploratory)
    pd = os.path.join(a.checks, 'paralog_effects')
    C = {(r['region'], r['identity'], r['copies']): r for r in rd(os.path.join(pd, 'paralog.counts.tsv'))}
    for reg, tag in (('full_its', 'ITS'), ('kruger', 'Kr'), ('full', 'Full')):
        for idn, t2 in (('98.5%', 'NE'), ('97%', 'NS')):
            ca, cm = C[(reg, idn, 'all copies')], C[(reg, idn, 'main-type copies')]
            add(f'Par{tag}{t2}All', ca['clusters'], f'{reg} at {idn}: clusters of all copies, of main-type copies, '
                'short-only and same-length-only clusters')
            add(f'Par{tag}{t2}Main', cm['clusters'])
            add(f'Par{tag}{t2}Short', ca['clusters_short_only'])
            add(f'Par{tag}{t2}Same', ca['clusters_same_length_only'])
            add(f'Par{tag}{t2}Join', ca['clusters_joining_species'])
    ps = [r for r in rd(os.path.join(pd, 'paralog.species.tsv')) if r['region'] == 'full_its' and r['identity'] == '98.5%']
    so = [int(r['of_which_short_only_clusters']) for r in ps]
    add('ParITSNESpMin', min(so), 'full ITS at 98.5%: short-only clusters per species that carries the paralog')
    add('ParITSNESpMax', max(so))
    add('ParITSNESpN', len(ps))
    add('ParITSNEShared', sum(int(r['short_only_clusters_shared_with_other_species']) for r in ps))
    T = {r['tree']: r for r in rd(os.path.join(pd, 'paralog.trees.tsv'))}
    add('ParShortN', T['ssu_flank']['short_copies'], 'short copies by 5.8S sequence (Glomerales)')
    add('ParSSUIdentical', T['ssu_flank']['short_in_tips_with_main_copies'], 'SSU flank: identical to a main-type copy')
    add('ParSSUPure', T['ssu_flank']['copies_in_short_only_clades_2plus'])
    add('ParSSUCrossN', T['ssu_flank']['short_only_clades_joining_taxa'])
    for tr, tag in (('full_its', 'ITS'), ('lsu_flank', 'LSU'), ('linked', 'Linked')):
        add(f'Par{tag}Pure', T[tr]['copies_in_short_only_clades_2plus'], f'{tr}: copies in short-only clades of 2+, '
            'clades joining taxa, their copies')
        add(f'Par{tag}CrossN', T[tr]['short_only_clades_joining_taxa'])
        add(f'Par{tag}Cross', T[tr]['copies_in_clades_joining_taxa'])
    tc = rd(os.path.join(pd, 'paralog.tree_clades.tsv'))
    ip = next(r for r in tc if r['tree'] == 'lsu_flank' and r['taxa'] == 'Rhizophagus intraradices; Rhizophagus prolifer')
    add('ParLSUIntraProlN', ip['copies'], 'LSU flank: the R. intraradices/R. prolifer short clade, its support, its sister')
    add('ParLSUIntraProlSup', f"{float(ip['sh_alrt']):g}/{float(ip['ufboot']):g}")
    add('ParLSUIntraProlSister', emph(ip['sister_main_taxa']))
    add('ParLSUIntraProlSisterSup', f"{float(ip['with_sister_sh_alrt']):g}/{float(ip['with_sister_ufboot']):g}")

    # ------------------------------------------------------------------ EUKARYOME long reads (exploratory)
    lr = os.path.join(a.checks, '..', 'eukaryome_longread')
    LS = {(r['region'], r['threshold']): r for r in rd(os.path.join(lr, 'longread.summary.tsv'))}
    sp_rows = rd(os.path.join(lr, 'longread.species.tsv'))
    add('LrReads', sum(int(r['reads']) for r in sp_rows), 'EUKARYOME long reads of named Glomeromycota species')
    add('LrSpecies', len(sp_rows))
    for reg, tag in (('its', 'ITS'), ('kruger', 'Kr'), ('read', 'Read')):
        for thr, t2 in (('97% identity', 'NS'), ('98.5% identity', 'NE'), ('best single cut-off', 'Best')):
            r = LS[(reg, thr)]
            add(f'Lr{tag}{t2}Id', f"{float(r['identity_pct']):.2f}\\%", f'{reg} at {thr}: identity, split, joined')
            add(f'Lr{tag}{t2}Split', f"{float(r['split_from_species_pct']):.1f}\\%")
            add(f'Lr{tag}{t2}Join', f"{float(r['joined_to_congener_pct']):.1f}\\%")
        add(f'Lr{tag}Spread', f"{100 - float(LS[(reg, 'within-species spread')]['identity_pct']):.1f}\\%")
    lp = [r for r in rd(os.path.join(lr, 'longread.pairs.tsv')) if r['region'] == 'read' and r['cutoff'] == 'best single cut-off']
    add('LrReadBestPairs', len(lp), 'whole read, best cut-off: species pairs behind the joins, and their reads')
    add('LrReadBestPairReads', sum(int(r['reads_joined']) for r in lp))
    lp97 = [r for r in rd(os.path.join(lr, 'longread.pairs.tsv')) if r['region'] == 'read' and r['cutoff'] == '97% identity']
    add('LrReadNSPairs', len(lp97))
    # whole read, best cut-off: reads split from their species or joined to a congener, either way (of all reads)
    near = collections.defaultdict(dict)
    for r in rd(os.path.join(lr, 'longread.nearest.tsv')):
        if r['region'] == 'read':
            near[r['read']][r['class']] = float(r['nearest_pct'])
    t = round(100 - float(LS[('read', 'best single cut-off')]['identity_pct']), 2)
    either = sum(1 for d in near.values()
                 if d.get('same species', -1) > t or ('same genus' in d and d['same genus'] <= t))
    add('LrReadBestEither', f"{100 * either / len(near):.1f}\\%",
        'whole read, best cut-off: reads split or joined, of all reads; and reads whose genus has another species')
    add('LrCongenerN', LS[('read', 'best single cut-off')]['reads_with_congener'])

    # Table S7: clusters per region and identity
    with open(os.path.join(a.paper, 'tables', 'table_clusters.tex'), 'w') as fh:
        fh.write('% written by analysis/copy_checks_numbers.py from benchmarks/copy_checks/paralog_effects\n')
        fh.write('\\begin{table}[htbp]\\centering\n\\small\\setlength{\\tabcolsep}{4pt}\n\\caption{Single-linkage clusters of the copies of the '
                 '\\NSpecies{} named species, within each order, at two identity thresholds: all retained copies, and '
                 'main-type copies only. Short only and same length only: clusters made only of short-5.8S or only of '
                 'same-length divergent copies; several species: clusters holding copies of more than one species. '
                 'Exploratory; no rule depends on it.}\\label{tab:clusters}\n')
        fh.write('\\begin{tabular}{llrrrrr}\n\\toprule\n & & \\multicolumn{2}{c}{Clusters} & '
                 '\\multicolumn{3}{c}{Clusters of all copies}\\\\\n\\cmidrule(lr){3-4}\\cmidrule(lr){5-7}\n'
                 'Region & Identity & all copies & main type & short only & same length only & several species'
                 '\\\\\n\\midrule\n')
        for reg, lab in (('full_its', 'Full ITS'), ('kruger', 'Kr\\"uger fragment'), ('full', 'Whole copy')):
            for idn in ('98.5%', '97%'):
                ca, cm = C[(reg, idn, 'all copies')], C[(reg, idn, 'main-type copies')]
                pc = idn.replace('%', '\\%')
                fh.write(f"{lab} & {pc} & {ca['clusters']} & {cm['clusters']} & "
                         f"{ca['clusters_short_only']} & {ca['clusters_same_length_only']} & "
                         f"{ca['clusters_joining_species']}\\\\\n")
        fh.write('\\bottomrule\n\\end{tabular}\n\\end{table}\n')

    with open(os.path.join(a.paper, 'numbers_copychecks.tex'), 'w') as fh:
        fh.write('% written by analysis/copy_checks_numbers.py from benchmarks/copy_checks\n' + '\n'.join(m) + '\n')
    print('\n'.join(x for x in m if not x.startswith('%')))
    print('species fixed by pruning (linked):', fixed)


if __name__ == '__main__':
    main()
