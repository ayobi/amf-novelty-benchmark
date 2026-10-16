#!/usr/bin/env python3
"""The evidence for each class of divergent copy, side by side (added after review, point 1.4).

Reads the outputs of the copy checks (benchmarks/copy_checks/artefacts and function) and the v4 reference; computes
nothing new except the grouping. Classes: the short-5.8S copies of Rhizophagus and [Rhizoglomus]; the short-5.8S
copies of Septoglomus viscosum; the same-length copies of Gigaspora rosea; the other same-length copies.

Evidence lines, each from the copy checks:
  * shared deletion: a 5.8S deletion of the same length, starting within 2 nt, in a copy under another culture label
    (copy_checks_numbers.py's post hoc reading), and whether one of those labels is another species;
  * whole ITS recurs: a copy within 1.0% in full ITS under another label (the check's 'recurring');
  * library replicate: a copy of the same type in another library of the same culture (12 cultures have 2-4);
  * public record: a UNITE+INSD record at >= 99% identity over >= 90% coverage, culture copies excluded;
  * degraded: share of broken base pairs above the 99th percentile of main-type copies in the SSU, 5.8S or LSU.

  python3 analysis/divergent_classes.py --checks ../benchmarks/copy_checks --ref ../reference/culture_ref_v4 --paper .
Writes tables/table_classes.tex and numbers_classes.tex.
"""
import argparse
import collections
import csv
import os
import re

UNNAMED = re.compile(r"\bsp\.?(\s|$)")   # as stage3/arms.py


def rd(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def lens(path):
    out, k = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            k = line[1:].split()[0]
            out[k] = 0
        elif k:
            out[k] += len(line)
    return out


CLASSES = [('rhizo', 'Short 5.8S, \\spn{Rhizophagus}'),
           ('septo', 'Short 5.8S, \\spn{Septoglomus viscosum}'),
           ('giga', 'Same length, \\spn{Gigaspora rosea}'),
           ('other', 'Same length, others')]


def cls(r):
    if r['div_kind'] == 'short_5.8S':
        return 'rhizo' if r['genus'] in ('Rhizophagus', 'Rhizoglomus') else 'septo'
    return 'giga' if r['genus'] == 'Gigaspora' else 'other'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--checks', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--paper', default='.')
    ap.add_argument('--comparison', default='../benchmarks/m4_external_holdout/frozen_m2b/comparison',
                    help='frozen M2b tables of the withholding tests (LSU placement, databases intact)')
    a = ap.parse_args()
    rec = {}
    for mode in ('culture', 'species'):
        for r in rd(os.path.join(a.comparison, f'epa_species_clades_physical_{mode}.copies.tsv')):
            rec[(mode, r['accession'])] = r['status']
    c = rd(os.path.join(a.checks, 'artefacts', 'artefacts.copies.tsv'))
    fn = {r['accession']: r for r in rd(os.path.join(a.checks, 'function', 'function.copies.tsv'))}
    thr = {r['region']: float(r['L_R']) for r in rd(os.path.join(a.checks, 'artefacts', 'artefacts.thresholds.tsv'))}
    # shared 5.8S deletions, exactly as copy_checks_numbers.py reads them
    ssu = lens(os.path.join(a.ref, 'fasta', 'ssu_flank.fasta'))
    its1 = lens(os.path.join(a.ref, 'fasta', 'its1.fasta'))
    dels = collections.defaultdict(list)
    for r in c:
        for i in r['indels'].split(';'):
            if i:
                k, reg, p, L = i.split(':')
                if k == 'del' and reg == 's58':
                    cm = r['comparator']
                    dels[r['accession']].append((int(L), int(p) - 1 - ssu[cm] - its1[cm]))
    by = {r['accession']: r for r in c}
    shared = {}
    for acc, ds in dels.items():
        hits = set()
        for L, o in ds:
            if L < thr['s58']:
                continue
            for acc2, ds2 in dels.items():
                if by[acc2]['culture_label'] != by[acc]['culture_label'] and any(L2 == L and abs(o2 - o) <= 2 for L2, o2 in ds2):
                    hits.add(acc2)
        shared[acc] = hits
    dv = [r for r in c if r['copy_class'] == 'divergent']
    g = collections.defaultdict(list)
    for r in dv:
        g[cls(r)].append(r)

    def degraded(r):
        f = fn.get(r['accession'], {})
        return any(f.get(f'{reg}_degraded') == 'yes' for reg in ('ssu_flank', 's58', 'lsu_flank'))

    stats = {}
    for k, _ in CLASSES:
        rs = g[k]
        named = {r['organism'] for r in rs if not UNNAMED.search(r['organism'])}
        sh = [r for r in rs if shared.get(r['accession'])]
        cross = [r for r in sh if any(by[x]['organism'] != r['organism'] for x in shared[r['accession']])]
        stats[k] = {
            'copies': len(rs), 'cultures': len({r['culture'] for r in rs}), 'species': len(named),
            'sp_only': len({r['organism'] for r in rs}) - len(named),
            'shared': len(sh), 'shared_cross': len(cross),
            'recur': sum(int(r['recur_labels_n'] or 0) > 0 for r in rs),
            'librep_yes': sum(r['library_replicate'] == 'yes' for r in rs),
            'librep_tested': sum(r['library_replicate'] in ('yes', 'no') for r in rs),
            'external': sum(r['external_record'] == 'yes' for r in rs),
            'degraded': sum(degraded(r) for r in rs),
            'deg_regions': ' / '.join(str(sum(fn.get(r['accession'], {}).get(f'{reg}_degraded') == 'yes' for r in rs))
                                      for reg in ('ssu_flank', 's58', 'lsu_flank')),
            'verdicts': collections.Counter(r['verdict'] for r in rs),
            'rec_culture': sum(rec[('culture', r['accession'])] == 'divergent class' for r in rs),
            'rec_species': sum(rec[('species', r['accession'])] == 'divergent class' for r in rs),
        }
    S = stats

    def cell(k, n, d=None):
        return f'{S[k][n]}' if d is None else f'{S[k][n]} of {S[k][d]}'
    rows = [
        ('Copies (cultures)', [f"{S[k]['copies']} ({S[k]['cultures']})" for k, _ in CLASSES]),
        ('Named species', [f"{S[k]['species']}" + (f" + sp." if S[k]['sp_only'] else '') for k, _ in CLASSES]),
        ('5.8S deletion shared with another label', [cell(k, 'shared') if k in ('rhizo', 'septo') else '--' for k, _ in CLASSES]),
        ('\\quad with another species', [cell(k, 'shared_cross') if k in ('rhizo', 'septo') else '--' for k, _ in CLASSES]),
        ('Recognized: culture / species withheld', [f"{S[k]['rec_culture']} / {S[k]['rec_species']}" if k in ('rhizo', 'septo') else '--' for k, _ in CLASSES]),
        ('Whole ITS recurs under another label', [cell(k, 'recur') for k, _ in CLASSES]),
        ('Same type in another library', [f"{S[k]['librep_yes']} of {S[k]['librep_tested']}" if S[k]['librep_tested'] else '--' for k, _ in CLASSES]),
        ('Public record ($\\ge$99\\%)', [cell(k, 'external') for k, _ in CLASSES]),
        ('Degraded structure (SSU / 5.8S / LSU)', [S[k]['deg_regions'] for k, _ in CLASSES]),
    ]
    body = '\n'.join(f"{lab} & " + ' & '.join(v) + '\\\\' for lab, v in rows)
    tex = r'''\begin{table}[htbp]\centering
\begin{threeparttable}
\caption{The evidence differs between classes of divergent copy.}
\label{tab:classes}
\footnotesize\setlength{\tabcolsep}{3pt}
\begin{tabular}{lcccc}
\toprule
 & Short 5.8S & Short 5.8S & Same length & Same length\\
 & \spn{Rhizophagus} & \spn{S.~viscosum} & \spn{G.~rosea} & others\\
\midrule
''' + body + r'''
\midrule
Reading & paralog family & genomic & co-dominant types & unresolved\\
\bottomrule
\end{tabular}
\begin{tablenotes}\footnotesize
\item \emph{Shared deletion}: a 5.8S deletion of the same length, starting within 2~nt, in a copy under another culture
label (read after the check). \emph{Whole ITS recurs}: a copy within 1.0\% in full ITS under another label, the
criterion of the artifact check. \emph{Recognized}: called a divergent copy by the decision steps in the tests
with withheld cultures and species (Table~\ref{tab:benchmark}); same-length copies are only flagged. \emph{Another library}: a copy of the same type from an independent library of the
same culture, among copies whose culture has two or more libraries. \emph{Public record}: a UNITE+INSD record at
$\ge$99\% identity over $\ge$90\% coverage, culture copies excluded. \emph{Degraded}: broken base pairs above the 99th
percentile of main-type copies in the SSU, 5.8S or LSU (Methods S8). The same-length copy of \spn{Rhizophagus
intraradices} is counted with the others. None of these lines shows that a copy is transcribed.
\end{tablenotes}
\end{threeparttable}
\end{table}
'''
    with open(os.path.join(a.paper, 'tables', 'table_classes.tex'), 'w') as fh:
        fh.write(tex)
    tags = {'rhizo': 'Rh', 'septo': 'Sv', 'giga': 'Gr', 'other': 'Ot'}
    with open(os.path.join(a.paper, 'numbers_classes.tex'), 'w') as fh:
        fh.write('% written by analysis/divergent_classes.py; do not edit by hand\n')
        for k, t in tags.items():
            for n in ('copies', 'cultures', 'species', 'shared', 'shared_cross', 'recur', 'librep_yes', 'librep_tested', 'external', 'degraded',
                      'rec_culture', 'rec_species'):
                name = 'Cl' + t + ''.join(p.capitalize() for p in n.split('_'))
                fh.write(f"\\newcommand{{\\{name}}}{{{S[k][n]}}}\n")
        fh.write(f"\\newcommand{{\\ClAllLibrepTested}}{{{sum(S[k]['librep_tested'] for k in tags)}}}\n")
        fh.write(f"\\newcommand{{\\ClAllCopies}}{{{sum(S[k]['copies'] for k in tags)}}}\n")
    for k, t in CLASSES:
        print(k, {n: v for n, v in S[k].items()})


if __name__ == '__main__':
    main()
