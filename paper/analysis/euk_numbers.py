#!/usr/bin/env python3
"""Generate numbers_euk.tex and tables/table_euk.tex from the EUKARYOME long-read test (PROTOCOL_EUK.md, note 1).

Two parts, like lgdb_numbers.py:
  test set   from the locked files in benchmarks/eukaryome (roles, predictions, removal plan, VT geometry,
             genus check) and the Stage 1 windows summary; available before the run.
  results    from the euk_test.py run (runs_euk: 06_scores, scores.json, 01_evidence/evidence_report.json) and
             runs_euk/report/report.json (euk_report.py). The main outcomes are recomputed here from the score
             tables with euk_report.py's definitions and must equal report.json; a difference stops the build.
A missing input does not stop the run: its macros become \\todo{...}, visible in the PDF, and the table is not
written.

Usage (from the paper directory)
  python3 analysis/euk_numbers.py --bench ../benchmarks/eukaryome [--runs PATH/amf_m4_step1/runs_euk]
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
import statistics
import sys

NOVEL = ('novel species, genus in v4', 'novel lineage, genus not in v4')
NOVEL_SP, NOVEL_LIN = NOVEL
KNOWN, GENUS_ONLY = 'known species', 'genus only'
STRATA = {'absent species, genus in v4': 'Sp', 'absent genus, named species': 'GenNamed',
          'absent genus, genus only': 'GenOnly', 'placeholder lineage': 'Placeholder'}
NO_GENUS_OUTCOME = {'Dentiscutata erythropus'}
NO_GENUS_OUTCOME_GENERA = {'Nanoglomus'}
MOTIF = re.compile(''.join({'Y': '[CT]', 'R': '[AG]'}.get(c, c) for c in 'YTYGGRGARTATG'))


def rd(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def pct(a, b, nd=1):
    return f'{100 * a / b:.{nd}f}\\%' if b else '--'


def mpct(x, nd=1):
    return f'{100 * x:.{nd}f}\\%' if x is not None else '--'


def num(n):
    return f'{n:,}'.replace(',', '{,}')


def spn(name):
    return f'\\spn{{{name}}}'


class Macros:
    def __init__(self):
        self.lines, self.todos = [], []

    def add(self, name, value, comment=None):
        if comment:
            self.lines.append(f'% {comment}')
        self.lines.append(f'\\newcommand{{\\{name}}}{{{value}}}')

    def todo(self, names, why):
        for n in names:
            self.lines.append(f'\\newcommand{{\\{n}}}{{\\todo{{{why}}}}}')
            self.todos.append(n)


RESULT_MACROS = (
    'EkFK', 'EkFKPct', 'EkFKTaxMacro', 'EkFKvOne', 'EkFKSec', 'EkFKSecPct', 'EkFKA', 'EkFKAvOne', 'EkFKOutsideLocked',
    'EkFKSecOutsideLocked', 'EkFKSp', 'EkFKGenNamed', 'EkFKGenOnly', 'EkFKPlaceholder', 'EkFKLockedCalled',
    'EkKept', 'EkKeptPct', 'EkKeptvOne', 'EkKeptvOnePct', 'EkNSixDemoted', 'EkPredInsideOwnKept', 'EkPredOutsideOwnKept',
    'EkKnownOther', 'EkGenusRight', 'EkGenusWrong', 'EkGenusNone', 'EkNovelSpAsLineage', 'EkLinScored', 'EkLinGenusGiven',
    'EkLinAsLineage', 'EkCallsGenusOnly', 'EkCallsGenusOnlyOutside', 'EkCallsKnownOutside', 'EkVTAnchor',
    'EkVTAnchorSpanned', 'EkSHAnchor', 'EkParityCopies', 'EkParityFolds', 'EkMotifDOne', 'EkMotifNovelLineage',
    'EkFKFlagged')


def is_known(r):
    return r['status'].startswith('known')


def called(r):
    return [s for s in r['placement'].split(';') if s]


def test_set(m, bench):
    """Macros that the locked test set fixes before the run."""
    asv = rd(bench / 'overlap.asvs.tsv')
    pred = {r['asv']: r for r in rd(bench / 'predictions.predictions.tsv')}
    test = [r for r in asv if r['role'] in NOVEL + (KNOWN, GENUS_ONLY)]
    ex = [r for r in asv if r['role'] == 'excluded']
    summ = json.loads((bench / 'stage1' / 'euk_test.summary.json').read_text())
    by = lambda role: [r for r in test if r['role'] == role]  # noqa: E731
    lin = by(NOVEL_LIN)
    m.add('EkReads', num(summ['reads']), 'EUKARYOME test set (eukaryome_overlap.py roles; PROTOCOL_EUK.md section 1)')
    m.add('EkWindows', num(len(asv)))
    idn = [r for r in ex if r['exclusion'] == 'identical in full ITS to a v4 copy']
    m.add('EkIdenticalN', num(len(idn)))
    m.add('EkIdenticalAllThree', num(sum(r['identical_v4_all3'] == 'yes' for r in idn)))
    m.add('EkIdenticalAssembly', num(sum(r['non_euk_id'] == 'yes' for r in idn)))
    m.add('EkOtherRemovedN', num(len(ex) - len(idn)))
    m.add('EkRemovedN', num(len(ex)))
    m.add('EkTestN', num(len(test)))
    m.add('EkKnownN', num(len(by(KNOWN))))
    m.add('EkKnownSpecies', num(len({r['organism'] for r in by(KNOWN)})))
    m.add('EkNovelSpN', num(len(by(NOVEL_SP))))
    m.add('EkNovelSpSpecies', num(len({r['organism'] for r in by(NOVEL_SP)})))
    m.add('EkNovelLinN', num(len(lin)))
    m.add('EkNovelLinTaxa', num(len({r['organism'] for r in lin})))
    named = [r for r in lin if r['stratum'] != 'placeholder lineage']
    m.add('EkNovelLinGenera', num(len({r['organism'].split()[0] for r in named})))
    m.add('EkNovelLinNamedN', num(len(named)))
    ph = [r for r in lin if r['stratum'] == 'placeholder lineage']
    m.add('EkPlaceholderN', num(len(ph)))
    m.add('EkPlaceholderTaxa', num(len({r['organism'] for r in ph})))
    m.add('EkNovelN', num(len(by(NOVEL_SP)) + len(lin)))
    go = by(GENUS_ONLY)
    m.add('EkGenusOnlyN', num(len(go)))
    m.add('EkGenusOnlyOpenN', num(sum(r['stratum'] == 'genus only, membership open' for r in go)))
    m.add('EkNonEukIDs', num(sum(r['non_euk_id'] == 'yes' for r in test)))
    m.add('EkVEighteenNear', num(sum(r['v18_best_pct'] not in ('', None) and float(r['v18_best_pct']) >= 99.5 for r in test)))
    m.add('EkLgdbIdentical', num(sum(bool(r['identical_lgdb']) for r in asv)))
    # predictions
    nov = [r for r in test if r['role'] in NOVEL]
    locked = [r for r in nov if pred[r['asv']]['inside']]
    m.add('EkLockedN', num(len(locked)), 'locked predictions (predict_longglodb.py, T at the paper\'s definition)')
    m.add('EkLockedCase', '; '.join(f"{spn(r['organism'])} inside {' and '.join(spn(s) for s in pred[r['asv']]['inside'].split('; '))}"
                                    for r in locked) or 'none')
    kn = by(KNOWN)
    own = [r for r in kn if pred[r['asv']]['own_species_in_v4'] in pred[r['asv']]['inside'].split('; ')]
    m.add('EkPredInsideOwn', num(len(own)))
    m.add('EkPredOutsideOwn', num(len(kn) - len(own)))
    m.add('EkGenusOnlyInside', num(sum(bool(pred[r['asv']]['inside']) for r in go)))
    lock = json.loads((bench / 'predictions.lock').read_text())
    for k, name in (('ITS', 'EkTITS'), ('LSU', 'EkTLSU'), ('SSU', 'EkTSSU')):
        m.add(name, f"{lock['T'][k]:.2f}\\%")
    # removal plan
    plan = rd(bench / 'overlap.v18_plan.tsv')
    apart = rd(bench / 'overlap.v18_apart.tsv')
    m.add('EkRoundsSp', num(sum(p['round'] == 'species' for p in plan)), 'condition B plan (overlap.v18_plan.tsv)')
    m.add('EkRoundsGen', num(sum(p['round'] == 'genus' for p in plan)))
    m.add('EkRounds', num(len(plan)))
    m.add('EkApartN', num(sum(int(p['reads']) for p in apart)))
    m.add('EkApartSpecies', num(sum(p['round'] == 'species' or (p['target'].count(' ') == 1 and p['round'] != 'none') for p in apart)))
    m.add('EkApartGenera', num(sum(p['round'] != 'none' and p['target'].count(' ') == 0 for p in apart)))
    # VT geometry
    off = rd(bench / 'vtgeom.offsets.tsv')
    top = max((r for r in off if r['offset_nt'] != 'not found'), key=lambda r: int(r['reads']))
    m.add('EkVTStart', f"{top['offset_nt']}~nt", 'VT geometry (eukaryome_overlap.py vtgeom)')
    m.add('EkVTStartN', num(int(top['reads'])))
    m.add('EkVTFull', num(sum(int(r['reads']) for r in off if r['offset_nt'] == '0')))
    mj = bench / 'vtgeom.maarjam.tsv'
    if mj.exists():
        g = {r['offset_nt']: r for r in rd(mj)}
        m.add('EkVTCoverPct', mpct(float(g[top['offset_nt']]['share'])))
    else:
        m.todo(['EkVTCoverPct'], 'vtgeom --maarjam')
    # genus check one rank up (note 1: distinct reads)
    gi = [r for r in rd(bench / 'genus_inside.reads.tsv') if r['inside_genus_only_culture']]
    m.add('EkNanoInsideN', num(len(gi)), 'genus_inside (PROTOCOL_EUK.md section 2, note 1)')
    m.add('EkNanoN', num(sum(r['organism'].split()[0] == 'Nanoglomus' for r in lin)))
    return test, pred


def motif_set(m, bench, wins_s58):
    """Short-5.8S motif in the test reads (descriptive, protocol section 6); needs the s58 test windows."""
    if not wins_s58 or not Path(wins_s58).exists():
        m.todo(['EkMotifN', 'EkMotifNano', 'EkMotifRhiz', 'EkMotifRhizN'], 'pass --s58 overlap.win.s58.fasta')
        return None
    seqs, name = {}, None
    for line in open(wins_s58):
        if line.startswith('>'):
            name = line[1:].strip()
            seqs[name] = ''
        else:
            seqs[name] += line.strip()
    asv = {r['asv']: r for r in rd(bench / 'overlap.asvs.tsv')}
    hit = {x for x, s in seqs.items() if MOTIF.search(s)}
    gen = lambda x: asv[x]['organism'].split()[0]  # noqa: E731
    rh = [x for x in seqs if gen(x) in ('Rhizophagus', 'Rhizoglomus')]
    m.add('EkMotifN', num(len(hit)), 'short-5.8S motif YTYGGRGARTATG in the test reads (descriptive)')
    m.add('EkMotifNano', num(sum(gen(x) == 'Nanoglomus' for x in hit)))
    m.add('EkMotifRhiz', num(sum(x in hit for x in rh)))
    m.add('EkMotifRhizN', num(len(rh)))
    return hit


def results(m, runs, test, pred, motif_hits, paper):
    """Macros and the table from the run; the main outcomes must equal report.json."""
    sc = runs / '06_scores'
    t = {p.name.replace('.copies.tsv', ''): rd(p) for p in sc.glob('*.copies.tsv')}
    meta = json.loads((sc / 'scores.json').read_text())
    rep = json.loads((runs / 'report' / 'report.json').read_text())
    ev = json.loads((runs / '01_evidence' / 'evidence_report.json').read_text())
    strat = {r['asv']: r['stratum'] for r in test}
    locked = {x for x, p in pred.items() if p['role'] in NOVEL and p['inside']}

    def primary_set(analysis, ver):
        a = {r['accession']: r for r in t[f'A_{analysis}_{ver}']}
        b = {r['accession']: r for r in t.get(f'B_{analysis}_{ver}', [])}
        return [b.get(x, a[x]) for x in sorted(a) if a[x]['test_role'] in NOVEL]

    ps = {(an, v): primary_set(an, v) for an in ('primary', 'secondary_vt') for v in ('v2', 'v1')}
    fk = {k: [r for r in v if is_known(r)] for k, v in ps.items()}
    p2 = ps[('primary', 'v2')]
    # cross-check with report.json
    want = rep['primary_primary_v2']
    if (len(p2), len(fk[('primary', 'v2')])) != (want['reads'], want['false_known']):
        sys.exit(f"primary outcome differs from report.json: {len(fk[('primary', 'v2')])}/{len(p2)} vs "
                 f"{want['false_known']}/{want['reads']}")
    by = defaultdict(list)
    for r in p2:
        by[r['organism']].append(is_known(r))
    macro = statistics.mean(statistics.mean(v) for v in by.values())
    m.add('EkFK', num(len(fk[('primary', 'v2')])), 'results (euk_test.py run, euk_report.py)')
    m.add('EkFKPct', pct(len(fk[('primary', 'v2')]), len(p2)))
    m.add('EkFKTaxMacro', mpct(macro))
    m.add('EkFKvOne', num(len(fk[('primary', 'v1')])))
    m.add('EkFKSec', num(len(fk[('secondary_vt', 'v2')])))
    m.add('EkFKSecPct', pct(len(fk[('secondary_vt', 'v2')]), len(ps[('secondary_vt', 'v2')])))
    a2 = [r for r in t['A_primary_v2'] if r['test_role'] in NOVEL]
    a1 = [r for r in t['A_primary_v1'] if r['test_role'] in NOVEL]
    m.add('EkFKA', num(sum(map(is_known, a2))))
    m.add('EkFKAvOne', num(sum(map(is_known, a1))))
    m.add('EkFKOutsideLocked', num(sum(r['accession'] not in locked for r in fk[('primary', 'v2')])))
    m.add('EkFKSecOutsideLocked', num(sum(r['accession'] not in locked for r in fk[('secondary_vt', 'v2')])))
    m.add('EkFKLockedCalled', num(sum(r['accession'] in locked for r in fk[('primary', 'v2')])))
    for s, short in STRATA.items():
        m.add(f'EkFK{short}', num(sum(strat[r['accession']] == s for r in fk[('primary', 'v2')])))
    flagged = {r['asv'] for r in test if r['non_euk_id'] == 'yes' or
               (r['v18_best_pct'] not in ('', None) and float(r['v18_best_pct']) >= 99.5)}
    m.add('EkFKFlagged', num(sum(r['accession'] in flagged for r in fk[('primary', 'v2')])))
    # known species kept (A)
    for ver, suffix in (('v2', ''), ('v1', 'vOne')):
        kn = [r for r in t[f'A_primary_{ver}'] if r['test_role'] == KNOWN]
        kept = [r for r in kn if is_known(r) and r['organism'] in called(r)]
        m.add(f'EkKept{suffix}', num(len(kept)))
        m.add(f'EkKept{suffix}Pct', pct(len(kept), len(kn)))
        if ver == 'v2':
            want = rep['known_retention_A']['v2']['all']
            if len(kept) != want['kept']:
                sys.exit(f"known kept differs from report.json: {len(kept)} vs {want['kept']}")
            ins = lambda r: pred[r['accession']]['own_species_in_v4'] in pred[r['accession']]['inside'].split('; ')  # noqa: E731
            m.add('EkPredInsideOwnKept', num(sum(ins(r) for r in kept)))
            m.add('EkPredOutsideOwnKept', num(sum(not ins(r) for r in kept)))
            m.add('EkKnownOther', num(sum(is_known(r) and r['organism'] not in called(r) for r in kn)))
            m.add('EkNSixDemoted', num(sum(r['its_rule_fired'] == '1' for r in kn)))
    # genus (primary set)
    sp = [r for r in p2 if strat[r['accession']] == 'absent species, genus in v4' and r['organism'] not in NO_GENUS_OUTCOME]
    m.add('EkGenusRight', num(sum(r['final_genus_call'] == '1' and r['final_genus_compatible'] == '1' for r in sp)))
    m.add('EkGenusWrong', num(sum(r['final_genus_call'] == '1' and r['final_genus_compatible'] != '1' for r in sp)))
    m.add('EkGenusNone', num(sum(r['final_genus_call'] != '1' for r in sp)))
    m.add('EkNovelSpAsLineage', num(sum(r['status'] == 'novel lineage' for r in p2 if r['test_role'] == NOVEL_SP)))
    lin = [r for r in p2 if r['test_role'] == NOVEL_LIN and r['organism'].split()[0] not in NO_GENUS_OUTCOME_GENERA]
    m.add('EkLinScored', num(len(lin)))
    m.add('EkLinGenusGiven', num(sum(r['final_genus_call'] == '1' for r in lin)))
    m.add('EkLinAsLineage', num(sum(r['status'] == 'novel lineage' for r in p2 if r['test_role'] == NOVEL_LIN)))
    # all known calls against the predictions (A, primary, N6)
    inside = {x: set(filter(None, p['inside'].split('; '))) for x, p in pred.items()}
    a_all = t['A_primary_v2']
    out = [r for r in a_all if is_known(r) and not set(called(r)) & inside[r['accession']]]
    m.add('EkCallsKnownOutside', num(len(out)))
    gok = [r for r in a_all if r['test_role'] == GENUS_ONLY and is_known(r)]
    m.add('EkCallsGenusOnly', num(len(gok)))
    m.add('EkCallsGenusOnlyOutside', num(sum(r in out for r in gok)))
    # evidence and integrity
    m.add('EkVTAnchor', num(ev['VT']['test_with_anchor']))
    m.add('EkVTAnchorSpanned', num(ev['VT_spanned']['test_with_anchor']))
    m.add('EkSHAnchor', num(ev['SH']['test_with_anchor']))
    m.add('EkParityCopies', num(meta['parity_m4']['copies']))
    m.add('EkParityFolds', num(meta['parity_m4']['folds']))
    # motif reads: divergent-copy calls
    if motif_hits is not None:
        mt = [r for r in a_all if r['accession'] in motif_hits]
        m.add('EkMotifDOne', num(sum(r['d1_short'] == '1' for r in mt)))
        m.add('EkMotifNovelLineage', num(sum(r['status'] == 'novel lineage' for r in mt)))
    else:
        m.todo(['EkMotifDOne', 'EkMotifNovelLineage'], 'pass --s58')
    # table, same columns as table_lgdb.tex
    nov_n = len(p2)
    kn_n = sum(1 for r in t['A_primary_v2'] if r['test_role'] == KNOWN)
    rows = []
    for analysis, label in (('primary', 'Main'), ('secondary_vt', 'Second')):
        for cond, rows_c in (('removed', ps[(analysis, 'v2')]), ('intact', [r for r in t[f'A_{analysis}_v2'] if r['test_role'] in NOVEL])):
            v1 = ps[(analysis, 'v1')] if cond == 'removed' else [r for r in t[f'A_{analysis}_v1'] if r['test_role'] in NOVEL]
            f2, f1 = sum(map(is_known, rows_c)), sum(map(is_known, v1))
            spc = [r for r in rows_c if strat[r['accession']] == 'absent species, genus in v4' and r['organism'] not in NO_GENUS_OUTCOME]
            gr = sum(r['final_genus_call'] == '1' and r['final_genus_compatible'] == '1' for r in spc)
            gw = sum(r['final_genus_call'] == '1' and r['final_genus_compatible'] != '1' for r in spc)
            gn = sum(r['final_genus_call'] != '1' for r in spc)
            lc = [r for r in rows_c if r['test_role'] == NOVEL_LIN and r['organism'].split()[0] not in NO_GENUS_OUTCOME_GENERA]
            lg = sum(r['final_genus_call'] == '1' for r in lc)
            if cond == 'intact':
                kn = [r for r in t[f'A_{analysis}_v2'] if r['test_role'] == KNOWN]
                kn1 = [r for r in t[f'A_{analysis}_v1'] if r['test_role'] == KNOWN]
                kk = f"{sum(is_known(r) and r['organism'] in called(r) for r in kn)} ({sum(is_known(r) and r['organism'] in called(r) for r in kn1)})"
            else:
                kk = '--'
            rows.append(f"{label if cond == 'removed' else ''} & {cond} & {num(f2)} ({num(f1)}) & {pct(f2, len(rows_c))} & "
                        f"{gr} / {gw} / {gn} & {num(lg)} & {kk}\\\\")
    tex = [r'\begin{table}[htbp]\centering', r'\begin{threeparttable}',
           r'\caption{Field long reads: the rules as fixed, on \EkTestN{} EUKARYOME reads.}', r'\label{tab:euk}',
           r'\footnotesize\setlength{\tabcolsep}{4pt}',
           r'\begin{tabular}{L{0.10\linewidth}L{0.12\linewidth}L{0.12\linewidth}L{0.07\linewidth}L{0.19\linewidth}L{0.14\linewidth}L{0.11\linewidth}}',
           r'\toprule',
           f"Analysis & LSU reference & False known (of {num(nov_n)}) & Share & Absent species: genus right/wrong/none (of {len(sp)}) & "
           f"Absent genera: genus given (of {num(len(lin))}) & Known kept (of {kn_n})\\\\", r'\midrule'] + rows + [
           r'\bottomrule', r'\end{tabular}', r'\begin{tablenotes}\footnotesize',
           r'\item With the ITS distance rule; in brackets, the same without it. \emph{LSU reference}: \emph{removed}, each absent '
           r'species or genus that V18 holds removed in its own round, and intact for the reads of taxa it does not name; '
           r'\emph{intact}, every read placed on the full reference. \emph{False known}: reads of taxa the reference does not '
           r'hold called a known species or species group. \emph{Absent species}: reads of species the reference does not hold, '
           r'in genera it holds (\spn{Dentiscutata erythropus}, whose genus is open, not counted). \emph{Absent genera}: reads of '
           r'genera the reference does not hold and of EUKARYOME placeholder lineages, so any genus given is wrong '
           r'(\spn{Nanoglomus} not counted, Section~\ref{sec:euk-methods}). \emph{Known kept}: reads of species the reference '
           r'holds whose call contains their species. \emph{Second}: VT coverage counted over the part of each MaarjAM '
           r'reference that the read spans.',
           r'\end{tablenotes}', r'\end{threeparttable}', r'\end{table}']
    (paper / 'tables' / 'table_euk.tex').write_text('\n'.join(tex) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bench', type=Path, default=Path('../benchmarks/eukaryome'))
    ap.add_argument('--runs', type=Path, help='runs_euk of euk_test.py (06_scores and report/report.json)')
    ap.add_argument('--s58', help='overlap.win.s58.fasta of the roles run (motif counts; sequences are not in git)')
    ap.add_argument('--paper', type=Path, default=Path('.'))
    a = ap.parse_args()
    m = Macros()
    m.lines.append('% GENERATED by analysis/euk_numbers.py from the EUKARYOME test; do not edit.')
    test, pred = test_set(m, a.bench)
    hits = motif_set(m, a.bench, a.s58)
    if a.runs and (a.runs / '06_scores' / 'scores.json').exists() and (a.runs / 'report' / 'report.json').exists():
        results(m, a.runs, test, pred, hits, a.paper)
    else:
        m.todo(RESULT_MACROS, 'EUKARYOME run')
    (a.paper / 'numbers_euk.tex').write_text('\n'.join(m.lines) + '\n')
    print(f"numbers_euk.tex: {sum(l.startswith(chr(92) + 'newcommand') for l in m.lines)} macros, {len(m.todos)} still TODO"
          + ('' if m.todos else '; tables/table_euk.tex written'))


if __name__ == '__main__':
    main()
