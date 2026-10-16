#!/usr/bin/env python3
"""Uncertainty of the headline rates, with species and cultures as the resampling units (added after review, 1.5 and C).

The 903 copies come from 151 cultures of 37 named species, so copies are not independent. Each rate is given with
two cluster-bootstrap intervals: species resampled with replacement, and cultures resampled with replacement (10,000
replicates, percentile intervals, seed fixed). A rate is pooled over the resampled clusters (events / copies). The
number of species and cultures that contribute events is given beside each rate.

Sources (all in the repository):
  * withheld cultures and species, external databases intact: the frozen M2b tables
    (benchmarks/m4_external_holdout/frozen_m2b/comparison/epa_species_clades_physical_{culture,species}.copies.tsv),
    to which the frozen ITS distance rule (its_rule.py, 99th percentile, refitted per round on the training cultures)
    is applied here exactly as compare_its_rule.py applied it; the counts are checked against the paper's macros;
  * LSU-removal test: per-species counts from tables/tableS1_m4_folds.tex (species level only; culture level needs
    the per-copy aggregate, --m4);
  * LongGloDB: benchmarks/m4_external_holdout/runs_lgdb/06_scores/{A,B}_primary_{v1,v2}.copies.tsv;
  * EUKARYOME (with --euk-runs once the run is in): runs_euk/06_scores, taxa as the units.

  python3 analysis/intervals.py --paper . [--m4 <aggregate_v2>]
Writes analysis/intervals.tsv, numbers_intervals.tex and tables/table_intervals.tex.
"""
import argparse
import collections
import csv
import os
import re
import statistics
import sys

import edlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'stage3'))
from arms import is_named, read_fasta  # noqa: E402

B, SEED = 10000, 20261015
SHORT = {'species kept, placement': 'kept', 'species kept, placement + ITS rule': 'kept, ITS rule',
         'demoted by the ITS rule': 'demoted by ITS rule', 'false known, placement': 'false known',
         'false known, placement + ITS rule': 'false known, ITS rule',
         'false known, LSU intact, ITS rule': 'false known, intact, ITS rule',
         'false known, species removed': 'false known, removed', 'false known, species removed, ITS rule': 'false known, removed, ITS rule',
         'species kept, no ITS rule': 'kept', 'species kept, ITS rule': 'kept, ITS rule'}
COMP = os.path.join(ROOT, 'benchmarks/m4_external_holdout/frozen_m2b/comparison')
LGDB = os.path.join(ROOT, 'benchmarks/m4_external_holdout/runs_lgdb/06_scores')
REF = os.path.join(ROOT, 'reference/culture_ref_v4')
NOVEL_ROLES = ('novel species, genus in v4', 'novel lineage, genus not in v4')
DEV = {'Glomus chinense', 'Glomus rugosae'}


def rd(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def macro(paper, name):
    for f in os.listdir(paper):
        if f.startswith('numbers_') and f.endswith('.tex') and f != 'numbers_intervals.tex':
            m = re.search(r'\\newcommand\{\\' + name + r'\}\{([^}]*)\}', open(os.path.join(paper, f)).read())
            if m:
                return m.group(1)
    raise KeyError(name)


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode='HW', task='distance')['editDistance'] / len(q)


class ItsRule:
    """its_rule.py (frozen), pool 'all', reading the v4 reference directly."""
    def __init__(self):
        self.cop = {r['accession']: r for r in rd(os.path.join(REF, 'copies.tsv')) if r['status'] == 'retained'}
        cul = rd(os.path.join(REF, 'cultures.tsv'))
        self.species = {r['culture']: r['organism'] for r in cul if is_named(r['organism']) and r.get('name_status', 'ok') == 'ok'}
        self.cultures = frozenset(r['culture'] for r in cul)
        self.by_species = collections.defaultdict(set)
        for c, sp in self.species.items():
            self.by_species[sp].add(c)
        self.by_culture = collections.defaultdict(list)
        for x, r in self.cop.items():
            self.by_culture[r['culture']].append(x)
        its = read_fasta(os.path.join(REF, 'fasta', 'full_its.fasta'))
        ok = sorted(x for x in self.cop if x in its)
        self.mains = {x for x in ok if self.cop[x]['copy_class'] == 'main'}
        self.dmin = {'all': collections.defaultdict(dict), 'main': collections.defaultdict(dict)}
        for i, x in enumerate(ok):
            for y in ok[i + 1:]:
                d = dist(its[x], its[y])
                for q, r in ((x, y), (y, x)):
                    c = self.cop[r]['culture']
                    if d < self.dmin['all'][q].get(c, float('inf')):
                        self.dmin['all'][q][c] = d
                    if r in self.mains and d < self.dmin['main'][q].get(c, float('inf')):
                        self.dmin['main'][q][c] = d
        self._thr = {}

    def threshold(self, training):
        if training in self._thr:
            return self._thr[training]
        vals = []
        for sp, cs in self.by_species.items():
            tc = cs & training
            if len(tc) < 2:
                continue
            for c in tc:
                for y in self.by_culture[c]:
                    if y in self.mains:
                        ds = [self.dmin['main'][y][o] for o in tc - {c} if o in self.dmin['main'][y]]
                        if ds:
                            vals.append(min(ds))
        self._thr[training] = statistics.quantiles(vals, n=100)[98]
        return self._thr[training]

    def apply(self, rows):
        out = []
        for r in rows:
            r = dict(r)
            training = self.cultures - frozenset(c for c in r['held_cultures'].split(';') if c)
            t = self.threshold(training)
            own = self.cop[r['accession']]['culture']
            ds = [d for c, d in self.dmin['all'].get(r['accession'], {}).items() if c in training and c != own]
            d = min(ds) if ds else None
            r['its_rule_fired'] = int(d is not None and d > t and r['status'].startswith('known'))
            if r['its_rule_fired']:
                r['status'] = 'novel species, genus known'
            out.append(r)
        return out


def known(r):
    return r['status'].startswith('known')


def kept(r):
    return known(r) and r['organism'] in r['placement'].split(';')


class Boot:
    def __init__(self):
        self.rng = np.random.default_rng(SEED)

    def ci(self, rows, event, unit):
        """rows: dicts; event(row)->bool; unit: key. Pooled rate with a cluster-bootstrap percentile interval."""
        by = collections.defaultdict(lambda: [0, 0])
        for r in rows:
            by[r[unit]][0] += bool(event(r))
            by[r[unit]][1] += 1
        e = np.array([v[0] for v in by.values()], float)
        n = np.array([v[1] for v in by.values()], float)
        idx = self.rng.integers(0, len(e), size=(B, len(e)))
        rates = e[idx].sum(1) / n[idx].sum(1)
        lo, hi = np.percentile(rates, [2.5, 97.5])
        return 100 * lo, 100 * hi, len(by), int((e > 0).sum())


def species_rows_from_table(path):
    rows = []
    for line in open(path):
        m = re.match(r"\\emph\{(.+?)\}(\$\^\{\\(dagger|ddagger)\}\$)? & (\d+) & (\d+) & (\d+) & (\d+) & (\d+)\\\\", line.strip())
        if m and not m.group(3):
            sp = m.group(1).replace('{[}', '[').replace('{]}', ']')
            n, c1, c2, v1, v2 = map(int, m.groups()[3:])
            for i in range(n):          # one pseudo-row per copy; events assigned in order (species level only)
                rows.append({'organism': sp, 'ctl_v1': i < c1, 'ctl_v2': i < c2, 'v1': i < v1, 'v2': i < v2})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--paper', default='.')
    ap.add_argument('--m4', help='aggregate_v2 of the LSU-removal run, for culture-level intervals')
    a = ap.parse_args()
    boot = Boot()
    out = []

    def add(test, outcome, rows, event, units=('organism', 'culture'), check=None, note=''):
        k, n = sum(bool(event(r)) for r in rows), len(rows)
        if check is not None and str(k) != str(check):
            raise AssertionError(f'{test} / {outcome}: {k} here, {check} in the paper')
        row = {'test': test, 'outcome': outcome, 'events': k, 'copies': n, 'pct': round(100 * k / n, 1)}
        for u, tag in (('organism', 'species'), ('culture', 'culture')):
            if u in units:
                lo, hi, nu, ne = boot.ci(rows, event, u)
                row.update({f'{tag}_lo': round(lo, 1), f'{tag}_hi': round(hi, 1), f'{tag}_units': nu, f'{tag}_with_events': ne})
            else:
                row.update({f'{tag}_lo': '', f'{tag}_hi': '', f'{tag}_units': '', f'{tag}_with_events': ''})
        row['note'] = note
        out.append(row)
        print(f"  {test:28s} {outcome:38s} {k:4d}/{n:4d} {row['pct']:5.1f}%  species {row['species_lo']}-{row['species_hi']} "
              f"({row['species_with_events']} of {row['species_units']})  culture {row['culture_lo']}-{row['culture_hi']} "
              f"({row['culture_with_events']} of {row['culture_units']})", flush=True)
        return row

    # ---- withheld cultures / species, databases intact
    rule = ItsRule()
    cul_clade = rd(os.path.join(COMP, 'epa_species_clades_physical_culture.copies.tsv'))
    spe_clade = rd(os.path.join(COMP, 'epa_species_clades_physical_species.copies.tsv'))
    cul_its, spe_its = rule.apply(cul_clade), rule.apply(spe_clade)
    kn_c = [r for r in cul_clade if r['copy_class'] == 'main' and r['group'] == 'named, species elsewhere']
    kn_i = [r for r in cul_its if r['copy_class'] == 'main' and r['group'] == 'named, species elsewhere']
    if len(kn_i) != int(macro(a.paper, 'KnownN')):
        raise AssertionError('known-species copies differ from KnownN')
    add('withheld cultures', 'species kept, placement', kn_c, kept, check=macro(a.paper, 'KnownKeptClade'))
    add('withheld cultures', 'species kept, placement + ITS rule', kn_i, kept, check=macro(a.paper, 'KnownKeptITS'))
    add('withheld cultures', 'demoted by the ITS rule', kn_i, lambda r: r['its_rule_fired'], check=macro(a.paper, 'ITSDemotedKnown'))
    held_c = [r for r in spe_clade if r['copy_class'] == 'main' and r['group'] == 'named, species held out']
    held_i = [r for r in spe_its if r['copy_class'] == 'main' and r['group'] == 'named, species held out']
    add('withheld species', 'false known, placement', held_c, known, check=macro(a.paper, 'FKClade'))
    add('withheld species', 'false known, placement + ITS rule', held_i, known, check=macro(a.paper, 'FKITS'))

    # ---- LSU-removal test (test species)
    if a.m4:
        for key, label, mac in (('control_v2', 'false known, LSU intact, ITS rule', 'MfourControlFKvTwo'),
                                ('v1', 'false known, species removed', 'MfourFKvOne'),
                                ('v2', 'false known, species removed, ITS rule', 'MfourFKvTwo')):
            rows = [r for r in rd(os.path.join(a.m4, f'species.{key}.copies.tsv'))
                    if r['copy_class'] == 'main' and r['group'] == 'named, species held out'
                    and r['development_fold'] == '0' and r['external_status'] == 'removed']
            add('LSU-removal test', label, rows, known, check=macro(a.paper, mac))
    else:
        rows = species_rows_from_table(os.path.join(a.paper, 'tables', 'tableS1_m4_folds.tex'))
        if len(rows) != int(macro(a.paper, 'MfourHeldN')):
            raise AssertionError('LSU-removal copies differ from MfourHeldN')
        add('LSU-removal test', 'false known, LSU intact, ITS rule', rows, lambda r: r['ctl_v2'], units=('organism',))
        add('LSU-removal test', 'false known, species removed', rows, lambda r: r['v1'], units=('organism',),
            check=macro(a.paper, 'MfourFKvOne'))
        add('LSU-removal test', 'false known, species removed, ITS rule', rows, lambda r: r['v2'], units=('organism',),
            check=macro(a.paper, 'MfourFKvTwo'))

    # ---- LongGloDB
    t = {n: rd(os.path.join(LGDB, f'{n}.copies.tsv')) for n in ('A_primary_v1', 'A_primary_v2', 'B_primary_v1', 'B_primary_v2')}
    nov = [r for r in t['B_primary_v2'] if r['test_role'] in NOVEL_ROLES]
    add('LongGloDB', 'false known, species removed, ITS rule', nov, known, check=macro(a.paper, 'LgFK'))
    k1 = [r for r in t['A_primary_v1'] if r['test_role'] == 'known species']
    k2 = [r for r in t['A_primary_v2'] if r['test_role'] == 'known species']
    add('LongGloDB', 'species kept, no ITS rule', k1, kept, check=macro(a.paper, 'LgKeptvOne'))
    add('LongGloDB', 'species kept, ITS rule', k2, kept, check=macro(a.paper, 'LgKept'))
    add('LongGloDB', 'demoted by the ITS rule', k2, lambda r: r['its_rule_fired'] == '1', check=macro(a.paper, 'LgNSixDemoted'))

    with open(os.path.join(HERE, 'intervals.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(out)

    # ---- macros and table
    def ci_text(r, tag):
        return f"{r[f'{tag}_lo']:.1f}--{r[f'{tag}_hi']:.1f}\\%" if r[f'{tag}_lo'] != '' else '--'
    names = {('withheld cultures', 'species kept, placement + ITS rule'): 'CiKnownKept',
             ('withheld cultures', 'demoted by the ITS rule'): 'CiDemoted',
             ('withheld species', 'false known, placement + ITS rule'): 'CiFKITS',
             ('LSU-removal test', 'false known, LSU intact, ITS rule'): 'CiMfourCtl',
             ('LSU-removal test', 'false known, species removed'): 'CiMfourVOne',
             ('LSU-removal test', 'false known, species removed, ITS rule'): 'CiMfourVTwo',
             ('LongGloDB', 'false known, species removed, ITS rule'): 'CiLgFK',
             ('LongGloDB', 'species kept, no ITS rule'): 'CiLgKeptvOne',
             ('LongGloDB', 'species kept, ITS rule'): 'CiLgKept',
             ('LongGloDB', 'demoted by the ITS rule'): 'CiLgDemoted'}
    with open(os.path.join(a.paper, 'numbers_intervals.tex'), 'w') as fh:
        fh.write('% written by analysis/intervals.py; do not edit by hand\n')
        for r in out:
            k = names.get((r['test'], r['outcome']))
            if not k:
                continue
            fh.write(f"\\newcommand{{\\{k}Pct}}{{{r['pct']:.1f}\\%}}\n")
            fh.write(f"\\newcommand{{\\{k}Sp}}{{{ci_text(r, 'species')}}}\n")
            fh.write(f"\\newcommand{{\\{k}Cu}}{{{ci_text(r, 'culture')}}}\n")
            fh.write(f"\\newcommand{{\\{k}SpN}}{{{r['species_with_events']}}}\n")
            fh.write(f"\\newcommand{{\\{k}CuN}}{{{r['culture_with_events']}}}\n")
    body = []
    last = None
    for r in out:
        test = r['test'][0].upper() + r['test'][1:] if r['test'] != last else ''
        last = r['test']
        cu = f"{ci_text(r, 'culture')} ({r['culture_with_events']}/{r['culture_units']})" if r['culture_units'] != '' else '--'
        body.append(f"{test} & {SHORT.get(r['outcome'], r['outcome'])} & {r['events']}/{r['copies']} ({r['pct']:.1f}\\%) & "
                    f"{ci_text(r, 'species')} ({r['species_with_events']}/{r['species_units']}) & {cu}\\\\")
    tex = r'''\begin{table}[htbp]\centering
\begin{threeparttable}
\caption{Main error and retention rates, with 95\% intervals that resample species or cultures.}
\label{tab:intervals}
\footnotesize\setlength{\tabcolsep}{2pt}
\begin{tabular}{llrrr}
\toprule
Test & Outcome & Events & \multicolumn{2}{c}{95\% interval (units with events)}\\
 & & & resampling species & resampling cultures\\
\midrule
''' + '\n'.join(body) + r'''
\bottomrule
\end{tabular}
\begin{tablenotes}\footnotesize
\item \emph{Events}: events over copies (ASVs for LongGloDB); \emph{kept}: called a known species that contains their own;
\emph{intact}, \emph{removed}: the LSU reference intact, or with the species removed from it. Pooled rates with percentile intervals from 10,000 cluster-bootstrap replicates
that resample species or cultures with replacement. In brackets: species or cultures with at least one event, of those tested. Withheld cultures and species: external databases intact, LSU placement. LSU-removal test: the
test species; intervals by culture need the per-copy results and are not given. LongGloDB: known species with the LSU
reference intact; absent species with their species removed from it.
\end{tablenotes}
\end{threeparttable}
\end{table}
'''
    with open(os.path.join(a.paper, 'tables', 'table_intervals.tex'), 'w') as fh:
        fh.write(tex)


if __name__ == '__main__':
    sys.exit(main())
