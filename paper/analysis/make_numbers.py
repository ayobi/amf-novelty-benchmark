#!/usr/bin/env python3
"""Generate numbers_generated.tex and the benchmark tables from the run outputs.

Every number that comes from a benchmark or from the M4 external holdout is written here, so a
number in the text cannot drift from the run that produced it. Numbers with no run output behind
them (curation counts, crosswalk counts, calibration) stay in numbers_static.tex.

Inputs
  --its-run  DIR   output of compare_its_rule.py (summary.json and the *_physical_*.copies.tsv files)
  --m4       DIR   the aggregate_v2 directory of an M4 run (species.summary.json, species.*.copies.tsv;
                   genus.* is read too when present). After `external_holdout_v2.py collect` this is
                   results/aggregate_v2 inside the ZIP; before, runs/aggregate_v2.
  --paper    DIR   the paper directory; writes numbers_generated.tex and tables/*.tex there

Missing pieces (an unfinished M4 run, no genus folds yet) do not stop the run: the affected macros
are written as \\todo{...}, so the paper still compiles and the gap is visible in the PDF.

Writes
  numbers_generated.tex
  tables/table4_benchmark.tex, tables/table5_external.tex, tables/tableS1_m4_folds.tex

Usage
  python3 analysis/make_numbers.py --its-run ../amf_m2b/its_rule_comparison \\
      --m4 ../amf_m4_step1/runs/aggregate_v2 --paper .
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

DEV = ('Glomus chinense', 'Glomus rugosae')          # PROTOCOL_v2.md, item 3


def rd(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def pct(a, b, nd=1):
    return f'{100 * a / b:.{nd}f}\\%' if b else '--'


class Macros:
    def __init__(self):
        self.lines, self.todos = [], []

    def add(self, name, value, comment=None):
        if comment:
            self.lines.append(f'% {comment}')
        self.lines.append(f'\\newcommand{{\\{name}}}{{{value}}}')

    def todo(self, name, why):
        self.todos.append(name)
        self.lines.append(f'\\newcommand{{\\{name}}}{{\\todo{{{why}}}}}')


def is_fk(r):
    return r['copy_class'] == 'main' and r['group'] == 'named, species held out' and r['status'].startswith('known')


def is_held(r):
    return r['copy_class'] == 'main' and r['group'] == 'named, species held out'


def short_name(full):
    """'Diversispora aestuarii' -> '\\emph{D.\\ aestuarii}'; brackets are braced for LaTeX."""
    parts = []
    for one in full.split(';'):
        one = one.strip()
        bracket = one.startswith('[')
        toks = one.replace('[', '').replace(']', '').split()
        if len(toks) >= 2:
            g = toks[0][0] + '.'
            g = '{[}' + g + '{]}' if bracket else g
            parts.append('\\emph{' + g + '\\ ' + ' '.join(toks[1:]) + '}')
        else:
            parts.append('\\emph{' + one + '}')
    return ' / '.join(parts)


def final_exact(rows):
    """Known-species copies whose final call is exactly their own species."""
    return sum(r['copy_class'] == 'main' and r['group'] == 'named, species elsewhere'
               and r['status'] == 'known species' and r['placement'] == r['organism'] for r in rows)


# ------------------------------------------------------------------ culture-bridge benchmark
def benchmark(m, its_dir):
    s = json.loads((its_dir / 'summary.json').read_text())
    names = {'Base': 'baseline', 'Clade': 'epa_species_clades', 'ITS': 'epa_species_clades+its'}
    cul = {k: s[f'{p}_physical_culture'] for k, p in names.items()}
    spe = {k: s[f'{p}_physical_species'] for k, p in names.items()}
    known_n, held_n, single_n = cul['Base']['known_n'], spe['Base']['held_n'], cul['Base']['held_n']
    m.add('KnownN', known_n, 'culture-bridge benchmark (compare_its_rule.py)')
    m.add('HeldN', held_n)
    m.add('SingletonN', single_n)
    rows = {}
    for k, p in names.items():
        rows[k] = rd(its_dir / f'{p}_physical_culture.copies.tsv')
        ex = final_exact(rows[k])
        m.add(f'KnownKept{k}', cul[k]['known_compatible'])
        m.add(f'Exact{k}', ex)
        m.add(f'FK{k}', spe[k]['false_known'])
        m.add(f'FK{k}Pct', pct(spe[k]['false_known'], held_n))
        m.add(f'SingletonFK{k}', cul[k]['false_known'])
    m.add('ExactITSPct', pct(final_exact(rows['ITS']), known_n))
    m.add('FKITSSpeciesMacro', f"{100 * spe['ITS']['false_known_species_macro']:.1f}\\%")
    m.add('FKITSCultureMacro', f"{100 * spe['ITS']['false_known_culture_macro']:.1f}\\%")
    div = [r for r in rows['ITS'] if r['copy_class'] == 'divergent' and r['group'] == 'named, species elsewhere']
    m.add('DivN', len(div), 'divergent-class copies of species held by another culture (culture holdout, placement + N6)')
    m.add('DivRecognised', sum(r['status'] == 'divergent class' for r in div))
    m.add('DivKnown', sum(r['status'].startswith('known') for r in div))
    m.add('DivNovel', sum(r['status'].startswith('novel') for r in div))
    m.add('ShortRecCulture', cul['ITS']['short_recognised'])
    m.add('ShortRecSpecies', spe['ITS']['short_recognised'])
    # what each arm and the reconciliation give before the decision tree (baseline configuration)
    known = [r for r in rows['Base'] if r['copy_class'] == 'main' and r['group'] == 'named, species elsewhere']
    for arm in ('VT', 'SH', 'LSU'):
        n = sum(r[f'{arm}_species'] == r['organism'] for r in known)
        m.add(f'ArmExact{arm}', n, 'exact species from one arm alone, known-species copies, culture holdout' if arm == 'VT' else None)
        m.add(f'ArmExact{arm}Pct', pct(n, known_n, 0))
    recon = [r for r in known if r['reconcile_rank'] == 'species' and r['call_species'] == r['organism']]
    m.add('ReconExact', len(recon))
    m.add('ReconExactPct', pct(len(recon), known_n, 0))
    fin = {r['accession'] for r in known if r['status'] == 'known species' and r['placement'] == r['organism']}
    if not fin <= {r['accession'] for r in recon}:
        raise ValueError('an exact final call that the reconciliation did not make')
    m.add('TreeDemoted', len(recon) - len(fin))
    # the ITS rule on known-species copies (culture holdout) and what is left in the species holdout
    fired = [r for r in rows['ITS'] if r['copy_class'] == 'main' and r['group'] == 'named, species elsewhere' and r['its_rule_fired'] == '1']
    m.add('ITSDemotedKnown', len(fired), 'the ITS rule on known-species copies (culture holdout, placement + N6)')
    m.add('ITSDemotedKnownPct', pct(len(fired), known_n))
    by = Counter(r['culture'] for r in fired)
    m.add('ITSDemotedText', ', '.join(f"{n} from {short_name(c.split(' | ')[0])} {c.split(' | ')[1].replace('_', chr(92) + '_')}"
                                      for c, n in by.most_common()))
    lit = [float(r['its_nearest_training_pct']) for r in fired if r['organism'] == 'Microdominikia litorea']
    if lit:
        m.add('LitoreaITS', f'{lit[0]:.1f}\\%')
    thr = [float(r['threshold_pct']) for r in rd(its_dir / 'its_thresholds.tsv') if r['mode'] == 'physical_culture']
    typical = Counter(round(t, 1) for t in thr).most_common(1)[0]
    m.add('ITSThrTypical', f'{typical[0]:.1f}\\%')
    m.add('ITSThrTypicalFolds', f'{typical[1]} of {len(thr)}')
    left = Counter(frozenset((r['organism'], r['placement'])) for r in rd(its_dir / 'epa_species_clades+its_physical_species.copies.tsv') if is_fk(r))
    if sum(left.values()) != spe['ITS']['false_known']:
        raise ValueError('species-holdout false known do not match the summary')
    want = {frozenset(('Diversispora varaderana', 'Diversispora aestuarii')): 'FKITSVarAes',
            frozenset(('Glomus chinense', 'Glomus rugosae')): 'FKITSChiRug'}
    if set(left) - set(want):
        raise ValueError(f'remaining species-holdout errors outside the two pairs named in the text: {set(left) - set(want)}')
    for k, name in want.items():
        m.add(name, left.get(k, 0))
    b = {r['accession']: r for r in rd(its_dir / 'baseline_physical_species.copies.tsv')}
    c = {r['accession']: r for r in rd(its_dir / 'epa_species_clades_physical_species.copies.tsv')}
    loss = sum(is_held(b[x]) and b[x]['status'] == 'novel species, genus known' and c[x]['status'] == 'novel lineage'
               for x in b)
    m.add('GenusLoss', loss, 'held-out copies moving from "novel species, genus known" to "novel lineage" (placement)')

    rowsT = [('Nearest-reference LSU', 'Base'), ('LSU species-clade placement', 'Clade'), ('Placement + ITS rule (N6)', 'ITS')]
    body = ''.join(f"{label} & {cul[k]['known_compatible']} & {final_exact(rows[k])} & "
                   f"{spe[k]['false_known']} ({pct(spe[k]['false_known'], held_n)}) & {cul[k]['false_known']}\\\\\n"
                   for label, k in rowsT)
    return f'''\\begin{{table}}[htbp]\\centering
\\begin{{threeparttable}}
\\caption{{Training-only benchmark with the external databases fixed.}}
\\label{{tab:benchmark}}
\\footnotesize\\setlength{{\\tabcolsep}}{{4pt}}
\\begin{{tabular}}{{lrrrr}}
\\toprule
Configuration & Known kept & Exact species & False known, & False known,\\\\
 & (of {known_n}) & & all held out (of {held_n}) & single-culture (of {single_n})\\\\
\\midrule
{body}\\bottomrule
\\end{{tabular}}
\\begin{{tablenotes}}\\footnotesize
\\item \\emph{{Known kept}}: main-type copies of species held by another training
culture that receive a known species or species-group call containing the true
species (culture holdout, 148 folds). \\emph{{False known}}: known calls for copies
whose species has no training culture (species holdout, 50 folds). Averaged over
species and over cultures, the final false-known rate is {100 * spe['ITS']['false_known_species_macro']:.1f}\\% and {100 * spe['ITS']['false_known_culture_macro']:.1f}\\%.
\\end{{tablenotes}}
\\end{{threeparttable}}
\\end{{table}}
'''


# ------------------------------------------------------------------ M4 external holdout
def per_target(m4_dir, rank):
    """{target: {status, dev, held, fk_v1, fk_v2, ctl_v1, ctl_v2, ...}} from the four copy tables."""
    out = {}
    for key in ('v1', 'v2', 'control_v1', 'control_v2'):
        p = m4_dir / f'{rank}.{key}.copies.tsv'
        if not p.exists():
            return {}
        for r in rd(p):
            t = out.setdefault(r['target'], {'status': r['external_status'], 'dev': int(r['development_fold']),
                                              'held': {}, 'fk': {}})
            t['held'][key] = t['held'].get(key, 0) + is_held(r)
            t['fk'][key] = t['fk'].get(key, 0) + is_fk(r)
    return out


def external(m, m4_dir):
    summ_path = m4_dir / 'species.summary.json'
    if not summ_path.exists():
        for n, why in (('MfourFolds', 'M4 species folds'), ('MfourHeldN', 'held-out copies'), ('MfourFKvOne', 'v1 false known'),
                       ('MfourFKvTwo', 'v2 false known'), ('MfourFKvOnePct', 'v1 rate'), ('MfourFKvTwoPct', 'v2 rate'),
                       ('MfourControlFKvOne', 'control false known'), ('MfourControlFKvTwo', 'control false known'),
                       ('MfourControlPct', 'control rate'), ('MfourGenusCalled', 'genus calls'), ('MfourGenusCompatible', 'genus calls'),
                       ('MfourGenusN', 'genus calls')):
            m.todo(n, f'{why}: no species.summary.json in {m4_dir}')
        return None
    summ = json.loads(summ_path.read_text())
    pt = per_target(m4_dir, 'species')
    done = sorted(pt)
    first = [t for t in done if not pt[t]['dev'] and pt[t]['status'] == 'removed']
    absent = [t for t in done if not pt[t]['dev'] and pt[t]['status'] != 'removed']
    dev = [t for t in DEV if t in pt]
    complete = bool(summ.get('complete'))
    m.add('MfourComplete', '' if complete else f"\\todo{{M4 species run incomplete: {len(done)} of {summ.get('expected_folds')} folds}}",
          'empty when the M4 species run is complete')
    m.add('MfourFolds', len(first))
    m.add('MfourAllFolds', len(done))
    # development folds
    for tgt, short in (('Glomus rugosae', 'Rug'), ('Glomus chinense', 'Chi')):
        if tgt in pt:
            t = pt[tgt]
            m.add(f'{short}Full', t['fk']['control_v1'])
            m.add(f'{short}FullITS', t['fk']['control_v2'])
            m.add(f'{short}Removed', t['fk']['v1'])
            m.add(f'{short}RemovedITS', t['fk']['v2'])
            m.add(f'{short}N', t['held']['v1'])
            m.add(f'{short}Rescued', t['fk']['v1'] - t['fk']['v2'])
        else:
            for suf in ('Full', 'FullITS', 'Removed', 'RemovedITS', 'N', 'Rescued'):
                m.todo(f'{short}{suf}', f'{tgt} fold missing')
    # first-outcome folds, target removed from V18
    if first:
        agg = {k: sum(pt[t]['fk'][k] for t in first) for k in ('v1', 'v2', 'control_v1', 'control_v2')}
        held = sum(pt[t]['held']['v1'] for t in first)
        m.add('MfourHeldN', held)
        m.add('MfourFKvOne', agg['v1'])
        m.add('MfourFKvTwo', agg['v2'])
        m.add('MfourFKvOnePct', pct(agg['v1'], held))
        m.add('MfourFKvTwoPct', pct(agg['v2'], held))
        m.add('MfourControlFKvOne', agg['control_v1'])
        m.add('MfourControlFKvTwo', agg['control_v2'])
        m.add('MfourControlPct', pct(agg['control_v1'], held))
        m.add('MfourControlTwoPct', pct(agg['control_v2'], held))
        stratum = summ.get('first_outcomes/removed', {})
        for k, name in (('v1', 'MfourFKvOneSpeciesMacro'), ('v2', 'MfourFKvTwoSpeciesMacro')):
            v = stratum.get(k, {}).get('false_known_species_macro')
            m.add(name, f'{100 * v:.1f}\\%') if v is not None else m.todo(name, 'species-averaged rate')
        g = stratum.get('v2', {}).get('all_main')
        if g:
            m.add('MfourGenusN', g['n'])
            m.add('MfourGenusCalled', g['genus_called'])
            m.add('MfourGenusCompatible', g['genus_compatible'])
            m.add('MfourGenusWrong', g['genus_wrong'])
        else:
            for n in ('MfourGenusN', 'MfourGenusCalled', 'MfourGenusCompatible', 'MfourGenusWrong'):
                m.todo(n, 'genus-call counts')
    else:
        for n in ('MfourHeldN', 'MfourFKvOne', 'MfourFKvTwo', 'MfourFKvOnePct', 'MfourFKvTwoPct', 'MfourControlFKvOne',
                  'MfourControlFKvTwo', 'MfourControlPct', 'MfourFKvOneSpeciesMacro', 'MfourFKvTwoSpeciesMacro',
                  'MfourGenusN', 'MfourGenusCalled', 'MfourGenusCompatible', 'MfourGenusWrong'):
            m.todo(n, 'no first-outcome fold scored yet')
    # the fold whose species V18 does not hold by label
    if absent:
        t = pt[absent[0]]
        m.add('MfourAbsentTarget', absent[0].replace('[', '{[}').replace(']', '{]}'))
        m.add('MfourAbsentFK', f"{t['fk']['v1']}/{t['held']['v1']}")
        m.add('MfourAbsentFKvTwo', f"{t['fk']['v2']}/{t['held']['v2']}")
    else:
        m.todo('MfourAbsentTarget', 'fold not scored yet')
        m.todo('MfourAbsentFK', 'fold not scored yet')
        m.todo('MfourAbsentFKvTwo', 'fold not scored yet')

    # what the remaining errors are, and how the held-out copies are called, in the first-outcome removed folds
    if first:
        v2rows = [r for r in rd(m4_dir / 'species.v2.copies.tsv') if r['target'] in first and is_held(r)]
        left = {}
        for r in v2rows:
            if r['status'].startswith('known'):
                left[(r['organism'], r['placement'])] = left.get((r['organism'], r['placement']), 0) + 1
        top = sorted(left.items(), key=lambda kv: -kv[1])
        text = ', '.join(f"{short_name(o)} called {short_name(p)} ({n})" for (o, p), n in top[:6])
        m.add('MfourRemaining', text or 'none', 'the largest groups of remaining false-known calls (first-outcome folds, v2)')
        m.add('MfourRemainingPairs', len(top))
        mix = {}
        for r in v2rows:
            k = ('known' if r['status'].startswith('known') else 'novelspecies' if r['status'] == 'novel species, genus known'
                 else 'novellineage' if r['status'] == 'novel lineage' else 'other')
            mix[k] = mix.get(k, 0) + 1
        for k, name in (('novelspecies', 'MfourNovelSpecies'), ('novellineage', 'MfourNovelLineage'), ('known', 'MfourKnownCalls'), ('other', 'MfourOtherCalls')):
            m.add(name, mix.get(k, 0))
        # why the novel-species calls were made (reason codes of the decision tree)
        why = {}
        for r in v2rows:
            if r['status'] == 'novel species, genus known':
                for code in ('N1', 'N2', 'N3', 'N6'):
                    if code in r['reasons']:
                        why[code] = why.get(code, 0) + 1
        for code, name in (('N1', 'MfourNOne'), ('N2', 'MfourNTwo'), ('N3', 'MfourNThree'), ('N6', 'MfourNSix')):
            m.add(name, why.get(code, 0), 'demotions to "novel species, genus known" by reason code (first-outcome folds, v2)')
        m.add('MfourNOnePct', pct(why.get('N1', 0), len(v2rows)))
    else:
        for n in ('MfourRemaining', 'MfourRemainingPairs', 'MfourNovelSpecies', 'MfourNovelLineage', 'MfourKnownCalls', 'MfourOtherCalls',
                  'MfourNOne', 'MfourNTwo', 'MfourNThree', 'MfourNSix', 'MfourNOnePct', 'MfourControlTwoPct'):
            m.todo(n, 'no first-outcome fold scored yet')

    # table 5
    def cell(k, t):
        return f"{t['fk'][k]}/{t['held'][k]}"
    lines = []
    for tgt in dev:
        t, name = pt[tgt], '\\emph{' + tgt.replace('Glomus', 'G.\\') + '}'
        lines.append(f"{name}$^{{\\dagger}}$ & full & {cell('control_v1', t)} & {cell('control_v2', t)}\\\\")
        lines.append(f"{name}$^{{\\dagger}}$ & target removed & {cell('v1', t)} & {cell('v2', t)}\\\\")
    lines.append('\\midrule')
    if first:
        held = sum(pt[t]['held']['v1'] for t in first)
        lines.append(f"{len(first)} first-outcome folds & full & {agg['control_v1']}/{held} ({pct(agg['control_v1'], held)}) & "
                     f"{agg['control_v2']}/{held} ({pct(agg['control_v2'], held)})\\\\")
        lines.append(f"{len(first)} first-outcome folds & target removed & {agg['v1']}/{held} ({pct(agg['v1'], held)}) & "
                     f"{agg['v2']}/{held} ({pct(agg['v2'], held)})\\\\")
    else:
        lines.append('first-outcome folds & target removed & \\todo{v1} & \\todo{v2}\\\\')
    table5 = f'''\\begin{{table}}[htbp]\\centering
\\begin{{threeparttable}}
\\caption{{Removing the target species from the LSU reference as well.}}
\\label{{tab:external}}
\\small
\\begin{{tabular}}{{llrr}}
\\toprule
Target & V18 & False known (v1) & With N6 (v2)\\\\
\\midrule
{chr(10).join(lines)}
\\bottomrule
\\end{{tabular}}
\\begin{{tablenotes}}\\footnotesize
\\item False-known calls among the main-type copies of each fold's target species.
$^{{\\dagger}}$\\,Development folds: these two species informed the N6 rule. \\emph{{Full}}: the fold's
copies placed on the intact V18; \\emph{{target removed}}: V18 realigned without the target's records, model
refitted, all copies re-placed. First-outcome folds pool the species whose V18 records were removed{f'; {len(absent)} further fold(s) whose species V18 does not hold by label are reported apart' if absent else ''}.
\\end{{tablenotes}}
\\end{{threeparttable}}
\\end{{table}}
'''
    # supplementary per-fold table
    rows = []
    for t in done:
        d = pt[t]
        nm = t.replace('[', '{[}').replace(']', '{]}')
        mark = '$^{\\dagger}$' if d['dev'] else ('$^{\\ddagger}$' if d['status'] != 'removed' else '')
        rows.append(f"\\emph{{{nm}}}{mark} & {d['held']['v1']} & {d['fk']['control_v1']} & {d['fk']['control_v2']} & {d['fk']['v1']} & {d['fk']['v2']}\\\\")
    tableS3 = f'''\\begin{{longtable}}{{lrrrrr}}
\\caption{{External holdout by species fold: false-known calls among the target species' main-type copies.}}\\label{{tab:m4folds}}\\\\
\\toprule
Target species & Copies & Intact & Intact + N6 & Removed & Removed + N6\\\\
\\midrule
\\endfirsthead
\\toprule
Target species & Copies & Intact & Intact + N6 & Removed & Removed + N6\\\\
\\midrule
\\endhead
\\bottomrule
\\endfoot
{chr(10).join(rows)}
\\end{{longtable}}
\\noindent{{\\footnotesize Columns: false-known calls with V18 intact, V18 intact plus the ITS rule (N6), the target removed from V18, and removed plus N6. $^{{\\dagger}}$ development fold; $^{{\\ddagger}}$ species not in V18 by label (removal changes nothing).}}
'''
    return table5, tableS3


def complexes(m, prefix, m4_dir, ref_dir):
    """Species complexes (species_groups.py output): counts, thresholds and group-level false known."""
    pairs = rd(Path(f'{prefix}.pairs.tsv'))
    cx = {}
    for r in rd(Path(f'{prefix}.complexes.tsv')):
        for sp in r['species'].split(';'):
            cx[sp] = r['complex']
    n_cx = len(set(cx.values()))
    m.add('NComplexes', n_cx, 'species complexes (species_groups.py, distance rule on the frozen reference)')
    m.add('NComplexSpecies', len(cx))
    lim = {r['marker']: float(r['T']) for r in rd(Path(f'{prefix}.limits.tsv'))}
    for mk in ('ITS', 'LSU', 'SSU'):
        m.add(f'{mk}Spread', f'{lim[mk]:.1f}\\%')
    # the pair closest to being a complex: fails on one marker only, by the smallest margin
    near = []
    for r in pairs:
        if r['inseparable'] == '1':
            continue
        over = [(float(r[f'{mk}_min']) - lim[mk], mk) for mk in ('ITS', 'LSU', 'SSU') if float(r[f'{mk}_min']) > lim[mk]]
        if len(over) == 1:
            near.append((over[0][0], over[0][1], r['species_a'], r['species_b']))
    near.sort()
    if near:
        d, mk, sa, sb = near[0]
        m.add('NearMissPair', f'{short_name(sa)} with {short_name(sb)}')
        m.add('NearMissMarker', mk)
        m.add('NearMissPoints', f'{d:.1f}')
        m.add('NearMissNext', f'{near[1][0]:.1f}' if len(near) > 1 else '--')
    if ref_dir:
        cul = {r['culture']: r for r in rd(Path(ref_dir) / 'cultures.tsv')}
        main = [r for r in rd(Path(ref_dir) / 'copies.tsv') if r['status'] == 'retained' and r['copy_class'] == 'main']
        inc = sum(cul.get(r['culture'], {}).get('organism') in cx for r in main)
        m.add('ComplexCopies', inc)
        m.add('ComplexCopiesPct', pct(inc, len(main)))
    sm = Path(m4_dir) / 'species.summary.json'
    if not sm.exists():
        return None
    def same(r):
        sp, pl = r['organism'], [p for p in r['placement'].split(';') if p]
        return sp in cx and any(cx.get(p) == cx[sp] for p in pl)
    for key, tag in (('v1', 'VOne'), ('v2', 'VTwo'), ('control_v1', 'CtlOne'), ('control_v2', 'CtlTwo')):
        rows = [r for r in rd(Path(m4_dir) / f'species.{key}.copies.tsv')
                if is_held(r) and r['development_fold'] == '0' and r['external_status'] == 'removed']
        fk = [r for r in rows if r['status'].startswith('known')]
        inside = sum(same(r) for r in fk)
        m.add(f'FKInside{tag}', inside)
        m.add(f'FKOutside{tag}', len(fk) - inside)
        m.add(f'FKOutside{tag}Pct', pct(len(fk) - inside, len(rows)))
    # table of the inseparable pairs
    body = ''
    for r in pairs:
        if r['inseparable'] == '1':
            a, b = (x.replace('[', '{[}').replace(']', '{]}') for x in (r['species_a'], r['species_b']))
            body += (f"\\emph{{{a}}} & \\emph{{{b}}} & {cx[r['species_a']]} & {r['ITS_min']} & {r['LSU_min']} & {r['SSU_min']}\\\\\n")
    return f'''\\begin{{table}}[htbp]\\centering
\\begin{{threeparttable}}
\\caption{{Species complexes: pairs of named species whose nearest copies lie within the within-species spread on all three markers.}}
\\label{{tab:complexes}}
\\footnotesize\\setlength{{\\tabcolsep}}{{4pt}}
\\begin{{tabular}}{{llcrrr}}
\\toprule
Species A & Species B & Complex & ITS & LSU & SSU\\\\
\\midrule
{body}\\bottomrule
\\end{{tabular}}
\\begin{{tablenotes}}\\footnotesize
\\item Columns give the distance (\\%) between the nearest copies of the two species. The limits are the 99th
percentile of the distance between cultures of one species: ITS {lim['ITS']:.1f}\\%, LSU {lim['LSU']:.1f}\\%, SSU {lim['SSU']:.1f}\\%.
\\end{{tablenotes}}
\\end{{threeparttable}}
\\end{{table}}
'''


def genus_recall(m, m4_dir, ref_dir):
    """Genus recall: held-out copies whose genus is held by another culture, and how many kept it."""
    p = Path(m4_dir) / 'species.v2.copies.tsv'
    if not p.exists() or not ref_dir:
        return
    syn = {'Rhizoglomus': 'Rhizophagus'}
    genus = lambda o: syn.get(o.split()[0].strip('[]'), o.split()[0].strip('[]'))
    orgs = [r['organism'] for r in rd(Path(ref_dir) / 'cultures.tsv')]
    rows = [r for r in rd(p) if is_held(r) and r['development_fold'] == '0' and r['external_status'] == 'removed']
    rep_ = [r for r in rows if any(genus(o) == genus(r['organism']) and o != r['organism'] for o in orgs)]
    right = sum(r.get('final_genus_call') == '1' and r.get('final_genus_compatible') == '1' for r in rep_)
    m.add('GenusRepN', len(rep_), 'held-out copies whose genus is held by another culture (first-outcome folds)')
    m.add('GenusRepRight', right)
    m.add('GenusRecallPct', pct(right, len(rep_)))
    m.add('GenusRepLost', sum(r['status'] == 'novel lineage' for r in rep_))
    m.add('GenusUnrepN', len(rows) - len(rep_))


def genus_rank(m, m4_dir):
    """M4 genus folds: the whole genus withheld from the cultures and from V18."""
    p = m4_dir / 'genus.summary.json'
    names = ('GenusFolds', 'GenusMainN', 'GenusHeldN', 'GenusFKvOne', 'GenusFKvTwo', 'GenusNovelLineage', 'GenusWrongN',
             'GenusWrongPct', 'GenusWrongText', 'GenusControlWrong', 'GenusDivN', 'GenusDivRecognised', 'GenusDivWrong', 'GenusBoundary')
    if not p.exists():
        for n in names:
            m.todo(n, 'genus folds not run yet')
        return
    summ = json.loads(p.read_text())
    rows = rd(m4_dir / 'genus.v2.copies.tsv'); rows_v1 = rd(m4_dir / 'genus.v1.copies.tsv'); ctl = rd(m4_dir / 'genus.control_v2.copies.tsv')
    main = [r for r in rows if r['copy_class'] == 'main']
    held = [r for r in main if r['group'] == 'named, species held out']
    wrong = [r for r in main if r['final_genus_call'] == '1' and r['final_genus_compatible'] != '1']
    div = [r for r in rows if r['copy_class'] == 'divergent']
    m.add('GenusFolds', len({r['target'] for r in rows}), 'M4 genus folds (genus withheld from the cultures and from V18)')
    m.add('GenusMainN', len(main))
    m.add('GenusHeldN', len(held))
    m.add('GenusFKvOne', sum(is_fk(r) for r in rows_v1))
    m.add('GenusFKvTwo', sum(is_fk(r) for r in rows))
    m.add('GenusNovelLineage', sum(r['status'] == 'novel lineage' for r in held))
    m.add('GenusWrongN', len(wrong))
    if any(r['group'] != 'named, species held out' for r in wrong):
        raise ValueError('a wrong genus call outside the named held-out copies; revise the genus paragraph')
    m.add('GenusWrongPct', pct(len(wrong), len(held)))
    pairs = Counter((r['target'], r['placement']) for r in wrong)
    m.add('GenusWrongText', ', '.join(f"\\spn{{{t}}} called \\spn{{{c}}} ({n})" for (t, c), n in pairs.most_common()) or 'none')
    m.add('GenusControlWrong', sum(r['copy_class'] == 'main' and r['final_genus_call'] == '1' and r['final_genus_compatible'] != '1' for r in ctl))
    m.add('GenusDivN', len(div))
    m.add('GenusDivRecognised', sum(r['status'] == 'divergent class' for r in div))
    m.add('GenusDivWrong', sum(r['final_genus_call'] == '1' and r['final_genus_compatible'] != '1' for r in div))
    m.add('GenusBoundary', sum(r['status'] == 'boundary' for r in main))
    if not summ.get('complete'):
        m.todo('GenusIncomplete', f"genus run incomplete: {len(summ.get('completed_folds', []))} of {summ.get('expected_folds')} folds")


def percentiles(m, path, m4_dir):
    """N6 and the complexes moved to the same percentile (percentile_sensitivity.py)."""
    r = json.loads(Path(path).read_text())
    for p, tag in ((95.0, 'NinetyFive'), (97.5, 'NinetySevenHalf'), (99.0, 'NinetyNine')):
        x, c, sp = r[f'm4@{p}'], r[f'culture@{p}'], r[f'species@{p}']
        m.add(f'Pct{tag}FK', x['false_known'], f'N6 and complexes at the {p}th percentile' if p == 95.0 else None)
        m.add(f'Pct{tag}FKPct', pct(x['false_known'], x['held']))
        m.add(f'Pct{tag}Outside', x['outside_complex'])
        m.add(f'Pct{tag}Complexes', x['complexes'])
        m.add(f'Pct{tag}KnownDemoted', c['known_demoted'])
        m.add(f'Pct{tag}KnownDemotedPct', pct(c['known_demoted'], c['n']))
        m.add(f'Pct{tag}SpeciesFK', sp['false_known'])
        m.add(f'Pct{tag}SpeciesFKPct', pct(sp['false_known'], sp['n']))
        allowed = {'Rhizophagus irregularis -> [Rhizoglomus] vesiculiferum'}
        if set(x['outside_pairs']) - allowed:
            raise ValueError(f"errors outside a complex other than R. irregularis -> vesiculiferum at {p}: {x['outside_pairs']}")
    rv = r['m4@99.0']['pairs'].get('Rhizophagus irregularis -> [Rhizoglomus] vesiculiferum', 0)
    m.add('RirrVesN', rv, 'R. irregularis copies called [R.] vesiculiferum (first-outcome folds, with N6)')
    body = ''
    for p, lab in ((95.0, '95th'), (97.5, '97.5th'), (99.0, '99th')):
        x, c, sp = r[f'm4@{p}'], r[f'culture@{p}'], r[f'species@{p}']
        body += (f"{lab} & {x['median_T']:.1f}\\% & {x['complexes']} ({x['complex_species']}) & "
                 f"{c['known_demoted']} ({pct(c['known_demoted'], c['n'])}) & {sp['false_known']} ({pct(sp['false_known'], sp['n'])}) & "
                 f"{x['false_known']} ({pct(x['false_known'], x['held'])}) & {x['outside_complex']}\\\\\n")
    n_c, n_s, n_m = r['culture@99.0']['n'], r['species@99.0']['n'], r['m4@99.0']['held']
    return f'''\\begin{{table}}[htbp]\\centering
\\begin{{threeparttable}}
\\caption{{The ITS rule and the species complexes moved to the same percentile.}}
\\label{{tab:percentile}}
\\footnotesize\\setlength{{\\tabcolsep}}{{4pt}}
\\begin{{tabular}}{{lrrrrrr}}
\\toprule
Percentile & ITS limit & Complexes & Known demoted & False known, & False known, & Outside\\\\
 & (median) & (species) & (of {n_c}) & species holdout (of {n_s}) & V18 removed (of {n_m}) & a complex\\\\
\\midrule
{body}\\bottomrule
\\end{{tabular}}
\\begin{{tablenotes}}\\footnotesize
\\item The ITS rule is re-applied to the stored decision-tree rows at each percentile, with the limit fitted
on each fold's training cultures, and the complexes are those of the distance rule at the same percentile
(\\texttt{{percentile\\_sensitivity.py}}). \\emph{{Known demoted}}: known-species copies the rule demotes in the
culture holdout. \\emph{{V18 removed}}: first-outcome M4 folds. The 99th-percentile row reproduces the reported
results exactly. Every error outside a complex is \\emph{{R.\\ irregularis}} called \\emph{{{{[}}R.{{]}}\\ vesiculiferum}}.
\\end{{tablenotes}}
\\end{{threeparttable}}
\\end{{table}}
'''


def deposit(m, ref_dir):
    """What the GenBank deposit itself says (organism and culture fields), against the reference."""
    import re
    rows = rd(Path(ref_dir) / 'copies.tsv')
    unnamed = lambda o: bool(re.search(r'\bsp\.', o))
    gb = {r['organism_gb'] for r in rows}
    m.add('DepositLabels', len({r['culture_gb'].split(' | ')[-1] for r in rows}), 'GenBank deposit (organism and culture fields)')
    m.add('DepositNames', len(gb))
    m.add('DepositSpecies', sum(not unnamed(o) for o in gb))
    m.add('DepositGenusOnly', sum(unnamed(o) for o in gb))
    fin = {r['organism'] for r in rows if r['status'] == 'retained'}
    m.add('RefNames', len(fin))
    m.add('RefGenusOnly', sum(unnamed(o) for o in fin))


# ------------------------------------------------------------------ crosswalk (three naming systems)
GENUS_SYN = {'Rhizoglomus': 'Rhizophagus', 'Claroideoglomus': 'Entrophospora'}   # arms.py GENUS_SYNONYMS


def species_key(name):
    """Genus (synonyms joined) + epithet without its Latin ending, so that culture names meet V18 labels:
    'Rhizophagus prolifer' = 'Rhizophagus proliferus', '[Rhizoglomus] vesiculiferum' = 'Rhizophagus
    vesiculiferus', 'Septoglomus viscosum' = V18's clade-first 'Septoglomus Viscospora viscosa'."""
    import re
    t = name.replace('[', '').replace(']', '').split()
    if len(t) < 2:
        return name
    return GENUS_SYN.get(t[0], t[0]) + ' ' + re.sub(r'(us|um|a|e)$', '', t[-1])


def v18_labels(path):
    """Species labels of the V18 FASTA as arms.parse_delavaux reads them (all tokens after the family)."""
    import re
    out = set()
    for line in open(path):
        if line.startswith('>'):
            m = re.match(r'^([A-Z]{1,2}_?\d+(?:\.\d+)?)_(.+)$', line[1:].split()[0])
            toks = m.group(2).split('_')
            out.add(' '.join(toks[1:] if toks[0].endswith('aceae') else toks))
    return out


def tex_taxon(name):
    """UNITE/MaarjAM record name for a table: binomials in italics, placeholders ('Glomeraceae sp') upright."""
    t = name.replace('_', ' ').split()
    if len(t) >= 2 and t[1] in ('sp', 'sp.'):
        return (f'\\emph{{{t[0]}}} sp.' if not t[0].endswith(('aceae', 'mycota')) else f'{t[0]} sp.')
    return f'\\emph{{{name}}}' if len(t) == 2 else name


def crosswalk(m, prefix, ref_dir, v18_path):
    """Counts of section 3.3 and the divergent and paralog tables, from crosswalk.py tables (here the replay
    of the v4 crosswalk on the final reference v4, analysis/crosswalk_replay.py)."""
    from statistics import median
    cop = rd(Path(f'{prefix}.copies.tsv'))
    cul = rd(Path(f'{prefix}.cultures.tsv'))
    kind = {r['accession']: r.get('div_kind', '') for r in rd(Path(ref_dir) / 'copies.tsv')} if ref_dir else {}
    counted = [r for r in cul if r['name_status'] == 'ok']
    named = [r for r in counted if r['named'] == 'True']
    m.add('XwCultures', len(cul), 'crosswalk counts (crosswalk.py report on the final reference v4; analysis/crosswalk_replay.py)')
    m.add('XwNamedCultures', len(named))
    m.add('XwNamedSpecies', len({r['organism'] for r in named}))
    m.add('VTAssigned', sum(r['VT_call'] == 'assigned' for r in cop))
    sh_assigned = sum(r['SH_call'] == 'assigned' for r in cop)
    m.add('SHAssignedFull', sh_assigned)
    m.add('SHAssignedFullPct', pct(sh_assigned, len(cop), 0))
    m.add('LSUAssigned', sum(r['LSU_call'] == 'assigned' for r in cop))
    # many-to-many over named cultures, by modal label (crosswalk.py section 3)
    for arm in ('VT', 'SH', 'LSU'):
        su, us = {}, {}
        for r in named:
            if r[f'{arm}_modal']:
                su.setdefault(r['organism'], set()).add(r[f'{arm}_modal'])
                us.setdefault(r[f'{arm}_modal'], set()).add(r['organism'])
        m.add(f'{arm}SpeciesMulti', sum(len(v) > 1 for v in su.values()))
        m.add(f'{arm}SpeciesN', len(su))
        m.add(f'{arm}UnitsMulti', sum(len(v) > 1 for v in us.values()))
        m.add(f'{arm}UnitsN', len(us))
    # cross-arm splits over counted cultures (section 5)
    for a1, a2 in (('SH', 'VT'), ('SH', 'LSU'), ('LSU', 'VT'), ('VT', 'SH')):
        mm = {}
        for r in counted:
            if r[f'{a1}_modal'] and r[f'{a2}_modal']:
                mm.setdefault(r[f'{a1}_modal'], set()).add(r[f'{a2}_modal'])
        m.add(f'{a1}SplitBy{a2}', sum(len(v) > 1 for v in mm.values()))
        m.add(f'{a1}SplitBy{a2}N', len(mm))
    # within-culture SH splits of main-type copies (section 4)
    by = {}
    for r in cop:
        if r['copy_class'] == 'main' and r['SH_call'] == 'assigned':
            by.setdefault(r['culture'], []).append(r['SH_label'])
    two = {c: set(v) for c, v in by.items() if len(v) >= 2}
    split = [len(v) for v in two.values() if len(v) > 1]
    m.add('SHSplitCultures', len(split))
    m.add('SHSplitDenom', len(two))
    m.add('SHSplitMax', max(split))
    # LSU nearest reference: own species among named cultures whose species V18 holds
    if v18_path:
        v18 = {species_key(x) for x in v18_labels(v18_path)}
        held = [r for r in named if species_key(r['organism']) in v18]
        own = sum(bool(r['LSU_modal']) and species_key(r['LSU_modal']) == species_key(r['organism']) for r in held)
        none = sum(not r['LSU_modal'] for r in held)
        g = lambda s: GENUS_SYN.get(s.replace('[', '').split()[0].strip(']'), s.replace('[', '').split()[0].strip(']'))
        gen_ok = sum(g(r['LSU_modal']) == g(r['organism']) for r in held if r['LSU_modal'])
        m.add('LSUOwnSpecies', own)
        m.add('LSUOwnDenom', len(held))
        m.add('LSUOwnPct', pct(own, len(held), 0))
        m.add('LSUCongener', len(held) - own - none)
        m.add('LSUNoCall', none)
        m.add('LSUGenusRight', 'all' if gen_ok == len(held) - none else f'{gen_ok} of {len(held) - none}')
    # divergent copies against the culture's main-type unit (section 6), by 5.8S kind
    modal = {r['culture']: r for r in cul}
    main_id = {arm: median(float(r[f'{arm}_best_pident']) for r in cop if r['copy_class'] == 'main' and r[f'{arm}_best_pident'])
               for arm in ('VT', 'SH', 'LSU')}
    col = {}
    for k in ('short_5.8S', 'same_length'):
        xs = [r for r in cop if r['copy_class'] == 'divergent' and kind.get(r['accession']) == k]
        c = {'n': len(xs), 'cultures': len({r['culture'] for r in xs})}
        for arm in ('VT', 'SH', 'LSU'):
            same = other = none = nomodal = 0
            for r in xs:
                mo = modal[r['culture']][f'{arm}_modal']
                if not mo:
                    nomodal += 1
                elif r[f'{arm}_call'] != 'assigned':
                    none += 1
                elif r[f'{arm}_label'] == mo:
                    same += 1
                else:
                    other += 1
            ids = [float(r[f'{arm}_best_pident']) for r in xs if r[f'{arm}_best_pident']]
            c[arm] = (same, other, none, nomodal, median(ids) if ids else float('nan'))
        col[k] = c
    s, d = col['short_5.8S'], col['same_length']
    m.add('DivSHMedian', f"{s['SH'][4]:.1f}\\%", 'divergent copies (crosswalk.py section 6)')
    m.add('MainSHMedian', f"{main_id['SH']:.1f}\\%")
    m.add('SameLenN', d['n'])
    m.add('SameLenCultures', d['cultures'])

    def row(label, f):
        return f'{label} & {f(s)} & {f(d)}\\\\\n'
    body = (row('Copies (cultures)', lambda c: f"{c['n']} ({c['cultures']})")
            + row('VT: same as the culture\'s main type / other / no call', lambda c: f"{c['VT'][0]} / {c['VT'][1]} / {c['VT'][2]}")
            + row('LSU: same species / another species / no call', lambda c: f"{c['LSU'][0]} / {c['LSU'][1]} / {c['LSU'][2]}")
            + row(f"LSU best-hit identity, median (main type {main_id['LSU']:.1f}\\%)", lambda c: f"{c['LSU'][4]:.1f}\\%")
            + row('SH: same / other SH / no call', lambda c: f"{c['SH'][0]} / {c['SH'][1]} / {c['SH'][2]}")
            + row(f"SH best-hit identity, median (main type {main_id['SH']:.1f}\\%)", lambda c: f"{c['SH'][4]:.1f}\\%"))
    nm = {k: [c[a][3] for a in ('VT', 'LSU', 'SH')] for k, c in col.items()}
    table_div = f'''\\begin{{table}}[htbp]\\centering
\\begin{{threeparttable}}
\\caption{{The two divergent copy classes mislead each marker differently.}}
\\label{{tab:divergent}}
\\small
\\begin{{tabular}}{{p{{0.50\\linewidth}}p{{0.18\\linewidth}}p{{0.18\\linewidth}}}}
\\toprule
 & Short 5.8S & Same length\\\\
\\midrule
{body}\\bottomrule
\\end{{tabular}}
\\begin{{tablenotes}}\\footnotesize
\\item Each divergent copy against the modal unit of its culture's main-type copies (crosswalk.py,
culture reference v4). Copies in cultures without a main-type unit on a marker are left out of that
marker's row (short 5.8S: VT {nm['short_5.8S'][0]}, LSU {nm['short_5.8S'][1]}, SH {nm['short_5.8S'][2]};
same length: VT {nm['same_length'][0]}, LSU {nm['same_length'][1]}, SH {nm['same_length'][2]}).
\\end{{tablenotes}}
\\end{{threeparttable}}
\\end{{table}}
'''
    # paralog units: units reached by divergent copies and by at most one main-type copy
    rows_p, singles, counts = [], [], {}
    for arm in ('SH', 'VT'):
        mainc = Counter(r[f'{arm}_label'] for r in cop if r['copy_class'] == 'main' and r[f'{arm}_call'] == 'assigned')
        reach = {}
        for r in cop:
            if r['copy_class'] == 'divergent' and r[f'{arm}_call'] == 'assigned':
                reach.setdefault(r[f'{arm}_label'], []).append(r)
        units = {u: xs for u, xs in reach.items() if mainc[u] <= 1}
        counts[arm] = (len(reach), sum(mainc[u] == 0 for u in reach), sum(mainc[u] == 1 for u in reach),
                       sum(mainc[u] == 0 and 'sp' in Counter(r[f'{arm}_best_taxon'] for r in xs).most_common(1)[0][0].replace('_', ' ').replace('.', '').split()
                           for u, xs in reach.items()))
        for u, xs in sorted(units.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            ids = sorted(float(r[f'{arm}_pident']) for r in xs)
            idt = f'{ids[0]:.1f}\\%' if ids[0] == ids[-1] else f'{ids[0]:.1f}--{ids[-1]:.1f}\\%'
            orgs = Counter(r['organism'] for r in xs)
            who = '; '.join(f"{short_name(o) if 'sp.' not in o else chr(92) + 'emph{' + o.split()[0] + '} sp.'} {n}" for o, n in orgs.most_common())
            ncul = len({r['culture'] for r in xs})
            name = Counter(r[f'{arm}_best_taxon'] for r in xs).most_common(1)[0][0]
            name = tex_taxon(name) if arm == 'SH' else 'environmental (MaarjAM)'
            if len(xs) >= 2:
                rows_p.append(f"{u} & {name} & {who} ({ncul} culture{'s' if ncul > 1 else ''}) & {mainc[u]} & {idt}\\\\")
            else:
                singles.append(u.replace('.10FU', ''))
    m.add('DivSHs', counts['SH'][0], 'SHs reached by divergent-class copies; of them, reached by no / one main-type copy')
    m.add('ParalogSH', counts['SH'][1])
    m.add('DivSHOneMain', counts['SH'][2])
    m.add('ParalogPlaceholder', counts['SH'][3])
    table_par = f'''\\begin{{table}}[htbp]\\centering
\\begin{{threeparttable}}
\\caption{{Paralog units: SHs and a VT reached by divergent-class copies and by at most one main-type copy.}}
\\label{{tab:paralog}}
\\footnotesize\\setlength{{\\tabcolsep}}{{4pt}}
\\begin{{tabular}}{{p{{0.15\\linewidth}}p{{0.20\\linewidth}}p{{0.34\\linewidth}}p{{0.06\\linewidth}}p{{0.13\\linewidth}}}}
\\toprule
Unit & Name of the best hit & Divergent copies & Main & Identity\\\\
\\midrule
{chr(10).join(rows_p)}
\\bottomrule
\\end{{tabular}}
\\begin{{tablenotes}}\\footnotesize
\\item \\emph{{Main}}: main-type copies of any culture in the unit. {len(singles)} further units take one divergent
copy each ({', '.join(singles)}). Names are those of the best UNITE record; the VT row is MaarjAM's.
\\end{{tablenotes}}
\\end{{threeparttable}}
\\end{{table}}
'''
    return table_div, table_par


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--its-run', type=Path, required=True)
    ap.add_argument('--m4', type=Path, required=True)
    ap.add_argument('--paper', type=Path, default=Path('.'))
    ap.add_argument('--complexes', help='prefix of species_groups.py output (optional)')
    ap.add_argument('--ref', help='frozen culture reference directory (for the share of copies in complexes)')
    ap.add_argument('--crosswalk', help='prefix of the crosswalk tables on the final reference (crosswalk_replay.py output)')
    ap.add_argument('--percentiles', help='percentile_sensitivity.py output (JSON)')
    ap.add_argument('--v18', help='V18_LSUDB_052025_AMFONLY.fasta (species V18 holds, for the LSU own-species count)')
    a = ap.parse_args()
    m = Macros()
    table4 = benchmark(m, a.its_run)
    ext = external(m, a.m4)
    genus_rank(m, a.m4)
    cx_table = complexes(m, a.complexes, a.m4, a.ref) if a.complexes else None
    genus_recall(m, a.m4, a.ref)
    if a.percentiles:
        (a.paper / 'tables').mkdir(parents=True, exist_ok=True)
        (a.paper / 'tables/table_percentile.tex').write_text(percentiles(m, a.percentiles, a.m4))
    if a.ref:
        deposit(m, a.ref)
    if a.crosswalk:
        t_div, t_par = crosswalk(m, a.crosswalk, a.ref, a.v18)
        (a.paper / 'tables').mkdir(parents=True, exist_ok=True)
        (a.paper / 'tables/table2_divergent.tex').write_text(t_div)
        (a.paper / 'tables/table_paralog.tex').write_text(t_par)
    if cx_table:
        (a.paper / 'tables').mkdir(parents=True, exist_ok=True)
        (a.paper / 'tables/table_complexes.tex').write_text(cx_table)
    (a.paper / 'tables').mkdir(parents=True, exist_ok=True)
    (a.paper / 'tables/table4_benchmark.tex').write_text(table4)
    if ext:
        (a.paper / 'tables/table5_external.tex').write_text(ext[0])
        (a.paper / 'tables/tableS1_m4_folds.tex').write_text(ext[1])
    (a.paper / 'numbers_generated.tex').write_text(
        '% GENERATED by analysis/make_numbers.py; do not edit. See numbers_static.tex for the rest.\n'
        + '\n'.join(m.lines) + '\n')
    print(f'wrote {a.paper}/numbers_generated.tex ({len(m.lines)} lines), tables/table4_benchmark.tex'
          + (', table5_external.tex, tableS1_m4_folds.tex' if ext else ''))
    if m.todos:
        print(f'{len(m.todos)} macros are \\todo placeholders (unfinished inputs): ' + ', '.join(sorted(set(m.todos))[:8]) + ' ...')


if __name__ == '__main__':
    sys.exit(main())
