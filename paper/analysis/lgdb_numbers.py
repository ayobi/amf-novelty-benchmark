#!/usr/bin/env python3
"""Generate numbers_lgdb.tex and tables/table_lgdb.tex from the LongGloDB test (PROTOCOL_LGDB.md,
amendments 1 and 2).

Every LongGloDB number in the paper is written here from the run outputs, as make_numbers.py does for
the benchmarks; nothing is fitted. Kept apart from make_numbers.py so the two can be rerun
independently. A missing input does not stop the run: its macros become \\todo{...}, visible in the PDF.

Inputs
  --runs   DIR   the lgdb_test.py output (runs_lgdb): 06_scores/*.copies.tsv, 06_scores/scores.json,
                 plan.json, 01_evidence/evidence_report.json, A/parity.json (optional)
  --bench  DIR   bench/longglodb: overlap.asvs.tsv and overlap.cultures.tsv (longglodb_overlap.py, current
                 version) and predictions.predictions.tsv (the locked predictions)
  --paper  DIR   the paper directory

Usage (from the paper directory)
  python3 analysis/lgdb_numbers.py --runs ../amf_m4_step1/runs_lgdb --bench ../bench/longglodb --paper .
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
STEFANI_LABS = ('STFP', 'STFI')
HDR = re.compile(r'^(?P<pre>.+)_(?P<lab>[^_]+)_(?P<asv>ASV\d+)_(?P<coll>[^_]+)$')


def rd(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def pct(a, b, nd=1):
    return f'{100 * a / b:.{nd}f}\\%' if b else '--'


def mpct(x, nd=1):
    return f'{100 * x:.{nd}f}\\%' if x is not None else '--'


def num(n):
    return f'{n:,}'.replace(',', '{,}')


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


def known(r):
    return r['status'].startswith('known')


def macro_rate(rows, key):
    by = defaultdict(list)
    for r in rows:
        by[r[key]].append(known(r))
    return statistics.mean(statistics.mean(v) for v in by.values()) if by else None


def kept(r):
    return known(r) and r['organism'] in r['placement'].split(';')


def tex_name(org):
    """'Gigaspora gigantea' -> \\spn{Gigaspora gigantea}; 'Racocetra sp.' -> \\spn{Racocetra} sp."""
    if org.endswith(' sp.'):
        return f'\\spn{{{org[:-4]}}} sp.'
    return f'\\spn{{{org}}}'


def called_name(subject, org):
    """Abbreviate a second name in the subject's genus: called(Funneliformis coronatus, Funneliformis mosseae)
    -> \\spn{F.\\ mosseae}."""
    g = subject.split()[0].strip('[]')
    if org.split()[0].strip('[]') == g and len(org.split()) > 1 and not org.endswith(' sp.'):
        return f'\\spn{{{g[0]}.\\ {" ".join(org.split()[1:])}}}'
    return tex_name(org)


def code_of(culture):
    return culture.split(' | ', 1)[1].replace('LGDB ', '').replace('_', '\\_') if ' | ' in culture else culture


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', type=Path, required=True)
    ap.add_argument('--bench', type=Path)
    ap.add_argument('--paper', type=Path, default=Path('.'))
    a = ap.parse_args()
    m = Macros()
    sc = a.runs / '06_scores'
    T = {p.name.replace('.copies.tsv', ''): rd(p) for p in sorted(sc.glob('*.copies.tsv'))}
    need = ['A_primary_v1', 'A_primary_v2', 'B_primary_v1', 'B_primary_v2',
            'A_secondary_vt_v1', 'A_secondary_vt_v2', 'B_secondary_vt_v1', 'B_secondary_vt_v2']
    missing = [n for n in need if n not in T]
    if missing:
        raise SystemExit(f'missing score tables in {sc}: {missing}')
    meta = json.loads((sc / 'scores.json').read_text())
    plan = json.loads((a.runs / 'plan.json').read_text())
    ev = json.loads((a.runs / '01_evidence/evidence_report.json').read_text())

    # ---- test set and integrity
    A1, A2, B1, B2 = T['A_primary_v1'], T['A_primary_v2'], T['B_primary_v1'], T['B_primary_v2']
    SA1, SA2, SB2 = T['A_secondary_vt_v1'], T['A_secondary_vt_v2'], T['B_secondary_vt_v2']
    role_n = Counter(r['test_role'] for r in A2)
    novel_a = [r for r in A2 if r['test_role'] in NOVEL]
    m.add('LgTestASVs', num(meta['test_copies']), 'LongGloDB test set (lgdb_test.py, 06_scores)')
    m.add('LgTestCultures', meta['test_cultures'])
    m.add('LgKnownN', role_n[KNOWN])
    m.add('LgKnownCultures', len({r['culture'] for r in A2 if r['test_role'] == KNOWN}))
    m.add('LgNovelN', len(novel_a))
    m.add('LgNovelSpN', role_n[NOVEL_SP])
    m.add('LgNovelLinN', role_n[NOVEL_LIN])
    m.add('LgNovelCultures', len({r['culture'] for r in novel_a}))
    m.add('LgNovelSpecies', len({r['organism'] for r in novel_a}))
    m.add('LgGenusOnlyN', role_n[GENUS_ONLY])
    m.add('LgFolds', len(plan['b_folds']), 'condition B folds; species V18 does not hold by label')
    m.add('LgAbsentByLabel', len(plan['absent_by_label']))
    m.add('LgParityCopies', meta['parity_m4']['copies'], 'M4 parity before scoring')
    m.add('LgParityFolds', meta['parity_m4']['folds'])
    pa = a.runs / 'A' / 'parity.json'
    if pa.exists():
        p = json.loads(pa.read_text())
        n = next((v for k, v in p.items() if isinstance(v, int) and 'alignment' in k), None)
        if n is not None:
            m.add('LgParityA', num(n))
        else:
            m.todo('LgParityA', 'no integer alignment count found in A/parity.json; check its fields')
    else:
        m.todo('LgParityA', 'A/parity.json not in the runs folder')
    m.add('LgVTAnchor', ev['VT']['test_with_anchor'], 'VT arm: test ASVs with a hit passing the frozen coverage')
    m.add('LgVTAnchorSpanned', ev['VT_spanned']['test_with_anchor'])
    m.add('LgSHAnchor', ev['SH']['test_with_anchor'])
    m.add('LgNSixThreshold', f"{float(meta['its_threshold_pct']):.1f}\\%")

    # ---- primary outcome
    def fk_block(rows, prefix):
        nov = [r for r in rows if r['test_role'] in NOVEL]
        fk = [r for r in nov if known(r)]
        m.add(f'{prefix}FK', len(fk))
        m.add(f'{prefix}FKPct', pct(len(fk), len(nov)))
        m.add(f'{prefix}FKSpMacro', mpct(macro_rate(nov, 'organism')))
        m.add(f'{prefix}FKCuMacro', mpct(macro_rate(nov, 'culture')))
        return nov, fk
    nov_b, fk_b = fk_block(B2, 'Lg')
    m.add('LgFKNovelSp', sum(r['test_role'] == NOVEL_SP for r in fk_b))
    m.add('LgFKNovelLin', sum(r['test_role'] == NOVEL_LIN for r in fk_b))
    m.add('LgFKvOne', sum(known(r) for r in B1 if r['test_role'] in NOVEL), 'same, without N6')
    fk_species = sorted({r['organism'] for r in fk_b})
    fk_cultures = sorted({r['culture'] for r in fk_b})
    fk_calls = sorted({r['placement'] for r in fk_b})
    if len(fk_cultures) == 1 and len(fk_calls) == 1:
        m.add('LgFKCase', f'{tex_name(fk_species[0])} {code_of(fk_cultures[0])} called {called_name(fk_species[0], fk_calls[0])}')
    else:
        m.add('LgFKCase', '; '.join(f'{tex_name(r["organism"])} {code_of(r["culture"])}' for r in fk_b))
    its = sorted(float(r['its_nearest_training_pct']) for r in fk_b if r['its_nearest_training_pct'])
    m.add('LgFKITSRange', f'{its[0]:.1f}--{its[-1]:.1f}\\%' if its else '--')
    fk_ids = {r['accession'] for r in fk_b}
    clades = Counter(r['lsu_external_species_clade'] for r in A2 if r['accession'] in fk_ids)
    m.add('LgFKLSUCladeA', tex_name(clades.most_common(1)[0][0]) if clades else '--',
          'V18 species clade of the false-known ASVs with V18 intact (condition A)')

    # A and B
    a_by = {r['accession']: r for r in A2}
    same = sum(a_by[r['accession']]['status'] == r['status'] and a_by[r['accession']]['placement'] == r['placement']
               for r in nov_b)
    m.add('LgABSame', same, 'novel ASVs with the same status and placement in A and B')
    m.add('LgLSUUnmatchedB', sum(r['LSU_state'] != 'matched' for r in nov_b))
    m.add('LgFKA', sum(known(r) for r in novel_a))

    # ---- secondary VT analysis
    nov_s, fk_s = fk_block(SB2, 'LgSec')
    m.add('LgSecVTChanged', sum(r['VT_state'] != a_by[r['accession']]['VT_state'] for r in SA2),
          'secondary VT: ASVs whose VT state changes (condition A)')
    m.add('LgSecStatusChanged', sum(r['status'] != a_by[r['accession']]['status'] for r in SA2))

    # ---- genus and rank
    def genus_counts(rows, role):
        s = [r for r in rows if r['test_role'] == role]
        called = sum(r['final_genus_call'] == '1' for r in s)
        right = sum(r['final_genus_compatible'] == '1' for r in s)
        return s, called, right
    s, called, right = genus_counts(B2, NOVEL_SP)
    m.add('LgGenusRight', right, 'novel species in held genera, primary B')
    m.add('LgGenusWrong', called - right)
    m.add('LgNovelSpAsLineage', sum(r['status'] == 'novel lineage' for r in s))
    m.add('LgBoundary', sum(r['status'] == 'boundary' for r in s))
    s, called, right = genus_counts(B2, NOVEL_LIN)
    m.add('LgLinGenusGiven', called)
    m.add('LgLinAsLineage', sum(r['status'] == 'novel lineage' for r in s))
    s, called, right = genus_counts(SB2, NOVEL_SP)
    m.add('LgSecGenusRight', right, 'secondary VT, B')
    m.add('LgSecGenusWrong', called - right)
    s, called, right = genus_counts(SB2, NOVEL_LIN)
    m.add('LgSecLinGenusGiven', called)
    wrong = [r for r in s if r['final_genus_call'] == '1']
    pairs = Counter((r['organism'], r['placement']) for r in wrong)
    m.add('LgSecWrongGenusDetail', '; '.join(f'the {n} ASVs of {tex_name(o)} were given \\spn{{{g}}}' for (o, g), n in pairs.most_common()) or 'none')
    units = Counter(r['VT_units'] for r in wrong)
    m.add('LgSecWrongGenusVT', units.most_common(1)[0][0] if units else '--')
    out_s = [r for r in fk_s]

    # ---- known species under A
    kn2 = [r for r in A2 if r['test_role'] == KNOWN]
    kn1 = [r for r in A1 if r['test_role'] == KNOWN]
    m.add('LgKept', sum(map(kept, kn2)), 'known species, condition A, primary, with N6')
    m.add('LgKeptPct', pct(sum(map(kept, kn2)), len(kn2)))
    m.add('LgKeptvOne', sum(map(kept, kn1)))
    m.add('LgKeptvOnePct', pct(sum(map(kept, kn1)), len(kn1)))
    m.add('LgNSixDemoted', sum(r.get('its_rule_fired') == '1' for r in kn2))
    wrong_sp = [r for r in kn2 if known(r) and not kept(r)]
    m.add('LgWrongSp', len(wrong_sp))
    ws = Counter((r['organism'], r['culture'], r['placement']) for r in wrong_sp)
    m.add('LgWrongSpDetail', '; '.join(f'{tex_name(o)} {code_of(c)} called {called_name(o, p)}' for (o, c, p), n in ws.most_common()) or 'none')
    wits = sorted(float(r['its_nearest_training_pct']) for r in wrong_sp if r['its_nearest_training_pct'])
    m.add('LgWrongSpITSRange', f'{wits[0]:.1f}--{wits[-1]:.1f}\\%' if wits else '--')
    kn2s = [r for r in SA2 if r['test_role'] == KNOWN]
    m.add('LgSecKept', sum(map(kept, kn2s)))
    m.add('LgSecKeptvOne', sum(map(kept, [r for r in SA1 if r['test_role'] == KNOWN])))
    dem = [r for r in kn2 if r.get('its_rule_fired') == '1']
    dits = sorted(float(r['its_nearest_training_pct']) for r in dem if r['its_nearest_training_pct'])
    m.add('LgNSixDemotedITSRange', f'{dits[0]:.1f}--{dits[-1]:.1f}\\%' if dits else '--')
    m.add('LgDOne', sum(r['status'] == 'divergent class' for r in A2), 'D1 calls among test ASVs, A primary')
    mx = sorted(float(r['its_nearest_training_pct']) for r in kn2 if 'MX982A' in r['culture'] and r['its_nearest_training_pct'])
    m.add('LgLeptoMXITS', (f'{mx[0]:.1f}\\%' if mx[0] == mx[-1] else f'{mx[0]:.1f}--{mx[-1]:.1f}\\%') if mx else '--',
          'A. leptoticha MX982A: full-ITS distance to the nearest reference culture')

    # ---- bench files: composition, identity, the lock
    asvs = cul = pred = None
    if a.bench:
        f = a.bench / 'overlap.asvs.tsv'
        asvs = rd(f) if f.exists() else None
        f = a.bench / 'overlap.cultures.tsv'
        cul = rd(f) if f.exists() else None
        f = a.bench / 'predictions.predictions.tsv'
        pred = {r['asv']: r for r in rd(f)} if f.exists() else None
    if asvs:
        stf = [r for r in asvs if HDR.match(r['asv']) and HDR.match(r['asv'])['lab'] in STEFANI_LABS]
        m.add('LgTotal', num(len(asvs)), 'LongGloDB composition (longglodb_overlap.py)')
        m.add('LgStefani', num(len(stf)))
        m.add('LgStefaniIdentical', num(sum(r['nearest_its_pct'] not in ('', None) and float(r['nearest_its_pct']) == 0 for r in stf)))
        near = {r['asv']: r for r in asvs}

        def identical(r):
            n = near.get(r['accession'])
            return bool(n) and n['nearest_its_pct'] not in ('', None) and float(n['nearest_its_pct']) == 0.0 \
                and n['nearest_v4'].split(' | ')[0] == r['organism']
        ident = [r for r in kn2 if identical(r)]
        other = [r for r in kn2 if not identical(r)]
        m.add('LgIdentN', len(ident), 'known-species ASVs identical in ITS to a reference copy of their species')
        m.add('LgIdentKept', sum(map(kept, ident)))
        m.add('LgNonIdentN', len(other))
        m.add('LgNonIdentKept', sum(map(kept, other)))
        m.add('LgNonIdentKeptPct', pct(sum(map(kept, other)), len(other)))
        # G. rosea: test cultures carrying the rDNA types shared by the four group-A reference cultures
        # (Stage 3 roadmap 2.3: identical types across 194757_703A, BEG9_770A, 4511_973A, CW0011_4010A)
        group_a = ('194757_703A', 'BEG9_770A', '4511_973A', 'CW0011_4010A')
        ro = [r for r in asvs if r['asv'].startswith('Gigaspora_rosea') and HDR.match(r['asv'])
              and HDR.match(r['asv'])['lab'] not in STEFANI_LABS and r['nearest_its_pct'] not in ('', None)
              and r['nearest_v4'].split(' | ')[-1] in group_a and float(r['nearest_its_pct']) <= 0.5]
        codes = sorted({HDR.match(r['asv'])['pre'].split('_')[2] for r in ro})
        if codes:
            m.add('LgRoseaCultures', ', '.join(codes[:-1]) + ' and ' + codes[-1] if len(codes) > 1 else codes[0],
                  'G. rosea test cultures with ASVs within 0.5% (ITS) of the group-A types')
            m.add('LgRoseaMaxITS', f"{max(float(r['nearest_its_pct']) for r in ro):.2f}\\%")
        else:
            m.todo('LgRoseaCultures', 'no G. rosea test ASV near the group-A types')
            m.todo('LgRoseaMaxITS', 'no G. rosea test ASV near the group-A types')
    else:
        for n in ('LgRoseaCultures', 'LgRoseaMaxITS'):
            m.todo(n, 'needs bench/longglodb/overlap.asvs.tsv')
        for n in ('LgTotal', 'LgStefani', 'LgStefaniIdentical', 'LgIdentN', 'LgIdentKept', 'LgNonIdentN',
                  'LgNonIdentKept', 'LgNonIdentKeptPct'):
            m.todo(n, 'needs bench/longglodb/overlap.asvs.tsv')
    if cul and any(r.get('role', '').startswith('code shared') for r in cul):
        m.add('LgCodeShared', sum(r['role'].startswith('code shared') for r in cul))
        m.add('LgCodeConflicts', sum(r['role'] == 'code shared with Stefani, names differ' for r in cul))
    else:
        m.todo('LgCodeShared', 'needs the current overlap.cultures.tsv')
        m.todo('LgCodeConflicts', 'needs the current overlap.cultures.tsv')
    if pred:
        locked = {x for x, r in pred.items() if r['role'] in NOVEL and r['inside']}
        m.add('LgLockedN', len(locked), 'the lock (predictions.predictions.tsv)')
        m.add('LgFKOutsideLocked', sum(r['accession'] not in locked for r in fk_b))
        m.add('LgSecFKOutsideLocked', sum(r['accession'] not in locked for r in fk_s))
        own_in = lambda x: pred[x]['own_species_in_v4'] and pred[x]['own_species_in_v4'] in pred[x]['inside'].split('; ')
        pin = [r for r in kn2 if own_in(r['accession'])]
        pout = [r for r in kn2 if not own_in(r['accession'])]
        m.add('LgPredInsideOwn', len(pin))
        m.add('LgPredInsideOwnKept', sum(map(kept, pin)))
        m.add('LgPredOutsideOwn', len(pout))
        m.add('LgPredOutsideOwnNotKept', sum(not kept(r) for r in pout))
        m.add('LgWrongSpPredicted', sum(r['placement'] in pred[r['accession']]['inside'].split('; ') for r in wrong_sp))
    else:
        for n in ('LgLockedN', 'LgFKOutsideLocked', 'LgSecFKOutsideLocked', 'LgPredInsideOwn', 'LgPredInsideOwnKept',
                  'LgPredOutsideOwn', 'LgPredOutsideOwnNotKept', 'LgWrongSpPredicted'):
            m.todo(n, 'needs bench/longglodb/predictions.predictions.tsv')
    sec_out = [r for r in fk_s if r['accession'] not in {x['accession'] for x in fk_b}]
    m.add('LgSecExtraFKDetail', '; '.join(f'{tex_name(r["organism"])} {code_of(r["culture"])} called '
                                          f'{" / ".join(called_name(r["organism"], p) for p in r["placement"].split(";"))}' for r in sec_out) or 'none')
    sits = [float(r['its_nearest_training_pct']) for r in sec_out if r['its_nearest_training_pct']]
    m.add('LgSecExtraFKITS', f'{min(sits):.1f}\\%' if sits else '--')

    # ---- table
    def row(name, cond, rows, rows_v1, known_rows=None, known_rows_v1=None):
        nov = [r for r in rows if r['test_role'] in NOVEL]
        fk = sum(known(r) for r in nov)
        fk1 = sum(known(r) for r in rows_v1 if r['test_role'] in NOVEL)
        sp = [r for r in rows if r['test_role'] == NOVEL_SP]
        right = sum(r['final_genus_compatible'] == '1' for r in sp)
        wrongg = sum(r['final_genus_call'] == '1' and r['final_genus_compatible'] != '1' for r in sp)
        lin = [r for r in rows if r['test_role'] == NOVEL_LIN]
        lin_g = sum(r['final_genus_call'] == '1' for r in lin)
        if known_rows is not None:
            k = [r for r in known_rows if r['test_role'] == KNOWN]
            k1 = [r for r in known_rows_v1 if r['test_role'] == KNOWN]
            kc = f'{sum(map(kept, k))} ({sum(map(kept, k1))})'
        else:
            kc = '--'
        return (f'{name} & {cond} & {fk} ({fk1}) & {pct(fk, len(nov))} & {right} / {wrongg} / {len(sp) - right - wrongg} '
                f'& {lin_g} & {kc}\\\\')
    body = '\n'.join([
        row('Primary', 'A: V18 intact', A2, A1, A2, A1),
        row('', 'B: target removed', B2, B1),
        row('VT spanned', 'A: V18 intact', SA2, SA1, SA2, SA1),
        row('', 'B: target removed', SB2, T['B_secondary_vt_v1']),
    ])
    table = f'''\\begin{{table}}[htbp]\\centering
\\begin{{threeparttable}}
\\caption{{Independent cultures: the frozen rules on {meta['test_cultures']} LongGloDB cultures from other collections.}}
\\label{{tab:lgdb}}
\\footnotesize\\setlength{{\\tabcolsep}}{{4pt}}
\\begin{{tabular}}{{p{{0.10\\linewidth}}p{{0.16\\linewidth}}p{{0.10\\linewidth}}p{{0.07\\linewidth}}p{{0.17\\linewidth}}p{{0.13\\linewidth}}p{{0.11\\linewidth}}}}
\\toprule
Analysis & Condition & False known (of {len(novel_a)}) & Share & Novel species: genus right / wrong / none (of {role_n[NOVEL_SP]}) & Novel genera: genus given (of {role_n[NOVEL_LIN]}) & Known kept (of {role_n[KNOWN]})\\\\
\\midrule
{body}
\\bottomrule
\\end{{tabular}}
\\begin{{tablenotes}}\\footnotesize
\\item Protocol v2 (with the ITS rule, N6); in brackets, the same without it. \\emph{{False known}}: ASVs of species the
reference does not hold called a known species or species group. \\emph{{Novel species}}: ASVs of species the reference
does not hold, in genera it holds. \\emph{{Novel genera}}: ASVs of genera it does not hold, so any genus given is wrong.
\\emph{{Known kept}}: ASVs of species the reference holds whose call contains their species; condition A only.
\\emph{{VT spanned}}: the secondary analysis, with VT coverage measured over the part of each MaarjAM reference the
query spans (amendment 2).
\\end{{tablenotes}}
\\end{{threeparttable}}
\\end{{table}}
'''
    (a.paper / 'tables').mkdir(parents=True, exist_ok=True)
    (a.paper / 'tables/table_lgdb.tex').write_text(table)
    (a.paper / 'numbers_lgdb.tex').write_text(
        '% GENERATED by analysis/lgdb_numbers.py from the LongGloDB test outputs; do not edit.\n'
        + '\n'.join(m.lines) + '\n')
    print(f'wrote {a.paper}/numbers_lgdb.tex ({len(m.lines)} lines) and tables/table_lgdb.tex')
    if m.todos:
        print(f'{len(m.todos)} macros are \\todo placeholders: ' + ', '.join(m.todos))


if __name__ == '__main__':
    sys.exit(main())
