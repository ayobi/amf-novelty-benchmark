#!/usr/bin/env python3
"""Artefact check on the divergent-class copies: long deletions per region, and recurrence across isolates.

PROTOCOL_COPY_CHECKS.md, section 2. Verdict (the rule of the comment):
                 intact        deleted
  recurring      paralog       recurrent deleted copy (a genomic copy; possibly a pseudogene)
  single         unresolved    probable artefact

Deletions. Each copy is aligned to the nearest main-type copy of its own culture (whole-copy distance, never
itself): global alignment, free end gaps, BLASTN default scores (match 2, mismatch -3, gap existence 5,
extension 2, so a gap of k nt costs 5 + 2k). Gap runs at either end of the alignment are trimming differences
and are ignored. Each internal deletion (comparator bases missing from the copy) goes to the comparator
region (SSU flank, ITS1, 5.8S, ITS2, LSU flank) holding most of its bases; insertions go to the region of the
comparator base before them. A deletion in region R is long when it is >= L_R nt, L_R = max(3, P99_R), P99_R
being the 99th percentile (exclusive method) of the longest deletion in R per main-type copy against its
nearest main-type culture-mate (0 when none), over main-type copies of cultures with two or more of them.

Recurrence. Same type: full-ITS distance <= 1.0% (edlib infix, % of the shorter sequence). A divergent copy
recurs across isolates when a copy of the same type, of any class, lies under another physical culture label.
Reported, not used in the verdict: library replicate, external record (best UNITE+INSD hit in the crosswalk at
>= 99.0% identity over >= 90% coverage), the class and species of the closest same-type copy under another
label, SH and VT paralog units (reached by at most one main-type copy), 5.8S covariance-model score.

  python3 stage3/copy_artefacts.py --ref reference/culture_ref_v4 \\
      --crosswalk paper/analysis/xw/crosswalk_v4final.copies.tsv \\
      --cmsearch benchmarks/stefani/s58.cmsearch.tbl --out benchmarks/copy_checks/artefacts

Writes {out}/artefacts.copies.tsv (all copies), artefacts.types.tsv, artefacts.summary.tsv,
artefacts.thresholds.tsv and artefacts.json.
"""
import argparse
import collections
import hashlib
import json
import os
import statistics
import sys
import time

import edlib
from Bio import Align
import Bio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arms import read_fasta, read_tsv, write_tsv  # noqa: E402

REGIONS = ['ssu_flank', 'its1', 's58', 'its2', 'lsu_flank']
CONSERVED = ['ssu_flank', 's58', 'lsu_flank']
SAME_TYPE = 1.0                 # % full ITS
SENS_TYPE = (0.5, 2.0)
FLOOR = 3                       # nt
STRICT = 10                     # nt, sensitivity
EXT_ID, EXT_COV = 99.0, 0.90    # external record
PCT = 99


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode='HW', task='distance')['editDistance'] / len(q)


def p99(vals):
    """statistics.quantiles(..., method='exclusive') at the 99th percentile."""
    if len(vals) < 2:
        return 0.0
    return statistics.quantiles(vals, n=100, method='exclusive')[PCT - 1]


def make_aligner():
    al = Align.PairwiseAligner()
    al.mode = 'global'
    al.match_score, al.mismatch_score = 2, -3
    al.open_gap_score, al.extend_gap_score = -7, -2   # BLAST existence 5 + extension 2 for the first base
    al.end_gap_score = 0                               # free end gaps (trimming differences)
    return al


def indels(aligner, query, target):
    """Internal indels of query against target: [('del'|'ins', target_pos, length)], target_pos 0-based.
    'del' = target bases [pos, pos+len) missing from the query; 'ins' = query bases inserted before target pos."""
    aln = aligner.align(target, query)[0]
    tb, qb = aln.aligned
    out = []
    for i in range(1, len(tb)):
        t_prev_end, t_start = tb[i - 1][1], tb[i][0]
        q_prev_end, q_start = qb[i - 1][1], qb[i][0]
        if t_start > t_prev_end:
            out.append(('del', int(t_prev_end), int(t_start - t_prev_end)))
        if q_start > q_prev_end:
            out.append(('ins', int(t_start), int(q_start - q_prev_end)))
    return out


def region_of(bounds, kind, pos, length):
    """bounds: [(region, start, end)] of the comparator, 0-based half-open."""
    if kind == 'del':
        best, bo = None, -1
        for r, s, e in bounds:
            o = min(e, pos + length) - max(s, pos)
            if o > bo:
                best, bo = r, o
        return best
    p = max(pos - 1, 0)
    for r, s, e in bounds:
        if s <= p < e:
            return r
    return bounds[-1][0]


def cm_scores(path):
    best = {}
    if not path:
        return best
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            c = line.split()
            acc, score = c[0], float(c[14])
            best[acc] = max(best.get(acc, -1e9), score)
    return best


def single_linkage(items, d, thr):
    parent = {x: x for x in items}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, x in enumerate(items):
        for y in items[i + 1:]:
            if d[(x, y)] <= thr:
                parent[find(x)] = find(y)
    groups = collections.defaultdict(list)
    for x in items:
        groups[find(x)].append(x)
    return sorted((sorted(v) for v in groups.values()), key=lambda v: (-len(v), v[0]))


def verdict(recurring, intact):
    if recurring:
        return 'paralog' if intact else 'recurrent deleted copy'
    return 'unresolved' if intact else 'probable artefact'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--crosswalk', required=True)
    ap.add_argument('--cmsearch', default='')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()

    rows = [r for r in read_tsv(os.path.join(a.ref, 'copies.tsv')) if r['status'] == 'retained']
    by_acc = {r['accession']: r for r in rows}
    seq = {m: read_fasta(os.path.join(a.ref, 'fasta', f'{m}.fasta')) for m in REGIONS + ['full', 'full_its']}
    for acc in by_acc:
        if seq['full'][acc] != ''.join(seq[m][acc] for m in REGIONS):
            sys.exit(f'{acc}: the regions do not tile the full copy')
    xw = {r['accession']: r for r in read_tsv(a.crosswalk)}
    missing = set(by_acc) - set(xw)
    if missing:
        sys.exit(f'{len(missing)} copies missing from the crosswalk, e.g. {sorted(missing)[:3]}')
    cms = cm_scores(a.cmsearch)
    by_cul = collections.defaultdict(list)
    for r in rows:
        by_cul[r['culture']].append(r['accession'])

    # 1. comparator: nearest main-type copy of the same culture, never the copy itself
    comp, comp_d = {}, {}
    for cul, accs in by_cul.items():
        mains = sorted(x for x in accs if by_acc[x]['copy_class'] == 'main')
        for acc in accs:
            cands = [(dist(seq['full'][acc], seq['full'][m]), m) for m in mains if m != acc]
            if cands:
                comp_d[acc], comp[acc] = min(cands)

    # 2. indels per region
    aligner = make_aligner()
    ind = {}
    for acc, c in comp.items():
        lens = [len(seq[m][c]) for m in REGIONS]
        bounds, s = [], 0
        for m, L in zip(REGIONS, lens):
            bounds.append((m, s, s + L))
            s += L
        ind[acc] = [(k, region_of(bounds, k, p, L), p, L) for k, p, L in indels(aligner, seq['full'][acc], seq['full'][c])]

    def longest(acc, kind, region):
        return max([L for k, r, p, L in ind.get(acc, []) if k == kind and r == region], default=0)

    # 3. thresholds from main-type copies with a main-type culture-mate
    base = [acc for acc in comp if by_acc[acc]['copy_class'] == 'main']
    thr = {}
    trows = []
    for m in REGIONS:
        vals = [longest(acc, 'del', m) for acc in base]
        p = p99(vals)
        thr[m] = max(float(FLOOR), p)
        trows.append({'region': m, 'main_copies': len(vals), 'with_deletion': sum(v > 0 for v in vals),
                      'max_deletion': max(vals), 'p99': f'{p:.2f}', 'L_R': f'{thr[m]:.2f}'})
    write_tsv(os.path.join(a.out, 'artefacts.thresholds.tsv'), trows,
              ['region', 'main_copies', 'with_deletion', 'max_deletion', 'p99', 'L_R'])

    def long_regions(acc, limit=None):
        out = []
        for m in REGIONS:
            L = longest(acc, 'del', m)
            if L >= (limit if limit is not None else thr[m]):
                out.append(m)
        return out

    # 4. recurrence of divergent copies (full ITS)
    divs = sorted(x for x in by_acc if by_acc[x]['copy_class'] == 'divergent')
    d_full = {}
    for x in divs:
        for y in by_acc:
            if x != y:
                d_full[(x, y)] = d_full[(y, x)] = dist(seq['full_its'][x], seq['full_its'][y])

    def recurrence(x, t):
        lab = by_acc[x]['culture_label']
        hits = [(d_full[(x, y)], y) for y in by_acc if y != x and by_acc[y]['culture_label'] != lab
                and d_full[(x, y)] <= t]
        return sorted(hits)

    # paralog units (paper: units reached by divergent copies and by at most one main-type copy)
    mainc = {arm: collections.Counter(r[f'{arm}_label'] for r in xw.values()
                                      if r['copy_class'] == 'main' and r[f'{arm}_call'] == 'assigned')
             for arm in ('SH', 'VT')}

    # types: single linkage among divergent copies
    types = single_linkage(divs, d_full, SAME_TYPE)
    type_of = {x: f'T{i + 1:02d}' for i, g in enumerate(types) for x in g}

    # culture main-type 5.8S CM score medians
    cm_med = {}
    for cul, accs in by_cul.items():
        v = [cms[x] for x in accs if by_acc[x]['copy_class'] == 'main' and x in cms]
        cm_med[cul] = statistics.median(v) if v else None

    out = []
    for acc in sorted(by_acc):
        r = by_acc[acc]
        o = {k: r[k] for k in ('accession', 'culture', 'culture_label', 'library', 'organism', 'genus', 'order',
                               'copy_class', 'div_kind', 's58_len', 's58_delta')}
        o['comparator'] = comp.get(acc, '')
        o['comparator_dist_pct'] = f'{comp_d[acc]:.2f}' if acc in comp_d else ''
        for m in REGIONS:
            o[f'del_{m}'] = longest(acc, 'del', m) if acc in comp else ''
            o[f'ins_{m}'] = longest(acc, 'ins', m) if acc in comp else ''
        o['indels'] = ';'.join(f'{k}:{rg}:{p + 1}:{L}' for k, rg, p, L in ind.get(acc, []))
        lr = long_regions(acc) if acc in comp else []
        o['long_deletion_regions'] = ','.join(lr)
        o['intact'] = ('yes' if not lr else 'no') if acc in comp else ''
        o['long_deletion_regions_10nt'] = ','.join(long_regions(acc, STRICT)) if acc in comp else ''
        sc = cms.get(acc)
        o['s58_cm_score'] = f'{sc:.1f}' if sc is not None else ''
        med = cm_med.get(r['culture'])
        o['s58_cm_delta'] = f'{sc - med:.1f}' if sc is not None and med is not None else ''
        x = xw[acc]
        for arm in ('SH', 'VT'):
            lab = x[f'{arm}_label'] if x[f'{arm}_call'] == 'assigned' else ''
            o[f'{arm}_call'] = x[f'{arm}_call']
            o[f'{arm}_label'] = lab
            o[f'{arm}_main_copies'] = mainc[arm][lab] if lab else ''
            o[f'{arm}_paralog_unit'] = ('yes' if mainc[arm][lab] <= 1 else '') if lab and r['copy_class'] == 'divergent' else ''
        if r['copy_class'] == 'divergent':
            hits = recurrence(acc, SAME_TYPE)
            labs = collections.OrderedDict()
            for dd, y in hits:
                labs.setdefault(by_acc[y]['culture_label'], (dd, y))
            o['recur_labels_n'] = len(labs)
            o['recur_labels'] = ','.join(f"{lb}:{by_acc[y]['organism']}:{by_acc[y]['copy_class']}:{dd:.2f}"
                                         for lb, (dd, y) in labs.items())
            for t in SENS_TYPE:
                o[f'recur_labels_n_{t:g}pct'] = len({by_acc[y]['culture_label'] for _, y in recurrence(acc, t)})
            others = sorted((d_full[(acc, y)], y) for y in by_acc
                            if y != acc and by_acc[y]['culture_label'] != r['culture_label'])
            dd, y = others[0]
            o['nearest_other_label'] = f"{y}|{by_acc[y]['culture']}"
            o['nearest_other_label_dist_pct'] = f'{dd:.2f}'
            o['nearest_other_label_class'] = by_acc[y]['copy_class']
            o['nearest_same_type_main_other_species'] = (
                'yes' if hits and by_acc[hits[0][1]]['copy_class'] == 'main'
                and by_acc[hits[0][1]]['organism'] != r['organism'] else '')
            libs = {by_acc[z]['library'] for z in by_cul[r['culture']]}
            if len(libs) < 2:
                o['library_replicate'] = 'single library'
            else:
                rep = sorted({by_acc[z]['library'] for z in by_cul[r['culture']] if z != acc
                              and by_acc[z]['library'] != r['library'] and d_full[(acc, z)] <= SAME_TYPE})
                o['library_replicate'] = 'yes' if rep else 'no'
                o['library_replicate_libs'] = ','.join(rep)
            bp = x['SH_best_pident']
            cov = x['SH_cov']
            o['external_best_pident'] = bp
            o['external_best_cov'] = cov
            o['external_best_taxon'] = x['SH_best_taxon']
            ext = bp != '' and cov != '' and float(bp) >= EXT_ID and float(cov) >= EXT_COV
            o['external_record'] = 'yes' if ext else ''
            intact = not lr
            o['type'] = type_of[acc]
            o['verdict'] = verdict(o['recur_labels_n'] > 0, intact)
            o['reproduced'] = ('yes' if o['recur_labels_n'] == 0 and (o['library_replicate'] == 'yes' or ext) else '')
            for t in SENS_TYPE:
                o[f'verdict_{t:g}pct'] = verdict(o[f'recur_labels_n_{t:g}pct'] > 0, intact)
            o['verdict_del10nt'] = verdict(o['recur_labels_n'] > 0, not long_regions(acc, STRICT))
        out.append(o)

    cols = (['accession', 'culture', 'culture_label', 'library', 'organism', 'genus', 'order', 'copy_class',
             'div_kind', 's58_len', 's58_delta', 'type', 'verdict', 'reproduced', 'intact',
             'long_deletion_regions', 'long_deletion_regions_10nt', 'recur_labels_n', 'recur_labels',
             'recur_labels_n_0.5pct', 'recur_labels_n_2pct', 'nearest_other_label', 'nearest_other_label_dist_pct',
             'nearest_other_label_class', 'nearest_same_type_main_other_species', 'library_replicate',
             'library_replicate_libs', 'external_record', 'external_best_pident', 'external_best_cov',
             'external_best_taxon', 'SH_call', 'SH_label', 'SH_main_copies', 'SH_paralog_unit', 'VT_call', 'VT_label',
             'VT_main_copies', 'VT_paralog_unit', 's58_cm_score', 's58_cm_delta', 'comparator', 'comparator_dist_pct']
            + [f'del_{m}' for m in REGIONS] + [f'ins_{m}' for m in REGIONS]
            + ['indels', 'verdict_0.5pct', 'verdict_2pct', 'verdict_del10nt'])
    write_tsv(os.path.join(a.out, 'artefacts.copies.tsv'), out, cols)
    D = {o['accession']: o for o in out}

    # types
    trow = []
    for i, g in enumerate(types):
        xs = [D[x] for x in g]
        trow.append({
            'type': f'T{i + 1:02d}', 'copies': len(g),
            'classes': ','.join(f'{k}:{v}' for k, v in collections.Counter(o['div_kind'] for o in xs).items()),
            'cultures': len({o['culture'] for o in xs}), 'labels': len({o['culture_label'] for o in xs}),
            'species': ','.join(sorted({o['organism'] for o in xs})),
            'libraries': len({o['library'] for o in xs}),
            's58_len': ','.join(sorted({o['s58_len'] for o in xs})),
            'long_deletion_regions': ','.join(sorted({rg for o in xs for rg in o['long_deletion_regions'].split(',') if rg})),
            'SH_units': ','.join(sorted({o['SH_label'] for o in xs if o['SH_label']})),
            'VT_units': ','.join(sorted({o['VT_label'] for o in xs if o['VT_label']})),
            'other_label_copies_same_type': max(o['recur_labels_n'] for o in xs),
            'verdicts': ','.join(f'{k}:{v}' for k, v in collections.Counter(o['verdict'] for o in xs).items()),
            'accessions': ','.join(g)})
    write_tsv(os.path.join(a.out, 'artefacts.types.tsv'), trow, list(trow[0]))

    # summary
    srow = []
    dv = [D[x] for x in divs]
    for col, label in (('verdict', f'primary (same type <= {SAME_TYPE:g}%, long deletion >= L_R)'),
                       ('verdict_0.5pct', 'same type <= 0.5%'), ('verdict_2pct', 'same type <= 2.0%'),
                       ('verdict_del10nt', 'long deletion >= 10 nt')):
        for cls in ('short_5.8S', 'same_length', 'all'):
            xs = [o for o in dv if cls == 'all' or o['div_kind'] == cls]
            c = collections.Counter(o[col] for o in xs)
            srow.append({'reading': label, 'class': cls, 'copies': len(xs),
                         'paralog': c['paralog'], 'recurrent deleted copy': c['recurrent deleted copy'],
                         'unresolved': c['unresolved'], 'probable artefact': c['probable artefact'],
                         'single but reproduced': sum(o['reproduced'] == 'yes' for o in xs) if col == 'verdict' else ''})
    write_tsv(os.path.join(a.out, 'artefacts.summary.tsv'), srow, list(srow[0]))
    flagged_main = [o for o in out if o['copy_class'] == 'main'
                    and any(rg in CONSERVED for rg in o['long_deletion_regions'].split(','))]
    meta = {'finished': time.strftime('%Y-%m-%d %H:%M:%S'), 'seconds': round(time.time() - t0, 1),
            'biopython': Bio.__version__, 'edlib_mode': 'HW (infix), % of shorter',
            'aligner': {'mode': 'global', 'match': 2, 'mismatch': -3, 'gap_first_base': -7, 'gap_extend': -2,
                        'end_gaps': 0},
            'same_type_pct': SAME_TYPE, 'sensitivity_same_type_pct': SENS_TYPE, 'floor_nt': FLOOR,
            'strict_nt': STRICT, 'external': {'identity': EXT_ID, 'coverage': EXT_COV},
            'thresholds': {m: thr[m] for m in REGIONS},
            'inputs': {'copies.tsv': sha256(os.path.join(a.ref, 'copies.tsv')),
                       'crosswalk': sha256(a.crosswalk), 'cmsearch': sha256(a.cmsearch) if a.cmsearch else '',
                       **{f'fasta/{m}.fasta': sha256(os.path.join(a.ref, 'fasta', f'{m}.fasta'))
                          for m in REGIONS + ['full', 'full_its']}},
            'main_copies_with_long_deletion_in_conserved_region': [o['accession'] for o in flagged_main],
            'outputs': {f: sha256(os.path.join(a.out, f)) for f in
                        ('artefacts.copies.tsv', 'artefacts.types.tsv', 'artefacts.summary.tsv',
                         'artefacts.thresholds.tsv')}}
    with open(os.path.join(a.out, 'artefacts.json'), 'w') as fh:
        json.dump(meta, fh, indent=1)

    print('long-deletion thresholds L_R (nt): ' + ', '.join(f'{m} {thr[m]:.2f}' for m in REGIONS))
    for s in srow:
        print('\t'.join(str(s[k]) for k in s))
    print(f'divergent types: {len(types)}; main-type copies with a long deletion in SSU, 5.8S or LSU: '
          f'{len(flagged_main)}')


if __name__ == '__main__':
    main()
