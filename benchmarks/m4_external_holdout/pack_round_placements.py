#!/usr/bin/env python3
"""Collect the per-round outputs of the LSU-removal test into the repository (release 0.3.2).

The full runs keep one folder per round under their --out directory (ignored by git as runs/). This
script finds each round of holdout_plan.json there, checks it, and with --copy writes a compact copy to
rounds/<key>/: the round's LSU placement table (05_summary/lsu_placements.copies.tsv, which Gate 0b of
the Rust port reads), its DONE and metrics files, and its four version 2 per-sequence tables. Alignments,
jplace files and trees stay on the PC, as in the collect step.

Checks, for every round:
  - DONE_v2.json's result hash matches 06_evaluation_v2/external_holdout_v2.copies.tsv;
  - a round with removed V18 tips has a placement table covering exactly the 903 culture copies;
  - species rounds only: every published row of aggregate_v2 (species.v1, v2, control_v1, control_v2)
    for that target equals the round's own row, so the packed placements are the ones behind the paper.

    python3 pack_round_placements.py RUNS [RUNS ...]          report only, changes nothing
    python3 pack_round_placements.py RUNS [RUNS ...] --copy   also writes rounds/ and rounds/MANIFEST.json

With several RUNS folders, the first that holds a complete round is used, and any later copy of the same
round is compared with it. --copy writes nothing unless all 52 rounds pass. Exit status 1 if a round is
missing, damaged or differs from the publication.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parent
REF = ROOT.parent.parent / 'reference/culture_ref_v4/copies.tsv'
EVAL_V2 = ('external_holdout_v1', 'external_holdout_v2', 'full_reference_v1', 'full_reference_v2')
PUBLISHED = {'external_holdout_v1': 'species.v1', 'external_holdout_v2': 'species.v2',
             'full_reference_v1': 'species.control_v1', 'full_reference_v2': 'species.control_v2'}
EXTRA = {'target', 'external_status', 'development_fold'}          # added by aggregate_v2


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def find_round(runs, key):
    """Folders holding this round with a version 2 result, in the order given."""
    return [r / key for r in runs if (r / key / 'DONE_v2.json').exists()]


def check_round(plan, fold, copies, published):
    problems = []
    done = json.loads((fold / 'DONE_v2.json').read_text())
    result = fold / '06_evaluation_v2/external_holdout_v2.copies.tsv'
    if not result.exists() or done.get('result_sha256') != sha(result):
        problems.append('DONE_v2 result hash does not match external_holdout_v2.copies.tsv')
    placement = fold / '05_summary/lsu_placements.copies.tsv'
    if plan['removed_ref_ids']:
        if not placement.exists():
            problems.append('no 05_summary/lsu_placements.copies.tsv')
        else:
            ids = {r['accession'] for r in read_tsv(placement)}
            if ids != copies:
                problems.append(f'placement table covers {len(ids)} accessions, not the {len(copies)} culture copies')
    for name in EVAL_V2:
        if not (fold / '06_evaluation_v2' / f'{name}.copies.tsv').exists():
            problems.append(f'missing 06_evaluation_v2/{name}.copies.tsv')
    if plan['rank'] == 'species' and not problems:
        for name, table in PUBLISHED.items():
            mine = {r['accession']: r for r in read_tsv(fold / '06_evaluation_v2' / f'{name}.copies.tsv')}
            theirs = [r for r in published[table] if r['target'] == plan['target']]
            if len(theirs) != len(mine):
                problems.append(f'{name}: {len(mine)} rows here, {len(theirs)} published')
                continue
            bad = 0
            for r in theirs:
                m = mine.get(r['accession'])
                if m is None or any(m.get(k, '') != v for k, v in r.items() if k not in EXTRA):
                    bad += 1
            if bad:
                problems.append(f'{name}: {bad} of {len(theirs)} rows differ from aggregate_v2/{table}.copies.tsv')
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs', nargs='+', type=Path, help='run output folders (the --out of external_holdout_v2.py)')
    ap.add_argument('--copy', action='store_true', help='write rounds/ and rounds/MANIFEST.json')
    a = ap.parse_args()
    runs = [r.expanduser().resolve() for r in a.runs]
    plans = json.loads((ROOT / 'holdout_plan.json').read_text())['folds']
    copies = {r['accession'] for r in read_tsv(REF) if r['status'] == 'retained'}
    published = {t: read_tsv(ROOT / 'aggregate_v2' / f'{t}.copies.tsv') for t in PUBLISHED.values()}
    manifest, failed, passed = {'source_runs': [str(r) for r in runs], 'rounds': {}}, 0, []
    for plan in plans:
        key = plan['key']
        found = find_round(runs, key)
        if not found:
            print(f'MISSING  {key}')
            failed += 1
            continue
        fold = found[0]
        problems = check_round(plan, fold, copies, published)
        for other in found[1:]:
            for rel in ('05_summary/lsu_placements.copies.tsv', '06_evaluation_v2/external_holdout_v2.copies.tsv'):
                p, q = fold / rel, other / rel
                if p.exists() and q.exists() and sha(p) != sha(q):
                    problems.append(f'{rel} differs in {other.parent}')
        status = 'OK' if not problems else 'PROBLEM'
        failed += bool(problems)
        own = 'own placement' if plan['removed_ref_ids'] else 'full-reference placement'
        print(f'{status:8} {key}  ({plan["rank"]}, {own}; from {fold.parent})')
        for p in problems:
            print(f'           {p}')
        if not problems:
            passed.append((plan, fold))
    for plan, fold in (passed if a.copy and not failed else []):
            key = plan['key']
            dest = ROOT / 'rounds' / key
            files = [f for f in ('DONE.json', 'DONE_v2.json', '05_summary/lsu_placements.copies.tsv',
                                 '06_evaluation/metrics.json', '06_evaluation_v2/metrics.json',
                                 *(f'06_evaluation_v2/{n}.copies.tsv' for n in EVAL_V2)) if (fold / f).exists()]
            for f in files:
                (dest / f).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(fold / f, dest / f)
            manifest['rounds'][key] = {'rank': plan['rank'], 'target': plan['target'],
                                       'own_placement': bool(plan['removed_ref_ids']),
                                       'source': str(fold), 'files': {f: sha(dest / f) for f in files}}
    species = sum(p['rank'] == 'species' and bool(p['removed_ref_ids']) for p in plans)
    genus = sum(p['rank'] == 'genus' for p in plans)
    print(f'\n{len(plans)} rounds planned ({species} species rounds with their own placement, '
          f'{genus} genus rounds); {failed} missing or with problems')
    if a.copy and failed:
        print('Nothing written: --copy needs every round to pass.')
    elif a.copy:
        (ROOT / 'rounds').mkdir(exist_ok=True)
        (ROOT / 'rounds' / 'MANIFEST.json').write_text(json.dumps(manifest, indent=1, sort_keys=True) + '\n')
        print(f'Wrote {ROOT / "rounds"}')
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
