#!/usr/bin/env python3
"""M4 external-LSU holdouts, protocol v2: version 1 unchanged plus the culture-ITS rule (N6).

Stages 01-05 and the version 1 evaluation (06_evaluation) run or are reused exactly as in
external_holdout.py, with the same stage keys, so completed folds are not recomputed. Version 2 adds a
separate stage, 06_evaluation_v2, that applies N6 to the version 1 calls of both the full-V18 control and
the external holdout, and writes <out>/aggregate_v2. Version 1 outputs are never modified.

Commands mirror external_holdout.py:
  python3 external_holdout_v2.py run --rank species --out runs --threads 8 [--jobs 2]
  python3 external_holdout_v2.py report --rank species --out runs
  python3 external_holdout_v2.py collect --out runs --zip ~/Downloads/amf_m4_v2_species.zip
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import external_holdout as v1          # sets up frozen_m2b imports
from its_rule import ItsRule
from lsu_placement import read_tsv, write_tsv, sha

DEVELOPMENT_FOLDS = {'Glomus chinense', 'Glomus rugosae'}   # PROTOCOL_v2.md, item 3
QUANTILE = 99
_CACHE = {}


def verify_v2():
    manifest = json.loads((ROOT / 'V2_MANIFEST.json').read_text())
    for name, expected in manifest['files'].items():
        if sha(ROOT / name) != expected:
            raise ValueError(f'Version 2 file changed: {name}')
    return v1.digest(manifest)


def rule_for(data):
    if 'rule' not in _CACHE:
        _CACHE['rule'] = ItsRule(data, QUANTILE)
    return _CACHE['rule']


def genus_fields(data, row):
    tokens = {s.split()[0].strip('[]') for s in row['placement'].split(';') if s}
    tokens = {'Rhizophagus' if g == 'Rhizoglomus' else 'Entrophospora' if g == 'Claroideoglomus' else g for g in tokens}
    row['final_genus_call'] = int(bool(tokens))
    row['final_genus_compatible'] = int(data.genus[row['culture']] in tokens)


def score_v2(plan, placement_path, path):
    """Version 1 calls for control and holdout, then N6 on copies of both. Writes into path."""
    data = v1.Data(v1.FROZEN / 'inputs/benchmark')
    ids, drop = plan['query_ids'], frozenset(plan['held_cultures'])
    if not set(ids) <= set(data.cop) or any(data.cop[x]['culture'] not in drop for x in ids):
        raise ValueError('Invalid fold queries or culture exclusions')
    old = {r['accession']: r for r in read_tsv(v1.FROZEN / 'placement/summary/lsu_placements.copies.tsv')}
    new = {r['accession']: r for r in read_tsv(placement_path)}
    rule, train = rule_for(data), data.cultures - drop
    threshold, n = rule.threshold(train)
    metrics, changes = {}, []
    for name, rows in (('full_reference', v1.infer(data, ids, drop, old)),
                       ('external_holdout', v1.infer(data, ids, drop, new))):
        before = [dict(r) for r in rows]
        after = [dict(r) for r in rows]
        rule.apply(after, {r['accession']: r for r in before}, train)
        for r in after:
            genus_fields(data, r)
        write_tsv(path / f'{name}_v1.copies.tsv', before)
        write_tsv(path / f'{name}_v2.copies.tsv', after)
        metrics[f'{name}_v1'] = v1.fold_metrics(before)
        metrics[f'{name}_v2'] = v1.fold_metrics(after)
        for b, r in zip(before, after):
            if b['status'] != r['status']:
                changes.append({'comparison': name, 'accession': r['accession'], 'culture': r['culture'],
                                'copy_class': r['copy_class'], 'v1_status': b['status'], 'v2_status': r['status'],
                                'v1_placement': b['placement'], 'v2_placement': r['placement'],
                                'its_nearest_training_pct': r['its_nearest_training_pct'], 'its_threshold_pct': r['its_threshold_pct']})
    write_tsv(path / 'v1_to_v2_changes.tsv', changes)
    v1.write_json(path / 'metrics.json', {'plan': plan, **metrics, 'its_threshold_pct': threshold, 'its_threshold_n': n,
                                          'its_quantile': QUANTILE, 'development_fold': plan['target'] in DEVELOPMENT_FOLDS,
                                          'independent_validation': False})


def run_fold_v2(plan, output, tools, threads, base_key, v2_key):
    v1.run_fold(plan, output, tools, threads, base_key)
    fold = Path(output) / plan['key']
    signature = v1.digest({'base': base_key, 'plan': plan})
    placement = (fold / '05_summary/lsu_placements.copies.tsv' if plan['removed_ref_ids']
                 else v1.FROZEN / 'placement/summary/lsu_placements.copies.tsv')

    def score(path):
        score_v2(plan, placement, path)
        for name in ('full_reference', 'external_holdout'):     # version 1 calls reproduced exactly
            mine = {r['accession']: r for r in read_tsv(path / f'{name}_v1.copies.tsv')}
            for r in read_tsv(fold / '06_evaluation' / f'{name}.copies.tsv'):
                for k in ('status', 'placement', 'reasons'):
                    if mine[r['accession']][k] != r[k]:
                        raise AssertionError(f'Version 1 parity failed: {name} {r["accession"]} {k}')
    v1.stage(fold / '06_evaluation_v2', signature + ':score_v2:' + v2_key, score)
    v1.write_json(fold / 'DONE_v2.json', {'signature': signature, 'v2_key': v2_key, 'complete': True,
                                          'result_sha256': sha(fold / '06_evaluation_v2/external_holdout_v2.copies.tsv')})
    return plan['key']


def aggregate_v2(output, rank, v2_key=None):
    plans = v1.load_plan(rank)
    rows = {k: [] for k in ('v1', 'v2', 'control_v1', 'control_v2')}
    done = []
    for plan in plans:
        fold = Path(output) / plan['key']
        if not (fold / 'DONE_v2.json').exists():
            continue
        state = json.loads((fold / 'DONE_v2.json').read_text())
        if v2_key and state['v2_key'] != v2_key:
            raise ValueError(f'{fold} was scored under different version 2 files')
        if state['result_sha256'] != sha(fold / '06_evaluation_v2/external_holdout_v2.copies.tsv'):
            raise ValueError(f'Completed version 2 result changed: {fold}')
        done.append(plan)
        for key, f in (('v1', 'external_holdout_v1'), ('v2', 'external_holdout_v2'),
                       ('control_v1', 'full_reference_v1'), ('control_v2', 'full_reference_v2')):
            for r in read_tsv(fold / '06_evaluation_v2' / f'{f}.copies.tsv'):
                r['target'], r['external_status'] = plan['target'], plan['external_status']
                r['development_fold'] = int(plan['target'] in DEVELOPMENT_FOLDS)
                rows[key].append(r)
    out = Path(output) / 'aggregate_v2'
    out.mkdir(parents=True, exist_ok=True)

    def subset(key, dev=None, status=None):
        return [r for r in rows[key] if (dev is None or r['development_fold'] == dev)
                and (status is None or r['external_status'] == status)]
    strata = {'all': {}, 'first_outcomes': {'dev': 0}, 'development': {'dev': 1}}
    statuses = sorted({p['external_status'] for p in done})
    result = {'rank': rank, 'complete': len(done) == len(plans), 'expected_folds': len(plans),
              'completed_folds': [p['target'] for p in done],
              'missing_folds': [p['target'] for p in plans if p not in done], 'independent_validation': False}
    for label, kw in strata.items():
        for st in [None] + statuses:
            name = label + ('' if st is None else f'/{st}')
            result[name] = {k: v1.fold_metrics(subset(k, kw.get('dev'), st)) for k in rows if subset(k, kw.get('dev'), st)}
    for k, rs in rows.items():
        write_tsv(out / f'{rank}.{k}.copies.tsv', rs) if rs else None
    v1.write_json(out / f'{rank}.summary.json', result)

    print(f"{rank}: {len(done)}/{len(plans)} folds scored under version 2")
    print(f"  {'stratum':34s}{'copies':>7s}{'false known v1':>16s}{'false known v2':>16s}{'full-V18 control v1 -> v2':>28s}")
    for label, kw in strata.items():
        for st in [None] + statuses:
            e1, e2 = subset('v1', kw.get('dev'), st), subset('v2', kw.get('dev'), st)
            if not e1:
                continue
            c1, c2 = subset('control_v1', kw.get('dev'), st), subset('control_v2', kw.get('dev'), st)
            fk = lambda rs: sum(r['copy_class'] == 'main' and r['group'] == 'named, species held out'
                                and r['status'].startswith('known') for r in rs)
            held = sum(r['copy_class'] == 'main' and r['group'] == 'named, species held out' for r in e1)
            name = label + ('' if st is None else f' / {st}')
            print(f"  {name:34s}{held:>7d}{fk(e1):>10d} ({100 * fk(e1) / max(held, 1):3.0f}%){fk(e2):>9d} ({100 * fk(e2) / max(held, 1):3.0f}%)"
                  f"{fk(c1):>14d} -> {fk(c2)}")
    left = Counter((r['organism'], r['placement']) for r in subset('v2', 0)
                   if r['copy_class'] == 'main' and r['group'] == 'named, species held out' and r['status'].startswith('known'))
    if left:
        print('  first-outcome false known under v2: ' + ', '.join(f'{o} -> {p} ({n})' for (o, p), n in left.most_common(12)))
    return result


def collect_v2(output, destination):
    output, destination = Path(output).resolve(), Path(destination).resolve()
    names = ('config.json', 'DONE.json', 'DONE_v2.json', 'stage.json', 'metrics.json', 'commands.json', 'validation.json',
             'full_reference.copies.tsv', 'external_holdout.copies.tsv', 'paired_changes.tsv', 'lsu_placements.copies.tsv',
             'placement_summary.json', 'fit.raxml.bestModel', 'fit.raxml.bestTree')
    files = sorted(p for p in output.rglob('*') if p.is_file() and p.resolve() != destination and (
        p.name in names or p.suffix == '.log' or '06_evaluation_v2' in p.parts or 'aggregate_v2' in p.parts
        or 'aggregate' in p.parts))
    with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in files:
            z.write(p, 'results/' + str(p.relative_to(output)))
        for name in ('PROTOCOL.md', 'FROZEN_M2B.json', 'holdout_plan.json', 'species_exclusions.tsv', 'external_holdout.py',
                     'PROTOCOL_v2.md', 'V2_MANIFEST.json', 'its_rule.py', 'external_holdout_v2.py'):
            z.write(ROOT / name, 'protocol/' + name)
    print(f'Created {destination} ({destination.stat().st_size:,} bytes). No upload performed.')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='command', required=True)
    run = sub.add_parser('run')
    run.add_argument('--rank', choices=['species', 'genus'], default='species')
    run.add_argument('--target', action='append')
    run.add_argument('--out', type=Path, default=ROOT / 'runs')
    run.add_argument('--threads', type=int, default=8)
    run.add_argument('--jobs', type=int, default=1)
    run.add_argument('--mafft', default='mafft')
    run.add_argument('--raxml-ng', default='raxml-ng')
    run.add_argument('--epa-ng', default='epa-ng')
    rep = sub.add_parser('report')
    rep.add_argument('--out', type=Path, default=ROOT / 'runs')
    rep.add_argument('--rank', choices=['species', 'genus'], default='species')
    col = sub.add_parser('collect')
    col.add_argument('--out', type=Path, default=ROOT / 'runs')
    col.add_argument('--zip', type=Path, required=True)
    a = ap.parse_args()
    if a.command == 'report':
        aggregate_v2(a.out, a.rank, verify_v2())
        return
    if a.command == 'collect':
        collect_v2(a.out, a.zip)
        return
    if a.threads < 1 or a.jobs < 1:
        ap.error('threads and jobs must be positive')
    frozen_hash = v1.verify_frozen()
    v2_key = verify_v2()
    tools = v1.tool_info(a)
    plans = v1.load_plan(a.rank, a.target)
    # identical to external_holdout.py, so completed version 1 stages are reused
    base_key = v1.digest({'frozen': frozen_hash, 'protocol': sha(ROOT / 'PROTOCOL.md'), 'plan': sha(ROOT / 'holdout_plan.json'),
                          'algorithm': 'fresh-reference-fixed-topology-v1', 'versions': {k: v['version'] for k, v in tools.items()}})
    a.out = a.out.resolve()
    a.out.mkdir(parents=True, exist_ok=True)
    v1.write_json(a.out / 'invocation_v2.json', {'argv': sys.argv, 'base_key': base_key, 'v2_key': v2_key,
                                                  'script_sha256': sha(__file__), 'tools': tools, 'threads_per_job': a.threads,
                                                  'jobs': a.jobs, 'python': sys.version, 'started': time.strftime('%Y-%m-%d %H:%M:%S')})
    if a.jobs == 1:
        for p in plans:
            run_fold_v2(p, a.out, tools, a.threads, base_key, v2_key)
    else:
        with ProcessPoolExecutor(max_workers=a.jobs) as pool:
            for f in as_completed([pool.submit(run_fold_v2, p, a.out, tools, a.threads, base_key, v2_key) for p in plans]):
                f.result()
    v1.aggregate(a.out, a.rank)
    aggregate_v2(a.out, a.rank, v2_key)


if __name__ == '__main__':
    main()
