#!/usr/bin/env python3
"""Resumable, frozen-rule external-LSU species/genus holdouts.

Run from any directory. Default output: this package's runs/ directory.
No classifier thresholds or labels are fitted by this program.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import time
import zipfile
from Bio import Phylo

ROOT = Path(__file__).resolve().parent
FROZEN = ROOT / 'frozen_m2b'
sys.path.insert(0, str(FROZEN))
from lsu_placement import read_fasta, write_fasta, read_tsv, write_tsv, sha, summarize
from compare_holdouts import Data, ArmTrainer, placement_fit, reconcile, evaluate, metrics


def write_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def verify_frozen():
    manifest = json.loads((ROOT / 'FROZEN_M2B.json').read_text())
    for name, expected in manifest['files'].items():
        if sha(FROZEN / name) != expected:
            raise ValueError(f'Frozen input changed: {name}')
    return digest(manifest)


def load_plan(rank, targets=None):
    plan = json.loads((ROOT / 'holdout_plan.json').read_text())['folds']
    selected = [p for p in plan if p['rank'] == rank and (not targets or p['target'] in targets)]
    if targets and set(targets) != {p['target'] for p in selected}:
        raise ValueError('Unknown target; use list to see the exact frozen names')
    return selected


def tool_info(args):
    result = {}
    for name, executable, flag, expected in [
        ('mafft', args.mafft, '--version', '7.526'),
        ('raxml_ng', args.raxml_ng, '--version', '1.2.2'),
        ('epa_ng', args.epa_ng, '--version', '0.3.8'),
    ]:
        found = shutil.which(executable)
        if not found:
            raise ValueError(f'Missing executable: {executable}')
        p = subprocess.run([found, flag], capture_output=True, text=True, timeout=30)
        output = (p.stdout + p.stderr).strip()
        if p.returncode or not re.search(r'(?<![\d.])' + re.escape(expected) + r'(?![\d.])', output):
            raise ValueError(f'{name}: expected {expected}; observed {output[:300]}')
        result[name] = {'executable': str(Path(found).resolve()), 'version': expected,
                        'version_output': output, 'executable_sha256': sha(found)}
    return result


def valid_stage(path, key):
    marker = path / 'stage.json'
    if not marker.exists():
        return False
    state = json.loads(marker.read_text())
    if state['key'] != key:
        raise ValueError(f'Configuration changed for completed stage {path}; choose a new output')
    for rel, expected in state['outputs'].items():
        if not (path / rel).is_file() or sha(path / rel) != expected:
            raise ValueError(f'Completed stage output changed: {path / rel}')
    return True


def stage(path, key, function):
    """Discard only an incomplete generated stage; verify completed stages."""
    if valid_stage(path, key):
        return
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    start = time.monotonic()
    function(path)
    files = {str(p.relative_to(path)): sha(p) for p in sorted(path.rglob('*')) if p.is_file()}
    write_json(path / 'stage.json', {'key': key, 'seconds': time.monotonic() - start, 'outputs': files})


def call(command, log, stdout=None):
    print(f'  {log.parent.name}: {Path(command[0]).name}', flush=True)
    with log.open('w') as err:
        if stdout:
            with stdout.open('w') as out:
                subprocess.run(command, stdout=out, stderr=err, check=True)
        else:
            subprocess.run(command, stdout=err, stderr=subprocess.STDOUT, check=True)


def pruned_tree(source, removed, expected):
    tree = Phylo.read(source, 'newick')
    tips = [n.name for n in tree.get_terminals()]
    if len(tips) != len(set(tips)) or not set(removed) <= set(tips):
        raise ValueError('Invalid pruning IDs or duplicate tips')
    for tip in sorted(removed):
        tree.prune(tip)
    if {n.name for n in tree.get_terminals()} != set(expected):
        raise ValueError('Pruning produced the wrong taxon set')
    return tree


def splits(path):
    t = Phylo.read(path, 'newick'); all_tips = {x.name for x in t.get_terminals()}
    out = set()
    for node in t.find_clades():
        side = {x.name for x in node.get_terminals()}; other = all_tips - side
        if len(side) > 1 and len(other) > 1:
            out.add(min(tuple(sorted(side)), tuple(sorted(other))))
    return out


def infer(data, ids, drop, placements):
    if set(placements) != set(data.cop):
        raise ValueError('Every training and query copy needs a fold-specific placement')
    train = data.cultures - drop
    fits = {n: ArmTrainer(data, n).fit(train) for n in ('VT', 'SH')}
    sets = {x: frozenset(filter(None, r['credible_species_clade_units'].split(',')))
            if r['bridge_alignment_eligible'] == '1' else frozenset()
            for x, r in placements.items()}
    if any(float(r['credible_mass_target']) != .95 for r in placements.values()):
        raise ValueError('Compatible mass differs from frozen 0.95')
    fits['LSU'] = placement_fit(data, train, sets)
    if any(set(v) & drop for v in fits['LSU'].holders.values()):
        raise AssertionError('Held-out culture leaked into LSU holders')
    rec = {x: reconcile(data, x, fits, drop) for x in ids}
    rows = []
    for n in evaluate(data, ids, fits, drop, rec):
        x = n['accession']; p = placements[x]
        if set(filter(None, n['d1_comparator_cultures'].split(';'))) & drop:
            raise AssertionError('D1 held-out comparator')
        if set(filter(None, n['nearest_ssu_cultures'].split(';'))) & drop:
            raise AssertionError('SSU held-out comparator')
        row = {**rec[x], **n}
        row.update({f'lsu_{k}': p[k] for k in ('best_edge_lwr', 'credible_edge_count',
                    'best_pendant_length', 'retained_fraction', 'external_species_clade',
                    'external_species_clade_mass', 'external_genus_clade', 'external_genus_clade_mass')})
        row['held_cultures'] = ';'.join(sorted(drop))
        row['training_genus_present'] = int(any(data.genus[c] == data.genus[row['culture']] for c in train))
        genus_tokens = {s.split()[0].strip('[]') for s in row['placement'].split(';') if s}
        genus_tokens = {'Rhizophagus' if g == 'Rhizoglomus' else 'Entrophospora' if g == 'Claroideoglomus' else g for g in genus_tokens}
        row['final_genus_call'] = int(bool(genus_tokens))
        row['final_genus_compatible'] = int(data.genus[row['culture']] in genus_tokens)
        rows.append(row)
    return sorted(rows, key=lambda r: r['accession'])


def fold_metrics(rows):
    m = metrics(rows)
    main = [r for r in rows if r['copy_class'] == 'main']
    for label, subset in [('all_main', main), ('training_genus_present', [r for r in main if int(r['training_genus_present'])]),
                          ('training_genus_absent', [r for r in main if not int(r['training_genus_present'])])]:
        m[label] = {'n': len(subset), 'genus_called': sum(int(r['final_genus_call']) for r in subset),
                    'genus_compatible': sum(int(r['final_genus_compatible']) for r in subset),
                    'genus_wrong': sum(int(r['final_genus_call']) and not int(r['final_genus_compatible']) for r in subset)}
    return m


def run_fold(plan, output, tools, threads, base_key):
    fold = Path(output) / plan['key']; fold.mkdir(parents=True, exist_ok=True)
    signature = digest({'base': base_key, 'plan': plan})
    config = fold / 'config.json'
    if config.exists() and json.loads(config.read_text())['signature'] != signature:
        raise ValueError(f'Changed configuration in {fold}; choose a new output')
    write_json(config, {'signature': signature, 'plan': plan, 'tools': tools, 'threads_this_invocation': threads})
    start = time.monotonic()
    print(f"Starting {plan['target']} ({len(plan['removed_ref_ids'])} external tips removed)", flush=True)
    removed = set(plan['removed_ref_ids'])
    if removed:
        def reference(path):
            meta = read_tsv(FROZEN / 'reference/reference_metadata.tsv')
            raw = read_fasta(FROZEN / 'reference/reference.unaligned.fasta')
            if not removed <= set(raw) or any(r['is_amf'] != '1' for r in meta if r['ref_id'] in removed):
                raise ValueError('Removal includes missing taxa or an outgroup')
            kept = {x: s for x, s in raw.items() if x not in removed}
            write_fasta(path / 'raw.fasta', kept)
            write_tsv(path / 'metadata.tsv', [r for r in meta if r['ref_id'] not in removed])
            tree = pruned_tree(FROZEN / 'reference/backbone.topology.newick', removed, kept)
            Phylo.write(tree, path / 'topology.newick', 'newick', format_branch_length='%1.15f')
            command = [tools['mafft']['executable'], '--localpair', '--maxiterate', '1000', '--thread', str(threads), '--threadit', '0', str(path / 'raw.fasta')]
            call(command, path / 'mafft.log', path / 'aligned.fasta')
            aligned = read_fasta(path / 'aligned.fasta')
            if set(aligned) != set(kept) or len({len(s) for s in aligned.values()}) != 1:
                raise ValueError('Invalid reference alignment')
            if any(aligned[x].replace('-', '').replace('.', '') != s.replace('-', '').replace('.', '') for x, s in kept.items()):
                raise ValueError('Reference residues changed')
            write_json(path / 'commands.json', [command])
        ref = fold / '01_reference'; stage(ref, signature + ':reference', reference)

        def model(path):
            command = [tools['raxml_ng']['executable'], '--evaluate', '--msa', str(ref / 'aligned.fasta'), '--tree', str(ref / 'topology.newick'), '--model', 'GTR+G4', '--threads', str(min(threads, 4)), '--seed', '20260929', '--prefix', str(path / 'fit'), '--force', 'perf_threads']
            call(command, path / 'raxml.log')
            if splits(ref / 'topology.newick') != splits(path / 'fit.raxml.bestTree'):
                raise ValueError('Fixed topology changed')
            write_json(path / 'commands.json', [command])
        model_dir = fold / '02_model'; stage(model_dir, signature + ':model', model)

        def align(path):
            shutil.copyfile(FROZEN / 'placement/queries.fasta', path / 'queries.fasta')
            command = [tools['mafft']['executable'], '--auto', '--addfragments', str(path / 'queries.fasta'), '--keeplength', '--mapout', '--thread', str(min(threads, 4)), str(ref / 'aligned.fasta')]
            call(command, path / 'mafft.log', path / 'combined.fasta')
            raw = read_fasta(path / 'queries.fasta'); base = read_fasta(ref / 'aligned.fasta'); full = read_fasta(path / 'combined.fasta')
            if set(raw) & set(base) or set(full) != set(base) | set(raw) or any(full[x] != s for x, s in base.items()):
                raise ValueError('Query addition altered the fixed reference or query IDs')
            write_fasta(path / 'queries.aligned.fasta', {x: full[x] for x in raw})
            write_json(path / 'commands.json', [command])
            write_json(path / 'validation.json', {'reference_columns_unchanged': True, 'queries': len(raw), 'columns': len(next(iter(base.values())))})
        alignment = fold / '03_queries'; stage(alignment, signature + ':query_alignment', align)

        def place(path):
            command = [tools['epa_ng']['executable'], '--ref-msa', str(ref / 'aligned.fasta'), '--tree', str(model_dir / 'fit.raxml.bestTree'), '--query', str(alignment / 'queries.aligned.fasta'), '--model', str(model_dir / 'fit.raxml.bestModel'), '--outdir', str(path), '--threads', str(threads), '--no-heur', '--no-pre-mask', '--filter-min-lwr', '0', '--filter-max', '2000', '--precision', '12']
            call(command, path / 'epa.log')
            write_json(path / 'commands.json', [command])
        epa = fold / '04_epa'; stage(epa, signature + ':epa', place)

        def interpret(path):
            summarize(epa / 'epa_result.jplace', ref / 'metadata.tsv', alignment / 'queries.fasta', alignment / 'queries.aligned.fasta', path, .95)
            # Full weights remain in the raw jplace; avoid a redundant 14 MB copy per fold.
            (path / 'placement_weights.json.gz').unlink()
        summary = fold / '05_summary'; stage(summary, signature + ':summary', interpret)
        placement_path = summary / 'lsu_placements.copies.tsv'
    else:
        placement_path = FROZEN / 'placement/summary/lsu_placements.copies.tsv'

    def score(path):
        data = Data(FROZEN / 'inputs/benchmark'); ids = plan['query_ids']; drop = frozenset(plan['held_cultures'])
        if not set(ids) <= set(data.cop) or any(data.cop[x]['culture'] not in drop for x in ids):
            raise ValueError('Invalid fold queries or culture exclusions')
        old = {r['accession']: r for r in read_tsv(FROZEN / 'placement/summary/lsu_placements.copies.tsv')}
        new = {r['accession']: r for r in read_tsv(placement_path)}
        baseline = infer(data, ids, drop, old); experimental = infer(data, ids, drop, new)
        if plan['rank'] == 'species':
            expected = {r['accession']: r for r in read_tsv(FROZEN / 'comparison/epa_species_clades_physical_species.copies.tsv')}
            for r in baseline:
                for k in ['status', 'placement', 'rank', 'call_species', 'call_genera', 'correct', 'reasons', 's58']:
                    if r[k] != expected[r['accession']][k]:
                        raise AssertionError(f'Full-reference control parity failed: {r["accession"]} {k}')
        byid = {r['accession']: r for r in baseline}; changes = []
        for r in experimental:
            b = byid[r['accession']]
            changes.append({'accession': r['accession'], 'culture': r['culture'], 'copy_class': r['copy_class'], 'div_kind': r['div_kind'], 'old_status': b['status'], 'new_status': r['status'], 'old_placement': b['placement'], 'new_placement': r['placement'], 'old_LSU_state': b['LSU_state'], 'new_LSU_state': r['LSU_state']})
        write_tsv(path / 'full_reference.copies.tsv', baseline); write_tsv(path / 'external_holdout.copies.tsv', experimental); write_tsv(path / 'paired_changes.tsv', changes)
        write_json(path / 'metrics.json', {'plan': plan, 'full_reference': fold_metrics(baseline), 'external_holdout': fold_metrics(experimental), 'control_parity': 'passed' if plan['rank'] == 'species' else 'same-fold full-V18 control computed', 'independent_validation': False})
    stage(fold / '06_evaluation', signature + ':score', score)
    write_json(fold / 'DONE.json', {'signature': signature, 'complete': True, 'seconds_this_invocation': time.monotonic() - start, 'external_status': plan['external_status'], 'result_sha256': sha(fold / '06_evaluation/external_holdout.copies.tsv')})
    print(f"Completed {plan['target']}", flush=True)
    return plan['key']


def aggregate(output, rank):
    plans = load_plan(rank); completed = []; all_new = []; all_old = []; strata = defaultdict(list)
    for plan in plans:
        folder = Path(output) / plan['key']
        if not (folder / 'DONE.json').exists():
            continue
        config = json.loads((folder / 'config.json').read_text())
        if not valid_stage(folder / '06_evaluation', config['signature'] + ':score'):
            raise ValueError(f'Missing completed scoring stage: {folder}')
        if json.loads((folder / 'DONE.json').read_text())['result_sha256'] != sha(folder / '06_evaluation/external_holdout.copies.tsv'):
            raise ValueError('Completed result digest changed')
        completed.append(plan['target'])
        new = read_tsv(folder / '06_evaluation/external_holdout.copies.tsv'); old = read_tsv(folder / '06_evaluation/full_reference.copies.tsv')
        for r in new + old:
            r['target'] = plan['target']; r['external_status'] = plan['external_status']
        all_new.extend(new); all_old.extend(old); strata[plan['external_status']].extend(new)
    out = Path(output) / 'aggregate'; out.mkdir(exist_ok=True, parents=True)
    result = {'rank': rank, 'complete': len(completed) == len(plans), 'completed_folds': completed,
              'expected_folds': len(plans), 'missing_folds': [p['target'] for p in plans if p['target'] not in completed],
              'external_holdout': fold_metrics(all_new), 'full_reference': fold_metrics(all_old),
              'external_strata': {k: fold_metrics(v) for k, v in strata.items()},
              'independent_validation': False}
    write_tsv(out / f'{rank}.copies.tsv', all_new); write_tsv(out / f'{rank}.control.tsv', all_old)
    write_json(out / f'{rank}.summary.json', result)
    print(json.dumps({k: result[k] for k in ('rank', 'complete', 'completed_folds', 'expected_folds', 'missing_folds')}, indent=2))
    return result


def collect(output, destination):
    # Compact review evidence. Full alignments and exhaustive jplace stay on the PC.
    output = Path(output).resolve(); destination = Path(destination).resolve()
    names = ('config.json', 'DONE.json', 'stage.json', 'metrics.json', 'commands.json', 'validation.json',
             'full_reference.copies.tsv', 'external_holdout.copies.tsv', 'paired_changes.tsv',
             'lsu_placements.copies.tsv', 'placement_summary.json', 'fit.raxml.bestModel', 'fit.raxml.bestTree')
    files = sorted(p for p in output.rglob('*') if p.is_file() and (p.name in names or p.suffix == '.log' or 'aggregate' in p.parts))
    files = [p for p in files if p.resolve() != destination]
    with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in files: z.write(p, 'results/' + str(p.relative_to(output)))
        for name in ('PROTOCOL.md', 'FROZEN_M2B.json', 'holdout_plan.json', 'species_exclusions.tsv', 'external_holdout.py'):
            z.write(ROOT / name, 'protocol/' + name)
    print(f'Created {destination} ({destination.stat().st_size:,} bytes). No upload performed.')


def main():
    ap = argparse.ArgumentParser(description=__doc__); sub = ap.add_subparsers(dest='command', required=True)
    ls = sub.add_parser('list'); ls.add_argument('--rank', choices=['species', 'genus'], default='species')
    run = sub.add_parser('run'); run.add_argument('--rank', choices=['species', 'genus'], default='species'); run.add_argument('--target', action='append')
    run.add_argument('--out', type=Path, default=ROOT / 'runs'); run.add_argument('--threads', type=int, default=8); run.add_argument('--jobs', type=int, default=1)
    run.add_argument('--mafft', default='mafft'); run.add_argument('--raxml-ng', default='raxml-ng'); run.add_argument('--epa-ng', default='epa-ng')
    report = sub.add_parser('report'); report.add_argument('--out', type=Path, default=ROOT / 'runs'); report.add_argument('--rank', choices=['species', 'genus'], default='species')
    exp = sub.add_parser('collect'); exp.add_argument('--out', type=Path, default=ROOT / 'runs'); exp.add_argument('--zip', type=Path, required=True)
    a = ap.parse_args()
    if a.command == 'list':
        for p in load_plan(a.rank): print(f"{p['target']}\t{len(p['removed_ref_ids'])} tips\t{len(p['query_ids'])} copies\t{p['external_status']}")
        return
    if a.command == 'report': aggregate(a.out, a.rank); return
    if a.command == 'collect': collect(a.out, a.zip); return
    if a.threads < 1 or a.jobs < 1: ap.error('threads and jobs must be positive')
    frozen_hash = verify_frozen(); ts = tool_info(a); plans = load_plan(a.rank, a.target)
    base_key = digest({'frozen': frozen_hash, 'protocol': sha(ROOT / 'PROTOCOL.md'), 'plan': sha(ROOT / 'holdout_plan.json'),
                       'algorithm': 'fresh-reference-fixed-topology-v1', 'versions': {k: v['version'] for k, v in ts.items()}})
    a.out = a.out.resolve(); a.out.mkdir(parents=True, exist_ok=True)
    write_json(a.out / 'invocation.json', {'argv': sys.argv, 'base_key': base_key, 'script_sha256': sha(__file__), 'tools': ts, 'threads_per_job': a.threads, 'jobs': a.jobs, 'python': sys.version})
    if a.jobs == 1:
        for p in plans: run_fold(p, a.out, ts, a.threads, base_key)
    else:
        with ProcessPoolExecutor(max_workers=a.jobs) as pool:
            fs = [pool.submit(run_fold, p, a.out, ts, a.threads, base_key) for p in plans]
            for f in as_completed(fs): f.result()
    aggregate(a.out, a.rank)


if __name__ == '__main__':
    main()
