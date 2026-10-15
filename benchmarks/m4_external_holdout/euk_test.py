#!/usr/bin/env python3
"""EUKARYOME long-read test (benchmarks/eukaryome/PROTOCOL_EUK.md): the LongGloDB runner on another test set.

Imports lgdb_test.py, which stays byte-identical, and replaces only what names the test set:
  EukInputs  the EUKARYOME folder: protocol record (PROTOCOL_EUK.md, predictions.lock, names.tsv,
             names_decisions.tsv, overlap.v18_plan.tsv and the four test-window files, all in
             protocol_hashes.txt), the lock, overlap.asvs.tsv and overlap.cultures.tsv, the test windows
             (overlap.win.*.fasta, written by stage3/eukaryome_overlap.py roles under the read IDs) and the
             Stage 1 summary (stage1/windows/euk_test.summary.json: the options of the deposited copies).
  b_plan     the fixed overlap.v18_plan.tsv: 30 species rounds (the species' V18 records out) and 21 genus
             rounds (every V18 record of the genus out); genus rounds score every read of the genus.
  test_set   culture IDs 'organism | EUKARYOME' (one per labelled taxon; EUKARYOME has no culture codes).
Everything else is lgdb_test.py's own code: evidence with parity against the frozen VT and SH evidence,
condition A, the B folds, the M4 parity check, scoring with and without N6, primary and secondary VT analyses.

Commands (from amf_m4_step1, tools as for lgdb_test.py):
  python3 euk_test.py check   --euk ../bench/eukaryome --maarjam ../ref/maarjam/maarjam.fasta \\
                              --unite ../ref/unite/UNITE_public_19.02.2025.fasta.gz
  python3 euk_test.py run     (same options) --threads 4 --jobs 8
  python3 euk_test.py collect --zip ~/Downloads/euk_results.zip
Default output: runs_euk/ beside this file. If this wrapper changes after scores exist, delete runs_euk/06_scores
before rerunning: the stage keys hash lgdb_test.py, not this file.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import lgdb_test as L                                     # noqa: E402

EUK_PROTOCOL_FILES = ('PROTOCOL_EUK.md', 'predictions.lock', 'names.tsv', 'names_decisions.tsv', 'overlap.v18_plan.tsv',
                      'overlap.win.ssu_flank.fasta', 'overlap.win.full_its.fasta', 'overlap.win.s58.fasta',
                      'overlap.win.LSU_LROR-FLR2.fasta')
STAGE1_SUMMARY = 'stage1/windows/euk_test.summary.json'
LAB = 'EUKARYOME'
_STATE = {}


class EukInputs(L.Inputs):
    """lgdb_test.Inputs with the EUKARYOME folder layout; every check of the original is kept."""

    def __init__(self, euk, maarjam=None, unite=None):            # noqa: C901 (mirrors L.Inputs.__init__)
        self.dir = Path(euk).resolve()
        d = self.dir
        # 1. protocol record
        record = {}
        for line in (d / 'protocol_hashes.txt').read_text().splitlines():
            if line.strip():
                h, name = line.split(None, 1)
                record[Path(name.strip().lstrip('*')).name] = h
        for name in EUK_PROTOCOL_FILES:
            if name not in record:
                raise ValueError(f'{name} is not in protocol_hashes.txt; record it before running')
            if L.sha(d / name) != record[name]:
                raise ValueError(f'{name} changed after it was recorded')
        self.protocol_hashes = {n: record[n] for n in EUK_PROTOCOL_FILES}
        # 2. the lock binds the predictions, the roles and the reference
        lock = json.loads((d / 'predictions.lock').read_text())
        if L.sha(d / 'predictions.predictions.tsv') != lock['predictions_sha256']:
            raise ValueError('predictions.predictions.tsv differs from the lock')
        if L.sha(d / 'predictions.cultures.tsv') != lock['cultures_sha256']:
            raise ValueError('predictions.cultures.tsv differs from the lock')
        by_name = {Path(k).name: v for k, v in lock['inputs'].items()}
        if L.sha(d / 'overlap.asvs.tsv') != by_name['overlap.asvs.tsv']:
            raise ValueError('overlap.asvs.tsv differs from the lock')
        ref = L.FROZEN / 'inputs/benchmark/ref'
        for name in ('copies.tsv', 'cultures.tsv'):
            if L.sha(ref / name) != by_name[name]:
                raise ValueError(f'frozen {name} is not the culture reference the lock used')
        for name in ('full_its.fasta', 'LSU_LROR-FLR2.fasta', 'ssu_flank.fasta'):
            if L.sha(ref / 'fasta' / name) != by_name[name]:
                raise ValueError(f'frozen {name} is not the one the lock used')
        if lock['ref_fingerprint'] != L.REF_FINGERPRINT:
            raise ValueError('lock fingerprint is not culture reference v4')
        self.lock = lock
        # 3. test set and roles
        self.asvs = [r for r in L.read_tsv(d / 'overlap.asvs.tsv') if r['role'] in L.TEST_ROLES]
        self.pred = {r['asv']: r for r in L.read_tsv(d / 'predictions.predictions.tsv')}
        self.cultures = {(r['organism'], r['code']): r for r in L.read_tsv(d / 'overlap.cultures.tsv')}
        if {r['asv'] for r in self.asvs} != set(self.pred):
            raise ValueError('test reads differ between overlap.asvs.tsv and the predictions')
        # 4. test windows: Stage 1 with the options of the deposited copies, under the read IDs
        summ = json.loads((d / STAGE1_SUMMARY).read_text())['params']
        want = {'version': '0.3.1', 'max_err_frac': 0.15, 'only_anchored': True, 'open_ends': True, 'keep_primers': False}
        if any(summ.get(k) != v for k, v in want.items()):
            raise ValueError(f'Stage 1 options differ from those of the deposited copies: {summ}')
        self.window_paths = {m: d / f'overlap.win.{f}.fasta' for m, f in L.WINDOWS.items()}
        self.win = {m: L.read_fasta(p) for m, p in self.window_paths.items()}
        ids = {r['asv'] for r in self.asvs}
        for m, seqs in self.win.items():
            if set(seqs) != ids:
                raise ValueError(f'test window {m} does not hold exactly the {len(ids)} test reads')
        # 5. external databases: the files of the v4 crosswalk
        prior = json.loads((L.FROZEN / 'inputs/benchmark/baseline/crosswalk_v4.run.json').read_text())['arms']
        self.ext = {}
        for arm, path in (('VT', maarjam), ('SH', unite)):
            if path is None:
                continue
            path = Path(path).resolve()
            if L.sha(path) != prior[arm]['sha256']:
                raise ValueError(f'{path.name} is not the {arm} reference of the v4 crosswalk')
            self.ext[arm] = path

    def file_hashes(self):
        out = {f'euk/{n}': L.sha(self.dir / n) for n in ('overlap.asvs.tsv', 'overlap.cultures.tsv', 'predictions.predictions.tsv',
                                                          'predictions.cultures.tsv', 'protocol_hashes.txt', STAGE1_SUMMARY,
                                                          *EUK_PROTOCOL_FILES)}
        out.update({f'window/{m}': L.sha(p) for m, p in self.window_paths.items()})
        return out


def euk_test_set(inputs, v4_species):
    """lgdb_test.test_set with culture IDs 'organism | EUKARYOME'; the same role checks."""
    cultures, copies = {}, {}
    for r in inputs.asvs:
        key = (r['organism'], r['code'])
        cid = f"{r['organism']} | {LAB}"
        role = r['role']
        if role == 'known species' and r['organism'] not in v4_species:
            raise ValueError(f'known-species read whose name is not a v4 species: {r["asv"]}')
        if role in L.NOVEL_ROLES and r['organism'] in v4_species:
            raise ValueError(f'novel read whose name is a v4 species: {r["asv"]}')
        c = cultures.setdefault(cid, {'culture': cid, 'organism': r['organism'], 'code': r['code'], 'labs': set(),
                                      'collections': set(), 'role': role, 'stratum': set(), 'meta': inputs.cultures.get(key, {})})
        c['labs'].add(r['lab']); c['collections'].add(r['collection']); c['stratum'].add(r['stratum'])
        if c['role'] != role:
            raise ValueError(f'culture {cid} has reads in two roles')
        copies[r['asv']] = {'asv': r['asv'], 'culture': cid, 'role': role, 'organism': r['organism'], 'row': r}
    for c in cultures.values():
        c['lab'] = ';'.join(sorted(c.pop('labs'))); c['collection'] = ';'.join(sorted(c.pop('collections')))
        c['stratum'] = ';'.join(sorted(c['stratum']))
    if len(cultures) != len({k for k in inputs.cultures if inputs.cultures[k]['role'] in L.TEST_ROLES}):
        raise ValueError('test cultures differ from overlap.cultures.tsv')
    return cultures, copies


def euk_b_plan(test_cultures, test_copies):
    """The fixed plan of overlap.v18_plan.tsv, checked against the frozen V18 metadata and the roles."""
    d = _STATE['dir']
    meta = {r['ref_id']: r for r in L.read_tsv(L.FROZEN / 'reference/reference_metadata.tsv')}
    plans, mapping = [], []
    for p in L.read_tsv(d / 'overlap.v18_plan.tsv'):
        ids = sorted(x for x in p['ref_ids'].split(';') if x)
        if not ids or any(x not in meta or meta[x]['is_amf'] != '1' for x in ids):
            raise ValueError(f'plan {p["target"]}: removal holds unknown records or an outgroup')
        if p['round'] == 'species':
            role = 'novel species, genus in v4'
            q = sorted(x for x, t in test_copies.items() if t['role'] == role and t['organism'] == p['target'])
            if {meta[x]['label'] for x in ids} != set(p['v18_labels'].split(';')):
                raise ValueError(f'plan {p["target"]}: records and labels disagree')
        elif p['round'] == 'genus':
            role = 'novel lineage, genus not in v4'
            q = sorted(x for x, t in test_copies.items() if t['role'] == role and t['row']['stratum'] != 'placeholder lineage'
                       and t['organism'].split()[0] == p['target'])
            if any(meta[x]['genus'] != p['target'] for x in ids) or \
                    set(ids) != {x for x, m in meta.items() if m['is_amf'] == '1' and m['genus'] == p['target']}:
                raise ValueError(f'plan {p["target"]}: a genus round must remove every V18 record of the genus')
        else:
            raise ValueError(f'unknown round {p["round"]}')
        if len(q) != int(p['reads']):
            raise ValueError(f'plan {p["target"]}: {len(q)} test reads, plan says {p["reads"]}')
        cult = sorted({test_copies[x]['culture'] for x in q})
        plans.append({'key': f"euk_{p['round']}_" + ''.join(ch if ch.isalnum() else '_' for ch in p['target'].lower()),
                      'target': p['target'], 'round': p['round'], 'role': role, 'removed_ref_ids': ids,
                      'v18_label': p['v18_labels'], 'mapping': f"{p['round']} round", 'target_cultures': cult, 'query_ids': q})
        mapping.append({'test_species': p['target'], 'v18_label': p['v18_labels'], 'ref_ids': ';'.join(ids),
                        'records': len(ids), 'mapping': f"{p['round']} round"})
    apart = [r['target'] for r in L.read_tsv(d / 'overlap.v18_apart.tsv')] if (d / 'overlap.v18_apart.tsv').exists() else []
    # every novel read is in exactly one round or reported apart
    in_rounds = [x for p in plans for x in p['query_ids']]
    if len(in_rounds) != len(set(in_rounds)):
        raise ValueError('a test read is in two rounds')
    return plans, apart, mapping


def install(euk_dir):
    """Point lgdb_test at the EUKARYOME test set."""
    _STATE['dir'] = Path(euk_dir).resolve()
    L.Inputs = EukInputs
    L.PROTOCOL_FILES = EUK_PROTOCOL_FILES
    L.test_set = euk_test_set
    L.b_plan = euk_b_plan


def cmd_collect(a):
    L.cmd_collect(a)
    with zipfile.ZipFile(a.zip.resolve(), 'a', zipfile.ZIP_DEFLATED) as z:
        z.write(__file__, 'protocol/euk_test.py')
        d = _STATE['dir']
        for name in ('PROTOCOL_EUK.md', 'protocol_hashes.txt', 'predictions.lock', 'names.tsv', 'names_decisions.tsv',
                     'overlap.v18_plan.tsv', 'overlap.v18_apart.tsv', 'overlap.log'):
            if (d / name).exists():
                z.write(d / name, f'protocol/{name}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='command', required=True)
    for name in ('check', 'run'):
        p = sub.add_parser(name)
        p.add_argument('--euk', type=Path, required=True)
        p.add_argument('--maarjam', type=Path)
        p.add_argument('--unite', type=Path)
        p.add_argument('--out', type=Path, default=ROOT / 'runs_euk')
        p.add_argument('--threads', type=int, default=4)
        p.add_argument('--jobs', type=int, default=1)
        p.add_argument('--mafft', default='mafft'); p.add_argument('--raxml-ng', default='raxml-ng'); p.add_argument('--epa-ng', default='epa-ng')
    c = sub.add_parser('collect')
    c.add_argument('--euk', type=Path, required=True)
    c.add_argument('--out', type=Path, default=ROOT / 'runs_euk')
    c.add_argument('--zip', type=Path, required=True)
    a = ap.parse_args()
    install(a.euk)
    a.lgdb = a.euk                      # lgdb_test's commands read a.lgdb
    if a.command == 'check':
        L.cmd_check(a)
        tc, tcp = L.test_set(L.Inputs(a.euk), set(L.v1.Data(L.FROZEN / 'inputs/benchmark').species.values()))
        print('reads by role: ' + '; '.join(f'{k} {v}' for k, v in Counter(t['role'] for t in tcp.values()).items()))
    elif a.command == 'run':
        if a.threads < 1 or a.jobs < 1:
            ap.error('threads and jobs must be positive')
        L.cmd_run(a)
        print('next: the report (euk_report.py) reads runs_euk/06_scores', flush=True)
    else:
        cmd_collect(a)


if __name__ == '__main__':
    main()
