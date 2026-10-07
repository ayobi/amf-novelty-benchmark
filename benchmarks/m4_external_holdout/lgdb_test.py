#!/usr/bin/env python3
"""LongGloDB independent-culture test (PROTOCOL_LGDB.md, amendments 1 and 2).

Lives next to external_holdout.py in amf_m4_step1 and reuses the frozen code unchanged:
  Data, ArmTrainer, reconcile, evaluate          frozen_m2b (M2b package)
  infer, stage, tool_info, verify_frozen         external_holdout.py (protocol v1)
  ItsRule, genus_fields, QUANTILE, verify_v2     its_rule.py, external_holdout_v2.py (protocol v2)
  load_reference, run_blast, vt_hits             frozen arms.py / arm_sets.py (VT and SH evidence)
  summarize                                      frozen lsu_placement.py
The test ASVs are added to the frozen data in memory as extra copies of extra cultures, and every
test culture is in the excluded set of every call, so no test sequence trains anything. Before any
test result is computed the program proves this: it rescores all 37 M4 species folds with the test
cultures present and excluded, and requires the stored M4 v2 results back exactly. It also requires
the frozen VT and SH evidence of all 903 reference copies to be reproduced exactly by its own BLAST
run before it trusts the evidence it makes for the test ASVs.

Commands (run from any directory; tools as in M4 plus makeblastdb and blastn on PATH):
  python3 lgdb_test.py check  --lgdb ../bench/longglodb --maarjam ../ref/maarjam/maarjam.fasta \\
                              --unite ../ref/unite/UNITE_public_19.02.2025.fasta.gz
  python3 lgdb_test.py run    (same options) --threads 4 --jobs 8
  python3 lgdb_test.py report --lgdb ../bench/longglodb
  python3 lgdb_test.py collect --zip ~/Downloads/lgdb_results.zip
Default output: runs_lgdb/ beside this file. Completed stages are verified and reused; a stage
whose inputs or settings changed stops the run instead of being overwritten.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import external_holdout_v2 as v2                       # noqa: E402  (imports the frozen modules)
v1 = v2.v1
from Bio import Phylo                                   # noqa: E402
from its_rule import ItsRule                            # noqa: E402
from lsu_placement import read_fasta, write_fasta, read_tsv, write_tsv, sha, summarize  # noqa: E402
import arms as frozen_arms                              # noqa: E402  (frozen_m2b/code/baseline)
import arm_sets as frozen_arm_sets                      # noqa: E402

FROZEN = v1.FROZEN
TEST_ROLES = ('known species', 'novel species, genus in v4', 'novel lineage, genus not in v4', 'genus only')
NOVEL_ROLES = ('novel species, genus in v4', 'novel lineage, genus not in v4')
PROTOCOL_FILES = ('PROTOCOL_LGDB.md', 'predictions.lock', 'PROTOCOL_LGDB_amendment1.md', 'PROTOCOL_LGDB_amendment2.md')
WINDOWS = {'ssu_flank': 'ssu_flank', 'full_its': 'full_its', 's58': 's58', 'LSU_LROR-FLR2': 'LSU_LROR-FLR2'}
REF_FINGERPRINT = '181938d695f3df40'
# V18 labels that hold a test species under another spelling (as in the M4 species_exclusions.tsv)
V18_VARIANTS = {'Dentiscutata erythropa': ('Dentiscutata erythropus', 'orthographic variant')}
# amendment 1: the ASV that the paper's percentile definition adds to the locked 'inside' set
FROZEN_DEFINITION_EXTRA = ('Gigaspora_decipiens_AU102_9_KU_ASV1303_INVAM',)
STEFANI_LIKE_KEYS = ('status', 'placement', 'rank', 'call_species', 'call_genera', 'correct', 'reasons', 's58')


# ---------------------------------------------------------------- small helpers
def write_json(path, value):
    v1.write_json(path, value)


def digest(value):
    return v1.digest(value)


def tool_versions(args):
    """M4 tools (version-checked by the frozen tool_info) plus BLAST (version recorded)."""
    tools = v1.tool_info(args)
    for name in ('makeblastdb', 'blastn'):
        found = shutil.which(name)
        if not found:
            raise ValueError(f'Missing executable: {name} (put the amfsplit environment bin on PATH after amf-m4)')
        out = subprocess.run([found, '-version'], capture_output=True, text=True, timeout=30)
        tools[name] = {'executable': str(Path(found).resolve()), 'version': (out.stdout + out.stderr).strip().splitlines()[0],
                       'executable_sha256': sha(found)}
    return tools


# ---------------------------------------------------------------- inputs and their record
class Inputs:
    """Everything the test reads, verified against the protocol record and the lock."""

    def __init__(self, lgdb, maarjam=None, unite=None):
        self.dir = Path(lgdb).resolve()
        d = self.dir
        # 1. protocol record: every protocol file listed in protocol_hashes.txt and unchanged
        record = {}
        for line in (d / 'protocol_hashes.txt').read_text().splitlines():
            if line.strip():
                h, name = line.split(None, 1)
                record[name.strip().lstrip('*')] = h
        for name in PROTOCOL_FILES:
            if name not in record:
                raise ValueError(f'{name} is not in protocol_hashes.txt; record it before running')
            if sha(d / name) != record[name]:
                raise ValueError(f'{name} changed after it was recorded')
        self.protocol_hashes = {n: record[n] for n in PROTOCOL_FILES}
        # 2. the lock binds the predictions, the roles and the reference
        lock = json.loads((d / 'predictions.lock').read_text())
        if sha(d / 'predictions.predictions.tsv') != lock['predictions_sha256']:
            raise ValueError('predictions.predictions.tsv differs from the lock')
        if sha(d / 'predictions.cultures.tsv') != lock['cultures_sha256']:
            raise ValueError('predictions.cultures.tsv differs from the lock')
        by_name = {Path(k).name: v for k, v in lock['inputs'].items()}
        if sha(d / 'overlap.asvs.tsv') != by_name['overlap.asvs.tsv']:
            raise ValueError('overlap.asvs.tsv differs from the lock')
        ref = FROZEN / 'inputs/benchmark/ref'
        for name in ('copies.tsv', 'cultures.tsv'):
            if sha(ref / name) != by_name[name]:
                raise ValueError(f'frozen {name} is not the culture reference the lock used')
        for name in ('full_its.fasta', 'LSU_LROR-FLR2.fasta', 'ssu_flank.fasta'):
            if sha(ref / 'fasta' / name) != by_name[name]:
                raise ValueError(f'frozen {name} is not the one the lock used')
        if lock['ref_fingerprint'] != REF_FINGERPRINT:
            raise ValueError('lock fingerprint is not culture reference v4')
        self.lock = lock
        # 3. test set and roles
        self.asvs = [r for r in read_tsv(d / 'overlap.asvs.tsv') if r['role'] in TEST_ROLES]
        self.pred = {r['asv']: r for r in read_tsv(d / 'predictions.predictions.tsv')}
        self.cultures = {(r['organism'], r['code']): r for r in read_tsv(d / 'overlap.cultures.tsv')}
        if {r['asv'] for r in self.asvs} != set(self.pred):
            raise ValueError('test ASVs differ between overlap.asvs.tsv and the predictions')
        # 4. Stage 1 windows of the test ASVs, made with the options of the deposited copies
        w = d / 'stage1/windows'
        summ = json.loads((w / 'lgdb_test.summary.json').read_text())['params']
        want = {'version': '0.3.1', 'max_err_frac': 0.15, 'only_anchored': True, 'open_ends': True, 'keep_primers': False}
        if any(summ.get(k) != v for k, v in want.items()):
            raise ValueError(f'Stage 1 options differ from those of the deposited copies: {summ}')
        self.window_paths = {m: w / f'lgdb_test.{f}.fasta' for m, f in WINDOWS.items()}
        self.win = {m: read_fasta(p) for m, p in self.window_paths.items()}
        ids = {r['asv'] for r in self.asvs}
        for m, seqs in self.win.items():
            if set(seqs) != ids:
                raise ValueError(f'Stage 1 window {m} does not hold exactly the {len(ids)} test ASVs')
        # 5. external databases: the files of the v4 crosswalk
        prior = json.loads((FROZEN / 'inputs/benchmark/baseline/crosswalk_v4.run.json').read_text())['arms']
        self.ext = {}
        for arm, path in (('VT', maarjam), ('SH', unite)):
            if path is None:
                continue
            path = Path(path).resolve()
            if sha(path) != prior[arm]['sha256']:
                raise ValueError(f'{path.name} is not the {arm} reference of the v4 crosswalk')
            self.ext[arm] = path

    def file_hashes(self):
        out = {f'lgdb/{n}': sha(self.dir / n) for n in ('overlap.asvs.tsv', 'overlap.cultures.tsv', 'predictions.predictions.tsv',
                                                         'predictions.cultures.tsv', 'protocol_hashes.txt', *PROTOCOL_FILES)}
        out.update({f'window/{m}': sha(p) for m, p in self.window_paths.items()})
        return out


def test_set(inputs, v4_species):
    """Cultures and copies of the test, with the checks the roles imply. A culture is an organism and a
    culture code, as in overlap.cultures.tsv: a code sequenced by two labs is one culture."""
    cultures, copies = {}, {}
    for r in inputs.asvs:
        key = (r['organism'], r['code'])
        cid = f"{r['organism']} | LGDB {r['code']}"
        role = r['role']
        if role == 'known species' and r['organism'] not in v4_species:
            raise ValueError(f'known-species ASV whose name is not a v4 species: {r["asv"]}')
        if role in NOVEL_ROLES and r['organism'] in v4_species:
            raise ValueError(f'novel ASV whose name is a v4 species: {r["asv"]}')
        c = cultures.setdefault(cid, {'culture': cid, 'organism': r['organism'], 'code': r['code'], 'labs': set(), 'collections': set(),
                                      'role': role, 'meta': inputs.cultures.get(key, {})})
        c['labs'].add(r['lab']); c['collections'].add(r['collection'])
        if c['role'] != role:
            raise ValueError(f'culture {cid} has ASVs in two roles')
        copies[r['asv']] = {'asv': r['asv'], 'culture': cid, 'role': role, 'organism': r['organism'], 'row': r}
    for c in cultures.values():
        c['lab'] = ';'.join(sorted(c.pop('labs'))); c['collection'] = ';'.join(sorted(c.pop('collections')))
    if len(cultures) != len({k for k in inputs.cultures if inputs.cultures[k]['role'] in TEST_ROLES}):
        raise ValueError('test cultures differ from overlap.cultures.tsv')
    return cultures, copies


# ---------------------------------------------------------------- condition B plan
def b_plan(test_cultures, test_copies):
    """One fold per novel species whose label V18 holds: its V18 records out, all copies re-placed."""
    meta = read_tsv(FROZEN / 'reference/reference_metadata.tsv')
    species = sorted({c['organism'] for c in test_cultures.values() if c['role'] in NOVEL_ROLES})
    plans, absent, mapping = [], [], []
    for sp in species:
        label, kind = (sp, 'exact frozen label')
        if sp in V18_VARIANTS:
            label, kind = V18_VARIANTS[sp]
        ids = sorted(r['ref_id'] for r in meta if r['is_amf'] == '1' and r['label'] == label)
        if not ids and sp not in V18_VARIANTS:
            absent.append(sp)
        mapping.append({'test_species': sp, 'v18_label': label if ids else '', 'ref_ids': ';'.join(ids), 'records': len(ids),
                        'mapping': kind if ids else 'absent_by_label'})
        if ids:
            cult = sorted(c for c, v in test_cultures.items() if v['organism'] == sp)
            q = sorted(x for x, v in test_copies.items() if v['culture'] in cult)
            plans.append({'key': 'lgdb_' + ''.join(ch if ch.isalnum() else '_' for ch in sp.lower()), 'target': sp,
                          'role': test_cultures[cult[0]]['role'], 'removed_ref_ids': ids, 'v18_label': label,
                          'mapping': kind, 'target_cultures': cult, 'query_ids': q})
    return plans, absent, mapping


# ---------------------------------------------------------------- VT / SH evidence
def blast_with_coordinates(query, ref_fasta, workdir, name, task, threads, evalue=1e-10, max_targets=5000):
    """Same search as the frozen run_blast, with alignment coordinates added (secondary VT analysis)."""
    os.makedirs(workdir, exist_ok=True)
    db = os.path.join(workdir, name)
    with open(db + '.fasta', 'w') as fh:
        fh.write(ref_fasta)
    subprocess.run(['makeblastdb', '-in', db + '.fasta', '-dbtype', 'nucl', '-out', db], check=True, stdout=subprocess.DEVNULL)
    cols = 'qseqid sseqid pident length qlen slen evalue bitscore qstart qend sstart send'
    res = subprocess.run(['blastn', '-task', task, '-query', query, '-db', db, '-dust', 'no', '-outfmt', f'6 {cols}',
                          '-evalue', str(evalue), '-max_target_seqs', str(max_targets), '-max_hsps', '1',
                          '-num_threads', str(threads)], check=True, capture_output=True, text=True)
    hits = defaultdict(list)
    for line in res.stdout.splitlines():
        q, s, pid, ln, ql, sl, ev, bs, qs, qe, ss, se = line.split('\t')
        hits[q].append({'ref': s, 'pident': float(pid), 'length': int(ln), 'qlen': int(ql), 'slen': int(sl), 'evalue': float(ev),
                        'bitscore': float(bs), 'qstart': int(qs), 'qend': int(qe), 'sstart': int(ss), 'send': int(se)})
    return hits


def spanned_coverage(h, n_N):
    """Coverage with the reference bases that lie beyond the ends of the query discounted, as reference
    Ns already are (amendment 2, secondary analysis)."""
    s_lo, s_hi = min(h['sstart'], h['send']), max(h['sstart'], h['send'])
    if h['sstart'] <= h['send']:
        left = max(0, (s_lo - 1) - (h['qstart'] - 1))
        right = max(0, (h['slen'] - s_hi) - (h['qlen'] - h['qend']))
    else:
        left = max(0, (s_lo - 1) - (h['qlen'] - h['qend']))
        right = max(0, (h['slen'] - s_hi) - (h['qstart'] - 1))
    return h['length'] / min(h['qlen'], max(1, h['slen'] - n_N - left - right))


def vt_hits_spanned(hs, refs, min_cov, max_ev, min_id, drop):
    """frozen arm_sets.vt_hits with spanned_coverage in place of the frozen coverage; nothing else differs."""
    ids, best = {}, None
    for h in hs:
        info = refs[h['ref']]
        if frozen_arms.base_acc(info['acc']) in drop:
            continue
        cov = spanned_coverage(h, info.get('n_N', 0))
        if cov < min_cov or h['evalue'] > max_ev:
            continue
        v = info['label']
        if h['pident'] > ids.get(v, -1):
            ids[v] = h['pident']
        if h['pident'] >= min_id and (best is None or (h['bitscore'], h['pident']) > best):
            best = (h['bitscore'], h['pident'])
    return ids, (best[1] if best else None)


def best_unit(ids, hs, refs):
    """Highest-identity unit; ties go to the unit of the higher-bit-score hit, then the smaller code."""
    if not ids:
        return ('', 0.)
    top = max(ids.values())
    tied = sorted(u for u, p in ids.items() if p == top)
    if len(tied) > 1:
        score = {u: max((h['bitscore'] for h in hs if refs[h['ref']]['label'] == u and h['pident'] == top), default=-1) for u in tied}
        tied.sort(key=lambda u: (-score[u], u))
    return (tied[0], top)


def make_evidence(path, inputs, threads):
    """BLAST the reference and test windows; require the frozen evidence back; keep the test evidence."""
    data = v1.Data(FROZEN / 'inputs/benchmark')
    drop = {frozen_arms.base_acc(r['accession']) for r in read_tsv(FROZEN / 'inputs/benchmark/ref/copies.tsv')}
    out, report = {}, {}
    for arm in ('VT', 'SH'):
        kind, region, task, min_id, min_cov, _ = frozen_arm_sets.ARMS[arm]
        e = data.evidence[arm]
        keep = 'Glomeromycota' if kind == 'unite' else None
        refs, fasta, nrec = frozen_arms.load_reference(str(inputs.ext[arm]), kind, keep)
        if len(refs) != e['reference_records']:
            raise ValueError(f'{arm}: {len(refs)} reference records kept, frozen evidence used {e["reference_records"]}')
        ref_q = read_fasta(FROZEN / 'inputs/benchmark/ref/fasta' / f'{region}.fasta')
        test_q = inputs.win[region]
        if set(ref_q) & set(test_q):
            raise ValueError('test and reference query IDs overlap')
        qpath = path / f'{arm}.queries.fasta'
        write_fasta(qpath, {**ref_q, **test_q})
        hits = frozen_arms.run_blast(str(qpath), fasta, str(path / f'{arm}.blastdb'), arm.lower(), task, threads)
        # parity with the frozen evidence for every reference copy
        bad = []
        for x in data.cop:
            ids, anchor = frozen_arm_sets.vt_hits(hits.get(x, []), refs, min_cov, 1e-50, min_id, drop)
            q = e['queries'][x]
            if ids != q['unit_pident'] or anchor != q['anchor_pident']:
                bad.append(x)
        if bad:
            raise AssertionError(f'{arm}: frozen evidence not reproduced for {len(bad)} reference copies (first {bad[:3]}); '
                                 'the BLAST installation or the reference file differs from the one behind the benchmark')
        arm_out = {}
        for x in test_q:
            ids, anchor = frozen_arm_sets.vt_hits(hits.get(x, []), refs, min_cov, 1e-50, min_id, drop)
            arm_out[x] = {'unit_pident': ids, 'anchor_pident': anchor, 'best': best_unit(ids, hits.get(x, []), refs),
                          'raw_hits': len(hits.get(x, []))}
        out[arm] = arm_out
        report[arm] = {'reference_copies_reproduced': len(data.cop), 'test_queries': len(test_q),
                       'test_with_any_unit': sum(bool(v['unit_pident']) for v in arm_out.values()),
                       'test_with_anchor': sum(v['anchor_pident'] is not None for v in arm_out.values()),
                       'reference_records': len(refs), 'blast_task': task, 'min_id': min_id, 'min_cov': min_cov}
        if arm == 'VT':
            # secondary analysis (amendment 2): coverage over the part of each reference the query spans
            hc = blast_with_coordinates(str(qpath), fasta, str(path / 'VT2.blastdb'), 'vt2', task, threads)
            for x in set(hits) | set(hc):
                a = sorted((h['ref'], h['pident'], h['length'], h['bitscore']) for h in hits.get(x, []))
                b = sorted((h['ref'], h['pident'], h['length'], h['bitscore']) for h in hc.get(x, []))
                if a != b:
                    raise AssertionError(f'coordinate BLAST differs from the frozen search for {x}')
            sec, changed_ref = {}, 0
            for x in data.cop:
                if vt_hits_spanned(hc.get(x, []), refs, min_cov, 1e-50, min_id, drop) != frozen_arm_sets.vt_hits(hits.get(x, []), refs, min_cov, 1e-50, min_id, drop):
                    changed_ref += 1
            for x in test_q:
                ids, anchor = vt_hits_spanned(hc.get(x, []), refs, min_cov, 1e-50, min_id, drop)
                sec[x] = {'unit_pident': ids, 'anchor_pident': anchor, 'best': best_unit(ids, hc.get(x, []), refs)}
            out['VT_spanned'] = sec
            report['VT_spanned'] = {'test_with_any_unit': sum(bool(v['unit_pident']) for v in sec.values()),
                                    'test_with_anchor': sum(v['anchor_pident'] is not None for v in sec.values()),
                                    'reference_copies_whose_evidence_would_change': changed_ref}
        shutil.rmtree(path / f'{arm}.blastdb', ignore_errors=True)
        shutil.rmtree(path / 'VT2.blastdb', ignore_errors=True)
    with gzip.open(path / 'test_evidence.json.gz', 'wt') as fh:
        json.dump(out, fh, sort_keys=True)
    write_json(path / 'evidence_report.json', report)


# ---------------------------------------------------------------- augmented frozen data
def augmented_data(test_cultures, test_copies, windows, evidence, vt_key='VT'):
    """Frozen Data plus the test copies, in memory. Every frozen check runs first, on the frozen files."""
    data = v1.Data(FROZEN / 'inputs/benchmark')
    syn = frozen_arms.GENUS_SYNONYMS
    for cid, c in test_cultures.items():
        if cid in data.cul:
            raise ValueError(f'test culture ID collides with a reference culture: {cid}')
        genus = frozen_arms.genus_of(c['organism'])
        data.cul[cid] = {'culture': cid, 'organism': c['organism'], 'culture_label': f"LGDB:{c['code']}",
                         'genus': genus, 'named': str(int(frozen_arms.is_named(c['organism']))), 'name_status': 'ok'}
        data.genus[cid] = syn.get(genus, genus)
        if frozen_arms.is_named(c['organism']):
            data.species[cid] = c['organism']
            data.by_species[c['organism']].add(cid)
        data.by_label[data.cul[cid]['culture_label']].add(cid)
    data.cultures = frozenset(data.cul)
    for x, t in test_copies.items():
        if x in data.cop:
            raise ValueError(f'test copy ID collides with a reference copy: {x}')
        data.cop[x] = {'accession': x, 'culture': t['culture'], 'organism': t['organism'], 'copy_class': 'main', 'div_kind': '',
                       'status': 'retained', 's58_len': str(len(windows['s58'][x]))}
        data.by_culture[t['culture']].append(x)
        for m in ('ssu_flank', 'full_its', 's58', 'LSU_LROR-FLR2'):
            data.fa[m][x] = windows[m][x]
        for arm, key in (('VT', vt_key), ('SH', 'SH')):
            ev = evidence[key][x]
            data.evidence[arm]['queries'][x] = {'unit_pident': ev['unit_pident'], 'anchor_pident': ev['anchor_pident']}
            data.best[arm][x] = tuple(ev['best'])
    # test copies are never comparators: they stay out of data.mains and data.main58 (and every call drops them)
    return data


def placements_for(data, *tables):
    out = {}
    for t in tables:
        for r in read_tsv(t):
            out.setdefault(r['accession'], r)
    missing = set(data.cop) - set(out)
    if missing:
        raise ValueError(f'{len(missing)} copies have no placement (first {sorted(missing)[:3]})')
    return {x: out[x] for x in data.cop}


def score(data, rule, ids, drop, placements):
    """Frozen infer(), then N6 and the genus fields exactly as protocol v2 scores them."""
    train = data.cultures - drop
    before = v1.infer(data, ids, drop, placements)
    after = [dict(r) for r in before]
    rule.apply(after, {r['accession']: r for r in before}, train)
    for r in after:
        v2.genus_fields(data, r)
    return before, after


def parity_m4(data, rule, test_drop, placements, runs):
    """Rescore every M4 species fold with the test cultures present and excluded; require the stored
    M4 v2 control (frozen placements, protocol v2) back exactly."""
    plans = v1.load_plan('species')
    checked = 0
    for plan in plans:
        stored = Path(runs) / plan['key'] / '06_evaluation_v2/full_reference_v2.copies.tsv'
        if not stored.exists():
            raise ValueError(f'M4 v2 control missing for parity: {stored}')
        expected = {r['accession']: r for r in read_tsv(stored)}
        drop = frozenset(plan['held_cultures']) | test_drop
        _, rows = score(data, rule, plan['query_ids'], drop, placements)
        for r in rows:
            e = expected[r['accession']]
            for k in STEFANI_LIKE_KEYS + ('its_rule_fired', 'final_genus_call', 'final_genus_compatible'):
                if str(r[k]) != str(e[k]):
                    raise AssertionError(f'M4 parity failed: {plan["target"]} {r["accession"]} {k}: {r[k]!r} vs {e[k]!r}')
            checked += 1
    return {'folds': len(plans), 'copies': checked}


# ---------------------------------------------------------------- placement (A and B)
def align_and_place(fold, ref_aligned, tree, model, metadata, queries, tools, threads, key):
    """Stages 03-05 of an M4 fold, with the combined query file (903 reference + test windows)."""
    def align(path):
        shutil.copyfile(queries, path / 'queries.fasta')
        command = [tools['mafft']['executable'], '--auto', '--addfragments', str(path / 'queries.fasta'), '--keeplength', '--mapout',
                   '--thread', str(min(threads, 4)), str(ref_aligned)]
        v1.call(command, path / 'mafft.log', path / 'combined.fasta')
        raw = read_fasta(path / 'queries.fasta'); base = read_fasta(ref_aligned); full = read_fasta(path / 'combined.fasta')
        if set(raw) & set(base) or set(full) != set(base) | set(raw) or any(full[x] != s for x, s in base.items()):
            raise ValueError('Query addition altered the fixed reference or query IDs')
        write_fasta(path / 'queries.aligned.fasta', {x: full[x] for x in raw})
        (path / 'combined.fasta').unlink()
        write_json(path / 'commands.json', [command])
    alignment = fold / '03_queries'; v1.stage(alignment, key + ':query_alignment', align)

    def place(path):
        command = [tools['epa_ng']['executable'], '--ref-msa', str(ref_aligned), '--tree', str(tree), '--query',
                   str(alignment / 'queries.aligned.fasta'), '--model', str(model), '--outdir', str(path), '--threads', str(threads),
                   '--no-heur', '--no-pre-mask', '--filter-min-lwr', '0', '--filter-max', '2000', '--precision', '12']
        v1.call(command, path / 'epa.log')
        write_json(path / 'commands.json', [command])
    epa = fold / '04_epa'; v1.stage(epa, key + ':epa', place)

    def interpret(path):
        summarize(epa / 'epa_result.jplace', metadata, alignment / 'queries.fasta', alignment / 'queries.aligned.fasta', path, .95)
        (path / 'placement_weights.json.gz').unlink()
    summary = fold / '05_summary'; v1.stage(summary, key + ':summary', interpret)
    return summary / 'lsu_placements.copies.tsv', alignment / 'queries.aligned.fasta'


def run_a(out, queries, tools, threads, base_key):
    fold = out / 'A'; fold.mkdir(parents=True, exist_ok=True)
    key = digest({'base': base_key, 'condition': 'A', 'queries': sha(queries)})
    ref = FROZEN / 'reference'
    table, aligned = align_and_place(fold, ref / 'reference.aligned.fasta', ref / 'fit.raxml.bestTree', ref / 'fit.raxml.bestModel',
                                     ref / 'reference_metadata.tsv', queries, tools, threads, key)
    # diagnostic: does adding the test windows change how the 903 reference windows align and place?
    frozen_aln = read_fasta(FROZEN / 'placement/queries.aligned.fasta'); new_aln = read_fasta(aligned)
    frozen_pl = {r['accession']: r for r in read_tsv(FROZEN / 'placement/summary/lsu_placements.copies.tsv')}
    new_pl = {r['accession']: r for r in read_tsv(table)}
    write_json(fold / 'parity.json', {
        'reference_windows': len(frozen_aln),
        'alignment_identical': sum(new_aln[x] == s for x, s in frozen_aln.items()),
        'credible_species_clade_units_identical': sum(new_pl[x]['credible_species_clade_units'] == r['credible_species_clade_units']
                                                     for x, r in frozen_pl.items()),
        'eligibility_identical': sum(new_pl[x]['bridge_alignment_eligible'] == r['bridge_alignment_eligible'] for x, r in frozen_pl.items()),
        'note': 'Scores use the frozen placements of the 903 reference copies; this only shows whether the combined run treats them alike.'})
    return table


def run_b(plan, out, queries, tools, threads, base_key):
    fold = out / 'B' / plan['key']; fold.mkdir(parents=True, exist_ok=True)
    key = digest({'base': base_key, 'condition': 'B', 'plan': plan, 'queries': sha(queries)})
    removed = set(plan['removed_ref_ids'])
    print(f"Starting B {plan['target']} ({len(removed)} V18 records removed)", flush=True)

    def reference(path):
        meta = read_tsv(FROZEN / 'reference/reference_metadata.tsv')
        raw = read_fasta(FROZEN / 'reference/reference.unaligned.fasta')
        if not removed <= set(raw) or any(r['is_amf'] != '1' for r in meta if r['ref_id'] in removed):
            raise ValueError('Removal includes missing taxa or an outgroup')
        kept = {x: s for x, s in raw.items() if x not in removed}
        write_fasta(path / 'raw.fasta', kept)
        write_tsv(path / 'metadata.tsv', [r for r in meta if r['ref_id'] not in removed])
        tree = v1.pruned_tree(FROZEN / 'reference/backbone.topology.newick', removed, kept)
        Phylo.write(tree, path / 'topology.newick', 'newick', format_branch_length='%1.15f')
        command = [tools['mafft']['executable'], '--localpair', '--maxiterate', '1000', '--thread', str(threads), '--threadit', '0', str(path / 'raw.fasta')]
        v1.call(command, path / 'mafft.log', path / 'aligned.fasta')
        aligned = read_fasta(path / 'aligned.fasta')
        if set(aligned) != set(kept) or len({len(s) for s in aligned.values()}) != 1:
            raise ValueError('Invalid reference alignment')
        if any(aligned[x].replace('-', '').replace('.', '') != s.replace('-', '').replace('.', '') for x, s in kept.items()):
            raise ValueError('Reference residues changed')
        write_json(path / 'commands.json', [command])
    ref = fold / '01_reference'; v1.stage(ref, key + ':reference', reference)

    def model(path):
        command = [tools['raxml_ng']['executable'], '--evaluate', '--msa', str(ref / 'aligned.fasta'), '--tree', str(ref / 'topology.newick'),
                   '--model', 'GTR+G4', '--threads', str(min(threads, 4)), '--seed', '20260929', '--prefix', str(path / 'fit'), '--force', 'perf_threads']
        v1.call(command, path / 'raxml.log')
        if v1.splits(ref / 'topology.newick') != v1.splits(path / 'fit.raxml.bestTree'):
            raise ValueError('Fixed topology changed')
        write_json(path / 'commands.json', [command])
    model_dir = fold / '02_model'; v1.stage(model_dir, key + ':model', model)

    table, _ = align_and_place(fold, ref / 'aligned.fasta', model_dir / 'fit.raxml.bestTree', model_dir / 'fit.raxml.bestModel',
                               ref / 'metadata.tsv', queries, tools, threads, key)
    write_json(fold / 'DONE.json', {'key': key, 'plan': plan, 'placements_sha256': sha(table)})
    print(f"Completed B {plan['target']}", flush=True)
    return plan['key']


# ---------------------------------------------------------------- scoring
def annotate(rows, condition, analysis, fold, test_cultures, test_copies):
    for r in rows:
        t = test_copies[r['accession']]; c = test_cultures[t['culture']]
        r.update(condition=condition, analysis=analysis, fold=fold, test_role=t['role'], test_code=c['code'], test_lab=c['lab'],
                 test_collection=c['collection'])
    return rows


def run_scores(path, out, inputs, plans, runs):
    """Parity first; then A and every B fold, primary and secondary VT analyses, before and after N6."""
    data_ref = v1.Data(FROZEN / 'inputs/benchmark')
    v4_species = set(data_ref.species.values())
    test_cultures, test_copies = test_set(inputs, v4_species)
    with gzip.open(out / '01_evidence/test_evidence.json.gz', 'rt') as fh:
        evidence = json.load(fh)
    for arm in evidence.values():
        for v in arm.values():
            v['best'] = tuple(v['best'])
    windows = inputs.win
    data = {'primary': augmented_data(test_cultures, test_copies, windows, evidence, 'VT'),
            'secondary_vt': augmented_data(test_cultures, test_copies, windows, evidence, 'VT_spanned')}
    rule = ItsRule(data['primary'], v2.QUANTILE)
    drop = frozenset(test_cultures)
    a_table = out / 'A/05_summary/lsu_placements.copies.tsv'
    frozen_table = FROZEN / 'placement/summary/lsu_placements.copies.tsv'
    pa = placements_for(data['primary'], frozen_table, a_table)
    # parity: the M4 v2 results back, with the test cultures present and excluded
    t0 = time.monotonic()
    parity = parity_m4(data['primary'], rule, drop, pa, runs)
    parity['seconds'] = round(time.monotonic() - t0, 1)
    print(f"M4 parity: {parity['copies']} copies in {parity['folds']} folds reproduced exactly", flush=True)
    ids_all = sorted(test_copies)
    tables = defaultdict(list)
    for analysis, d in data.items():
        rule_d = rule if analysis == 'primary' else ItsRule(d, v2.QUANTILE)
        before, after = score(d, rule_d, ids_all, drop, pa)
        tables[f'A_{analysis}_v1'] += annotate(before, 'A', analysis, 'A', test_cultures, test_copies)
        tables[f'A_{analysis}_v2'] += annotate(after, 'A', analysis, 'A', test_cultures, test_copies)
        for plan in plans:
            fold_table = out / 'B' / plan['key'] / '05_summary/lsu_placements.copies.tsv'
            pb = placements_for(d, fold_table)
            before, after = score(d, rule_d, plan['query_ids'], drop, pb)
            tables[f'B_{analysis}_v1'] += annotate(before, 'B', analysis, plan['key'], test_cultures, test_copies)
            tables[f'B_{analysis}_v2'] += annotate(after, 'B', analysis, plan['key'], test_cultures, test_copies)
    for name, rows in tables.items():
        write_tsv(path / f'{name}.copies.tsv', sorted(rows, key=lambda r: (r['fold'], r['accession'])))
    t, n = rule.threshold(data['primary'].cultures - drop)
    write_json(path / 'scores.json', {'parity_m4': parity, 'its_threshold_pct': t, 'its_threshold_n': n, 'test_cultures': len(test_cultures),
                                      'test_copies': len(test_copies), 'b_folds': [p['key'] for p in plans],
                                      'tables': {k: len(v) for k, v in tables.items()}})


# ---------------------------------------------------------------- commands
def prepare(a):
    inputs = Inputs(a.lgdb, a.maarjam, a.unite)
    frozen_hash = v1.verify_frozen(); v2_key = v2.verify_v2()
    data_ref = v1.Data(FROZEN / 'inputs/benchmark')
    test_cultures, test_copies = test_set(inputs, set(data_ref.species.values()))
    plans, absent, mapping = b_plan(test_cultures, test_copies)
    return inputs, frozen_hash, v2_key, test_cultures, test_copies, plans, absent, mapping


def cmd_check(a):
    inputs, frozen_hash, v2_key, tc, tcp, plans, absent, mapping = prepare(a)
    roles = Counter(v['role'] for v in tcp.values())
    print('frozen code and protocol v2: verified')
    print('protocol record: ' + ', '.join(f'{n} {h[:8]}' for n, h in inputs.protocol_hashes.items()))
    print(f'test set: {len(tcp)} ASVs in {len(tc)} cultures; ' + '; '.join(f'{k} {v}' for k, v in roles.items()))
    print(f'condition B: {len(plans)} folds; V18 holds no label for: {", ".join(absent) or "none"}')
    for m in mapping:
        print(f"  {m['test_species']:28s} {m['records']} records  {m['mapping']}{'  (' + m['v18_label'] + ')' if m['v18_label'] and m['v18_label'] != m['test_species'] else ''}")
    if a.maarjam and a.unite:
        print('external databases: the v4 crosswalk files')
    try:
        tools = tool_versions(a)
        print('tools: ' + ', '.join(f"{k} {v['version'][:40]}" for k, v in tools.items()))
    except ValueError as err:
        print(f'tools: {err}')
    runs = ROOT / 'runs'
    have = sum((runs / p['key'] / '06_evaluation_v2/full_reference_v2.copies.tsv').exists() for p in v1.load_plan('species'))
    print(f'M4 v2 species controls for the parity check: {have} of 37')


def cmd_run(a):
    inputs, frozen_hash, v2_key, tc, tcp, plans, absent, mapping = prepare(a)
    if not (a.maarjam and a.unite):
        raise SystemExit('run needs --maarjam and --unite')
    tools = tool_versions(a)
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    base_key = digest({'frozen': frozen_hash, 'v2': v2_key, 'protocol': inputs.protocol_hashes,
                       'versions': {k: v['version'] for k, v in tools.items()}})
    write_json(out / 'invocation.json', {'argv': sys.argv, 'base_key': base_key, 'script_sha256': sha(__file__), 'tools': tools,
                                         'inputs': inputs.file_hashes(), 'threads_per_job': a.threads, 'jobs': a.jobs,
                                         'python': sys.version, 'started': time.strftime('%Y-%m-%d %H:%M:%S')})
    write_json(out / 'plan.json', {'b_folds': plans, 'absent_by_label': absent})
    write_tsv(out / 'v18_mapping.tsv', mapping)
    # 1. evidence, with parity against the frozen evidence
    ev_key = digest({'base': base_key, 'inputs': inputs.file_hashes(),
                     'ext': {k: sha(p) for k, p in inputs.ext.items()}, 'code': sha(__file__)})
    v1.stage(out / '01_evidence', ev_key + ':evidence', lambda p: make_evidence(p, inputs, a.threads * a.jobs))
    print('evidence: frozen VT and SH evidence reproduced for all reference copies', flush=True)
    # 2. combined query file: the 903 frozen reference windows and the test windows
    qdir = out / '02_queries'; qdir.mkdir(exist_ok=True)
    queries = qdir / 'queries.fasta'
    ref_q = read_fasta(FROZEN / 'placement/queries.fasta'); test_q = inputs.win['LSU_LROR-FLR2']
    combined = {**ref_q, **test_q}
    if queries.exists() and read_fasta(queries) != combined:
        raise ValueError('existing combined query file differs; choose a new output')
    if not queries.exists():
        write_fasta(queries, combined)
    # 3. condition A, then the B folds in parallel
    run_a(out, queries, tools, a.threads, base_key)
    print('condition A placed', flush=True)
    if a.jobs == 1:
        for p in plans:
            run_b(p, out, queries, tools, a.threads, base_key)
    else:
        with ProcessPoolExecutor(max_workers=a.jobs) as pool:
            for f in as_completed([pool.submit(run_b, p, out, queries, tools, a.threads, base_key) for p in plans]):
                f.result()
    # 4. scores (parity with the M4 v2 results first)
    sc_key = digest({'base': base_key, 'evidence': json.loads((out / '01_evidence/stage.json').read_text())['outputs'],
                     'A': sha(out / 'A/05_summary/lsu_placements.copies.tsv'),
                     'B': {p['key']: sha(out / 'B' / p['key'] / '05_summary/lsu_placements.copies.tsv') for p in plans},
                     'code': sha(__file__)})
    v1.stage(out / '06_scores', sc_key + ':scores', lambda p: run_scores(p, out, inputs, plans, ROOT / 'runs'))
    print(f'scores written to {out / "06_scores"}; next: python3 {Path(__file__).name} report --lgdb {a.lgdb}', flush=True)


def cmd_collect(a):
    out = a.out.resolve(); dest = a.zip.resolve()
    names = ('invocation.json', 'plan.json', 'v18_mapping.tsv', 'stage.json', 'commands.json', 'evidence_report.json', 'parity.json',
             'DONE.json', 'placement_summary.json', 'lsu_placements.copies.tsv', 'scores.json')
    files = sorted(p for p in out.rglob('*') if p.is_file() and p.resolve() != dest and
                   (p.name in names or p.suffix == '.log' or '06_scores' in p.parts or 'report' in p.parts or p.name == 'test_evidence.json.gz'))
    with zipfile.ZipFile(dest, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in files:
            z.write(p, 'results/' + str(p.relative_to(out)))
        z.write(__file__, 'protocol/lgdb_test.py')
    print(f'Created {dest} ({dest.stat().st_size:,} bytes). No upload performed.')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='command', required=True)
    for name in ('check', 'run'):
        p = sub.add_parser(name)
        p.add_argument('--lgdb', type=Path, required=True)
        p.add_argument('--maarjam', type=Path)
        p.add_argument('--unite', type=Path)
        p.add_argument('--out', type=Path, default=ROOT / 'runs_lgdb')
        p.add_argument('--threads', type=int, default=4)
        p.add_argument('--jobs', type=int, default=1)
        p.add_argument('--mafft', default='mafft'); p.add_argument('--raxml-ng', default='raxml-ng'); p.add_argument('--epa-ng', default='epa-ng')
    r = sub.add_parser('report')
    r.add_argument('--lgdb', type=Path, required=True)
    r.add_argument('--out', type=Path, default=ROOT / 'runs_lgdb')
    c = sub.add_parser('collect')
    c.add_argument('--out', type=Path, default=ROOT / 'runs_lgdb')
    c.add_argument('--zip', type=Path, required=True)
    a = ap.parse_args()
    if a.command == 'check':
        cmd_check(a)
    elif a.command == 'run':
        if a.threads < 1 or a.jobs < 1:
            ap.error('threads and jobs must be positive')
        cmd_run(a)
    elif a.command == 'report':
        import lgdb_report
        lgdb_report.report(Inputs(a.lgdb), a.out.resolve())
    else:
        cmd_collect(a)


if __name__ == '__main__':
    main()
