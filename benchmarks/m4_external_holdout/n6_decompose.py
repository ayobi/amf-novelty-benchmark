#!/usr/bin/env python3
"""Where does the one-copy-per-culture sensitivity of the ITS distance rule (N6) come from?

Reimplements the frozen rule (benchmarks/m4_external_holdout/its_rule.py, release 0.3.1) on one
pairwise full-ITS distance matrix, checks it against the published per-row values, then varies one
thing at a time:

  V0   frozen rule: copy-weighted exclusive 99th percentile, every copy in calibration and query pool
  V1c  same calibration values, inverse-CDF 99th percentile (estimator change only)
  V1   culture-weighted inverse-CDF 99th percentile (each calibration culture weighs 1)
  V2   one main copy per training culture in the calibration only (query pool keeps every copy)
  V3   frozen threshold, one main copy per training culture in the query pool only
  V4   one copy per culture in both (the 31 October experiment)

Scored on the three published designs: culture holdout (cost: correct known calls demoted), species
holdout with the LSU reference intact, and the LSU-removal test (benefit: false-known calls removed).
"""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import random
import statistics
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parent        # benchmarks/m4_external_holdout
sys.path.insert(0, str(ROOT))
import external_holdout as v1          # noqa: E402  (sets up the frozen_m2b imports)
from arm_sets import dist              # noqa: E402  (the frozen distance)

OUT = ROOT / 'audit_outputs'
OUT.mkdir(exist_ok=True)
REPS, SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 100, 20261031


def load(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


from lsu_placement import read_fasta, read_tsv   # noqa: E402
from arms import is_named                         # noqa: E402
import hashlib                                    # noqa: E402


FASTA = ROOT / 'frozen_m2b/inputs/benchmark/ref/fasta'
if not FASTA.exists():
    sys.exit(f'{FASTA} is missing: copy reference/culture_ref_v4/fasta there (the hash check below confirms it)')


class LiteData:
    """The parts of the frozen bridge.Data the ITS rule reads (the arm evidence is not needed),
    with the same hash check on the reference."""
    def __init__(self, inputs):
        ref = Path(inputs) / 'ref'
        manifest = json.loads((ref / 'MANIFEST.json').read_text())
        for name, expected in manifest['outputs'].items():
            if hashlib.sha256((ref / name).read_bytes()).hexdigest() != expected:
                raise ValueError(f'Frozen reference hash mismatch: {name}')
        self.cop = {r['accession']: r for r in read_tsv(str(ref / 'copies.tsv')) if r['status'] == 'retained'}
        self.cul = {r['culture']: r for r in read_tsv(str(ref / 'cultures.tsv'))}
        self.cultures = frozenset(self.cul)
        self.species = {c: r['organism'] for c, r in self.cul.items()
                        if is_named(r['organism']) and r.get('name_status', 'ok') == 'ok'}
        self.by_species = defaultdict(set)
        self.by_culture = defaultdict(list)
        for c, sp in self.species.items():
            self.by_species[sp].add(c)
        for x, r in self.cop.items():
            self.by_culture[r['culture']].append(x)
        self.fa = {'full_its': read_fasta(str(ref / 'fasta' / 'full_its.fasta'))}


data = LiteData(v1.FROZEN / 'inputs/benchmark')
its = data.fa['full_its']
acc = sorted(x for x in data.cop if x in its)
ix = {a: i for i, a in enumerate(acc)}
N = len(acc)
mains = {x for x in acc if data.cop[x]['copy_class'] == 'main'}
culture_of = {x: data.cop[x]['culture'] for x in acc}
copies_of = defaultdict(list)            # culture -> copies with full ITS
for x in acc:
    copies_of[culture_of[x]].append(x)
mains_of = {c: [x for x in xs if x in mains] for c, xs in copies_of.items()}

dpath = OUT / 'D.npy'
if dpath.exists():
    D = np.load(dpath)
else:
    t0 = time.time()
    D = np.zeros((N, N))
    for i in range(N):
        a = its[acc[i]]
        for j in range(i + 1, N):
            D[i, j] = D[j, i] = dist(a, its[acc[j]])     # sorted order, first argument on ties, as frozen
    np.save(dpath, D)
    print(f'distance matrix {N}x{N} in {time.time() - t0:.0f}s', flush=True)

cultures = sorted(copies_of)
ci = {c: k for k, c in enumerate(cultures)}
# copy -> culture minima, frozen semantics: pool 'all' (every copy) and 'main' (main-type copies)
DALL = np.full((N, len(cultures)), np.inf)
DMAIN = np.full((N, len(cultures)), np.inf)
for c, xs in copies_of.items():
    cols = [ix[x] for x in xs]
    DALL[:, ci[c]] = D[:, cols].min(axis=1)
    mcols = [ix[x] for x in mains_of[c]]
    if mcols:
        DMAIN[:, ci[c]] = D[:, mcols].min(axis=1)
for x in acc:                            # a copy is never its own nearest
    DALL[ix[x], ci[culture_of[x]]] = DMAIN[ix[x], ci[culture_of[x]]] = np.inf
    # (own culture is excluded from both the threshold and the query statistic below)

species_sets = {sp: set(cs) for sp, cs in data.by_species.items()}


def calib(train, sample=None):
    """[(value, culture)] for main-type copies of named species with >= 2 training cultures."""
    vals = []
    for sp, cs in species_sets.items():
        tc = [c for c in cs & train]
        if len(tc) < 2:
            continue
        for c in tc:
            others = [ci[o] for o in tc if o != c]
            ys = ([sample[c]] if sample and c in sample else []) if sample is not None else mains_of.get(c, [])
            for y in ys:
                if sample is not None:
                    ds = [D[ix[y], ix[sample[o]]] for o in tc if o != c and o in sample]
                else:
                    ds = [d for d in DMAIN[ix[y], others] if np.isfinite(d)]
                if ds:
                    vals.append((min(ds), c))
    return vals


def q_exclusive(vals, q=99):
    v = [a for a, _ in vals]
    return statistics.quantiles(v, n=100)[q - 1] if len(v) >= 2 else None


def q_weighted(vals, q=0.99, by_culture=True):
    if len(vals) < 2:
        return None
    per = Counter(c for _, c in vals)
    w = [(a, (1 / per[c]) if by_culture else 1.0) for a, c in vals]
    w.sort()
    tot, acc_w = sum(b for _, b in w), 0.0
    for a, b in w:
        acc_w += b
        if acc_w >= q * tot - 1e-12:
            return a
    return w[-1][0]


def nearest(x, train, sample=None):
    own = culture_of[x]
    if sample is None:
        cols = [ci[c] for c in train if c != own]
        d = DALL[ix[x], cols].min() if cols else np.inf
    else:
        ds = [D[ix[x], ix[sample[c]]] for c in train if c != own and c in sample]
        d = min(ds) if ds else np.inf
    return None if not np.isfinite(d) else d


def held_of(row):
    return frozenset(h for h in row['held_cultures'].split(';') if h) if row.get('held_cultures') else None


# --- the three designs -------------------------------------------------------------------------
plan = {f['target']: frozenset(f['held_cultures']) for f in
        json.loads((ROOT / 'holdout_plan.json').read_text())['folds'] if f['rank'] == 'species'}
cult = load(ROOT / 'frozen_m2b/comparison/epa_species_clades_physical_culture.copies.tsv')
spec = load(ROOT / 'frozen_m2b/comparison/epa_species_clades_physical_species.copies.tsv')
m4v1 = load(ROOT / 'aggregate_v2/species.v1.copies.tsv')
m4v2 = {r['accession'] + '|' + r['target']: r for r in load(ROOT / 'aggregate_v2/species.v2.copies.tsv')}
for r in m4v1:
    r['_held'] = plan[r['target']]
for rows in (cult, spec):
    for r in rows:
        r['_held'] = held_of(r)
all_cultures = frozenset(cultures) | frozenset(data.cultures)


def folds(rows):
    g = defaultdict(list)
    for r in rows:
        g[r['_held']].append(r)
    return g


DESIGNS = {'culture': folds(cult), 'species': folds(spec), 'm4': folds(m4v1)}
DEV = {'Glomus chinense', 'Glomus rugosae'}


def score(variant_threshold, variant_sample_pool, rng=None):
    """variant_threshold(train, sample) -> T; variant_sample_pool: use the sample for the query pool."""
    out = {}
    sample = None
    if rng is not None:
        sample = {c: rng.choice(ms) for c, ms in sorted(mains_of.items()) if ms}
    for name, g in DESIGNS.items():
        fired, res = [], Counter()
        for held, rows in g.items():
            train = frozenset(data.cultures) - held
            T = variant_threshold(train, sample)
            if variant_sample_pool:
                tc = [c for c in sorted(train) if c in sample]
                cols, names = np.array([ix[sample[c]] for c in tc]), np.array(tc, dtype=object)
                M = D
            else:
                tc = [c for c in sorted(train) if c in ci]
                cols, names = np.array([ci[c] for c in tc]), np.array(tc, dtype=object)
                M = DALL
            for r in rows:
                x = r['accession']
                d = None
                if x in ix and len(cols):
                    keep = names != culture_of[x]
                    if keep.any():
                        v = M[ix[x], cols[keep]].min()
                        d = v if np.isfinite(v) else None
                known = r['status'].startswith('known')
                fire = bool(known and d is not None and T is not None and d > T)
                main = r['copy_class'] == 'main'
                if name == 'culture':
                    ok = known and r['correct'] in ('yes', 'contains') and main
                    res['kept_before'] += ok
                    res['kept_after'] += ok and not fire
                else:
                    named_main = main and is_named(r['organism'])      # the paper's 803 named-species copies
                    res['fk_before'] += known and named_main
                    res['fk_after'] += known and named_main and not fire
                    if name == 'm4':
                        tag = ('dev' if r['target'] in DEV else 'callosa' if r['target'].startswith('Ambispora callosa')
                               else 'test34')
                        res[f'{tag}_before'] += known and main
                        res[f'{tag}_after'] += known and main and not fire
                        if fire:
                            fired.append(x + '|' + r['target'])
        out[name] = res
        if name == 'm4':
            out['m4_fired'] = fired
    return out


def thr_frozen(train, sample):
    return q_exclusive(calib(train))


_cache = {}


def cached(key, fn):
    if key not in _cache:
        _cache[key] = fn()
    return _cache[key]


def V0(train, sample):
    return cached(('V0', train), lambda: q_exclusive(calib(train)))


def V1c(train, sample):
    return cached(('V1c', train), lambda: q_weighted(calib(train), by_culture=False))


def V1(train, sample):
    return cached(('V1', train), lambda: q_weighted(calib(train), by_culture=True))


def onecopy(train, sample):
    return q_exclusive(calib(train, sample))


if __name__ == '__main__':
    t0 = time.time()
    # parity with the published LSU-removal rows (threshold, nearest distance, firing)
    bad, n81 = 0, 0
    for held, rows in DESIGNS['m4'].items():
        train = frozenset(data.cultures) - held
        T = V0(train, None)
        for r in rows:
            p = m4v2[r['accession'] + '|' + r['target']]
            d = nearest(r['accession'], train)
            fire = r['status'].startswith('known') and d is not None and d > T
            n81 += fire
            if (round(T, 3) != float(p['its_threshold_pct']) or
                    (p['its_nearest_training_pct'] and round(d, 3) != float(p['its_nearest_training_pct'])) or
                    int(fire) != int(p['its_rule_fired'])):
                bad += 1
    print(f'parity vs published v2 rows: {len(m4v1)} rows, {bad} mismatches, {n81} firings', flush=True)

    res = {}
    for name, fn in (('V0 frozen', V0), ('V1c inverse-CDF q99', V1c), ('V1 culture-weighted q99', V1)):
        res[name] = score(fn, False)
        print(name, {k: dict(v) if isinstance(v, Counter) else len(v) for k, v in res[name].items()}, flush=True)
    orig81 = set(res['V0 frozen']['m4_fired'])

    def reps(fn, pool):
        rng = random.Random(SEED)
        rows = []
        for _ in range(REPS):
            s = score(fn, pool, rng)
            rows.append({'kept': s['culture']['kept_after'], 'fk_species': s['species']['fk_after'],
                         'fk_m4': s['m4']['fk_after'], 'fk_test34': s['m4']['test34_after'],
                         'refire': len(orig81 & set(s['m4_fired'])), 'fired': len(s['m4_fired'])})
        return rows

    summary = {}
    for name, fn, pool in (('V2 one-copy calibration only', onecopy, False),
                           ('V3 one-copy query pool only', V0, True),
                           ('V4 one-copy both', onecopy, True)):
        rows = reps(fn, pool)
        summary[name] = {k: (min(r[k] for r in rows), statistics.median(r[k] for r in rows), max(r[k] for r in rows))
                         for k in rows[0]}
        print(name, summary[name], flush=True)

    # how the one-copy threshold is made: calibration size and position against the largest values
    full = frozenset(data.cultures)
    allc = calib(full)
    rng = random.Random(SEED)
    pos = []
    for _ in range(REPS):
        smp = {c: rng.choice(ms) for c, ms in sorted(mains_of.items()) if ms}
        vs = sorted(a for a, _ in calib(full, smp))
        T = statistics.quantiles(vs, n=100)[98]
        pos.append((len(vs), T, vs[-1], vs[-2]))
    jack = []
    for c in sorted({c for _, c in allc}):
        vs = [a for a, cc in allc if cc != c]
        jack.append(statistics.quantiles(vs, n=100)[98])
    extra = {
        'all_training': {'n': len(allc), 'T': q_exclusive(allc), 'culture_weighted_T': q_weighted(allc),
                         'cultures': len({c for _, c in allc}),
                         'values_above_T': sum(a > q_exclusive(allc) for a, _ in allc)},
        'one_copy_all_training': {'n': pos[0][0], 'T_min': min(p[1] for p in pos), 'T_median': statistics.median(p[1] for p in pos),
                                  'T_max': max(p[1] for p in pos),
                                  'share_T_above_second_largest': sum(p[1] >= p[3] for p in pos) / len(pos)},
        'jackknife_one_culture_out': {'T_min': min(jack), 'T_max': max(jack)},
    }
    print(json.dumps(extra, indent=1, default=float), flush=True)
    (OUT / 'n6_decompose_results.json').write_text(json.dumps(
        {'fixed': {k: {kk: dict(vv) if isinstance(vv, Counter) else sorted(vv) for kk, vv in v.items()} for k, v in res.items()},
         'replicates': summary, 'threshold': extra, 'reps': REPS, 'seed': SEED}, indent=1, default=float))
    table = [('Species clades, no ITS rule', 751, 71, 201, 78, '')]
    f = res['V0 frozen']; table.append(('Frozen ITS rule (paper, V0)', f['culture']['kept_after'], f['species']['fk_after'],
                                        f['m4']['fk_after'], f['m4']['test34_after'], '81 of 81'))
    for name in ('V1 culture-weighted q99', 'V1c inverse-CDF q99'):
        f = res[name]
        table.append((name, f['culture']['kept_after'], f['species']['fk_after'], f['m4']['fk_after'], f['m4']['test34_after'],
                      f"{len(orig81 & set(f['m4_fired']))} of 81 ({len(f['m4_fired'])} fired)"))
    for name, v in summary.items():
        med = lambda k: v[k][1]
        table.append((name + f' (median of {REPS})', med('kept'), med('fk_species'), med('fk_m4'), med('fk_test34'),
                      f"{med('refire')} of 81 (range {v['refire'][0]} to {v['refire'][2]})"))
    with open(OUT / 'n6_decompose_table.tsv', 'w') as fh:
        fh.write('variant\tculture_holdout_kept_of_751\tspecies_holdout_false_known_of_71\t'
                 'lsu_removal_false_known_of_201\tlsu_removal_34_test_of_78\toriginal_firings_kept\n')
        for row in table:
            fh.write('\t'.join(map(str, row)) + '\n')
    print(f'done in {time.time() - t0:.0f}s')
