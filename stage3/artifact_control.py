#!/usr/bin/env python3
"""Positive and negative controls for the artifact screen (step 1 of the novelty steps), added after review.

The screen flags a copy when its SSU flank differs by >= 1% from the nearest main-type copy of another culture
and, for every SSU-tied nearest copy, by at least half as much as the full ITS (frozen decision_tree.py,
artifact_min_ssu = 1, artifact_ratio = 0.5). The paper showed it fires for no real copy; this script asks what it
does fire for. Exploratory and written after the results; it changes no rule.

Queries, each scored with its parents' culture labels withheld (as in the culture-withheld test):
  real            every main-type copy, unchanged (negative control)
  errors_<r>      one main-type copy per culture label, with random errors at rate r (% of positions; 60%
                  substitutions, 20% insertions, 20% deletions), applied to the SSU flank and the full ITS alike
  chim_junction   SSU flank of copy A + full ITS of copy B (breakpoint between SSU and ITS)
  chim_ssu        SSU flank of A up to a random point and of B after it, + full ITS of B (breakpoint in the SSU)
Chimera parents are main-type copies of named species, from different genera or from congeneric species.
Distances are the frozen dist() (edlib infix, % of the shorter sequence). Seeded; deterministic.

  python3 stage3/artifact_control.py --ref reference/culture_ref_v4 --out benchmarks/copy_checks/artifact_control
"""
import argparse
import collections
import csv
import math
import os
import random
import sys
import time
from multiprocessing import Pool

import edlib

RATES = (0.5, 1, 2, 3, 5, 8)
N_CHIM = 150
MIN_SSU, RATIO = 1.0, 0.5


def read_fasta(path):
    out, k = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            k = line[1:].split()[0]
            out[k] = []
        elif k:
            out[k].append(line.upper())
    return {k: ''.join(v) for k, v in out.items()}


def dist(a, b, k=-1):                     # frozen novelty_baseline_consistent.dist
    if not a or not b:
        raise ValueError('Cannot compare empty marker sequences')
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    d = edlib.align(q, t, mode='HW', task='distance', k=k)['editDistance']
    return None if d < 0 else 100 * d / len(q)


G = {}


def init(state):
    G.update(state)


def screen(job):
    """Frozen A1 screen: SSU-nearest training main copies (all ties), then the ratio test against each."""
    name, q_ssu, q_its, drop, genus_hint = job
    ssu, its, mains, cop = G['ssu'], G['its'], G['mains'], G['cop']
    cand = [y for y in mains if cop[y]['culture_label'] not in drop]
    cand.sort(key=lambda y: (cop[y]['genus'] != genus_hint, y))          # order only; the result does not depend on it
    near_d, nearest = None, []
    for y in cand:
        k = -1 if near_d is None else int(near_d * min(len(q_ssu), len(ssu[y])) / 100) + 1
        d = dist(q_ssu, ssu[y], k)
        if d is None:
            continue
        if near_d is None or d < near_d - 1e-12:
            near_d, nearest = d, [y]
        elif math.isclose(d, near_d, rel_tol=0, abs_tol=1e-12):
            nearest.append(y)
    flag, its_d = False, []
    if near_d is not None and near_d >= MIN_SSU:
        its_d = [dist(q_its, its[y]) for y in nearest if y in its]
        votes = [near_d >= RATIO * d for d in its_d]
        flag = len(votes) == len(nearest) and bool(votes) and all(votes)
    return {'query': name, 'ssu_pct': round(near_d, 3) if near_d is not None else '',
            'its_pct': ';'.join(f'{d:.2f}' for d in its_d), 'flag': int(flag)}


def mutate(seq, rate, rng):
    out = []
    for ch in seq:
        if rng.random() < rate:
            u = rng.random()
            if u < 0.6:
                out.append(rng.choice([b for b in 'ACGT' if b != ch]))
            elif u < 0.8:
                out.append(ch + rng.choice('ACGT'))
            # else: deletion
        else:
            out.append(ch)
    return ''.join(out)


def wilson(k, n, z=1.96):
    if n == 0:
        return ('', '')
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (round(100 * max(0, c - h), 1), round(100 * min(1, c + h), 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=20261015)
    ap.add_argument('--procs', type=int, default=os.cpu_count())
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(a.out, exist_ok=True)
    cop = {r['accession']: r for r in csv.DictReader(open(os.path.join(a.ref, 'copies.tsv')), delimiter='\t')
           if r['status'] == 'retained'}
    ssu = read_fasta(os.path.join(a.ref, 'fasta', 'ssu_flank.fasta'))
    its = read_fasta(os.path.join(a.ref, 'fasta', 'full_its.fasta'))
    mains = sorted(x for x, r in cop.items() if r['copy_class'] == 'main' and x in ssu and x in its)
    named = [x for x in mains if ' sp.' not in cop[x]['organism'] and not cop[x]['organism'].endswith(' sp')]
    rng = random.Random(a.seed)
    jobs, kinds = [], {}

    def add(name, kind, q_ssu, q_its, drop, genus):
        jobs.append((name, q_ssu, q_its, frozenset(drop), genus))
        kinds[name] = kind

    for x in mains:                                                   # negative control
        add(f'real|{x}', 'real', ssu[x], its[x], {cop[x]['culture_label']}, cop[x]['genus'])
    by_label = collections.defaultdict(list)
    for x in mains:
        by_label[cop[x]['culture_label']].append(x)
    picks = [rng.choice(sorted(v)) for _, v in sorted(by_label.items())]
    for rate in RATES:                                                # random errors
        for x in picks:
            add(f'errors_{rate}|{x}', f'errors_{rate}', mutate(ssu[x], rate / 100, rng), mutate(its[x], rate / 100, rng),
                {cop[x]['culture_label']}, cop[x]['genus'])
    for kind in ('chim_junction', 'chim_ssu'):                        # chimeras
        for rel in ('other genus', 'same genus'):
            n = 0
            while n < N_CHIM:
                A, B = rng.choice(named), rng.choice(named)
                ga, gb = cop[A]['genus'], cop[B]['genus']
                if cop[A]['organism'] == cop[B]['organism'] or (rel == 'other genus') != (ga != gb):
                    continue
                if kind == 'chim_junction':
                    q_ssu = ssu[A]
                else:
                    f = rng.uniform(0.1, 0.9)
                    q_ssu = ssu[A][:int(f * len(ssu[A]))] + ssu[B][int(f * len(ssu[B])):]
                add(f'{kind}|{rel}|{A}+{B}|{n}', f'{kind} ({rel})', q_ssu, its[B],
                    {cop[A]['culture_label'], cop[B]['culture_label']}, ga)
                n += 1
    print(f'{len(jobs)} queries; {len(mains)} main-type copies as candidates; {a.procs} processes', flush=True)
    with Pool(a.procs, initializer=init, initargs=({'ssu': ssu, 'its': its, 'mains': mains, 'cop': cop},)) as pool:
        res = pool.map(screen, jobs, chunksize=8)
    rows = [dict(r, kind=kinds[r['query']]) for r in res]
    with open(os.path.join(a.out, 'artifact_control.queries.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['kind', 'query', 'ssu_pct', 'its_pct', 'flag'], delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    summ = []
    for kind in dict.fromkeys(kinds.values()):
        rs = [r for r in rows if r['kind'] == kind]
        k = sum(r['flag'] for r in rs)
        lo, hi = wilson(k, len(rs))
        ssu_med = sorted(r['ssu_pct'] for r in rs if r['ssu_pct'] != '')
        summ.append({'kind': kind, 'queries': len(rs), 'flagged': k, 'flagged_pct': round(100 * k / len(rs), 1),
                     'ci95_low': lo, 'ci95_high': hi, 'median_ssu_pct': ssu_med[len(ssu_med) // 2] if ssu_med else ''})
    with open(os.path.join(a.out, 'artifact_control.summary.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(summ)
    for s in summ:
        print(f"  {s['kind']:28s} {s['flagged']:4d} of {s['queries']:4d} flagged ({s['flagged_pct']}%, 95% CI "
              f"{s['ci95_low']}-{s['ci95_high']}); median SSU distance to nearest {s['median_ssu_pct']}%")
    print(f'({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    sys.exit(main())
