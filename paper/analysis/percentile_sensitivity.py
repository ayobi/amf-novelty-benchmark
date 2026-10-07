#!/usr/bin/env python3
"""Move the ITS rule and the species complexes to the same percentile, and see what changes.

The ITS rule (N6) and the species complexes both use the 99th percentile of within-species distance.
This script re-applies N6 at other percentiles to the stored decision-tree rows (no BLAST, placement
or tree rerun is needed: N6 only post-processes the rows, and each row keeps the copy's full-ITS
distance to the nearest training culture) and counts the remaining false-known calls against the
complexes that species_groups.py gives at the same percentile.

  M4 species folds (first-outcome folds, target removed from V18): stored v1 rows (before N6) plus the
    distance from the v2 rows; T per fold from the fold's training cultures.
  Culture-bridge benchmark: placement configuration before N6 (epa_species_clades) plus the distance
    from the +its rows; culture holdout (known-species copies demoted) and species holdout (false known).

T comes from the frozen ItsRule on the frozen data (exclusive quantiles, generalised to fractional
percentiles). At the 99th percentile every stored v2 / +its status must be reproduced exactly.

Usage:
  python3 percentile_sensitivity.py --m4root <amf_m4_step1> --m4 <aggregate_v2> --its-run <its_rule_rerun> \\
      --complexes 95=analysis/pct/v4_p95 97.5=analysis/pct/v4_p97.5 99=analysis/v4 --out analysis/pct/n6_percentile.json
"""
import argparse
import collections
import csv
import json
import statistics
import sys
from pathlib import Path


def rd(p):
    with open(p, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def q_excl(vals, p):
    """statistics.quantiles(method='exclusive') at any percentile p (identical at integer p)."""
    d = sorted(vals)
    n = len(d)
    h = (n + 1) * p / 100.0
    j = min(max(int(h), 1), n - 1)
    return d[j - 1] + (d[j] - d[j - 1]) * (h - j)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--m4root', type=Path, required=True, help='amf_m4_step1 (external_holdout.py and frozen_m2b)')
    ap.add_argument('--m4', type=Path, required=True)
    ap.add_argument('--its-run', type=Path, required=True)
    ap.add_argument('--complexes', nargs='+', required=True, help='percentile=prefix of species_groups.py output')
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()

    sys.path.insert(0, str(a.m4root))
    import external_holdout as v1       # sets up the frozen_m2b imports
    from its_rule import ItsRule

    class Rule(ItsRule):
        def threshold_p(self, training, p):
            data, vals = self.data, []
            training = frozenset(training)
            for sp, cs in data.by_species.items():
                tc = cs & training
                if len(tc) < 2:
                    continue
                for c in tc:
                    for y in data.by_culture[c]:
                        if y in self.mains:
                            ds = [self.dmin['main'][y][o] for o in tc - {c} if o in self.dmin['main'][y]]
                            if ds:
                                vals.append(min(ds))
            for ip in (95, 99):
                assert abs(q_excl(vals, ip) - statistics.quantiles(vals, n=100)[ip - 1]) < 1e-9
            return q_excl(vals, p)

    data = v1.Data(v1.FROZEN / 'inputs/benchmark')
    rule = Rule(data, 99)
    allc = set(data.cul)
    thr = {}

    def T(held, p):
        if (held, p) not in thr:
            thr[(held, p)] = rule.threshold_p(allc - held, p)
        return thr[(held, p)]

    def demote(status, d, t):
        return 'novel species, genus known' if d != '' and float(d) > t and status.startswith('known') else status

    cx = {}
    for item in a.complexes:
        p, prefix = item.split('=')
        m = {}
        for r in rd(f'{prefix}.complexes.tsv'):
            for s in r['species'].split(';'):
                m[s] = r['complex']
        cx[float(p)] = m
    res = {'percentiles': sorted(cx)}
    # M4 species folds, first outcome
    rows = [r for r in rd(a.m4 / 'species.v1.copies.tsv') if r['copy_class'] == 'main' and r['group'] == 'named, species held out'
            and r['development_fold'] == '0' and r['external_status'] == 'removed']
    v2 = {(r['target'], r['accession']): r for r in rd(a.m4 / 'species.v2.copies.tsv')}
    for p, c in sorted(cx.items()):
        fk, outside, ts = collections.Counter(), collections.Counter(), {}
        for r in rows:
            held = frozenset(r['held_cultures'].split(';'))
            ts[held] = T(held, p)
            st = demote(r['status'], v2[(r['target'], r['accession'])]['its_nearest_training_pct'], ts[held])
            if p == 99.0 and st != v2[(r['target'], r['accession'])]['status']:
                raise SystemExit(f'99th percentile does not reproduce the stored v2 status of {r["accession"]}')
            if st.startswith('known'):
                pair = f"{r['organism']} -> {r['placement']}"
                fk[pair] += 1
                sp, pl = r['organism'], [x for x in r['placement'].split(';') if x]
                if not (sp in c and any(c.get(x) == c[sp] for x in pl)):
                    outside[pair] += 1
        res[f'm4@{p}'] = {'held': len(rows), 'false_known': sum(fk.values()), 'outside_complex': sum(outside.values()),
                          'outside_pairs': dict(outside), 'pairs': dict(fk), 'complexes': len(set(c.values())),
                          'complex_species': len(c), 'median_T': statistics.median(ts.values())}
    # culture-bridge benchmark, placement configuration
    for mode, grp in (('culture', 'named, species elsewhere'), ('species', 'named, species held out')):
        pre = [r for r in rd(a.its_run / f'epa_species_clades_physical_{mode}.copies.tsv') if r['copy_class'] == 'main' and r['group'] == grp]
        post = {r['accession']: r for r in rd(a.its_run / f'epa_species_clades+its_physical_{mode}.copies.tsv')}
        for p in sorted(cx):
            n = 0
            for r in pre:
                held = frozenset(r['held_cultures'].split(';'))
                st = demote(r['status'], post[r['accession']]['its_nearest_training_pct'], T(held, p))
                if p == 99.0 and st != post[r['accession']]['status']:
                    raise SystemExit(f'99th percentile does not reproduce the stored +its status of {r["accession"]}')
                n += (r['status'].startswith('known') and not st.startswith('known')) if mode == 'culture' else st.startswith('known')
            res[f'{mode}@{p}'] = {'n': len(pre), ('known_demoted' if mode == 'culture' else 'false_known'): n}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + '\n')
    for k, v in res.items():
        print(k, v)


if __name__ == '__main__':
    sys.exit(main())
