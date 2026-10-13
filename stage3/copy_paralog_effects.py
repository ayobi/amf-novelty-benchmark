#!/usr/bin/env python3
"""What the short-5.8S paralog does to species counts and to trees (exploratory, no decision rule): L. Tedersoo's
follow-up to II.1 (13 October 2026). He notes that if short copies of different species group together, a species
carries more than one structural variant of the marker, which would raise species numbers counted from the ITS and
give conflicting signals between SSU, ITS and LSU trees.

  counts   Single-linkage clusters of all retained copies of the 37 named species, within each order, on the full
           ITS, the Kruger fragment and the whole copy (as in copy_long_fragments.py), at 98.5% and 97% identity;
           once with every copy and once with main-type copies only. Reports clusters against species, and the
           clusters made only of divergent copies.
  trees    Where the short-5.8S copies sit in the Glomerales tree of each marker (SSU flank, full ITS, LSU flank,
           linked copy; benchmarks/copy_checks/trees), rooted on a Funneliformis tip: clades made only of short
           copies, the taxa they join, their support, and short copies whose sequence is identical to a main-type
           copy (one tip) in that marker.

Short copies are the short-5.8S divergent class plus any main-type copy whose 5.8S is identical to a short-class 5.8S:
one copy, PX214776.1 (R. irregularis 4569_1185B), sits in a culture of two copies, where the 5.8S length rule
(3 nt or more below the culture's main-type median) cannot flag it.

  PYTHONPATH=<benchmark repo>/stage3 python3 stage3/copy_paralog_effects.py --ref reference/culture_ref_v4 \\
      --trees benchmarks/copy_checks/trees --out benchmarks/copy_checks/paralog_effects

Writes paralog.counts.tsv, paralog.species.tsv, paralog.trees.tsv and paralog.tree_clades.tsv.
"""
import argparse
import collections
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arms import is_named, read_fasta, read_tsv, write_tsv  # noqa: E402
from copy_long_fragments import load  # noqa: E402
from copy_thresholds import dist  # noqa: E402

REGIONS = [('full_its', 'Full ITS'), ('kruger', 'Kruger fragment'), ('full', 'Whole copy')]
CUTS = [(1.5, '98.5%'), (3.0, '97%')]
TREES = ['ssu_flank', 'full_its', 'lsu_flank', 'linked']


def find(par, x):
    while par[x] != x:
        par[x] = par[par[x]]
        x = par[x]
    return x


def cmd_counts(a, rows, seq, short):
    out, per_sp = [], []
    by_order = collections.defaultdict(list)
    for r in rows:
        by_order[r['order']].append(r)
    for m, lab in REGIONS:
        cache = {}
        pairs = []                                            # (distance, acc1, acc2) within an order
        for rs in by_order.values():
            for i, r1 in enumerate(rs):
                for r2 in rs[i + 1:]:
                    pairs.append((dist(seq[m][r1['accession']], seq[m][r2['accession']], cache),
                                  r1['accession'], r2['accession']))
        info = {r['accession']: r for r in rows}
        for t, tl in CUTS:
            for cset in ('all copies', 'main-type copies'):
                keep = {x for x, r in info.items() if cset == 'all copies' or (r['copy_class'] == 'main'
                                                                                 and x not in short)}
                par = {x: x for x in keep}
                for d, x, y in pairs:
                    if d <= t and x in keep and y in keep:
                        par[find(par, x)] = find(par, y)
                cl = collections.defaultdict(list)
                for x in keep:
                    cl[find(par, x)].append(info[x])
                cls = list(cl.values())
                sp_cl = collections.defaultdict(set)
                for k, c in enumerate(cls):
                    for r in c:
                        sp_cl[r['organism']].add(k)
                kinds = collections.Counter()
                for c in cls:
                    dk = {'short_5.8S' if r['accession'] in short else
                          (r['div_kind'] if r['copy_class'] == 'divergent' else 'main') for r in c}
                    if dk == {'short_5.8S'}:
                        kinds['short-5.8S only'] += 1
                    elif dk == {'same_length'}:
                        kinds['same-length only'] += 1
                    elif 'main' not in dk:
                        kinds['divergent only, mixed'] += 1
                out.append({'region': m, 'identity': tl, 'copies': cset, 'n_copies': len(keep),
                            'species': len(sp_cl), 'clusters': len(cls),
                            'clusters_per_species': round(len(cls) / len(sp_cl), 2),
                            'species_in_one_cluster': sum(len(v) == 1 for v in sp_cl.values()),
                            'species_in_2plus_clusters': sum(len(v) > 1 for v in sp_cl.values()),
                            'max_clusters_one_species': max(len(v) for v in sp_cl.values()),
                            'clusters_joining_species': sum(len({r['organism'] for r in c}) > 1 for c in cls),
                            'clusters_short_only': kinds['short-5.8S only'],
                            'clusters_same_length_only': kinds['same-length only'],
                            'clusters_divergent_mixed_only': kinds['divergent only, mixed']})
                if cset == 'all copies':
                    for sp, ks in sorted(sp_cl.items()):
                        rs = [r for r in rows if r['organism'] == sp]
                        nshort = sum(r['accession'] in short for r in rs)
                        if not nshort:
                            continue
                        own = [k for k in ks if all(r['accession'] in short for r in cls[k] if r['organism'] == sp)
                               and any(r['organism'] == sp for r in cls[k])]
                        short_only = [k for k in own if all(r['accession'] in short for r in cls[k])]
                        per_sp.append({'region': m, 'identity': tl, 'species': sp, 'copies': len(rs),
                                       'short_copies': nshort, 'clusters': len(ks),
                                       'clusters_holding_only_its_short_copies': len(own),
                                       'of_which_short_only_clusters': len(short_only),
                                       'short_only_clusters_shared_with_other_species':
                                           sum(len({r['organism'] for r in cls[k]}) > 1 for k in short_only)})
    write_tsv(os.path.join(a.out, 'paralog.counts.tsv'), out, list(out[0]))
    write_tsv(os.path.join(a.out, 'paralog.species.tsv'), per_sp, list(per_sp[0]))
    for r in out:
        print(f"{r['region']:9s} {r['identity']:6s} {r['copies']:17s} {r['clusters']:3d} clusters for {r['species']} "
              f"species ({r['clusters_per_species']}/species; split {r['species_in_2plus_clusters']}, "
              f"joining {r['clusters_joining_species']}); short-only {r['clusters_short_only']}, "
              f"same-length-only {r['clusters_same_length_only']}")


def cmd_trees(a, ref, short):
    from Bio import Phylo
    out, clades = [], []
    for t in TREES:
        src = os.path.join(a.trees, 'Glomerales', t)
        mem = {r['tip']: r['accessions'].split(',') for r in read_tsv(src + '.members.tsv')}
        tree = Phylo.read(src + '.treefile', 'newick')
        info = lambda n: [ref[x] for x in mem[n]]  # noqa: E731
        outg = next(tp for tp in tree.get_terminals() if all(r['genus'] == 'Funneliformis' for r in info(tp.name)))
        tree.root_with_outgroup(outg)
        short_tip = lambda n: all(x in short for x in mem[n])  # noqa: E731
        shorts = [x for n in mem for x in mem[n] if x in short]
        shared = [x for n in mem if not short_tip(n) for x in mem[n] if x in short]
        parents = {c: p for p in tree.find_clades(order='level') for c in p.clades}
        pure = []

        def walk(cl):
            if all(short_tip(tp.name) for tp in cl.get_terminals()):
                pure.append(cl)
                return
            for ch in cl.clades:
                walk(ch)
        walk(tree.root)
        n_multi = n_cross = n_cross_clades = 0
        sister_own = sister_other = 0
        for cl in pure:
            rs = [r for tp in cl.get_terminals() for r in info(tp.name)]
            taxa = sorted({r['organism'] for r in rs})
            lab = cl.name or (str(cl.confidence) if cl.confidence is not None else '')
            mm = re.match(r'([\d.]+)/([\d.]+)', lab or '')
            # the main-type copies of the sister group (the parent clade minus this one)
            par = parents.get(cl)
            sis = [r for ch in (par.clades if par else []) if ch is not cl for tp in ch.get_terminals()
                   for r in info(tp.name) if r['copy_class'] == 'main' and r['accession'] not in short]
            sis_taxa = sorted({r['organism'] for r in sis})
            own = bool(sis_taxa) and set(sis_taxa) <= set(taxa)
            sister_own += len(rs) if own else 0
            sister_other += len(rs) if not own else 0
            if len(rs) >= 2:
                n_multi += len(rs)
                if len(taxa) > 1:
                    n_cross += len(rs)
                    n_cross_clades += 1
            plab = (par.name or (str(par.confidence) if par.confidence is not None else '')) if par else ''
            pm = re.match(r'([\d.]+)/([\d.]+)', plab or '')
            clades.append({'tree': t, 'copies': len(rs), 'tips': len(cl.get_terminals()), 'taxa': '; '.join(taxa),
                           'sh_alrt': mm.group(1) if mm else '', 'ufboot': mm.group(2) if mm else '',
                           'sister_main_taxa': '; '.join(sis_taxa),
                           'with_sister_sh_alrt': pm.group(1) if pm else '', 'with_sister_ufboot': pm.group(2) if pm else '',
                           'sister_is_own_species_only': 'yes' if own else 'no'})
        out.append({'tree': t, 'short_copies': len(shorts), 'short_in_tips_with_main_copies': len(shared),
                    'short_only_clades': len(pure), 'short_only_clades_2plus': sum(
                        1 for cl in pure if sum(len(mem[tp.name]) for tp in cl.get_terminals()) >= 2),
                    'copies_in_short_only_clades_2plus': n_multi, 'short_only_clades_joining_taxa': n_cross_clades,
                    'copies_in_clades_joining_taxa': n_cross,
                    'short_copies_whose_sister_is_own_species_only': sister_own,
                    'short_copies_whose_sister_has_other_species': sister_other})
    write_tsv(os.path.join(a.out, 'paralog.trees.tsv'), out, list(out[0]))
    write_tsv(os.path.join(a.out, 'paralog.tree_clades.tsv'), clades, list(clades[0]))
    for r in out:
        print(r)
    for c in clades:
        if c['copies'] >= 2:
            print(c)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--trees', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--only', choices=['counts', 'trees'])
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cul = {r['culture']: r for r in read_tsv(os.path.join(a.ref, 'cultures.tsv'))}
    named = {c for c, r in cul.items() if is_named(r['organism']) and r.get('name_status', 'ok') == 'ok'}
    allr = [r for r in read_tsv(os.path.join(a.ref, 'copies.tsv')) if r['status'] == 'retained']
    ref = {r['accession']: r for r in allr}
    s58 = read_fasta(os.path.join(a.ref, 'fasta', 's58.fasta'))
    cls58 = {s58[x] for x, r in ref.items() if r['div_kind'] == 'short_5.8S'}
    short = {x for x in ref if ref[x]['div_kind'] == 'short_5.8S' or s58.get(x) in cls58}
    print(f'{len(short)} short-5.8S copies ({len(short) - sum(r["div_kind"] == "short_5.8S" for r in allr)} '
          f'classed main type: {sorted(x for x in short if ref[x]["div_kind"] != "short_5.8S")})')
    if a.only != 'counts':
        cmd_trees(a, ref, short)
    if a.only != 'trees':
        seq, _ = load(a.ref)
        rows = [r for r in allr if r['culture'] in named and all(r['accession'] in seq[m] for m, _ in REGIONS)]
        cmd_counts(a, rows, seq, short)


if __name__ == '__main__':
    main()
