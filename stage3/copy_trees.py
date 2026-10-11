#!/usr/bin/env python3
"""Per-marker trees of the culture copies: are cultures and species monophyletic?

PROTOCOL_COPY_CHECKS.md, section 1. One tree set per order (ITS cannot be aligned across orders) on the
SSU flank, the full ITS, the LSU flank and the linked copy (the three concatenated, one partition each).
Identical sequences are collapsed before alignment; each tree tip stands for every copy that carries its
sequence. Monophyly is read from the bipartitions of the unrooted tree, so no rooting is needed.

Subcommands, in run order (run_copy_trees.sh calls them between MAFFT and IQ-TREE):
  prep    write {out}/{order}/{marker}.unique.fasta and {marker}.members.tsv for the three markers
  concat  after MAFFT: build {order}/linked.aln.fasta, linked.members.tsv and linked.partitions.nex
  score   after IQ-TREE: test every culture, named species and complex on every tree, with all copies
          and with main-type copies only; write monophyly.groups.tsv, monophyly.summary.tsv, run.json

  python3 stage3/copy_trees.py prep   --ref reference/culture_ref_v4 --out benchmarks/copy_checks/trees
  python3 stage3/copy_trees.py concat --out benchmarks/copy_checks/trees
  python3 stage3/copy_trees.py score  --ref reference/culture_ref_v4 --out benchmarks/copy_checks/trees \\
      --complexes paper/analysis/v4.complexes.tsv

Outcome per group and tree:
  monophyletic      a branch separates exactly the group's copies from all others
  unresolved        no branch does, and none contradicts the group
  non-monophyletic  a branch contradicts it; 'strong' when a contradicting branch has SH-aLRT >= 80 and
                    UFBoot >= 95 (a branch that only joins copies with the same sequence has no support)
escapees: the group's copies outside the largest clade made only of its copies.
intruders: other copies inside the smallest clade holding all of the group's copies.
"""
import argparse
import collections
import hashlib
import json
import os
import re
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arms import is_named, read_fasta, read_tsv, write_tsv  # noqa: E402

MARKERS = ['ssu_flank', 'full_its', 'lsu_flank']
TREES = MARKERS + ['linked']
ORDERS = ['Archaeosporales', 'Paraglomerales', 'Diversisporales', 'Glomerales']
ALRT, UFB = 80.0, 95.0          # strong support (PROTOCOL_COPY_CHECKS.md, section 1)
SEED = 20261011


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def retained(ref):
    return [r for r in read_tsv(os.path.join(ref, 'copies.tsv')) if r['status'] == 'retained']


def write_members(path, groups):
    """groups: {tip: [accessions]}"""
    with open(path, 'w') as fh:
        fh.write('tip\tn\taccessions\n')
        for tip in sorted(groups):
            fh.write(f'{tip}\t{len(groups[tip])}\t{",".join(groups[tip])}\n')


def read_members(path):
    return {r['tip']: r['accessions'].split(',') for r in read_tsv(path)}


# ------------------------------------------------------------------ prep
def cmd_prep(a):
    rows = retained(a.ref)
    seqs = {m: read_fasta(os.path.join(a.ref, 'fasta', f'{m}.fasta')) for m in MARKERS}
    summary = []
    for order in ORDERS:
        accs = sorted(r['accession'] for r in rows if r['order'] == order)
        os.makedirs(os.path.join(a.out, order), exist_ok=True)
        for m in MARKERS:
            by_seq = collections.defaultdict(list)
            for acc in accs:
                by_seq[seqs[m][acc]].append(acc)
            groups = {members[0]: members for members in by_seq.values()}
            with open(os.path.join(a.out, order, f'{m}.unique.fasta'), 'w') as fh:
                for tip in sorted(groups):
                    fh.write(f'>{tip}\n{seqs[m][tip]}\n')
            write_members(os.path.join(a.out, order, f'{m}.members.tsv'), groups)
            summary.append(f'{order} {m}: {len(accs)} copies, {len(groups)} unique')
    print('\n'.join(summary))


# ------------------------------------------------------------------ concat
def cmd_concat(a):
    for order in ORDERS:
        d = os.path.join(a.out, order)
        aln, mem = {}, {}
        for m in MARKERS:
            aln[m] = read_fasta(os.path.join(d, f'{m}.aln.fasta'))
            mem[m] = read_members(os.path.join(d, f'{m}.members.tsv'))
            if set(aln[m]) != set(mem[m]):
                sys.exit(f'{d}/{m}.aln.fasta does not hold the tips of {m}.members.tsv')
        tip_of = {m: {acc: tip for tip, accs in mem[m].items() for acc in accs} for m in MARKERS}
        accs = sorted(tip_of[MARKERS[0]])
        width = {m: len(next(iter(aln[m].values()))) for m in MARKERS}
        by_row = collections.defaultdict(list)
        for acc in accs:
            by_row[''.join(aln[m][tip_of[m][acc]] for m in MARKERS)].append(acc)
        groups = {members[0]: members for members in by_row.values()}
        row_of = {members[0]: row for row, members in by_row.items()}
        with open(os.path.join(d, 'linked.aln.fasta'), 'w') as fh:
            for tip in sorted(groups):
                fh.write(f'>{tip}\n{row_of[tip]}\n')
        write_members(os.path.join(d, 'linked.members.tsv'), groups)
        start, parts = 1, []
        for m in MARKERS:
            parts.append(f'  charset {m} = {start}-{start + width[m] - 1};')
            start += width[m]
        with open(os.path.join(d, 'linked.partitions.nex'), 'w') as fh:
            fh.write('#nexus\nbegin sets;\n' + '\n'.join(parts) + '\nend;\n')
        print(f'{order} linked: {len(accs)} copies, {len(groups)} unique rows, '
              + ', '.join(f'{m} {width[m]}' for m in MARKERS) + ' columns')


# ------------------------------------------------------------------ Newick
def parse_newick(text):
    """Minimal Newick parser. Returns nested nodes: {'name', 'label', 'children'}; tips have no children.
    Internal labels are kept as written (IQ-TREE: 'SH-aLRT/UFBoot')."""
    text = text.strip().rstrip(';')
    i, n = 0, len(text)

    def read_token():
        nonlocal i
        if i < n and text[i] == "'":
            j = text.index("'", i + 1)
            tok = text[i + 1:j]
            i = j + 1
            return tok
        j = i
        while j < n and text[j] not in ',():;':
            j += 1
        tok = text[i:j]
        i = j
        return tok.strip()

    def node():
        nonlocal i
        nd = {'children': [], 'name': '', 'label': ''}
        if text[i] == '(':
            i += 1
            while True:
                nd['children'].append(node())
                if text[i] == ',':
                    i += 1
                    continue
                if text[i] == ')':
                    i += 1
                    break
                raise ValueError(f'unexpected {text[i]!r} at {i}')
            tok = read_token()
            nd['label'] = tok.split(':')[0]
        else:
            tok = read_token()
            nd['name'] = tok.split(':')[0]
        if i < n and text[i] == ':':
            i += 1
            while i < n and text[i] not in ',();':
                i += 1
        return nd

    root = node()
    return root


def support(label):
    """'87.5/99' -> (87.5, 99.0); anything else -> (None, None)."""
    parts = label.split('/')
    try:
        vals = [float(p) for p in parts]
    except ValueError:
        return None, None
    if len(vals) == 2:
        return vals[0], vals[1]
    return None, None


def tree_splits(root, members, index):
    """Bipartitions over copies. Returns (full mask, {mask: (supported, aLRT, UFBoot)}, tip masks).
    Each internal edge gives a split; each tip edge gives the split {its copies} | rest (identity, no support).
    Masks are canonical: the side without the lowest-indexed copy."""
    tips = {}
    raw = []   # (mask, supported, alrt, ufb)

    def walk(nd, is_root):
        if not nd['children']:
            if nd['name'] not in members:
                sys.exit(f'tree tip {nd["name"]} not in the member map')
            m = 0
            for acc in members[nd['name']]:
                m |= 1 << index[acc]
            tips[nd['name']] = m
            raw.append((m, False, None, None))
            return m
        m = 0
        for c in nd['children']:
            m |= walk(c, False)
        if not is_root:
            al, uf = support(nd['label'])
            ok = al is not None and al >= ALRT and uf >= UFB
            raw.append((m, ok, al, uf))
        return m

    full = walk(root, True)
    return full, raw, tips


def restrict(raw, keep, full):
    """Restrict raw splits to the copies in keep; drop trivial ones; merge duplicates (support: any edge)."""
    full_k = full & keep
    low = full_k & -full_k
    out = {}
    for m, ok, al, uf in raw:
        a = m & keep
        b = full_k & ~a
        if bin(a).count('1') < 2 or bin(b).count('1') < 2:
            continue
        c = b if a & low else a
        prev = out.get(c)
        if prev is None:
            out[c] = (ok, al, uf)
        else:
            # keep the best-supported edge that implies this bipartition
            best = prev
            if (ok and not prev[0]) or (ok == prev[0] and (uf or -1) > (prev[2] or -1)):
                best = (ok, al, uf)
            out[c] = best
    return full_k, out


def popcount(x):
    return bin(x).count('1')


def test_group(g, full, splits, tip_masks):
    """Outcome of group mask g on a restricted split set."""
    rest = full & ~g
    low = full & -full
    canon = rest if g & low else g
    res = {'status': '', 'strong': False, 'support': '', 'n_contra': 0, 'n_contra_strong': 0,
           'identical_shared': any(t & g and t & rest for t in tip_masks), 'esc': 0, 'intr': 0}
    if canon in splits:
        ok, al, uf = splits[canon]
        res['status'] = 'monophyletic'
        res['support'] = f'{al:g}/{uf:g}' if al is not None else 'identity'
        return res, 0, 0
    best_pure, best_cover, best_cover_mask, best_pure_mask = 1, None, None, 0
    for a, (ok, al, uf) in splits.items():
        b = full & ~a
        if (a & g) and (a & rest) and (b & g) and (b & rest):
            res['n_contra'] += 1
            if ok:
                res['n_contra_strong'] += 1
        for side in (a, b):
            if side & rest == 0:            # pure side: only group copies
                k = popcount(side)
                if k > best_pure:
                    best_pure, best_pure_mask = k, side
            elif side & g == g:              # side covering the whole group
                k = popcount(side)
                if best_cover is None or k < best_cover:
                    best_cover, best_cover_mask = k, side
    if res['n_contra'] == 0:
        res['status'] = 'unresolved'
        return res, 0, 0
    res['status'] = 'non-monophyletic'
    res['strong'] = res['n_contra_strong'] > 0
    if best_pure_mask == 0:                  # largest pure clade is a single copy
        best_pure_mask = g & -g
    esc_mask = g & ~best_pure_mask
    intr_mask = (best_cover_mask & ~g) if best_cover_mask is not None else rest
    return res, esc_mask, intr_mask


def best_model(prefix):
    p = prefix + '.iqtree'
    if not os.path.exists(p):
        return ''
    txt = open(p).read()
    m = re.search(r'Best-fit model according to BIC:\s*(\S+)', txt)
    if m:
        return m.group(1)
    m = re.search(r'Model of substitution:\s*(\S+)', txt)
    return m.group(1) if m else 'partitioned'


def cmd_score(a):
    rows = retained(a.ref)
    by_acc = {r['accession']: r for r in rows}
    cul = {r['culture']: r for r in read_tsv(os.path.join(a.ref, 'cultures.tsv'))}
    named = {c for c, r in cul.items() if is_named(r['organism']) and r.get('name_status', 'ok') == 'ok'}
    complexes = {r['complex']: r['species'].split(';') for r in read_tsv(a.complexes)} if a.complexes else {}
    n_cult = collections.Counter(r['organism'] for c, r in cul.items() if c in named)

    groups_out, trees_info = [], {}
    for order in ORDERS:
        accs = sorted(r['accession'] for r in rows if r['order'] == order)
        index = {acc: i for i, acc in enumerate(accs)}
        main_mask = 0
        for acc in accs:
            if by_acc[acc]['copy_class'] == 'main':
                main_mask |= 1 << index[acc]
        defs = []   # (level, name, mask)
        by = collections.defaultdict(int)
        for acc in accs:
            r = by_acc[acc]
            by[('culture', r['culture'])] |= 1 << index[acc]
            if r['culture'] in named:
                by[('species', r['organism'])] |= 1 << index[acc]
        for cx, sps in complexes.items():
            m = 0
            for acc in accs:
                if by_acc[acc]['organism'] in sps and by_acc[acc]['culture'] in named:
                    m |= 1 << index[acc]
            if m:
                by[('complex', f'{cx}: ' + ' + '.join(sps))] = m
        defs = sorted((lv, nm, m) for (lv, nm), m in by.items())
        for t in TREES:
            prefix = os.path.join(a.out, order, t)
            tf = prefix + '.treefile'
            log = prefix + '.log'
            if not (os.path.exists(tf) and os.path.exists(log)
                    and any(x.startswith('Total wall-clock time used') for x in open(log))):
                print(f'{prefix}: no finished IQ-TREE run; skipped')
                continue
            members = read_members(os.path.join(a.out, order, f'{t}.members.tsv'))
            root = parse_newick(open(tf).read())
            full, raw, tips = tree_splits(root, members, index)
            if full != (1 << len(accs)) - 1:
                sys.exit(f'{tf}: tree does not hold every copy of {order}')
            trees_info[f'{order}/{t}'] = {'tips': len(tips), 'copies': len(accs), 'model': best_model(prefix),
                                         'treefile_sha256': sha256(tf)}
            for cset, keep in (('all', full), ('main', main_mask)):
                fk, splits = restrict(raw, keep, full)
                tmask = [m & keep for m in tips.values()]
                for lv, nm, gm in defs:
                    g = gm & keep
                    if popcount(g) < 2 or g == fk:
                        continue
                    res, esc, intr = test_group(g, fk, splits, tmask)
                    g_accs = [acc for acc in accs if g >> index[acc] & 1]
                    esc_accs = [acc for acc in accs if esc >> index[acc] & 1]
                    intr_accs = [acc for acc in accs if intr >> index[acc] & 1]
                    sp = nm if lv == 'species' else ''
                    groups_out.append({
                        'order': order, 'tree': t, 'copy_set': cset, 'level': lv, 'group': nm,
                        'n_copies': len(g_accs), 'n_cultures': len({by_acc[x]['culture'] for x in g_accs}),
                        'species_cultures': n_cult.get(sp, '') if sp else '',
                        'n_divergent': sum(by_acc[x]['copy_class'] == 'divergent' for x in g_accs),
                        'status': res['status'], 'strong': 'yes' if res['strong'] else '',
                        'split_support': res['support'], 'contradicting_branches': res['n_contra'],
                        'contradicting_strong': res['n_contra_strong'],
                        'identical_shared': 'yes' if res['identical_shared'] else '',
                        'escapees_n': len(esc_accs),
                        'escapees_divergent_n': sum(by_acc[x]['copy_class'] == 'divergent' for x in esc_accs),
                        'escapees': ','.join(f"{x}:{by_acc[x]['copy_class']}" for x in esc_accs),
                        'intruders_n': len(intr_accs),
                        'intruders': ','.join(f"{x}:{by_acc[x]['organism']}|{by_acc[x]['culture_label']}:"
                                              f"{by_acc[x]['copy_class']}" for x in intr_accs),
                    })
    cols = ['order', 'tree', 'copy_set', 'level', 'group', 'n_copies', 'n_cultures', 'species_cultures',
            'n_divergent', 'status', 'strong', 'split_support', 'contradicting_branches', 'contradicting_strong',
            'identical_shared', 'escapees_n', 'escapees_divergent_n', 'escapees', 'intruders_n', 'intruders']
    write_tsv(os.path.join(a.out, 'monophyly.groups.tsv'), groups_out, cols)

    # summary: per tree type, copy set and level; pooled over orders ('all') and per order
    summ = []
    levels = [('culture', lambda r: r['level'] == 'culture'),
              ('species', lambda r: r['level'] == 'species'),
              ('species, >= 2 cultures', lambda r: r['level'] == 'species' and int(r['species_cultures']) >= 2),
              ('complex', lambda r: r['level'] == 'complex')]
    for t in TREES:
        for cset in ('all', 'main'):
            for lv, f in levels:
                for order in ['all'] + ORDERS:
                    xs = [r for r in groups_out if r['tree'] == t and r['copy_set'] == cset and f(r)
                          and (order == 'all' or r['order'] == order)]
                    if not xs:
                        continue
                    c = collections.Counter(r['status'] for r in xs)
                    strong = sum(r['strong'] == 'yes' for r in xs)
                    summ.append({'tree': t, 'copy_set': cset, 'level': lv, 'order': order, 'tested': len(xs),
                                 'monophyletic': c['monophyletic'], 'unresolved': c['unresolved'],
                                 'non_monophyletic': c['non-monophyletic'], 'strongly_non_monophyletic': strong,
                                 'pct_monophyletic': f"{100 * c['monophyletic'] / len(xs):.1f}",
                                 'pct_not_strongly_rejected': f"{100 * (len(xs) - strong) / len(xs):.1f}"})
    scols = ['tree', 'copy_set', 'level', 'order', 'tested', 'monophyletic', 'unresolved', 'non_monophyletic',
             'strongly_non_monophyletic', 'pct_monophyletic', 'pct_not_strongly_rejected']
    write_tsv(os.path.join(a.out, 'monophyly.summary.tsv'), summ, scols)

    versions = {}
    vf = os.path.join(a.out, 'versions.txt')
    if os.path.exists(vf):
        versions = dict(line.rstrip('\n').split('\t', 1) for line in open(vf) if '\t' in line)
    run = {'finished': time.strftime('%Y-%m-%d %H:%M:%S'), 'reference': os.path.abspath(a.ref),
           'reference_manifest_sha256': sha256(os.path.join(a.ref, 'MANIFEST.json'))
           if os.path.exists(os.path.join(a.ref, 'MANIFEST.json')) else '',
           'copies_tsv_sha256': sha256(os.path.join(a.ref, 'copies.tsv')),
           'complexes_sha256': sha256(a.complexes) if a.complexes else '',
           'strong_support': {'SH-aLRT': ALRT, 'UFBoot': UFB}, 'versions': versions, 'trees': trees_info,
           'outputs': {f: sha256(os.path.join(a.out, f)) for f in ('monophyly.groups.tsv', 'monophyly.summary.tsv')}}
    with open(os.path.join(a.out, 'run.json'), 'w') as fh:
        json.dump(run, fh, indent=1)

    # console: the headline table
    print('tree       set   level                      tested  mono  unres  non  strong  %mono')
    for r in summ:
        if r['order'] == 'all':
            print(f"{r['tree']:<10} {r['copy_set']:<5} {r['level']:<26} {r['tested']:>6} {r['monophyletic']:>5} "
                  f"{r['unresolved']:>6} {r['non_monophyletic']:>4} {r['strongly_non_monophyletic']:>7} "
                  f"{r['pct_monophyletic']:>6}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prep')
    p.add_argument('--ref', required=True)
    p.add_argument('--out', required=True)
    p = sub.add_parser('concat')
    p.add_argument('--out', required=True)
    p = sub.add_parser('score')
    p.add_argument('--ref', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--complexes', default='')
    a = ap.parse_args()
    {'prep': cmd_prep, 'concat': cmd_concat, 'score': cmd_score}[a.cmd](a)


if __name__ == '__main__':
    main()
