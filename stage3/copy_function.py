#!/usr/bin/env python3
"""Are the divergent copies functional? (ADDENDUM1_FUNCTION.md)

(a) Irregularities in SSU and LSU: substitutions per 100 aligned positions and internal indel events per region,
    each copy against the nearest main-type copy of its culture (the comparator of copy_artefacts.py, same
    alignment scores). Irregular = above the 99th percentile of the baseline main-type copies.
(b) Secondary structure: cmalign (Infernal, truncated alignment allowed) of the SSU flank, 5.8S and LSU flank to
    Rfam RF01960, RF00002 and RF02543. A consensus base pair the copy covers on both sides is broken when it is
    not Watson-Crick or G-U, or one side is deleted. Degraded = broken-pair fraction above the baseline's 99th
    percentile. A region whose model is not given is skipped and reported as such.
(c) ITS2 folding (descriptive): RNAfold MFE per nt, and the base-pair distance to the comparator's structure
    after mapping the comparator's pairs onto the copy through the pairwise alignment.

  python3 stage3/copy_function.py --ref reference/culture_ref_v4 \\
      --artefacts benchmarks/copy_checks/artefacts --cm-58s benchmarks/stefani/RF00002.cm \\
      --cm-ssu RF01960.cm --cm-lsu RF02543.cm --out benchmarks/copy_checks/function --cpu 8
"""
import argparse
import collections
import hashlib
import json
import os
import statistics
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arms import read_fasta, read_tsv, write_tsv  # noqa: E402
from copy_artefacts import make_aligner, REGIONS  # noqa: E402

STRUCT = [('ssu_flank', 'cm_ssu', 'RF01960'), ('s58', 'cm_58s', 'RF00002'), ('lsu_flank', 'cm_lsu', 'RF02543')]
CANON = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}
OPEN, CLOSE = '<([{', '>)]}'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def p99(vals):
    return statistics.quantiles(vals, n=100, method='exclusive')[98] if len(vals) >= 2 else 0.0


# ------------------------------------------------------------------ (a) substitutions and indels per region
def region_stats(aligner, query, target, bounds):
    """Per comparator region: aligned positions, substitutions, internal indel events."""
    aln = aligner.align(target, query)[0]
    tb, qb = aln.aligned
    st = {r: [0, 0, 0] for r, _, _ in bounds}

    def reg(pos):
        for r, s, e in bounds:
            if s <= pos < e:
                return r
        return bounds[-1][0]
    mapping = {}
    for (ts, te), (qs, qe) in zip(tb, qb):
        for k in range(te - ts):
            t, q = target[ts + k], query[qs + k]
            mapping[ts + k] = qs + k
            r = reg(ts + k)
            st[r][0] += 1
            if t != q and t in 'ACGT' and q in 'ACGT':
                st[r][1] += 1
    for i in range(1, len(tb)):
        if tb[i][0] > tb[i - 1][1]:
            st[reg(tb[i - 1][1])][2] += 1
        if qb[i][0] > qb[i - 1][1]:
            st[reg(max(tb[i][0] - 1, 0))][2] += 1
    return st, mapping


# ------------------------------------------------------------------ (b) cmalign and broken pairs
def pairs_from_ss(ss):
    stacks, pairs = collections.defaultdict(list), {}
    for i, c in enumerate(ss):
        if c in OPEN:
            stacks[c].append(i)
        elif c in CLOSE:
            j = stacks[OPEN[CLOSE.index(c)]].pop()
            pairs[j] = i
        elif c.isalpha() and c.isupper():
            stacks['pk' + c].append(i)
        elif c.isalpha() and c.islower():
            j = stacks['pk' + c.upper()].pop()
            pairs[j] = i
    return pairs


def run_cmalign(cm, seqs, workdir, tag, cpu):
    fa = os.path.join(workdir, f'{tag}.fasta')
    with open(fa, 'w') as fh:
        for k, s in seqs.items():
            fh.write(f'>{k}\n{s}\n')
    out, sf = os.path.join(workdir, f'{tag}.pfam'), os.path.join(workdir, f'{tag}.scores')
    cmd = ['cmalign', '--cpu', str(cpu), '--outformat', 'Pfam', '--sfile', sf, '-o', out, cm, fa]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    rows, ss, rf = {}, '', ''
    for line in open(out):
        if line.startswith('#=GC SS_cons'):
            ss = line.split()[-1]
        elif line.startswith('#=GC RF'):
            rf = line.split()[-1]
        elif line.strip() and not line.startswith('#') and not line.startswith('//'):
            k, s = line.split()
            rows[k] = s
    scores = {}
    for line in open(sf):
        if line.startswith('#') or not line.strip():
            continue
        c = line.split()
        scores[c[1]] = (float(c[6]), int(c[2]))          # bit score, sequence length
    return rows, ss, rf, scores, ' '.join(cmd)


def broken_pairs(row, ss, rf):
    cons = [i for i, c in enumerate(rf) if c not in '.~']
    res = [i for i in cons if row[i] not in '-.~']
    if not res:
        return 0, 0, {}
    first, last = res[0], res[-1]
    pairs = pairs_from_ss(ss)
    n = br = 0
    for i, j in pairs.items():
        if i < first or j > last or row[i] == '~' or row[j] == '~':
            continue
        a, b = row[i].upper().replace('T', 'U'), row[j].upper().replace('T', 'U')
        n += 1
        if a in '-.' or b in '-.' or (a + b) not in CANON:
            br += 1
    return br, n, pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--artefacts', required=True)
    ap.add_argument('--cm-ssu', default='')
    ap.add_argument('--cm-58s', default='')
    ap.add_argument('--cm-lsu', default='')
    ap.add_argument('--out', required=True)
    ap.add_argument('--cpu', type=int, default=4)
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, 'cmalign'), exist_ok=True)
    t0 = time.time()
    rows = {r['accession']: r for r in read_tsv(os.path.join(a.ref, 'copies.tsv')) if r['status'] == 'retained'}
    art = {r['accession']: r for r in read_tsv(os.path.join(a.artefacts, 'artefacts.copies.tsv'))}
    flagged = set(json.load(open(os.path.join(a.artefacts, 'artefacts.json')))
                  ['main_copies_with_long_deletion_in_conserved_region'])
    seq = {m: read_fasta(os.path.join(a.ref, 'fasta', f'{m}.fasta')) for m in REGIONS + ['full']}

    def group(acc):
        r = rows[acc]
        if r['copy_class'] == 'divergent':
            return 'short 5.8S' if r['div_kind'] == 'short_5.8S' else 'same length'
        return 'flagged main' if acc in flagged else 'baseline'

    # (a)
    aligner = make_aligner()
    out = {}
    for acc, r in sorted(rows.items()):
        comp = art[acc]['comparator']
        o = {'accession': acc, 'culture': r['culture'], 'organism': r['organism'], 'group': group(acc),
             'comparator': comp}
        if comp:
            bounds, s = [], 0
            for m in REGIONS:
                L = len(seq[m][comp])
                bounds.append((m, s, s + L))
                s += L
            st, mapping = region_stats(aligner, seq['full'][acc], seq['full'][comp], bounds)
            for m in ('ssu_flank', 's58', 'lsu_flank'):
                al, sub, ind = st[m]
                o[f'{m}_aligned'] = al
                o[f'{m}_subs_per100'] = round(100 * sub / al, 3) if al else ''
                o[f'{m}_indels'] = ind
            o['_map'], o['_bounds'] = mapping, bounds
        out[acc] = o
    base = [o for o in out.values() if o['group'] == 'baseline' and o['comparator']]
    thr = {}
    for m in ('ssu_flank', 's58', 'lsu_flank'):
        thr[f'{m}_subs_per100'] = p99([o[f'{m}_subs_per100'] for o in base])
        thr[f'{m}_indels'] = p99([o[f'{m}_indels'] for o in base])
    for o in out.values():
        for m in ('ssu_flank', 's58', 'lsu_flank'):
            if o['comparator']:
                o[f'{m}_irregular'] = 'yes' if (o[f'{m}_subs_per100'] > thr[f'{m}_subs_per100']
                                                or o[f'{m}_indels'] > thr[f'{m}_indels']) else ''

    # (b)
    meta_cm, cm_rows, cm_ss = {}, {}, {}
    for m, opt, acc_rf in STRUCT:
        cm = getattr(a, opt)
        if not cm:
            meta_cm[m] = 'not run (no model given)'
            continue
        uniq = {}
        rep = {}
        for acc in sorted(rows):
            s = seq[m][acc]
            rep[acc] = uniq.setdefault(s, acc)
        useqs = {v: k for k, v in uniq.items()}
        r_, ss, rf, scores, cmd = run_cmalign(cm, useqs, os.path.join(a.out, 'cmalign'), m, a.cpu)
        meta_cm[m] = {'model': acc_rf, 'sha256': sha256(cm), 'command': cmd, 'unique_sequences': len(useqs)}
        cm_rows[m], cm_ss[m] = (r_, rep), (ss, rf)
        for acc, o in out.items():
            row = r_[rep[acc]]
            br, n, _ = broken_pairs(row, ss, rf)
            o[f'{m}_pairs'] = n
            o[f'{m}_broken'] = br
            o[f'{m}_broken_frac'] = round(br / n, 4) if n else ''
            sc, L = scores[rep[acc]]
            o[f'{m}_bits_per_nt'] = round(sc / L, 4) if L else ''
        vals = [o[f'{m}_broken_frac'] for o in out.values() if o['group'] == 'baseline' and o['comparator'] and o[f'{m}_broken_frac'] != '']
        thr[f'{m}_broken_frac'] = p99(vals)
        for o in out.values():
            o[f'{m}_degraded'] = 'yes' if o[f'{m}_broken_frac'] != '' and o[f'{m}_broken_frac'] > thr[f'{m}_broken_frac'] else ''

    # 5.8S deletions on the RF00002 consensus
    drows = []
    if 's58' in cm_rows:
        (r_, rep), (ss, rf) = cm_rows['s58'], cm_ss['s58']
        cons_idx = {}
        k = 0
        for i, c in enumerate(rf):
            if c not in '.~':
                k += 1
                cons_idx[i] = k
        pairs = pairs_from_ss(ss)
        paired = set(pairs) | set(pairs.values())
        for acc, o in sorted(out.items()):
            if o['group'] == 'baseline' or not o['comparator']:
                continue
            row, crow = r_[rep[acc]], r_[rep[o['comparator']]]
            dels = [i for i in cons_idx if row[i] == '-' and crow[i] not in '-.~']
            if not dels:
                continue
            drows.append({'accession': acc, 'group': o['group'], 'organism': o['organism'],
                          'deleted_consensus_positions': ','.join(str(cons_idx[i]) for i in dels),
                          'n_deleted': len(dels), 'n_paired': sum(i in paired for i in dels),
                          'ss_cons': ''.join(ss[i] for i in dels)})
        write_tsv(os.path.join(a.out, 'function.deletions.tsv'), drows,
                  ['accession', 'group', 'organism', 'n_deleted', 'n_paired', 'deleted_consensus_positions', 'ss_cons'])

    # (c) ITS2 folding
    import RNA
    fold = {}
    for acc in rows:
        s = seq['its2'][acc]
        if s not in fold:
            st, mfe = RNA.fold(s.replace('T', 'U'))
            fold[s] = (st, mfe)
    for acc, o in out.items():
        s = seq['its2'][acc]
        st, mfe = fold[s]
        o['its2_mfe_per_nt'] = round(mfe / len(s), 4) if s else ''
        comp = o['comparator']
        if comp:
            off_c = len(seq['ssu_flank'][comp]) + len(seq['its1'][comp]) + len(seq['s58'][comp])
            off_q = len(seq['ssu_flank'][acc]) + len(seq['its1'][acc]) + len(seq['s58'][acc])
            cst = fold[seq['its2'][comp]][0]
            cp = {(i, j) for i, j in pairs_from_ss(cst).items()}
            qp = {(i, j) for i, j in pairs_from_ss(st).items()}
            mp = o['_map']
            mapped = set()
            for i, j in cp:
                qi, qj = mp.get(off_c + i), mp.get(off_c + j)
                if qi is not None and qj is not None:
                    mapped.add((qi - off_q, qj - off_q))
            o['its2_bp_distance'] = len(mapped ^ qp)
            o['its2_bp_distance_per_nt'] = round(len(mapped ^ qp) / len(s), 4)

    # summary
    groups = ['short 5.8S', 'same length', 'flagged main', 'baseline']
    from scipy.stats import wilcoxon
    srows = []
    for m in ('ssu_flank', 's58', 'lsu_flank'):
        for g in groups:
            xs = [o for o in out.values() if o['group'] == g and o['comparator']]
            row = {'region': m, 'group': g, 'copies': len(xs),
                   'irregular': sum(o.get(f'{m}_irregular') == 'yes' for o in xs),
                   'median_subs_per100': round(statistics.median([o[f'{m}_subs_per100'] for o in xs]), 3) if xs else ''}
            if f'{m}_broken_frac' in thr:
                bf = [o[f'{m}_broken_frac'] for o in xs if o[f'{m}_broken_frac'] != '']
                row.update(degraded=sum(o[f'{m}_degraded'] == 'yes' for o in xs),
                           median_broken_frac=round(statistics.median(bf), 4) if bf else '')
                if g in ('short 5.8S', 'same length'):
                    d = []
                    for o in xs:
                        mains = [p[f'{m}_broken_frac'] for p in out.values() if p['culture'] == o['culture']
                                 and rows[p['accession']]['copy_class'] == 'main' and p[f'{m}_broken_frac'] != '']
                        if mains and o[f'{m}_broken_frac'] != '':
                            d.append(o[f'{m}_broken_frac'] - statistics.mean(mains))
                    nz = [x for x in d if x != 0]
                    row['paired_median_diff'] = round(statistics.median(d), 4) if d else ''
                    row['paired_wilcoxon_p'] = f'{wilcoxon(nz).pvalue:.3g}' if len(nz) >= 5 else ''
            srows.append(row)
    write_tsv(os.path.join(a.out, 'function.summary.tsv'), srows,
              ['region', 'group', 'copies', 'irregular', 'median_subs_per100', 'degraded', 'median_broken_frac',
               'paired_median_diff', 'paired_wilcoxon_p'])

    # verdicts (ADDENDUM1, section b)
    verdicts = {}
    for g in ('short 5.8S', 'same length'):
        regs = [m for m, _, _ in STRUCT if f'{m}_broken_frac' in thr]
        if len(regs) < 3:
            verdicts[g] = f'incomplete: structure run for {",".join(regs) or "no region"}'
            continue
        xs = [o for o in out.values() if o['group'] == g and o['comparator']]
        share = {m: sum(o[f'{m}_degraded'] == 'yes' for o in xs) / len(xs) for m in regs}
        dl = [d for d in drows if d['group'] == g]
        paired_del = any(d['n_paired'] > 0 for d in dl)
        if all(v <= 0.10 for v in share.values()) and not paired_del:
            verdicts[g] = 'structurally intact'
        elif any(v >= 0.50 for v in share.values()) or paired_del:
            verdicts[g] = 'structurally degraded'
        else:
            verdicts[g] = 'mixed'
        verdicts[g + ' (shares degraded)'] = {m: round(v, 3) for m, v in share.items()}

    cols = ['accession', 'culture', 'organism', 'group', 'comparator']
    for m in ('ssu_flank', 's58', 'lsu_flank'):
        cols += [f'{m}_aligned', f'{m}_subs_per100', f'{m}_indels', f'{m}_irregular']
        if f'{m}_broken_frac' in thr:
            cols += [f'{m}_pairs', f'{m}_broken', f'{m}_broken_frac', f'{m}_degraded', f'{m}_bits_per_nt']
    cols += ['its2_mfe_per_nt', 'its2_bp_distance', 'its2_bp_distance_per_nt']
    write_tsv(os.path.join(a.out, 'function.copies.tsv'), list(out.values()), cols)
    import Bio
    meta = {'finished': time.strftime('%Y-%m-%d %H:%M:%S'), 'seconds': round(time.time() - t0, 1),
            'infernal': subprocess.run(['cmalign', '-h'], capture_output=True, text=True).stdout.split('\n')[1].strip(),
            'viennarna': RNA.__version__, 'biopython': Bio.__version__, 'models': meta_cm,
            'thresholds_p99_baseline': thr, 'verdicts': verdicts,
            'outputs': {f: sha256(os.path.join(a.out, f)) for f in os.listdir(a.out) if f.endswith('.tsv')}}
    with open(os.path.join(a.out, 'function.json'), 'w') as fh:
        json.dump(meta, fh, indent=1)
    for r in srows:
        print('\t'.join(str(r.get(k, '')) for k in ('region', 'group', 'copies', 'irregular', 'median_subs_per100',
                                                    'degraded', 'median_broken_frac', 'paired_median_diff',
                                                    'paired_wilcoxon_p')))
    print('verdicts:', json.dumps(verdicts))


if __name__ == '__main__':
    main()
