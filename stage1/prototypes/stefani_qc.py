#!/usr/bin/env python3
"""Reference hygiene checks on the Stefani et al. (2025) rDNA copies before they seed Stage 2.

1. SSU outliers: copies whose SSU flank sits far from their culture's medoid copy. Intragenomic
   SSU variation is ~0.1%, so a copy several % away is more likely a contaminant, chimera or
   mislabel than a genuine rDNA variant.
2. Divergent 5.8S copies: copies whose 5.8S sits far from the culture medoid. The 5.8S is under
   strong functional constraint, so a divergent 5.8S in a subset of copies is a classic sign of
   non-functional rDNA. Reports whether these sequences are shared across species and whether
   they drive each culture's ITS2 maximum.
3. The recurring V4 difference: where the edits sit in the most divergent AMV4.5NF-AMDGR pair of
   cultures whose V4 maximum falls in a given range (default 2.0-2.6%, about 5 edits).

All distances are infix (shorter sequence aligned inside the longer, end gaps free).

Usage (from amfsplit_stage1):
  python3 stefani_qc.py            # also writes ../bench/stefani/stefani_qc.flags.tsv
  python3 stefani_qc.py --groups ../bench/stefani/stefani_groups.tsv --prefix ../bench/stefani/out/windows/stefani
"""
import argparse
import collections
import csv
import itertools
import math
import statistics

import edlib

V4 = "SSU_V4_AMV4.5NF-AMDGR"


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks).upper()
    return seqs


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode="HW", task="distance")["editDistance"] / len(q)


def medoid_dists(seqs):
    """Distance of every copy to the culture medoid (the copy with the smallest summed distance)."""
    accs = list(seqs)
    pair = {frozenset(p): dist(seqs[p[0]], seqs[p[1]]) for p in itertools.combinations(accs, 2)}
    d = lambda a, b: 0.0 if a == b else pair[frozenset((a, b))]
    med = min(accs, key=lambda a: sum(d(a, b) for b in accs))
    return {a: d(a, med) for a in accs}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", default="../bench/stefani/stefani_groups.tsv")
    ap.add_argument("--prefix", default="../bench/stefani/out/windows/stefani")
    ap.add_argument("--min-copies", type=int, default=3)
    ap.add_argument("--ssu-cut", type=float, default=2.0, help="SSU flank %% from medoid to flag (default 2)")
    ap.add_argument("--s58-cut", type=float, default=3.0, help="5.8S %% from medoid to flag (default 3)")
    ap.add_argument("--v4-range", type=float, nargs=2, default=[2.0, 2.6])
    ap.add_argument("--show", type=int, default=6, help="V4 alignments to print")
    ap.add_argument("--out", default="../bench/stefani/stefani_qc",
                    help="prefix for <out>.flags.tsv (accession, culture, check, pct_from_medoid)")
    a = ap.parse_args()

    grp = {r["accession"]: r for r in csv.DictReader(open(a.groups), delimiter="\t")}
    rec = {r["read_id"]: r for r in csv.DictReader(open(f"{a.prefix}.records.tsv"), delimiter="\t")}
    try:
        tri = {r["read_id"]: r for r in csv.DictReader(open(f"{a.prefix}.triage.tsv"), delimiter="\t")}
    except FileNotFoundError:
        tri = {}
    seq = {r: read_fasta(f"{a.prefix}.{r}.fasta") for r in ("ssu_flank", "s58", "its2", V4)}
    cultures = collections.defaultdict(list)
    for acc, r in grp.items():
        cultures[r["culture"]].append(acc)
    cultures = {c: v for c, v in sorted(cultures.items()) if len(v) >= a.min_copies}

    flags = []  # (accession, culture, check, pct)

    # 1. SSU outliers
    print(f"1. SSU flank outliers (> {a.ssu_cut}% from the culture medoid)")
    n_out = 0
    for c, accs in cultures.items():
        s = {x: seq["ssu_flank"][x] for x in accs if x in seq["ssu_flank"]}
        if len(s) < a.min_copies:
            continue
        for x, v in sorted(medoid_dists(s).items(), key=lambda kv: -kv[1]):
            if v > a.ssu_cut:
                n_out += 1
                flags.append((x, c, "ssu_outlier", v))
                t = tri.get(x, {})
                print(f"  {c[:46]:46s} {x:11s} {v:5.2f}%  LSU best hit {t.get('best_target', '-')[:38]} "
                      f"({t.get('identity', '-')})")
    print(f"  {n_out} copies flagged")

    # 2. divergent 5.8S copies
    print(f"\n2. Divergent 5.8S copies (> {a.s58_cut}% from the culture medoid)")
    flagged, drives, expected = {}, [], 0.0
    for c, accs in cultures.items():
        s = {x: seq["s58"][x] for x in accs if x in seq["s58"]}
        if len(s) < a.min_copies:
            continue
        f = {x: v for x, v in medoid_dists(s).items() if v > a.s58_cut}
        if not f:
            continue
        flagged.update(f)
        flags += [(x, c, "divergent_5.8S", v) for x, v in f.items()]
        s2 = {x: seq["its2"][x] for x in accs if x in seq["its2"]}
        if len(s2) >= 2:
            p = max(itertools.combinations(s2, 2), key=lambda pq: dist(s2[pq[0]], s2[pq[1]]))
            drives.append(any(x in f for x in p))
            n, k = len(s2), sum(x in f for x in s2)
            expected += 1 - math.comb(n - k, 2) / math.comb(n, 2)
    if flagged:
        by_genus = collections.Counter(grp[x]["genus"] for x in flagged)
        tot = collections.Counter(r["genus"] for r in grp.values())
        print(f"  {len(flagged)} copies in {len(drives)} cultures; by genus (flagged/all copies): " +
              ", ".join(f"{g} {n}/{tot[g]}" for g, n in by_genus.most_common()))
        med_len = lambda xs: statistics.median(len(seq["s58"][x]) for x in xs)
        rest = [x for x in grp if x in seq["s58"] and x not in flagged]
        print(f"  median 5.8S length: flagged {med_len(flagged)} bp, all others {med_len(rest)} bp")
        var = collections.defaultdict(set)
        for x in flagged:
            var[seq["s58"][x]].add(grp[x]["organism"])
        multi = sorted((v for v in var.values() if len(v) > 1), key=len, reverse=True)
        print(f"  {len(var)} distinct divergent 5.8S sequences; {len(multi)} found in more than one species"
              + (f", widest: {sorted(multi[0])}" if multi else ""))
        print(f"  the culture's most divergent ITS2 pair includes a flagged copy in {sum(drives)} of "
              f"{len(drives)} cultures ({expected:.1f} expected if flagged copies were random)")
    else:
        print("  none")

    # 3. V4 edits
    lo, hi = a.v4_range
    print(f"\n3. Most divergent {V4} pair in cultures with a V4 max of {lo}-{hi}% (complete windows)")
    shown = 0
    for c, accs in cultures.items():
        s = {x: seq[V4][x] for x in accs if x in seq[V4] and rec.get(x, {}).get(f"{V4}_status") == "ok"}
        if len(s) < 2:
            continue
        d, p, q = max((dist(s[p], s[q]), p, q) for p, q in itertools.combinations(s, 2))
        if not lo <= d <= hi:
            continue
        qs, ts = (s[p], s[q]) if len(s[p]) <= len(s[q]) else (s[q], s[p])
        aln = edlib.getNiceAlignment(edlib.align(qs, ts, mode="HW", task="path"), qs, ts)
        m = aln["matched_aligned"]
        pos = [i for i, ch in enumerate(m) if ch != "|"]
        print(f"  {c[:52]}  {p} vs {q}  {d:.2f}%, edits at aligned columns {pos}")
        i, j = max(0, pos[0] - 6), pos[-1] + 7
        if j - i <= 90:
            for key in ("query_aligned", "matched_aligned", "target_aligned"):
                print(f"      {aln[key][i:j]}")
        shown += 1
        if shown >= a.show:
            break
    if not shown:
        print("  none in range")

    with open(f"{a.out}.flags.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "culture", "check", "pct_from_medoid"])
        w.writerows((x, c, k, f"{v:.2f}") for x, c, k, v in flags)
    print(f"\nflags: {a.out}.flags.tsv ({len(flags)} rows)")


if __name__ == "__main__":
    main()
