#!/usr/bin/env python3
"""Copy types within and across cultures: mixed cultures, or rDNA types inside one genome?

Selects cultures (by modal VT and/or a culture regex), clusters their copies on one marker
(single linkage, default LROR-FLR2 at <= 1% distance), then reports:
  1. the types: size, cultures, copy classes, 5.8S length, nearest V18 reference (from crosswalk)
  2. types per culture, and which types turn up in every selected culture
  3. distances between types per marker (SSU flank, ITS1, 5.8S, ITS2, LROR-FLR2), within a
     culture and between cultures
  4. each culture's minority types: distance to its own majority type against distance to the
     same type in other cultures

Reading it:
  intragenomic types  every culture carries every type (in varying ratios); a type is near-identical
                      across cultures; two types differ by the same amount inside a culture as between
  mixed culture       a type turns up in some cultures only, usually as a minority, and matches another
                      culture's majority type on every marker, SSU included
A type present in every culture is the strongest sign of intragenomic types: independent
contamination would not recur in all of them. Nothing here proves either case.

Distances are % edit distance with the shorter sequence aligned inside the longer (Stage 1 infix).

Usage (from amfsplit_stage3):
  python3 copy_types.py --ref ../ref/culture_ref/v3 --crosswalk ../bench/stefani/crosswalk_v3.copies.tsv \\
      --vt VTX00064 --out ../bench/stefani/types_VTX00064
"""
import argparse
import collections
import itertools
import os
import re
import statistics
import sys

import edlib

from arms import read_fasta, read_tsv, write_tsv

MARKERS = [("ssu_flank", "SSU flank"), ("its1", "ITS1"), ("s58", "5.8S"), ("its2", "ITS2"),
           ("LSU_LROR-FLR2", "LROR-FLR2")]


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode="HW", task="distance")["editDistance"] / len(q)


def med(xs):
    return statistics.median(xs) if xs else None


def f1(x):
    return f"{x:.1f}" if x is not None else "-"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="frozen culture reference directory")
    ap.add_argument("--crosswalk", help="crosswalk copies table, for each copy's nearest V18 reference")
    ap.add_argument("--vt", nargs="*", default=["VTX00064"], help="select cultures whose modal VT is one of these")
    ap.add_argument("--culture", help="also select cultures matching this regex (e.g. '4292_708F|194757')")
    ap.add_argument("--cluster-on", default="LSU_LROR-FLR2", choices=[m for m, _ in MARKERS])
    ap.add_argument("--cut", type=float, default=1.0, help="single-linkage distance cut, %% (default 1)")
    ap.add_argument("--out", help="prefix for <out>.types.tsv and <out>.pairs.tsv")
    a = ap.parse_args()

    cul = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    pat = re.compile(a.culture) if a.culture else None
    chosen = {cu for cu, r in cul.items() if r["vt_modal"] in (a.vt or []) or (pat and pat.search(cu))}
    if not chosen:
        raise SystemExit("no culture matches the selection")
    cop = {r["accession"]: r for r in read_tsv(os.path.join(a.ref, "copies.tsv"))
           if r["status"] == "retained" and r["culture"] in chosen}
    lsu = {}
    if a.crosswalk:
        lsu = {r["accession"]: r for r in read_tsv(a.crosswalk) if r["accession"] in cop}
    seqs = {}
    for m, _ in MARKERS:
        p = os.path.join(a.ref, "fasta", f"{m}.fasta")
        seqs[m] = {k: v for k, v in read_fasta(p).items() if k in cop} if os.path.exists(p) else {}

    # ---- single-linkage types on one marker
    s = seqs[a.cluster_on]
    parent = {x: x for x in cop}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    have = sorted(x for x in cop if x in s)
    for x, y in itertools.combinations(have, 2):
        if dist(s[x], s[y]) <= a.cut:
            parent[find(x)] = find(y)
    groups = collections.defaultdict(list)
    for x in have:
        groups[find(x)].append(x)
    ordered = sorted(groups.values(), key=lambda g: (-len(g), min(g)))
    typ = {x: f"T{i + 1}" for i, g in enumerate(ordered) for x in g}
    for x in cop:
        typ.setdefault(x, "none")
    types = [f"T{i + 1}" for i in range(len(ordered))] + (["none"] if "none" in typ.values() else [])

    def cls(x):
        c = cop[x]
        return c["div_kind"] if c["copy_class"] == "divergent" and c.get("div_kind") else c["copy_class"]

    cultures = sorted(chosen)
    vts = collections.Counter(cul[cu]["vt_modal"] for cu in cultures)
    print(f"{len(cop)} copies in {len(cultures)} cultures (modal VT: {', '.join(f'{k} {v}' for k, v in vts.items())}); "
          f"types by single linkage on {a.cluster_on} at <= {a.cut:g}%")

    # ---- 1. types
    print(f"\n1. types\n  {'type':6s}{'copies':>7s}{'cultures':>9s}   {'nearest V18 reference (median id)':38s}"
          f"{'classes':28s}5.8S bp")
    for t in types:
        xs = [x for x in cop if typ[x] == t]
        near = collections.Counter(lsu[x]["LSU_best_label"] for x in xs if x in lsu and lsu[x]["LSU_best_label"])
        top = near.most_common(1)[0][0] if near else ""
        ids = [float(lsu[x]["LSU_best_pident"]) for x in xs if x in lsu and lsu[x]["LSU_best_label"] == top]
        nearest = f"{top} ({f1(med(ids))})" + (f" +{len(near) - 1} other" if len(near) > 1 else "") if top else "-"
        classes = ", ".join(f"{k} {v}" for k, v in collections.Counter(cls(x) for x in xs).most_common())
        lens = collections.Counter(cop[x]["s58_len"] for x in xs if cop[x].get("s58_len"))
        print(f"  {t:6s}{len(xs):>7d}{len({cop[x]['culture'] for x in xs}):>9d}   {nearest[:37]:38s}{classes[:27]:28s}"
              + ",".join(sorted(lens)))

    # ---- 2. types per culture
    print("\n2. types per culture (m = main, d = divergent)")
    print(f"  {'culture':44s}" + "".join(f"{t:>9s}" for t in types))
    everywhere = []
    for t in types:
        if all(any(typ[x] == t and cop[x]["culture"] == cu for x in cop) for cu in cultures):
            everywhere.append(t)
    for cu in cultures:
        cells = []
        for t in types:
            xs = [x for x in cop if cop[x]["culture"] == cu and typ[x] == t]
            nm = sum(cop[x]["copy_class"] == "main" for x in xs)
            cells.append(f"{nm}m{len(xs) - nm}d" if xs else ".")
        print(f"  {cu[:43]:44s}" + "".join(f"{c:>9s}" for c in cells))
    print(f"  types found in every selected culture: {', '.join(everywhere) or 'none'}")

    # ---- 3. distances between types, within and between cultures
    pairs = []
    for x, y in itertools.combinations(sorted(cop), 2):
        row = {"a": x, "b": y, "culture_a": cop[x]["culture"], "culture_b": cop[y]["culture"],
               "type_a": typ[x], "type_b": typ[y], "same_culture": cop[x]["culture"] == cop[y]["culture"]}
        for m, _ in MARKERS:
            row[m] = round(dist(seqs[m][x], seqs[m][y]), 2) if x in seqs[m] and y in seqs[m] else ""
        pairs.append(row)
    print(f"\n3. distance between types, median % over copy pairs")
    print(f"  {'types':9s}{'where':9s}{'pairs':>6s}" + "".join(f"{lab:>11s}" for _, lab in MARKERS))
    real = [t for t in types if t != "none"]
    for t1, t2 in itertools.combinations_with_replacement(real, 2):
        for where, same in (("within", True), ("between", False)):
            ps = [p for p in pairs if {p["type_a"], p["type_b"]} == {t1, t2} and p["same_culture"] == same
                  and (t1 != t2 or p["type_a"] == p["type_b"])]
            if not ps:
                continue
            cells = [f1(med([p[m] for p in ps if p[m] != ""])) for m, _ in MARKERS]
            print(f"  {t1 + '-' + t2:9s}{where:9s}{len(ps):>6d}" + "".join(f"{c:>11s}" for c in cells))

    # ---- 4. minority types against own majority and against the same type elsewhere
    print("\n4. minority types: median % to the culture's own majority type vs to the same type in other cultures")
    for cu in cultures:
        xs = [x for x in cop if cop[x]["culture"] == cu and typ[x] != "none"]
        cnt = collections.Counter(typ[x] for x in xs)
        if len(cnt) < 2:
            continue
        major = cnt.most_common(1)[0][0]
        for t in [t for t in cnt if t != major]:
            mine = [x for x in xs if typ[x] == t]
            own = [x for x in xs if typ[x] == major]
            other = [x for x in cop if typ[x] == t and cop[x]["culture"] != cu]
            print(f"  {cu}: {t} ({len(mine)}) in a culture of mostly {major} ({cnt[major]})"
                  + ("" if other else f"; {t} occurs in no other selected culture"))
            for label, ref in ((f"to own {major}", own), (f"to {t} elsewhere ({len(other)})", other)):
                if not ref:
                    continue
                cells = []
                for m, lab in MARKERS:
                    ds = [dist(seqs[m][p], seqs[m][q]) for p in mine for q in ref if p in seqs[m] and q in seqs[m]]
                    cells.append(f"{lab} {f1(med(ds))}")
                print(f"      {label:24s}" + " | ".join(cells))

    if a.out:
        trows = [{"accession": x, "culture": cop[x]["culture"], "organism": cop[x]["organism"],
                  "copy_class": cop[x]["copy_class"], "div_kind": cop[x].get("div_kind", ""),
                  "s58_len": cop[x].get("s58_len", ""), "type": typ[x],
                  "LSU_best_label": lsu.get(x, {}).get("LSU_best_label", ""),
                  "LSU_best_pident": lsu.get(x, {}).get("LSU_best_pident", "")} for x in sorted(cop)]
        write_tsv(f"{a.out}.types.tsv", trows, list(trows[0]))
        write_tsv(f"{a.out}.pairs.tsv", pairs, list(pairs[0]))
        print(f"\ntables: {a.out}.types.tsv, {a.out}.pairs.tsv")


if __name__ == "__main__":
    sys.exit(main())
