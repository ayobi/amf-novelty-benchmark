#!/usr/bin/env python3
"""Phase 2, SSU arm: compatible VT sets with a calibrated margin δ.

Stage 1 showed the best-hit VT is fragile (71% of assigned copies have the next VT within 1 point).
Here each copy gets a set instead. The anchor is the copy's Stage 1 best hit (highest bit score among
hits passing >= 97% identity, >= 95% coverage of the shorter sequence, e-value <= 1e-50; Ns in the
reference not counted, hits to culture-reference accessions dropped). The set holds every VT whose
best identity is at least (anchor identity - δ), and >= 97%. The Stage 1 VT is therefore always in
the set, and so is any VT with a higher identity over a shorter alignment. A copy with no passing hit
gets an empty set (unassigned).

Merged units (--merge auto, the default): VTs that split the main-type copies of one culture at δ = 0
are joined into one unit (e.g. VTX00113+VTX00115) when the closest copies across the split are no
further apart in SSU than --merge-max-dist (default 1.0%, about the 99th percentile of within-culture
SSU spread in the Stefani set, 0.95%; a fixed value, so a mixed culture cannot raise it).
These are the inseparable pairs the culture reference shows directly; merging them keeps δ from being
widened for everyone to cover a few pairs. A split with the two sides further apart (a mixed culture
such as an unsplit BR105_2435A) is reported and not merged.

Calibration, on main-type copies. A culture with >= 2 main-type copies that have a set is consistent
when those sets share a unit. For each δ on a grid the report shows consistency and the cost:
copies left with one unit, mean set size, and named species still resolved (no other named species
shares any of their units). The chosen δ is the larger of
  (a) the smallest δ at which >= --target of cultures are consistent (default 0.99), and
  (b) the --spread-quantile percentile of within-culture SSU spread (default 95th), rounded up to
      the grid: two copies d% apart can differ by up to d points in identity to any reference,
      so a single read needs that much room even when whole genomes are consistent at a smaller δ.
--delta fixes δ instead.

Reports:
  1. within-culture SSU spread of main-type copies
  2. merged units
  3. calibration table and the chosen δ
  4. at the chosen δ: multi-unit culture sets, cultures still inconsistent, units shared by species
  5. divergent-class copies: does their set contain their culture's unit?
Writes {out}.copies.tsv, {out}.cultures.tsv and {out}.calibration.tsv.

Usage (from amfsplit_stage3):
  python3 vt_sets.py --ref ../ref/culture_ref/v4 --maarjam ../ref/maarjam/maarjam.fasta \\
      --out ../bench/stefani/vt_sets_v4 --threads 8
"""
import argparse
import collections
import itertools
import math
import os
import shutil
import statistics
import sys

import edlib

from arms import base_acc, is_named, load_reference, read_fasta, read_tsv, run_blast, write_tsv


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode="HW", task="distance")["editDistance"] / len(q)


def vt_hits(hs, refs, min_cov, max_ev, min_id, drop):
    """({VT: best identity}, anchor identity): identities over hits passing coverage and e-value; the
    anchor is the identity of the highest-bit-score hit that also passes min_id (None if none does)."""
    ids, best = {}, None
    for h in hs:
        info = refs[h["ref"]]
        if base_acc(info["acc"]) in drop:
            continue
        cov = h["length"] / min(h["qlen"], max(1, h["slen"] - info.get("n_N", 0)))
        if cov < min_cov or h["evalue"] > max_ev:
            continue
        v = info["label"]
        if h["pident"] > ids.get(v, -1):
            ids[v] = h["pident"]
        if h["pident"] >= min_id and (best is None or (h["bitscore"], h["pident"]) > best):
            best = (h["bitscore"], h["pident"])
    return ids, (best[1] if best else None)


def vt_set(ids, anchor, delta, min_id, unit):
    if anchor is None:
        return frozenset()
    thr = max(min_id, anchor - delta) - 1e-9
    return frozenset(unit.get(v, v) for v, p in ids.items() if p >= thr)


def union_find(groups):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for g in groups:
        g = sorted(g)
        for v in g:
            find(v)
        for v in g[1:]:
            parent[find(v)] = find(g[0])
    out = collections.defaultdict(set)
    for v in parent:
        out[find(v)].add(v)
    return list(out.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="frozen culture reference directory")
    ap.add_argument("--maarjam", required=True)
    ap.add_argument("--min-id", type=float, default=97.0)
    ap.add_argument("--min-cov", type=float, default=0.95)
    ap.add_argument("--max-evalue", type=float, default=1e-50)
    ap.add_argument("--merge", choices=["auto", "none"], default="auto",
                    help="join VTs that split one culture's main-type copies at δ = 0 (default auto)")
    ap.add_argument("--merge-max-dist", type=float, default=1.0,
                    help="max SSU %% distance between the closest copies across a split for it to merge (default 1.0)")
    ap.add_argument("--target", type=float, default=0.99, help="share of cultures that must be consistent")
    ap.add_argument("--spread-quantile", type=int, default=95,
                    help="percentile of within-culture SSU spread that δ must cover (0 = ignore)")
    ap.add_argument("--max-delta", type=float, default=3.0)
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--delta", type=float, help="use this δ instead of calibrating")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = read_tsv(os.path.join(a.ref, "copies.tsv"))
    cop = {r["accession"]: r for r in rows if r["status"] == "retained"}
    drop = {base_acc(r["accession"]) for r in rows}
    cul_meta = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    qpath = os.path.join(a.ref, "fasta", "ssu_flank.fasta")
    ssu = read_fasta(qpath)
    refs, fasta, _ = load_reference(a.maarjam, "maarjam")
    hits = run_blast(qpath, fasta, a.out + ".blastdb", "vt", "megablast", a.threads)
    shutil.rmtree(a.out + ".blastdb", ignore_errors=True)
    hv = {x: vt_hits(hits.get(x, []), refs, a.min_cov, a.max_evalue, a.min_id, drop) for x in cop}

    by_cul = collections.defaultdict(list)
    for x, r in cop.items():
        by_cul[r["culture"]].append(x)
    main_of = {cu: [x for x in xs if cop[x]["copy_class"] == "main" and x in ssu] for cu, xs in by_cul.items()}
    names = {cu: cul_meta[cu]["organism"] for cu in by_cul}
    named = {cu for cu in by_cul if is_named(names[cu]) and cul_meta[cu].get("name_status", "ok") == "ok"}
    print(f"{len(cop)} copies, {len(by_cul)} cultures; MaarjAM {len(refs)} sequences, "
          f"{len({v['label'] for v in refs.values()})} VT")

    # ---- 1. within-culture SSU spread
    spread = []
    for cu, xs in main_of.items():
        spread += [dist(ssu[p], ssu[q]) for p, q in itertools.combinations(xs, 2)]
    qs = statistics.quantiles(spread, n=100) if len(spread) > 1 else [0.0] * 99
    print(f"\n1. within-culture SSU flank spread, main-type copy pairs (n = {len(spread)}): median "
          f"{statistics.median(spread):.2f}%, 90th {qs[89]:.2f}%, 95th {qs[94]:.2f}%, 99th {qs[98]:.2f}%, "
          f"max {max(spread):.2f}%")

    # ---- 2. merged units from within-culture splits at δ = 0
    unit, why, kept_split = {}, collections.defaultdict(list), []
    merge_max = a.merge_max_dist
    if a.merge == "auto":
        splits = []
        for cu, xs in main_of.items():
            side = collections.defaultdict(list)
            for x in xs:
                s0 = vt_set(*hv[x], 0.0, a.min_id, {})
                if s0:
                    side[s0].append(x)
            if len(side) < 2 or frozenset.intersection(*side):
                continue
            gap = max(min(dist(ssu[p], ssu[q]) for p in side[s1] for q in side[s2])
                      for s1, s2 in itertools.combinations(side, 2))
            if gap <= merge_max:
                splits.append((cu, frozenset().union(*side), gap))
            else:
                kept_split.append((cu, sorted("/".join(sorted(k)) for k in side), gap))
        for g in union_find([s for _, s, _ in splits]):
            if len(g) > 1:
                name = "+".join(sorted(g))
                for v in g:
                    unit[v] = name
        for cu, s, gap in splits:
            why[unit[next(iter(s))]].append(f"{cu} ({gap:.2f}%)")
    print(f"\n2. merged units: VTs that split one culture's main-type copies at δ = 0, closest copies across "
          f"the split <= {merge_max:.2f}% apart in SSU" + (": none" if not why else ""))
    for u, cus in sorted(why.items()):
        print(f"      {u}: {', '.join(sorted(cus))}")
    for cu, sides, gap in kept_split:
        print(f"      not merged, {gap:.2f}% apart (mixed culture?): {cu}: {' vs '.join(sides)}")

    # ---- 3. calibration
    grid = [round(i * a.step, 3) for i in range(int(round(a.max_delta / a.step)) + 1)]
    if a.delta is not None and a.delta not in grid:
        grid = sorted(grid + [a.delta])
    all_named_species = {names[cu] for cu in named}
    cal, state = [], {}
    for d in grid:
        sets = {x: vt_set(*hv[x], d, a.min_id, unit) for x in cop}
        testable = consistent = 0
        culset = {}
        for cu, xs in main_of.items():
            ss = [sets[x] for x in xs if sets[x]]
            if not ss:
                culset[cu] = frozenset()
                continue
            inter = frozenset.intersection(*ss)
            culset[cu] = inter if inter else frozenset().union(*ss)
            if len(ss) >= 2:
                testable += 1
                consistent += bool(inter)
        filled = [s for s in sets.values() if s]
        groups = union_find([s for s in culset.values() if s])
        group_of = {u: i for i, g in enumerate(groups) for u in g}
        g_species = collections.defaultdict(set)
        for cu in named:
            for u in culset[cu]:
                g_species[group_of[u]].add(names[cu])
        sp_groups = collections.defaultdict(set)
        for cu in named:
            for u in culset[cu]:
                sp_groups[names[cu]].add(group_of[u])
        resolved = sum(1 for sp, gs in sp_groups.items() if all(len(g_species[g]) == 1 for g in gs))
        row = {"delta": d, "cultures_tested": testable, "consistent": consistent,
               "consistent_share": round(consistent / testable, 4) if testable else "",
               "copies_with_set": len(filled),
               "single_unit_share": round(sum(len(s) == 1 for s in filled) / len(filled), 4) if filled else "",
               "mean_set_size": round(statistics.mean(len(s) for s in filled), 3) if filled else "",
               "unit_groups": len(groups), "groups_with_several_species": sum(len(s) > 1 for s in g_species.values()),
               "species_resolved": resolved, "named_species": len(all_named_species)}
        cal.append(row)
        state[d] = (sets, culset, groups, g_species, group_of)
    ok = [r["delta"] for r in cal if r["consistent_share"] != "" and r["consistent_share"] >= a.target]
    d_cons = ok[0] if ok else grid[-1]
    d_spread = (math.ceil(round(qs[a.spread_quantile - 1] / a.step, 6)) * a.step) if a.spread_quantile else 0.0
    d_spread = round(min(d_spread, grid[-1]), 3)
    chosen = a.delta if a.delta is not None else max(d_cons, d_spread)
    print(f"\n3. calibration: a culture is consistent when its main-type copies' sets share a unit")
    print(f"  {'δ':>5s}{'consistent':>14s}{'one unit':>10s}{'mean set':>10s}{'species resolved':>18s}")
    show = {0.0, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, d_cons, d_spread, chosen}
    for r in cal:
        if r["delta"] in show:
            tags = [t for t, v in (("consistency", d_cons), ("spread", d_spread), ("chosen", chosen)) if r["delta"] == v]
            share = 100 * r["consistent_share"] if r["consistent_share"] != "" else float("nan")
            print(f"  {r['delta']:>5.1f}{r['consistent']:>6d}/{r['cultures_tested']:<3d}{share:>4.0f}%"
                  f"{100 * r['single_unit_share']:>9.0f}%{r['mean_set_size']:>10.2f}"
                  f"{r['species_resolved']:>11d} of {r['named_species']:<3d}" + (f"  <- {', '.join(tags)}" if tags else ""))
    if a.delta is not None:
        print(f"  chosen δ = {chosen:g} (fixed by --delta)")
    else:
        print(f"  consistency: smallest δ with >= {100 * a.target:.0f}% of cultures consistent = {d_cons:g}; "
              f"spread: {a.spread_quantile}th percentile of within-culture SSU spread, rounded up = {d_spread:g}")
        print(f"  chosen δ = {chosen:g} (the larger)")
    sets0 = state[grid[0]][0]
    inset = sum(cop[x]["vt_call"] == "assigned" and unit.get(cop[x]["vt"], cop[x]["vt"]) in sets0[x] for x in cop)
    print(f"  at δ = 0 the Stage 1 best-hit VT is in the set for {inset} of "
          f"{sum(cop[x]['vt_call'] == 'assigned' for x in cop)} assigned copies")

    # ---- 4. at the chosen δ
    sets, culset, groups, g_species, group_of = state[chosen]
    print(f"\n4. at δ = {chosen:g}")
    multi = collections.defaultdict(list)
    for cu, s in culset.items():
        if len(s) > 1:
            multi[s].append(cu)
    print(f"  culture sets with > 1 unit: {sum(len(v) for v in multi.values())} cultures in {len(multi)} distinct sets")
    for s, cus in sorted(multi.items(), key=lambda kv: -len(kv[1]))[:10]:
        orgs = collections.Counter(names[c] for c in cus)
        print(f"      {'/'.join(sorted(s))}: {', '.join(f'{o} ({n})' for o, n in orgs.most_common())}")
    bad = [cu for cu, xs in main_of.items() if len([x for x in xs if sets[x]]) >= 2
           and not frozenset.intersection(*[sets[x] for x in xs if sets[x]])]
    print(f"  cultures still inconsistent: {len(bad)}")
    for cu in sorted(bad)[:10]:
        cnt = collections.Counter("/".join(sorted(sets[x])) or "-" for x in main_of[cu])
        print(f"      {cu}: " + ", ".join(f"{k} x{n}" for k, n in cnt.most_common()))
    shared = [(i, sp) for i, sp in g_species.items() if len(sp) > 1]
    print(f"  unit groups shared by > 1 named species: {len(shared)} of {len(groups)}")
    for i, sp in sorted(shared, key=lambda t: -len(t[1]))[:10]:
        print(f"      {'/'.join(sorted(groups[i]))}: {', '.join(sorted(sp))}")

    # ---- 5. divergent-class copies
    print(f"\n5. non-main copies at δ = {chosen:g}: does the copy's set contain a unit of its culture's set?")
    for kind in ("short_5.8S", "same_length", "longer", ""):
        xs = [x for x in cop if cop[x]["copy_class"] == "divergent" and cop[x].get("div_kind", "") == kind]
        if not xs:
            continue
        yes = sum(bool(sets[x] & culset[cop[x]["culture"]]) for x in xs)
        empty = sum(not sets[x] for x in xs)
        nocul = sum(bool(sets[x]) and not culset[cop[x]["culture"]] for x in xs)
        other = [x for x in xs if sets[x] and culset[cop[x]["culture"]] and not sets[x] & culset[cop[x]["culture"]]]
        print(f"  divergent ({kind or 'no 5.8S length'}): {len(xs)} copies; shares a culture unit {yes}, "
              f"other units only {len(other)}, no VT {empty}, culture has no VT {nocul}")
        for x in other[:6]:
            print(f"      {x} {cop[x]['culture']}: {'/'.join(sorted(sets[x]))} (culture {'/'.join(sorted(culset[cop[x]['culture']]))})")

    # ---- tables
    crow = []
    for x in sorted(cop):
        ids, anchor = hv[x]
        top = sorted(ids.items(), key=lambda kv: -kv[1])[:5]
        crow.append({"accession": x, "culture": cop[x]["culture"], "organism": cop[x]["organism"],
                     "copy_class": cop[x]["copy_class"], "div_kind": cop[x].get("div_kind", ""),
                     "stage1_call": cop[x]["vt_call"], "stage1_vt": cop[x]["vt"],
                     "anchor_pident": anchor if anchor is not None else "",
                     "vt_set": ",".join(sorted(sets[x])), "set_size": len(sets[x]),
                     "culture_set": ",".join(sorted(culset[cop[x]["culture"]])),
                     "top_vts": ",".join(f"{v}:{p:.2f}" for v, p in top)})
    write_tsv(f"{a.out}.copies.tsv", crow, list(crow[0]))
    urow = []
    for cu in sorted(by_cul):
        xs = [x for x in main_of[cu] if sets[x]]
        urow.append({"culture": cu, "organism": names[cu], "named": cu in named, "main_with_set": len(xs),
                     "culture_set": ",".join(sorted(culset[cu])),
                     "consistent": (bool(frozenset.intersection(*[sets[x] for x in xs])) if len(xs) >= 2 else ""),
                     "copy_sets": ";".join(f"{k}:{n}" for k, n in collections.Counter(
                         "/".join(sorted(sets[x])) for x in xs).most_common())})
    write_tsv(f"{a.out}.cultures.tsv", urow, list(urow[0]))
    write_tsv(f"{a.out}.calibration.tsv", cal, list(cal[0]))
    print(f"\ntables: {a.out}.copies.tsv, .cultures.tsv, .calibration.tsv (δ = {chosen:g}, merge {a.merge})")


if __name__ == "__main__":
    sys.exit(main())
