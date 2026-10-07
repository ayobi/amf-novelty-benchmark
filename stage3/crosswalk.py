#!/usr/bin/env python3
"""Phase 1: VT / SH / LSU crosswalk from the frozen culture reference.

Every retained copy carries all three markers, so each copy gets three labels from physically
linked sequences rather than from name matching across databases:
  VT   SSU flank         vs MaarjAM   >= 97% identity, >= 95% coverage   (vt_check.py criteria)
  SH   full ITS          vs UNITE     >= 98.5%, >= 90%                   (the 1.5% SH distance)
  LSU  LROR-FLR2 window  vs Delavaux  >= 97%, >= 90%, nearest reference  (placement comes in Phase 2)
Coverage is of the shorter sequence (Ns in the reference not counted) and e-value <= 1e-50 on
every arm. Hits to any copy of the
culture reference (every accession in copies.tsv, excluded ones included) are dropped: newer UNITE
and INSD-derived releases may already hold the Stefani deposits, and a copy must not take its label
from itself or from a culture-mate (--keep-refset-hits turns this off). Arms without a reference
are skipped, so VT and LSU can run before UNITE is downloaded.

Counting unit is the culture from freeze_ref.py. A culture's label on an arm is the modal
assigned label of its main-type copies (ties go to the higher median identity); divergent-class
and SSU-outlier copies are reported against it. Cultures on name hold stay in the copy and
culture tables but are left out of species tallies, and 'sp.' organisms are not counted as named
species.

Reports (stdout):
  1. assignment per arm and order, including how many cultures UNITE gives an SH at all
  2. the VT arm against the Stage 1 vt_check calls carried in the reference (regression)
  3. many-to-many per arm: named species in > 1 unit, units holding > 1 named species
  4. within-culture splits: main-type copies of one culture in two units (inseparable pairs)
  5. cross-arm: units of one arm split by another
  6. divergent-class and SSU-outlier copies against their culture's main-type unit, per arm; with a
     reference frozen with 5.8S kinds, short-5.8S and same-length divergent copies are reported apart
  7. genus of the culture vs genus of the best reference hit (SH, LSU)
SH codes are shown with the SH's UNITE name, as 'SHcode (name)'.
Writes {out}.copies.tsv, {out}.cultures.tsv, {out}.species.tsv, {out}.units.tsv, and {out}.run.json
(arguments, culture_ref version and fingerprint, sha256 and record counts of every reference).

Usage (from amfsplit_stage3):
  python3 crosswalk.py --ref ../ref/culture_ref/v1 --maarjam ../ref/maarjam/maarjam.fasta \\
      --unite ../ref/unite/<UNITE general release>.fasta \\
      --lsu ../ref/delavaux/repo/2024_AMFPipeline/2024_AMFPipeline_ASV/V18_LSUDB_052025_AMFONLY.fasta \\
      --out ../bench/stefani/crosswalk_v1 --threads 8
"""
import argparse
import collections
import itertools
import json
import os
import shutil
import statistics
import sys

import edlib

from arms import (GENUS_SYNONYMS, assign, base_acc, is_named, load_reference, read_fasta, read_tsv,
                  run_blast, sha256, write_tsv)

ARMS = [  # name, loader, query region in <ref>/fasta, BLAST task, min id, min coverage
    ("VT", "maarjam", "ssu_flank", "megablast", 97.0, 0.95),
    ("SH", "unite", "full_its", "blastn", 98.5, 0.90),
    ("LSU", "delavaux", "LSU_LROR-FLR2", "blastn", 97.0, 0.90),
]
FIELDS = ["call", "label", "pident", "cov", "margin", "tied", "best_label", "best_pident",
          "best_taxon", "best_genus", "dropped"]


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode="HW", task="distance")["editDistance"] / len(q)


def fmt_counts(c, n=None, show=str):
    return ", ".join(f"{show(k)} ({v})" for k, v in c.most_common(n))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="frozen reference directory (freeze_ref.py)")
    ap.add_argument("--maarjam", help="MaarjAM FASTA (gDAT snapshot)")
    ap.add_argument("--unite", help="UNITE FASTA whose headers carry SH codes")
    ap.add_argument("--unite-keep", default="Glomeromycota",
                    help="keep only UNITE records whose header matches this regex ('' keeps all; default Glomeromycota)")
    ap.add_argument("--lsu", help="Delavaux LSU database FASTA (e.g. V18_LSUDB_052025_AMFONLY.fasta)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--max-evalue", type=float, default=1e-50)
    for name, _, _, task, mid, mcov in ARMS:
        ap.add_argument(f"--{name.lower()}-min-id", type=float, default=mid)
        ap.add_argument(f"--{name.lower()}-min-cov", type=float, default=mcov)
        ap.add_argument(f"--{name.lower()}-task", default=task,
                        choices=["megablast", "dc-megablast", "blastn"])
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--keep-refset-hits", action="store_true",
                    help="keep hits to accessions of the culture reference itself (default: dropped)")
    ap.add_argument("--keep-db", action="store_true", help="keep the BLAST databases")
    a = ap.parse_args()

    for tool in ("blastn", "makeblastdb"):
        if not shutil.which(tool):
            raise SystemExit(f"{tool} not found; install BLAST+ into the amfsplit env")
    man = json.load(open(os.path.join(a.ref, "MANIFEST.json")))
    all_rows = read_tsv(os.path.join(a.ref, "copies.tsv"))
    cop = {r["accession"]: r for r in all_rows if r["status"] == "retained"}
    drop = set() if a.keep_refset_hits else {base_acc(r["accession"]) for r in all_rows}
    cul_meta = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    print(f"culture_ref {man['version']} (fingerprint {man['fingerprint']}): {len(cop)} copies, "
          f"{len(cul_meta)} cultures, {sum(r['name_status'] == 'hold' for r in cul_meta.values())} on name hold")

    paths = {"maarjam": a.maarjam, "unite": a.unite, "delavaux": a.lsu}
    arms, refs_by, seqs_by, show, used = [], {}, {}, {}, {}
    print(f"\n{'arm':5s}{'reference':46s}{'seqs':>8s}{'units':>8s}{'query':>16s}{'copies':>8s}"
          f"{'ref-set seqs':>14s}{'best hit dropped':>18s}")
    for name, kind, region, _, _, _ in ARMS:
        path = paths[kind]
        if not path:
            print(f"{name:5s}(no reference given; arm skipped)")
            continue
        refs, fasta, n_all = load_reference(path, kind, a.unite_keep if kind == "unite" else None)
        qpath = os.path.join(a.ref, "fasta", f"{region}.fasta")
        if not os.path.exists(qpath):
            raise SystemExit(f"{qpath} missing from the reference")
        seqs_by[name] = read_fasta(qpath)
        hits = run_blast(qpath, fasta, a.out + ".blastdb", name.lower(), getattr(a, f"{name.lower()}_task"),
                         a.threads)
        mid, mcov = getattr(a, f"{name.lower()}_min_id"), getattr(a, f"{name.lower()}_min_cov")
        best_dropped = 0
        for acc in cop:
            hs = hits.get(acc, [])
            if hs and drop:
                top = max(hs, key=lambda h: (h["bitscore"], h["pident"]))
                best_dropped += base_acc(refs[top["ref"]]["acc"]) in drop
            call = assign(hits.get(acc, []), refs, mid, mcov, a.max_evalue, drop_accs=drop) \
                if acc in seqs_by[name] else {**{f: "" for f in FIELDS}, "call": "no_query", "dropped": 0}
            for f in FIELDS:
                cop[acc][f"{name}_{f}"] = call[f]
        refs_by[name] = refs
        arms.append(name)
        used[name] = {"path": os.path.abspath(path), "sha256": sha256(path), "records": n_all,
                      "kept": len(refs), "units": len({v["label"] for v in refs.values()}), "query": region,
                      "task": getattr(a, f"{name.lower()}_task"), "min_id": mid, "min_cov": mcov}
        n_units = len({v["label"] for v in refs.values()})
        n_refset = sum(base_acc(v["acc"]) in drop for v in refs.values())
        kept = f" [{len(refs)} of {n_all} kept]" if len(refs) < n_all else ""
        base = os.path.basename(path)
        shown = (base if len(base) + len(kept) <= 45 else base[:44 - len(kept)] + "~") + kept
        print(f"{name:5s}{shown[:45]:46s}{len(refs):>8d}{n_units:>8d}{region:>16s}{len(seqs_by[name]):>8d}"
              f"{n_refset:>14d}{best_dropped:>18d}")
        taxa = collections.defaultdict(collections.Counter)
        for v in refs.values():
            taxa[v["label"]][v["taxon"]] += 1
        names = {u: c.most_common(1)[0][0] for u, c in taxa.items()}
        show[name] = (lambda u, names=names: f"{u} ({names[u]})" if u in names else u) if name == "SH" else str
    if not arms:
        raise SystemExit("no arm could run: give at least one of --maarjam, --unite, --lsu")
    if a.keep_refset_hits:
        print("  ref-set hits kept (--keep-refset-hits): copies can take labels from themselves or culture-mates")
    else:
        print("  ref-set seqs: reference records that are copies of the culture reference. Hits to them are dropped;\n"
              "  best hit dropped: copies whose top hit was one of them (their call comes from the next hit)")
    if not a.keep_db:
        shutil.rmtree(a.out + ".blastdb", ignore_errors=True)

    # ---- culture level: modal label of main-type copies per arm
    by_cul = collections.defaultdict(list)
    for x, r in cop.items():
        by_cul[r["culture"]].append(r)
    cultures = []
    for cu, rs in sorted(by_cul.items()):
        m = cul_meta[cu]
        row = {k: m.get(k, "") for k in ("culture", "organism", "order", "genus", "name_status", "n_copies",
                                         "n_main", "n_divergent", "n_outlier", "n_libraries")}
        row["named"] = is_named(m["organism"])
        for arm in arms:
            main = [r for r in rs if r["copy_class"] == "main"]
            ok = [r for r in main if r[f"{arm}_call"] == "assigned"]
            cnt = collections.Counter(r[f"{arm}_label"] for r in ok)
            if cnt:
                top = max(cnt.values())
                tied = [u for u, n in cnt.items() if n == top]
                med = {u: statistics.median(float(r[f"{arm}_pident"]) for r in ok if r[f"{arm}_label"] == u)
                       for u in tied}
                modal = max(sorted(tied), key=lambda u: med[u])
            else:
                modal, tied = "", []
            row[f"{arm}_modal"] = modal
            row[f"{arm}_modal_tie"] = int(len(tied) > 1)
            row[f"{arm}_main_assigned"] = f"{len(ok)}/{len(main)}"
            allc = collections.Counter(r[f"{arm}_label"] if r[f"{arm}_call"] == "assigned" else f"({r[f'{arm}_call']})"
                                       for r in rs)
            row[f"{arm}_labels"] = ",".join(f"{k}:{v}" for k, v in allc.most_common())
            div = collections.Counter(r[f"{arm}_label"] if r[f"{arm}_call"] == "assigned" else f"({r[f'{arm}_call']})"
                                      for r in rs if r["copy_class"] == "divergent")
            row[f"{arm}_divergent"] = ",".join(f"{k}:{v}" for k, v in div.most_common())
        cultures.append(row)
    cul = {r["culture"]: r for r in cultures}
    for r in cop.values():
        for arm in arms:
            r[f"{arm}_culture_modal"] = cul[r["culture"]][f"{arm}_modal"]
    counted = [r for r in cultures if r["name_status"] == "ok"]
    named = [r for r in counted if r["named"]]

    # ---- 1. assignment by arm and order
    print(f"\n1. assignment: copies assigned (%) | cultures with a label (%) | median best-hit identity")
    orders = sorted({r["order"] for r in cultures})
    print(f"{'order':<17}{'cultures':>9}" + "".join(f"{arm:>26s}" for arm in arms))
    for o in orders + ["all"]:
        cs = [r for r in cultures if o == "all" or r["order"] == o]
        xs = [r for r in cop.values() if o == "all" or cul[r["culture"]]["order"] == o]
        cells = []
        for arm in arms:
            ids = [float(r[f"{arm}_best_pident"]) for r in xs if r[f"{arm}_best_pident"] != ""]
            cp = 100 * sum(r[f"{arm}_call"] == "assigned" for r in xs) / len(xs)
            cc = 100 * sum(bool(r[f"{arm}_modal"]) for r in cs) / len(cs)
            cells.append(f"{cp:5.0f} | {cc:3.0f} | {statistics.median(ids) if ids else float('nan'):5.1f}")
        print(f"{o:<17}{len(cs):>9}" + "".join(f"{c:>26s}" for c in cells))
    for arm in arms:
        c = collections.Counter(r[f"{arm}_call"] for r in cop.values())
        print(f"  {arm}: " + ", ".join(f"{k} {v}" for k, v in c.most_common()))

    # ---- 2. VT regression against Stage 1
    if "VT" in arms:
        same = sum((r["VT_call"], r["VT_label"]) == (r["vt_call"], r["vt"]) for r in cop.values())
        print(f"\n2. VT arm vs Stage 1 vt_check calls in the reference: {same}/{len(cop)} identical")
        for r in [r for r in cop.values() if (r["VT_call"], r["VT_label"]) != (r["vt_call"], r["vt"])][:8]:
            print(f"    {r['accession']}: now {r['VT_call']} {r['VT_label']}, Stage 1 {r['vt_call']} {r['vt']}")

    # ---- 3. many-to-many, counting cultures
    print(f"\n3. many-to-many, counting cultures by their modal label ({len(named)} cultures of "
          f"{len({r['organism'] for r in named})} named species; name holds and 'sp.' left out)")
    sp_units, unit_sp = {}, {}
    for arm in arms:
        su, us = collections.defaultdict(collections.Counter), collections.defaultdict(collections.Counter)
        for r in named:
            if r[f"{arm}_modal"]:
                su[r["organism"]][r[f"{arm}_modal"]] += 1
                us[r[f"{arm}_modal"]][r["organism"]] += 1
        sp_units[arm], unit_sp[arm] = su, us
        multi_s = sorted((s for s in su if len(su[s]) > 1), key=lambda s: (-len(su[s]), s))
        multi_u = sorted((u for u in us if len(us[u]) > 1), key=lambda u: (-len(us[u]), u))
        print(f"  {arm}: {len(multi_s)} of {len(su)} species span > 1 {arm}; "
              f"{len(multi_u)} of {len(us)} {arm}s hold > 1 species")
        for s in multi_s[:4]:
            print(f"      {s}: {fmt_counts(su[s], show=show[arm])}")
        for u in multi_u[:4]:
            print(f"      {show[arm](u)}: {fmt_counts(us[u])}")

    # ---- 4. within-culture splits of main-type copies
    print("\n4. within-culture splits: main-type copies of one culture assigned to two units")
    for arm in arms:
        pairs = collections.defaultdict(list)
        for cu, rs in by_cul.items():
            byu = collections.defaultdict(list)
            for r in rs:
                if r["copy_class"] == "main" and r[f"{arm}_call"] == "assigned":
                    byu[r[f"{arm}_label"]].append(r["accession"])
            for u1, u2 in itertools.combinations(sorted(byu), 2):
                s = seqs_by[arm]
                d = [dist(s[p], s[q]) for p in byu[u1] for q in byu[u2] if p in s and q in s]
                pairs[(u1, u2)].append((cu, min(d) if d else None))
        n_cul = len({cu for v in pairs.values() for cu, _ in v})
        print(f"  {arm}: {len(pairs)} unit pairs in {n_cul} cultures")
        for (u1, u2), v in sorted(pairs.items(), key=lambda kv: -len(kv[1]))[:6]:
            ds = [d for _, d in v if d is not None]
            orgs = collections.Counter(cul[cu]["organism"] for cu, _ in v)
            pair = f"{show[arm](u1)} / {show[arm](u2)}"
            print(f"      {pair}: {len(v)} culture(s) ({fmt_counts(orgs, 2)}); closest pair across the split "
                  f"{statistics.median(ds):.2f}% median, {max(ds):.2f}% max ({ARMS[[x[0] for x in ARMS].index(arm)][2]})"
                  if ds else f"      {pair}: {len(v)} culture(s)")

    # ---- 5. cross-arm splits
    if len(arms) > 1:
        print("\n5. cross-arm: units of one arm whose cultures carry several modal labels of another")
        for a1, a2 in itertools.permutations(arms, 2):
            m = collections.defaultdict(collections.Counter)
            for r in counted:
                if r[f"{a1}_modal"] and r[f"{a2}_modal"]:
                    m[r[f"{a1}_modal"]][r[f"{a2}_modal"]] += 1
            split = sorted((u for u in m if len(m[u]) > 1), key=lambda u: (-len(m[u]), u))
            print(f"  {a1} split by {a2}: {len(split)} of {len(m)} {a1}s")
            for u in split[:3]:
                print(f"      {show[a1](u)}: {fmt_counts(m[u], 4, show[a2])}")

    # ---- 6. divergent-class and SSU-outlier copies
    print("\n6. divergent-class and SSU-outlier copies against their culture's main-type unit")
    groups6 = [("divergent", None, "divergent-class"), ("outlier", None, "SSU-outlier")]
    if any(r.get("div_kind") for r in cop.values()):  # references frozen with 5.8S kinds
        groups6[:1] = [("divergent", "short_5.8S", "divergent-class (short 5.8S)"),
                       ("divergent", "same_length", "divergent 5.8S, same length (substitutions)"),
                       ("divergent", "longer", "divergent 5.8S, longer than the main type"),
                       ("divergent", "", "divergent-class (no 5.8S length)")]
    for cls, kind, what in groups6:
        xs = [r for r in cop.values() if r["copy_class"] == cls and (kind is None or r.get("div_kind", "") == kind)]
        if not xs:
            if kind is None:
                print(f"  no {what} copies")
            continue
        print(f"  {len(xs)} {what} copies in {len({r['culture'] for r in xs})} cultures")
        print(f"    {'arm':6s}{'same unit':>11s}{'other unit':>12s}{'no call':>9s}{'no culture unit':>17s}"
              f"   median best-hit identity ({cls} vs main-type)")
        for arm in arms:
            same = other = none = nomodal = 0
            for r in xs:
                m = r[f"{arm}_culture_modal"]
                if not m:
                    nomodal += 1
                elif r[f"{arm}_call"] != "assigned":
                    none += 1
                elif r[f"{arm}_label"] == m:
                    same += 1
                else:
                    other += 1
            idd = [float(r[f"{arm}_best_pident"]) for r in xs if r[f"{arm}_best_pident"] != ""]
            idm = [float(r[f"{arm}_best_pident"]) for r in cop.values()
                   if r["copy_class"] == "main" and r[f"{arm}_best_pident"] != ""]
            print(f"    {arm:6s}{same:>11d}{other:>12d}{none:>9d}{nomodal:>17d}   "
                  f"{statistics.median(idd) if idd else float('nan'):.1f} vs "
                  f"{statistics.median(idm) if idm else float('nan'):.1f}")

    # ---- 7. genus agreement
    gen_arms = [x for x in arms if x in ("SH", "LSU")]
    if gen_arms:
        print("\n7. culture genus vs genus of the best reference hit (main-type modal hit, counting cultures)")
        syn = lambda g: GENUS_SYNONYMS.get(g, g)
        for arm in gen_arms:
            agree, dis, via, nogenus = 0, collections.Counter(), collections.Counter(), 0
            for cu, r in cul.items():
                if r["name_status"] != "ok" or not r[f"{arm}_modal"]:
                    continue
                hits = [x for x in by_cul[cu] if x["copy_class"] == "main" and x[f"{arm}_call"] == "assigned"
                        and x[f"{arm}_label"] == r[f"{arm}_modal"]]
                g = collections.Counter(x[f"{arm}_best_genus"] for x in hits).most_common(1)[0][0]
                if not g:  # reference carries no genus (unidentified, gen_Incertae_sedis)
                    nogenus += 1
                elif g == r["genus"]:
                    agree += 1
                elif g and syn(g) == syn(r["genus"]):
                    agree += 1
                    via[f"{r['genus']} = {g}"] += 1
                else:
                    dis[f"{r['genus']} -> {g or '?'}"] += 1
            tot = agree + sum(dis.values())
            print(f"  {arm}: {agree}/{tot} agree" + (f" ({fmt_counts(via)} as synonyms)" if via else "")
                  + (f"; differ: {fmt_counts(dis, 6)}" if dis else "")
                  + (f"; {nogenus} more where the reference has no genus" if nogenus else ""))

    # ---- tables
    ccols = (["accession", "culture", "organism", "order", "copy_class", "flags", "name_status"]
             + [f"{arm}_{f}" for arm in arms for f in FIELDS + ["culture_modal"]])
    write_tsv(f"{a.out}.copies.tsv", list(cop.values()), ccols)
    write_tsv(f"{a.out}.cultures.tsv", cultures, list(cultures[0]))
    species = []
    for org in sorted({r["organism"] for r in counted}):
        rs = [r for r in counted if r["organism"] == org]
        row = {"organism": org, "named": is_named(org), "order": rs[0]["order"], "cultures": len(rs),
               "copies": sum(int(r["n_copies"]) for r in rs)}
        for arm in arms:
            c = collections.Counter(r[f"{arm}_modal"] or "(none)" for r in rs)
            row[f"{arm}_units"] = len([k for k in c if k != "(none)"])
            row[f"{arm}_cultures"] = fmt_counts(c, show=show[arm])
        species.append(row)
    write_tsv(f"{a.out}.species.tsv", species, list(species[0]))
    units = []
    for arm in arms:
        taxa = collections.defaultdict(collections.Counter)
        for v in refs_by[arm].values():
            taxa[v["label"]][v["taxon"]] += 1
        for u in sorted({r[f"{arm}_modal"] for r in counted if r[f"{arm}_modal"]}):
            rs = [r for r in counted if r[f"{arm}_modal"] == u]
            row = {"arm": arm, "unit": u, "reference_taxa": fmt_counts(taxa[u], 3), "cultures": len(rs),
                   "named_species": fmt_counts(collections.Counter(r["organism"] for r in rs if r["named"])),
                   "unnamed": fmt_counts(collections.Counter(r["organism"] for r in rs if not r["named"]))}
            for other in arms:
                if other != arm:
                    row[f"{other}_modal"] = fmt_counts(collections.Counter(r[f"{other}_modal"] or "(none)" for r in rs),
                                                       show=show[other])
            units.append(row)
    ucols = ["arm", "unit", "reference_taxa", "cultures", "named_species", "unnamed"] + [f"{x}_modal" for x in arms]
    write_tsv(f"{a.out}.units.tsv", units, ucols)
    with open(f"{a.out}.run.json", "w") as fh:
        json.dump({"culture_ref": {"path": os.path.abspath(a.ref), "version": man["version"],
                                   "fingerprint": man["fingerprint"]},
                   "arms": used, "max_evalue": a.max_evalue, "refset_hits": "kept" if a.keep_refset_hits else "dropped",
                   "script_sha256": sha256(os.path.abspath(__file__)), "args": vars(a)}, fh, indent=1)
        fh.write("\n")
    print(f"\ntables: {a.out}.copies.tsv, .cultures.tsv, .species.tsv, .units.tsv, .run.json "
          f"(reference {man['version']}, fingerprint {man['fingerprint']})")


if __name__ == "__main__":
    sys.exit(main())
