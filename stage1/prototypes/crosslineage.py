#!/usr/bin/env python3
"""crosslineage.py: how the Stage 1 region model holds up across AMF lineages.

Input: GenBank rDNA copies (e.g. Stefani et al. 2025, PX207096-PX207274 and
PX214440-PX215164) run through run_stage1.sh, plus the original FASTA for the
organism/isolate in each header.

Groups come from gb_groups.py (--groups), at culture level (strain, else shared isolate
label) or library level (isolate minus the per-copy _s<N> / _s<N>_<M> suffix). Without
--groups they are parsed from the FASTA headers, which pools every copy whose header
lacks an isolate by species; keep that only for old runs.

Per group and per order it reports:
  * anchored: share of copies with a full ITSxRust chain
  * AMF%: share of copies triaged as AMF (all inputs are AMF, so the rest are false negatives)
  * window recovery: share of copies with a usable window (ok or open-ended)
  * within-group divergence: mean and max pairwise edit distance / mean length,
    per region, between copies of the same group (the no-barcode-gap problem,
    measured with this region model)

Usage:
  crosslineage.py --fasta stefani_rdna.fasta --prefix out/windows/stefani \
      --groups stefani_groups.tsv --level culture --out crosslineage.culture
"""
import argparse
import csv
import itertools
import re
import statistics
from collections import defaultdict

import edlib

ORDER = {  # genus -> order (Glomeromycota); extend as needed
    **dict.fromkeys(["Rhizophagus", "Glomus", "Funneliformis", "Septoglomus", "Sclerocystis",
                     "Dominikia", "Microdominikia", "Kamienskia", "Nanoglomus", "Oehlia",
                     "Rhizoglomus", "Simiglomus", "Halonatospora", "Microkamienskia",
                     "Epigeocarpum", "Albahypha", "Orientoglomus"], "Glomerales"),
    **dict.fromkeys(["Claroideoglomus", "Entrophospora", "Viscospora"], "Entrophosporales"),
    **dict.fromkeys(["Diversispora", "Redeckera", "Otospora", "Acaulospora", "Kuklospora",
                     "Pacispora", "Gigaspora", "Scutellospora", "Racocetra", "Cetraspora",
                     "Dentiscutata", "Fuscutata", "Quatunica", "Intraornatospora",
                     "Paradentiscutata", "Bulbospora", "Corymbiglomus", "Sieverdingia"],
                    "Diversisporales"),
    **dict.fromkeys(["Archaeospora", "Ambispora", "Geosiphon", "Palaeospora"], "Archaeosporales"),
    **dict.fromkeys(["Paraglomus", "Innospora"], "Paraglomerales"),
}
REGIONS = ["ssu_flank", "its1", "s58", "its2", "SSU_VT_NS31-AM1", "SSU_V4_AMV4.5NF-AMDGR",
           "LSU_LROR-FLR2", "LSU_D2_FLR3-FLR4"]
WINDOWS = ["SSU_VT_NS31-AM1", "SSU_V4_AMV4.5NF-AMDGR", "LSU_LROR-FLR2", "LSU_D2_FLR3-FLR4"]
USABLE = {"ok", "open_5p", "open_3p"}
CUT = re.compile(r"\s(?:clone|small subunit|18S|internal transcribed|5\.8S|large subunit|28S|"
                 r"genomic|ribosomal|partial|complete|rRNA)\b.*$", re.I)


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


def culture_of(header):
    """'PX1.1 Rhizophagus irregularis isolate DAOM197198 clone 3 18S ...' -> (culture, genus)"""
    desc = header.split(None, 1)[1] if " " in header else header
    desc = CUT.sub("", desc).strip().rstrip(",;")
    genus = desc.split()[0].strip("[]") if desc else "?"  # '[Rhizoglomus]' -> 'Rhizoglomus'
    return desc or "?", genus


def pair_distance(a, b, metric):
    """Pairwise distance in %.

    nw:    global edit distance / mean length; any length difference counts, including
           truncated read ends and a few bases of region-boundary jitter
    infix: the shorter sequence aligned inside the longer (end gaps free) / shorter length;
           ignores terminal overhangs, so it is a lower bound where real terminal indels exist
    """
    if metric == "nw":
        ed = edlib.align(a, b, mode="NW", task="distance")["editDistance"]
        return 100 * ed / ((len(a) + len(b)) / 2)
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    ed = edlib.align(q, t, mode="HW", task="distance")["editDistance"]
    return 100 * ed / len(q)


def divergence(seqs, metric="infix"):
    """Mean and max pairwise distance (%), capped at 400 pairs."""
    pairs = [(a, b) for a, b in itertools.combinations(seqs, 2) if a and b][:400]
    d = [pair_distance(a, b, metric) for a, b in pairs]
    return (round(statistics.mean(d), 2), round(max(d), 2)) if d else (None, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="GenBank FASTA that was run through run_stage1.sh")
    ap.add_argument("--prefix", required=True, help="amfsplit output prefix (out/windows/<sample>)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--groups", help="gb_groups.py TSV (accession, genus, library, culture)")
    ap.add_argument("--level", choices=["culture", "library"], default="culture",
                    help="grouping level when --groups is given (default culture)")
    ap.add_argument("--min-copies", type=int, default=3, help="copies per group for divergence (default 3)")
    ap.add_argument("--metric", choices=["infix", "nw"], default="infix",
                    help="pairwise distance: infix (end gaps free, default) or nw (global)")
    ap.add_argument("--focus", default="Rhizophagus irregularis",
                    help="print per-group divergence for groups whose name starts with this ('' to skip)")
    args = ap.parse_args()

    headers = {}
    for line in open(args.fasta):
        if line.startswith(">"):
            headers[line[1:].split()[0]] = line[1:].strip()
    culture, genus = {}, {}
    if args.groups:
        with open(args.groups) as fh:
            grp = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
        missing = [a for a in headers if not grp.get(a, {}).get(args.level)]
        if missing:
            raise SystemExit(f"{len(missing)} accessions lack a {args.level} label in {args.groups}, "
                             f"e.g. {missing[:3]}")
        for acc in headers:
            culture[acc], genus[acc] = grp[acc][args.level], grp[acc]["genus"]
        level = args.level
    else:
        for acc, h in headers.items():
            culture[acc], genus[acc] = culture_of(h)
        level = "header"

    records = {}
    with open(f"{args.prefix}.records.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            records[row["read_id"]] = row
    triage = {}
    try:
        with open(f"{args.prefix}.triage.tsv") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                triage[row["read_id"]] = row["call"]
    except FileNotFoundError:
        pass
    region_seqs = {}
    for r in REGIONS:
        try:
            region_seqs[r] = read_fasta(f"{args.prefix}.{r}.fasta")
        except FileNotFoundError:
            region_seqs[r] = {}

    by_culture = defaultdict(list)
    for acc in headers:
        by_culture[culture[acc]].append(acc)

    rows = []
    for cult, accs in sorted(by_culture.items()):
        g = genus[accs[0]]
        row = {"culture": cult, "genus": g, "order": ORDER.get(g, "unassigned"), "copies": len(accs),
               "anchored_pct": round(100 * sum(a in records for a in accs) / len(accs), 1)}
        row["triage_AMF_pct"] = (round(100 * sum(triage.get(a) == "AMF" for a in accs) / len(accs), 1)
                                 if triage else None)
        for w in WINDOWS:
            ok = sum(records.get(a, {}).get(f"{w}_status") in USABLE for a in accs)
            row[f"{w}_pct"] = round(100 * ok / len(accs), 1)
        for r in REGIONS:
            # primer-cut windows: complete ones only, since open-ended windows run to the read end
            seqs = [region_seqs[r][a] for a in accs if a in region_seqs[r]
                    and (r not in WINDOWS or records.get(a, {}).get(f"{r}_status") == "ok")]
            mean_d, max_d = divergence(seqs, args.metric) if len(seqs) >= args.min_copies else (None, None)
            row[f"{r}_div_n"] = len(seqs)
            row[f"{r}_div_mean"], row[f"{r}_div_max"] = mean_d, max_d
        rows.append(row)

    cols = list(rows[0].keys()) if rows else []
    with open(f"{args.out}.cultures.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # order-level summary (copy-weighted recovery; median of culture divergences)
    by_order = defaultdict(list)
    for r in rows:
        by_order[r["order"]].append(r)
    print(f"{len(headers)} copies, {len(rows)} groups ({level} level)\n")
    head = f"{'order':<17}{'cult':>5}{'copies':>7}{'anch%':>7}{'AMF%':>7}" + "".join(
        f"{w.split('_', 1)[1][:13]:>15}" for w in WINDOWS)
    print("window recovery (% of copies)\n" + head)
    for o, rs in sorted(by_order.items()):
        n = sum(r["copies"] for r in rs)
        wavg = lambda k: sum(r[k] * r["copies"] for r in rs) / n
        amf = f"{wavg('triage_AMF_pct'):>7.1f}" if rs[0]["triage_AMF_pct"] is not None else f"{'-':>7}"
        print(f"{o:<17}{len(rs):>5}{n:>7}{wavg('anchored_pct'):>7.1f}{amf}" +
              "".join(f"{wavg(f'{w}_pct'):>15.1f}" for w in WINDOWS))
    print(f"\nwithin-group divergence ({level} level, {args.metric} distance; complete windows only), median over "
          f"groups of mean (max over groups of max) pairwise % per region; n = groups with >= {args.min_copies} copies")
    print(f"{'order':<17}{'n':>4}" + "".join(f"{r[:13]:>16}" for r in REGIONS))
    for o, rs in sorted(by_order.items()):
        cells = []
        for r in REGIONS:
            means = [x[f"{r}_div_mean"] for x in rs if x[f"{r}_div_mean"] is not None]
            maxes = [x[f"{r}_div_max"] for x in rs if x[f"{r}_div_max"] is not None]
            cells.append(f"{statistics.median(means):.1f} ({max(maxes):.1f})" if means else "-")
        n_div = sum(x["copies"] >= args.min_copies for x in rs)
        print(f"{o:<17}{n_div:>4}" + "".join(f"{c:>16}" for c in cells))

    focus = [r for r in rows if args.focus and r["culture"].startswith(args.focus)]
    if focus:
        show = ["ssu_flank", "its1", "s58", "its2", "LSU_LROR-FLR2"]
        print(f"\n{args.focus}, per {level}: mean (max) pairwise %")
        print(f"{'group':<40}{'n':>4}" + "".join(f"{s[:13]:>16}" for s in show))
        for r in focus:
            name = r["culture"].split(" | ")[-1] if " | " in r["culture"] else r["culture"][len(args.focus):].strip()
            cells = [f"{r[f'{s}_div_mean']} ({r[f'{s}_div_max']})" if r[f"{s}_div_mean"] is not None else "-"
                     for s in show]
            print(f"{(name or '(species pool)')[:39]:<40}{r['copies']:>4}" + "".join(f"{c:>16}" for c in cells))
        # identical sequences shared across groups: one variant catalog for the species, or per-group copy sets?
        print(f"\n{args.focus}: identical sequences across {level}s")
        for s in ["its1", "its2", "LSU_LROR-FLR2"]:
            var = defaultdict(set)
            for acc in headers:
                if culture[acc].startswith(args.focus) and acc in region_seqs[s] and (
                        s not in WINDOWS or records.get(acc, {}).get(f"{s}_status") == "ok"):
                    var[region_seqs[s][acc]].add(culture[acc])
            n_copies = sum(1 for acc in headers if culture[acc].startswith(args.focus) and acc in region_seqs[s])
            shared = [g for g in var.values() if len(g) > 1]
            print(f"  {s:<14} {n_copies:>4} copies -> {len(var):>3} distinct; {len(shared)} shared by >1 {level} "
                  f"(widest in {max(map(len, var.values()), default=0)} {level}s)")
    print(f"\nper-group table: {args.out}.cultures.tsv")


if __name__ == "__main__":
    main()
