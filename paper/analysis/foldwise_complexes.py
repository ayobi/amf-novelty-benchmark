#!/usr/bin/env python3
"""The complex rule applied prospectively inside each round of the LSU-removal test (added after review, point 1.1).

The paper's complexes (species_groups.py) were derived from all cultures, the withheld species included. Here the
rule is applied as it would be to a species nobody has cultured: in each species round of the LSU-removal test, the
withheld cultures (holdout_plan.json, 'held_cultures') are removed, the within-species spread T is recomputed on the
remaining cultures only (99th percentile, statistics.quantiles' exclusive method, as in species_groups.py), and each
main-type copy of the withheld species is tested against every remaining named species. A copy is INSIDE a species
when, on all three markers (full ITS, LROR-FLR2, SSU flank), its distance to the nearest main-type copy of that
species is at most T, the rule used to lock the LongGloDB and EUKARYOME predictions. Under the paper's claim, a copy
can be called a known species only if it lies inside one; a false-known call outside that set contradicts it.

The withheld species' own copies are used only as queries, as a field read would be; they inform neither T nor the
reference. Nothing here reads a test result. Written after the results were known, so it is a retrospective
analysis of a prospective rule, not a locked prediction.

Species level (always): observed false-known counts per round are read from tables/tableS1_m4_folds.tex
(make_numbers.py writes it from the M4 aggregate). Copy level (with --m4 <aggregate_v2>): every false-known call is
matched to its copy and the species it was called.

  python3 analysis/foldwise_complexes.py --ref ../reference/culture_ref_v4 \\
      --plan ../benchmarks/m4_external_holdout/holdout_plan.json --out analysis/foldwise [--m4 <aggregate_v2>] --paper .
Writes {out}.limits.tsv (T per round), {out}.copies.tsv (one row per withheld main-type copy), {out}.species.tsv
(one row per round), numbers_foldwise.tex.
"""
import argparse
import collections
import csv
import json
import os
import re
import statistics
import sys
import time

import edlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'stage3'))
from arms import is_named, read_fasta, read_tsv, write_tsv  # noqa: E402

MARKERS = [("full_its", "ITS"), ("LSU_LROR-FLR2", "LSU"), ("ssu_flank", "SSU")]
DEV = {"Glomus chinense", "Glomus rugosae"}


def percentile(vals, p):
    d = sorted(vals)
    n = len(d)
    h = (n + 1) * p / 100.0
    j = min(max(int(h), 1), n - 1)
    t = d[j - 1] + (d[j] - d[j - 1]) * (h - j)
    if float(p).is_integer() and abs(t - statistics.quantiles(vals, n=100)[int(p) - 1]) > 1e-9:
        raise AssertionError("percentile() departs from statistics.quantiles")
    return t


def dist(a, b, cap=None):
    """edlib infix distance as % of the shorter sequence; None when above cap (%)."""
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    k = -1 if cap is None else int(cap * len(q) / 100) + 1
    d = edlib.align(q, t, mode="HW", task="distance", k=k)["editDistance"]
    if d < 0:
        return None
    d = 100 * d / len(q)
    return None if cap is not None and d > cap else d


def table_s1(path):
    """Per-species false-known counts from the generated supplementary table."""
    out = {}
    for line in open(path):
        m = re.match(r"\\emph\{(.+?)\}(\$\^\{\\(dagger|ddagger)\}\$)? & (\d+) & (\d+) & (\d+) & (\d+) & (\d+)\\\\", line.strip())
        if m:
            sp = m.group(1).replace("{[}", "[").replace("{]}", "]")
            out[sp] = {"mark": m.group(3) or "", "copies": int(m.group(4)), "ctl_v1": int(m.group(5)),
                       "ctl_v2": int(m.group(6)), "v1": int(m.group(7)), "v2": int(m.group(8))}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--paper", default=".")
    ap.add_argument("--m4", help="aggregate_v2 directory of the LSU-removal run (species.{v1,v2,control_v1,control_v2}.copies.tsv)")
    ap.add_argument("--quantile", type=float, default=99.0)
    a = ap.parse_args()
    t0 = time.time()
    cul = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    sp_of = {c: m["organism"] for c, m in cul.items() if is_named(m["organism"]) and m.get("name_status", "ok") == "ok"}
    copies = [r for r in read_tsv(os.path.join(a.ref, "copies.tsv")) if r["status"] == "retained"]
    main_named = [r for r in copies if r["copy_class"] == "main" and r["culture"] in sp_of]
    seq = {m: read_fasta(os.path.join(a.ref, "fasta", f"{m}.fasta")) for m, _ in MARKERS}
    species = sorted({sp_of[r["culture"]] for r in main_named})
    items = {}                                  # (marker, species) -> sorted unique (sequence, culture)
    for m, _ in MARKERS:
        for sp in species:
            items[(m, sp)] = sorted({(seq[m][r["accession"]], r["culture"]) for r in main_named
                                     if sp_of[r["culture"]] == sp and r["accession"] in seq[m]})
    # all within-species distances once; each round reuses them on its remaining cultures
    pair = {}
    for m, _ in MARKERS:
        for sp in species:
            it = items[(m, sp)]
            for i, (s1, c1) in enumerate(it):
                for s2, c2 in it[i + 1:]:
                    if c1 != c2 and (s1, s2) not in pair:
                        pair[(s1, s2)] = pair[(s2, s1)] = dist(s1, s2)

    def spread(drop):
        T = {}
        for m, lab in MARKERS:
            vals = []
            for sp in species:
                it = [x for x in items[(m, sp)] if x[1] not in drop]
                if len({c for _, c in it}) < 2:
                    continue
                for s1, c1 in it:
                    vals.append(min(0.0 if s1 == s2 else pair[(s1, s2)] for s2, c2 in it if c2 != c1))
            T[lab] = percentile(vals, a.quantile)
        return T

    T_all = spread(frozenset())
    ref_T = {r["marker"]: float(r["T"]) for r in read_tsv(os.path.join(HERE, "v4.limits.tsv"))}
    for lab in T_all:
        if abs(T_all[lab] - ref_T[lab]) > 1e-3:
            raise AssertionError(f"T with nothing withheld departs from v4.limits.tsv: {lab} {T_all[lab]} vs {ref_T[lab]}")
    print("T, nothing withheld: " + ", ".join(f"{k} {v:.4f}%" for k, v in T_all.items()) + " (equals v4.limits.tsv)")
    complexes = {}
    for r in read_tsv(os.path.join(HERE, "v4.complexes.tsv")):
        for sp in r["species"].split(";"):
            complexes[sp] = r["complex"]
    plan = [p for p in json.load(open(a.plan))["folds"] if p["rank"] == "species"]
    acc = {r["accession"]: r for r in copies}
    lim_rows, copy_rows, sp_rows = [], [], []
    obs = table_s1(os.path.join(a.paper, "tables", "tableS1_m4_folds.tex"))
    for p in plan:
        target, drop = p["target"], frozenset(p["held_cultures"])
        T = spread(drop)
        lim_rows.append({"target": target, **{lab: round(T[lab], 4) for lab in T}})
        remaining = [sp for sp in species if any(x[1] not in drop for m, _ in MARKERS for x in items[(m, sp)])]
        queries = [x for x in p["query_ids"] if acc[x]["copy_class"] == "main"]
        for q in queries:
            d_in = {}
            for sp in remaining:
                ok = True
                for m, lab in MARKERS:
                    s = seq[m].get(q)
                    refs = [x[0] for x in items[(m, sp)] if x[1] not in drop]
                    ds = [dist(s, r, T[lab]) for r in refs] if s and refs else []
                    ds = [d for d in ds if d is not None]
                    if not ds:
                        ok = False
                        break
                    d_in[(sp, lab)] = min(ds)
                if not ok:
                    for _, lab in MARKERS:
                        d_in.pop((sp, lab), None)
            inside = sorted({sp for sp, _ in d_in})
            copy_rows.append({"target": target, "accession": q, "culture": acc[q]["culture"],
                              "complex_retrospective": complexes.get(target, ""),
                              "inside_prospective": ";".join(inside),
                              "inside_detail": ";".join(f"{sp} (" + ", ".join(f"{lab} {d_in[(sp, lab)]:.2f}" for _, lab in MARKERS) + ")"
                                                        for sp in inside)})
        rs = [r for r in copy_rows if r["target"] == target]
        o = obs.get(target, {})
        n_in = sum(bool(r["inside_prospective"]) for r in rs)
        role = "development" if target in DEV else ("not in V18" if p.get("external_status") != "removed" else "test")
        sp_rows.append({"target": target, "role": role, "main_copies": len(rs),
                        "complex_retrospective": complexes.get(target, ""),
                        "inside_prospective": n_in,
                        "inside_species": "; ".join(f"{s} ({n})" for s, n in collections.Counter(
                            s for r in rs for s in r["inside_prospective"].split(";") if s).most_common()),
                        "fk_intact_v1": o.get("ctl_v1", ""), "fk_intact_v2": o.get("ctl_v2", ""),
                        "fk_removed_v1": o.get("v1", ""), "fk_removed_v2": o.get("v2", "")})
        if o and o["copies"] != len(rs):
            raise AssertionError(f"{target}: {len(rs)} main-type queries but table S1 has {o['copies']} copies")
        print(f"  {target:32s} T ITS {T['ITS']:.2f} LSU {T['LSU']:.2f} SSU {T['SSU']:.2f}; inside {n_in:3d}/{len(rs):3d}; "
              f"FK removed v1 {o.get('v1', '-')}, v2 {o.get('v2', '-')}", flush=True)
    write_tsv(f"{a.out}.limits.tsv", lim_rows, list(lim_rows[0]))
    write_tsv(f"{a.out}.copies.tsv", copy_rows, list(copy_rows[0]))
    write_tsv(f"{a.out}.species.tsv", sp_rows, list(sp_rows[0]))

    # ---- summary over the test species
    test = [r for r in sp_rows if r["role"] == "test"]
    S = lambda k, rows=test: sum(int(r[k]) for r in rows)  # noqa: E731
    mac = {}
    mac["FwTestSpecies"] = len(test)
    mac["FwCopies"] = S("main_copies")
    mac["FwInsideCopies"] = S("inside_prospective")
    mac["FwInsideSpecies"] = sum(int(r["inside_prospective"]) > 0 for r in test)
    mac["FwRetroCopies"] = sum(int(r["main_copies"]) for r in test if r["complex_retrospective"])
    mac["FwRetroSpecies"] = sum(bool(r["complex_retrospective"]) for r in test)
    for k, tag in (("fk_removed_v2", "VTwo"), ("fk_removed_v1", "VOne"), ("fk_intact_v2", "CtlTwo"), ("fk_intact_v1", "CtlOne")):
        mac[f"FwFK{tag}"] = S(k)
        mac[f"FwFKSpecies{tag}"] = sum(int(r[k]) > 0 for r in test)
        # copies that cannot all be inside: false known in excess of the copies predicted inside, per species
        mac[f"FwExcess{tag}"] = sum(max(0, int(r[k]) - int(r["inside_prospective"])) for r in test)
        mac[f"FwFKSpNotInside{tag}"] = sum(int(r[k]) > 0 and int(r["inside_prospective"]) == 0 for r in test)
    mac["FwInsideNoFK"] = sum(int(r["inside_prospective"]) > 0 and int(r["fk_removed_v2"]) == 0 for r in test)
    mac["FwNoInsideSpecies"] = sum(int(r["inside_prospective"]) == 0 for r in test)
    mac["FwMaxAnticipated"] = sum(min(int(r["fk_removed_v2"]), int(r["inside_prospective"])) for r in test)
    # the species whose false-known calls outnumber its copies inside: which marker keeps the rest outside
    ex = [r for r in test if int(r["fk_removed_v2"]) > int(r["inside_prospective"])]
    if len(ex) > 1:
        raise AssertionError("more than one species with errors beyond the prediction; the text names one")
    if ex:
        r = ex[0]
        partner = r["inside_species"].split(" (")[0]
        p = next(p for p in plan if p["target"] == r["target"])
        drop = frozenset(p["held_cultures"])
        T = next(x for x in lim_rows if x["target"] == r["target"])
        refs = {lab: [x[0] for x in items[(m, partner)] if x[1] not in drop] for m, lab in MARKERS}
        fails = collections.Counter()
        far = collections.defaultdict(list)
        for c in [c for c in copy_rows if c["target"] == r["target"] and partner not in c["inside_prospective"].split(";")]:
            d = {lab: min(dist(seq[m][c["accession"]], s) for s in refs[lab]) for m, lab in MARKERS}
            bad = tuple(lab for _, lab in MARKERS if d[lab] > T[lab])
            fails[bad] += 1
            for lab in bad:
                far[lab].append(d[lab])
        main_lab = max(far, key=lambda k: len(far[k]))
        def it(sp):
            return "\\spn{" + sp + "}"
        mac.update({"FwExSpecies": it(r["target"]), "FwExPartner": it(partner), "FwExInside": int(r["inside_prospective"]),
                    "FwExCopies": int(r["main_copies"]), "FwExFK": int(r["fk_removed_v2"]),
                    "FwExMarker": {"LSU": "LROR--FLR2", "ITS": "the full ITS", "SSU": "the SSU flank"}[main_lab],
                    "FwExFarN": len(far[main_lab]), "FwExRange": f"{min(far[main_lab]):.1f}--{max(far[main_lab]):.1f}\\%",
                    "FwExT": f"{T[main_lab]:.2f}\\%", "FwExOtherN": sum(n for k, n in fails.items() if main_lab not in k)})
        print("outside", r["target"], "vs", partner, dict(fails))
    lt = [r for r in lim_rows if any(s["target"] == r["target"] and s["role"] == "test" for s in sp_rows)]
    for lab in ("ITS", "LSU", "SSU"):
        v = [r[lab] for r in lt]
        mac[f"FwT{lab}Range"] = f"{min(v):.2f}--{max(v):.2f}\\%"
    copy_level = {}
    if a.m4:
        cx_of = {r["accession"]: r for r in copy_rows}
        dev_rows = {}
        for key, tag in (("v2", "VTwo"), ("v1", "VOne"), ("control_v2", "CtlTwo"), ("control_v1", "CtlOne")):
            rows = [r for r in read_tsv(os.path.join(a.m4, f"species.{key}.copies.tsv"))
                    if r["copy_class"] == "main" and r["group"] == "named, species held out"
                    and r["development_fold"] == "0" and r["external_status"] == "removed"]
            fk = [r for r in rows if r["status"].startswith("known")]
            inside = [r for r in fk if cx_of[r["accession"]]["inside_prospective"]]
            same = [r for r in fk if set(p for p in r["placement"].split(";") if p)
                    & set(cx_of[r["accession"]]["inside_prospective"].split(";"))]
            copy_level[key] = {"false_known": len(fk), "inside_any": len(inside), "called_species_inside": len(same),
                               "outside": [r["accession"] + " " + r["organism"] + " -> " + r["placement"] for r in fk
                                           if not cx_of[r["accession"]]["inside_prospective"]]}
            mac[f"FwCopyInside{tag}"] = len(inside)
            mac[f"FwCopyOutside{tag}"] = len(fk) - len(inside)
            mac[f"FwCopySame{tag}"] = len(same)
            if key == "v2" and len(fk) != mac["FwFKVTwo"]:
                raise AssertionError("aggregate and table S1 disagree on false-known calls")
        json.dump(copy_level, open(f"{a.out}.copy_level.json", "w"), indent=1)
    # ---- supplementary table: rounds where anything lies inside or anything was called known
    def tex(sp):
        return "\\emph{" + sp.replace("[", "{[}").replace("]", "{]}") + "}"
    def short(sp):                      # epithet only: every species inside is a congener of the withheld one
        return "\\emph{" + sp.split(" ", 1)[1] + "}"
    body, quiet = [], [r for r in test if not int(r["inside_prospective"]) and not int(r["fk_removed_v2"]) and not int(r["fk_removed_v1"])]
    for r in sorted(sp_rows, key=lambda r: (r["role"] != "test", r["target"])):
        if r in quiet or r["role"] == "not in V18":
            continue
        ins = "; ".join(f"{short(s.rsplit(' (', 1)[0])} ({s.rsplit(' (', 1)[1]}" for s in r["inside_species"].split("; ") if s) or "--"
        mark = "$^{\\dagger}$" if r["role"] == "development" else ""
        body.append(f"{tex(r['target'])}{mark} & {r['main_copies']} & {r['complex_retrospective'] or '--'} & {r['inside_prospective']} & "
                    f"{ins} & {r['fk_removed_v1']} & {r['fk_removed_v2']}\\\\")
    body.append(f"{len(quiet)} other test species & {sum(int(r['main_copies']) for r in quiet)} & -- & 0 & -- & 0 & 0\\\\")
    with open(os.path.join(a.paper, "tables", "table_foldwise.tex"), "w") as fh:
        fh.write(r"""\begin{table}[htbp]\centering
\begin{threeparttable}
\caption{The complex rule applied inside each round of the LSU-removal test.}
\label{tab:foldwise}
\footnotesize\setlength{\tabcolsep}{3pt}
\begin{tabular}{lrcrlrr}
\toprule
Withheld species & Copies & Complex & \multicolumn{2}{c}{Inside a remaining congener} & \multicolumn{2}{c}{False known}\\
 & & (all cultures) & copies & epithet (copies) & no rule & ITS rule\\
\midrule
""" + "\n".join(body) + r"""
\bottomrule
\end{tabular}
\begin{tablenotes}\footnotesize
\item Main-type copies of each withheld species. \emph{Complex}: the species' complex when the rule is computed from all
cultures (Table~\ref{tab:complexes}). \emph{Inside}: with the withheld cultures removed and the within-species spread
recomputed without them, copies that lie within the spread of a remaining species on all three markers; only these
can be called a known species if the rule holds. False-known calls with the species removed from the LSU reference
(Table~\ref{tab:m4folds}). $^{\dagger}$\,Development species.
\end{tablenotes}
\end{threeparttable}
\end{table}
""")
    with open(os.path.join(a.paper, "numbers_foldwise.tex"), "w") as fh:
        fh.write("% written by analysis/foldwise_complexes.py; do not edit by hand\n")
        for k, v in mac.items():
            v = f"{v:,}" if isinstance(v, int) else v
            fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
    print(json.dumps(mac, indent=1))
    if copy_level:
        print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "outside"} | {"n_outside": len(v["outside"])}
                          for k, v in copy_level.items()}, indent=1))
    print(f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    sys.exit(main())
