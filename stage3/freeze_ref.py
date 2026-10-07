#!/usr/bin/env python3
"""Phase 0: freeze a versioned culture reference from the Stage 1 outputs.

The reference is a function of the Stage 1 tables plus an edits file, so every change is written
down, checked against the data before it is applied, and recorded in the manifest.

Edits TSV (lines starting with # are comments), columns: action, target, value, expect, reason
  drop          <accessions>                         copies leave the reference
  relabel       <accessions>   <organism>            copies move to another organism (a new culture
                                                     if their culture holds two lineages)
  copy_class    <accessions>   main|divergent|outlier  sets the copy class (overrides the Stage 1 flags)
  name_culture  <culture label> <organism>           every copy of the culture gets this organism; a
                                                     label filed under several organisms becomes one culture
  check_name    <culture label> <species,species>    culture is held out of species-level tallies and
                                                     its ITS is compared against each candidate species
<accessions> is one accession or a comma-separated list. Edits apply in the order drop, name_culture,
relabel, copy_class, check_name: culture-wide names first, then copy-level exceptions to them. A culture label matches the part of the gb_groups culture
after ' | ', exactly or as a suffix after one of ' _-:/'. expect holds ';'-separated conditions
tested before the edit is applied:
  organism=X  every matched copy is currently X     label=X   the copy's culture label ends with X
  n=N         N copies matched                      groups=N  N gb_groups cultures matched
  vt=VTX...   every matched copy with an assigned VT has it (and at least one does)
  class=X     every matched copy is currently in copy class X
  s58=N, s58=A-B  every matched copy has a 5.8S of N bp (or A to B bp)
A failed condition stops the run.

Culture = 'curated organism | culture label', the gb_groups convention, so Stage 1 tools read the
curated grouping unchanged (groups.tsv). Copy class: 'divergent' for divergent_5.8S copies, 'outlier'
for other ssu_outlier copies, else 'main', unless a copy_class edit says otherwise; only main-type
copies define a culture's labels. stefani_qc flags copies whose 5.8S is > 3% from the culture
medoid, which catches two things: the short-5.8S class (5-6 bp shorter than the main type) and
same-length copies with several 5.8S substitutions. Each copy gets s58_delta (its 5.8S length minus
the median of its culture's main-type copies) and divergent copies get div_kind (short_5.8S,
same_length or longer). Where the medoid falls in the short-5.8S class the flags come out inverted;
the copy-class check lists every culture whose divergent class is longer than its main type.

Writes <out>/ (refuses to overwrite a frozen version without --force):
  copies.tsv    every input copy: status, curated organism and culture, class, flags, VT, window status
  cultures.tsv  one row per retained culture
  groups.tsv    gb_groups.py format, retained copies, curated labels
  fasta/        full copies and every region for retained copies; primer windows only when complete
  MANIFEST.json inputs and outputs with sha256, edits with the copies they touched, counts, fingerprint

Usage (from amfsplit_stage3):
  python3 freeze_ref.py --edits ref_edits_v1.tsv --version v1 --out ../ref/culture_ref/v1
"""
import argparse
import collections
import datetime
import json
import os
import re
import shutil
import statistics
import sys

import edlib

from arms import (ORDER, genus_of, is_named, read_fasta, read_tsv, sha256, write_fasta,
                  write_tsv)

ITSX_REGIONS = ["ssu_flank", "its1", "s58", "its2", "full_its", "lsu_flank"]
ACTIONS = ("drop", "name_culture", "relabel", "copy_class", "check_name")
CLASSES = ("main", "divergent", "outlier")
ACC_ACTIONS = ("drop", "relabel", "copy_class")
NAME_REGIONS = [("its1", "ITS1"), ("its2", "ITS2"), ("LSU_LROR-FLR2", "LROR-FLR2"), ("ssu_flank", "SSU flank")]


def vt_key(c):
    """VT as tallied in tables: the VT when assigned, else the call with any tied VTs."""
    if c["vt_call"] == "assigned":
        return c["vt"]
    return f"({c['vt_call']} {c['vt_tied']})" if c.get("vt_tied") else f"({c['vt_call']})"


def label_of(culture):
    return culture.split(" | ", 1)[1] if " | " in culture else culture


def label_matches(label, target):
    return label == target or (label.endswith(target) and label[-len(target) - 1] in " _-:/")


def dist(a, b):
    """% edit distance, shorter sequence aligned inside the longer (end gaps free), as in Stage 1."""
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode="HW", task="distance")["editDistance"] / len(q)


def read_edits(path):
    rows, header = [], None
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            f = line.split("\t")
            if header is None:
                header = [x.strip() for x in f]
                if header[:2] != ["action", "target"]:
                    raise SystemExit(f"{path}: first non-comment line must be the header "
                                     "(action, target, value, expect, reason)")
                continue
            r = dict(zip(header, (x.strip() for x in f)))
            r["line"] = n
            if r.get("action") not in ACTIONS:
                raise SystemExit(f"{path} line {n}: unknown action {r.get('action')!r}; use one of {ACTIONS}")
            if r["action"] == "copy_class" and r.get("value") not in CLASSES:
                raise SystemExit(f"{path} line {n}: copy_class value must be one of {CLASSES}")
            rows.append(r)
    return rows


def check_expect(edit, accs, copies):
    """Return a list of failed conditions (empty when all hold)."""
    fails = []
    for cond in filter(None, (c.strip() for c in edit.get("expect", "").split(";"))):
        key, _, val = cond.partition("=")
        key, val = key.strip(), val.strip()
        if key == "organism":
            bad = sorted({copies[a]["organism"] for a in accs} - {val})
            if bad:
                fails.append(f"{cond} (found {bad})")
        elif key == "label":
            bad = sorted({copies[a]["culture_label"] for a in accs if not label_matches(copies[a]["culture_label"], val)})
            if bad:
                fails.append(f"{cond} (found {bad})")
        elif key == "n":
            if len(accs) != int(val):
                fails.append(f"{cond} (found {len(accs)})")
        elif key == "groups":
            g = {copies[a]["culture_gb"] for a in accs}
            if len(g) != int(val):
                fails.append(f"{cond} (found {len(g)}: {sorted(g)})")
        elif key == "vt":
            vts = collections.Counter(copies[a]["vt"] for a in accs if copies[a]["vt_call"] == "assigned")
            if not vts or set(vts) != {val}:
                fails.append(f"{cond} (found {dict(vts) or 'no assigned copies'})")
        elif key == "class":
            bad = sorted({copies[a]["copy_class"] for a in accs} - {val})
            if bad:
                fails.append(f"{cond} (found {bad})")
        elif key == "s58":
            lo, _, hi = val.partition("-")
            lo, hi = int(lo), int(hi or lo)
            found = collections.Counter(copies[a]["s58_len"] for a in accs)
            if any(n == "" or not lo <= n <= hi for n in found):
                fails.append(f"{cond} (found {dict(found)})")
        else:
            fails.append(f"{cond} (unknown condition)")
    return fails


def match_targets(edit, copies):
    t = edit["target"]
    if edit["action"] in ACC_ACTIONS:
        accs = [x.strip() for x in t.split(",") if x.strip()]
        absent = [x for x in accs if x not in copies]
        if absent or not accs:
            raise SystemExit(f"edits line {edit['line']}: accession(s) {absent or t!r} not in the input")
        return accs
    accs = [a for a, c in copies.items() if c["status"] == "retained" and label_matches(c["culture_label"], t)]
    if not accs:
        raise SystemExit(f"edits line {edit['line']}: no retained copy has a culture label matching {t!r}")
    return accs


def name_check(culture, cand, copies, seqs, held, gap):
    """Compare the main-type copies of a held culture with the other cultures of each candidate species."""
    mine = [a for a, c in copies.items() if c["culture"] == culture and c["status"] == "retained"]
    main = [a for a in mine if copies[a]["copy_class"] == "main"]
    vts = collections.Counter(vt_key(copies[a]) for a in mine)
    print(f"\n  {culture}: {len(mine)} copies, {len(main)} main-type; VT {dict(vts.most_common())}")
    if not main:
        print("    no main-type copies to compare; stays on hold")
        return None
    print(f"    median % from each main-type copy to the nearest main-type copy of the candidate's other cultures")
    print(f"    {'candidate':32s}{'cultures':>9s}" + "".join(f"{lab:>12s}" for _, lab in NAME_REGIONS) + "   modal VT of its cultures")
    med = {}
    for sp in cand:
        refs = collections.defaultdict(list)
        for a, c in copies.items():
            if c["status"] == "retained" and c["organism"] == sp and c["culture"] not in held and c["copy_class"] == "main":
                refs[c["culture"]].append(a)
        pool = [a for v in refs.values() for a in v]
        cells, med[sp] = [], {}
        for key, _ in NAME_REGIONS:
            ds = []
            for a in main:
                if a in seqs[key] and any(b in seqs[key] for b in pool):
                    ds.append(min(dist(seqs[key][a], seqs[key][b]) for b in pool if b in seqs[key]))
            med[sp][key] = statistics.median(ds) if ds else None
            cells.append(f"{med[sp][key]:.1f}" if ds else "-")
        modal = collections.Counter()
        for accs in refs.values():
            v = collections.Counter(vt_key(copies[a]) for a in accs)
            if v:
                modal[v.most_common(1)[0][0]] += 1
        print(f"    {sp[:31]:32s}{len(refs):>9d}" + "".join(f"{x:>12s}" for x in cells)
              + "   " + (", ".join(f"{k} ({n})" for k, n in modal.most_common(3)) or "-"))
    scored = [sp for sp in cand if all(med[sp].get(k) is not None for k in ("its1", "its2"))]
    if len(scored) < 2:
        print("    fewer than two candidates with ITS references; stays on hold")
        return None
    best = min(scored, key=lambda sp: med[sp]["its1"] + med[sp]["its2"])
    clear = all(med[o][k] - med[best][k] >= gap for o in scored if o != best for k in ("its1", "its2"))
    if clear:
        print(f"    -> closer to {best} by >= {gap:g} points in ITS1 and ITS2")
        return best
    print(f"    -> no candidate is closer by >= {gap:g} points in both ITS1 and ITS2; stays on hold")
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", default="../bench/stefani/stefani_groups.tsv")
    ap.add_argument("--flags", default="../bench/stefani/stefani_qc.flags.tsv")
    ap.add_argument("--vt", default="../bench/stefani/vt_check.copies.tsv", help="vt_check.py per-copy table")
    ap.add_argument("--prefix", default="../bench/stefani/out/windows/stefani", help="amfsplit window prefix")
    ap.add_argument("--fasta", default="../bench/stefani/stefani_rdna.fasta", help="full copies")
    ap.add_argument("--edits", required=True)
    ap.add_argument("--version", required=True, help="e.g. v1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--name-gap", type=float, default=1.0,
                    help="points by which a candidate must be closer in ITS1 and ITS2 to suggest it (default 1)")
    ap.add_argument("--force", action="store_true", help="overwrite an existing version directory")
    a = ap.parse_args()

    if os.path.exists(a.out) and os.listdir(a.out) and not a.force:
        raise SystemExit(f"{a.out} exists; a frozen version is not overwritten (use a new --version, or --force)")

    groups = {r["accession"]: r for r in read_tsv(a.groups)}
    vt = {r["accession"]: r for r in read_tsv(a.vt)}
    flags = collections.defaultdict(set)
    if a.flags and os.path.exists(a.flags):
        for r in read_tsv(a.flags):
            flags[r["accession"]].add(r["check"])
    records = {r["read_id"]: r for r in read_tsv(f"{a.prefix}.records.tsv")}
    full = read_fasta(a.fasta)
    s58_path = f"{a.prefix}.s58.fasta"
    s58_len = {k: len(v) for k, v in read_fasta(s58_path).items()} if os.path.exists(s58_path) else {}
    windows = sorted(k[:-len("_status")] for k in next(iter(records.values())) if k.endswith("_status"))

    missing = {name: sorted(set(full) - set(d))[:3] for name, d in
               (("groups", groups), ("vt table", vt), ("records", records)) if set(full) - set(d)}
    if missing:
        raise SystemExit(f"accessions in {a.fasta} missing from: {missing}")

    copies = {}
    for acc in full:
        g, v = groups[acc], vt[acc]
        fl = ",".join(sorted(flags.get(acc, set()) or set(filter(None, v.get("flags", "").split(",")))))
        copies[acc] = {
            "accession": acc, "status": "retained", "reason": "", "organism": g["organism"],
            "organism_gb": g["organism"], "culture": g["culture"], "culture_gb": g["culture"],
            "culture_label": label_of(g["culture"]), "library": g.get("library", ""),
            "library_gb": g.get("library", ""),
            "copy_class": ("divergent" if "divergent_5.8S" in fl else "outlier" if "ssu_outlier" in fl
                           else "main"), "flags": fl, "s58_len": s58_len.get(acc, ""),
            "name_status": "ok", "edits": "",
            "vt_call": v.get("call", ""), "vt": v.get("vt", ""), "vt_pident": v.get("pident", ""),
            "vt_margin": v.get("margin", ""), "vt_tied": v.get("tied_vts", ""), "vt_best": v.get("best_vt", ""),
            "vt_best_pident": v.get("best_pident", ""),
            **{f"{w}_status": records[acc].get(f"{w}_status", "") for w in windows},
        }

    # ---- apply edits by action (drop, name_culture, relabel, copy_class, check_name), each checked first
    edits = read_edits(a.edits)
    print(f"culture reference {a.version}: {len(copies)} copies in, {len(edits)} edits from {a.edits}")
    held, candidates, log = set(), {}, []
    for e in sorted(edits, key=lambda e: ACTIONS.index(e["action"])):
        accs = match_targets(e, copies)
        fails = check_expect(e, accs, copies)
        if fails:
            raise SystemExit(f"edits line {e['line']} ({e['action']} {e['target']}): expectation failed: "
                             + "; ".join(fails) + "\nNothing was written. Check the data or the edit.")
        tag = f"{e['action']}:{e['target']}"
        touched = sorted({copies[x]["culture"] for x in accs})
        for x in accs:
            c = copies[x]
            c["edits"] = ",".join(filter(None, [c["edits"], tag]))
            if e["action"] == "drop":
                c["status"], c["reason"] = "excluded", e.get("reason", "")
            elif e["action"] == "copy_class":
                c["copy_class"] = e["value"]
            elif e["action"] in ("relabel", "name_culture"):
                c["organism"] = e["value"]
                c["culture"] = f"{e['value']} | {c['culture_label']}"
                if c["library"]:
                    c["library"] = f"{e['value']} | {label_of(c['library'])}"
        if e["action"] == "check_name":
            now = sorted({copies[x]["culture"] for x in accs})
            held.update(now)
            for cu in now:
                candidates[cu] = [s.strip() for s in e["value"].split(",") if s.strip()]
            for x in accs:
                copies[x]["name_status"] = "hold"
        after = sorted({copies[x]["culture"] for x in accs})
        log.append({**{k: e.get(k, "") for k in ("action", "target", "value", "expect", "reason")},
                    "copies": accs, "cultures_before": touched, "cultures_after": after})
        what = {"drop": "dropped", "relabel": f"-> {e['value']}", "name_culture": f"-> {e['value']}",
                "copy_class": f"class -> {e['value']}", "check_name": "held for name check"}[e["action"]]
        shown = e["target"] if len(e["target"]) <= 17 else e["target"].split(",")[0] + ",.."
        print(f"  {e['action']:13s}{shown:18s}{len(accs):>3d} copies, {len(touched)} culture(s) {what}"
              + ("" if e["action"] in ("drop", "copy_class", "check_name") else f"; now {len(after)} culture(s)"))

    kept = {x: c for x, c in copies.items() if c["status"] == "retained"}
    for c in kept.values():
        c["genus"] = genus_of(c["organism"])
        c["order"] = ORDER.get(c["genus"], "unassigned")
    unknown = sorted({c["genus"] for c in kept.values() if c["order"] == "unassigned"})
    if unknown:
        print(f"  note: genera without an order in arms.ORDER: {unknown}")

    # copy classes against the 5.8S length signature: s58_delta per copy, div_kind per divergent copy
    by_c = collections.defaultdict(lambda: collections.defaultdict(list))
    for c in kept.values():
        if c["s58_len"] != "":
            by_c[c["culture"]][c["copy_class"]].append(c["s58_len"])
    for c in copies.values():
        mains = by_c.get(c["culture"], {}).get("main") if c["status"] == "retained" else None
        c["s58_delta"] = (round(c["s58_len"] - statistics.median(mains))
                          if mains and c["s58_len"] != "" else "")
        d = c["s58_delta"]
        c["div_kind"] = ("" if c["copy_class"] != "divergent" or d == "" else
                         "short_5.8S" if d <= -3 else "longer" if d >= 3 else "same_length")
    class_check = {"inverted": [], "same_length": [], "div_kinds": {}}
    both = {cu: d for cu, d in by_c.items() if d["main"] and d["divergent"]}
    if both:
        allm = [n for d in both.values() for n in d["main"]]
        alld = [n for d in both.values() for n in d["divergent"]]
        kinds = collections.Counter(c["div_kind"] or "no 5.8S" for c in kept.values() if c["copy_class"] == "divergent")
        class_check["div_kinds"] = dict(kinds)
        print(f"\ncopy classes vs 5.8S length ({len(both)} cultures hold both): main-type median "
              f"{statistics.median(allm):g} bp, divergent-class median {statistics.median(alld):g} bp")
        print("  divergent copies by 5.8S length against their culture's main type: "
              + ", ".join(f"{k} {v}" for k, v in kinds.most_common()))
        for cu in sorted(both):
            ks = [c["div_kind"] for c in kept.values() if c["culture"] == cu and c["copy_class"] == "divergent"]
            row = {"culture": cu, "main": sorted(by_c[cu]["main"]), "divergent": sorted(by_c[cu]["divergent"])}
            if statistics.median(by_c[cu]["divergent"]) - statistics.median(by_c[cu]["main"]) >= 3:
                class_check["inverted"].append(row)
            elif ks and all(k == "same_length" for k in ks):
                class_check["same_length"].append(row)
        for r in class_check["inverted"]:
            print(f"  {r['culture']}: divergent class LONGER than the main type (main {r['main']}, divergent "
                  f"{r['divergent']} bp)  <- classes likely inverted; see copy_class")
        if not class_check["inverted"]:
            print("  no culture has its divergent class longer than its main type")
        if class_check["same_length"]:
            print(f"  same-length divergent copies only (5.8S substitutions, not the short-5.8S class): "
                  + ", ".join(r["culture"] for r in class_check["same_length"]))

    # labels filed under several organisms that no edit merged: review list
    by_label = collections.defaultdict(set)
    for c in kept.values():
        by_label[c["culture_label"]].add(c["organism"])
    shared = {lab: sorted(orgs) for lab, orgs in by_label.items() if len(orgs) > 1}
    print(f"\nculture labels still under more than one organism: {len(shared)}")
    for lab, orgs in sorted(shared.items()):
        n = collections.Counter(c["organism"] for c in kept.values() if c["culture_label"] == lab)
        print(f"  {lab}: " + ", ".join(f"{o} x{n[o]}" for o in orgs)
              + ("  (mixed culture, by edit)" if any("relabel" in c["edits"] for c in kept.values()
                                                   if c["culture_label"] == lab) else "  <- review"))

    # ---- sequences for retained copies
    seqs = {}
    for r in ITSX_REGIONS:
        p = f"{a.prefix}.{r}.fasta"
        seqs[r] = {k: s for k, s in read_fasta(p).items() if k in kept} if os.path.exists(p) else {}
    for w in windows:
        p = f"{a.prefix}.{w}.fasta"
        seqs[w] = ({k: s for k, s in read_fasta(p).items() if k in kept and kept[k][f"{w}_status"] == "ok"}
                   if os.path.exists(p) else {})
    if not seqs["full_its"] and all(seqs[r] for r in ("its1", "s58", "its2")):
        seqs["full_its"] = {k: seqs["its1"][k] + seqs["s58"][k] + seqs["its2"][k]
                            for k in seqs["its1"] if k in seqs["s58"] and k in seqs["its2"]}
        print("  note: no full_its.fasta; built full ITS as its1 + s58 + its2")
    for need in ("ssu_flank", "full_its", "LSU_LROR-FLR2"):
        lack = len(kept) - len(seqs.get(need, {}))
        if lack:
            print(f"  note: {lack} retained copies have no {need} sequence")

    # ---- name checks
    suggestions = {}
    if held:
        print(f"\nname check ({len(held)} culture(s) on hold)")
        for cu in sorted(held):
            best = name_check(cu, candidates[cu], copies, seqs, held, a.name_gap)
            if best:
                suggestions[cu] = best
        if suggestions:
            print("\n  for the next version:")
            for cu, sp in sorted(suggestions.items()):
                n = sum(c["culture"] == cu for c in kept.values())
                named_by_edit = any("name_culture:" in c["edits"] for c in kept.values() if c["culture"] == cu)
                if sp == cu.split(" | ")[0]:
                    print(f"    {cu}: ITS supports the current name; delete its check_name line")
                elif named_by_edit:
                    print(f"    {cu}: set its name_culture value to {sp} and delete its check_name line")
                else:
                    print(f"    {cu}: replace its check_name line with\n"
                          f"      name_culture\t{label_of(cu)}\t{sp}\tn={n}\tname check {a.version}: ITS nearest to {sp}")

    # ---- cultures table
    by_cul = collections.defaultdict(list)
    for x, c in kept.items():
        by_cul[c["culture"]].append(c)
    cultures = []
    for cu, cs in sorted(by_cul.items()):
        main_vt = collections.Counter(c["vt"] for c in cs if c["copy_class"] == "main" and c["vt_call"] == "assigned")
        all_vt = collections.Counter(vt_key(c) for c in cs)
        libs = sorted({label_of(c["library"]) for c in cs if c["library"]})
        cultures.append({
            "culture": cu, "organism": cs[0]["organism"], "culture_label": cs[0]["culture_label"],
            "genus": cs[0]["genus"], "order": cs[0]["order"], "named": is_named(cs[0]["organism"]),
            "name_status": "hold" if cu in held else "ok", "n_copies": len(cs),
            "n_main": sum(c["copy_class"] == "main" for c in cs),
            "n_divergent": sum(c["copy_class"] == "divergent" for c in cs),
            "n_outlier": sum(c["copy_class"] == "outlier" for c in cs),
            "n_libraries": len(libs), "libraries": ",".join(libs),
            "vt_modal": main_vt.most_common(1)[0][0] if main_vt else "",
            "vt_counts": ",".join(f"{k}:{n}" for k, n in all_vt.most_common()),
            "organisms_gb": ",".join(sorted({c["organism_gb"] for c in cs})),
            "cultures_gb": ",".join(sorted({c["culture_gb"] for c in cs})),
            "edits": ",".join(sorted({e for c in cs for e in c["edits"].split(",") if e})),
        })

    # ---- write
    if os.path.exists(a.out) and a.force:
        shutil.rmtree(a.out)
    os.makedirs(os.path.join(a.out, "fasta"))
    ccols = (["accession", "status", "reason", "organism", "organism_gb", "culture", "culture_gb",
              "culture_label", "library", "library_gb", "genus", "order", "copy_class", "flags", "s58_len", "s58_delta", "div_kind", "name_status", "edits",
              "vt_call", "vt", "vt_pident", "vt_margin", "vt_tied", "vt_best", "vt_best_pident"]
             + [f"{w}_status" for w in windows])
    for c in copies.values():
        c.setdefault("genus", genus_of(c["organism"]))
        c.setdefault("order", ORDER.get(c["genus"], "unassigned"))
    write_tsv(os.path.join(a.out, "copies.tsv"), list(copies.values()), ccols)
    write_tsv(os.path.join(a.out, "cultures.tsv"), cultures, list(cultures[0]))
    gcols = list(next(iter(groups.values())))
    grows = [{**groups[x], "organism": c["organism"], "genus": c["genus"], "culture": c["culture"],
              "library": c["library"]} for x, c in kept.items()]
    write_tsv(os.path.join(a.out, "groups.tsv"), grows, gcols)
    write_fasta(os.path.join(a.out, "fasta", "full.fasta"), {x: full[x] for x in kept})
    for r, s in seqs.items():
        if s:
            write_fasta(os.path.join(a.out, "fasta", f"{r}.fasta"), s)

    named = {c["organism"] for c in cultures if c["named"] and c["name_status"] == "ok"}
    by_order = collections.defaultdict(lambda: [0, 0])
    for c in cultures:
        by_order[c["order"]][0] += 1
        by_order[c["order"]][1] += c["n_copies"]
    outputs = {os.path.relpath(os.path.join(d, f), a.out): sha256(os.path.join(d, f))
               for d, _, fs in os.walk(a.out) for f in sorted(fs)}
    counts = {"copies_in": len(copies), "retained": len(kept), "excluded": len(copies) - len(kept),
              "relabelled": sum(c["organism"] != c["organism_gb"] for c in kept.values()),
              "divergent_class": sum(c["copy_class"] == "divergent" for c in kept.values()),
              "outlier_class": sum(c["copy_class"] == "outlier" for c in kept.values()),
              "cultures": len(cultures), "named_species": len(named),
              "cultures_on_name_hold": sorted(held),
              "by_order": {o: {"cultures": v[0], "copies": v[1]} for o, v in sorted(by_order.items())}}
    here = os.path.abspath(__file__)
    manifest = {
        "name": "culture_ref", "version": a.version,
        "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "source": "Stefani et al. (2025) rDNA copies, GenBank PX207096-PX207274 and PX214440-PX215164",
        "fingerprint": outputs["copies.tsv"][:16],
        "script": {"file": os.path.basename(here), "sha256": sha256(here)},
        "inputs": {k: {"path": p, "sha256": sha256(p), "bytes": os.path.getsize(p)} for k, p in
                   [("groups", a.groups), ("flags", a.flags), ("vt", a.vt), ("records", f"{a.prefix}.records.tsv"),
                    ("fasta", a.fasta), ("edits", a.edits)] if p and os.path.exists(p)},
        "edits": log, "name_check_suggestions": suggestions, "shared_labels": shared,
        "copy_class_check": class_check,
        "counts": counts, "outputs": outputs,
    }
    with open(os.path.join(a.out, "MANIFEST.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
        fh.write("\n")

    print(f"\n{a.version}: {counts['retained']} copies retained, {counts['excluded']} excluded, "
          f"{counts['relabelled']} relabelled; {counts['cultures']} cultures "
          f"({len(held)} on name hold), {counts['named_species']} named species; "
          f"{counts['divergent_class']} divergent-class and {counts['outlier_class']} SSU-outlier copies")
    for o, v in sorted(by_order.items()):
        print(f"  {o:<17}{v[0]:>4} cultures {v[1]:>5} copies")
    print(f"fingerprint {manifest['fingerprint']}  ->  {a.out}/")


if __name__ == "__main__":
    sys.exit(main())
