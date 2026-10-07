#!/usr/bin/env python3
"""LongGloDB against the v4 culture reference: what is left for an independent test?

LongGloDB headers are  Genus_epithet_CODE..._LAB_ASVnnnn_COLLECTION, one ASV per sequence, several ASVs per
culture. The species is the first two fields and the culture code everything between it and the lab.

  0. Sources. The share of each lab's ASVs that is identical in ITS to a Stefani copy. The labs given in
     --stefani-labs (STFP, STFI: collection CCAMF) are the Stefani cultures reprocessed and are not
     independent.
  1. Culture codes shared with Stefani labels, outside that part: tokens with letters and digits (ON205A,
     BEG9), a letter token joined to the number after it (BEG_35 = BEG35), leading zeros ignored. Bare
     numbers and short tokens (11, 9, 4A) are replicate numbers, not codes, and are not matched.
  2. Roles for the rest: known species (named in v4), novel species with the genus in v4, novel lineage
     (genus not in v4), genus only. Two flags per culture: an ASV identical in ITS to a Stefani copy of the
     same species (possibly the same isolate held elsewhere), and a V18 record >= 99.5% identical in LSU
     (the LSU arm may see the culture's own record).
Names: Rhizoglomus is read as Rhizophagus, and the spelling variants in RENAME (extend with --rename) are
mapped to the v4 names; both are printed. Distances: edlib infix edit distance, % of the shorter sequence.

Writes {out}.cultures.tsv, {out}.asvs.tsv (with the role of every ASV) and the test set as
{out}.test_full.fasta, .test_ssu, .test_its and .test_lsu (known species, novel species, novel lineage and
genus-only cultures).

Usage (from amfsplit_stage3):
  python3 longglodb_overlap.py --longglodb ../ref/glodb --ref ../ref/culture_ref/v4 \
      --v18 ../ref/delavaux/repo/2024_AMFPipeline/2024_AMFPipeline_ASV/V18_LSUDB_052025_AMFONLY.fasta \
      --out ../bench/longglodb/overlap
"""
import argparse
import collections
import csv
import glob
import os
import re
import statistics
import sys
import time

import edlib

HDR = re.compile(r"^(?P<pre>.+)_(?P<lab>[^_]+)_(?P<asv>ASV\d+)_(?P<coll>[^_]+)$")
GENUS_SYN = {"rhizoglomus": "rhizophagus"}
# spelling and gender variants seen in LongGloDB (printed in the log; extend with --rename)
RENAME = {"Funneliformis coronatum": "Funneliformis coronatus", "Sclerocystis sinuos": "Sclerocystis sinuosa"}


def parse_header(h, rename):
    m = HDR.match(h)
    if not m:
        return None
    t = m["pre"].split("_")
    if len(t) < 2:
        return None
    if t[1].rstrip(".").lower() in ("sp", "spp"):
        t = [t[0], "sp."] + t[2:]
    org = " ".join(t[:2])
    return {"org": rename.get(org, org), "org_raw": org, "code": "_".join(t[2:]), "code_tokens": t[2:],
            "lab": m["lab"], "asv": m["asv"], "coll": m["coll"]}


def is_code(t):
    return len(t) >= 3 and bool(re.search(r"[A-Za-z]", t)) and bool(re.search(r"\d", t)) and not re.fullmatch(r"s\d+", t)


def code_set(tokens):
    """Culture-code keys from a label: tokens with letters and digits (ON205A, BEG9, FL879N), plus a letter
    token joined to the number after it (BEG_35 -> BEG35). Bare numbers (11, 9) and short tokens (4A) are
    left out: they are sample or replicate numbers, not culture codes."""
    out = {code_keys(t)[1] for t in tokens if is_code(t)}
    for a, b in zip(tokens, tokens[1:]):
        if re.fullmatch(r"[A-Za-z]+", a) and re.fullmatch(r"\d+[A-Za-z]?", b):
            out.add(code_keys(a + b)[1])
    return out


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def dist(a, b, k=-1):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    d = edlib.align(q, t, mode="HW", task="distance", k=k)["editDistance"]
    return None if d < 0 else 100 * d / len(q)


def nearest(q, cands, seqs):
    best, bd = None, None
    for y in cands:
        s = seqs.get(y)
        if not s:
            continue
        k = -1 if bd is None else int(bd * min(len(q), len(s)) / 100) + 1
        d = dist(q, s, k)
        if d is not None and (bd is None or d < bd):
            best, bd = y, d
    return best, bd


def norm_species(name):
    t = name.replace("_", " ").split()
    if not t:
        return "", ""
    g = GENUS_SYN.get(t[0].strip("[]").lower(), t[0].strip("[]").lower())
    ep = t[1].lower().rstrip(".") if len(t) > 1 else ""
    named = bool(ep) and ep not in ("sp", "spp", "cf", "aff")
    return (f"{g} {ep}" if named else ""), g


def code_keys(s):
    exact = re.sub(r"[^A-Z0-9]", "", s.upper())
    loose = "".join(p.lstrip("0") or "0" if p.isdigit() else p for p in re.findall(r"[A-Z]+|\d+", exact))
    return exact, loose


def find_file(d, pattern):
    hits = sorted(glob.glob(os.path.join(d, pattern)))
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--longglodb", required=True, help="folder with the four LongGloDB FASTA files")
    ap.add_argument("--ref", required=True, help="culture reference folder (v4)")
    ap.add_argument("--v18", help="V18 AMF-only FASTA (optional: the self-match check)")
    ap.add_argument("--stefani-labs", default="STFP,STFI",
                    help="LongGloDB lab fields that are the Stefani cultures reprocessed")
    ap.add_argument("--rename", help="TSV of LongGloDB name -> name used in v4 (adds to the built-in list)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    rename = dict(RENAME)
    if a.rename:
        for line in open(a.rename):
            p = line.rstrip("\n").split("\t")
            if len(p) == 2 and p[0] and not p[0].startswith("#"):
                rename[p[0]] = p[1]
    stf_labs = set(split_list(a.stefani_labs))

    # ---- LongGloDB
    f_full = find_file(a.longglodb, "LongGloDB_July*.fasta") or find_file(a.longglodb, "LongGloDB_[!SIL]*.fasta")
    f_its = find_file(a.longglodb, "LongGloDB_ITS_*.fasta")
    f_lsu = find_file(a.longglodb, "LongGloDB_LSU_*.fasta")
    f_ssu = find_file(a.longglodb, "LongGloDB_SSU_*.fasta")
    if not f_full or not f_its:
        raise SystemExit(f"need the full and ITS LongGloDB files in {a.longglodb}")
    full, its = read_fasta(f_full), read_fasta(f_its)
    lsu = read_fasta(f_lsu) if f_lsu else {}
    ssu = read_fasta(f_ssu) if f_ssu else {}
    asvs, bad = {}, []
    for h in full:
        m = parse_header(h, rename)
        (asvs.__setitem__(h, m) if m else bad.append(h))
    lens = sorted(len(full[h]) for h in asvs) or [0]
    print(f"LongGloDB: {len(full)} sequences (length median {statistics.median(lens):.0f}, {lens[0]}-{lens[-1]} bp); "
          f"in ITS file {sum(h in its for h in full)}, LSU {sum(h in lsu for h in full)}, SSU {sum(h in ssu for h in full)}")
    if bad:
        print(f"  {len(bad)} headers not parsed (left out), e.g. {bad[:3]}")
    used = sorted({(m['org_raw'], m['org']) for m in asvs.values() if m['org_raw'] != m['org']})
    if used:
        print("  names changed to match v4: " + "; ".join(f"{x} -> {y}" for x, y in used))

    # ---- v4 reference
    rows = list(csv.DictReader(open(os.path.join(a.ref, "copies.tsv")), delimiter="\t"))
    ref = {r["accession"]: r for r in rows if r.get("status", "retained") == "retained"}
    mains = [x for x, r in ref.items() if r.get("copy_class", "main") == "main"]
    v4_its = read_fasta(os.path.join(a.ref, "fasta", "full_its.fasta"))
    v4_cultures = sorted({r["culture"] for r in ref.values()})
    v4_species, v4_genera = set(), set()
    stef_codes = collections.defaultdict(set)
    for cu in v4_cultures:
        org, label = (cu.split(" | ", 1) + [""])[:2]
        sp, g = norm_species(org)
        if sp:
            v4_species.add(sp)
        v4_genera.add(g)
        for c in code_set(re.split(r"[_\-\s]+", label)):
            stef_codes[c].add(cu)
    print(f"v4: {len(ref)} copies, {len(v4_cultures)} cultures, {len(v4_species)} named species, {len(v4_genera)} genera")

    # ---- per ASV: nearest Stefani main-type copy (full ITS), nearest V18 record (LSU)
    v18 = {h: s.strip("N") for h, s in read_fasta(a.v18).items()} if a.v18 else {}
    v18_sp = {h: " ".join(h.split("_")[2:4]) if len(h.split("_")) >= 4 else h for h in v18}
    genus_of_copy = {x: norm_species(ref[x]["culture"].split(" | ")[0])[1] for x in mains}
    asv_rows = []
    for i, (h, m) in enumerate(asvs.items()):
        sp, g = norm_species(m["org"])
        row = {"asv": h, "organism": m["org"], "code": m["code"], "lab": m["lab"], "collection": m["coll"],
               "length": len(full[h]), "nearest_v4": "", "nearest_its_pct": "", "v18_best": "", "v18_best_pct": ""}
        if h in its:
            cands = sorted((x for x in mains if x in v4_its), key=lambda x: genus_of_copy[x] != g)
            y, d = nearest(its[h], cands, v4_its)
            if y:
                row.update(nearest_v4=ref[y]["culture"], nearest_its_pct=round(d, 2))
        if v18 and h in lsu:
            y, d = nearest(lsu[h], sorted(v18, key=lambda r: norm_species(v18_sp[r])[1] != g), v18)
            if y:
                row.update(v18_best=v18_sp[y], v18_best_pct=round(100 - d, 2))
        asv_rows.append(row)
        if (i + 1) % 250 == 0:
            print(f"  {i + 1} ASVs compared ({time.time() - t0:.0f} s)")

    # ---- 0. sources
    print("\n0. sources: share of ASVs identical in ITS to a Stefani copy")
    src = collections.defaultdict(list)
    for r in asv_rows:
        src[(r["lab"], r["collection"])].append(r)
    for (lab, coll), rs in sorted(src.items(), key=lambda kv: -len(kv[1])):
        ds = [r["nearest_its_pct"] for r in rs if r["nearest_its_pct"] != ""]
        tag = "  <- Stefani cultures reprocessed" if lab in stf_labs else ""
        print(f"  {lab:6s} {coll:14s}{len(rs):>5d} ASVs   0.0%: {sum(d == 0 for d in ds):>4d}   "
              f"median {statistics.median(ds) if ds else float('nan'):.2f}%{tag}")

    # ---- per culture (organism + full code), roles
    cult = collections.defaultdict(list)
    for r in asv_rows:
        cult[(r["organism"], r["code"])].append(r)
    out = []
    for (org, code), rs in sorted(cult.items()):
        sp, g = norm_species(org)
        lab = rs[0]["lab"]
        hits = sorted({cu for c in code_set(asvs[rs[0]["asv"]]["code_tokens"]) for cu in stef_codes.get(c, ())})
        agree = any(norm_species(c.split(" | ")[0])[0] == sp for c in hits) if hits and sp else False
        ds = [(r["nearest_its_pct"], r["nearest_v4"]) for r in rs if r["nearest_its_pct"] != ""]
        best_d, best_cu = min(ds) if ds else ("", "")
        identical = any(r["nearest_its_pct"] == 0 and sp and norm_species(r["nearest_v4"].split(" | ")[0])[0] == sp
                        for r in rs)
        v = [r["v18_best_pct"] for r in rs if r["v18_best_pct"] != ""]
        if lab in stf_labs:
            role = "Stefani reprocessed"
        elif hits:
            role = "code shared with Stefani" + ("" if agree else ", names differ")
        elif not sp:
            role = "genus only"
        elif sp in v4_species:
            role = "known species"
        elif g in v4_genera:
            role = "novel species, genus in v4"
        else:
            role = "novel lineage, genus not in v4"
        for r in rs:
            r["role"] = role
        out.append({"organism": org, "code": code, "lab": lab, "collection": rs[0]["collection"], "asvs": len(rs),
                    "role": role, "stefani_code_match": "; ".join(hits),
                    "identical_to_stefani_same_species": "yes" if identical else "",
                    "nearest_v4_culture": best_cu, "nearest_its_pct": best_d,
                    "v18_best_pct": max(v) if v else "", "v18_best_species":
                    next((r["v18_best"] for r in rs if v and r["v18_best_pct"] == max(v)), "")})

    with open(f"{a.out}.cultures.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    acols = ["asv", "organism", "code", "lab", "collection", "role", "length", "nearest_v4", "nearest_its_pct",
             "v18_best", "v18_best_pct"]
    with open(f"{a.out}.asvs.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=acols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(asv_rows)
    test_roles = ("known species", "novel species, genus in v4", "novel lineage, genus not in v4", "genus only")
    test = [r["asv"] for r in asv_rows if r["role"] in test_roles]
    for name, seqs in (("full", full), ("ssu", ssu), ("its", its), ("lsu", lsu)):
        if seqs:
            with open(f"{a.out}.test_{name}.fasta", "w", newline="\n") as fh:
                for h in test:
                    if h in seqs:
                        fh.write(f">{h}\n{seqs[h]}\n")

    # ---- report
    print("\n1. culture codes shared with the Stefani set (outside the reprocessed part)")
    for r in [r for r in out if r["role"].startswith("code shared")]:
        print(f"  {r['organism']} {r['code']} [{r['lab']}/{r['collection']}] -> {r['stefani_code_match']} "
              f"(nearest ITS {r['nearest_its_pct']}% to {r['nearest_v4_culture']})")

    print("\n2. what is left for a test")
    print(f"  {'role':40s}{'cultures':>9s}{'ASVs':>6s}{'species':>8s}{'identical ASV':>15s}{'V18 >=99.5%':>13s}")
    for k in ("Stefani reprocessed", "code shared with Stefani", "code shared with Stefani, names differ") + test_roles:
        rs = [r for r in out if r["role"] == k]
        if not rs:
            continue
        hi = sum(r["v18_best_pct"] != "" and r["v18_best_pct"] >= 99.5 for r in rs)
        print(f"  {k:40s}{len(rs):>9d}{sum(r['asvs'] for r in rs):>6d}{len({r['organism'] for r in rs}):>8d}"
              f"{sum(bool(r['identical_to_stefani_same_species']) for r in rs):>15d}{hi:>13d}")
    for k in ("novel species, genus in v4", "novel lineage, genus not in v4", "genus only"):
        rs = [r for r in out if r["role"] == k]
        if rs:
            print(f"    {k}: " + "; ".join(f"{r['organism']} {r['code']}" for r in rs))
    print("  identical ASV: an ASV 0.0% in ITS from a Stefani copy of the same species (perhaps the same isolate held"
          " in another collection, or a species with uniform ITS)")
    print("  V18 >=99.5%: the LSU arm may see a record from the same culture; drop or report it in the test")

    unknown = sorted({r["organism"] for r in out if r["role"].startswith(("novel species", "novel lineage"))})
    hints = []
    for o in unknown:
        sp, _ = norm_species(o)
        ep = sp.split()[1] if sp else ""
        for vsp in sorted(v4_species):
            vg, vep = vsp.split()
            if len(ep) >= 6 and vep[:6] == ep[:6] and vep != ep:
                hints.append(f"{o} ~ {vg.capitalize()} {vep}")
    if hints:
        print("  possible synonyms left (decide before the test; add to --rename): " + "; ".join(hints))
    print(f"\ntables: {a.out}.cultures.tsv, {a.out}.asvs.tsv; test sequences: {a.out}.test_{{full,ssu,its,lsu}}.fasta "
          f"({len(test)} ASVs)  ({time.time() - t0:.0f} s)")


def split_list(s):
    return [x for x in (s or "").split(",") if x]


if __name__ == "__main__":
    sys.exit(main())
