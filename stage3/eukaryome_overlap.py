#!/usr/bin/env python3
"""EUKARYOME v2.1 Glomeromycota long reads against the v4 culture reference: names, roles and removals.

The test set of PROTOCOL_EUK.md. EUKARYOME reads carry no culture codes, only EUKARYOME's classification
(k__;p__;c__;o__;f__;g__;s__ in the header), so a read's role comes from its labels, mapped to the names of
culture reference v4 by a fixed table (names.tsv). Nothing here searches MaarjAM, UNITE or V18, places a
read or scores a call.

Subcommands:
  names   Lists every EUKARYOME label of the test reads with a suggested mapping to v4 and to V18 labels,
          and marks the cases that need a decision. Its output is reviewed and saved as names.tsv, which
          is fixed (hashed) before the predictions.
  roles   Applies names.tsv: one role per read in the vocabulary of predict_longglodb.py and lgdb_test.py
          (known species; novel species, genus in v4; novel lineage, genus not in v4; genus only), plus a
          stratum for reporting. Reads whose full ITS is identical to a v4 copy are removed (role
          'excluded'). Flags: IDs that are not EUKARYOME's own (INSDC, UNITE or assembly accessions), and,
          with --v18, the best V18 record in LSU (>= 99.5% is reported apart, as in the LongGloDB test);
          with --longglodb, reads identical in full ITS to a public LongGloDB v1 ASV (counts only).
          Also writes the V18 removal plan: one round per absent species whose label V18 holds and one per
          absent genus that V18 holds (all its records out).
  vtgeom  Where each read starts inside the NS31-AM1 fragment, from the Stage 1 windows and the v4 NS31-AM1
          copies; with --maarjam, the share of MaarjAM references a read starting there could cover at the
          VT arm's 95%. Describes the reads and the database only; no read is searched against MaarjAM.

Names: Rhizoglomus is read as Rhizophagus, as in the paper. An epithet matches when its stem matches
(Latin endings -us/-um/-a/-is/-e removed: irregulare = irregularis). 'cf.' names count as genus only.
EUKARYOME placeholder taxa (Glomeraceae.gen25, Archaeosporomycetes.ord01, ...) are lineages without a named
member and count as genera not in v4. Distances: edlib infix edit distance, % of the shorter sequence.

Usage (from the repository root):
  python3 stage3/eukaryome_overlap.py names --windows STAGE1/windows/euk_test \\
      --ref reference/culture_ref_v4 --out bench/eukaryome/names_candidates.tsv
  python3 stage3/eukaryome_overlap.py roles --windows STAGE1/windows/euk_test --ref reference/culture_ref_v4 \\
      --names benchmarks/eukaryome/names.tsv [--v18 V18_LSUDB_052025_AMFONLY.fasta] [--longglodb ref/glodb] \\
      --out bench/eukaryome/overlap
  python3 stage3/eukaryome_overlap.py vtgeom --windows STAGE1/windows/euk_test --ref reference/culture_ref_v4 \\
      [--maarjam ref/maarjam/maarjam.fasta] --out bench/eukaryome/vtgeom
"""
import argparse
import collections
import csv
import glob
import hashlib
import os
import re
import sys
import time

import edlib

HERE = os.path.dirname(os.path.abspath(__file__))
V18_META = os.path.join(HERE, "..", "benchmarks", "m4_external_holdout", "frozen_m2b", "reference",
                        "reference_metadata.tsv")
GENUS_SYN = {"Rhizoglomus": "Rhizophagus"}
ROLES = ("known species", "novel species, genus in v4", "novel lineage, genus not in v4", "genus only")
METHOD_WINDOWS = ("ssu_flank", "full_its", "s58", "LSU_LROR-FLR2")   # what lgdb_test.py reads (WINDOWS)
MARKERS = (("full_its", "ITS"), ("LSU_LROR-FLR2", "LSU"), ("ssu_flank", "SSU"))
PLACEHOLDER = re.compile(r"\.(gen|fam|ord)\d+$")
ENDINGS = ("ensis", "us", "um", "is", "ae", "a", "e", "i")
LAB, COLLECTION = "EUK", "EUKARYOME v2.1"


# ---------------------------------------------------------------- io
def read_fasta(path, full_header=False):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:] if full_header else line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def dist(a, b, k=-1):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    d = edlib.align(q, t, mode="HW", task="distance", k=k)["editDistance"]
    return None if d < 0 else 100 * d / len(q)


# ---------------------------------------------------------------- names
def stem(ep):
    """Epithet without its Latin ending and trailing vowels: irregulare = irregularis, bareae = bareai,
    drummondii = drummondi, viscosa = viscosum."""
    ep = ep.lower().rstrip(".")
    for e in ("ensis", "us", "um", "is"):
        if ep.endswith(e) and len(ep) - len(e) >= 4:
            ep = ep[: -len(e)]
            break
    s = ep.rstrip("aeiou")
    return s if len(s) >= 3 else ep


def v4_genus(g):
    g = g.strip("[]")
    return GENUS_SYN.get(g, g)


def parse_header(h):
    """'EUK1124192;k__Fungi;...;g__Dominikia;s__minuta region=...' -> (read id, {rank: value})."""
    parts = h.split()[0].split(";")
    tax = {}
    for x in parts[1:]:
        if "__" in x:
            k, v = x.split("__", 1)
            tax[k] = v
    return parts[0], tax


def label_key(tax):
    return (tax.get("o", "unclassified"), tax.get("f", "unclassified"), tax.get("g", "unclassified"),
            tax.get("s", "unclassified"))


def load_v4(ref):
    cul = read_tsv(os.path.join(ref, "cultures.tsv"))
    species = sorted({r["organism"] for r in cul if r["named"] == "True" and r.get("name_status", "ok") == "ok"})
    genera = sorted({v4_genus(r["organism"].split()[0]) for r in cul})
    by_stem = collections.defaultdict(list)          # epithet stem -> v4 species
    for sp in species:
        g, ep = sp.split()[0], sp.split()[1]
        by_stem[stem(ep)].append(sp)
    return species, genera, by_stem


def load_v18(meta_path):
    meta = [r for r in read_tsv(meta_path) if r["is_amf"] == "1"]
    by_label = collections.defaultdict(list)
    by_genus = collections.defaultdict(list)
    for r in meta:
        by_label[r["label"]].append(r["ref_id"])
        by_genus[r["genus"]].append(r["ref_id"])
    return meta, by_label, by_genus


def v18_species_labels(genus, ep, by_label):
    """V18 labels for a species: same genus (or synonym) and same epithet stem; compound labels
    ('Entrophospora Albahypha drummondi') match on their first and last words."""
    gs = {genus, GENUS_SYN.get(genus, genus)} | {k for k, v in GENUS_SYN.items() if v == genus}
    out = []
    for lab in by_label:
        w = lab.split()
        if w[0] in gs and stem(w[-1]) == stem(ep):
            out.append(lab)
    return sorted(out)


def suggest(key, n, v4, v18):
    """One names.tsv row with a suggested mapping; 'review' says whether a decision is needed."""
    species, genera, by_stem = v4
    _, by_label, by_genus = v18
    o, f, g, s = key
    row = {"euk_order": o, "euk_family": f, "euk_genus": g, "euk_species": s, "reads": n,
           "level": "", "v4_organism": "", "v4_status": "", "v4_genus": "", "genus_resolved": "yes",
           "v18_labels": "", "v18_genus_records": "", "basis": "", "review": "no", "note": ""}
    cf = s.startswith("cf.")
    if g == "unclassified" or PLACEHOLDER.search(g):
        named_family = f != "unclassified" and not PLACEHOLDER.search(f)
        if g == "unclassified" and (o == "unclassified" or (named_family and f not in ("",))):
            # genus unknown inside a named family (could be a v4 genus), or nothing known at all
            row.update(level="unplaced", v4_status="excluded", basis="genus unknown",
                       note="no genus; the family or order holds named genera")
            if o != "unclassified" and (PLACEHOLDER.search(f) or PLACEHOLDER.search(o)):
                row.update(v4_status="genus not in v4", level="placeholder", basis="placeholder taxon", note="")
            return row
        deepest = g if g != "unclassified" else (f if f != "unclassified" else o)
        row.update(level="placeholder", v4_organism=f"{deepest} sp.", v4_status="genus not in v4",
                   basis="placeholder taxon")
        return row
    vg = v4_genus(g)
    held = vg in genera
    row["v4_genus"] = vg if held else ""
    row["v18_genus_records"] = len(by_genus.get(vg, [])) + (len(by_genus.get(g, [])) if g != vg else 0)
    if s == "unclassified" or cf:
        row.update(level="genus", v4_organism=f"{g} sp.", v4_status="genus in v4" if held else "genus not in v4",
                   basis="genus name" + (" (cf. name read as genus only)" if cf else ""))
        return row
    row["level"] = "species"
    hits = by_stem.get(stem(s), [])
    same_genus = [sp for sp in hits if v4_genus(sp.split()[0]) == vg]
    other_genus = [sp for sp in hits if v4_genus(sp.split()[0]) != vg]
    if same_genus:
        row.update(v4_organism=same_genus[0], v4_status="species in v4",
                   basis="exact" if same_genus[0] == f"{g} {s}" else "genus synonym or epithet ending")
    else:
        labs = v18_species_labels(vg, s, by_label) or v18_species_labels(g, s, by_label)
        # name the species as V18 spells it when V18 holds exactly one label for it (Rhizophagus invermaius)
        org = labs[0] if len(labs) == 1 and len(labs[0].split()) == 2 else f"{vg} {s}"
        row.update(v4_organism=org, v4_status="species not in v4, genus in v4" if held else "genus not in v4",
                   basis="genus name", v18_labels=";".join(labs))
    if other_genus:
        row.update(review="yes", note=f"epithet also in v4 under another genus: {', '.join(other_genus)}")
    return row


def apply_decisions(rows, path):
    """Decisions (names_decisions.tsv) replace the suggestion for a label; every decision must hit a label."""
    dec = {(d["euk_genus"], d["euk_species"]): d for d in read_tsv(path)}
    used = set()
    for r in rows:
        d = dec.get((r["euk_genus"], r["euk_species"]))
        if not d:
            continue
        used.add((r["euk_genus"], r["euk_species"]))
        for k in ("level", "v4_organism", "v4_status", "v4_genus", "genus_resolved", "v18_labels"):
            if d.get(k, "") != "":
                r[k] = d[k] if d[k] != "-" else ""
        r.update(basis="decision", review="decided", note=d["reason"])
    unused = set(dec) - used
    if unused:
        sys.exit(f"decisions that match no label: {sorted(unused)}")
    left = [r for r in rows if r["review"] == "yes"]
    if left:
        sys.exit(f"labels still to decide: {[(r['euk_genus'], r['euk_species']) for r in left]}")
    return rows


def cmd_names(a):
    hdrs = read_fasta(f"{a.windows}.full_its.fasta", full_header=True)
    counts = collections.Counter(label_key(parse_header(h)[1]) for h in hdrs)
    v4 = load_v4(a.ref)
    v18 = load_v18(a.v18_meta)
    rows = [suggest(k, n, v4, v18) for k, n in sorted(counts.items())]
    # genera that hold a v4 species only under another genus name are decisions, not lookups; so are
    # families that hold a v4 genus but whose EUKARYOME genus is not a v4 genus name (Dentiscutataceae: Fuscutata)
    v4_families = {}
    for r in rows:
        if r["level"] == "genus" and r["v4_status"] == "genus not in v4":
            sibs = [x for x in rows if x["euk_genus"] == r["euk_genus"] and x["review"] == "yes"]
            if sibs:
                r.update(review="yes", note="a species of this genus matches a v4 species of another genus")
    fam_genera = collections.defaultdict(set)
    for r in rows:
        if r["level"] in ("genus", "species"):
            fam_genera[r["euk_family"]].add(r["euk_genus"])
    v4_genera = set(v4[1])
    for r in rows:
        if r["level"] in ("genus", "species") and r["v4_status"] == "genus not in v4" and r["review"] == "no":
            fam = r["euk_family"]
            stem_fam = fam[:-6] if fam.endswith("aceae") else fam   # Dentiscutataceae -> Dentiscut...
            if any(g.startswith(stem_fam[:6]) for g in v4_genera) and not any(
                    v4_genus(g) in v4_genera for g in fam_genera[fam]):
                r.update(review="yes", note=f"family {fam} is named for a v4 genus that EUKARYOME does not use here")
    if a.decisions:
        rows = apply_decisions(rows, a.decisions)
    cols = ["euk_order", "euk_family", "euk_genus", "euk_species", "reads", "level", "v4_organism", "v4_status",
            "v4_genus", "genus_resolved", "v18_labels", "v18_genus_records", "basis", "review", "note"]
    write_tsv(a.out, rows, cols)
    st = collections.Counter((r["level"], r["v4_status"]) for r in rows)
    rd = collections.Counter()
    for r in rows:
        rd[(r["level"], r["v4_status"])] += int(r["reads"])
    print(f"{len(hdrs)} reads, {len(rows)} labels -> {a.out}")
    for k in sorted(st):
        print(f"  {k[0]:12s} {k[1]:32s} labels {st[k]:4d}  reads {rd[k]:5d}")
    print(f"  to review: {sum(r['review'] == 'yes' for r in rows)} labels")


# ---------------------------------------------------------------- roles
STATUS_ROLE = {"species in v4": "known species", "species not in v4, genus in v4": "novel species, genus in v4",
               "genus not in v4": "novel lineage, genus not in v4", "genus in v4": "genus only",
               "genus open": "genus only"}   # membership of the genus in v4 unresolved: scored, no outcome


def stratum(nm):
    lv, st = nm["level"], nm["v4_status"]
    if st == "excluded":
        return "unplaced"
    if lv == "placeholder":
        return "placeholder lineage"
    if lv == "species":
        return {"species in v4": "known species", "species not in v4, genus in v4": "absent species, genus in v4",
                "genus not in v4": "absent genus, named species"}[st]
    if st == "genus in v4":
        return "genus only, genus in v4" if nm.get("genus_resolved", "yes") != "no" else "genus only, membership open"
    if st == "genus open":
        return "genus only, membership open"
    return "absent genus, genus only"


def cmd_roles(a):
    t0 = time.time()
    names = {(r["euk_order"], r["euk_family"], r["euk_genus"], r["euk_species"]): r for r in read_tsv(a.names)}
    rec = {parse_header(r["read_id"])[0]: r for r in read_tsv(f"{a.windows}.records.tsv")}
    its = {parse_header(h)[0]: (s, h) for h, s in read_fasta(f"{a.windows}.full_its.fasta", full_header=True).items()}
    lsu = {parse_header(h)[0]: s for h, s in read_fasta(f"{a.windows}.LSU_LROR-FLR2.fasta", full_header=True).items()}
    ssu = {parse_header(h)[0]: s for h, s in read_fasta(f"{a.windows}.ssu_flank.fasta", full_header=True).items()}
    v4_its = read_fasta(os.path.join(a.ref, "fasta", "full_its.fasta"))
    v4_lsu = read_fasta(os.path.join(a.ref, "fasta", "LSU_LROR-FLR2.fasta"))
    v4_ssu = read_fasta(os.path.join(a.ref, "fasta", "ssu_flank.fasta"))
    copies = {r["accession"]: r for r in read_tsv(os.path.join(a.ref, "copies.tsv"))}
    v4_species = set(load_v4(a.ref)[0])
    its_index = collections.defaultdict(list)
    for acc, s in v4_its.items():
        its_index[s].append(acc)

    missing = [h for h in rec if h not in its]
    unknown = sorted({label_key(parse_header(its[h][1])[1]) for h in its} - set(names))
    if unknown:
        sys.exit(f"labels missing from {a.names}: {unknown[:5]} ...")
    for k, nm in names.items():
        if nm["v4_status"] in ("species in v4",) and nm["v4_organism"] not in v4_species:
            sys.exit(f"names.tsv maps {k} to {nm['v4_organism']}, which is not a v4 species")

    # the four windows the frozen method reads; a read without one of them cannot be scored
    wins = {w: {parse_header(hh)[0]: s for hh, s in read_fasta(f"{a.windows}.{w}.fasta", full_header=True).items()}
            for w in METHOD_WINDOWS}
    rows = []
    for h, (seq, header) in sorted(its.items()):
        tax = parse_header(header)[1]
        nm = names[label_key(tax)]
        role = STATUS_ROLE.get(nm["v4_status"], "excluded")
        row = {"asv": h, "organism": nm["v4_organism"], "code": LAB, "lab": LAB, "collection": COLLECTION,
               "role": role, "stratum": stratum(nm), "euk_label": ";".join(f"{k}__{tax.get(k, '')}" for k in "ofgs"),
               "length": len(seq), "non_euk_id": "yes" if not h.startswith("EUK") else "no",
               "identical_v4": "", "identical_v4_all3": "", "v18_best": "", "v18_best_pct": "",
               "identical_lgdb": "", "exclusion": ""}
        if role == "excluded":
            row["exclusion"] = nm["note"] or "unplaced"
        lacking = [w for w in METHOD_WINDOWS if not wins[w].get(h)]
        if lacking and role != "excluded":
            row.update(exclusion=f"no Stage 1 window: {', '.join(lacking)}", role_before_exclusion=role,
                       role="excluded")
            rows.append(row)
            continue
        same = its_index.get(seq, [])
        if same and role != "excluded":
            all3 = [x for x in same
                    if h in lsu and x in v4_lsu and dist(lsu[h], v4_lsu[x], 0) == 0
                    and h in ssu and x in v4_ssu and dist(ssu[h], v4_ssu[x], 0) == 0]
            row["identical_v4"] = ";".join(sorted({copies[x]["culture"] for x in same if x in copies}))
            row["identical_v4_all3"] = "yes" if all3 else "no"
            row["exclusion"] = "identical in full ITS to a v4 copy"
            row["role_before_exclusion"] = role
            row["role"] = "excluded"
        rows.append(row)

    # V18: best record in LSU (the LSU arm may see the read's own record)
    if a.v18:
        v18 = {hh: s.strip("N") for hh, s in read_fasta(a.v18).items()}
        v18_sp = {hh: " ".join(hh.split("_")[2:4]) if len(hh.split("_")) >= 4 else hh for hh in v18}
        for i, row in enumerate(rows):
            q = lsu.get(row["asv"])
            if not q or row["role"] not in ROLES:      # excluded reads (some without a name) are not flagged
                continue
            g = row["organism"].split()[0] if row["organism"] else ""
            order = sorted(v18, key=lambda r: v18_sp[r].split()[0] != g)
            best, bd = None, None
            for y in order:
                k = -1 if bd is None else int(bd * min(len(q), len(v18[y])) / 100) + 1
                d = dist(q, v18[y], k)
                if d is not None and (bd is None or d < bd):
                    best, bd = y, d
            if best:
                row.update(v18_best=v18_sp[best], v18_best_pct=round(100 - bd, 2))
            if (i + 1) % 1000 == 0:
                print(f"  V18: {i + 1} reads ({time.time() - t0:.0f} s)", flush=True)

    # LongGloDB v1 (public files): reads identical in full ITS to an ASV, counted by lab
    lg_counts = collections.Counter()
    if a.longglodb:
        f = sorted(glob.glob(os.path.join(a.longglodb, "LongGloDB_ITS_*.fasta")))
        if not f:
            sys.exit(f"no LongGloDB_ITS_*.fasta in {a.longglodb}")
        lg = collections.defaultdict(set)
        for hh, s in read_fasta(f[0]).items():
            m = re.match(r"^.+_(?P<lab>[^_]+)_ASV\d+_[^_]+$", hh)
            lg[s].add(m["lab"] if m else "?")
        for row in rows:
            labs = lg.get(its[row["asv"]][0])
            if labs:
                row["identical_lgdb"] = ";".join(sorted(labs))
                for lab in labs:
                    lg_counts[lab] += 1

    acols = ["asv", "organism", "code", "lab", "collection", "role", "stratum", "euk_label", "length", "non_euk_id",
             "identical_v4", "identical_v4_all3", "v18_best", "v18_best_pct", "identical_lgdb", "exclusion",
             "role_before_exclusion"]
    write_tsv(f"{a.out}.asvs.tsv", rows, acols)
    # the test reads' windows under their read IDs, for predict_longglodb.py (--test-fasta {out}.win) and the
    # runner; EUKARYOME sequences stay out of git, like LongGloDB's
    test_ids = [r["asv"] for r in rows if r["role"] in ROLES]
    vt_win = {parse_header(hh)[0]: s for hh, s in
              read_fasta(f"{a.windows}.SSU_VT_NS31-AM1.fasta", full_header=True).items()}
    for w, seqs in list(wins.items()) + [("SSU_VT_NS31-AM1", vt_win)]:
        with open(f"{a.out}.win.{w}.fasta", "w") as fh:
            for h in test_ids:
                if seqs.get(h):
                    fh.write(f">{h}\n{seqs[h]}\n")

    # one 'culture' per organism (EUKARYOME has no culture codes), as lgdb_test.py groups them
    cul = collections.OrderedDict()
    for row in rows:
        if row["role"] not in ROLES:
            continue
        c = cul.setdefault((row["organism"], row["code"]), {"organism": row["organism"], "code": row["code"],
                                                            "role": row["role"], "asvs": 0, "strata": set(),
                                                            "v18_best_pct": ""})
        if c["role"] != row["role"]:
            sys.exit(f"{row['organism']} has reads in two roles")
        c["asvs"] += 1
        c["strata"].add(row["stratum"])
        if row["v18_best_pct"] != "" and (c["v18_best_pct"] == "" or row["v18_best_pct"] > c["v18_best_pct"]):
            c["v18_best_pct"] = row["v18_best_pct"]
    crow = [{**c, "strata": ";".join(sorted(c["strata"]))} for c in cul.values()]
    write_tsv(f"{a.out}.cultures.tsv", crow, ["organism", "code", "role", "asvs", "strata", "v18_best_pct"])

    # V18 removal plan: species rounds (absent species, genus in v4) and genus rounds (absent named genera)
    meta, by_label, by_genus = load_v18(a.v18_meta)
    plan, apart = [], []
    test = [r for r in rows if r["role"] in ROLES]
    for org in sorted({r["organism"] for r in test if r["role"] == "novel species, genus in v4"}):
        nm = next(n for n in names.values() if n["v4_organism"] == org)
        labs = [x for x in nm["v18_labels"].split(";") if x]
        ids = sorted(i for lab in labs for i in by_label.get(lab, []))
        n = sum(r["organism"] == org for r in test)
        (plan if ids else apart).append({"round": "species", "target": org, "v18_labels": ";".join(labs),
                                         "records": len(ids), "ref_ids": ";".join(ids), "reads": n})
    genus_reads = collections.Counter(r["organism"].split()[0] for r in test
                                      if r["role"] == "novel lineage, genus not in v4"
                                      and r["stratum"] != "placeholder lineage")
    for g, n in sorted(genus_reads.items()):
        ids = sorted(by_genus.get(g, []))
        (plan if ids else apart).append({"round": "genus", "target": g, "v18_labels": ";".join(
            sorted({m["label"] for m in meta if m["genus"] == g})), "records": len(ids), "ref_ids": ";".join(ids),
            "reads": n})
    n_ph = sum(r["stratum"] == "placeholder lineage" for r in test)
    apart.append({"round": "none", "target": "placeholder lineages", "v18_labels": "", "records": 0, "ref_ids": "",
                  "reads": n_ph})
    pcols = ["round", "target", "v18_labels", "records", "reads", "ref_ids"]
    write_tsv(f"{a.out}.v18_plan.tsv", plan, pcols)
    write_tsv(f"{a.out}.v18_apart.tsv", apart, pcols)

    # ---- log
    print(f"reads with Stage 1 windows: {len(its)}; in records.tsv without windows: {len(missing)} {missing}")
    ex = collections.Counter(r["exclusion"] for r in rows if r["role"] == "excluded")
    print(f"excluded: {sum(ex.values())}  " + "; ".join(f"{k} {v}" for k, v in ex.most_common()))
    idn = [r for r in rows if r["identical_v4"]]
    print(f"  identical in full ITS to a v4 copy: {len(idn)} reads ({sum(r['identical_v4_all3'] == 'yes' for r in idn)} "
          f"also identical in LROR-FLR2 and SSU flank); by role before removal: "
          + "; ".join(f"{k} {v}" for k, v in collections.Counter(r['role_before_exclusion'] for r in idn).most_common()))
    for r in idn:
        print(f"    {r['asv']}  {r['organism']}  ->  {r['identical_v4']}")
    print("test set by role and stratum (reads / taxa):")
    for role in ROLES:
        rs = [r for r in test if r["role"] == role]
        print(f"  {role:34s} {len(rs):5d} / {len({r['organism'] for r in rs})}")
        for s_, n in collections.Counter(r["stratum"] for r in rs).most_common():
            print(f"      {s_:38s} {n:5d} / {len({r['organism'] for r in rs if r['stratum'] == s_})}")
    print(f"non-EUKARYOME IDs in the test: {sum(r['non_euk_id'] == 'yes' for r in test)}")
    if a.v18:
        hi = [r for r in test if r["v18_best_pct"] != "" and r["v18_best_pct"] >= 99.5]
        print(f"V18 record >= 99.5% in LSU: {len(hi)} reads; by role: "
              + "; ".join(f"{k} {v}" for k, v in collections.Counter(r['role'] for r in hi).most_common()))
    if a.longglodb:
        print(f"identical in full ITS to a LongGloDB v1 ASV: {sum(bool(r['identical_lgdb']) for r in rows)} reads; "
              f"by LongGloDB lab: {dict(lg_counts)}")
    print(f"V18 removal plan: {sum(p['round'] == 'species' for p in plan)} species rounds, "
          f"{sum(p['round'] == 'genus' for p in plan)} genus rounds ({sum(p['reads'] for p in plan)} reads); "
          f"reported under the intact LSU reference only: "
          + "; ".join(f"{p['target']} ({p['reads']})" for p in apart))
    print(f"inputs: names.tsv {sha(a.names)[:16]}, records {sha(a.windows + '.records.tsv')[:16]}, "
          f"full_its {sha(a.windows + '.full_its.fasta')[:16]}")
    print(f"files: {a.out}.asvs.tsv, .cultures.tsv, .v18_plan.tsv, .v18_apart.tsv ({time.time() - t0:.0f} s)")


# ---------------------------------------------------------------- genus-level prediction
def cmd_genus_inside(a):
    """The complex rule one rank up. A read of a genus v4 does not hold can still be given a v4 genus through
    a v4 culture named only to genus (Glomus sp., Septoglomus sp., ...), whose genus may be EUKARYOME's
    absent genus under an older name. A read is 'inside' such a culture when its nearest main-type copy of
    the culture lies within T on all three markers (T as in the predictions). Genera with reads inside a
    culture named only to genus get no genus outcome; their false-known outcome is unchanged."""
    T = {k: float(v) for k, v in (x.split("=") for x in a.thresholds.split(","))}
    asv = [r for r in read_tsv(f"{a.overlap}.asvs.tsv") if r["role"] == "novel lineage, genus not in v4"]
    win = {lab: read_fasta(f"{a.overlap}.win.{m}.fasta") for m, lab in MARKERS}
    cul = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    unnamed = {c for c, r in cul.items() if r["named"] != "True"}
    copies = [r for r in read_tsv(os.path.join(a.ref, "copies.tsv"))
              if r["status"] == "retained" and r["copy_class"] == "main" and r["culture"] in unnamed]
    seq = {m: read_fasta(os.path.join(a.ref, "fasta", f"{m}.fasta")) for m, _ in MARKERS}
    uniq = {(m, c): sorted({seq[m][r["accession"]] for r in copies if r["culture"] == c and r["accession"] in seq[m]})
            for m, _ in MARKERS for c in unnamed}
    rows = []
    for r in asv:
        inside, near = [], {}
        for c in sorted(unnamed):
            ds = {}
            for m, lab in MARKERS:
                q = win[lab].get(r["asv"])
                vals = [dist(q, s, int(30 * min(len(q), len(s)) / 100) + 1) for s in uniq[(m, c)]] if q else []
                vals = [v for v in vals if v is not None]
                ds[lab] = min(vals) if vals else None
            if all(ds[lab] is not None and ds[lab] <= T[lab] for _, lab in MARKERS):
                inside.append(c)
            if ds["ITS"] is not None and (not near or ds["ITS"] < near["ITS"]):
                near = {"culture": c, **ds}
        rows.append({"asv": r["asv"], "organism": r["organism"], "stratum": r["stratum"],
                     "inside_genus_only_culture": ";".join(inside),
                     "nearest_by_ITS": near.get("culture", ""),
                     **{f"{lab}_pct": "" if near.get(lab) is None else round(near[lab], 2) for _, lab in MARKERS}})
    write_tsv(f"{a.out}.reads.tsv", rows, ["asv", "organism", "stratum", "inside_genus_only_culture",
                                           "nearest_by_ITS", "ITS_pct", "LSU_pct", "SSU_pct"])
    by_genus = collections.defaultdict(collections.Counter)
    for r in rows:
        g = r["organism"].split()[0]
        by_genus[g]["reads"] += 1
        for c in r["inside_genus_only_culture"].split(";"):
            if c:
                by_genus[g][c] += 1
    out = []
    for g, cn in sorted(by_genus.items()):
        inside = {k: v for k, v in cn.items() if k != "reads"}
        if inside:
            out.append({"euk_genus": g, "reads": cn["reads"],
                        "inside": "; ".join(f"{k} {v}" for k, v in sorted(inside.items(), key=lambda x: -x[1])),
                        "genus_outcome": "none (a v4 culture named only to genus may hold this genus)"})
    write_tsv(f"{a.out}.genera.tsv", out, ["euk_genus", "reads", "inside", "genus_outcome"])
    print(f"T: {T}; {len(rows)} absent-lineage reads against {len(unnamed)} v4 cultures named only to genus")
    for o in out:
        print(f"  {o['euk_genus']:28s} {o['reads']:5d} reads; inside: {o['inside']}")
    print(f"files: {a.out}.reads.tsv, {a.out}.genera.tsv")


# ---------------------------------------------------------------- VT geometry
def cmd_vtgeom(a):
    rec = {parse_header(r["read_id"])[0]: r for r in read_tsv(f"{a.windows}.records.tsv")}
    ssu = {parse_header(h)[0]: s for h, s in read_fasta(f"{a.windows}.ssu_flank.fasta", full_header=True).items()}
    vt = read_fasta(os.path.join(a.ref, "fasta", "SSU_VT_NS31-AM1.fasta"))
    prefixes = sorted({s[:80] for s in vt.values()})
    offs = {}
    for h, s in ssu.items():
        st = rec[h]["SSU_VT_NS31-AM1_status"]
        if st == "ok":
            offs[h] = 0
            continue
        best = None
        for p in prefixes:
            al = edlib.align(s[:30], p, mode="HW", task="locations", k=6)
            if al["editDistance"] >= 0 and (best is None or al["editDistance"] < best[0]):
                best = (al["editDistance"], al["locations"][0][0])
        offs[h] = best[1] if best else None
    c = collections.Counter(offs.values())
    rows = [{"offset_nt": k if k is not None else "not found", "reads": v} for k, v in
            sorted(c.items(), key=lambda x: (x[0] is None, x[0] if x[0] is not None else 0))]
    write_tsv(f"{a.out}.offsets.tsv", rows, ["offset_nt", "reads"])
    print(f"{len(offs)} reads: offset of the read start from the 3' end of the NS31 site (0 = whole NS31-AM1 "
          f"fragment present):")
    for r in rows:
        print(f"  {r['offset_nt']:>9}  {r['reads']}")
    if a.maarjam:
        # Where the NS31-AM1 fragment begins in each MaarjAM reference: align the first 40 nt of the v4
        # fragments (just after NS31) to the reference; u = reference bases before that point.
        mj = read_fasta(a.maarjam)
        heads = sorted({s[:40] for s in vt.values()})
        geo = []
        for name, s in mj.items():
            best = None
            for p in heads:
                al = edlib.align(p, s, mode="HW", task="locations", k=8)
                if al["editDistance"] >= 0 and (best is None or al["editDistance"] < best[0]):
                    best = (al["editDistance"], al["locations"][0][0])
            n_n = s.count("N")
            geo.append((name, len(s), best[1] if best else None, n_n))
        print(f"\nMaarjAM: {len(geo)} references; fragment start found in {sum(g[2] is not None for g in geo)}")
        lens = sorted(g[1] for g in geo)
        print(f"  length median {lens[len(lens) // 2]}, 5-95% {lens[len(lens) // 20]}-{lens[19 * len(lens) // 20]}")
        out = []
        for off in sorted({o for o in offs.values() if o is not None}):
            ok = 0
            for _, ln, u, n_n in geo:
                if u is None:
                    continue
                cover = (ln - u - off) / max(1, ln - n_n)   # best case: aligned to the reference's 3' end
                ok += cover >= 0.95
            n = sum(1 for g in geo if g[2] is not None)
            out.append({"offset_nt": off, "reads": c[off], "maarjam_refs_coverable_95": ok,
                        "share": round(ok / n, 4) if n else ""})
            print(f"  read offset {off:3d} nt ({c[off]} reads): {ok} of {n} references ({100 * ok / n:.1f}%) "
                  f"coverable at >= 95%")
        write_tsv(f"{a.out}.maarjam.tsv", out, ["offset_nt", "reads", "maarjam_refs_coverable_95", "share"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("names", "roles", "vtgeom", "genus_inside"):
        p = sub.add_parser(name)
        if name == "genus_inside":
            p.add_argument("--overlap", required=True, help="prefix of the roles run")
            p.add_argument("--thresholds", required=True, help="T as passed to predict_longglodb.py")
        else:
            p.add_argument("--windows", required=True, help="Stage 1 window prefix (.../windows/euk_test)")
        p.add_argument("--ref", required=True, help="culture reference v4 folder")
        p.add_argument("--out", required=True)
        if name in ("names", "roles"):
            p.add_argument("--v18-meta", default=V18_META, help="V18 reference_metadata.tsv of the frozen package")
        if name == "names":
            p.add_argument("--decisions", help="names_decisions.tsv: reviewed decisions; with it, the output is names.tsv")
        if name == "roles":
            p.add_argument("--names", required=True, help="the fixed names.tsv")
            p.add_argument("--v18", help="V18 AMF-only FASTA (flag: best V18 record in LSU)")
            p.add_argument("--longglodb", help="LongGloDB v1 folder (count reads identical in ITS to an ASV)")
        if name == "vtgeom":
            p.add_argument("--maarjam", help="MaarjAM FASTA of the v4 crosswalk")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    {"names": cmd_names, "roles": cmd_roles, "vtgeom": cmd_vtgeom, "genus_inside": cmd_genus_inside}[a.cmd](a)


if __name__ == "__main__":
    main()
