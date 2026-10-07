#!/usr/bin/env python3
"""Predictions for the LongGloDB test, fixed before the frozen pipeline is run on it.

The paper's complex rule (species_groups.py) says two named species are inseparable when, on every
marker (full ITS, LROR-FLR2, SSU flank), their nearest copies lie within the within-species spread T
(the 99th percentile of the distance from a main-type copy to the nearest main-type copy of another
culture of the same species). The paper's claim is that false-known calls arise only there.

Applied to cultures the reference has never seen, the rule makes a prediction that can fail:
  * a test ASV can be wrongly called a known species only if it lies within T of that species on all
    three markers ("inside"). The set of such ASVs is the PREDICTED set of possible false knowns.
  * any false-known call outside that set, once the pipeline is run, contradicts the claim.
For cultures of species v4 already holds, the same distances say whether the ASV lies inside its own
species' spread (a known call is possible) and whether it also lies inside another species (a group call).

T is recomputed from v4 exactly as species_groups.py does and printed next to the paper's values
(ITS 6.0, LSU 3.0, SSU 0.8); --thresholds fixes them instead. Distances: edlib infix edit distance as %
of the shorter sequence, against the deduplicated main-type copies of every named v4 species.

Region check first. LongGloDB's SSU, ITS and LSU files are cut by its authors, not by Stage 1. The
Stefani cultures reprocessed inside LongGloDB are the same molecules as v4 copies, so their distance to
their own v4 culture on each marker shows whether the files cover comparable stretches (it should be
about 0%). If it isn't, run Stage 1 on the test sequences and pass its windows with --test-fasta.

Writes {out}.predictions.tsv (one row per test ASV), {out}.cultures.tsv (per culture), and
{out}.lock (sha256 of the predictions with T, the v4 fingerprint and the input hashes) so the
predictions can be shown to predate the run.

Usage (from amfsplit_stage3, after longglodb_overlap.py):
  python3 predict_longglodb.py --ref ../ref/culture_ref/v4 --overlap ../bench/longglodb/overlap \\
      --longglodb ../ref/glodb --out ../bench/longglodb/predictions
"""
import argparse
import collections
import datetime
import glob
import hashlib
import json
import os
import re
import statistics
import sys
import time

import edlib

from arms import is_named, read_fasta, read_tsv, write_tsv

MARKERS = [("full_its", "ITS", "its"), ("LSU_LROR-FLR2", "LSU", "lsu"), ("ssu_flank", "SSU", "ssu")]
PAPER_T = {"ITS": 6.0, "LSU": 3.0, "SSU": 0.8}
CAP = 30.0  # distances above this are recorded as > CAP (only the inside test and the nearest species need them)


def dist(a, b, cap=None):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    k = -1 if cap is None else int(cap * len(q) / 100) + 1
    d = edlib.align(q, t, mode="HW", task="distance", k=k)["editDistance"]
    return None if d < 0 else 100 * d / len(q)


def norm(name):
    t = name.replace("_", " ").split()
    if not t:
        return ""
    g = t[0].strip("[]")
    g = {"Rhizoglomus": "Rhizophagus"}.get(g, g)
    return " ".join([g] + t[1:2])


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def quantile(vals, q):
    v = sorted(vals)
    pos = (len(v) - 1) * q / 100.0
    lo = int(pos)
    return v[lo] + (v[min(lo + 1, len(v) - 1)] - v[lo]) * (pos - lo)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="culture reference (v4)")
    ap.add_argument("--overlap", required=True, help="prefix of the longglodb_overlap.py run")
    ap.add_argument("--longglodb", help="LongGloDB folder, for the region check on the reprocessed Stefani part")
    ap.add_argument("--test-fasta", help="prefix of Stage 1 windows for the test ASVs ({p}.full_its.fasta, ...) "
                                         "instead of the LongGloDB SSU/ITS/LSU files")
    ap.add_argument("--quantile", type=float, default=99.0)
    ap.add_argument("--thresholds", help="fix T instead of recomputing it, e.g. ITS=6.0,LSU=3.0,SSU=0.8")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()

    # ---- v4: named species, deduplicated main-type copies per marker
    cul = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    sp_of = {c: m["organism"] for c, m in cul.items() if is_named(m["organism"]) and m.get("name_status", "ok") == "ok"}
    copies = read_tsv(os.path.join(a.ref, "copies.tsv"))
    mains = [r for r in copies if r["status"] == "retained" and r["copy_class"] == "main" and r["culture"] in sp_of]
    by_sp = collections.defaultdict(list)
    for r in mains:
        by_sp[sp_of[r["culture"]]].append(r)
    species = sorted(by_sp)
    seq = {m: read_fasta(os.path.join(a.ref, "fasta", f"{m}.fasta")) for m, _, _ in MARKERS}
    uniq = {(m, sp): sorted({(seq[m][r["accession"]], r["culture"]) for r in by_sp[sp] if r["accession"] in seq[m]})
            for m, _, _ in MARKERS for sp in species}
    norm_sp = {norm(sp): sp for sp in species}
    genera = {norm(sp).split()[0] for sp in species} | {norm(c.split(" | ")[0]).split()[0] for c in cul}

    # ---- T, as species_groups.py
    if a.thresholds:
        T = {k: float(v) for k, v in (x.split("=") for x in a.thresholds.split(","))}
        print("within-species spread T fixed: " + ", ".join(f"{k} {v:.2f}%" for k, v in T.items()))
    else:
        T = {}
        for m, lab, _ in MARKERS:
            vals = []
            for sp in species:
                items = uniq[(m, sp)]
                if len({c for _, c in items}) < 2:
                    continue
                for s1, c1 in items:
                    vals.append(min(dist(s1, s2) for s2, c2 in items if c2 != c1))
            T[lab] = quantile(vals, a.quantile)
        print(f"within-species spread T ({a.quantile:g}th pct, recomputed from v4): "
              + ", ".join(f"{lab} {T[lab]:.2f}% (paper {PAPER_T[lab]})" for _, lab, _ in MARKERS))
        off = [lab for lab in T if abs(T[lab] - PAPER_T[lab]) > 0.05]
        if off:
            print(f"  WARNING: {', '.join(off)} differ from the paper's values; check the reference version, "
                  f"or pass --thresholds to use the paper's")

    # ---- test ASVs and their regions
    arows = read_tsv(f"{a.overlap}.asvs.tsv")
    roles = {r["asv"]: r for r in arows}
    test_roles = {"known species", "novel species, genus in v4", "novel lineage, genus not in v4", "genus only"}
    if a.test_fasta:
        tseq = {lab: read_fasta(f"{a.test_fasta}.{m}.fasta") for m, lab, _ in MARKERS}
    else:
        tseq = {lab: read_fasta(f"{a.overlap}.test_{short}.fasta") for _, lab, short in MARKERS}
    test = [r["asv"] for r in arows if r["role"] in test_roles]
    print(f"test: {len(test)} ASVs in {len({(roles[h]['organism'], roles[h]['code']) for h in test})} cultures; "
          + ", ".join(f"{lab} {sum(h in tseq[lab] for h in test)}" for _, lab, _ in MARKERS))

    # ---- region check on the reprocessed Stefani part
    if a.longglodb and not a.test_fasta:
        lg = {}
        for _, lab, short in MARKERS:
            f = sorted(glob.glob(os.path.join(a.longglodb, f"LongGloDB_{short.upper()}_*.fasta")))
            lg[lab] = read_fasta(f[0]) if f else {}
        acc_by_cul = collections.defaultdict(list)
        for r in copies:
            if r["status"] == "retained" and r["copy_class"] == "main":
                acc_by_cul[r["culture"]].append(r["accession"])
        stf = [r for r in arows if r["role"] == "Stefani reprocessed" and r["nearest_its_pct"] in ("0", "0.0")]
        print(f"\nregion check: {len(stf)} reprocessed Stefani ASVs identical in ITS to a v4 copy, against that culture")
        for m, lab, _ in MARKERS:
            ds = []
            for r in stf:
                q = lg[lab].get(r["asv"])
                accs = [x for x in acc_by_cul.get(r["nearest_v4"], []) if x in seq[m]]
                if q and accs:
                    dd = [dist(q, seq[m][x], CAP) for x in accs]
                    ds.append(min(CAP if d is None else d for d in dd))
            if ds:
                print(f"  {lab}: median {statistics.median(ds):.2f}%, 95th {quantile(ds, 95):.2f}%, "
                      f"share <= 0.5%: {sum(d <= 0.5 for d in ds) / len(ds):.0%} (n {len(ds)})")
        print("  near 0% means LongGloDB's files and the v4 regions cover comparable stretches")

    # ---- distances from every test ASV to every named v4 species, per marker
    rows = []
    for i, h in enumerate(test):
        r = roles[h]
        own = norm_sp.get(norm(r["organism"]))
        d_by = {}
        for m, lab, _ in MARKERS:
            q = tseq[lab].get(h)
            for sp in species:
                if not q or not uniq[(m, sp)]:
                    d_by[(lab, sp)] = None
                    continue
                ds = [dist(q, s, CAP) for s, _ in uniq[(m, sp)]]
                ds = [d for d in ds if d is not None]
                d_by[(lab, sp)] = min(ds) if ds else None   # None = above CAP
        inside = [sp for sp in species
                  if all(d_by[(lab, sp)] is not None and d_by[(lab, sp)] <= T[lab] for _, lab, _ in MARKERS)]
        near = {}
        for _, lab, _ in MARKERS:
            cand = [(d_by[(lab, sp)], sp) for sp in species if d_by[(lab, sp)] is not None]
            near[lab] = min(cand) if cand else (None, "")
        role = r["role"]
        if role in ("novel species, genus in v4", "novel lineage, genus not in v4"):
            pred = "possible false known: inside " + "; ".join(inside) if inside else "not known (outside every v4 species)"
        elif role == "known species":
            others = [sp for sp in inside if sp != own]
            if own in inside:
                pred = "known possible" + (f" (group with {'; '.join(others)})" if others else "")
            else:
                pred = "outside own species' spread" + (f"; inside {'; '.join(others)}" if others else "")
        else:
            pred = "genus only: " + ("inside " + "; ".join(inside) if inside else "outside every v4 species")
        row = {"asv": h, "organism": r["organism"], "code": r["code"], "lab": r["lab"], "role": role,
               "own_species_in_v4": own or "", "prediction": pred, "inside": "; ".join(inside)}
        for _, lab, _ in MARKERS:
            d, sp = near[lab]
            row[f"{lab}_nearest_species"] = sp
            row[f"{lab}_nearest_pct"] = "" if d is None else round(d, 2)
            if own:
                row[f"{lab}_own_pct"] = "" if d_by[(lab, own)] is None else round(d_by[(lab, own)], 2)
        rows.append(row)
        if (i + 1) % 50 == 0:
            print(f"  {i + 1} ASVs ({time.time() - t0:.0f} s)")

    cols = ["asv", "organism", "code", "lab", "role", "own_species_in_v4", "prediction", "inside"]
    for _, lab, _ in MARKERS:
        cols += [f"{lab}_nearest_species", f"{lab}_nearest_pct", f"{lab}_own_pct"]
    for r in rows:
        for c in cols:
            r.setdefault(c, "")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    write_tsv(f"{a.out}.predictions.tsv", rows, cols)

    # ---- per culture
    cul_rows = []
    by_c = collections.defaultdict(list)
    for r in rows:
        by_c[(r["organism"], r["code"])].append(r)
    for (org, code), rs in sorted(by_c.items()):
        ins = collections.Counter(sp for r in rs for sp in r["inside"].split("; ") if sp)
        cul_rows.append({"organism": org, "code": code, "role": rs[0]["role"], "asvs": len(rs),
                         "asvs_inside_any": sum(bool(r["inside"]) for r in rs),
                         "inside_species": "; ".join(f"{sp} ({n})" for sp, n in ins.most_common()),
                         "own_inside": (sum(r["own_species_in_v4"] in r["inside"].split("; ") for r in rs)
                                        if rs[0]["own_species_in_v4"] else "")})
    write_tsv(f"{a.out}.cultures.tsv", cul_rows, list(cul_rows[0]))

    # ---- lock: the predictions, T and the inputs, hashed
    lock = {"created": datetime.datetime.now().isoformat(timespec="seconds"), "T": T, "quantile": a.quantile,
            "predictions_sha256": sha(f"{a.out}.predictions.tsv"), "cultures_sha256": sha(f"{a.out}.cultures.tsv"),
            "inputs": {p: sha(p) for p in [os.path.join(a.ref, "copies.tsv"), os.path.join(a.ref, "cultures.tsv"),
                                           f"{a.overlap}.asvs.tsv"] + [os.path.join(a.ref, "fasta", f"{m}.fasta")
                                                                        for m, _, _ in MARKERS]}}
    manifest = os.path.join(a.ref, "MANIFEST.json")
    if os.path.exists(manifest):
        try:
            lock["ref_fingerprint"] = json.load(open(manifest)).get("fingerprint", "")
        except ValueError:
            pass
    with open(f"{a.out}.lock", "w") as fh:
        json.dump(lock, fh, indent=1)

    # ---- report
    print("\npredictions (per ASV, then per culture)")
    for role in ("novel species, genus in v4", "novel lineage, genus not in v4", "known species", "genus only"):
        rs = [r for r in rows if r["role"] == role]
        cs = [c for c in cul_rows if c["role"] == role]
        if not rs:
            continue
        if role.startswith("novel"):
            print(f"  {role}: {len(rs)} ASVs; possible false known (inside a v4 species) {sum(bool(r['inside']) for r in rs)}; "
                  f"cultures with any such ASV {sum(c['asvs_inside_any'] > 0 for c in cs)} of {len(cs)}")
            for c in cs:
                if c["asvs_inside_any"]:
                    print(f"      {c['organism']} {c['code']}: {c['asvs_inside_any']}/{c['asvs']} ASVs inside {c['inside_species']}")
        elif role == "known species":
            own_in = sum(r["own_species_in_v4"] in r["inside"].split("; ") for r in rs)
            grp = sum(len([s for s in r["inside"].split("; ") if s]) > 1 for r in rs)
            print(f"  known species: {len(rs)} ASVs; inside own species {own_in}; inside another species too {grp}")
            for r in [r for r in rs if r["own_species_in_v4"] not in r["inside"].split("; ")][:8]:
                print(f"      outside own: {r['organism']} {r['code']} (ITS {r.get('ITS_own_pct', '')}%, "
                      f"LSU {r.get('LSU_own_pct', '')}%, SSU {r.get('SSU_own_pct', '')}% to own species)")
        else:
            print(f"  genus only: {len(rs)} ASVs; inside a v4 species {sum(bool(r['inside']) for r in rs)}")
    print(f"\nfalsifiable claim: no novel-species or novel-lineage ASV outside the 'inside' set is called a known "
          f"species once the frozen pipeline runs.\nfiles: {a.out}.predictions.tsv, {a.out}.cultures.tsv, "
          f"{a.out}.lock ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    sys.exit(main())
