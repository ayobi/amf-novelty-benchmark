#!/usr/bin/env python3
"""Phase 4: novelty decision tree on top of the Phase 3 reconciliation.

Every copy gets one status and a reason code, testing the cheapest explanation first:

  1. artifact     selection signature: against the nearest main-type copy of another culture (by SSU
                  flank), the SSU differs by >= --artifact-min-ssu % and by at least --artifact-ratio of
                  the full-ITS difference. Functional copies differ about ten times less in SSU than in
                  ITS; errors and pseudogenes diverge at a similar rate everywhere (the PX215117.1 case).
                  Reason A1.
  2. divergent    short-5.8S class, read from the copy itself: its 5.8S is >= --short-58s bp shorter than
                  the main-type 5.8S of the nearest culture (D1), or it sits in a paralog unit, an SH or VT
                  that other cultures reach only through divergent-class copies (D2). Placed by the VT arm
                  alone, since its SSU is the genome's.
                  A same-length 5.8S that differs by >= 3% from the nearest culture's main type is only
                  flagged (D3): such copies turned out to be second lineages in Septoglomus and co-dominant
                  types in Gigaspora, so they are not an explanation and go on through the tree.
  3. boundary     no VT call, but the best VT is within --boundary points of the VT threshold (B1). Near a
                  known VT, not a new lineage (Septoglomus sp. C4 at 96.93%).
  4. placement / novelty, from the Phase 3 call:
       known species / known species group   all matched arms agree and neither VT nor LSU is unmatched (K1)
       novel species, genus known            a species call demoted because VT or LSU is unmatched (N1, N2),
                                             or no species survives while a genus does (N3)
       discordant                            matched arms share no genus (X1)
       novel lineage                         no arm matched (N4; N5 when the best VT is also below the
                                             boundary, i.e. off every VT)

Scoring against what the culture reference knows (the reconcile run's holdout decides the truth groups):
  main-type copies, species held by another culture    expected known species or group (the M3 control)
  main-type copies, species held out                   expected novel species / novel lineage, not a
                                                       known species (false species)
  divergent-class copies (reference labels)            expected divergent (short 5.8S) or, for same-length
                                                       copies, anything but a false call of novelty
Writes {out}.copies.tsv with status, reason codes and placement for every copy.

Usage (from amfsplit_stage3):
  python3 novelty.py --ref ../ref/culture_ref/v4 --reconcile ../bench/stefani/reconcile_v4 \\
      --vt ../bench/stefani/vt_sets_v4 --sh ../bench/stefani/sets_SH_v4 --lsu ../bench/stefani/sets_LSU_v4 \\
      --out ../bench/stefani/novelty_v4
  (and the same with --reconcile ../bench/stefani/reconcile_v4_species --holdout species)
"""
import argparse
import collections
import os
import statistics
import sys

import edlib

from arms import is_named, read_fasta, read_tsv, write_tsv


def dist(a, b, k=-1):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    d = edlib.align(q, t, mode="HW", task="distance", k=k)["editDistance"]
    return None if d < 0 else 100 * d / len(q)


def split(s, sep=","):
    return [u for u in s.split(sep) if u] if s else []


def load_sets(prefix):
    rows = read_tsv(f"{prefix}.copies.tsv")
    col = "unit_set" if "unit_set" in rows[0] else "vt_set"
    top = "top_units" if "top_units" in rows[0] else "top_vts"
    units = {r["accession"]: set(split(r[col])) for r in rows}
    best = {}
    for r in rows:
        t = split(r[top])
        best[r["accession"]] = (t[0].rsplit(":", 1)[0], float(t[0].rsplit(":", 1)[1])) if t else ("", 0.0)
    return units, best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--reconcile", required=True, help="prefix of a reconcile.py run")
    ap.add_argument("--vt", required=True)
    ap.add_argument("--sh")
    ap.add_argument("--lsu")
    ap.add_argument("--holdout", choices=["culture", "species"], default="culture",
                    help="must match the reconcile run")
    ap.add_argument("--artifact-min-ssu", type=float, default=1.0)
    ap.add_argument("--artifact-ratio", type=float, default=0.5)
    ap.add_argument("--short-58s", type=int, default=3)
    ap.add_argument("--subst-58s", type=float, default=3.0)
    ap.add_argument("--vt-min-id", type=float, default=97.0)
    ap.add_argument("--boundary", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ref_rows = [r for r in read_tsv(os.path.join(a.ref, "copies.tsv")) if r["status"] == "retained"]
    cop = {r["accession"]: r for r in ref_rows}
    cul_meta = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    rec = {r["accession"]: r for r in read_tsv(f"{a.reconcile}.copies.tsv")}
    groups = collections.Counter(r["group"] for r in rec.values())
    if (a.holdout == "species") != (groups.get("named, species elsewhere", 0) == 0):
        raise SystemExit(f"--holdout {a.holdout} does not match the reconcile run ({dict(groups)})")
    fa = {m: read_fasta(os.path.join(a.ref, "fasta", f"{m}.fasta")) for m in ("ssu_flank", "full_its", "s58")}
    vt_units, vt_best = load_sets(a.vt)
    arm_units = {"VT": vt_units}
    for name, p in (("SH", a.sh), ("LSU", a.lsu)):
        if p:
            arm_units[name] = load_sets(p)[0]

    species = {c: m["organism"] for c, m in cul_meta.items()
               if is_named(m["organism"]) and m.get("name_status", "ok") == "ok"}
    by_species = collections.defaultdict(set)
    for c, sp in species.items():
        by_species[sp].add(c)

    def held(c):
        return by_species[species[c]] if a.holdout == "species" and c in species else {c}

    mains = [x for x, r in cop.items() if r["copy_class"] == "main" and x in fa["ssu_flank"]]
    main_58 = collections.defaultdict(list)
    for x in mains:
        if cop[x].get("s58_len"):
            main_58[cop[x]["culture"]].append(x)

    # units reached by main-type vs divergent copies, per culture, for the paralog-unit test
    reach = {n: (collections.defaultdict(set), collections.defaultdict(set)) for n in arm_units}
    for n, units in arm_units.items():
        for x, us in units.items():
            if x not in cop:
                continue
            side = 0 if cop[x]["copy_class"] == "main" else 1
            for u in us:
                reach[n][side][u].add(cop[x]["culture"])

    out = []
    for x, r in cop.items():
        c = r["culture"]
        drop = held(c)
        reasons, flags = [], []
        # nearest main-type copy of another culture, by SSU flank (bounded edlib, same-VT candidates first)
        near, near_d = None, None
        if x in fa["ssu_flank"]:
            q = fa["ssu_flank"][x]
            mine = vt_units.get(x, set())
            cands = [y for y in mains if cop[y]["culture"] not in drop]
            cands.sort(key=lambda y: not (vt_units.get(y, set()) & mine))
            for y in cands:
                k = -1 if near_d is None else int(near_d * min(len(q), len(fa["ssu_flank"][y])) / 100) + 1
                d = dist(q, fa["ssu_flank"][y], k)
                if d is not None and (near_d is None or d < near_d):
                    near, near_d = y, d
        its_d = dist(fa["full_its"][x], fa["full_its"][near]) if near and x in fa["full_its"] and near in fa["full_its"] else None

        # 1. artifact
        artifact = (near_d is not None and its_d is not None and near_d >= a.artifact_min_ssu
                    and near_d >= a.artifact_ratio * its_d)
        if artifact:
            reasons.append(f"A1 SSU {near_d:.1f}% vs ITS {its_d:.1f}% to nearest culture")

        # 2. divergent class from the copy's own sequence
        short = paralog = subst = False
        ncul = cop[near]["culture"] if near else ""
        if r.get("s58_len") and main_58.get(ncul):
            ref_len = statistics.median(int(cop[y]["s58_len"]) for y in main_58[ncul])
            delta = int(r["s58_len"]) - ref_len
            if delta <= -a.short_58s:
                short = True
                reasons.append(f"D1 5.8S {r['s58_len']} bp, {-delta:g} bp shorter than {ncul.split(' | ')[-1]}")
            elif abs(delta) < a.short_58s and x in fa["s58"]:
                d58 = min(dist(fa["s58"][x], fa["s58"][y]) for y in main_58[ncul] if y in fa["s58"])
                if d58 >= a.subst_58s:
                    subst = True
                    flags.append(f"D3 5.8S {d58:.1f}% from {ncul.split(' | ')[-1]}, same length")
        for n, units in arm_units.items():
            for u in units.get(x, ()):
                m_c = reach[n][0].get(u, set()) - drop
                d_c = reach[n][1].get(u, set()) - drop
                if d_c and not m_c:
                    paralog = True
                    reasons.append(f"D2 paralog unit {u} ({n})")
        divergent = short or paralog

        # 3. boundary
        rc = rec[x]
        vt_state = rc["VT_state"]
        best_vt, best_id = vt_best.get(x, ("", 0.0))
        boundary = vt_state == "no call" and best_id >= a.vt_min_id - a.boundary
        off_vt = vt_state == "no call" and best_id < a.vt_min_id - a.boundary

        # 4. placement and novelty
        rank, um = rc["rank"], set(split(rc["unmatched_arms"], ";"))
        if artifact:
            status = "artifact"
        elif divergent:
            status = "divergent class"
            place = rc["VT_species"] or rc["VT_genera"]
            reasons.append(f"placed by VT: {place or 'no VT match'}")
        elif boundary:
            status = "boundary"
            reasons.append(f"B1 best VT {best_vt} {best_id:.2f}%")
        elif rank in ("species", "species group"):
            if um & {"VT", "LSU"}:
                status = "novel species, genus known"
                reasons += [f"{'N1' if n == 'VT' else 'N2'} {n} unmatched" for n in ("VT", "LSU") if n in um]
            else:
                status = "known " + rank
                reasons.append("K1")
        elif rank == "genus":
            status = "novel species, genus known"
            reasons.append("N3 no species shared by the arms")
        elif rank == "discordant":
            status = "discordant"
            reasons.append("X1 arms share no genus")
        else:
            status = "novel lineage"
            reasons.append("N5 off every VT" if off_vt else "N4 no arm matched")
        placement = (rc["call_species"] if status.startswith("known") else
                     (rc["VT_species"] or rc["VT_genera"]) if status == "divergent class" else
                     rc["call_genera"] if status == "novel species, genus known" else "")
        out.append({"accession": x, "culture": c, "organism": r["organism"], "copy_class": r["copy_class"],
                    "div_kind": r.get("div_kind", ""), "group": rc["group"], "status": status,
                    "placement": placement, "reasons": "; ".join(reasons), "flags": "; ".join(flags),
                    "nearest_culture": ncul, "nearest_ssu_pct": round(near_d, 2) if near_d is not None else "",
                    "nearest_its_pct": round(its_d, 2) if its_d is not None else "",
                    "best_vt": best_vt, "best_vt_pident": best_id, "reconcile_rank": rank,
                    "unmatched_arms": rc["unmatched_arms"]})

    # ---- report
    def truth(r):
        if r["copy_class"] == "divergent":
            return f"divergent ({r['div_kind'] or '?'})"
        if r["copy_class"] != "main":
            return r["copy_class"]
        return {"named, species elsewhere": "main, species known", "named, species held out": "main, species held out",
                "unnamed ('sp.')": "main, 'sp.' culture"}[r["group"]]
    statuses = ["known species", "known species group", "novel species, genus known", "novel lineage",
                "boundary", "divergent class", "discordant", "artifact"]
    heads = ["known sp.", "known grp", "novel sp.", "novel lin.", "boundary", "divergent", "discord.", "artifact"]
    print(f"{len(out)} copies; holdout {a.holdout}; arms {', '.join(arm_units)}")
    print(f"\n1. status by what the reference knows (rows) ")
    print(f"  {'':30s}{'copies':>7s}" + "".join(f"{h:>11s}" for h in heads))
    rows_order = ["main, species known", "main, species held out", "main, 'sp.' culture",
                  "divergent (short_5.8S)", "divergent (same_length)", "divergent (longer)", "divergent (?)", "outlier"]
    for t in rows_order:
        rs = [r for r in out if truth(r) == t]
        if rs:
            cnt = collections.Counter(r["status"] for r in rs)
            print(f"  {t:30s}{len(rs):>7d}" + "".join(f"{cnt[s]:>11d}" for s in statuses))

    known = [r for r in out if truth(r) == "main, species known"]
    heldo = [r for r in out if truth(r) == "main, species held out"]
    short = [r for r in out if truth(r) == "divergent (short_5.8S)"]
    print("\n2. checks")
    if known:
        bad = [r for r in known if not r["status"].startswith("known")]
        print(f"  M3 negative control: main-type copies of known species not called known: {len(bad)} of {len(known)}")
        for cu, n in collections.Counter(r["culture"] for r in bad).most_common(10):
            r0 = next(r for r in bad if r["culture"] == cu)
            print(f"      {cu} ({n}): {r0['status']} | {r0['reasons']}"[:220])
    if heldo:
        fs = [r for r in heldo if r["status"].startswith("known")]
        print(f"  false species: held-out species called a known species or group: {len(fs)} of {len(heldo)}"
              f" ({100 * len(fs) / len(heldo):.0f}%)")
        for cu, n in collections.Counter(r["culture"] for r in fs).most_common(8):
            r0 = next(r for r in fs if r["culture"] == cu)
            print(f"      {cu} ({n}): {r0['status']} {r0['placement']}"[:200])
    if short:
        det = sum(r["status"] == "divergent class" for r in short)
        how = collections.Counter(code for r in short if r["status"] == "divergent class"
                                  for code in ("D1", "D2") if code in r["reasons"])
        print(f"  short-5.8S copies recognised from their own sequence: {det} of {len(short)} "
              f"(5.8S length {how['D1']}, paralog unit {how['D2']}); called novel: "
              f"{sum(r['status'].startswith('novel') for r in short)}")
    fp = [r for r in out if r["copy_class"] == "main" and r["status"] == "divergent class"]
    print(f"  main-type copies called divergent: {len(fp)}")
    for r in fp[:6]:
        print(f"      {r['accession']} {r['culture']}: {r['reasons']}"[:200])
    art = [r for r in out if r["status"] == "artifact"]
    print(f"  artifacts: {len(art)}" + (" (the culture reference should have none left)" if not art else ""))
    for r in art[:6]:
        print(f"      {r['accession']} {r['culture']}: {r['reasons']}"[:200])
    d3 = [r for r in out if "D3" in r["flags"]]
    print(f"  same-length 5.8S flags (D3): {len(d3)} copies "
          f"({sum(r['copy_class'] == 'divergent' for r in d3)} of them divergent-class in the reference)")
    ratio = [r["nearest_ssu_pct"] / r["nearest_its_pct"] for r in out
             if r["copy_class"] == "main" and r["nearest_its_pct"] not in ("", 0, 0.0) and r["nearest_ssu_pct"] != ""]
    if ratio:
        q = statistics.quantiles(ratio, n=100)
        print(f"  SSU/ITS divergence ratio to the nearest culture, main-type copies: median {statistics.median(ratio):.3f}, "
              f"95th {q[94]:.3f}, max {max(ratio):.3f} (artifact needs >= {a.artifact_ratio} with SSU >= "
              f"{a.artifact_min_ssu}%)")

    write_tsv(f"{a.out}.copies.tsv", out, list(out[0]))
    print(f"\ntable: {a.out}.copies.tsv")


if __name__ == "__main__":
    sys.exit(main())
