#!/usr/bin/env python3
"""Stage 1 VT check: assign MaarjAM virtual taxa (VT) to rDNA copies from their SSU flank.

Each copy's SSU flank (the ITSxRust region) is searched with blastn against the MaarjAM
SSU set. The references span NS31-AML2 (~520 bp), so they align inside the flank, and
copies with no NS31-AM1 window (Paraglomerales, where AM1 fails) stay in the check.

Assignment follows the usual MaarjAM criteria: best hit by bit score with identity >= 97%,
an alignment covering >= 95% of the shorter sequence (here the reference), and
e-value <= 1e-50. A copy whose top-scoring hits span several VT is called ambiguous.
The margin column is best identity minus the best identity to any other VT, so a small
margin means the copy sits near a VT boundary.

Reports:
  * assignment per order (order of the best hit)
  * cultures whose assigned copies fall in more than one VT
  * per QC flag class (from stefani_qc.flags.tsv): does a flagged copy get the same VT
    as the unflagged copies of its culture? (prediction 1 uses divergent_5.8S)
  * --focus-copy (default PX214583.1) against its culture-mates (prediction 2)
Writes {out}.copies.tsv, one row per copy.

Usage:
  vt_check.py --query ../bench/stefani/out/windows/stefani.ssu_flank.fasta \\
      --maarjam ../ref/maarjam/maarjam.fasta --groups ../bench/stefani/stefani_groups.tsv \\
      --flags ../bench/stefani/stefani_qc.flags.tsv --out ../bench/stefani/vt_check
"""
import argparse
import collections
import csv
import os
import re
import shutil
import statistics
import subprocess
import sys


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def load_maarjam(path, taxa_path):
    """Return {ref_id: info} and write a copy of the set with simple unique IDs.

    MaarjAM headers: >gi|MAARJAM_ID|gb|ACCESSION| Family Genus species VTX00000
    MaarjAM IDs are not unique in the gDAT snapshot, so references get ref00001-style IDs.
    """
    order_of = {}
    if taxa_path and os.path.exists(taxa_path):
        for line in open(taxa_path):
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                lin = [x.strip() for x in parts[1].split(";")]
                if len(lin) >= 5:
                    order_of[lin[4]] = lin[3]  # family -> order
    refs, out, seq = {}, [], []
    header = None

    def flush():
        if header is None:
            return
        idpart, _, rest = header.partition(" ")
        p, toks = idpart.split("|"), rest.split()
        if len(p) < 4 or not toks or not re.fullmatch(r"VTX\d+", toks[-1]):
            raise SystemExit(f"unexpected MaarjAM header: {header}")
        rid = f"ref{len(refs) + 1:05d}"
        fam = toks[0]
        refs[rid] = {"maarjam_id": p[1], "accession": p[3], "family": fam,
                     "order": order_of.get(fam, fam), "taxon": " ".join(toks[1:-1]),
                     "vt": toks[-1]}
        out.append(f">{rid}\n{''.join(seq)}\n")

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                flush()
                header, seq = line[1:], []
            elif line:
                seq.append(line)
    flush()
    return refs, "".join(out)


def run_blast(query, ref_fasta, workdir, threads, evalue):
    os.makedirs(workdir, exist_ok=True)
    db = os.path.join(workdir, "maarjam")
    with open(db + ".fasta", "w") as fh:
        fh.write(ref_fasta)
    subprocess.run(["makeblastdb", "-in", db + ".fasta", "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL)
    cols = "qseqid sseqid pident length qlen slen evalue bitscore"
    res = subprocess.run(
        ["blastn", "-task", "megablast", "-query", query, "-db", db, "-dust", "no",
         "-outfmt", f"6 {cols}", "-evalue", str(evalue), "-max_target_seqs", "5000",
         "-max_hsps", "1", "-num_threads", str(threads)],
        check=True, capture_output=True, text=True)
    hits = collections.defaultdict(list)
    for line in res.stdout.splitlines():
        q, s, pid, ln, ql, sl, ev, bs = line.split("\t")
        hits[q].append({"ref": s, "pident": float(pid), "length": int(ln), "qlen": int(ql),
                        "slen": int(sl), "evalue": float(ev), "bitscore": float(bs)})
    return hits


def assign(hs, refs, min_id, min_cov, max_ev):
    """Best-hit VT call for one copy."""
    row = {"call": "no_hit", "vt": "", "pident": "", "cov": "", "best_vt": "", "best_pident": "",
           "best_taxon": "", "order": "", "margin": "", "tied_vts": ""}
    if not hs:
        return row
    for h in hs:
        h["cov"] = h["length"] / min(h["qlen"], h["slen"])
        h["vt"] = refs[h["ref"]]["vt"]
    best = max(hs, key=lambda h: (h["bitscore"], h["pident"]))
    info = refs[best["ref"]]
    best_by_vt = {}
    for h in hs:
        if h["pident"] > best_by_vt.get(h["vt"], -1):
            best_by_vt[h["vt"]] = h["pident"]
    others = [p for v, p in best_by_vt.items() if v != best["vt"]]
    row.update(best_vt=best["vt"], best_pident=round(best["pident"], 2), best_taxon=info["taxon"],
               order=info["order"], cov=round(best["cov"], 3),
               margin=round(best["pident"] - max(others), 2) if others else "")
    passing = [h for h in hs if h["pident"] >= min_id and h["cov"] >= min_cov and h["evalue"] <= max_ev]
    if not passing:
        row["call"] = "unassigned"
        return row
    top = max(h["bitscore"] for h in passing)
    tied = sorted({h["vt"] for h in passing if h["bitscore"] == top})
    pbest = max((h for h in passing if h["bitscore"] == top), key=lambda h: h["pident"])
    row.update(pident=round(pbest["pident"], 2), vt=tied[0] if len(tied) == 1 else "",
               call="assigned" if len(tied) == 1 else "ambiguous",
               tied_vts=",".join(tied) if len(tied) > 1 else "")
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", required=True, help="SSU flank FASTA, IDs = accessions")
    ap.add_argument("--maarjam", required=True, help="MaarjAM SSU FASTA (gDAT snapshot format)")
    ap.add_argument("--taxa", help="maarjam.taxa.VT.txt (default: next to --maarjam)")
    ap.add_argument("--groups", required=True, help="gb_groups.py TSV")
    ap.add_argument("--level", choices=["culture", "library"], default="culture")
    ap.add_argument("--flags", help="stefani_qc.py flags TSV (accession, check)")
    ap.add_argument("--focus-copy", default="PX214583.1")
    ap.add_argument("--min-id", type=float, default=97.0)
    ap.add_argument("--min-cov", type=float, default=0.95)
    ap.add_argument("--max-evalue", type=float, default=1e-50)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()

    for tool in ("blastn", "makeblastdb"):
        if not shutil.which(tool):
            raise SystemExit(f"{tool} not found; install BLAST+: "
                             "mamba install -n amfsplit -c bioconda -c conda-forge blast")

    taxa = args.taxa or os.path.join(os.path.dirname(os.path.abspath(args.maarjam)), "maarjam.taxa.VT.txt")
    refs, ref_fasta = load_maarjam(args.maarjam, taxa)
    queries = read_fasta(args.query)
    with open(args.groups) as fh:
        grp = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    flags = collections.defaultdict(set)
    if args.flags:
        with open(args.flags) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                flags[r["accession"]].add(r["check"])

    missing = [a for a in queries if a not in grp]
    if missing:
        raise SystemExit(f"{len(missing)} query IDs are not in {args.groups}, e.g. {missing[:3]}")
    n_vt = len({r["vt"] for r in refs.values()})
    print(f"{len(queries)} copies vs {len(refs)} MaarjAM references ({n_vt} VT)")

    hits = run_blast(args.query, ref_fasta, args.out + ".blastdb", args.threads, 1e-10)
    rows = []
    for acc in queries:
        g = grp[acc]
        group = g[args.level]
        organism = g.get("organism") or group.split(" | ")[0]
        r = {"accession": acc, "organism": organism, "group": group, "ssu_len": len(queries[acc]),
             "flags": ",".join(sorted(flags.get(acc, ())))}
        r.update(assign(hits.get(acc, []), refs, args.min_id, args.min_cov, args.max_evalue))
        rows.append(r)

    cols = ["accession", "organism", "group", "ssu_len", "flags", "call", "vt", "pident", "cov",
            "tied_vts", "best_vt", "best_pident", "best_taxon", "order", "margin"]
    with open(args.out + ".copies.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # 1. assignment per order
    print(f"\nassignment (>= {args.min_id:g}% identity, >= {args.min_cov:.0%} coverage), by order of best hit")
    print(f"{'order':<17}{'copies':>7}{'assigned':>10}{'ambiguous':>11}{'unassigned':>12}"
          f"{'no hit':>8}{'median id':>11}")
    by_order = collections.defaultdict(list)
    for r in rows:
        by_order[r["order"] or "(no hit)"].append(r)
    for o, rs in sorted(by_order.items()):
        c = collections.Counter(r["call"] for r in rs)
        ids = [r["best_pident"] for r in rs if r["best_pident"] != ""]
        med = f"{statistics.median(ids):.1f}" if ids else "-"
        print(f"{o:<17}{len(rs):>7}{c['assigned']:>10}{c['ambiguous']:>11}{c['unassigned']:>12}"
              f"{c['no_hit']:>8}{med:>11}")

    # 2. groups split across VT
    by_group = collections.defaultdict(list)
    for r in rows:
        by_group[r["group"]].append(r)
    multi = []
    for gname, rs in by_group.items():
        vts = collections.Counter(r["vt"] for r in rs if r["call"] == "assigned")
        if len(vts) > 1:
            multi.append((gname, vts, rs))
    n_eval = sum(sum(r["call"] == "assigned" for r in rs) >= 2 for rs in by_group.values())
    print(f"\n{args.level} groups with >= 2 assigned copies: {n_eval}; split across VT: {len(multi)}")
    for gname, vts, rs in sorted(multi, key=lambda x: x[0]):
        minority = [r for r in rs if r["call"] == "assigned" and r["vt"] != vts.most_common(1)[0][0]]
        det = "; ".join(f"{r['accession']} {r['vt']} margin {r['margin']}"
                        + (f" [{r['flags']}]" if r["flags"] else "") for r in minority)
        print(f"  {gname[:48]:<48} {dict(vts.most_common())}  minority: {det}")

    # 3. flagged copies vs unflagged culture-mates
    def main_vt(r):
        mates = collections.Counter(m["vt"] for m in by_group[r["group"]]
                                    if m is not r and not m["flags"] and m["call"] == "assigned")
        return mates.most_common(1)[0][0] if mates else None

    classes = sorted({c for s in flags.values() for c in s})
    if classes:
        print("\nflagged copies vs unflagged copies of the same group")
        print(f"{'flag class':<22}{'copies':>7}{'same VT':>9}{'other VT':>10}{'not comparable':>16}")
        for cls in classes:
            fr = [r for r in rows if cls in r["flags"].split(",")]
            same = diff = nc = 0
            odd = []
            for r in fr:
                mv = main_vt(r)
                if r["call"] != "assigned" or mv is None:
                    nc += 1
                elif r["vt"] == mv:
                    same += 1
                else:
                    diff += 1
                    odd.append(f"{r['accession']} {r['vt']} vs {mv}")
            print(f"{cls:<22}{len(fr):>7}{same:>9}{diff:>10}{nc:>16}")
            for o in odd[:8]:
                print(f"    {o}")

    # 4. focus copy
    fc = next((r for r in rows if r["accession"] == args.focus_copy), None)
    if fc:
        print(f"\n{args.focus_copy} ({fc['group']}): {fc['call']} {fc['vt'] or fc['best_vt']} "
              f"at {fc['best_pident']}% to {fc['best_taxon']}, margin {fc['margin']}"
              + (f", flags {fc['flags']}" if fc["flags"] else ""))
        top = {}
        for h in sorted(hits.get(args.focus_copy, []), key=lambda h: -h["bitscore"]):
            v = refs[h["ref"]]["vt"]
            if v not in top:
                top[v] = f"{v} {h['pident']:.2f}% ({refs[h['ref']]['taxon']})"
            if len(top) == 3:
                break
        print("  closest VT: " + "; ".join(top.values()))
        mates = [m for m in by_group[fc["group"]] if m is not fc]
        mv = collections.Counter(m["vt"] or m["call"] for m in mates)
        ids = [m["best_pident"] for m in mates if m["best_pident"] != ""]
        print(f"  {len(mates)} group-mates: {dict(mv.most_common())}, best-hit identity "
              f"{min(ids):.2f} to {max(ids):.2f}%" if ids else f"  {len(mates)} group-mates")
    elif args.focus_copy:
        print(f"\n{args.focus_copy} not in the query set")

    # 5. organism <-> VT
    org_vts = collections.defaultdict(set)
    vt_orgs = collections.defaultdict(set)
    for r in rows:
        if r["call"] == "assigned" and not re.search(r"\bsp\.?$", r["organism"]):
            org_vts[r["organism"]].add(r["vt"])
            vt_orgs[r["vt"]].add(r["organism"])
    many = sorted((o for o, v in org_vts.items() if len(v) > 1), key=lambda o: -len(org_vts[o]))
    shared = sorted((v for v, o in vt_orgs.items() if len(o) > 1), key=lambda v: -len(vt_orgs[v]))
    print(f"\nnamed species in more than one VT: {len(many)} of {len(org_vts)}"
          + (f", e.g. {', '.join(f'{o} ({len(org_vts[o])})' for o in many[:4])}" if many else ""))
    print(f"VT holding more than one named species: {len(shared)} of {len(vt_orgs)}"
          + (f", e.g. {', '.join(f'{v} ({len(vt_orgs[v])})' for v in shared[:4])}" if shared else ""))
    print(f"\nper-copy table: {args.out}.copies.tsv")


if __name__ == "__main__":
    sys.exit(main())
