#!/usr/bin/env python3
"""amfsplit v0.2 (Stage 1 prototype).

Turns each long AMF rDNA read (SSU-ITS-LSU amplicon, PacBio HiFi or ONT) into a
multi-marker record:

  * ITS regions (ITS1, 5.8S, ITS2, full ITS) and the SSU / LSU flanks, taken
    from ITSxRust's per-read anchor table (--anchors-tsv), and
  * "legacy barcode windows": the exact short-read regions that the MaarjAM,
    SSU-V4, LSU-D2 and Delavaux LSU literature was built on, cut from the long
    read at the same primer sites (IUPAC-aware, mismatch-tolerant search).

Cutting windows at the original primer sites keeps every window directly
comparable with its legacy reference database, which is what the Stage 3
VT <-> SH <-> LSU-clade bridge needs.

Coordinates in the output are 1-based, inclusive, in the rRNA-sense
orientation of the read (the same convention ITSxRust uses for minus-strand
reads).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass

import edlib

# ----------------------------------------------------------------------------- sequence helpers
IUPAC = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
EQUALITIES = [(code, base) for code, bases in IUPAC.items() for base in bases]
_COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_sequences(path: str):
    """Yield (read_id, sequence) from FASTA or FASTQ, gzipped or not."""
    with open_text(path) as fh:
        first = fh.readline()
        if not first:
            return
        if first.startswith("@"):
            header = first
            while header:
                seq = fh.readline().strip()
                fh.readline()
                fh.readline()
                yield header[1:].split()[0], seq.upper()
                header = fh.readline()
        elif first.startswith(">"):
            rid, chunks = first[1:].split()[0], []
            for line in fh:
                if line.startswith(">"):
                    yield rid, "".join(chunks).upper()
                    rid, chunks = line[1:].split()[0], []
                else:
                    chunks.append(line.strip())
            yield rid, "".join(chunks).upper()
        else:
            sys.exit(f"error: {path} is neither FASTA nor FASTQ")


# ----------------------------------------------------------------------------- panel
@dataclass
class Window:
    name: str
    target: str
    fwd_name: str
    fwd: str
    rev_name: str
    rev: str
    min_len: int
    max_len: int

    @property
    def rev_site(self) -> str:
        """Reverse primer as it appears in the rRNA-sense read."""
        return revcomp(self.rev)


def load_panel(path: str) -> list[Window]:
    windows = []
    with open(path) as fh:
        rows = (line for line in fh if line.strip() and not line.startswith("#"))
        for row in csv.DictReader(rows, delimiter="\t"):
            windows.append(Window(row["window"], row["target"], row["fwd_name"],
                                  row["fwd_seq"].upper(), row["rev_name"],
                                  row["rev_seq"].upper(), int(row["min_len"]),
                                  int(row["max_len"])))
    return windows


def find_site(site: str, seq: str, max_err: int):
    """Best infix match of a primer site in seq -> (start0, end0, edits, n_best) or None."""
    r = edlib.align(site, seq, mode="HW", task="locations", k=max_err,
                    additionalEqualities=EQUALITIES)
    if r["editDistance"] < 0:
        return None
    start, end = r["locations"][0]
    return start, end, r["editDistance"], len(r["locations"])


def max_edits(primer: str, frac: float) -> int:
    return max(1, int(len(primer) * frac))


# ----------------------------------------------------------------------------- ITSxRust anchors
def load_anchors(path: str | None) -> dict[str, dict]:
    if not path:
        return {}
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["read_id"]] = row
    return out


def _int(v: str | None):
    return int(v) if v not in (None, "") else None


def itsx_regions(row: dict, read_len: int) -> dict[str, tuple[int, int] | None]:
    """ITS sub-regions and flanks (1-based inclusive, oriented) from an anchors row."""
    reg: dict[str, tuple[int, int] | None] = {}
    i1 = (_int(row.get("its1_start")), _int(row.get("its1_end")))
    i2 = (_int(row.get("its2_start")), _int(row.get("its2_end")))
    fl = (_int(row.get("full_start")), _int(row.get("full_end")))
    reg["its1"] = i1 if None not in i1 else None
    reg["its2"] = i2 if None not in i2 else None
    reg["full_its"] = fl if None not in fl else None
    reg["s58"] = (i1[1] + 1, i2[0] - 1) if reg["its1"] and reg["its2"] else None
    reg["ssu_flank"] = (1, fl[0] - 1) if reg["full_its"] and fl[0] > 1 else None
    reg["lsu_flank"] = (fl[1] + 1, read_len) if reg["full_its"] and fl[1] < read_len else None
    return reg


def orientation_from_anchors(row: dict) -> str | None:
    for key in ("ssu_strand", "s58s_strand", "s58e_strand", "lsu_strand"):
        if row.get(key) in ("+", "-"):
            return row[key]
    return None


# ----------------------------------------------------------------------------- core
ITSX_REGIONS = ["ssu_flank", "its1", "s58", "its2", "full_its", "lsu_flank"]
USABLE = ("ok", "open_5p", "open_3p")


def orient_by_primers(seq: str, windows: list[Window], frac: float) -> str:
    """Fallback orientation when ITSxRust gave no anchors: majority of primer sites."""
    rc = revcomp(seq)
    score = {"+": 0, "-": 0}
    for w in windows:
        for site in (w.fwd, w.rev_site):
            k = max_edits(site, frac)
            if find_site(site, seq, k):
                score["+"] += 1
            if find_site(site, rc, k):
                score["-"] += 1
    if score["+"] == score["-"] == 0:
        return "?"
    return "+" if score["+"] >= score["-"] else "-"


def cut_window(w: Window, seq: str, frac: float, keep_primers: bool, open_ends: bool = False):
    f = find_site(w.fwd, seq, max_edits(w.fwd, frac))
    r = find_site(w.rev_site, seq, max_edits(w.rev_site, frac))
    rec = {"fwd_ed": f[2] if f else "", "rev_ed": r[2] if r else "",
           "multi": int(bool(f and f[3] > 1) or bool(r and r[3] > 1))}
    if not f and not r:
        return {**rec, "status": "both_missing"}
    if (not f or not r) and not open_ends:
        return {**rec, "status": "fwd_missing" if not f else "rev_missing"}
    open_tag = ""
    if not f:
        # forward site lies outside the amplicon: window runs from the read's 5' end
        s0, e0, open_tag = 0, (r[1] if keep_primers else r[0] - 1), "open_5p"
    elif not r:
        # reverse site lies outside the amplicon: window runs to the read's 3' end
        s0, e0, open_tag = (f[0] if keep_primers else f[1] + 1), len(seq) - 1, "open_3p"
    elif keep_primers:
        s0, e0 = f[0], r[1]
    else:
        s0, e0 = f[1] + 1, r[0] - 1
    if e0 < s0:
        return {**rec, "status": "out_of_order"}
    length = e0 - s0 + 1
    min_len = w.min_len // 2 if open_tag else w.min_len
    status = open_tag or "ok"
    if length < min_len:
        status = "too_short"
    elif length > w.max_len:
        status = "too_long"
    return {**rec, "status": status, "start": s0 + 1, "end": e0 + 1, "len": length}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--reads", required=True, help="reads (FASTA/FASTQ, .gz ok)")
    ap.add_argument("-a", "--anchors", help="ITSxRust --anchors-tsv for the same reads")
    ap.add_argument("-p", "--panel", required=True, help="primer-window panel TSV")
    ap.add_argument("-o", "--out", required=True, help="output prefix")
    ap.add_argument("--max-err-frac", type=float, default=0.15,
                    help="max edits per primer as a fraction of its length (default 0.15)")
    ap.add_argument("--keep-primers", action="store_true",
                    help="include primer sites in windows (default: primer-trimmed)")
    ap.add_argument("--open-ends", action="store_true",
                    help="when one primer site lies outside the amplicon, cut the window to the read "
                         "end and flag it open_5p/open_3p (e.g. NS31 windows on nested amplicons "
                         "that start downstream of NS31)")
    ap.add_argument("--only-anchored", action="store_true",
                    help="skip reads ITSxRust could not anchor (e.g. off-target host rDNA)")
    args = ap.parse_args(argv)

    windows = load_panel(args.panel)
    anchors = load_anchors(args.anchors)

    region_names = ITSX_REGIONS + [w.name for w in windows]
    fasta = {name: open(f"{args.out}.{name}.fasta", "w") for name in region_names}
    rec_fh = open(f"{args.out}.records.tsv", "w", newline="")
    cols = ["read_id", "read_len", "orient", "orient_source", "itsx_confidence"]
    cols += [f"{r}_coords" for r in ITSX_REGIONS]
    for w in windows:
        cols += [f"{w.name}_{k}" for k in ("status", "start", "end", "len", "fwd_ed", "rev_ed", "multi")]
    cols += ["windows_ok", "flank_consistent"]
    writer = csv.DictWriter(rec_fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    writer.writeheader()

    n_reads = 0
    n_processed = 0
    status_counts = {w.name: Counter() for w in windows}
    lengths = defaultdict(list)
    orient_src = Counter()
    region_counts = Counter()
    inconsistent = Counter()

    for rid, raw in read_sequences(args.reads):
        n_reads += 1
        row = anchors.get(rid)
        if row is None and args.only_anchored:
            orient_src["skipped_unanchored"] += 1
            continue
        n_processed += 1
        orient = orientation_from_anchors(row) if row else None
        source = "itsxrust" if orient else "primers"
        if orient is None:
            orient = orient_by_primers(raw, windows, args.max_err_frac)
        orient_src[source if orient != "?" else "unresolved"] += 1
        seq = revcomp(raw) if orient == "-" else raw

        rec = {"read_id": rid, "read_len": len(seq), "orient": orient,
               "orient_source": source if orient != "?" else "unresolved",
               "itsx_confidence": (row or {}).get("confidence", "unanchored")}

        regions = itsx_regions(row, len(seq)) if row else {r: None for r in ITSX_REGIONS}
        for name in ITSX_REGIONS:
            coords = regions.get(name)
            rec[f"{name}_coords"] = f"{coords[0]}-{coords[1]}" if coords else ""
            if coords and coords[1] >= coords[0]:
                region_counts[name] += 1
                lengths[name].append(coords[1] - coords[0] + 1)
                fasta[name].write(f">{rid} region={name} coords={coords[0]}-{coords[1]} "
                                  f"orient={orient}\n{seq[coords[0] - 1:coords[1]]}\n")

        n_ok = 0
        flank_ok = True
        for w in windows:
            res = cut_window(w, seq, args.max_err_frac, args.keep_primers, args.open_ends) \
                if orient != "?" \
                else {"status": "unoriented"}
            status_counts[w.name][res["status"]] += 1
            for k in ("status", "start", "end", "len", "fwd_ed", "rev_ed", "multi"):
                rec[f"{w.name}_{k}"] = res.get(k, "")
            if res["status"] in USABLE:
                n_ok += 1
                lengths[w.name].append(res["len"])
                fasta[w.name].write(f">{rid} window={w.name} coords={res['start']}-{res['end']} "
                                    f"status={res['status']} orient={orient} "
                                    f"fwd_ed={res['fwd_ed']} rev_ed={res['rev_ed']}\n"
                                    f"{seq[res['start'] - 1:res['end']]}\n")
                # cross-check primer-derived windows against HMM anchors
                full = regions.get("full_its")
                if full:
                    if w.name.startswith("SSU") and res["end"] >= full[0]:
                        flank_ok = False
                        inconsistent[w.name] += 1
                    if w.name.startswith("LSU") and res["start"] <= full[1]:
                        flank_ok = False
                        inconsistent[w.name] += 1
        rec["windows_ok"] = n_ok
        rec["flank_consistent"] = "yes" if flank_ok else "no"
        writer.writerow(rec)

    for fh in fasta.values():
        fh.close()
    rec_fh.close()

    def dist(v):
        if not v:
            return None
        q = statistics.quantiles(v, n=4) if len(v) > 1 else [v[0]] * 3
        return {"n": len(v), "median": statistics.median(v), "q1": q[0], "q3": q[2],
                "min": min(v), "max": max(v)}

    summary = {
        "reads": n_reads,
        "reads_processed": n_processed,
        "orientation_source": dict(orient_src),
        "itsxrust_regions": {r: {"reads": region_counts[r], "length": dist(lengths[r])}
                             for r in ITSX_REGIONS},
        "windows": {w.name: {"target": w.target,
                             "primers": f"{w.fwd_name}/{w.rev_name}",
                             "ok": status_counts[w.name]["ok"],
                             "open": status_counts[w.name]["open_5p"] + status_counts[w.name]["open_3p"],
                             "usable": sum(status_counts[w.name][k] for k in USABLE),
                             "ok_pct": round(100 * status_counts[w.name]["ok"] / max(1, n_processed), 1),
                             "usable_pct": round(100 * sum(status_counts[w.name][k] for k in USABLE)
                                                 / max(1, n_processed), 1),
                             "status": dict(status_counts[w.name]),
                             "length": dist(lengths[w.name]),
                             "crosses_its_anchor": inconsistent[w.name]}
                    for w in windows},
        "params": {"max_err_frac": args.max_err_frac, "keep_primers": args.keep_primers,
                   "open_ends": args.open_ends,
                   "only_anchored": args.only_anchored},
    }
    with open(f"{args.out}.summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)

    # short console report
    print(f"{n_reads} reads ({n_processed} processed)  orientation: {dict(orient_src)}",
          file=sys.stderr)
    for r in ITSX_REGIONS:
        d = summary["itsxrust_regions"][r]
        med = d["length"]["median"] if d["length"] else "-"
        print(f"  {r:<24} {d['reads']:>6} reads  median len {med}", file=sys.stderr)
    for w in windows:
        d = summary["windows"][w.name]
        med = d["length"]["median"] if d["length"] else "-"
        print(f"  {w.name:<24} {d['usable']:>6} usable ({d['usable_pct']}%; {d['open']} open)  "
              f"median len {med}  "
              f"{dict(status_counts[w.name])}", file=sys.stderr)


if __name__ == "__main__":
    main()
