#!/usr/bin/env python3
"""Independent check of the concatemer screen from read lengths alone.

A concatemer of two amplicon copies is about twice as long as a single copy. This script does not
use primers: it compares the read-length distribution of flagged and unflagged reads in an amfsplit
v0.3 records.tsv, and lists unflagged reads that are long enough to be candidate misses.

  * flagged reads should sit near 2x (or more) the median length of unflagged reads;
  * unflagged reads at >= --factor x the median are the candidate misses. Their number, over the
    number flagged, is a rough upper bound on how many concatemers the screen missed (long reads
    can also be genuine long amplicons, so it is an upper bound, not a count).

Usage:
  python3 check_concat_lengths.py LIB1.records.tsv [LIB2.records.tsv ...] [--factor 1.7] [--show 5]
"""
import argparse
import csv
import statistics
import sys


def load(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh.read().replace("\r\n", "\n").splitlines(), delimiter="\t"))


def pct(v, p):
    return statistics.quantiles(v, n=100)[p - 1] if len(v) > 1 else (v[0] if v else float("nan"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("records", nargs="+")
    ap.add_argument("--factor", type=float, default=1.7, help="length factor for candidate misses (default 1.7)")
    ap.add_argument("--show", type=int, default=5)
    a = ap.parse_args()
    for path in a.records:
        rows = load(path)
        if "concatemer" not in rows[0]:
            print(f"{path}: no concatemer column (records from v0.2?)")
            continue
        single = [int(r["read_len"]) for r in rows if r["concatemer"] == "no"]
        flagged = [r for r in rows if r["concatemer"] == "yes"]
        med = statistics.median(single)
        fl = [int(r["read_len"]) for r in flagged]
        cand = [r for r in rows if r["concatemer"] == "no" and int(r["read_len"]) >= a.factor * med]
        print(f"\n{path}\n  {len(rows)} anchored reads; unflagged length: median {med:.0f}, 5th {pct(single, 5):.0f}, 95th {pct(single, 95):.0f}, max {max(single)}")
        if flagged:
            near2 = sum(1 for x in fl if x >= 1.6 * med)
            print(f"  flagged: {len(flagged)} ({100 * len(flagged) / len(rows):.2f}%); median length {statistics.median(fl):.0f} "
                  f"= {statistics.median(fl) / med:.2f}x the unflagged median; {near2} of {len(flagged)} are >= 1.6x")
            if "concat_kind" in rows[0]:
                import collections as _c
                kc = _c.Counter(r["concat_kind"] for r in flagged)
                print(f"  kinds: {kc['full']} full, {kc['partial']} partial")
            if "concat_evidence" in rows[0]:
                import collections
                pat = collections.Counter(r["concat_evidence"].split("@")[0][:70] for r in flagged)
                print("  evidence of flagged reads (top): " + "; ".join(f"{k} x{v}" for k, v in pat.most_common(4)))
        else:
            print("  flagged: 0")
        print(f"  unflagged reads >= {a.factor}x the median: {len(cand)} ({100 * len(cand) / len(rows):.2f}% of reads); "
              f"candidate misses per flagged read: {len(cand) / max(1, len(flagged)):.1f}")
        for r in sorted(cand, key=lambda r: -int(r["read_len"]))[:a.show]:
            print(f"     {r['read_id']}  length {r['read_len']}  ({int(r['read_len']) / med:.2f}x)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
