#!/usr/bin/env python3
"""Compare two amfsplit records.tsv files for the same reads (e.g. v0.2 output against v0.3).

Answers the question that matters before the new version is used: do windows that were usable
before keep exactly the same coordinates, and where do statuses change? Handles CRLF (v0.2) and LF.

Usage:
  python3 compare_records.py OLD.records.tsv NEW.records.tsv [--show 8]
"""
import argparse
import collections
import csv
import sys

USABLE = ("ok", "open_5p", "open_3p")


def load(path):
    with open(path, newline="") as fh:
        text = fh.read().replace("\r\n", "\n")
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--show", type=int, default=8)
    a = ap.parse_args()
    old = {r["read_id"]: r for r in load(a.old)}
    new = {r["read_id"]: r for r in load(a.new)}
    only_old, only_new = set(old) - set(new), set(new) - set(old)
    common = sorted(set(old) & set(new))
    cols = [c for c in next(iter(old.values())) if c.endswith("_status")]
    windows = [c[:-len("_status")] for c in cols]
    print(f"{len(old)} reads in old, {len(new)} in new, {len(common)} shared; only old {len(only_old)}, only new {len(only_new)}")
    moved = []
    print(f"\n{'window':28s}{'usable old':>11s}{'usable new':>11s}{'lost':>6s}{'gained':>8s}{'coord change':>14s}   transitions (old -> new, changed statuses only)")
    for w in windows:
        uo = un = lost = gained = coord = 0
        trans = collections.Counter()
        for rid in common:
            so, sn = old[rid][f"{w}_status"], new[rid].get(f"{w}_status", "")
            uo += so in USABLE
            un += sn in USABLE
            if so in USABLE and sn not in USABLE:
                lost += 1
            if so not in USABLE and sn in USABLE:
                gained += 1
            if so in USABLE and sn in USABLE and (old[rid][f"{w}_start"], old[rid][f"{w}_end"]) != (new[rid][f"{w}_start"], new[rid][f"{w}_end"]):
                coord += 1
                moved.append((rid, w, old[rid][f"{w}_start"], old[rid][f"{w}_end"], new[rid][f"{w}_start"], new[rid][f"{w}_end"]))
            if so != sn:
                trans[f"{so}->{sn}"] += 1
        print(f"{w:28s}{uo:>11d}{un:>11d}{lost:>6d}{gained:>8d}{coord:>14d}   " + ", ".join(f"{k} {v}" for k, v in trans.most_common(6)))
    conc = [r for r in common if new[r].get("concatemer") == "yes"]
    print(f"\nconcatemers flagged in new: {len(conc)} of {len(common)} ({100 * len(conc) / max(1, len(common)):.2f}%)")
    if "flank_consistent" in old[common[0]]:
        bo = sum(old[r]["flank_consistent"] == "no" for r in common)
        bn = sum(new[r]["flank_consistent"] == "no" for r in common)
        print(f"windows crossing an ITS anchor (flank_consistent = no): old {bo}, new {bn}")
    if moved:
        print(f"\n{len(moved)} usable windows changed coordinates; first {a.show}:")
        for row in moved[:a.show]:
            print("   ", *row)
    else:
        print("\nno usable window changed coordinates")
    return 1 if moved else 0


if __name__ == "__main__":
    sys.exit(main())
