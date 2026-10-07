#!/usr/bin/env python3
"""Species complexes from the culture reference: pairs of named species that the markers cannot separate.

For each marker (full ITS, LROR-FLR2, SSU flank) the within-species spread T is the --quantile
(default 99th) percentile of the distance from a main-type copy to the nearest main-type copy of
another culture of the same species (pooled over species with >= 2 cultures). The percentile is
statistics.quantiles' exclusive method, the definition the frozen ITS rule (its_rule.py) uses;
fractional percentiles (97.5) interpolate on the same definition. Two named species are
INSEPARABLE when, on every marker, the nearest pair of copies across the two species is no further
apart than T: a copy of one can lie inside the spread of the other. Connected pairs form a complex.
Nothing here uses any test result: the rule reads only the frozen reference, so it can be fixed
before a benchmark and applied to it.

Usage: python3 species_groups.py --ref ../ref/culture_ref/v4 --out ../bench/stefani/species_groups_v4
Writes {out}.pairs.tsv (every species pair with its nearest distances), {out}.complexes.tsv and
{out}.limits.tsv (T per marker).
"""
import argparse, csv, collections, itertools, os, statistics, sys
import edlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arms import is_named, read_fasta, read_tsv, write_tsv

MARKERS = [("full_its", "ITS"), ("LSU_LROR-FLR2", "LSU"), ("ssu_flank", "SSU")]


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode="HW", task="distance")["editDistance"] / len(q)


def percentile(vals, p):
    """statistics.quantiles(vals, method='exclusive') at percentile p; any real p, identical at integers."""
    d = sorted(vals)
    n = len(d)
    h = (n + 1) * p / 100.0
    j = min(max(int(h), 1), n - 1)
    t = d[j - 1] + (d[j] - d[j - 1]) * (h - j)
    if float(p).is_integer() and abs(t - statistics.quantiles(vals, n=100)[int(p) - 1]) > 1e-9:
        raise AssertionError("percentile() departs from statistics.quantiles")
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--quantile", type=float, default=99.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cul = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    sp_of = {c: m["organism"] for c, m in cul.items() if is_named(m["organism"]) and m.get("name_status", "ok") == "ok"}
    cop = [r for r in read_tsv(os.path.join(a.ref, "copies.tsv"))
           if r["status"] == "retained" and r["copy_class"] == "main" and r["culture"] in sp_of]
    by_sp = collections.defaultdict(list)
    for r in cop:
        by_sp[sp_of[r["culture"]]].append(r)
    species = sorted(by_sp)
    # one representative set per culture keeps the work small: all main copies, deduplicated by sequence per marker
    seq = {m: {k: v for k, v in read_fasta(os.path.join(a.ref, "fasta", f"{m}.fasta")).items()} for m, _ in MARKERS}
    uniq = {}
    for m, _ in MARKERS:
        for sp in species:
            uniq[(m, sp)] = sorted({(seq[m][r["accession"]], r["culture"]) for r in by_sp[sp] if r["accession"] in seq[m]})
    # within-species spread
    T = {}
    for m, lab in MARKERS:
        vals = []
        for sp in species:
            items = uniq[(m, sp)]
            if len({c for _, c in items}) < 2:
                continue
            for s1, c1 in items:
                vals.append(min(dist(s1, s2) for s2, c2 in items if c2 != c1))
        T[m] = percentile(vals, a.quantile)
    print("within-species spread T (" + str(a.quantile) + "th pct): " + ", ".join(f"{lab} {T[m]:.2f}%" for m, lab in MARKERS))
    rows = []
    for s1, s2 in itertools.combinations(species, 2):
        row = {"species_a": s1, "species_b": s2}
        ok = True
        for m, lab in MARKERS:
            A, B = uniq[(m, s1)], uniq[(m, s2)]
            d = min(dist(x, y) for x, _ in A for y, _ in B) if A and B else None
            row[f"{lab}_min"] = round(d, 2) if d is not None else ""
            ok &= d is not None and d <= T[m]
        row["inseparable"] = int(ok)
        rows.append(row)
    write_tsv(f"{a.out}.pairs.tsv", rows, list(rows[0]))
    write_tsv(f"{a.out}.limits.tsv", [{"marker": lab, "quantile": a.quantile, "T": round(T[m], 4)} for m, lab in MARKERS],
              ["marker", "quantile", "T"])
    # complexes = connected components of inseparable pairs
    parent = {s: s for s in species}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for r in rows:
        if r["inseparable"]:
            parent[find(r["species_a"])] = find(r["species_b"])
    comp = collections.defaultdict(list)
    for s in species:
        comp[find(s)].append(s)
    cx = [sorted(v) for v in comp.values() if len(v) > 1]
    out = [{"complex": f"C{i + 1}", "species": ";".join(v), "n_species": len(v),
            "n_cultures": len({r['culture'] for s in v for r in by_sp[s]})} for i, v in enumerate(sorted(cx, key=lambda v: -len(v)))]
    write_tsv(f"{a.out}.complexes.tsv", out, ["complex", "species", "n_species", "n_cultures"]) if out else None
    print(f"{len(species)} named species; {sum(r['inseparable'] for r in rows)} inseparable pairs; {len(cx)} complexes covering {sum(len(v) for v in cx)} species")
    for o in out:
        print(f"  {o['complex']}: {o['species'].replace(';', ', ')}  ({o['n_cultures']} cultures)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
