#!/usr/bin/env python3
"""Stage 2 prototype: copy-aware abundance estimation over culture genomes.

A read is scored against every retained copy of the culture reference on the ITS and LSU windows.
Hard assignment ('nearest') sends each read to the genome of its single nearest copy. The copy-aware
model instead treats each genome (culture) as a MIXTURE of its own copies, including the divergent
classes, so a read from a short-5.8S paralog is explained by its own genome and not by an unrelated
genome whose main copy happens to lie closer:

    L(r | g) = sum over copies c of g of  w_gc * exp(-beta * d(r, c))         (w_gc = 1 / n_g)
    L(r | novel) = exp(-beta * d0)
    P(r) = sum_g pi_g L(r | g) + pi_0 L(r | novel);   EM over (pi_g, pi_0)

d(r, c) is the edit distance per base (infix) averaged over the two windows. pi_0 is a novelty
component: reads farther than d0 from every copy are absorbed by it instead of being forced onto
the nearest genome. Abundances are aggregated to species and to species complexes
(species_groups.py). Abundance here means share of rDNA reads; genome copy number is not estimable.

Two modes:
  simulate  draw mock communities from the reference itself and score methods against the truth,
            leaving each read's source culture out of the reference ('new isolate of a known
            species'; a species with one culture is then genuinely unknown to the model).
  reads     estimate abundance from Stage 1 window FASTAs (ITS and LSU) of a real library.

Usage (from amfsplit_stage3):
  python3 copy_em.py simulate --ref ../ref/culture_ref/v4 --groups ../bench/stefani/species_groups_v4 \\
      --out ../bench/stefani/stage2_sim --threads 8
  python3 copy_em.py reads --ref ../ref/culture_ref/v4 --groups ../bench/stefani/species_groups_v4 \\
      --its LIB.full_its.amf.fasta --lsu LIB.LSU_LROR-FLR2.amf.fasta --out lib_stage2 --threads 8
"""
import argparse
import collections
import math
import os
import random
import sys
from multiprocessing import Pool

import edlib

from arms import is_named, read_fasta, read_tsv, write_tsv

MARKERS = (("full_its", 0.5), ("LSU_LROR-FLR2", 0.5))
COMP = str.maketrans("ACGT", "TGCA")


def dist(query: str, target: str) -> float:
    if not query or not target:
        return 1.0
    q, t = (query, target) if len(query) <= len(target) else (target, query)
    return edlib.align(q, t, mode="HW", task="distance")["editDistance"] / len(q)


class Reference:
    def __init__(self, ref_dir, groups_prefix=None):
        rows = [r for r in read_tsv(os.path.join(ref_dir, "copies.tsv")) if r["status"] == "retained"]
        cul = {r["culture"]: r for r in read_tsv(os.path.join(ref_dir, "cultures.tsv"))}
        self.seq = {m: read_fasta(os.path.join(ref_dir, "fasta", f"{m}.fasta")) for m, _ in MARKERS}
        self.copies = [r for r in rows if all(r["accession"] in self.seq[m] for m, _ in MARKERS)]
        self.genome = {r["accession"]: r["culture"] for r in self.copies}
        self.cultures = sorted({r["culture"] for r in self.copies})
        self.species = {c: (cul[c]["organism"] if is_named(cul[c]["organism"]) and cul[c].get("name_status", "ok") == "ok"
                            else f"{cul[c]['organism']} [{c.split(' | ')[-1]}]") for c in self.cultures}
        self.complex_of = {}
        if groups_prefix and os.path.exists(f"{groups_prefix}.complexes.tsv"):
            for r in read_tsv(f"{groups_prefix}.complexes.tsv"):
                for s in r["species"].split(";"):
                    self.complex_of[s] = r["complex"] + ": " + " / ".join(x.split()[-1] for x in r["species"].split(";"))
        # unique sequences per marker, and which copies share them
        self.uniq = {}
        for m, _ in MARKERS:
            seqs = sorted({self.seq[m][r["accession"]] for r in self.copies})
            self.uniq[m] = seqs
        self.index = {m: {s: i for i, s in enumerate(self.uniq[m])} for m, _ in MARKERS}

    def group_of(self, culture):
        sp = self.species[culture]
        return self.complex_of.get(sp, sp)


_REF = None


def _init(ref_dir, groups_prefix):
    global _REF
    _REF = Reference(ref_dir, groups_prefix)


def _read_dists(item):
    """Distances of one read to every unique reference sequence, per marker: {marker: [d...]}"""
    rid, seqs = item
    out = {}
    for m, _ in MARKERS:
        s = seqs.get(m, "")
        out[m] = [dist(s, t) for t in _REF.uniq[m]] if s else None
    return rid, out


def copy_distance_matrix(ref, reads, threads):
    """reads: {rid: {marker: seq}} -> {rid: [d(r, copy) for copy in ref.copies]} (weighted marker mean)."""
    if threads > 1:
        with Pool(threads, initializer=_init, initargs=(ref.ref_dir, ref.groups_prefix)) as p:
            res = dict(p.imap_unordered(_read_dists, list(reads.items()), chunksize=8))
    else:
        _init(ref.ref_dir, ref.groups_prefix)
        res = dict(map(_read_dists, reads.items()))
    mats = {}
    for rid, per in res.items():
        vec = []
        for r in ref.copies:
            num = den = 0.0
            for m, w in MARKERS:
                if per[m] is not None:
                    num += w * per[m][ref.index[m][ref.seq[m][r["accession"]]]]
                    den += w
            vec.append(num / den if den else 1.0)
        mats[rid] = vec
    return mats


def likelihoods(ref, dvec, beta, d0, exclude=None, nearest=False):
    """(dict culture -> L(r|g), L_novel, nearest culture). exclude: culture left out of the reference."""
    per = collections.defaultdict(list)
    best, best_c = 9.0, None
    for r, d in zip(ref.copies, dvec):
        c = ref.genome[r["accession"]]
        if c == exclude:
            continue
        per[c].append(d)
        if d < best:
            best, best_c = d, c
    L = {}
    for c, ds in per.items():
        if nearest:
            L[c] = math.exp(-beta * min(ds))
        else:
            L[c] = sum(math.exp(-beta * d) for d in ds) / len(ds)
    return L, math.exp(-beta * d0), best_c, best


def auto_d0(dmins, quantile=0.25, spread=0.045):
    """Novelty distance that adapts to the library: the noise floor (a low quantile of the
    nearest-copy distances, which is the read error plus isolate divergence) plus the within-species
    spread of the reference (mean of the 99th-percentile ITS and LSU spreads, 6.0% and 3.0%)."""
    v = sorted(dmins)
    return v[int(quantile * (len(v) - 1))] + spread


def em(ref, liks, novel_l, iters=300, prior=0.05, use_novel=True):
    cult = ref.cultures
    pi = {c: 1.0 / (len(cult) + 1) for c in cult}
    pi0 = 1.0 / (len(cult) + 1)
    n = len(liks)
    for _ in range(iters):
        acc = collections.defaultdict(float)
        acc0 = 0.0
        for L in liks:
            tot = sum(pi[c] * v for c, v in L.items()) + (pi0 * novel_l if use_novel else 0.0)
            if tot <= 0:
                continue
            for c, v in L.items():
                acc[c] += pi[c] * v / tot
            if use_novel:
                acc0 += pi0 * novel_l / tot
        new = {c: (acc[c] + prior / len(cult)) / (n + prior) for c in cult}
        new0 = (acc0 + prior / (len(cult) + 1)) / (n + prior) if use_novel else 0.0
        z = sum(new.values()) + new0
        new = {c: v / z for c, v in new.items()}
        new0 /= z
        delta = max(abs(new[c] - pi[c]) for c in cult)
        pi, pi0 = new, new0
        if delta < 1e-7:
            break
    return pi, pi0


def aggregate(ref, pi, pi0, level):
    out = collections.Counter()
    for c, v in pi.items():
        out[ref.species[c] if level == "species" else ref.group_of(c)] += v
    out["NOVEL"] += pi0
    return out


def tvd(p, q):
    keys = set(p) | set(q)
    return 0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in keys)


# ---------------------------------------------------------------------------- simulation
def noisy(seq, rate, rng):
    out = []
    for ch in seq:
        if rng.random() < rate:
            r = rng.random()
            if r < 0.6:
                out.append(rng.choice([b for b in "ACGT" if b != ch]))
            elif r < 0.8:
                out.append(ch + rng.choice("ACGT"))
        else:
            out.append(ch)
    return "".join(out)


def simulate(args):
    rng = random.Random(args.seed)
    d0_list = [x if x == "auto" else float(x) for x in (args.d0_scan or args.d0).split(",")]
    ref = Reference(args.ref, args.groups)
    ref.ref_dir, ref.groups_prefix = args.ref, args.groups
    by_cult = collections.defaultdict(list)
    for r in ref.copies:
        by_cult[r["culture"]].append(r)
    n_cult_of_species = collections.Counter(ref.species[c] for c in ref.cultures)
    named = [c for c in ref.cultures if is_named(ref.species[c].split(" [")[0])]
    results = []
    for k in range(args.communities):
        # communities of cultures from distinct species, log-normal abundances
        pool = rng.sample(named, len(named))
        chosen, seen = [], set()
        for c in pool:
            if ref.species[c] not in seen:
                chosen.append(c); seen.add(ref.species[c])
            if len(chosen) == args.species:
                break
        w = [math.exp(rng.gauss(0, 1.0)) for _ in chosen]
        z = sum(w)
        truth_c = {c: v / z for c, v in zip(chosen, w)}
        reads, source = {}, {}
        for i in range(args.reads):
            c = rng.choices(chosen, weights=w)[0]
            cp = rng.choice(by_cult[c])
            rid = f"c{k}_r{i}"
            reads[rid] = {m: noisy(ref.seq[m][cp["accession"]], args.error, rng) for m, _ in MARKERS}
            source[rid] = c
        mats = copy_distance_matrix(ref, reads, args.threads)
        truth = {"species": collections.Counter(), "group": collections.Counter()}
        for c, v in truth_c.items():
            truth["species"][ref.species[c]] += v
            truth["group"][ref.group_of(c)] += v
        # species present in the reference only through the left-out culture are unknown to the model
        unknown = {c for c in chosen if n_cult_of_species[ref.species[c]] == 1}
        dmins = [min(d for r, d in zip(ref.copies, dvec) if ref.genome[r["accession"]] != source[rid]) for rid, dvec in mats.items()]
        for d0_spec in d0_list:
          d0 = auto_d0(dmins, spread=args.spread) if d0_spec == "auto" else d0_spec
          for method in ("nearest", "copy_aware", "copy_aware_novel"):
            if method != "copy_aware_novel" and d0_spec != d0_list[0]:
                continue                                   # d0 only matters for the novelty component
            liks, novel_l, hard = [], None, collections.Counter()
            for rid, dvec in mats.items():
                L, nl, bc, bd = likelihoods(ref, dvec, args.beta, d0, exclude=source[rid], nearest=(method == "nearest"))
                liks.append(L); novel_l = nl
                hard[bc] += 1
            if method == "nearest":
                pi = {c: hard[c] / len(mats) for c in ref.cultures}; pi0 = 0.0
            else:
                pi, pi0 = em(ref, liks, novel_l, use_novel=(method == "copy_aware_novel"))
            for level in ("species", "group"):
                est = aggregate(ref, pi, pi0, level)
                t = dict(truth[level])
                for c in unknown:                      # truth for unknown species is 'NOVEL'
                    key = ref.species[c] if level == "species" else ref.group_of(c)
                    t["NOVEL"] = t.get("NOVEL", 0.0) + t.pop(key, 0.0)
                results.append({"community": k, "method": method, "d0": d0_spec if method == "copy_aware_novel" else "", "d0_used": round(d0, 3) if method == "copy_aware_novel" else "", "level": level,
                                "tvd": round(tvd(t, est), 4),
                                "novel_true": round(sum(truth_c[c] for c in unknown), 3), "novel_est": round(est.get("NOVEL", 0.0), 3),
                                "n_species": len(chosen), "n_reads": len(mats)})
        print(f"community {k}: {len(chosen)} species, {len(unknown)} unknown to the model", file=sys.stderr, flush=True)
    write_tsv(f"{args.out}.tsv", results, list(results[0]))
    print("\nmean total-variation distance to the true composition (0 = perfect), "
          f"{args.communities} communities x {args.reads} reads, {int(100 * args.error)}% read error")
    mean = lambda xs: sum(xs) / len(xs)
    known = {r["community"] for r in results if float(r["novel_true"]) == 0.0}
    print(f"  {'method':18s}{'d0':>13s}{'species':>9s}{'group':>8s}{'all known':>11s}{'with unknown':>14s}{'novel est (true)':>20s}")
    combos = [("nearest", ""), ("copy_aware", "")] + [("copy_aware_novel", d) for d in d0_list]
    used = {d: [r["d0_used"] for r in results if r["method"] == "copy_aware_novel" and r["d0"] == d and r["level"] == "species"] for d in d0_list}
    for method, d0 in combos:
        sel = [r for r in results if r["method"] == method and r["d0"] == d0]
        sp = [r for r in sel if r["level"] == "species"]
        gr = [r for r in sel if r["level"] == "group"]
        kn = [r["tvd"] for r in sp if r["community"] in known]
        un = [r["tvd"] for r in sp if r["community"] not in known]
        label = str(d0) + (f" ({mean(used[d0]):.3f})" if d0 == "auto" else "")
        print(f"  {method:18s}{label:>13s}{mean([r['tvd'] for r in sp]):>9.3f}{mean([r['tvd'] for r in gr]):>8.3f}"
              f"{(mean(kn) if kn else float('nan')):>11.3f}{(mean(un) if un else float('nan')):>14.3f}"
              f"{mean([r['novel_est'] for r in sp]):>13.3f} ({mean([r['novel_true'] for r in sp]):.3f})")
    return 0


# ---------------------------------------------------------------------------- real reads
def reads_mode(args):
    ref = Reference(args.ref, args.groups)
    ref.ref_dir, ref.groups_prefix = args.ref, args.groups
    its, lsu = read_fasta(args.its), read_fasta(args.lsu)
    ids = sorted(set(its) & set(lsu))
    if args.max_reads and len(ids) > args.max_reads:
        ids = random.Random(1).sample(ids, args.max_reads)
    reads = {i: {"full_its": its[i], "LSU_LROR-FLR2": lsu[i]} for i in ids}
    print(f"{len(ids)} reads with both windows (ITS {len(its)}, LSU {len(lsu)})", file=sys.stderr)
    mats = copy_distance_matrix(ref, reads, args.threads)
    dmins = [min(mats[rid]) for rid in ids]
    d0 = auto_d0(dmins, spread=args.spread) if args.d0 == "auto" else float(args.d0)
    print(f"novelty distance d0 = {d0:.3f} ({'adaptive: noise floor + spread' if args.d0 == 'auto' else 'fixed'}); "
          f"median nearest-copy distance {sorted(dmins)[len(dmins) // 2]:.3f}", file=sys.stderr)
    args.d0 = d0
    liks, per_read, novel_l = [], [], None
    for rid in ids:
        L, novel_l, bc, bd = likelihoods(ref, mats[rid], args.beta, d0)
        liks.append(L)
        per_read.append({"read_id": rid, "nearest_genome": bc, "nearest_species": ref.species[bc], "nearest_group": ref.group_of(bc),
                         "nearest_distance": round(bd, 4), "beyond_d0": int(bd > args.d0)})
    pi, pi0 = em(ref, liks, novel_l)
    write_tsv(f"{args.out}.reads.tsv", per_read, list(per_read[0]))
    for level in ("species", "group"):
        est = aggregate(ref, pi, pi0, level)
        rows = [{"unit": k, "share": round(v, 5), "reads_equiv": round(v * len(ids), 1)} for k, v in est.most_common() if v > 1e-4]
        write_tsv(f"{args.out}.{level}.tsv", rows, ["unit", "share", "reads_equiv"])
        print(f"\nabundance by {level} ({len(ids)} reads):")
        for r in rows[:12]:
            print(f"  {r['share']:7.3f}  {r['unit']}")
    hard = collections.Counter(r["nearest_group"] for r in per_read)
    print(f"\nreads beyond d0 = {args.d0} of every copy: {sum(r['beyond_d0'] for r in per_read)} ({100 * sum(r['beyond_d0'] for r in per_read) / len(ids):.1f}%)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    for name in ("simulate", "reads"):
        p = sub.add_parser(name)
        p.add_argument("--ref", required=True)
        p.add_argument("--groups", help="prefix of species_groups.py output")
        p.add_argument("--out", required=True)
        p.add_argument("--threads", type=int, default=4)
        p.add_argument("--beta", type=float, default=150.0, help="likelihood sharpness per unit edit distance (default 150)")
        p.add_argument("--d0", default="auto", help="novelty distance per base: a number, or 'auto' (default): "
                                                    "the library's noise floor plus --spread")
        p.add_argument("--spread", type=float, default=0.045, help="within-species spread added to the noise floor (default 0.045)")
        p.add_argument("--d0-scan", help="simulate: comma-separated d0 values to compare, e.g. 0.07,0.09,0.12")
    s = sub.choices["simulate"]
    s.add_argument("--communities", type=int, default=6)
    s.add_argument("--species", type=int, default=8)
    s.add_argument("--reads", type=int, default=300)
    s.add_argument("--error", type=float, default=0.05)
    s.add_argument("--seed", type=int, default=3)
    r = sub.choices["reads"]
    r.add_argument("--its", required=True)
    r.add_argument("--lsu", required=True)
    r.add_argument("--max-reads", type=int, default=0)
    args = ap.parse_args()
    return simulate(args) if args.mode == "simulate" else reads_mode(args)


if __name__ == "__main__":
    sys.exit(main())
