#!/usr/bin/env python3
"""Re-run crosswalk.py's report on a newer culture reference without re-running BLAST.

crosswalk.py's per-copy calls (VT, SH, LSU) depend only on the query window, the external
reference and the set of culture-reference accessions whose hits are dropped. Culture names,
copy classes and name holds enter only after the calls, in the culture-level report. When a
reference edit changes only labels and classes, the stored per-copy calls of the earlier run
are therefore the calls of the new reference, and the report can be replayed on them.

This script checks that premise and replays the report:
  1. the retained accessions are the same in both references, and every copy whose label
     or class changed is listed;
  2. with --evidence (the frozen per-unit evidence captured on the new reference, M2b
     inputs/benchmark/evidence), the query windows hash to the evidence, and on every arm and
     copy the stored call's identity equals the evidence anchor (identity of the highest-bit-
     score hit passing the threshold) and the stored best unit is in the evidence at an
     identity at least that of the stored best hit; any disagreement stops the replay;
  3. crosswalk.py's main() is run unchanged with BLAST replaced by the stored calls.
     Reference record names are rebuilt from the stored best-hit taxa (display only).

Usage:
  python3 crosswalk_replay.py --code <dir with crosswalk.py, arms.py> --ref <culture_ref/v4> \\
      --calls crosswalk_v4_copies.tsv [--evidence <M2b inputs/benchmark/evidence>] --out crosswalk_v4final
Writes the crosswalk tables with --out as prefix, {out}.log and {out}.replay.json.
"""
import argparse
import collections
import contextlib
import gzip
import hashlib
import io
import json
import os
import sys
import types

ARM_KIND = {"VT": "maarjam", "SH": "unite", "LSU": "delavaux"}
KIND_ARM = {v: k for k, v in ARM_KIND.items()}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--code", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--calls", required=True, help="{out}.copies.tsv of the earlier crosswalk run")
    ap.add_argument("--evidence", help="frozen per-unit evidence captured on --ref (optional check)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    sys.path.insert(0, os.path.abspath(a.code))
    import arms as armsmod
    import crosswalk as cw

    calls = {r["accession"]: r for r in armsmod.read_tsv(a.calls)}
    ref_rows = armsmod.read_tsv(os.path.join(a.ref, "copies.tsv"))
    retained = {r["accession"]: r for r in ref_rows if r["status"] == "retained"}
    if set(calls) != set(retained):
        raise SystemExit(f"accessions differ: {len(set(calls) ^ set(retained))} not in both")
    man = json.load(open(os.path.join(a.ref, "MANIFEST.json")))

    # what the edits changed between the two runs (labels and classes only)
    changed = []
    for x, c in calls.items():
        r = retained[x]
        diff = {k: (c[k], r[k]) for k in ("culture", "organism", "order", "copy_class", "name_status")
                if c[k] != r[k]}
        if diff:
            changed.append({"accession": x, **{k: f"{o} -> {n}" for k, (o, n) in diff.items()}})

    # premise 2: stored best hits against the frozen evidence of the new reference
    ev_check = {}
    if a.evidence:
        for arm in ("VT", "SH", "LSU"):
            e = json.load(gzip.open(os.path.join(a.evidence, f"{arm}.evidence.json.gz"), "rt"))
            if e["copies_sha256"] != sha256(os.path.join(a.ref, "copies.tsv")):
                raise SystemExit(f"{arm} evidence was not captured on this reference")
            qpath = os.path.join(a.ref, "fasta", f"{e['region']}.fasta")
            if e["query_sha256"] != sha256(qpath):
                raise SystemExit(f"{arm}: query windows differ from the evidence")
            # anchor = identity of the highest-bit-score hit passing the identity threshold (arm_sets.py),
            # i.e. the identity of the stored call; the stored best unit must be in the evidence at
            # an identity at least that of the stored best hit (rounded to two decimals there)
            anchor_ok = best_ok = 0
            bad_ex = []
            for x, c in calls.items():
                q = e["queries"][x]
                a_, p = q["anchor_pident"], c[f"{arm}_pident"]
                if (a_ is None and p == "") or (a_ is not None and p != "" and abs(a_ - float(p)) < 0.006):
                    anchor_ok += 1
                else:
                    bad_ex.append((x, "anchor", a_, p))
                b, bp = c[f"{arm}_best_label"], c[f"{arm}_best_pident"]
                if not b or q["unit_pident"].get(b, -1) >= float(bp) - 0.006:
                    best_ok += 1
                else:
                    bad_ex.append((x, "best", b, bp))
            ev_check[arm] = {"copies": len(calls), "call_identity_equals_anchor": anchor_ok,
                             "best_unit_in_evidence": best_ok, "examples_of_disagreement": bad_ex[:5]}
            if bad_ex:
                raise SystemExit(f"{arm}: stored calls disagree with the evidence of this reference: {bad_ex[:5]}")

    # stand-ins for the BLAST step: the stored calls, keyed through a marker hit
    state = {"arm": None}

    def load_reference(path, kind, keep=None):
        arm = KIND_ARM[kind]
        refs = {}
        for c in calls.values():
            for lab, tax, gen in ((c[f"{arm}_best_label"], c[f"{arm}_best_taxon"], c[f"{arm}_best_genus"]),
                                  (c[f"{arm}_label"], "", "")):
                if lab and (f"{arm}:{lab}" not in refs or not refs[f"{arm}:{lab}"]["taxon"]):
                    refs[f"{arm}:{lab}"] = {"acc": f"replay_{arm}", "label": lab, "taxon": tax, "genus": gen}
        refs["__marker__"] = {"acc": "replay_marker", "label": "__marker__", "taxon": "", "genus": ""}
        return refs, "", len(refs)

    def run_blast(query, ref_fasta, workdir, name, task="megablast", threads=4, **kw):
        state["arm"] = name.upper()
        return {x: [{"ref": "__marker__", "bitscore": 0.0, "pident": 0.0, "q": x}] for x in calls}

    def assign(hs, refs, min_id, min_cov, max_ev, drop_accs=None):
        if not hs:
            raise SystemExit("replay: a copy without a stored call")
        c = calls[hs[0]["q"]]
        return {f: c[f"{state['arm']}_{f}"] for f in cw.FIELDS}

    cw.load_reference, cw.run_blast, cw.assign = load_reference, run_blast, assign
    cw.sha256 = lambda path, n=None: "replay"
    real_shutil = cw.shutil
    cw.shutil = types.SimpleNamespace(which=lambda t: "/replay", rmtree=real_shutil.rmtree)
    argv = ["crosswalk.py", "--ref", a.ref, "--maarjam", a.calls, "--unite", a.calls, "--lsu", a.calls,
            "--out", a.out]
    buf = io.StringIO()
    old = sys.argv
    sys.argv = argv
    try:
        with contextlib.redirect_stdout(buf):
            cw.main()
    finally:
        sys.argv = old
    log = buf.getvalue()
    head = (f"REPLAY of crosswalk.py (sha256 {sha256(os.path.join(a.code, 'crosswalk.py'))}) on stored per-copy "
            f"calls {os.path.basename(a.calls)} (sha256 {sha256(a.calls)})\n"
            f"reference {man['version']} fingerprint {man['fingerprint']}; {len(changed)} copies changed label or class\n"
            "The arm table's record and unit counts and the 'dropped' columns describe the replay stand-in, not the "
            "external references; everything from section 1 on is the crosswalk report.\n\n")
    with open(a.out + ".log", "w") as fh:
        fh.write(head + log)
    run = json.load(open(a.out + ".run.json"))
    run["replay"] = {"calls": os.path.abspath(a.calls), "calls_sha256": sha256(a.calls),
                     "crosswalk_sha256": sha256(os.path.join(a.code, "crosswalk.py")),
                     "arms_sha256": sha256(os.path.join(a.code, "arms.py")),
                     "changed_copies": changed, "evidence_check": ev_check}
    run["arms"] = "replayed from stored calls; see replay"
    with open(a.out + ".replay.json", "w") as fh:
        json.dump(run, fh, indent=1)
    os.remove(a.out + ".run.json")
    print(head + log)
    print("changed copies:", len(changed))
    for c in changed:
        print("  ", c)
    for arm, v in ev_check.items():
        print(f"evidence check {arm}: {v}")


if __name__ == "__main__":
    main()
