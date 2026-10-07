#!/usr/bin/env python3
"""Synthetic tests for amfsplit v0.3 against the v0.2 script (amfsplit_v02_orig.py).

Reads are built from random sequence with made-up primers, so the tests check logic, not primer
biology. Run: python3 tests/test_v03.py   (from the amfsplit folder; needs edlib)
"""
import csv
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import edlib  # noqa: F401
except ImportError:
    sys.exit("edlib is not installed for this python (" + sys.executable + ").\n"
             "Activate the environment that has it first, e.g.:  conda activate amfsplit")

HERE = Path(__file__).resolve().parent.parent
V02, V03 = HERE / "amfsplit_v02_orig.py", HERE / "amfsplit.py"
rng = random.Random(7)
COMP = str.maketrans("ACGT", "TGCA")


def rnd(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(COMP)[::-1]


def mutate(s, n_edits):
    s = list(s)
    for i in rng.sample(range(len(s)), n_edits):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


P = {name: rnd(20) for name in ("AML1", "NS31", "AMV45", "AM1", "AMDGR", "LROR", "FLR2", "FLR3", "FLR4", "LSUMBR")}
PANEL = [
    ("AMPLICON_AML1-LSUmBr", "AML1", P["AML1"], "LSUmBr", P["LSUMBR"], 1800, 3500),
    ("SSU_VT_NS31-AM1", "NS31", P["NS31"], "AM1", P["AM1"], 400, 650),
    ("SSU_V4_AMV4.5NF-AMDGR", "AMV4.5NF", P["AMV45"], "AMDGR", P["AMDGR"], 150, 300),
    ("LSU_LROR-FLR2", "LROR", P["LROR"], "FLR2", P["FLR2"], 500, 900),
    ("LSU_D2_FLR3-FLR4", "FLR3", P["FLR3"], "FLR4", P["FLR4"], 250, 450),
]


def put(seq, pos, site):
    return seq[:pos] + site + seq[pos + len(site):]


def make_copy():
    """SSU 1600 + ITS 465 + LSU 1100, sites placed in sense orientation."""
    ssu, its, lsu = rnd(1600), rnd(465), rnd(1100)
    ssu = put(ssu, 20, P["AML1"])
    ssu = put(ssu, 1000, P["NS31"])
    ssu = put(ssu, 1100, P["AMV45"])
    ssu = put(ssu, 1300, rc(P["AMDGR"]))
    ssu = put(ssu, 1500, rc(P["AM1"]))
    lsu = put(lsu, 40, P["LROR"])
    lsu = put(lsu, 500, P["FLR3"])
    lsu = put(lsu, 780, rc(P["FLR2"]))
    lsu = put(lsu, 900, rc(P["FLR4"]))
    lsu = put(lsu, 1050, rc(P["LSUMBR"]))
    return ssu, its, lsu


def anchors_row(rid, n_ssu=1600, n_its=465):
    a, b = n_ssu + 1, n_ssu + n_its
    return {"read_id": rid, "its1_start": a, "its1_end": a + 108, "its2_start": a + 265, "its2_end": b,
            "full_start": a, "full_end": b, "confidence": "confident", "ssu_strand": "+", "s58s_strand": "+",
            "s58e_strand": "+", "lsu_strand": "+"}


def build():
    ssu, its, lsu = make_copy()
    full = ssu + its + lsu
    reads, anchors = {}, {}

    def add(name, seq, anchored=True, anc=None):
        reads[name] = seq
        if anchored:
            anchors[name] = anc or anchors_row(name)

    add("A_clean", full)
    # B: AM1 site destroyed, read covers where it should be
    add("B_am1_absent", put(full, 1500, mutate(rc(P["AM1"]), 9)))
    # C: as B plus an exact AM1 decoy inside the LSU, which the v0.2 whole-read search prefers
    add("C_am1_absent_decoy", put(put(full, 1500, mutate(rc(P["AM1"]), 9)), 1600 + 465 + 600, rc(P["AM1"])))
    # D: nested-amplicon read starting inside the NS31 site (site beyond the read)
    cut = 1010
    add("D_nested_5p", full[cut:], anc=anchors_row("D_nested_5p", 1600 - cut))
    # E: read ends before FLR2 (site beyond the read at the 3' end)
    add("E_trunc_3p", full[:1600 + 465 + 700])
    # F: LROR destroyed while the read continues past where it should lie
    add("F_lror_absent", put(full, 1600 + 465 + 40, mutate(P["LROR"], 9)))
    # G: two copies in the same orientation
    # a concatemer holds copies of the SAME molecule: the second copy is the first with ~1% differences
    add("G_concat_tandem", full + mutate(full, 30), anc=anchors_row("G_concat_tandem"))
    # H: head-to-head inverted concatemer
    add("H_concat_inverted", full + rc(full), anc=anchors_row("H_concat_inverted"))
    # I: unanchored read
    add("I_unanchored", full, anchored=False)
    # J: off-target
    add("J_offtarget", rnd(3200), anchored=False)
    # K: single copy with spurious repeat sites of four primers at unrelated spacings (offsets to the
    #    real site: 850, 665, 1280 and 2445 nt), which v0.3.0 called a concatemer: repeated primers
    #    but no shared copy spacing
    noisy = put(full, 150, mutate(P["NS31"], 3))
    noisy = put(noisy, 1900, mutate(P["FLR3"], 3))
    noisy = put(put(noisy, 1300, mutate(P["AML1"], 3)), 2300, mutate(P["AML1"], 3))
    noisy = put(noisy, 400, mutate(rc(P["FLR2"]), 3))
    add("K_spurious_repeats", noisy)
    # N: a concatemer of two DIFFERENT molecules (other taxon: flanks 5-8% diverged, unrelated ITS)
    ssu_, its_, lsu_ = full[:1600], full[1600:2065], full[2065:]
    other = mutate(ssu_, 80) + rnd(len(its_)) + mutate(lsu_, 88)
    add("N_concat_different_molecule", full + other, anc=anchors_row("N_concat_different_molecule"))
    # L / M: a single copy with its last 330 nt repeated (tandem) or folded back (inverted): partial repeats
    tail = full[-330:]
    add("L_partial_tandem", full + tail, anc=anchors_row("L_partial_tandem"))
    add("M_partial_foldback", full + rc(tail), anc=anchors_row("M_partial_foldback"))
    return reads, anchors


def write_inputs(d, reads, anchors):
    with open(d / "reads.fa", "w") as fh:
        for k, v in reads.items():
            fh.write(f">{k}\n{v}\n")
    cols = ["read_id", "its1_start", "its1_end", "its2_start", "its2_end", "full_start", "full_end",
            "confidence", "ssu_strand", "s58s_strand", "s58e_strand", "lsu_strand"]
    with open(d / "anchors.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(anchors.values())
    with open(d / "panel.tsv", "w") as fh:
        fh.write("window\ttarget\tfwd_name\tfwd_seq\trev_name\trev_seq\tmin_len\tmax_len\n")
        for row in PANEL:
            fh.write("\t".join(map(str, row[:1] + ("test",) + row[1:])) + "\n")


def run(script, d, tag, *extra):
    out = d / tag
    r = subprocess.run([sys.executable, str(script), "-i", str(d / "reads.fa"), "-a", str(d / "anchors.tsv"),
                        "-p", str(d / "panel.tsv"), "-o", str(out), "--open-ends", "--max-err-frac", "0.15", *extra],
                       capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"{script.name} failed:\n{r.stderr[-1500:]}")
    raw = (d / f"{tag}.records.tsv").read_bytes()
    rows = {r["read_id"]: r for r in csv.DictReader(raw.decode().splitlines(), delimiter="\t")}
    return rows, raw, json.load(open(d / f"{tag}.summary.json"))


def main():
    reads, anchors = build()
    fails = []

    def check(cond, msg):
        print(("ok   " if cond else "FAIL ") + msg)
        if not cond:
            fails.append(msg)

    with tempfile.TemporaryDirectory() as t:
        d = Path(t)
        write_inputs(d, reads, anchors)
        old, old_raw, _ = run(V02, d, "v02")
        new, new_raw, summ = run(V03, d, "v03")
        st = lambda rows, r, w: rows[r][f"{w}_status"]

        check(b"\r" in old_raw and b"\r" not in new_raw, "v0.2 writes CRLF, v0.3 writes LF")
        # clean read: identical windows
        same = all(old["A_clean"][f"{w[0]}_{k}"] == new["A_clean"][f"{w[0]}_{k}"]
                   for w in PANEL for k in ("status", "start", "end", "len"))
        check(same and all(st(new, "A_clean", w[0]) == "ok" for w in PANEL), "clean read: all windows ok and identical to v0.2")
        check(new["A_clean"]["concatemer"] == "no" and new["A_clean"]["region_limited"] == "yes", "clean read: not a concatemer, region limited")
        # B: relabelling
        check(st(old, "B_am1_absent", "SSU_VT_NS31-AM1") in ("too_long", "open_3p"), f"v0.2 B: {st(old, 'B_am1_absent', 'SSU_VT_NS31-AM1')} (misleading)")
        check(st(new, "B_am1_absent", "SSU_VT_NS31-AM1") == "missing_3p", "v0.3 B: AM1 absent on a covering read is missing_3p")
        # C: decoy
        check(st(new, "C_am1_absent_decoy", "SSU_VT_NS31-AM1") == "missing_3p", "v0.3 C: AM1 decoy in the LSU is ignored")
        check(new["C_am1_absent_decoy"]["flank_consistent"] == "yes", "v0.3 C: no window crosses the ITS")
        check(old["C_am1_absent_decoy"]["flank_consistent"] == "no" or st(old, "C_am1_absent_decoy", "SSU_VT_NS31-AM1") != "missing_3p",
              "v0.2 C: the decoy corrupted the window (motivation for the region limit)")
        # D: nested amplicon: NS31 beyond the read stays open_5p
        check(st(new, "D_nested_5p", "SSU_VT_NS31-AM1") == "open_5p", "v0.3 D: site beyond the read start stays open_5p")
        check(st(new, "D_nested_5p", "AMPLICON_AML1-LSUmBr") == "missing_5p" or st(new, "D_nested_5p", "AMPLICON_AML1-LSUmBr") == "open_5p",
              f"v0.3 D: amplicon window {st(new, 'D_nested_5p', 'AMPLICON_AML1-LSUmBr')}")
        # E: truncated 3'
        check(st(new, "E_trunc_3p", "LSU_LROR-FLR2") == "open_3p" and st(new, "E_trunc_3p", "LSU_D2_FLR3-FLR4") == "open_3p",
              "v0.3 E: read ends before the site: open_3p")
        # F: LROR absent
        check(st(new, "F_lror_absent", "LSU_LROR-FLR2") == "missing_5p", "v0.3 F: LROR absent on a covering read is missing_5p")
        check(st(old, "F_lror_absent", "LSU_LROR-FLR2") in ("too_long", "open_5p"), f"v0.2 F: {st(old, 'F_lror_absent', 'LSU_LROR-FLR2')} (window started at read base 1)")
        # concatemers
        for r in ("G_concat_tandem", "H_concat_inverted", "N_concat_different_molecule", "L_partial_tandem", "M_partial_foldback"):
            check(new[r]["concatemer"] == "yes" and all(st(new, r, w[0]) == "concatemer" for w in PANEL),
                  f"v0.3 {r}: flagged ({new[r]['concat_evidence'][:60]})")
            check(old[r]["windows_ok"] not in ("0", ""), f"v0.2 {r}: no flag, {old[r]['windows_ok']} windows written")
        check(all(new[r]["concatemer"] == "no" for r in ("A_clean", "B_am1_absent", "C_am1_absent_decoy", "D_nested_5p", "E_trunc_3p", "F_lror_absent", "I_unanchored", "J_offtarget", "K_spurious_repeats")),
              "v0.3: no false concatemer among the single-copy reads")
        for r, kind in (("G_concat_tandem", "full"), ("H_concat_inverted", "full"), ("N_concat_different_molecule", "full"),
                        ("L_partial_tandem", "partial"), ("M_partial_foldback", "partial")):
            check(new[r]["concatemer"] == "yes" and new[r]["concat_kind"] == kind,
                  f"v0.3.1 {r}: flagged as {new[r]['concat_kind'] or 'nothing'} repeat ({new[r]['concat_evidence'][-40:]})")
        check(new["K_spurious_repeats"]["concatemer"] == "no",
              f"v0.3.1 K: spurious repeats at unrelated spacings are not a concatemer ({new['K_spurious_repeats']['concat_evidence'] or 'no evidence'})")
        check(summ["concatemers"]["flagged"] == 5 and summ["concatemers"]["full"] == 3 and summ["concatemers"]["partial"] == 2, "summary counts 3 full and 2 partial repeats")
        # flagged reads stay in records but not in the FASTA
        fasta_ids = {l[1:].split()[0] for l in open(d / "v03.SSU_VT_NS31-AM1.fasta") if l.startswith(">")}
        check("G_concat_tandem" not in fasta_ids and "A_clean" in fasta_ids, "flagged concatemers are left out of the window FASTA")
        # unanchored / off-target
        check(new["I_unanchored"]["region_limited"] == "no" and st(new, "I_unanchored", "SSU_VT_NS31-AM1") == "ok",
              "unanchored read: whole-read search, windows found, region_limited = no")
        check(st(new, "J_offtarget", "SSU_VT_NS31-AM1") in ("both_missing", "unoriented"), "off-target read: no window")
        # --keep-concatemers and --no-region-limit
        keep, _, _ = run(V03, d, "keep", "--keep-concatemers")
        check(keep["G_concat_tandem"]["concatemer"] == "yes" and st(keep, "G_concat_tandem", "SSU_VT_NS31-AM1") != "concatemer",
              "--keep-concatemers writes the windows and keeps the flag")
        nolim, _, _ = run(V03, d, "nolim", "--no-region-limit", "--no-concat-screen")
        same_as_old = all(nolim[r][f"{w[0]}_status"] == old[r][f"{w[0]}_status"] and nolim[r][f"{w[0]}_start"] == old[r][f"{w[0]}_start"]
                          for r in ("A_clean", "C_am1_absent_decoy", "D_nested_5p", "E_trunc_3p") for w in PANEL
                          if old[r][f"{w[0]}_status"] not in ("too_long", "fwd_missing", "rev_missing"))
        check(same_as_old, "--no-region-limit --no-concat-screen reproduces v0.2 windows (except relabelled statuses)")

    print(f"\n{len(fails)} failure(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
