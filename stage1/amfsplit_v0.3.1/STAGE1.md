# amfsplit — Stage 1 (reuse the stack)

## v0.3 / v0.3.1: the queued fixes (30 September 2026)

All four were identified in the Stefani and ONT checks and are now in `amfsplit.py`. `tests/test_v03.py` builds synthetic reads and checks each one against the v0.2 script (`amfsplit_v02_orig.py`), and `compare_records.py` compares two real `records.tsv` files.

- **Region-limited primer search.** Each primer is searched only where it can lie: SSU primers before the ITSxRust full-ITS start, LSU primers after its end, both widened by `--region-margin` (default 50 nt). Reads without a full-ITS anchor fall back to the whole read (`region_limited = no`). Regions follow the window name (`SSU_*`, `LSU_*`, `AMPLICON_*`); optional panel columns `fwd_region` and `rev_region` (`ssu`, `lsu`, `any`) override them. Before this, a spurious FLR4 hit inside the SSU flank produced a window that crossed the ITS. `--no-region-limit` restores v0.2.
- **`missing_5p` / `missing_3p`.** A primer that is absent while the read covers the place where it should lie now gets its own status. `open_5p` / `open_3p` are kept for sites beyond the read: the search region is cut by the read end and the window would still be shorter than the panel's `max_len`. Before, an AM1 failure on a full-length copy came out as `too_long`, and an absent LROR gave a window that started at base 1.
- **Repeat (concatemer) screen.** A read is flagged when at least `--concat-min-primers` different primers (default 2) each occur twice, at least `--concat-min-sep` apart (default 150), with the *same spacing* (a tandem repeat unit, agreeing within `--concat-tol`, 5% and at least 40 nt), or each occur on the opposite strand with the same mirror sum (an inverted repeat). Chance primer hits can agree on a spacing by accident, so a cluster of only two primers must also pass a primer-free check: the two supposed copies (up to 800 nt) must resemble each other, at most `--concat-sim` edits per base (default 0.40; unrelated DNA is about 0.5). A cluster of `--concat-strong` (default 3) or more primers agreeing on one spacing is accepted without it, because the two units of a real concatemer are often different molecules (other taxa or rDNA variants) whose ITS does not resemble each other. `concat_kind` says whether two comparable copies are present (`full`: spacing at least 0.4 of the read, or mirror centre in the middle third) or only a segment is repeated (`partial`: typically an end duplicated or folded back, which makes the read 10-15% longer). It searches the whole read and does not depend on the region limit, so nothing simply disappears. Flagged reads stay in `records.tsv` (`concatemer = yes`, `concat_evidence` with the spacing and copy similarity, `concat_kind`, every window `concatemer`) and are left out of the FASTA files unless `--keep-concatemers` is set. Synthetic reads with 3-15% error: 100% of tandem and inverted repeats flagged, 0% of single copies; single copies carrying 2-6 chance primer sites: 0 of 3,000 flagged (16.7% without the copy-similarity check).
  - **Why v0.3.1, and what the real data show.** v0.3.0 flagged any read in which two primers each occurred twice. On four benchmark libraries that flagged 0.01-0.02% of HiFi and 0.11-0.22% of ONT reads, but the flagged ONT reads were single-copy length (median 1.11x the unflagged median, none at 1.6x or more): chance primer hits. v0.3.1 requires the spacing to agree across primers and the two supposed copies to resemble each other. On the same four libraries it flags 1, 3, 6 and 8 reads (0.00-0.03%). The ONT `.sub` benchmark files are length-selected (nothing below 2.5 kb or above 3.2 kb), so they hold almost no concatemers by construction. Positive control on the unfiltered SRR27457425 library: 1,480 of 244,012 reads are 4.5 kb or longer (0.61%); run without anchors, the screen flagged 811 of them (54.8%; 790 full, 21 partial), which is 0.33% of the whole library and matches the 0.3-0.4% estimated earlier from unfiltered reads. Of the 669 unflagged long reads, 562 (84%) yielded two or more usable windows and had a median length of 5.4 kb, about twice a single amplicon: AMF concatemers that the similarity check had rejected (simulated concatemers of two different molecules were detected only 46-82% of the time with the similarity check alone, close to the 59% seen here). Version 0.3.1 as first released therefore missed about 41% of real ONT concatemers; the strong-cluster rule above restores 100% on the simulation (0.8% false flags in a stress test with 2-6 chance primer sites per read; 0% without decoys). `full` and `partial` describe the geometry of the repeat (comparable units versus a repeated segment), not the completeness of the amplicon: a read followed by the reverse complement of a truncated copy counts as `full` although it is only about 1.2x the median length. `check_concat_lengths.py` reruns the length check on any records.tsv.
- **LF line endings** in `records.tsv` (the csv default is CRLF). The other writers already used `lineterminator`.

`records.tsv` gains three columns at the end (`region_limited`, `concatemer`, `concat_evidence`); existing readers that select columns by name are unaffected. The summary JSON gains `concatemers`, `region_limited_reads` and the new parameters.

**Before using v0.3 on anything that feeds Stage 3**, rerun it on a library that already has v0.2 records and run `python3 compare_records.py OLD NEW`. Usable windows should keep identical coordinates; only statuses (`too_long` to `missing_*`, flagged concatemers) and the two new flags should change.

## v0.2: changes after the first benchmark run (PRJDB40262, PRJNA1062293)

**Triage is now its own step.** In the first run, the best LSU hit for most anchored reads was the non-AMF *Mortierella* outgroup, in three of five soil libraries:

| Library | Best hit Mortierella | Best hit AMF (≥90%) |
|---|---|---|
| DRR907033 | 88.7% | 8.9% |
| SRR27457425 (ASS1) | 93.0% | 0.1% |
| SRR27457429 (N1) | 95.5% | 2.6% |

These fungi anchor perfectly well in ITSxRust, so anchoring can delimit regions but cannot decide taxonomy. `triage.py` labels every read (AMF / unresolved / non_AMF / no_hit / no_window) from its LSU window's best hit, then writes `*.amf.fasta` subsets. By default it uses the Delavaux database's four outgroups; for real use, supply a broad LSU reference with lineage in its headers (`TRIAGE_DB`, e.g. EUKARYOME LSU).

**Tolerant anchoring is the default.** In DRR907034, 57% of HiFi reads got only partial chains under the strict HiFi thresholds, and 94% of those became full-chain with tolerant thresholds (E ≤ 1e-3, score ≥ 15). Most of them best-hit nothing closer than 76% in the Delavaux database, so they are a divergent non-AMF lineage. The strict thresholds were acting as an accidental filter, and triage now does that job explicitly. Set `ANCHOR_MODE=strict` to restore ITSxRust's preset thresholds.

**Open-ended windows (`--open-ends`).** The nested amplicons in both benchmark projects start about 350 nt downstream of AML1, just before AMV4.5NF, so NS31 and WANDA lie outside the read. amfsplit now cuts those windows to the read end and flags them `open_5p` / `open_3p`. On AAFC reads trimmed to mimic this, the SSU flank matched the benchmark (1,194 bp), and NS31–AM1 came back as a 460 bp open window, 91% of the full region.

**New defaults in `run_stage1.sh`.** Tolerant anchoring, `--only-anchored --open-ends`, and triage right after windowing. `MIN_ID` (the identity cut-off for a confident AMF call) defaults to 90 for HiFi and 85 for ONT R9.4.1. Validation against Delavaux and MaarjAM now runs on the AMF-only subsets.


Stage 1 turns each long AMF rDNA read into a **multi-marker record**:

- **ITS regions and flanks from ITSxRust.** ITS1, 5.8S, ITS2, full ITS, and the SSU and LSU flanks come from `--anchors-tsv`, so there is only one nhmmer pass.
- **Legacy barcode windows.** These are the exact short-read regions that the MaarjAM, SSU-V4, LSU-D2 and Delavaux LSU literature was built on, cut from the long read at the original primer sites. The search is IUPAC-aware and tolerates mismatches (edlib infix alignment). Because each window matches the region its legacy database was built from, it stays directly comparable with that database. That is what the Stage 3 VT ↔ SH ↔ LSU-clade bridge needs.

Coordinates are 1-based, inclusive, and given in the read's rRNA-sense orientation. This is the convention ITSxRust uses for minus-strand reads; it was verified against ITSxRust's own ITS1 output.

## Files

| File | Purpose |
|---|---|
| `amfsplit.py` | Per-read records (`.records.tsv`), one FASTA per region/window, `.summary.json` |
| `panels/amf_legacy_windows.tsv` | Primer-window panel: window, target database, primers, loose length range |
| `validate_windows.py` | Best vsearch hit of each window against a legacy reference: identity, query/target coverage |
| `triage.py` | AMF / non-AMF call per read from window best hits; writes AMF-only window FASTAs |
| `run_stage1.sh` | One-sample driver: ITSxRust (all + ssu + lsu from one tblout), amfsplit, triage, then validation on the AMF subsets |
| `results/` | Summaries from the AAFC test run below |

## Quick start

```bash
HMM=/path/to/F.hmm DELAVAUX_DB=/path/to/LRORFLR2.V10seqs_8.4.20.fasta THREADS=8 \
  ./run_stage1.sh sample.fq.gz sample hifi out/
# or, with existing ITSxRust anchors:
python3 amfsplit.py -i sample.fq.gz -a sample.anchors.tsv -p panels/amf_legacy_windows.tsv \
  -o out/sample --only-anchored --open-ends
```

Requirements: `itsxrust` ≥ 0.3.0, HMMER 3 (`nhmmer`), `vsearch`, Python 3 with `edlib`.

## v0.1 results: AAFC PacBio HiFi test reads

The input was the two example datasets shipped with the AAFC AMF-rDNA-PacBio pipeline: `bc1013` and `bc1026`, AML1–LSUmBr amplicons from *Rhizophagus irregularis* spores. The tools were ITSxRust v0.3.0 (release binary, `--preset hifi`), the ITSx fungal `F.hmm`, HMMER 3.4, vsearch 2.27 and edlib 1.3.9.

| | bc1013 | bc1026 |
|---|---|---|
| Reads | 1,707 | 2,137 |
| ITSxRust confident / partial / skipped | 1,275 / 48 / 384 | 2,019 / 118 / 0 |
| ITS1 / 5.8S / ITS2 median length (bp) | 109 / 156 / 202 | 109 / 156 / 201 |
| Legacy windows recovered (anchored reads) | 98.7–99.9% | 99.2–99.9% |
| Windows crossing the HMM ITS anchors | 0 | 0 |
| LROR–FLR2 vs Delavaux DB: hit / median id / target cov | 100% / 97.2% / 100% | 100% / 97.1% / 100% |

The table below gives median window lengths in bp, primer-trimmed, for bc1026. Adding the primers back gives about 549 bp for NS31–AM1, in line with the ~550 bp fragment behind MaarjAM VTs.

| Window | Median length |
|---|---|
| AML1–LSUmBr | 2,766 |
| NS31–AM1 | 508 |
| NS31–AML2 | 522 |
| WANDA–AML2 | 500 |
| AMV4.5NF–AMDGR | 216 |
| LROR–FLR2 | 739 |
| FLR3–FLR4 | 332 |

### What the test run showed

1. **ITSxRust anchors Glomeromycota HiFi reads without tuning.** Every bc1026 read was kept, and in bc1013 the only rejections were a separate ~2.9 kb read population.
2. **Those 384 bc1013 rejections are plant-derived.** Their best LSU hit is the *Citrus limon* outgroup in the Delavaux DB, at about 86% identity and 21% coverage, which fits host-root co-amplification.
   - These reads still carry sites for the less specific SSU primers (NS31–AML2, WANDA–AML2, AMV4.5NF–AMDGR). Without `--only-anchored`, plant SSU windows would flow into MaarjAM classification.
   - The AMF-specific sites (AM1, FLR2, FLR3/FLR4) are absent from them, as the primer designs predict. `--only-anchored` should be the production default until a Glomeromycota triage step exists.
3. **SSUplex is not the right triage for AMF-specific amplicons.** Host 18S and AMF 18S both come out as *eukaryota*; the fungal ITS anchors are what separate them. SSUplex stays relevant for universal-primer soil and root libraries with organelle or bacterial carry-over.
4. **Primer-cut windows match the legacy database region better than HMM flanks.** LROR–FLR2 windows reach 93% median query coverage with full target coverage, against 77% for the whole LSU flank.
5. **Stage 2 preview.** Within a single *R. irregularis* culture, the LSU window's identity to its best reference has quartiles of 95.5 / 97.1 / 97.4% (bc1026), and about 40% of reads fall below 97%. Part of this is strain-versus-reference divergence, but a 97% threshold would still split this one organism.
6. **Flag for Stage 3.** 20 bc1026 LSU windows (1%) best-match *Scutellospora dipurpurascens* at 94–97%. This is either a second lineage or chimeras, and exactly what a cross-region concordance check should catch.

## Remaining Stage 1 work (on the workstation)

1. **Benchmark data**

   ```bash
   for acc in PRJDB40262 PRJNA1062293; do
     curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=run_accession,fastq_ftp&format=tsv" > $acc.runs.tsv
     tail -n +2 $acc.runs.tsv | cut -f2 | tr ';' '\n' | sed 's#^#https://#' | xargs -n1 -P4 curl -sO
   done
   ```

   PRJDB40262 is PacBio Revio HiFi (`hifi`). PRJNA1062293 is ONT R9.4.1 (`ont`); `run_stage1.sh` raises the primer edit budget to 0.20 for ONT. Report window recovery as a function of `--max-err-frac`.

2. **Validation against references not reachable from here**
   - **MaarjAM.** Compare the `SSU_VT_NS31-AM1` windows against the VT reference sequences by setting `MAARJAM_DB`. Expect near-full target coverage and clean VT assignment.
   - **Krüger et al. 2012 reference alignment.** Run ITSxRust and amfsplit on the reference sequences themselves. This tests primer sites and anchors across lineages, especially Archaeosporales and Paraglomerales, where AM1 and `F.hmm` are most likely to fail.
   - **EUKARYOME long-read SSU–ITS–LSU subset (Glomeromycota).** Build a boundary-error distribution against reference annotations, reusing ITSxRust's `bench/sim/04_eval_boundaries.py`.

3. **Fallback for primer-site failures.** Project window coordinates by aligning to a reference operon wherever the primer search misses.

4. **Rust port.** Once the panel and window definitions settle, port to Rust as a separate crate that consumes ITSxRust anchors, which keeps ITSxRust's published scope intact. rust-bio's Myers matcher with ambiguity codes is the likely replacement for edlib.

**Proposed Stage 1 exit criteria**
- At least 95% of Glomeromycota HiFi reads anchored.
- At least 95% window recovery wherever the primer is expected to bind.
- Zero windows crossing the HMM anchors.
- LSU and VT windows covering at least 95% of their reference entries.

Primer sequences in the panel were entered from the primary literature. Confirm each against its source paper before publication; the near-complete site recovery on *R. irregularis* is a good sign but does not replace that check.
