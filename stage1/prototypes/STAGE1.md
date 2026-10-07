# amfsplit — Stage 1 (reuse the stack)

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
