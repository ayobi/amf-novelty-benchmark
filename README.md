# Detecting unknown arbuscular mycorrhizal fungi from linked long-read rDNA copies

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23219691.svg)](https://doi.org/10.5281/zenodo.23219691)

Code, curated reference and benchmark outputs for the paper of the same title (O'Brien, in preparation).

Each long-read rDNA copy from an identified AMF culture spans SSU, ITS and LSU. Cut into the barcode
windows of MaarjAM (SSU virtual taxa), UNITE (ITS species hypotheses) and the Delavaux LSU database,
the copies link the three naming systems on the same molecules and give a test of how reliably a
combined call names known species and flags unknown ones.

## Contents

| Folder | What it holds |
|---|---|
| `stage1/amfsplit_v0.3.1/` | The Stage 1 window caller used to cut each copy into barcode windows (amfsplit v0.3.1) |
| `stage1/prototypes/` | Stage 1 prototypes: GenBank grouping, QC flags, VT checks, primer panels |
| `stage3/` | Stage 3 code (`arms.py`, `freeze_ref.py`, `crosswalk.py`, `reconcile.py`, `novelty.py`, `longglodb_overlap.py`, `predict_longglodb.py`) and the reference edits files with their reasons |
| `reference/culture_ref_v4/` | Culture reference v4 (fingerprint `181938d695f3df40`): 903 copies, 151 cultures, 37 named species, with its manifest |
| `benchmarks/stefani/` | Crosswalk and holdout outputs on the Stefani et al. (2025) copies |
| `benchmarks/m4_external_holdout/` | External holdout (protocols v1 and v2), the frozen M2b code with its LSU placement, and the LongGloDB runner and report (`lgdb_test.py`, `lgdb_report.py`) with their per-fold outputs |
| `benchmarks/longglodb/` | Independent-culture test: protocol, amendments, locked predictions and their hashes, score tables. No LongGloDB sequences (see below) |
| `paper/` | LaTeX source; `paper/analysis/` regenerates every number and generated table from the run outputs |
| `external_refs/` | Checksums of the third-party references used (not redistributed) |
| `env/` | Conda environments (`amfsplit`, `amf-m4`) |
| `MANIFEST.tsv` | Every file with its size and SHA-256 |

## Data

- **Culture copies:** GenBank PX207096–PX207274 and PX214440–PX215164 (Stefani et al. 2025,
  *New Phytologist* 248: 1501–1515, doi:10.1111/nph.70557). The culture reference in this repository
  is derived from those records; the edits files record every change and its reason.
- **Third-party references, not included.** Download them from their sources and check them against
  `external_refs/CHECKSUMS.sha256`:
  - MaarjAM: gDAT snapshot, `github.com/ut-planteco/gDAT` (`db/`), trimmed to NS31–AML2
    (sha256 `2b2c57e664e76347…`).
  - UNITE general FASTA release 10.0 of 19 February 2025, doi:10.15156/BIO/3301229.
  - Delavaux LSU database V18 (AMF only), `github.com/c383d893/AMF-LSU-Database-and-Pipeline2`,
    `2024_AMFPipeline/2024_AMFPipeline_ASV/V18_LSUDB_052025_AMFONLY.fasta`.
  - LongGloDB v1 (July 2026), from `glodb.org` (Databases page). The test set is rebuilt from the
    public files by `stage3/longglodb_overlap.py`; no LongGloDB sequence is stored here.

## Reproducing

Outline, in run order; each script takes its inputs as command-line options:

1. **Stage 1:** cut copies into windows with `stage1/amfsplit_v0.3.1/run_stage1.sh` (needs the ITSxRust
   HMM profiles and the primer panel in `stage1/prototypes/panels/`).
2. **Reference:** `stage3/freeze_ref.py` builds the culture reference from the Stage 1 tables and an
   edits file; v4 is the version used in the paper.
3. **Crosswalk:** `stage3/crosswalk.py` assigns VT, SH and LSU labels per copy against the three references.
4. **Holdouts:** culture, species and external holdouts, run under the protocols in
   `benchmarks/m4_external_holdout/`.
5. **Independent cultures:** `stage3/longglodb_overlap.py`, the locked prediction (`stage3/predict_longglodb.py`), then
   `benchmarks/m4_external_holdout/lgdb_test.py` and `lgdb_report.py`,
   following `benchmarks/longglodb/PROTOCOL_LGDB.md` and its amendments.
6. **Paper:** `paper/analysis/make_numbers.py`, `lgdb_numbers.py` and `make_figures.py` regenerate the
   numbers and generated tables; then `pdflatex main.tex` twice.

## Licence

Code: MIT (`LICENSE`). Tables and other outputs produced here: CC BY 4.0 (`DATA_LICENSE.md`).
Third-party data keep their own terms.

## Citation

O'Brien, A. Detecting unknown arbuscular mycorrhizal fungi from linked long-read rDNA copies.
In preparation. Code and data: v0.1.0, doi:10.5281/zenodo.23219691.
