# Detecting unknown arbuscular mycorrhizal fungi from linked long-read rRNA gene sequences

Code, curated reference and benchmark outputs for the paper of the same title (O'Brien, Marín and Tedersoo,
preprint). Release 0.3.1 goes with the preprint; the paper itself is in `paper/` (`main.tex`, `main.pdf`).

Each long-read rRNA gene sequence from an identified AMF culture spans the SSU, ITS and LSU. Cut into the
regions of MaarjAM (SSU virtual taxa), UNITE (ITS species hypotheses) and the Delavaux LSU database, the
sequences link the three naming systems on the same molecules. They test how reliably the three names,
combined, identify known species and flag species missing from the reference.

## Contents

| Folder | What it holds |
|---|---|
| `stage1/amfsplit_v0.3.1/` | The Stage 1 code that cuts each sequence into marker regions (v0.3.1) |
| `stage1/prototypes/` | Stage 1 prototypes: GenBank grouping, QC flags, VT checks, primer panels |
| `stage3/` | Stage 3 code (`arms.py`, `freeze_ref.py`, `crosswalk.py`, `reconcile.py`, `novelty.py`, `longglodb_overlap.py`, `predict_longglodb.py`) and the reference edits files with their reasons |
| `reference/culture_ref_v4/` | Culture reference v4 (fingerprint `181938d695f3df40`): 903 sequences, 151 cultures, 37 named species, with its manifest |
| `benchmarks/stefani/` | Crosswalk and holdout outputs on the Stefani et al. (2025) sequences |
| `benchmarks/m4_external_holdout/` | LSU-removal test (protocols v1 and v2), the frozen M2b code with its LSU placement, the LongGloDB runner and report (`lgdb_test.py`, `lgdb_report.py`) and the EUKARYOME wrappers (`euk_test.py`, `euk_report.py`), with their outputs |
| `benchmarks/longglodb/` | Independent-culture test: protocol, amendments, locked predictions and their hashes, score tables. No LongGloDB sequences (see below) |
| `benchmarks/eukaryome/` | EUKARYOME long-read test: protocol and note with their hashes, name table and decisions, roles, locked predictions and removal plan. Read IDs only, no sequences |
| `benchmarks/copy_checks/` | Checks of divergent sequences (locked protocol and addendum, with their hashes): per-marker trees (Supplementary Data S1), artifact check, rRNA structure, distance thresholds on short and long regions, the effect of the short-5.8S paralog on clusters and trees, the motif that marks it, and the controls of the SSU-discordance flag (`artifact_control/`). Scripts: `stage3/copy_*.py`, `stage3/short58_motif*.py`, `stage3/artifact_control.py` |
| `benchmarks/eukaryome_longread/` | Exploratory, not used in the paper: distance thresholds on whole EUKARYOME v2.1 long reads of named Glomeromycota (`stage3/longread_thresholds.py`). Read IDs and distances only, no sequences |
| `paper/` | LaTeX source and PDF; `paper/analysis/` regenerates every number, generated table and figure from the run outputs; `paper/review/` holds the response to the first review |
| `external_refs/` | Checksums of the third-party references used (not redistributed) |
| `MANIFEST.tsv` | Every file of release v0.1.0 with its size and SHA-256 |

## Data

- **Culture sequences:** GenBank PX207096–PX207274 and PX214440–PX215164 (Stefani et al. 2025,
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

1. **Stage 1:** cut the sequences into marker regions with `stage1/amfsplit_v0.3.1/run_stage1.sh` (needs the ITSxRust
   HMM profiles and the primer panel in `stage1/prototypes/panels/`).
2. **Reference:** `stage3/freeze_ref.py` builds the culture reference from the Stage 1 tables and an
   edits file; v4 is the version used in the paper.
3. **Crosswalk:** `stage3/crosswalk.py` names each sequence by VT, SH and LSU against the three references.
4. **Withholding tests:** cultures, species, and species removed from the LSU reference, run under the protocols in
   `benchmarks/m4_external_holdout/`.
5. **Independent cultures:** `stage3/longglodb_overlap.py`, the locked prediction (`stage3/predict_longglodb.py`), then
   `benchmarks/m4_external_holdout/lgdb_test.py` and `lgdb_report.py`,
   following `benchmarks/longglodb/PROTOCOL_LGDB.md` and its amendments.
6. **EUKARYOME:** `stage3/eukaryome_overlap.py`, then `euk_test.py` and `euk_report.py` in
   `benchmarks/m4_external_holdout/`, following `benchmarks/eukaryome/PROTOCOL_EUK.md`.
7. **Paper:** `paper/analysis/make_numbers.py`, `lgdb_numbers.py`, `euk_numbers.py`, `foldwise_complexes.py`,
   `intervals.py`, `divergent_classes.py`, `make_figures.py` and `make_fig_predicted.py` regenerate the
   numbers, generated tables and figures; `paper/analysis/display_names.py` then prints the names and terms used in the
   paper (it changes no number); then `pdflatex main.tex` twice.

## Licence

Code: MIT (`LICENSE`). Tables and other outputs produced here: CC BY 4.0 (`DATA_LICENSE.md`).
Third-party data keep their own terms.

## Citation

O'Brien A, Marín C, Tedersoo L. Detecting unknown arbuscular mycorrhizal fungi from linked long-read rRNA gene
sequences. Preprint.
Code and data: release 0.3.1 (Zenodo); earlier releases 0.3.0, doi:10.5281/zenodo.23422391 (paper v19, before the
EUKARYOME results), 0.2.0, doi:10.5281/zenodo.23374725, and 0.1.0, doi:10.5281/zenodo.23219691. `CITATION.cff` holds the same in machine-readable form.
