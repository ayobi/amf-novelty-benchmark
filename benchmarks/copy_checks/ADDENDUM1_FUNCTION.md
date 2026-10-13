# Addendum 1: are the divergent copies functional?

Written 2026-10-12, before any result of these analyses was computed, in answer to L. Tedersoo's comment II.1
(are the short-5.8S copies still functional; irregularities in SSU and LSU; secondary structure). It adds
analyses to PROTOCOL_COPY_CHECKS.md and changes nothing in it. Same reference (culture reference v4), same
comparators (the nearest main-type copy of the same culture, whole-copy distance) and the same percentile
method (exclusive, 99th) as section 2 of the protocol.

## Groups

- **Short 5.8S**: the 43 short-5.8S copies.
- **Same length**: the 8 same-length divergent copies.
- **Flagged main**: main-type copies with a long deletion in the SSU flank, 5.8S or LSU flank in the artefact
  check (listed in `artefacts/artefacts.json`).
- **Baseline**: all other main-type copies with a main-type culture-mate.

## (a) Irregularities in SSU and LSU

From the alignment of the artefact check (each copy against its comparator; BLASTN scores, free end gaps),
per region (SSU flank, 5.8S, LSU flank): substitutions per 100 aligned positions and the number of internal
indel events. A copy is **irregular** in a region when either exceeds the 99th percentile of the baseline.

## (b) Secondary structure

Each region is aligned with Infernal 1.1.5 `cmalign` (truncated alignment allowed) to its Rfam model:
RF01960 (eukaryotic SSU) for the SSU flank, RF00002 (5.8S) for the 5.8S and RF02543 (eukaryotic LSU) for the
LSU flank. For every consensus base pair the copy covers on both sides, the pair is **broken** when it is not
Watson-Crick (A-U, G-C) or G-U, or when one side is deleted. Per copy and region: the broken-pair fraction and
the bit score per nucleotide.

A copy is **degraded** in a region when its broken-pair fraction exceeds the 99th percentile of the baseline.
For each group and region: the share of degraded copies, and a paired comparison of each divergent copy with
the mean of its culture's main-type copies (two-sided Wilcoxon signed-rank test).

The 5.8S deletions of the short copies are located on the RF00002 consensus: which consensus positions they
remove and whether those positions are paired or unpaired in the consensus structure.

**Verdict for the short-5.8S class:**
- *structurally intact* if at most 10% of its copies are degraded in each region and its shared deletion
  removes only unpaired consensus positions;
- *structurally degraded* (consistent with pseudogenes) if at least 50% are degraded in some region, or the
  shared deletion removes paired consensus positions;
- *mixed* otherwise.
The same rule is applied to the same-length group, descriptively.

## (c) ITS2 folding (descriptive)

ITS2 of each copy folded with ViennaRNA 2.7.2 `RNAfold` (default parameters, 37 C): minimum free energy per
nucleotide, and the base-pair distance between the copy's structure and its comparator's. ITS2 is a spacer, so
no verdict is drawn from it.

## Outputs (`benchmarks/copy_checks/function/`)

`function.copies.tsv` (every copy: per-region substitutions, indels, broken pairs, bit score, ITS2 folding),
`function.summary.tsv` (per group and region), `function.deletions.tsv` (the 5.8S deletions on the RF00002
consensus), `function.json` (versions, model checksums, thresholds, hashes).
