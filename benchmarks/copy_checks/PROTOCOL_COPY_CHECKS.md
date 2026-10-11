# Copy checks: per-marker monophyly and divergent-copy artefacts

Written 2026-10-11, before any tree or artefact result was computed. It answers points 1 and 2 of
L. Tedersoo's comments on the paper. Input: culture reference v4 (fingerprint `181938d695f3df40`;
903 retained copies, 151 cultures, 148 physical culture labels, 37 named species) and the crosswalk
on v4 (`paper/analysis/xw/crosswalk_v4final.copies.tsv`). Nothing here changes the reference, the
crosswalk or any result of the paper. Any change to this file after its hash is recorded in
`LOCK_COPY_CHECKS.sha256` is a dated, hashed amendment.

Distances are the paper's: edlib infix edit distance as a percentage of the shorter sequence.
Percentiles are Python's exclusive method (`statistics.quantiles`), as in `species_groups.py`.

## 1. Per-marker trees and monophyly

**Copies and markers.** All 903 retained copies, one tree set per order (Archaeosporales 30 copies,
Paraglomerales 28, Diversisporales 245, Glomerales 600). ITS cannot be aligned across orders, and every
culture and species lies inside one order. Markers: SSU flank (`ssu_flank`), full ITS (`full_its`),
LSU flank (`lsu_flank`, which contains the LROR-FLR2 window) and the linked copy (the three
concatenated). The three regions tile each copy exactly.

**Alignment.** Identical sequences are collapsed per order and marker. MAFFT L-INS-i
(`--localpair --maxiterate 1000`), no trimming. The linked alignment concatenates the three marker
alignments, one partition per marker.

**Trees.** IQ-TREE with ModelFinder (`-m MFP`; for the linked alignment, one model per partition and
edge-linked proportional branch lengths), 1,000 ultrafast bootstraps (`-B 1000`), 1,000 SH-aLRT
replicates (`--alrt 1000`), near-zero branches collapsed (`-czb`), seed 20261011. Each tip is one
unique sequence and stands for every copy that carries it, as a polytomy.

**Groups.** Culture: the 151 curated cultures (the three mixed labels are already split by lineage).
Species: the 37 named species, all cultures pooled. Complex: the four complexes of
`paper/analysis/v4.complexes.tsv`, each as one unit. A group is tested on a tree when it has at least
two copies in the copy set. Two copy sets: all copies, and main-type copies only (the 51
divergent-class copies pruned from the same tree).

**Outcome per group and tree,** read from the bipartitions of the unrooted tree:
- *monophyletic*: a branch separates exactly the group's copies from all others;
- *unresolved*: no branch does, and none contradicts the group (it sits in a polytomy, or shares an
  identical sequence with another group's copy);
- *non-monophyletic*: at least one branch contradicts the group; *strongly* non-monophyletic when a
  contradicting branch has SH-aLRT >= 80 and UFBoot >= 95. A branch that only joins copies with the
  same sequence has no support value and never counts as strong.

For a non-monophyletic group: *escapees* are its copies outside the largest clade made only of its
copies, and *intruders* are other copies inside the smallest clade that holds all of its copies;
both are listed with their copy class.

**Reported,** per marker, pooled over orders and per order: the number and share of cultures and of
named species that are monophyletic, unresolved, non-monophyletic and strongly non-monophyletic, with
all copies and with main-type copies only. Species are also reported for the 30 species with two or
more cultures, the informative case. The answer to the comment is the share of cultures and of named
species that are monophyletic per marker, with all copies.

**Reading.** A culture that is non-monophyletic with all copies and monophyletic with main-type copies
only is split by its divergent-class copies. A species that is non-monophyletic with main-type copies
is split by its cultures, not by its divergent copies.

## 2. Artefact check on the divergent copies

**Targets.** The 51 divergent-class copies (43 short-5.8S, 8 same-length). The 852 main-type copies
give the baseline.

**Deletions.** Each copy is aligned to the nearest main-type copy of its own culture (whole-copy
distance; never itself). Global alignment with free end gaps, Biopython `PairwiseAligner`, match 2,
mismatch -3, gap open -5, gap extend -2 (the BLASTN defaults). Gap runs that touch either end of the
alignment are trimming differences and are ignored. Each internal deletion (comparator bases missing
from the copy) is assigned to the region of the comparator (SSU flank, ITS1, 5.8S, ITS2, LSU flank)
that holds most of its bases; insertions are recorded the same way, descriptively.

A deletion in region R is *long* when it is at least L_R nt, with L_R = max(3, P99_R), where P99_R is
the 99th percentile of the longest deletion in R per main-type copy against its nearest main-type
culture-mate (0 when it has none), over the main-type copies of cultures with two or more of them. The
floor of 3 nt keeps homopolymer-length errors out. A copy is *intact* when it has no long deletion in
any region and *deleted* otherwise.

**Recurrence.** Two copies are the *same type* when their full-ITS distance is <= 1.0%. A divergent copy
*recurs across isolates* when a copy of the same type, of any class, occurs under another physical
culture label (the two lineages of a mixed label count as one isolate).

**Verdict** (the rule of the comment):

| | intact | deleted |
|---|---|---|
| recurring | paralog | recurrent deleted copy (a genomic copy; possibly a pseudogene) |
| single | unresolved | probable artefact |

**Reported with each copy, not used in the verdict:** a *library replicate* (a copy of the same type
from another library of the same culture; 12 cultures have 2-4 libraries); an *external record* (best
UNITE+INSD hit in the crosswalk, copies of the culture reference excluded, at >= 99.0% identity over
>= 90% coverage); whether the closest same-type copy under another label is a main-type copy of
another species (a mixed culture or contamination would look like this); the copy's SH and VT and
whether either is a paralog unit (reached by at most one main-type copy, the paper's Table definition);
the 5.8S covariance-model score (Infernal, RF00002, `benchmarks/stefani/s58.cmsearch.tbl`) minus the
median of its culture's main-type copies. A single copy with a library replicate or an external record
is marked *reproduced*; its verdict does not change.

**Types** (descriptive). Divergent copies are grouped by single linkage at <= 1.0% full ITS. Per type:
copies, cultures, labels, species, libraries, regions with long deletions, SH and VT units, verdicts.

**Sensitivity, fixed now.** Verdict counts with same type at <= 0.5% and <= 2.0%, and with long
deletion read as >= 10 nt in every region.

**By-product.** Main-type copies with a long deletion in the SSU flank, 5.8S or LSU flank are listed.
Any change to the reference that follows is a dated edit of the edits file, not part of this check.

**Caveat stated in advance.** One isolate can be held by several collections under different codes, so
recurrence across labels can overstate independence. The labels behind every recurrence are listed so
that this can be checked.

## Outputs (`benchmarks/copy_checks/`)

- `trees/`: per order and marker, the unique-sequence FASTA and member map, MAFFT alignment, IQ-TREE
  files; `trees/monophyly.groups.tsv` (every group on every tree), `trees/monophyly.summary.tsv`,
  `trees/run.json` (tool versions, commands, hashes).
- `artefacts/`: `artefacts.copies.tsv` (all 903 copies), `artefacts.types.tsv`,
  `artefacts.summary.tsv`, `artefacts.thresholds.tsv` (L_R per region), `artefacts.json`.
