# AMF M2b: LSU placement experiment

Completed 29 September 2026. Frozen culture reference: v4, fingerprint `181938d695f3df40`.

## Result

LSU phylogenetic placement is implemented and has been run on all 903 LROR–FLR2 copies. The method that groups compatible placements within pure reference species clades is a promising experimental replacement for the current LSU matching arm. It improves species resolution and reduces false-known calls in this cohort, with measurable remaining errors and regressions.

The comparison below uses the **training-only Phase 4 baseline**, including physical-culture exclusions and the earlier D1 consistency changes. It therefore differs from the original global v4 logs with 200/803 false-known calls. This run reproduced the frozen training-only baseline exactly before evaluating the new LSU methods. The original v4 reference and previous results were not modified.

| Measure | Frozen baseline | Compatible tree edges | Compatible species clades |
|---|---:|---:|---:|
| Represented-species copies retaining a compatible final known call | 748/755 (99.1%) | 746/755 (98.8%) | **751/755 (99.5%)** |
| Final exact-species calls among those 755 copies | 603 (79.9%) | 696 (92.2%) | **701 (92.8%)** |
| Final containing-species-group calls | 145 | 50 | 50 |
| Wrong final known names among represented-species copies | 0 | 0 | 0 |
| False-known calls, seven single-culture species | 3/48 (6.3%) | 0/48 | **0/48** |
| False-known calls, all named species held out | 215/803 (26.8%) | 71/803 (8.8%) | **71/803 (8.8%)** |
| False-known rate, each species weighted equally | 23.5% | 8.5% | 8.5% |
| False-known rate, each culture weighted equally | 24.7% | 7.3% | 7.3% |

“False known” includes both a known-species call and a known-species-group call when the true species has no training culture. Copies within a culture are not independent biological replicates; the macro averages help show the effect of unequal copy counts. These are descriptive results, without a claim of independent validation.

Before the Phase 4 novelty rules, the species-clade method reconciles 704/755 copies to their exact species and 51/755 to a group containing it: all 755 are compatible. After Phase 4, three become genus-level novelty candidates and one is flagged divergent. Thus the reconciliation negative control passes on this run, but the complete Phase 4 classifier still has known-species demotions.

## What changed, and what did not

The external V18 backbone contains 495 reference sequences: 460 AMF sequences and 35 outgroups. Its topology was retained while branch lengths and the GTR+G4 model were fitted to a reference-only MAFFT alignment. All 903 culture queries were then aligned into those fixed columns and placed with EPA-ng.

Each query contributes a compatible set of tree edges containing at least 95% of its likelihood weight, including every tie at the cutoff. Two LSU interpretations were specified before inspecting the placement outcomes:

1. Keep each compatible edge as a distinct unit.
2. Collapse compatible edges inside maximal clades bearing a single V18 species label. Different components of a non-monophyletic label remain separate; mixed-species edges remain individual units.

Only training cultures supply the links from those units to culture species names. A training culture uses the intersection of its eligible main-copy sets if nonempty, otherwise their union, matching the established bridge rule. No held-out culture contributes to those links. VT/SH training, D1/D2, and the Phase 4 decision rules were kept fixed. The experiment changes the LSU units, not the entire classifier.

There are 148 physical-culture folds and 50 species/unnamed-culture folds. A physical-culture exclusion also removes other lineages carrying the same culture label. A named-species exclusion removes every culture of that species and its co-cultured lineages. The external V18 taxa remain in the tree: this is **culture-bridge exclusion, not full external-species removal**.

## Gains and regressions

All six earlier known-species demotions are recovered by the species-clade method:

| Copies | Culture species | Final culture-holdout result |
|---|---|---|
| PX207141.1 | Gigaspora rosea | Known species |
| PX214440.1, PX214441.1 | Funneliformis caledonium | Known species |
| PX214953.1 | Diversispora epigaea | Known species |
| PX214979.1, PX214980.1 | Diversispora varaderana | Known species |

The exact-edge method recovers four of the six; the two F. caledonium copies need the clade grouping.

The clade method introduces three new known-species demotions, all because the LSU unit has no training-culture holder under the existing culture-set rule:

| Copy | Culture | New final status |
|---|---|---|
| PX214751.1 | Dominikia duoreactiva, 4392_1156C | Novel species, genus known |
| PX214760.1 | Dominikia bonfanteae, C10660_2264B | Novel species, genus known |
| PX214942.1 | Diversispora epigaea, BEG47_1140C | Novel species, genus known |

The earlier PX214776.1 main/divergent annotation conflict remains a divergent-class call. The reference labels were not edited to make the score improve. The exact-edge method introduces six known-species demotions, listed in `comparison/audit/known_species_regressions.tsv`.

In species holdout, both placement methods resolve 150 of the baseline's 215 false-known calls, retain 65, and introduce **six new false-known calls**, leaving 71. Those new errors are PX214652.1, PX214658.1, PX214661.1 and PX214663.1 (G. chinense called G. rugosae), plus PX214979.1 and PX214980.1 (D. varaderana called D. aestuarii). The latter two improve when another culture of their species is available, but become falsely named when their entire species is withheld.

All 71 remaining false-known calls involve six species in 17 cultures:

| Held-out species | Incorrect known call | Copies |
|---|---|---:|
| Glomus chinense | Glomus rugosae | 28 |
| Glomus rugosae | Glomus chinense | 4 |
| Ambispora callosa | Ambispora leptoticha | 20 |
| Ambispora leptoticha | Ambispora callosa | 10 |
| Diversispora varaderana | Diversispora aestuarii | 8 |
| Diversispora aestuarii | Diversispora varaderana | 1 |

Five false-known copies have best-edge likelihood weight ratios of at least 0.95. High placement support alone is therefore insufficient to establish a correct species name. No new LWR or pendant-length cutoff was fitted to these errors.

## The resolution cost

There is also a loss of genus resolution among held-out species:

| Final status among 803 named main copies, species holdout | Baseline | Species-clade LSU |
|---|---:|---:|
| Known species or species group: false-known calls | 215 | 71 |
| Novel species, genus known | 368 | 366 |
| Novel lineage | 219 | 365 |
| Divergent class | 1 | 1 |

In the paired comparison, **146 copies change from “novel species, genus known” to “novel lineage.”** This is loss of a genus call from the culture bridge; it does not prove that the sequences represent new genera or deeper lineages. External clade annotations are provided separately and are not substituted for validated culture predictions.

Short-5.8S recognition is unchanged: 36/43 in culture holdout and 29/43 in species holdout. Two main-class copies are flagged divergent in each mode, as in the frozen baseline. No copy is flagged as an artifact. LSU placement does not, by itself, settle the same-length divergent-copy interpretation or validate environmental novelty.

## Verification and reproducibility

- The pinned V18 repository commit is `f4a0014336f49be6aea48e8374dae052360a88bd`; the AMF-only FASTA digest matches the previously captured evidence. All 495 tree tips match the full reference FASTA, with no duplicate IDs.
- The reference alignment has 1,844 columns. Reference residues are preserved, optimized-tree bipartitions match the supplied topology, and query addition leaves every reference column unchanged.
- All 903 queries placed successfully, each evaluated on all 987 tree edges. Saved LWR sums range from 0.999999999997 to 1.000000000002. Exhaustive output is retained at 12-decimal precision in `placement/epa/epa_result.jplace`.
- All queries pass the prespecified 90% retention requirement. Minimum retention is 99.46%; median retention is 100%; 241 query insertion bases are omitted in total. Three independently aligned check queries match their full-batch alignments exactly.
- The median compatible set contains two edges; the largest contains seven. The median best-edge LWR is 0.886. Neither statistic is a calibrated species-accuracy estimate.
- Both holdout baselines reproduce the previous status, placement, rank, species/genus calls, correctness, reason and 5.8S fields exactly. Training holders and D1/nearest-SSU pools exclude held-out cultures.
- Five targeted unit tests pass: cutoff ties, truncated-weight rejection, separate components for repeated species labels, nested-clade collapse, and exclusion of query-only units from training support.
- Source code, frozen inputs, fitted model, complete placement output, per-copy comparisons, audit tables, commands, software versions and checksums are included. The redundant derived `placement_weights.json.gz` is omitted from the archive; it is regenerated by `summarize` from the included full jplace file.

Software used: MAFFT 7.526, RAxML-NG 1.2.2, EPA-ng 0.3.8, Python 3.12.14, Biopython 1.88 and edlib 1.3.9.post1. See `SOURCES.md` for primary documentation and `software_provenance.json` for versions, source commits and binary hashes. Third-party executables are not bundled.

## Interpretation and next stage

Use the species-clade implementation as the **leading experimental M2b option**, keeping the baseline available. It gives better known-species sensitivity and resolution than exact edges on this cohort while achieving the same false-known count. This recommendation is based on development results, not a new independent test.

The next validation should freeze this configuration and remove the target species from the external reference as well as from the culture bridge. It should also assess genus coverage and divergent-copy behavior with independent material. The three new known-species demotions and the three unresolved sister-species pairs are the concrete cases to investigate before changing production defaults. Any resulting rule changes need a fresh evaluation, not retuning against these same scores.

The reference itself limits interpretation: 51 of its 222 multi-record species labels do not form exclusive clades on the supplied backbone. That is a label/topology audit, not proof that those species are biologically invalid. No new backbone bootstrap analysis was run, and the original bootstrap comments are not represented as support for the newly fitted alignment/model. LWRs describe conditional relative attachment support, not probabilities of correct species identification. This completed run covers LROR–FLR2; D2 placement, environmental library replication, and the outstanding Stage 1 flags remain separate tasks.

To inspect or reproduce this run, start with `README.md`. Re-running the comparison does not require rerunning the alignment or downloading the large external databases.
