# What the short-5.8S paralog does to species counts and trees (exploratory)

For L. Tedersoo's follow-up to II.1 (13 October 2026): if short copies of different species group together, a
species carries more than one structural variant of the marker, which would raise species numbers counted from the
ITS and give conflicting signals between SSU, ITS and LSU trees. It fixes no rule and changes nothing locked.

Built by `stage3/copy_paralog_effects.py` from culture reference v4 and the Glomerales trees in `../trees`.

**Short copies.** The 43 short-5.8S divergent copies, plus PX214776.1 (*R. irregularis*, culture 4569_1185B):
its 5.8S is identical to a short-class 5.8S, but its culture has only two copies (5.8S of 152 and 156 nt), so the
length rule (3 nt or more below the culture's main-type median) cannot flag it, and it is classed main type in v4.
That makes 44 short copies, 41 of them in *Rhizophagus*.

## Species counts

Single-linkage clusters of all 850 retained copies of the 37 named species (within each order), and of their
main-type copies only (the short copy above left out):

| Region | Identity | Copies | Clusters | Per species | Species in 2+ clusters | Clusters joining species | Short-5.8S-only clusters | Same-length-only clusters |
|---|---|---|---|---|---|---|---|---|
| Full ITS | 98.5% | all copies | 240 | 6.49 | 34 | 2 | 18 | 2 |
| Full ITS | 98.5% | main-type copies | 220 | 5.95 | 33 | 2 | 0 | 0 |
| Full ITS | 97% | all copies | 110 | 2.97 | 24 | 2 | 12 | 2 |
| Full ITS | 97% | main-type copies | 96 | 2.59 | 21 | 2 | 0 | 0 |
| Krüger fragment | 98.5% | all copies | 110 | 2.97 | 24 | 1 | 13 | 2 |
| Krüger fragment | 98.5% | main-type copies | 95 | 2.57 | 21 | 1 | 0 | 0 |
| Krüger fragment | 97% | all copies | 56 | 1.51 | 14 | 4 | 8 | 2 |
| Krüger fragment | 97% | main-type copies | 46 | 1.24 | 9 | 4 | 0 | 0 |
| Whole copy | 98.5% | all copies | 71 | 1.92 | 16 | 3 | 11 | 2 |
| Whole copy | 98.5% | main-type copies | 58 | 1.57 | 12 | 3 | 0 | 0 |
| Whole copy | 97% | all copies | 36 | 0.97 | 8 | 7 | 5 | 1 |
| Whole copy | 97% | main-type copies | 30 | 0.81 | 3 | 7 | 0 | 0 |

Per species (`paralog.species.tsv`): on the full ITS at 98.5%, each species that carries the short paralog gains
1-7 clusters made only of its short copies; no short-only cluster joins two species.

## Trees

Where the short copies sit in the Glomerales tree of each marker, rooted on a *Funneliformis* tip
(`paralog.trees.tsv`; every short-only clade, its support, the main-type taxa of its sister group and the support
of the clade it forms with them in
`paralog.tree_clades.tsv`):

| Tree | Short copies | Identical to a main-type copy | In short-only clades of 2+ copies | Clades joining taxa (copies) |
|---|---|---|---|---|
| SSU flank | 44 | 21 | 12 | 0 (0) |
| Full ITS | 44 | 0 | 43 | 3 (37) |
| LSU flank | 44 | 1 | 42 | 2 (10) |
| Linked copy | 44 | 0 | 43 | 3 (37) |

- **SSU:** 21 of the 44 short copies are identical to a main-type copy of their species; the rest sit with their
  own species. The SSU carries no paralog signal.
- **Full ITS and linked copy:** 43 of 44 sit in clades of short copies only; three of them join two taxa:
  *R. irregularis* with *[Rhizoglomus] vesiculiferum* (27 copies; the pair is complex C4), *R. intraradices* with
  *R. prolifer* (8), and *R. irregularis* with *Rhizophagus* sp. (2).
- **LSU:** 42 of 44 sit in short-only clades, but most split by species (the *R. irregularis* and *[Rhizoglomus]
  vesiculiferum* short copies form separate clades). The *R. intraradices*/*R. prolifer* short clade (99/100) is
  sister to *Sclerocystis sinuosa* (100/100) on the LSU tree, but to the main-type copies of its own species on the
  ITS and linked trees: a supported conflict between markers.

  PYTHONPATH=<benchmark repo>/stage3 python3 stage3/copy_paralog_effects.py --ref reference/culture_ref_v4 \
      --trees benchmarks/copy_checks/trees --out benchmarks/copy_checks/paralog_effects
