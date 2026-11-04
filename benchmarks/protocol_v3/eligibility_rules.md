# Protocol v3: eligibility rules (extracted from PROTOCOL_v3.md)

**Confirmation set, by eligibility rule.**

1. Sources, in this order: Glomeromycota rRNA sequences deposited in INSD after the EUKARYOME v2.1 release of 9 September 2026; cultures added by a LongGloDB release after v1; reads added by a EUKARYOME release after v2.1; new culture sequencing.
2. Each sequence spans the SSU region, full ITS and LR0R–FLR2 with at least 300 nt each, as in the EUKARYOME test, and is a PacBio HiFi read or a long-read ASV.
3. Labels are culture names where they exist. A database classification forms its own stratum, because EUKARYOME labels put 14.6% (99th percentile) between reads of one species in full ITS, against 6.0% between cultures.
4. Excluded: sequences identical in full ITS to one in v4, LongGloDB v1 or EUKARYOME v2.1; sequences from a culture whose code matches a v4 or LongGloDB v1 culture. Records of MaarjAM, UNITE or V18 themselves are flagged and reported apart, as in the EUKARYOME test.
5. Strata: K, species v4 holds; A, absent species in v4 genera; L, absent genera and lineages; G, sequences labelled only to a genus.
6. Withholding: each K sequence is scored twice, with the reference intact (kept, demotions) and with its species withheld from the culture reference and from V18 in its own round (false known among absent close relatives). Genus rounds do the same at genus rank for R3; MaarjAM and UNITE stay intact, as in every earlier test.

For K, the independence is new sequences and new collections, not new species; A and L add new species and genera where the sources have them.

**Triage sets.** Positives: the 2,629 EUKARYOME absent-lineage reads, already measured as development (section 7). Negatives: non-Glomeromycota long reads spanning the same three regions, taken from EUKARYOME v2.1 by a taxon list fixed in section 10 and never run in this project. They measure how many non-AMF sequences reach naming under R6, and how many are then called new lineage.
