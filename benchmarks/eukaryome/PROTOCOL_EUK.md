# EUKARYOME long-read test, protocol

Written 15 October 2026, after the predictions were locked (`predictions.lock`) and before any test read was
searched against MaarjAM, UNITE or V18, placed, scored or classified. The frozen M2b code and protocol v2 (N6)
are used unchanged, as in the LongGloDB test (`benchmarks/longglodb/PROTOCOL_LGDB.md` and its two amendments,
whose provisions are written into this document from the start). This document adds no rule and no threshold
to the method.

Suggested by L. Tedersoo (point 3 of his comments on version 15): the EUKARYOME long-read Glomeromycota as a
further external test, after removing copies identical to the culture references.

## 1. Test set

- **Source.** EUKARYOME v2.1 long-read set (`General_EUK_longread_v2.1.fasta`, 400,217 reads, SHA-256
  536e3e10…, CC BY 4.0; Tedersoo et al. 2024, *Database* baae043). Reads: the 5,546 Glomeromycota reads whose
  SSU, ITS and LSU windows are all at least 300 nt (the Glomeromycota queries of rdnalink Phase 0,
  `search/queries.tsv`), extracted with their headers to `euk21_glo_query.fasta` (SHA-256 aaab407b…).
- **Windows.** Stage 1 (amfsplit v0.3.1 with ITSxRust 0.3.0, `F.hmm` SHA-256 a8c38c5c…, HiFi preset, the
  options of the deposited copies: max_err_frac 0.15, only anchored, open ends, primers removed). 5,540 reads
  have windows: 2 had no ITSxRust anchors (FJ461846 and AJ133706, *Sclerocystis sinuosa*), 3 partial reads were
  not anchored, and 1 read is a concatemer (PQFF01000047). Stage 1 records SHA-256 efd27f2d…, full ITS
  87dccbda…. In all 5,541 anchored reads the SSU flank has the same length as rdnalink's SSU window, and in all
  but 4 (1–2 nt) the full ITS has the same length as its ITS window.
- **Names.** Each read's EUKARYOME labels are mapped to the names of culture reference v4 by `names.tsv`,
  built by `stage3/eukaryome_overlap.py names` from its suggestions and the seven decisions of
  `names_decisions.tsv`, each with its reason. Rules: *Rhizoglomus* is read as *Rhizophagus*; epithets match on
  their stem (Latin ending and trailing vowels removed: *irregulare* = *irregularis*, *bareae* = *bareai*);
  "cf." names count as genus only; EUKARYOME placeholder taxa (Glomeraceae.gen25, Archaeosporomycetes.ord01, …)
  are lineages without a named member and count as genera not in v4. EUKARYOME's labels are taken as they stand
  and are not changed in the light of any result.
- **Removed (29).** 19 reads whose full ITS is identical to a v4 copy (12 with known-species labels, 5 genus
  only, 2 *Nanoglomus*; 3 of the 19 are also identical in LROR–FLR2 and SSU flank; 4 are rDNA copies of genome
  assemblies of *R. irregularis* and *G. rosea*); 6 reads without a genus inside a named family or order;
  3 reads without an LROR–FLR2 window; 1 read of *Appendiculaspora leptoticha* (`names_decisions.tsv`).
- **Test.** 5,511 reads, with roles in the vocabulary of the LongGloDB test (`overlap.asvs.tsv`):

  | Role | Reads | Taxa | Strata |
  |---|---|---|---|
  | known species | 199 | 20 species | |
  | novel species, genus in v4 | 77 | 33 species | absent species in genera v4 holds |
  | novel lineage, genus not in v4 | 2,629 | 82 taxa | 2,017 genus-only reads of 26 named genera; 82 reads of 20 named species; 530 reads of 36 placeholder taxa |
  | genus only | 2,606 | 17 genera | 2,568 reads of 15 genera v4 holds; 38 reads of *Fuscutata* and *Appendiculaspora*, whose membership is open |

  EUKARYOME has no culture codes; each labelled taxon is one "culture" in the runner's sense (code `EUK`).
- **Flags, reported apart (as in the LongGloDB test).** 23 test reads carry INSDC, UNITE or assembly IDs rather
  than EUKARYOME's own, and may be records of the external databases themselves. 57 reads have a V18 record at
  least 99.5% identical in LSU (33 genus only, 20 known species, 3 novel lineage, 1 novel species); the LSU arm
  may see their own record under condition A. No read is identical in full ITS to an ASV of the public
  LongGloDB v1, so this test and the LongGloDB test share no sequence.

## 2. Predictions (locked before this protocol)

- The complex rule of `species_groups.py` applied to the test reads by `predict_longglodb.py`, unchanged, with
  the within-species spread T passed at the paper's definition (99th percentile, `statistics.quantiles`
  exclusive method, as in `its_rule.py`): ITS 5.98803%, LSU 3.02698%, SSU 0.80652% (full precision in the
  lock). This is the definition that LongGloDB amendment 1 identified as the paper's, so there is no second set.
  A read is "inside" a v4 species when its nearest copies of that species lie within T on all three markers.
- **Predicted:** of the 2,706 reads of taxa v4 does not hold (77 novel species, 2,629 novel lineage), only
  EUK1266127 (*Diversispora bareai*), inside *D. arenaria* and *D. jakucsiae* (complex C3), can be called a
  known species. Of the 199 known-species reads, 158 lie inside their own species (13 of them inside another
  species too) and 41 do not.
- **One rank up** (`eukaryome_overlap.py genus_inside`, same T): 50 *Nanoglomus* reads lie inside a v4 culture
  named only to genus, *Glomus* sp. 5089_2351A (42) or *Glomus* sp. 4406 (8), and two *Nanoglomus* reads are
  identical in full ITS to 5089_2351A (removed above). These cultures may be *Nanoglomus*, so a genus given to a
  *Nanoglomus* read through them is not counted as wrong: *Nanoglomus* has no genus outcome. Its false-known
  outcome is unchanged. No other genus has a read inside a v4 culture named only to genus.

## 3. Method (frozen)

As in `PROTOCOL_LGDB.md` section 3, run by `benchmarks/m4_external_holdout/euk_test.py`, which imports
`lgdb_test.py` unchanged and replaces only the reader of the test set, the plan of condition B (section 5) and
the culture IDs:

- Training: all 151 cultures of culture reference v4. Every test read belongs to a test culture that is always
  in the excluded set, so no test read trains merges, margins, unit links, D1/D2 comparators, nearest-copy
  pools or the N6 threshold. Before any test result, the 37 M4 species folds are rescored with the test reads
  present and excluded and must equal the stored M4 v2 results.
- Windows: the Stage 1 SSU flank, 5.8S, full ITS and LROR–FLR2 (section 1). EUKARYOME's own sequences outside
  these windows are not used.
- VT and SH evidence: the frozen `vt_hits` on BLAST against the MaarjAM and UNITE files of the v4 crosswalk
  (SHA-256 2b2c57e6… and abe98f43…), dropping hits to culture-reference accessions. The frozen VT and SH evidence
  of all 903 reference copies must be reproduced first.
- LSU: M2b species-clade units from EPA-ng placement on V18 (SHA-256 e6794cbb…), 0.95 credible mass, unchanged.
- Decision tree and N6 (protocol v2) unchanged. The copy class of a test read is unknown; every read is scored
  as a main-type copy and the D1 test applies as usual.
- Hits of a read to its own record in UNITE (UNITE and INSDC IDs among the reads) are not removed: an SH that
  no reference culture carries does not change a call.

## 4. VT coverage (fixed here, as LongGloDB amendment 2)

- Where the reads start (`eukaryome_overlap.py vtgeom`, from the Stage 1 windows and the v4 NS31–AM1
  copies): 5,083 reads start 21 nt after the 3′ end of the NS31 site and 70 more within 2 nt of that (19–23 nt),
  147 at 2 nt, 16 at 51 nt and 2 at 8 or 15 nt; 172 hold the whole NS31–AM1 fragment; in 50, mostly
  Paraglomerales, the site was not located.
- MaarjAM (3,066 references; the fragment start located in 2,933; median length 520 nt): a read starting 21 nt
  inside can cover at least 95% of 99.0% of these references at best, that is, when the alignment reaches the
  reference's 3′ end. Unlike the LongGloDB amplicons (about 60 nt missing, VT search nearly silent), the VT arm
  is therefore expected to work for most reads, but close to its coverage threshold. This count describes
  the database only; no read was searched against MaarjAM.
- **Primary analysis:** the frozen VT arm as it is.
- **Secondary analysis:** the VT arm with coverage measured over the part of each MaarjAM reference that the read
  spans (`vt_hits_spanned` of `lgdb_test.py`, unchanged); identity (97%), e-value, merged units, δ and every
  other rule unchanged.
- Reported: the VT state of the reads (matched, unmatched, silent) under each analysis and condition, the
  number whose VT state changes, B1 "boundary" calls and B2 flags.

## 5. Conditions

- **A. V18 intact.** Test windows added to the frozen full-V18 alignment and placed on the frozen tree and
  model; the 903 reference placements are the frozen ones.
- **B. Target removed from V18** (`overlap.v18_plan.tsv`), each round as an M4 fold (records out, alignment and
  model refitted, the 903 copies and all test windows re-placed):
  - 30 species rounds, one per absent species of a v4 genus whose label V18 holds (stem match on genus and
    epithet; *Dentiscutata erythropus* for EUKARYOME's *Fuscutata erythropus*); its reads are scored in its round.
  - 21 genus rounds, one per absent named genus that V18 holds: every V18 record of the genus out, compound
    labels included (*Entrophospora Alborhynchus walkeri*). All reads of the genus, named or not, are scored in
    its round.
- **Under A only, reported apart:** absent species whose label V18 does not hold (*Ambispora gerdemannii*,
  *Dominikia indica*, *Rhizophagus custos*; 5 reads); absent genera V18 does not hold (*Antiquispora*,
  *Durabilispora*, *Geosiphon*, *Palaeospora*, *Parvocarpum*, *Pseudoentrophospora*; 48 reads); and the 530
  reads of placeholder taxa. These are taxa that no curated reference names, the case the LSU-removal test
  approximates.

## 6. Outcomes

- **Primary:** false-known calls (a known species or species group) among the 2,706 reads of taxa v4 does not
  hold, under B (under A for the taxa of section 5 reported apart), with N6: by read, averaged over labelled
  taxa (species where named, otherwise the genus or placeholder taxon), and by stratum (absent species in v4
  genera; absent named genera; placeholder taxa).
- **Test of the complex claim:** every false-known call under A or B is checked against the locked predictions.
  A false-known call for any read other than EUK1266127 contradicts the claim and is reported as such.
- **Calls outside the predicted set, all reads:** every known-species call on any test read, the 2,606
  genus-only reads included, is checked against the predictions. For reads of absent taxa this is the test
  above; for the other reads, whose species may be the one called, it shows how far the rule describes the
  method's calls, and it does not depend on EUKARYOME's species labels. Reported by role.
- **Secondary:**
  - Known species kept under A (199 reads), with and without N6; split by the prediction (inside or outside
    their own species); the reads N6 demotes, with their ITS distance to the nearest v4 culture of their
    species, set beside the LongGloDB result.
  - Genus: for absent species in v4 genera, genus right, wrong or none (no genus outcome for *Dentiscutata
    erythropus*); for absent named genera and placeholder taxa, any genus given is wrong (no genus outcome for
    *Nanoglomus*, section 2); the share called "new lineage"; under A and B.
  - Descriptive: for genus-only reads of v4 genera, whether a genus given is the label's v4 genus
    (`names.tsv`, v4_genus); divergent-copy calls (D1, D2) and likely-artifact calls by role; reads whose
    5.8S carries the motif YTYGGRGARTATG (`short58_motif.py`), by label, and how they were called.
  - Every outcome with and without the flagged reads of section 1.

## 7. Limits and no retuning

- All rules, thresholds and pools stay fixed whatever the outcome; a change needs a new protocol version and a
  fresh evaluation. MaarjAM and UNITE stay intact, as in protocol v2.
- EUKARYOME's labels are its classification of environmental and reference reads, not identified cultures, and
  are noisy at species level: in the exploratory run (`benchmarks/eukaryome_longread`), the 99th percentile of
  the distance from a read to the nearest read with the same species label was 14.6% in full ITS, against 6.0%
  between cultures of one species in v4. A false-known call can therefore reflect a label. Label conflicts are
  listed for EUKARYOME's maintainers as disagreements between the call and the label, not as errors of either.
- EUKARYOME is curated by one of the authors (L.T.). Its labels are independent of the method, not of the
  author list.
- There are no genomes or cultures: copy classes are unknown and nothing is learned from the test reads.

## 8. Record

The SHA-256 of this file, `predictions.lock`, `names.tsv`, `names_decisions.tsv`, `overlap.asvs.tsv` and
`overlap.v18_plan.tsv` go into `protocol_hashes.txt`, which is committed and pushed, with the tables (read IDs
only, no sequences), before the run starts.
