# LongGloDB independent-culture test, protocol

Written 6 October 2026, after the predictions were locked (`bench/longglodb/predictions.lock`) and
before any test sequence was placed, scored or classified. Frozen M2b code and protocol v2 (N6) are
used unchanged; this document adds no rule and no threshold.

## 1. Test set

- Source: LongGloDB v1 (July 2026). Its 810 sequences from labs STFP/STFI (collection CCAMF) are the
  Stefani cultures reprocessed and are excluded. So are four cultures sharing a culture code with the
  Stefani set: BEG35 (same species) and three with conflicting names (MT106, ON205A, BEG9), which are
  reported apart as curation findings.
- Test: 280 ASVs in 52 cultures from labs KU, ZU and ZUII, with the roles fixed by
  `longglodb_overlap.py` (`overlap.asvs.tsv`, hashed in the lock):
  known species 22 cultures (164 ASVs), novel species with the genus in v4 17 (59), novel lineage with
  the genus absent from v4 9 (43), genus only 4 (14).
- Names: Rhizoglomus read as Rhizophagus; Funneliformis coronatum and Sclerocystis sinuos mapped to
  the v4 names. No other renaming.

## 2. Predictions (locked before this protocol)

The complex rule of `species_groups.py`, applied to the test ASVs with T recomputed from v4 by the
same method (ITS 5.98%, LSU 2.91%, SSU 0.75%; 99th percentile). An ASV is "inside" a v4 species when
its nearest copies lie within T of it on all three markers. Predicted: of 102 novel-species and
novel-lineage ASVs, only the 4 of *Gigaspora gigantea* KS302 (inside *G. rosea*) can be called a known
species.

## 3. Method (frozen)

- Training: all 151 cultures of culture reference v4. Every test culture is always in the excluded set,
  so no test sequence trains merges, margins, unit links, D1/D2 comparators, nearest-copy pools or the
  N6 threshold.
- Test sequences: Stage 1 windows (amfsplit v0.3.1, the options used for the deposited copies): SSU
  flank, 5.8S, full ITS, LROR-FLR2. LongGloDB's own SSU/ITS/LSU files are not used.
- VT and SH evidence: the frozen `vt_hits` on BLAST against the same MaarjAM and UNITE files as the v4
  crosswalk (SHA-256 checked against `crosswalk_v4.run.json`), dropping hits to culture-reference
  accessions, as for the reference copies.
- LSU: M2b species-clade units from EPA-ng placement, 0.95 credible mass, unchanged.
- Decision tree and N6 (protocol v2) unchanged. The copy class of a test ASV is unknown; every test ASV
  is scored as a main-type copy and the D1 test applies as usual.

## 4. Conditions

- **A. V18 intact.** Test windows added to the frozen full-V18 reference alignment and placed on the
  frozen tree and model; the 903 reference placements are the frozen ones.
- **B. Target removed from V18.** For each novel test species whose label V18 holds, one fold as in M4
  stages 01-05: its V18 records removed, alignment and model refitted, the 903 copies and the test ASVs
  re-placed. Novel species V18 does not hold by label are reported under A only, apart.

## 5. Outcomes

- **Primary:** false-known calls (known species or species group) among novel-species and
  novel-lineage ASVs under B, by ASV and averaged over species and over cultures.
- **Test of the complex claim:** every false known under A or B is checked against the locked
  predictions. A false known outside the predicted set contradicts the claim and is reported as such.
- **Secondary:** known-species retention under A for the 22 known cultures, split by whether an ASV is
  identical in ITS to a Stefani copy; genus calls (right, wrong, absent) and novelty rank against the
  expected rank; D1 calls among test ASVs; cultures whose LSU matches a V18 record at >= 99.5%, reported
  apart under A because V18 may hold their own record.

## 6. Limits and no retuning

All rules, thresholds and pools stay fixed whatever the outcome; a change needs a new protocol
version and a fresh evaluation. MaarjAM and UNITE stay intact, as in protocol v2. The test cultures
come from other collections but some may be the same isolates under other codes (flagged in
`overlap.cultures.tsv`); results are reported with and without them.
