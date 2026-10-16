# Response to the referee report (MEE format, major revision): working draft

Draft of 15 October 2026, for the authors. It maps each point to the change in `paper/main.tex` (v18 draft), and
it can become the response letter. **Pending** marks work that needs the EUKARYOME run or files from the
workstation. Section numbers refer to the v18 draft.

## Summary of changes

1. The complex claim is split into a retrospective statement (errors fall in complexes defined from all cultures)
   and a prospective test. For the prospective test, the rule is applied within each round of the LSU-removal test,
   with the spread recomputed without the withheld cultures (new Table 10).
2. The rate for an absent species now leads with the LSU-removal test (7.9%, 95% interval over species 1.4–17.3%).
   "Unknown" is defined by the reference condition (Introduction).
3. EUKARYOME results: **pending** (the run is in progress; placeholder subsection 3.9).
4. Intervals resample species and cultures (new Table 12). The ITS rule is presented as a chosen threshold with a
   trade-off.
5. The paralog evidence is split by class (new Table 5). The title of Section 3.3, the Abstract and the Conclusions
   now claim a paralog family only for *Rhizophagus*.
6. The first decision step is renamed an SSU-discordance flag (text and Figure 1), with new in silico controls.

## Point by point

### 1.1 The species-complex analysis is not independent

Agreed. We now separate the two statements the referee distinguishes.

- Methods 2.5 (new paragraph) gives the rule's two readings. Applied to reference species, it gives the complexes
  (retrospective). Applied to a query, it gives "inside" a species (prospective, needs no test result), and this
  reading fixed the LongGloDB and EUKARYOME predictions before those runs.
- Results 3.7 (retitled) labels the original finding retrospective. It adds the analysis the referee requested:
  in each LSU-removal round we removed the withheld cultures, recomputed the spread from the remaining cultures
  (ITS 5.98–6.08%, LSU 2.31–3.11%, SSU 0.74–0.81% across rounds), and tested each withheld copy against every
  remaining species (`paper/analysis/foldwise_complexes.py`; Table 10; Methods S6).
- **Species level, the rule anticipates the errors.** Each of the 6 test species with false-known calls had copies
  inside a congener. None of the 27 species without such copies had a false-known call. One species predicted at
  risk (*[Rhizoglomus] vesiculiferum*) had none.
- **Copy level, it falls short.** Only 13 of 56 *Rhizophagus irregularis* copies lie inside *[R.] vesiculiferum*,
  yet 18 were called that species. The others lie outside it in LROR–FLR2 (3.5–5.4% against a recomputed spread of
  2.95%). So at least 5 of the 53 errors fall outside the prediction, and at most 48 inside it. All of them are in
  the one pair that is not a complex at the 95th and 97.5th percentiles.
- Discussion 4.4 is retitled "Where the markers set a floor" and states the observed limit and the prospective
  result separately.
- **Pending:** exact copy-level counts. These need `amf_m4_step1/runs/aggregate_v2/species.{v1,v2,control_v1,control_v2}.copies.tsv`;
  then `foldwise_complexes.py --m4 <dir>` writes the counts and the `\todo` in 3.7 is removed.

### 1.2 The headline rate depends on what was removed

Agreed on the framing.

- The Abstract now leads with the LSU-removal rate and its interval, with the intact-reference rate beside it.
- The Introduction defines an unknown species as one absent from the reference against which a call is made, and
  states that no test measures the error for a species never sequenced.
- Methods 2.6 states that the LSU-removal test is the closest to such a species, but not equal to it.
- Results 3.6 says this rate, not the 2.2% of the databases-intact test, is the better guide for an absent species
  with close relatives.
- Discussion 4.5 says each rate belongs to its reference condition.
- Not yet done: removal from UNITE as well, planned before journal submission. Removal from MaarjAM is ill-defined,
  because VTs are not species. The existing bound stays: the VT accounted for the novelty of 18 of 673 copies
  (2.7%), so for a species never sequenced the rate may be higher by up to that much.

### 1.3 EUKARYOME results missing

Agreed; the reviewer saw the draft while the pre-registered run was still going.

- Placeholder Section 3.9 lists what will be reported:
  - false-known calls for absent species and lineages, averaged over taxa as well as reads;
  - calls outside the locked prediction;
  - the genus given to reads of absent genera;
  - the effect of LSU removal and of the second VT analysis;
  - retention and ITS-rule demotions, as a third collection;
  - sensitivity to uncertain EUKARYOME labels.
- **Pending** the run.

### 1.4 The paralog interpretation is too broad

Agreed.

- Section 3.3 is retitled "Divergent copies, and a paralog family in *Rhizophagus*"; "They are not artifacts" now
  reads "Most are not artifacts".
- New paragraph and Table 5 (`divergent_classes.py`) set the evidence side by side for four classes, from the
  existing copy checks:
  - **Short-5.8S *Rhizophagus*** (40 copies, 22 cultures). All share the deletion with other species. 22 recur
    in whole ITS under another label. 10 of the 10 from cultures with replicate libraries recur in another library
    of the same culture. 26 match a public record at ≥99%.
  - ***Septoglomus viscosum* short copies** (3 copies, 2 cultures of one species). They share their own deletion,
    2 of 2 recur in a second library and match public records, and their SSU is degraded. Reading: genomic, origin
    beyond this species unknown.
  - ***Gigaspora rosea* same-length copies.** Co-dominant types of one genome.
  - **Other same-length copies** (5). Unresolved.
- Replicate amplifications exist only for cultures with 2–4 libraries, which covers 16 of the 51 divergent copies;
  stated in Results and Limitations. No independent genomic reads were used.
- Discussion 4.2 is retitled and adds that predicted folding shows a structural change, not loss of function, and
  that "pseudogene" would go beyond the evidence.
- Abstract and Conclusions now claim the paralog family for *Rhizophagus* only.

### 1.5 The distance rule is not portable

Agreed.

- New Table 12 (`intervals.py`) gives cluster-bootstrap intervals that resample species and cultures, for false
  known and for retention, with the number of species and cultures with events. The withholding tests were
  rescored from the frozen tables with the frozen ITS rule, and every count was checked against the paper.
- The ITS rule's demotions: 0.8% in the culture test (0.0–2.1% by culture), against 9.1% in LongGloDB (2.0–21.4%),
  where it prevented no error (Results 3.8).
- Discussion 4.7 now calls the rule an empirically chosen decision threshold with a trade-off, not a calibrated
  error rate.
- **Pending:**
  - culture-level intervals for the LSU-removal test (same files as 1.1);
  - EUKARYOME as a third collection for the transfer question. Its labels are a classification, not culture names,
    so a leave-one-collection-out calibration across three collections needs care.

### 1.6 The artifact detector has no positive control

Agreed on both counts.

- The step is renamed "SSU discordance (possible artifact)" in the text, in Methods S5 and in Figure 1.
- In silico controls (`stage3/artifact_control.py`; Methods S8; Results 3.4):
  - none of 852 real copies is flagged, although 26 lie ≥1% from every other culture in the SSU;
  - copies with even errors: 0.7% flagged at 0.5% errors, 29% at 2%, 92% at 8%;
  - chimeras that join one species' SSU to another's ITS: 0 of 300 flagged;
  - chimeras with the breakpoint inside the SSU: 10 of 300 flagged.
- Discussion 4.3 now says the flag catches noisy reads, is not a chimera screen, and could be tripped by a genuine
  lineage with no cultured relative.

### Section 2: claims to narrow

| Claim | Change | Where |
|---|---|---|
| SH counts inflate richness | "can be inflated by variation within genomes"; field magnitude not tested | Discussion 4.1 |
| Floor for novelty | Observed limit vs prospective performance, kept separate | Results 3.7, Discussion 4.4, Conclusions |
| Divergent copies are paralogs | By class; *Rhizophagus* only for "paralog family" | Results 3.3, Table 5, Discussion 4.2, Abstract, Conclusions |
| Protein-coding genes don't separate complexes | Limited to pairs with records; other genes or criteria may separate them | Discussion 4.6, Conclusions ("available protein-coding sequences") |
| "Unknown" | Defined by reference condition; each rate reported with its condition | Introduction, Methods 2.6, Discussion 4.5, 4.9 |
| 13-nt motif | Applies to full-ITS and linked reads, and to ITS2 before ITSx trimming; motif written from the copies it was tested on | Discussion 4.2 |

### Section 3: additional analyses

- **A. Fold-wise complexes.** Done at species level (Table 10); copy level **pending** the four aggregate files.
- **B. All-systems removal.** Framing done (1.2). UNITE removal is planned for journal submission. MaarjAM removal
  is ill-defined, and the limitation is stated.
- **C. Intervals and counts of species and cultures.** Done (Table 12). LSU-removal culture level and EUKARYOME
  are **pending**.
- **D. V18 label sensitivity** (51 of 220 multi-record labels not exclusive). Not done. It needs placement re-runs
  on the workstation. Planned for journal submission.
- **E. Divergent-copy detection by class.**
  - Short-5.8S copies were recognized only in *Rhizophagus*: 36 of 40 with the culture withheld, 29 with the
    species withheld. None of the 3 *S. viscosum* copies was recognized.
  - In the genus test no divergent copy was recognized, because both tests need a related culture.
  - Only the *Rhizophagus* motif works without a related culture, and it was not one of the fixed steps.
  - Where: Table 5, Results 3.5 and 3.6.

### Section 4: additional issues

- **Culture identities as ground truth.** Limitations now name the cases that touch the results:
  - one of the six known-species copies demoted by the ITS rule comes from *M. litorea* 5274;
  - the reference's *A. leptoticha* (ON205A) is *A. gerdemannii* in INVAM;
  - agreement with a deposited name is not proof of species.

  A rescoring without disputed labels is not done yet.
- **Genus output (*Sclerocystis*).**
  - Results 3.6 now says these copies get a definite genus rather than a flag, which is why the failure is easy to
    miss.
  - New paragraph in Discussion 4.5 proposes a flag for a genus that rests only on the nearest represented clade,
    for testing; it is not part of the fixed steps.
  - Implications says to treat such genus calls as provisional.
- **Tool claims.** The Introduction now states the evidence is a crosswalk and a novelty benchmark, not a
  validation of community profiles or richness.

### Section 5: requested revisions

1. EUKARYOME: **pending** the run (Section 3.9).
2. Complexes without held-out information: done at species level; copy level **pending**.
3. Rates framed by reference condition, leading with LSU removal: done.
4. Uncertainty and transfer: intervals done; transfer across three collections **pending** EUKARYOME.
5. Narrower biology with the function caveat: done.
6. Ecological implications revised: done (Discussion 4.1, 4.9; Introduction).

## For the authors

- **Title.** Kept, with "unknown" now defined. A narrower alternative, if wanted: "Linked long-read rDNA copies show
  where AMF naming systems disagree and when novelty can be called".
- **Abstract length.** About 440 words; MEE allows 350 (v17 was already about 370). Trim before journal
  submission.
- **Exploratory analyses.** All analyses added after review are exploratory and change no rule or prediction.
  They are labelled "added after review" in the Methods.
