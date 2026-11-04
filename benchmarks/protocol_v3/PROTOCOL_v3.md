# Protocol v3: locked test of the amfsplit rule changes

Oct 19, 2026 · @Bioinformatics

## Status and purpose

This protocol decides which of nine rule changes become amfsplit defaults, from one locked run on a fresh labelled set. It was locked on 4 November 2026 (section 9, step 6), and any change after that is an amendment with its reason.

The changes come from the paper, the design brief of 16 October 2026 and the audits of 19 to 22 October: unnamed cultures as taxa (R1), the inside gate (R2), a provisional-genus flag (R3), limits fitted on two collections (R4), the *Rhizoglomus* motif step (R5), triage that labels instead of filtering (R6), VT coverage over the spanned part (V), a genus-clade fallback for calls without a genus (R7), and a flag on novelty that rests on one marker (R8). Results already seen suggested each one, so the culture folds, LongGloDB v1 and EUKARYOME v2.1 give development estimates only.

It also locks one measurement on development data, E: each test species removed from all three external databases, with the paper's own bound as the prediction (section 5).

A change that meets its adoption rule (section 8) becomes a default; otherwise it stays an opt-in flag. The protocol makes no claim about undescribed species, community profiles or richness, or about inputs other than PacBio HiFi reads and ASVs; ONT consensus input needs its own error-series test.

## Prerequisites

The locked run starts only when these four steps are done, in this order.

1. **Gate 0 passed (0a, 0b and Gate 3 done).** amfsplit v0.1 reproduces release 0.3.1 assignment for assignment, in two checks. Gate 0a, the nearest-reference arm in the 148 culture and 50 species rounds, passed byte for byte on 20 October with all 34 workspace tests. Gate 0b, every other per-call table of release 0.3.1, passed byte for byte by 1 November. Gate 3 then reproduced the paper's calls from v0.1's own BLAST evidence and EPA-ng placements on 2 November, and v0.1.0-rc.1 is tagged.
   - Exact, from the stored evidence and placements. Culture and species rounds in all three LSU arms: nearest reference 748/755 kept, 603 exact, 215/803 false known; species clade 751, 701, 71; species clade with N6 745, 697, 18; single-culture species 3, 0, 0 of 48. Then the 37 LSU-removal rounds (53/673 with N6), the 15 genus rounds, LongGloDB (4/102 false known, 136/164 kept) and EUKARYOME (13/2,706, 152/199), both VT analyses, conditions A and B.
   - Rebuild, done on 20 October. Placement refitted from the V18 source at f4a0014 matched every frozen artifact and all 903 rows with the pinned RAxML-NG 1.2.2 binary (SHA-256 b5038949…); an unpinned build moved PX207259.1 to another edge. MAFFT 7.526 and EPA-ng 0.3.8 matched by output, and every amfsplit db manifest records their executable hashes.
   - Parity details. Distances equal edlib's infix (HW) distance with the shorter sequence as query and the first argument on equal lengths, because HW distance is not symmetric for equal-length pairs. Percentiles follow `statistics.quantiles` (exclusive method). No extra equalities, so N and other IUPAC codes count as mismatches. Evidence floats are parsed with round-trip precision. N6 and the complex rule pool different copies for their limits (section 7), and each keeps its own pool.
2. **Arm code written and development estimates computed** (section 7). The arm code computes R1 to R8, V and S from v0.1's output and its refits; it is checked against the development estimates, and both are hashed into the lock.
3. **Name checks for R4 written as curation edits**, each with its check, before any limit is fitted. Done on 3 November 2026: benchmarks/protocol\_v3/r4\_name\_checks.tsv leaves out the four code-sharing cultures (three swaps, one duplicate) and MR101 and MX982A, whose sequences confirm neither label; MR101 lies inside F. mosseae on every marker, which the sensitivity arm uses.
4. **Fresh set assembled by the rule in section 5, at the set lock (section 9, step 7)**: identifiers, labels and stratum counts only. No fresh sequence is cut, searched, placed, triaged or scored before the lock, not even to smoke-test v0.1 or the arm code.

## Rules under test

Each rule is fixed here exactly as it will run; everything else stays the frozen method (rules v2: steps A1 to X1 with N6). Every quantity a rule learns from cultures is fitted on the training cultures of each round only.

**R1, unnamed cultures as taxa (fix 1).** Each culture not identified to species is a taxon of its own, labelled by genus and culture code. Two variants run as separate arms.

- R1a, inside. After the decision steps, a sequence inside at least one unnamed culture on all three markers, and inside no named species, is called "unnamed culture" (code U1): the cultures' genus and codes, no species. Divergent, artifact and discordant calls keep their call and gain flag U.
- R1b, block. A matched VT, SH or LSU name held only by unnamed training cultures blocks a species: a known call becomes genus only (as N3), flagged with the cultures' codes. This is the paper's alternative.
- A sequence inside both a named species and an unnamed culture keeps its named call, written as a group that includes the unnamed culture.

**R2, inside gate (fix 2).** After N6, and after R1a when both run, a known species or species group stands only if the sequence lies inside the called species, or one member of the called group, on every marker it holds. Inside means within the marker's within-species limit of the nearest main-type sequence of that species, in full ITS, LR0R–FLR2 and the SSU region. Limits are the complex rule's (99th percentile, exclusive method, unique sequence–culture pairs; ITS 5.98803%, LSU 3.02698%, SSU 0.80652% with all cultures), refitted per round; a failing call becomes "novel species, genus known" (code N7) with the failing markers recorded.

**R3, provisional genus (fix 3).** A flag (G1) that never changes a call, on any call that gives a genus. Two statistics run as separate arms.

- R3a, distance. Within-genus limit per marker: 99th percentile of the distance from each main-type sequence of a named species to the nearest main-type sequence of another named species of its genus, over genera with two or more named training species. Flag when the sequence lies beyond that limit from every training culture of the given genus (on any marker; on all markers as a sensitivity).
- R3b, placement. Flag when the pendant length of the best LSU placement (the smallest among the edges tied for the top weight, because EPA-ng lists tied placements in run-dependent order) exceeds the 99th percentile of that length among culture main-type sequences in the 37 LSU-removal rounds, where the species is absent from V18 and its genus present.

**R4, limits on two collections (fix 4).** The within-species limits of N6 and of the complex rule, which R2 uses, are refitted, each with its own pool, on culture reference v4 plus the LongGloDB v1 cultures of other laboratories, with nearest other-culture distances allowed to cross collections. LongGloDB ASVs get copy classes by v4's 5.8S rule. Left out of the fit: the four cultures sharing a code with a v4 label, and MX982A and MR101 unless the name checks clear them; the percentile (95th, 97.5th or 99th) is chosen before the lock by section 8, from cross-fitted estimates (fit on one collection, score the other).

**R5, *Rhizoglomus* motif step (fix 5).** A new step D4, after D1 and D2 and before B1. A sequence is a divergent sequence, named by its VT alone, when its 5.8S as ITSxRust delimits it contains YTYGGRGARTATG exactly and every matched VT name is linked only to *Rhizoglomus* training cultures; a motif without that SSU evidence is flag M1 and changes nothing. On the culture reference PX214776.1 (*R. irregulare* 4569\_1185B, classed main type, carrying the deletion) is already divergent through D1 in every design, and R5 leaves it, the 29 *Oehlia diaphana* sequences and PX214775.1 (*Glomus* sp. 5089\_2351A) unchanged.

**R6, triage labels instead of filtering.** Triage (`triage.py`: LR0R–FLR2 against V18 with its four outgroups, vsearch, `--min-id 90`) labels each sequence AMF, unresolved, non\_AMF or no\_hit. AMF and unresolved sequences go on to naming, and the other two are reported but not named. A new lineage (N4, N5) on an unresolved sequence is reported as "new lineage, AMF membership unconfirmed".

**V, VT coverage over the spanned part.** The VT arm counts coverage over the part of each MaarjAM reference the sequence spans (`vt_hits_spanned`, as in LongGloDB amendment 2); identity, δ and everything else unchanged. It runs only with R3 on, because its genus error in the LongGloDB test (*Racocetra* given *Cetraspora* through VTX00041, LSU name unmatched) is the case R3 should flag.

**R7, genus-clade fallback.** A call without a genus (new lineage, N4 or N5) takes the genus of the V18 genus clade its LSU places in, when that clade's unit is held by two or more training cultures of one genus that agrees with the clade's label (code G2). It never gives a species, and its two-culture minimum was settled in section 10. It runs only with R3 on, because its known failure is a sister genus standing in for an absent one (section 7).

**R8, novelty on one marker.** A new lineage resting on one unmatched arm, with the other two no call, keeps its call and carries flag N8. The variant R8d demotes it to "insufficient evidence". R8d is the rule audited as N6\_MULTI\_UNMATCHED\_v0.1, renamed because N6 is the paper's ITS rule.

## Arms

Fourteen arms answer the questions the rules raise, all from one execution: one set of fragments, evidence and placements per condition, with arms differing only in post-processing and fitted limits. Every arm is reported, whatever its result.

| Arm | Rules on | Question | Read for |
| --- | --- | --- | --- |
| F | none (v0.1, rules v2) | Baseline; the complex claim on fresh data | Every outcome |
| F+R1a | R1a | Does labelling unnamed-culture matches remove false species? | False known, U1 calls, kept |
| F+R1b | R1b | The paper's block alternative | False known, kept |
| F+R2 | R2 | The gate alone | Kept, false known |
| F+R1a+R2 | R1a, R2 | The gate with unnamed cultures, as the brief proposes | Kept, false known, U1 calls |
| F+R4 | R4 | Two-collection limits alone (N6 and complexes) | Kept, false known |
| F+R2+R4 | R2, R4 | The gate under two-collection limits | Kept, false known |
| F+R5 | R5 | The motif step | Paralogs recognized, wrong reclassifications |
| F+R3a | R3a | Distance-based genus flag | Flags on wrong and right genera |
| F+R3b | R3b | Placement-based genus flag | Flags on wrong and right genera |
| F+V+R3 | V, and the R3 variant chosen | VT coverage over the spanned part | Genus errors, flags, VT states |
| F+R7 | R7, with both R3 flags computed | Does the fallback restore genera without adding unflagged wrong ones? | Genus right, wrong and flagged |
| F+R8d | R8d | Does demoting single-marker novelty spare known species a false "new"? | New lineage and insufficient evidence, by stratum |
| S | The rules adopted under section 8, in S's order: R5, R4, R1, R2, R6, R7, R8, R3 | The candidate default | Every outcome |

R6 changes which sequences reach naming, not how they are named, so it is read on the triage sets of section 5 rather than as an arm. R8's flag changes no call either, so its counts come from F.

## Data

Everything already seen is development data; the confirmation set is defined by a rule and counted, but not run, at the set lock.

**Development data (estimates only).** Culture reference v4 (fingerprint 181938d695f3df40) in its 148 culture, 50 species, 37 LSU-removal and 15 genus rounds; LongGloDB v1 (280 test ASVs of 52 cultures, local public files, never redistributed); EUKARYOME v2.1 (5,511 reads). Because R4 fits on LongGloDB, R4's estimates are cross-fitted: fit on v4 and score LongGloDB, then the reverse.

**Confirmation set, by eligibility rule.**

1. Sources, in this order: Glomeromycota rRNA sequences deposited in INSD after the EUKARYOME v2.1 release of 9 September 2026; cultures added by a LongGloDB release after v1; reads added by a EUKARYOME release after v2.1; new culture sequencing.
2. Each sequence spans the SSU region, full ITS and LR0R–FLR2 with at least 300 nt each, as in the EUKARYOME test, and is a PacBio HiFi read or a long-read ASV.
3. Labels are culture names where they exist. A database classification forms its own stratum, because EUKARYOME labels put 14.6% (99th percentile) between reads of one species in full ITS, against 6.0% between cultures.
4. Excluded: sequences identical in full ITS to one in v4, LongGloDB v1 or EUKARYOME v2.1; sequences from a culture whose code matches a v4 or LongGloDB v1 culture. Records of MaarjAM, UNITE or V18 themselves are flagged and reported apart, as in the EUKARYOME test.
5. Strata: K, species v4 holds; A, absent species in v4 genera; L, absent genera and lineages; G, sequences labelled only to a genus.
6. Withholding: each K sequence is scored twice, with the reference intact (kept, demotions) and with its species withheld from the culture reference and from V18 in its own round (false known among absent close relatives). Genus rounds do the same at genus rank for R3; MaarjAM and UNITE stay intact, as in every earlier test.

For K, the independence is new sequences and new collections, not new species; A and L add new species and genera where the sources have them.

**Triage sets.** Positives: the 2,629 EUKARYOME absent-lineage reads, already measured as development (section 7). Negatives: non-Glomeromycota long reads spanning the same three regions, taken from EUKARYOME v2.1 by a taxon list fixed in section 10 and never run in this project. They measure how many non-AMF sequences reach naming under R6, and how many are then called new lineage.

**Measurement E, removal from all three databases.** E reruns the 37 species rounds of the paper's M4 harness rather than a parallel one. In each round V18 loses the species' records as in M4. MaarjAM and UNITE lose every record whose nearest named v4 species is the withheld one and that lies within that species' limit (NS31–AM1 1.71795%, full ITS 5.98803%), and a VT whose v4 holders all belong to that species goes whole. Records nearer a congener stay, because they would exist even if the species had never been sequenced.

The prediction, locked before any round beyond the pilot, is the paper's bound. With N6, false known among the 34 test species rises from 53/673 to about 71/673 at most: the 18 sequences whose novelty rested on an unmatched VT. The 34 calls that rested on names sharing no species (N3) could also flip, so the bound can fail, which is what makes it a test.

The *F. mosseae* pilot is E's first round and counts as development. Its C1 removes AY639162.1 and AY639270.1, which are M4's R0040 and R0041, so it must reproduce M4's 122 new lineages and 8 new species in *Funneliformis*. Its C3 has already run: no false known under the paper's arm, and the genus returned for all 130 copies.

## Outcomes

Scoring follows the paper and the EUKARYOME protocol, with four additions for the call types the new rules create.

- **Kept** (K, reference intact): a known species or group that contains the sequence's own species. **Demoted**: kept under F, not kept under the arm.
- **False known** (A, L, and K with its species withheld): a known species or group although the sequence's species is absent from the reference it was scored against.
- **Genus** (A, L, and K in genus rounds): right, wrong or none; for L any genus given is wrong. *Nanoglomus* has no genus outcome, as in the EUKARYOME test, because v4's *Glomus* sp. 5089\_2351A and 4406 may be *Nanoglomus*.
- **Paralogs** (R5): among motif carriers, the share called divergent, by label; where the culture gives the copy class (v4's 5.8S rule), recall of short-5.8S sequences and wrong reclassifications of main types.
- **Genus flags** (R3): flags on wrong genera (recall) and on right genera (false flags).
- **Triage** (R6): positives reaching naming; negatives reaching naming; negatives called new lineage.
- **The complex claim**, under F only: every false-known call is checked against inside predictions computed before the run. Under arms with R2 the claim holds by construction and is not an outcome.

The six additions:

1. U1 (unnamed culture) is not a known species: it is never kept and never false known, and is counted by stratum.
2. A group that includes an unnamed culture counts as a known group when it holds a named species.
3. N7 (the gate) counts as novel species, genus known, like N6.
4. A new lineage with AMF membership unconfirmed counts as a new lineage; on a negative it counts as a non-AMF sequence called new lineage.
5. G2 (R7) is a genus outcome like any other: right, wrong or none, reported with the R3 flags it carries.
6. Insufficient evidence (R8d) is neither known nor new: it is never false known and is counted by stratum.

**Uncertainty.** Rates by sequence and averaged over species and over cultures, with 95% intervals from 10,000 cluster-bootstrap replicates that resample species, and cultures where the source has them (fixed seed). Differences between an arm and F use the same replicates for both, so each interval is paired. Every outcome is reported with and without the flagged database records.

## Development estimates

On seen data, R2 removes every false-known call that lies outside its called species (EUKARYOME 13 to 1, LSU removal 53 to 48) at a cost of 1 to 7 kept sequences per test. R7 restores every genus the method loses where the genus is present, and R8 changes no culture call on the main arm. R3b flags 40 of 69 wrong genera in the genus rounds but only about a fifth in strata A and L, and R5 with an SSU fallback recognizes 93 of 95 *Rhizoglomus* motif carriers.

Every row was first computed from stored outputs of release 0.3.1 (commit 02c0844) by October scripts that were never committed; each first reproduced the frozen counts (13, 152, 4, 136, 53, 9, 71, 186, 487) and the paper's limits. Every row the arm code reproduces is now pinned in its tests, and the rows it does not reproduce are marked superseded in the paragraphs below the table.

| Rule | Data (seen) | Outcome | Frozen | With the rule | Note |
| --- | --- | --- | --- | --- | --- |
| R2 | v4, withheld cultures | Kept | 745/755 | 738/755 | The brief's estimate: 7 more demoted, in 4 cultures |
| R2 | v4, withheld species, V18 intact | False known | 18/803 | 18/803 | The gate alone also leaves 18 of the 71 placement calls, so it subsumes N6 here |
| R2 | v4, LSU removal, V18 removed | False known | 53/673 | 48/673 | The 5 removed lie outside the called species |
| R2 | v4, LSU removal, V18 intact | False known | 9/673 | 9/673 |  |
| R2 | LongGloDB, absent taxa | False known | 4/102 | 4/102 | All inside *G. rosea*, as predicted |
| R2 | LongGloDB, known species | Kept | 136/164 | 135/164 | The one lost, *F. mosseae* BEG12 (ASV1250), has all three windows and lies outside the species' SSU limit: a real cost of R2 |
| R2 | EUKARYOME, absent taxa | False known | 13/2,706 | 1/2,706 | EUK1266127 (*D. bareai*) stays, as predicted |
| R2 | EUKARYOME, known species | Kept | 152/199 | 149/199 |  |
| R2 | EUKARYOME, genus-only reads | Known calls outside the called species | 593 calls | 113 demoted | Species unknown; descriptive only |
| R1a | EUKARYOME, absent taxa | False known | 13 | 3; with R2, 1 | 10 *Septoglomus* reads lie inside only 4404\_2057A or FS\_9\_2328B and become U1. Computed by the arm code: known species kept 152 to 151 (one P. laccatum read inside BR105\_2435A's unnamed lineage, which R2 also removes), and no change on any v4 design or LongGloDB |
| R1a | EUKARYOME, *Nanoglomus* inside unnamed cultures | Calls | 42: 33 divergent, 8 new species, 1 boundary | 9 become U1 | Inside *Glomus* sp. 5089\_2351A (8 also inside 4406, which PROTOCOL\_EUK.md counts twice to reach 50; the paper's 42 is right); divergent calls keep their call |
| R3b | EUKARYOME, absent genera given a genus | Flagged | no flag | 181/268 | *Silvaspora* 87/87, *Orientoglomus* 50/114, *Racocetra* 0/11 |
| R3b | EUKARYOME, absent species, right genus | False flags | no flag | 1/49 |  |
| R3b | EUKARYOME, known species kept | False flags | no flag | 0/152 | Cut: pendant length 0.07488, 99th percentile of 803 main-type copies in the LSU-removal rounds |
| R5 | EUKARYOME, *Rhizoglomus*-labelled motif carriers | Called divergent | 46/95 | 66/95 by VT; 93/95 with the SSU fallback | Fallback: VT not matched and nearest SSU culture *Rhizoglomus* within the within-genus SSU limit: 1.970% on all of v4 by R3a's definition (the 2.108% recorded earlier does not follow from it; no development call depends on the difference) |
| R5 | EUKARYOME, other motif carriers | Reclassified | none | 1/438 | 408 *Nanoglomus* (1 reclassified), 17 *Entrophospora*, 7 *Oehlia* (none) |
| R6 | EUKARYOME new-lineage calls, *Nanoglomus* aside | Dropped by filtering triage at `--min-id` 85, 88, 90, 92 | none | 79, 607, 787, 907 of 1,259 | No known (199) or absent-species (77) read falls below any of these cuts |
| R7 | v4, species withheld, V18 intact | Genus for new-lineage calls | None of 365 | 179 right, 0 wrong, 186 abstain | Unit level; R3b flags none; the 186 are in genera with no other culture |
| R7 | v4, LSU removal, V18 removed | Genus, of 487 whose genus another culture holds | 326/487 | 487/487, 0 wrong | Label level; R3b flags none (in-sample cut) |
| R7 | LongGloDB, absent species in v4 genera | Genus, of 59 | 6/59 | 50/59, 0 wrong | Label level; R3b flags 3 of the 44 added |
| R7 | LongGloDB, absent genera | Wrong genus, of 43 | 0 | 0 | All 43 abstain |
| R7 | EUKARYOME, absent species in v4 genera | Genus for new-lineage calls | None | 24 right, 0 wrong | Label level; R3b flags none |
| R7 | EUKARYOME, absent lineages, Nanoglomus aside | Wrong genus, of 1,259 new-lineage calls | 0 | 258 | R3b flags 212 (82%); most given Rhizoglomus (106) or Ambispora (76) |
| R8 | v4, species clades: culture, species and LSU-removal rounds | Single-marker new lineage | 28, 379 and 360 new lineages | 0 flagged | The baseline arm's 3 changes do not occur on this arm |
| R8 | EUKARYOME and LongGloDB | Single-marker new lineage | 1,973 and 101 new lineages | 76 and 31 flagged | All absent taxa or genus-only labels, no known species; LongGloDB's VT is silent by amplicon design |

**R1b, computed by the arm code.** R1b changes no call on any v4 design or on LongGloDB. On EUKARYOME it takes false known among absent taxa from 13 to 3 of 2,706 (1 with R2), and known species kept from 152 to 150 of 199 (148 with R2). All 35 blocked calls come from three LSU placement units, each held only by one unnamed culture: E218 (4404\_2057A), E975 (BR105\_2435A) and SC209 (FS\_9\_2328B). They include the 10 *Septoglomus* false knowns that R1a also removes, and 2 correct *P. laccatum* calls that share edge E975 with BR105\_2435A's unnamed lineage, one more than R1a costs.

**R5, computed by the arm code.** R5 reproduces its EUKARYOME estimates and changes no main-type call on any design. On culture data it turns only true divergent copies divergent, all of them previously called new species or new lineages. Short-5.8S copies recognised under species holdout rise from 29 to 32 of 43 by VT, and to 40 with the SSU fallback; under culture holdout from 36 to 36 and 40. PX214776.1 is already divergent through D1 in every design, so R5 has nothing to move there.

**R3, computed by the arm code.** On the frozen genus calls, R3b flags 40 of 69 wrong genera in the 15 genus rounds and 199 of 1,063 absent-genus calls on EUKARYOME. It falsely flags at most one right-genus call: 1 of 326 in the LSU-removal rounds (the paper's first-outcome folds, A. callosa apart), none in either holdout, and none of 47 absent species or 152 known species on EUKARYOME. Its cut is 0.074880850445, itself an observed pendant length, so the eight sequences at exactly that value are not flagged. R3a on any marker flags all 69 wrong genera, but also 19 of 487 right-genus calls in the species holdout and 6 of 47 on EUKARYOME; the species-holdout flags are all *Septoglomus viscosum*, whose SSU lies beyond the genus limit. R3a on all markers flags 2 of 69. The R3b rows on EUKARYOME in the table (181 of 268) were measured with R7 on, and are rechecked once R7 exists.

**R7, computed by the arm code as written (unit level).** In the species holdout, 179 of 365 new lineages get the right genus and none a wrong one. In the LSU-removal rounds every genus another culture holds is restored, 326 to 487 of 487. Both match the rows above. On LongGloDB, 51 of 59 absent species in v4 genera end with the right genus (7 frozen; one read more than the row above in both columns, not traced), and none of the 43 absent genera gets one. On EUKARYOME, 27 of 29 new lineages of absent species in v4 genera get the right genus, and only 24 of 1,468 absent lineages (*Nanoglomus* aside) get a wrong one, 17 of them flagged by R3b. The label-level EUKARYOME rows above (258 of 1,259, and 24 right) came from an approximation that does not reproduce from the published data, and are superseded. In the 15 genus rounds, where every genus is wrong, R7 gives 12 genera and R3b flags none of them: the known failure exists, and R3b does not catch it there.

**R8, computed by the arm code.** It reproduces the rows above on v4 (none of 28, 379 and 360 new lineages flagged) and on LongGloDB (31 of 101). On EUKARYOME it flags 77 of the 2,186 new lineages in the paper's primary table: 70 absent genera and 7 genus-only labels, none a known species. The recorded 76 of 1,973 comes from no published table (the secondary analysis gives 74 of 2,183) and is superseded. R8d demotes the same calls to insufficient evidence and changes no false-known count.

**V, computed from the published secondary analyses.** The external tests already ran V as their secondary VT analysis. Against the primary analysis, V keeps three more known species on LongGloDB (136 to 139 of 164), but adds one false known there (4 to 5 of 102) and one on EUKARYOME (13 to 14 of 2,706, kept unchanged at 152 of 199). Its genus error is larger than one case: all 18 LongGloDB wrong-genus calls under V are *Racocetra gregaria* called *Cetraspora*, new lineages under the primary analysis, and neither R3a nor R3b flags any of them. On development data, the pairing with R3 does not contain V's known failure. V's arm code is `amfsplit name --vt-spanned`: on the LongGloDB windows it reproduces the test's spanned evidence for all 280 sequences, best units included, and every shared column of the published secondary analysis.

**R4, cross-fitted for N6 (prototype, 3 November 2026).** The pool's other laboratories are LongGloDB's KU, ZU and ZUII cultures: 56, of which the six name-check exclusions leave 50, all among the test cultures. Fitted on v4 (4.441, 5.020 and 5.993% at the 95th, 97.5th and 99th percentiles), LongGloDB keeps 4 of 102 false known and 136 of 164 kept at every percentile. Fitted on LongGloDB alone (7.258, 8.222 and 19.923%, from only 187 distances), v4 keeps 750, 750 and 751 of 755, and its species-holdout false known becomes 21, 33 and 71 of 803. Decision 6 therefore picks the 95th percentile. Pooled on v4 and the 50 cultures, with 185 of 977 nearest distances crossing collections, N6's threshold is 5.031% at the 95th percentile (6.048% at the 97.5th, 7.404% at the 99th), against the frozen 5.993%. R4's arm code verifies the name checks, then reproduces the frozen N6 exactly when fitted on v4 at the 99th percentile: no call changes in 198 holdout rounds. At the decided 95th percentile, pooled, it lowers false known (species holdout 18 to 13 of 803; EUKARYOME 13 to 9 of 2,706) but demotes more known calls than F (culture holdout kept 745 to 734 of 755, interval -3.21 to -0.39 points; EUKARYOME 152 to 144 of 199). As decided, R4 therefore fails its own adoption rule, which requires demotion of K to fall against F's limits. Whether to keep decision 6 or amend it before the lock is open. With R2 (F+R2+R4), the pooled inside limits at the 95th percentile cost 48 correct calls (culture holdout kept 697 of 755, interval -12.25 to -2.86 points; species-holdout false known 8). At the 97.5th percentile R4 does what it is for: alone it keeps 747 (false known 18; EUKARYOME 14 and 152), and with R2 741 against R2's own 738 (false known 18; EUKARYOME 1 and 148). At the 99th it keeps more (750 alone, 749 with R2) but false known rises (22 and 21; EUKARYOME 18 and 5). An amended decision 6, the lowest false known among percentiles whose pooled limits demote no more than F's, would select the 97.5th.

**R6, positives (4 November 2026).** Triage of the 5,511 EUKARYOME overlap reads (protocol\_v3/r6\_triage\_positives.tsv) labels 4,351 AMF, 1,159 unresolved and 1 non\_AMF. All 199 known-species reads and all 77 absent-species reads are AMF, so known species kept (152) and false known (13) do not change. Of the 2,629 absent-lineage positives, filtering names 1,496; R6 names 2,628 (interval +25.6 to +72.6 points) and reports 847 of them as new lineages with AMF membership unconfirmed. The 1,132 unresolved positives are real Glomeromycota that V18 represents poorly (*Polonospora* 281, *Nanoglomus* 156, Glomeraceae genera 25 and 14). The one non\_AMF read's best hit is *Geosiphon*, a member of Glomeromycota that V18's AMF-only file leaves out. The negatives measure R6's cost. Of the 394,660 non-Glomeromycota reads of EUKARYOME v2.1, 19,402 have full ITS, LR0R-FLR2 and SSU flank windows of at least 300 nt; triage labels 18,768 non\_AMF, 495 unresolved and none AMF (139 have no hit). Filtering names none of them; R6 names the 495 unresolved, all fungi (Basidiomycota 152, Chytridiomycota 131, Rozellomycota 108, Mortierellomycota 34). Of these, 485 are called new lineages, which R6 reports with AMF membership unconfirmed, and 10 artifacts; none gets a species or a genus. R6's cost is therefore 2.6% of non-AMF fungi reaching naming, every one labelled or set aside.

**S, computed by the arm code (4 November 2026).** S applies its rules in a fixed order: R5 (variant), R4, R1a, R2, R6, R7, R8, R3 (variant). R6 runs before R7, so a sequence whose AMF membership triage could not confirm never takes an AMF genus from R7; on EUKARYOME this cuts R7's wrong genera on absent lineages from 24 to 7. With the development variants (R5 with the SSU fallback, R3b), S keeps 741 of 755 in the culture holdout and leaves species-holdout false known at 18 of 803. It recognises 40 of 43 short-5.8S copies in both holdouts, and gives 179 of 365 new lineages the right genus and none a wrong one. On EUKARYOME, false known falls from 13 to 1 of 2,706 and kept from 152 to 148 of 199; 847 absent lineages are reported with AMF membership unconfirmed. On LongGloDB, with R4 fitted on v4, false known stays at 4 of 102, kept falls from 136 to 135, and 51 of 59 absent species get the right genus. In the LSU-removal rounds, false known with V18 removed falls from 53 to 36 of 673, exactly R4 with R2 (R2 alone gives 48). In the genus rounds, R3b flags 40 of 75 wrong genera with V18 removed, but none of R7's 6 with V18 intact.

**Dress rehearsal, first findings (4 November 2026).** With LongGloDB v1 and EUKARYOME v2.1 standing in, the strata are K 164 and 199, A 59 and 77, L 43 and 2,629, G 14 and 2,606, plus the 19,402 negatives. Three gaps must close before step 9 can run as written. First, withheld K does not exist for non-culture sequences: condition B of both external tests removes the absent taxa's relatives from V18 (102 and 2,123 rows), but no K sequence has been scored with its own species withheld, and v0.1 cannot yet prune the culture reference and V18 per round and re-place the queries; the same holds at genus rank for R3. Second, the scoring module scores each source apart, while section 8 speaks of pooled strata. Third, the R5 carrier and V outcomes exist only in tests, not in the scoring module. Applying section 8 to the estimates already in hand, R6 (485 of 19,402 negatives called new lineage, 2.5%), R7 (1 of 7 added wrong genera flagged on EUKARYOME, 0 of 6 in the genus rounds) and V (0 of 18 VT-only wrong genera flagged) fail as written, and R1a fails on EUKARYOME's K alone (added demotion bound 1.56 points) pending the pooled K. These are rehearsal results and decide nothing.

**Withheld K, the first gap closed (4 November 2026).** `amfsplit withhold` rebuilds the LSU-removal rounds as the frozen runs did. Across all 51 rounds, all 102 evaluation tables are byte-identical; the placement summaries match up to tied best edges, except likelihood weights differing in their last digits in 4 rows of 3 rounds (EPA-ng's floating-point summation), which change no call. All 24 K species of the stand-in sets have a species fold, so `amfsplit withheld-k` names each K sequence with its species withheld in its own round (condition W). The frozen method then calls 18 of 363 K sequences known. On LongGloDB all 8 of 164 are MR101, labelled *F. coronatus* but *F. mosseae* by sequence (the R4 name checks), so its true figure is 0 of 156. On EUKARYOME they are 10 of 199: *R. irregularis* 5 of 22, *G. macrocarpum* 3 of 7 and *F. caledonium* 2 of 16, all species with close relatives in the reference. Pooled strata scoring and the report remain.

**Pooled strata, the second gap closed (4 November 2026).** `--score strata` pools LongGloDB and EUKARYOME by each row's stratum and condition. The frozen method keeps 288 of 363 K sequences with the reference intact, calls 18 of 363 withheld K known, and makes all 17 of its absent-taxa false knowns in A (17 of 136), none in L (0 of 2,672). S keeps 283 of 363 K (interval -3.02 to -0.45 points), cuts withheld-K false known to 12 and A to 5, and recognises 108 of 110 *Rhizoglomus* carriers. Read against section 8 before each rule is run alone, three things stand out: S's demotion bound (3.02 points) sits just over R2's 3-point ceiling; R3b flags 200 of 1,070 wrong genera in A and L (19%, against 50%); and R7's added wrong genera are flagged in 2 of 8 (against 80%). `amfsplit adopt` will apply section 8 rule by rule.

**Dress rehearsal done: procedure step 5 (4 November 2026).** `amfsplit adopt` applied section 8 to the development stand-ins, each rule alone on the pooled strata and S last; the report and every score are in `protocol_v3/rehearsal/`. It decides nothing. Mechanically, three rules would become defaults: R2 (demotion bound 2.34 points; false known in A, L and withheld K 35 to 22), R5 with the SSU fallback (98.2% of carriers; 0.23% of others reclassified) and the R8 flag (R8d not needed: no single-marker new lineage among K). S, built from these three, keeps the same 284 of 363 and passes. R1a fails narrowly (bound 1.05 against 1 point) and R1b clearly (2.10). R4 at the 97.5th percentile fails, since it leaves K's demotion unchanged and adds one false known. R5 by VT alone fails (72.7%), and R6 fails its cap (2.50% of negatives, bound 4.55%), so filtering stays. No R3 variant qualifies: R3a catches 61.6% of wrong genera but flags 2.17% of right ones, against 2%; R3b catches 21.1%. R7 and V therefore stay opt-in, and would fail their own checks besides (R7 flags 18 of 25 of its wrong genera; V flags 5 of 46 of its VT-only ones).

**Limits fitted on all v4 cultures.** Release 0.3.1 has two pools for the within-species limit, and both are frozen. N6 (`its_rule.py`) pools every main-type copy; the complex rule and the inside test (`species_groups.py`) pool unique sequence–culture pairs. The within-genus limit of R3a uses unique sequence–species pairs over the 8 genera with two or more named species.

| Marker | Within species, N6 pool | Within species, complex pool | Within genus, 95th / 97.5th / 99th |
| --- | --- | --- | --- |
| Full ITS | 5.99255% | 5.98803% | 18.913 / 20.288 / 20.868% |
| LR0R–FLR2 | not used | 3.02698% | 7.659 / 10.526 / 11.504% |
| SSU region | not used | 0.80652% | 1.970 / 2.038 / 2.108% |

**Still to compute, and where.** Measurement E remains to run, against its prediction locked on 3 November. The per-stratum predictions for the confirmation set are the rates of the dress rehearsal (protocol\_v3/predictions\_v3.tsv); they become counts at the set lock (section 9).

## Adoption rules

Each rule is judged on the confirmation set by its benefit against its cost, using paired 95% intervals against F. Development estimates set expectations but decide nothing. Its numbers were settled on 2 November 2026 (section 10).

- **R1a, R1b, R2.** These only remove species calls, so false known cannot rise. Each is adopted if its added demotion of K, reference intact, has an upper bound of at most 3 points for R2 and 1 point for R1. If both R1 variants pass, R1a is preferred, because it names the matching culture.
- **R4.** False known can rise here, so the rule is two-sided. Adopted if demotion of K falls against F's limits, and false known in A, L and withheld K rises by at most 1 point with an upper bound under 3 points.
- **R5.** Adopted if at least 80% of *Rhizoglomus*-labelled motif carriers are called divergent and at most 1% of other carriers are reclassified. Of the two genus checks, VT only or with the SSU fallback, the one with higher recall that meets the cap is adopted.
- **R3.** The flag is on by default if it catches at least 50% of wrong genera, in genus rounds and strata A and L, while flagging at most 2% of right genera. If both variants qualify, the one with higher recall is used; otherwise flags stay opt-in.
- **R7.** Adopted with the chosen R3 variant if it gives the right genus to at least 50% of F's genus-less calls whose genus another culture holds, and at least 80% of the wrong genera it adds carry the provisional flag. Without an adopted R3, R7 stays opt-in.
- **R8.** The flag is on by default, because it changes no call. R8d is adopted only if single-marker new lineages turn up among K sequences with the reference intact; on development data every one was an absent taxon or a genus-only label.
- **V.** Adopted with R3 only if every wrong genus given through a VT alone is flagged, and no false-known call arises through a VT shared by species.
- **R6.** Labelling replaces filtering if at most 1% of triage negatives are called new lineage, with an upper bound under 2%. Otherwise filtering stays, and every run reports what triage dropped, with the identities.
- **S.** S checks that the adopted rules work together. If S demotes more of K than the adopted rules' largest ceiling allows, the interacting rules are named and the later one in the order R1, R2, R4, R5, R7 stays opt-in.

**Uninformative outcome.** A rule whose data fall below these minimums is reported as untested and stays opt-in, whatever its point estimates:

- K: 100 sequences of 10 species from 2 collections.
- R1 and R2: at least one false-known call of the kind they target.
- R5: 20 *Rhizoglomus* motif carriers.
- R3: 30 wrong-genus calls.
- R7: 30 genus-less calls whose genus another culture holds, and 30 from genera the reference lacks.
- R6: 200 negatives.

**Power check.** When the confirmation set is assembled (section 9, step 7), each rule's power is simulated on its stratum counts. Power is the chance of meeting the rule's cap at its development rate, with the cluster bootstrap of section 6 and a fixed seed. A rule under 80% power is declared underpowered before the lock: it still runs and is reported, but stays opt-in, and the report says the test could not decide, not that the rule failed. The simulation and its result are hashed into the lock.

## Procedure and locking

The confirmation set is processed once, after steps 1 to 7 are done, hashed and published.

1. Settle the decisions in section 10 and replace each proposal with its final value. Done on 2 November 2026.
2. Done: Gate 0b passed and v0.1.0-rc.1 is tagged (prerequisite 1); the lock records its commit.
3. Write the arm code and compute the pending development estimates of section 7 with it, locking measurement E's prediction before its rounds run; commit the scripts and their outputs.
4. Done: the R4 name checks are curation edits, each with its check (r4\_name\_checks.tsv).
5. Dress rehearsal: run steps 9 and 10 on development data, with LongGloDB v1 and EUKARYOME v2.1 standing in for the confirmation set, and report it as a rehearsal. It tests the arm code, the scoring and the report, and decides nothing.
6. Protocol lock. Put the SHA-256 of this protocol (exported to Markdown as `PROTOCOL_v3.md`), the eligibility rules, the per-stratum predictions as rates, the name-check edits, the arm code and the v0.1 commit into `protocol_hashes.txt`. Commit, push and archive it.
7. Set lock, when the confirmation sources are available (at the latest, the first LongGloDB release after v1). Assemble the confirmation set and the triage negatives by their rules (identifiers, labels, stratum counts), turn the predictions into counts and run the power check of section 8. Add both set manifests, the predictions and the power check to `protocol_hashes.txt`, then commit, push and archive it before any confirmation sequence is cut.
8. Rescore the stored development results with the confirmation set present and excluded; they must not change.
9. Run v0.1 with its arm code once, on every condition and withholding round of the confirmation set and on the negatives, and compute all arms from that output.
10. Report every arm and outcome, apply section 8, and record which rules become rules v3 defaults and which stay opt-in.

**Amendments and retuning.** A departure from this protocol is written, dated and hashed before the result it affects is seen. Nothing is retuned after the run; as protocol v2 put it, a change needs a new protocol version and fresh data.

## Decisions before locking

Fifteen choices, all decided: two by 20 October and the other thirteen on 2 November 2026. Their values replace the proposals throughout.

| Decision | Options | Recommended | Reason | Status |
| --- | --- | --- | --- | --- |
| Missing regions in R2 and R1a | Check the markers present and flag the rest; or fail on any missing marker | Check the markers present, flag the rest, fail when two are missing | *F. mosseae* BEG12 was thought to lack ITS, but it has all three windows and fails R2 on SSU; no development sequence lacks a marker, so the rule guards window failures in new data | Decided |
| R3a, "beyond the limit" | On any marker; on all markers | Any marker, with all markers reported as a sensitivity | Within-genus limits are wide (ITS 20.9%, LSU 11.5%), so one marker beyond them is informative | Decided |
| Benefit needed for R1 and R2 | At least one targeted error removed on the confirmation set; or cost alone | At least one | A default should rest on benefit seen on unseen data; cost alone is a weaker claim | Decided |
| Ceilings and caps marked (p) | As proposed in section 8; other values | As proposed, with the power check of section 8 | R2's 3 points leaves room above its 0.9-point culture cost for spread across collections, where N6 lost 9.1%. In simulation, R2 at its development rate met that cap in only 55 to 75% of draws at 100 to 500 sequences, hence the power check | Decided |
| Minimum sizes | As proposed in section 8; other values | As proposed, with the power check of section 8 | Small enough to reach, large enough for an interval | Decided |
| R4 percentile | 95th, 97.5th or 99th, by a rule | Amended 3 November 2026, from the development estimates and before the lock: the lowest cross-fitted false known among percentiles whose pooled limits demote no more known calls than F's, alone and with R2 (culture holdout). This selects the 97.5th. The original rule (cross-fitted demotion at most 3%) selected the 95th, whose pooled limits demote more than F's (kept 745 to 734 alone, 738 to 697 with R2), against R4's own adoption rule | Uses both directions of the cross-fit | Decided |
| R4 fit exclusions | The four code-sharing cultures, MX982A and MR101 | Exclude all six; report MR101 relabelled *F. mosseae* as a sensitivity | Each shared code is a duplicate or a swap, and MX982A and MR101 conflict with v4 names | Decided |
| Triage negatives | A taxon list | All non-Glomeromycota EUKARYOME v2.1 long reads with three regions of at least 300 nt, by phylum; Endogonomycetes and Mortierellomycota apart | The closest relatives, and the Stage 1 contaminant | Decided |
| Triage reference for R6 | V18 with its four outgroups; a broad `TRIAGE_DB` | V18 with its outgroups | No new dependency; a broad reference needs its own leave-out test | Decided |
| Gate 0 rebuild tolerance | Exact; or differences only at the 0.95 weight cut, each listed | Exact, with pinned binaries | Settled 20 October: the pinned RAxML-NG reproduced every artifact, and an unpinned build moved PX207259.1 | Decided |
| Confirmation sources | INSD after EUKARYOME v2.1, LongGloDB after v1, EUKARYOME after v2.1, new sequencing | Search all four; the INSD cutoff is 9 September 2026, the EUKARYOME v2.1 release | EUKARYOME announced Release 2.1 on 9 September 2026 ([news](https://sisu.ut.ee/eukaryome/news-0/)) | Decided |
| Taxonomy convention | The frozen `normalise_genus` and the EUKARYOME `names.tsv` rules; or an external authority | The frozen convention (*Rhizoglomus* = *Rhizophagus*, *Claroideoglomus* = *Entrophospora*) as one versioned crosswalk file | The paper and the Rust port already apply it; the 23 genus "mismatches" came from comparing raw strings | Decided |
| R7 support minimum | 1, 2, 3 or 5 training cultures | 2 | In the 22 October screen, 1 and 2 gave identical calls, while 3 and 5 only lowered coverage | Decided |
| Measurement E removal rule | By name; by whole VT unit; by sequence (nearest named species, within its limit) | By sequence, with a whole VT only when its v4 holders all belong to the species | Environmental records rarely carry the species name, and a congener's records would exist anyway | Decided |
| R8 default | Flag only; or demote (R8d) | Flag only | On the main arm, every single-marker new lineage in development data was an absent taxon or a genus-only label | Decided |

## What this protocol cannot show

- **Undescribed species.** Outside measurement E, MaarjAM and UNITE stay intact, so a withheld species keeps any VT or SH of its own. E removes records by their nearest named species, which approximates a species nobody has sequenced without equalling one.
- **New species from K.** Withheld K sequences are new sequences from new collections, not new species; only strata A and L test unseen species.
- **The complex claim under R2.** With the gate on, every remaining false-known call lies inside its called species by construction, so only arm F tests the claim.
- **Copy classes outside cultures.** Sequences without a culture are scored as main type, as in the EUKARYOME test, so divergent copies among them are counted as they are called.
- **Out-of-scope inputs and outputs.** ONT consensus input, a broad triage reference, community profiles, richness and abundance.
