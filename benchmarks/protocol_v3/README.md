# Protocol v3 locks

Protocol v3 is the locked test of the amfsplit rule changes R1 to R8. It stays a draft until its
protocol lock, when it is exported here as `PROTOCOL_v3.md` and hashed into `protocol_hashes.txt`
with the arm code and the predictions. The set lock adds the confirmation set's manifests and the
power check, before any confirmation sequence is cut.

| File | Locked | What |
|---|---|---|
| `E_PREDICTION.md`, `E_prediction.lock` | 2 November 2026 | Measurement E's prediction, before its rounds run |
| `r4_name_checks.tsv` | 3 November 2026 | R4's name checks: the six LongGloDB cultures left out of the fit, each with a check R4 verifies before fitting (inside means within the complex rule's limits on ITS, LSU and SSU) |
| `r6_triage_positives.tsv`, `.summary.json` | 4 November 2026 | R6's development input: `amfsplit triage` of the 5,511 EUKARYOME overlap reads (LR0R-FLR2 against V18_LSUDB_052025 with its outgroups, vsearch, --min-id 90; AMF reference SHA-256 e6794cbb...): 4,351 AMF, 1,159 unresolved, 1 non_AMF |
| `r6_negatives_triage.tsv.gz`, `.summary.json`, `r6_negatives_named_rows.tsv.gz` | 4 November 2026 | R6's negatives: the 394,660 non-Glomeromycota reads of EUKARYOME v2.1 (`General_EUK_longread_v2.1.fasta`) through Stage 1 (ITSxRust 0.3.0, HiFi preset), the 19,402 with full ITS, LR0R-FLR2 and SSU flank windows of at least 300 nt triaged as the positives were (19,263 with a hit: 18,768 non_AMF, 495 unresolved, none AMF), and the named rows of the 495 that reach naming under R6 |
