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
