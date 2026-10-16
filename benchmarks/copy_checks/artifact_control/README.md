# What the artifact screen catches (controls added after review)

Exploratory, written after review (point 1.6), on 15 October 2026. It changes no rule and no result.

The screen (novelty step 1, frozen `decision_tree.py`) flags a copy when its SSU flank differs by at least 1% from
the nearest main-type copy of another culture and, for every SSU-tied nearest copy, by at least half as much as
the full ITS. The paper reported that it fires for no real copy. This asks what it does fire for.

    python3 stage3/artifact_control.py --ref reference/culture_ref_v4 --out benchmarks/copy_checks/artifact_control

Seed 20261015; 30 s on 2 cores. Each query is scored with its parents' culture labels withheld, against the 852
main-type copies of v4 (`artifact_control.queries.tsv`, one row per query; `artifact_control.summary.tsv`).

| Queries | n | Flagged | 95% CI (Wilson) | Median SSU distance to nearest |
|---|---:|---:|---|---:|
| Real main-type copies (negative control) | 852 | 0 | 0–0.4% | 0.0% |
| Random errors, 0.5% of positions | 148 | 1 (0.7%) | 0.1–3.7% | 0.54% |
| Random errors, 1% | 148 | 18 (12.2%) | 7.8–18.4% | 1.02% |
| Random errors, 2% | 148 | 43 (29.1%) | 22.3–36.8% | 2.04% |
| Random errors, 3% | 148 | 68 (45.9%) | 38.1–54.0% | 2.99% |
| Random errors, 5% | 148 | 106 (71.6%) | 63.9–78.3% | 5.04% |
| Random errors, 8% | 148 | 136 (91.9%) | 86.4–95.3% | 7.92% |
| Chimera, SSU of A + ITS of B, other genus | 150 | 0 | 0–2.5% | 0.0% |
| Chimera, SSU of A + ITS of B, same genus | 150 | 0 | 0–2.5% | 0.0% |
| Chimera, breakpoint inside the SSU, other genus | 150 | 10 (6.7%) | 3.7–11.8% | 2.17% |
| Chimera, breakpoint inside the SSU, same genus | 150 | 0 | 0–2.5% | 0.27% |

Random errors: one main-type copy per culture label (148), 60% substitutions, 20% insertions, 20% deletions,
at the same rate in the SSU flank and the full ITS. Chimeras: parents are main-type copies of different named
species, drawn at random.

## Reading

- The screen flags copies whose SSU is about as divergent as their ITS. Errors spread evenly along a read do that,
  once they reach a few percent: 29% of copies flagged at 2%, 92% at 8%. At 0.5% (still far above HiFi error
  rates) it flags almost none.
- It does not catch chimeras. A copy whose SSU comes intact from one culture and whose ITS comes from another sits
  at 0% SSU distance from its SSU parent's species, below the 1% floor (0 of 300). With the breakpoint inside the
  SSU, 136 of 300 pass the 1% floor, but the nearest SSU copy usually belongs to a genus whose ITS is far more
  divergent, so the ratio test passes them (10 of 300 flagged).
- So the step is a screen for SSU–ITS discordance, sensitive to noisy reads, and not a chimera screen. Of the real
  main-type copies, 26 are at least 1% from every other culture in the SSU flank; none is flagged.
