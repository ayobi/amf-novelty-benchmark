# LSU-removal test: per-sequence results (aggregate_v2)

One row per culture sequence for every species round of the LSU-removal test (protocol v2), as aggregated by
the run: `species.v1` and `species.v2` with the species removed from the LSU reference (v1 without, v2 with the
ITS distance rule), and `species.control_v1` and `species.control_v2` with the LSU reference intact.
`species.summary.json` holds the per-round summaries behind Table S4.

Used by `paper/analysis/foldwise_complexes.py` and `paper/analysis/intervals.py` with
`--m4 ../benchmarks/m4_external_holdout/aggregate_v2` (sequence-level matching of the complex rule, and
intervals that resample cultures).
