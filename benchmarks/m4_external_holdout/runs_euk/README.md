# EUKARYOME long-read test: run outputs

Outputs of `euk_test.py run` (protocol `benchmarks/eukaryome/PROTOCOL_EUK.md`), collected with
`euk_test.py collect` and reported by `euk_report.py` (`report/report.md`). No sequences are stored here:
the evidence, placements and score tables refer to EUKARYOME v2.1 sequence identifiers only.

- `invocation.json`: command, inputs and their SHA-256 (they equal the locked files in `benchmarks/eukaryome`).
- `01_evidence/`: VT and SH evidence for the test sequences, with the parity check against the frozen evidence.
- `A/`, `B/`: LSU placement with the reference intact (condition A) and in each removal round (condition B).
- `06_scores/`: one row per sequence for each condition, analysis (primary, secondary VT) and rule version
  (v1 without, v2 with the ITS distance rule).
- `report/`: `euk_report.py` output; `report.md` is the readable version.

Tables over 200 kB are gzip-compressed. Before rerunning `euk_report.py` or the paper scripts
(`paper/analysis/euk_numbers.py`, `euk_followup.py`, `make_fig_predicted.py`), decompress them in place:

    find benchmarks/m4_external_holdout/runs_euk -name '*.tsv.gz' -exec gunzip -k {} \;
