# Step 1: one external-reference holdout

> **Note, 1 November 2026.** This README was written for the pilot. The full runs' results are in `aggregate_v2/` (species rounds, as published) and `rounds/` (every round, release 0.3.2); see the two sections at the end.

Run only **Glomus chinense**. Its two V18 reference records and all its culture links are excluded. The existing M2b species-clade rules and thresholds remain frozen.

The fresh 493-reference alignment and fitted model are already completed and included. The remaining work aligns all 903 culture copies into that reference, performs EPA placements, and scores the 35 G. chinense main copies against the matched full-reference control. Training copies must also be re-placed because their LSU units depend on the new tree.

## Run

```bash
unzip -n "$HOME/Downloads/AMF_M4_Step1.zip" -d "$HOME/Downloads/AMF_test"
cd "$HOME/Downloads/AMF_test/amf_m4_step1"
conda env create -f environment.yml
conda run --no-capture-output -n amf-m4 bash run_pilot.sh
```

If `amf-m4` already exists, skip environment creation. This uses a separate environment from amfsplit; it does not depend on your existing `.venv`. The environment pins MAFFT 7.526, RAxML-NG 1.2.2 and EPA-ng 0.3.8. The runner verifies versions before starting.

Expect the remaining alignment and placements to take several minutes; the complete PC stage has not been timed. Default is eight threads, with query alignment using up to four. Keep the PC awake. Only one species is run.

When it finishes, **upload `amf_m4_step1_results.zip` from this directory**. Then we will inspect the result before running more species, genus holdouts or independent material. The detailed results are also in `runs/species_glomus_chinense/06_evaluation/`.

## Resume or report a failure

Repeat the same run command to resume. Completed stages are hash-checked; an interrupted generated stage restarts. Do not edit the delivered code or inputs and reuse cached results.

To collect logs if it fails:

```bash
conda run -n amf-m4 python3 external_holdout.py collect \
  --out runs --zip amf_m4_step1_results.zip
```

Full placements remain on your PC; the upload ZIP contains compact evidence and logs. Nothing is uploaded automatically.

## What is already verified

- Classifier code and biological inputs are frozen against the delivered M2b package.
- Both G. chinense external tips are removed; 493 reference taxa remain.
- The new reference alignment preserves retained sequences, and the fitted tree preserves the pruned topology.
- Six tests cover exclusions, aliases, co-cultured lineages, pruning and interrupted-stage recovery.
- The pilot wrapper checks saved reference/model hashes before continuing.

The full species/genus protocol is included for provenance. This delivery runs just this pilot and **does not contain its final classification result yet**. The earlier no-deletion control and case audits informed preparation but are not substituted for the unfinished pilot.

This remains a development test on the existing culture collection and a pruned published backbone. VT/SH external databases are unchanged. It does not establish independent species accuracy or make likelihood weights probabilities of correct species names.

## Per-round outputs (`rounds/`, release 0.3.2)

`rounds/<key>/` holds, for each of the 52 rounds of `holdout_plan.json` (37 species, 15 genus), the round's LSU placement table (`05_summary/lsu_placements.copies.tsv`), its DONE and metrics files, and its four version 2 per-sequence tables. The *Ambispora callosa* round reuses the full-reference placement, so it has no table of its own. `pack_round_placements.py` copied them from the original runs after checking, for every round, the result hash and that the placement table covers the 903 culture sequences, and for the 37 species rounds that every row equals the published `aggregate_v2` tables. `rounds/MANIFEST.json` records each file's SHA-256 and source. Alignments, jplace files and trees stay on the PC.

A rerun of the 36 species rounds with their own placement, with the same inputs and the same MAFFT, RAxML-NG and EPA-ng executables but 8 threads per job instead of 4, reproduced every call. Placement statistics agreed to about 1e-12 (2 cells differed), and the reported best edge changed for 50 of 32,508 sequence-rounds with its likelihood weight unchanged, consistent with ties broken in a different order.

## Robustness of the ITS distance rule

- `n6_decompose.py` reruns the frozen rule, with no mismatch against the 850 published rows of `aggregate_v2`, and changes one thing at a time: equal weight per culture, and one main-type copy per training culture in the calibration, the query pool, or both (100 runs, seed 20261031). It writes `audit_outputs/n6_decompose_table.tsv` and `n6_decompose_results.json` here. Thinning only the calibration switches the rule off, because with one copy per culture its 99th percentile rests on 131 values; thinning only the query pool keeps all 81 firings in every run; equal culture weights give the published calls.
- `audit_n6_nearest_references.py` and `audit_n6_pool_comparison.py`, run from the repository root, write `audit_outputs/n6_nearest_references.tsv` and `audit_outputs/n6_all_vs_main_pool.tsv` there. `audit_n6_distance_sensitivity.py` is the check of whether the edit distance depends on argument order.
