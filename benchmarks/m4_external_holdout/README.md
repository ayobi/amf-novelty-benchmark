# Step 1: one external-reference holdout

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
