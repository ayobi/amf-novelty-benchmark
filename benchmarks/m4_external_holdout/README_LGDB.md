# LongGloDB independent-culture test: runner

Implements `PROTOCOL_LGDB.md` with amendments 1 and 2. Unzip into `amf_m4_step1`, next to
`external_holdout.py`. Nothing in the frozen package changes; `verify_code.py` still passes.

| File | What it does |
|---|---|
| `lgdb_test.py` | `check`, `run` (evidence, condition A, condition B folds, scores) and `collect` |
| `lgdb_report.py` | `lgdb_test.py report`: the protocol's outcomes, compared with the lock |
| `tests/preflight_parity.py` | no external tools; adds the test ASVs to the frozen data and requires the 37 M4 species results back exactly |

## What the run does

1. **Evidence.** BLASTs the 903 reference windows and the 280 test windows against the MaarjAM and
   UNITE files of the v4 crosswalk, with the frozen loaders, search and filters. It stops unless
   the frozen VT and SH evidence of all 903 reference copies comes back exactly. The secondary VT
   analysis (amendment 2) uses a second search with coordinates, checked hit for hit against the first.
2. **Condition A.** The 903 reference windows and the test windows are added to the frozen V18
   alignment and placed on the frozen tree and model. Scores use the frozen placements of the
   reference copies; `A/parity.json` shows whether the combined run treated them alike.
3. **Condition B.** One fold per novel test species whose label V18 holds (19 folds; *Dentiscutata
   erythropa* as V18's *D. erythropus*): its V18 records removed, alignment and model refitted, all
   windows re-placed, exactly as in the M4 folds.
4. **Scores.** The test ASVs join the frozen data in memory, and every test culture is excluded
   from every call. Before any test result, all 37 M4 species folds are rescored and must equal
   the stored M4 v2 results. Then A and every B fold are scored with the frozen `infer()`, before
   and after N6, in the primary and the secondary VT analysis.

Completed stages are verified by hash and reused; a changed input stops the run.

## Commands

```bash
cd ~/Downloads/AMF_test/amf_m4_step1
export PATH=$HOME/miniforge3/envs/amf-m4/bin:$PATH:$HOME/miniforge3/envs/amfsplit/bin
LG="--lgdb ../bench/longglodb --maarjam ../ref/maarjam/maarjam.fasta --unite ../ref/unite/UNITE_public_19.02.2025.fasta.gz"
python3 verify_code.py
python3 tests/preflight_parity.py ../bench/longglodb      # about 2 minutes
python3 lgdb_test.py check $LG
nohup python3 lgdb_test.py run $LG --threads 4 --jobs 8 > runs_lgdb.log 2>&1 &
# when it ends:
python3 lgdb_test.py report --lgdb ../bench/longglodb
python3 lgdb_test.py collect --zip ~/Downloads/lgdb_results.zip
```
