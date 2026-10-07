# M4 external-LSU holdout protocol, version 2

Written 29 September 2026. Adds one rule to the frozen version 1 protocol (`PROTOCOL.md`); everything else in version 1 stands.

1. **Unchanged.** Stages 01–05 (fold reference realignment, model fit, query addition, EPA-ng placement, placement summary) and the version 1 evaluation (`06_evaluation`) are run or reused exactly as specified in `PROTOCOL.md`, with the same stage keys. Version 1 outputs are never modified.

2. **Added rule N6 (culture ITS).** Applied after the unchanged decision tree, to both the full-V18 control and the external holdout of every fold, in a separate stage `06_evaluation_v2`:
   - Threshold T, per fold, from training cultures only. For every main-type copy of a named species with at least two training cultures, take the full-ITS distance to the nearest main-type copy of another training culture of the same species. T is the 99th percentile of those distances.
   - Query statistic: the full-ITS distance to the nearest copy, main-type or divergent-class, of any training culture.
   - Distance: edlib infix edit distance as a percentage of the shorter sequence (the project's standard).
   - A `known species` or `known species group` call becomes `novel species, genus known`, placed at the call's genera, when the query statistic exceeds T. No other status, placement or threshold changes.
   - The rule reads only the culture reference. It uses no external ITS database.

3. **Development folds.** The rule and its 99th percentile were chosen after inspecting the training-only culture-bridge benchmark, and after applying the rule post hoc to the two completed pilot folds (*Glomus chinense*, *Glomus rugosae*). These two are development folds and are reported separately. The other species folds, and all genus folds, are the first version 2 outcomes.

4. **Reporting.** For every fold, report version 1 and version 2 calls side by side, for both the full-V18 control and the external holdout, with:
   - false-known rates by copy, and averaged over species and over cultures;
   - genus-call correctness and coverage;
   - novelty rank;
   - behaviour of divergent-class copies.

   Keep development folds, first-outcome folds, deleted and already-absent strata, and species and genus ranks apart. Known-species retention is not measured by these folds; it comes from the culture-holdout benchmark (745/755 under species clades plus N6).

5. **No retuning.** The threshold, percentile, comparator pool and every other rule stay fixed whatever the version 2 outcomes are. A change needs a version 3 protocol and a fresh evaluation.

6. **Limits.** As in version 1: VT and SH external databases stay fixed, the cohort has already been used for development, and independent validation needs new verified cultures.
