# Thresholds on long fragments (exploratory)

For L. Tedersoo's follow-up to II.2-3 (13 October 2026): "How would the entire Krüger fragment, taken together, work?
And how would the entire 4.5-kb fragment work here?" It fixes no rule and changes nothing in the locked protocol or
its addendum.

Built by `stage3/copy_long_fragments.py` from culture reference v4 (main-type copies only), with the distances,
classes and within-species spread of `copy_thresholds.py`:

- **Krüger fragment:** from the SSUmAf1 site (Krüger et al. 2009) to the end of the copy at LSUmBr; median
  1577 nt (1451-1678). All 903 copies carry the site with at most one mismatch. The copies stop
  at LSUmBr, so the last stretch of the SSUmAf-LSUmAr amplicon is missing.
- **Whole copy:** AML1 to LSUmBr, primers trimmed; median 2758 nt (2630-2856). This is the longest linked
  fragment the cultures give. It is not the 4.5-kb fragment, which only long reads that cover it (such as
  EUKARYOME's) can test, and those carry no culture, so they show how close species come but not how far one
  genome spreads.
- **Within-species spread** (99th percentile, as in the paper): full ITS 5.99%, Krüger fragment
  3.30%, whole copy 1.93%.

For a distance cut-off t, a main-type copy is **split from its culture** if it lies more than t from every other
copy of its culture (839 copies in cultures with two or more main-type copies), **split from its species** if more
than t from every other copy of its species (803 copies of named species), and **joined to another species** if a
copy of another species of its genus lies within t (617 copies with a congener in the reference). The column
"outside the four complexes" leaves out the partners of the paper's complexes (`paper/analysis/v4.complexes.tsv`),
which no marker separates. "Best" is the cut-off that minimises the larger of the culture split and the join.

| Region | Cut-off | Identity | Split from its culture | Split from species | Joined to another species | ... outside the four complexes | Joined to another genus |
|---|---|---|---|---|---|---|---|
| Full ITS | 98.5% identity | 98.50% | 24.4% | 13.7% | 1.8% | 0.0% | 0.0% |
| Full ITS | 97% identity | 97.00% | 10.6% | 3.4% | 6.5% | 0.0% | 0.0% |
| Full ITS | within-species spread | 94.01% | 2.7% | 0.5% | 39.5% | 11.8% | 0.0% |
| Full ITS | best single threshold | 96.55% | 8.7% | 2.1% | 7.8% | 0.0% | 0.0% |
| Full ITS | best single threshold, complexes aside | 94.95% | 3.6% | 0.7% | 32.4% | 3.1% | 0.0% |
| Krüger fragment | 98.5% identity | 98.50% | 9.3% | 2.9% | 1.1% | 0.0% | 0.0% |
| Krüger fragment | 97% identity | 97.00% | 4.1% | 0.6% | 19.1% | 0.3% | 0.0% |
| Krüger fragment | within-species spread | 96.70% | 3.3% | 0.6% | 25.6% | 2.1% | 0.0% |
| Krüger fragment | best single threshold | 98.15% | 6.7% | 1.7% | 7.0% | 0.0% | 0.0% |
| Krüger fragment | best single threshold, complexes aside | 96.55% | 3.0% | 0.5% | 32.1% | 2.6% | 0.0% |
| Whole copy | 98.5% identity | 98.50% | 5.1% | 1.2% | 15.6% | 0.0% | 0.0% |
| Whole copy | 97% identity | 97.00% | 1.1% | 0.1% | 77.1% | 36.0% | 0.0% |
| Whole copy | within-species spread | 98.07% | 3.9% | 0.6% | 22.2% | 1.9% | 0.0% |
| Whole copy | best single threshold | 98.90% | 7.0% | 2.0% | 7.0% | 0.0% | 0.0% |
| Whole copy | best single threshold, complexes aside | 97.90% | 3.3% | 0.4% | 30.0% | 2.9% | 0.0% |

Files: `long.limits.tsv` (lengths and spreads), `long.nearest.tsv` (nearest copy per class and region),
`long.quantiles.tsv`, `long.sweep.tsv` (every cut-off from 0 to 12% in steps of 0.05), `long.summary.tsv` (the
table above) and `fig_thresholds_long.{pdf,png}` (split and join against the cut-off, per region).

  PYTHONPATH=<benchmark repo>/stage3 python3 stage3/copy_long_fragments.py --ref reference/culture_ref_v4 \
      --complexes paper/analysis/v4.complexes.tsv --out benchmarks/copy_checks/thresholds_long
