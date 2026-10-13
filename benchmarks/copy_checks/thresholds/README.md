# Distances at each level of relatedness (exploratory)

Exploratory analysis for L. Tedersoo's comments II.2 and II.3 (similarity thresholds; ITS vs SSU vs LSU
divergence). It fixes no rule and changes nothing in the locked protocol or its addendum.

Built by `stage3/copy_thresholds.py` from culture reference v4 (main-type copies only):

- `thresholds.nearest.tsv`: for every main-type copy and region, the distance to its nearest main-type copy in
  four classes (same culture; same species, other culture; same genus, other species; other genus, same order).
- `thresholds.quantiles.tsv`: median, 95th and 99th percentile per region and class, and the share of copies
  within the within-species spread (`below_spread`) and within the database threshold (`below_db`).
- `thresholds.pairs.tsv`: distances between the medoid copies (by full ITS) of every pair of cultures in an order.
- `thresholds.limits.tsv`: the within-species spread per region (99th percentile, computed as in
  `species_groups.py`; the script checks the three paper markers against the paper's limits).
- `fig_thresholds_nearest.*` (Figure S1) and `fig_thresholds_regions.*` (Figure S2).

Regions: the SSU NS31-AM1 fragment (on which MaarjAM VTs are defined; the Paraglomerales copies lack it), the whole
SSU flank, full ITS, and LSU LROR-FLR2. Database thresholds are drawn only on the fragments they apply to: VT 97% on
NS31-AM1, SH 98.5% on full ITS, V18 97% on LROR-FLR2.

Revision of 12 October 2026 (evening): the first version drew the VT 97% line on the whole SSU flank, which varies
less than NS31-AM1; the NS31-AM1 fragment was added as its own region and the VT line moved there. The
within-species spreads are now computed exactly as in `species_groups.py` (0.8065, 5.988, 3.027) instead of the
rounded values 0.8, 6.0, 3.0; no reported share changed.
