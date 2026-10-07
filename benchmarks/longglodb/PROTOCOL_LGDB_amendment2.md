# LongGloDB independent-culture test, amendment 2

Written 6 October 2026, after the Stage 1 windows of the 280 test ASVs were made
(`bench/longglodb/stage1`) and before any test sequence was searched against MaarjAM, UNITE or
V18, placed, scored or classified. Amendment 1 (sha256 f2672f7d…) stands.

## 1. The test amplicons begin inside the VT region

The test sequences start just downstream of the NS31Glo3 site, about 60 bp inside the NS31–AM1
region on which MaarjAM virtual taxa are defined. The NS31 site is absent from 275 of 280 ASVs
(window `open_5p`, median 446 bp against 508 bp in the reference copies), and the SSU flank has a
median of 1,180 bp against 1,473 bp. The frozen VT arm counts a MaarjAM hit only when it covers at
least 95% of the shorter sequence, which here is the MaarjAM reference (median 520 bp). A test SSU
flank can cover a median of about 88% of a MaarjAM reference, and at least 95% of about 1% of them.
The VT arm will therefore make no call for nearly all test ASVs.

## 2. Consequences under the frozen decision tree

- N1 (VT unmatched) cannot fire.
- VT drops out of the species and genus intersections; calls rest on SH and LSU.
- The VT route of D2 cannot fire; the SH route can.
- D1 uses its fallback comparators: cultures whose nearest copy is within 1% in the SSU flank.
- B1 labels an ASV "boundary" when its recorded best VT identity lies within one point below 97%.

## 3. Decision

- **Primary analysis: unchanged.** The frozen VT arm is applied as it is. This is what the method
  does with amplicons that begin at this primer, and it is the harsher test of the novelty claim,
  because the N1 demotion is unavailable.
- **Secondary analysis, fixed here.** The VT arm with coverage measured over the part of each MaarjAM
  reference that the query spans: reference bases upstream of the query's 5' end are discounted, as
  reference Ns already are. Identity (97%), e-value, merged units, δ and every other rule are
  unchanged. Reported apart, with the number of ASVs whose VT state changes.
- **Added to the report:** the number of test ASVs with status "boundary" (B1), with flag B2, and the
  VT state (matched, unmatched, silent) under each condition.

## 4. Reading of the lock

The locked table (`predictions.predictions.tsv`), not the printed log, is authoritative; the log
listed at most eight known-species ASVs outside their own species. The table records that all 8
ASVs of *Funneliformis coronatus* MR101 lie inside *F. mosseae* and not inside *F. coronatus*, and
that 17 known-species ASVs lie inside no v4 species (*Ambispora leptoticha* MX982A 8,
*Dentiscutata heterogama* IL203A 2 and FL728 1, *Acaulospora morrowiae* BR983A 2, *Gigaspora rosea*
BR155B 2, *F. mosseae* BEG12 1, *A. leptoticha* CR312 1). Nothing in the lock changes.

## 5. Record

The SHA-256 of this file is added to `protocol_hashes.txt` before the run starts.
