# LongGloDB independent-culture test, amendment 1

Written 6 October 2026, after `PROTOCOL_LGDB.md` (sha256 bf16ed1d…) and `predictions.lock`
(sha256 c1aaed69…) and before any test sequence was windowed, placed, scored or classified.
It changes no rule, threshold, role, condition or outcome of the method.

## 1. Two definitions of the 99th percentile

The predictions were computed with the within-species spread T by linear interpolation:
ITS 5.975%, LSU 2.913%, SSU 0.747%. The frozen method (`its_rule.py`, protocol v2) and the
paper's species complexes use the default of Python's `statistics.quantiles` (exclusive
method): ITS 5.988%, LSU 3.027%, SSU 0.807%. The paper states 6.0%, 3.0% and 0.8%. The
difference came from a copy of `species_groups.py` changed for a percentile sensitivity
analysis. It does not affect the decision tree, whose only percentile threshold (N6) is
computed inside the frozen `its_rule.py`. Section 2 of the protocol, which says T was
recomputed "by the same method", should read: same percentile, different interpolation.

## 2. Effect on the predictions

The predictions were recomputed under both definitions by an independent reimplementation,
which reproduces the locked own-species distances exactly. One novel-species ASV changes:
*Gigaspora decipiens* AU102_9 ASV1303 (ITS 5.94%, LSU 2.95%, SSU 0.44% to *G. rosea*) lies
inside *G. rosea* under the frozen definition and outside it under the locked one, because its
LSU distance falls between the two LSU thresholds. Nothing else in the outcome roles changes:
no novel-lineage ASV is inside a v4 species under either definition, and 139 of 164
known-species ASVs are inside their own species under both. Genus-only ASVs, which have no
outcome, change from 3 to 6 inside a v4 species.

## 3. Decision

The locked set (the 4 ASVs of *G. gigantea* KS302 inside *G. rosea*) stays the primary test of
the complex claim, because it is the stricter set and was fixed first. The frozen-definition
set (those 4 plus *G. decipiens* AU102_9 ASV1303) is reported alongside it. A false-known call
for AU102_9 ASV1303 therefore counts as a contradiction in the primary test, and the report
states that the paper's own definition would have predicted it.

## 4. Record

The SHA-256 of this file is recorded with those of `PROTOCOL_LGDB.md` and `predictions.lock`
before the run starts.
