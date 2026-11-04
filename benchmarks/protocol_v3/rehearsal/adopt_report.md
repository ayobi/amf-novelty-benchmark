# Protocol v3, section 8 applied

Development stand-ins (LongGloDB v1, EUKARYOME v2.1): a rehearsal, which decides nothing. The power check runs at the set lock.

| Rule | Criterion | Measured | Verdict |
| --- | --- | --- | --- |
| R1A | added demotion of K, reference intact, bounded at 1 point(s); removes at least one targeted error | K kept 287 of 363 (frozen 288; difference -1.05 to +0.00 points); demotion bound 1.05 points; false known in A, L and withheld K 35 to 25 | fails |
| R1B | added demotion of K, reference intact, bounded at 1 point(s); removes at least one targeted error | K kept 286 of 363 (frozen 288; difference -2.10 to +0.00 points); demotion bound 2.10 points; false known in A, L and withheld K 35 to 25 | fails |
| R2 | added demotion of K, reference intact, bounded at 3 point(s); removes at least one targeted error | K kept 284 of 363 (frozen 288; difference -2.34 to -0.32 points); demotion bound 2.34 points; false known in A, L and withheld K 35 to 22 | passes |
| R4 | demotion of K falls; false known in A, L and withheld K rises at most 1 point, upper bound under 3 | K kept 288 of 363 (frozen 288; difference +0.00 to +0.00 points); false known 36 of 3171 (frozen 35; difference +0.00 to +0.13 points) | fails |
| R5 | at least 80% of Rhizoglomus-labelled carriers called divergent; at most 1% of other carriers reclassified | 72.7% (80 of 110); 0.23% (1 of 438) | fails |
| R5S | at least 80% of Rhizoglomus-labelled carriers called divergent; at most 1% of other carriers reclassified | 98.2% (108 of 110); 0.23% (1 of 438) | passes |
| R3A | at least 50% of wrong genera flagged in the genus rounds and strata A and L; at most 2% of right genera flagged | 61.6% (697 of 1132); 2.17% (9 of 415) | fails |
| R3A_ALL | at least 50% of wrong genera flagged in the genus rounds and strata A and L; at most 2% of right genera flagged | 1.7% (19 of 1132); 0.00% (0 of 415) | fails |
| R3B | at least 50% of wrong genera flagged in the genus rounds and strata A and L; at most 2% of right genera flagged | 21.1% (239 of 1132); 0.24% (1 of 415) | fails |
| R7 | with an adopted R3: the right genus for at least 50% of F's genus-less calls whose genus another culture holds; at least 80% of its wrong genera flagged | right genus 71 of 73; wrong genera flagged 18 of 25 (with r3b, no R3 adopted) | stays opt-in (no R3 adopted) |
| R8 | flag on by default; R8d adopted only if single-marker new lineages turn up among K with the reference intact | 0 of 363 K sequences | flag on; R8d not adopted |
| R6 | at most 1% of triage negatives called new lineage, upper bound under 2% | 2.50% (485 of 19402; upper bound 4.55%; 495 named under R6) | filtering stays (each run reports what triage dropped) |
| V | with an adopted R3: every wrong genus given through a VT alone flagged; no false known through a VT shared by species | VT-alone wrong genera flagged 5 of 46; new false knowns through a shared VT 2 (with r3b) | stays opt-in (no R3 adopted) |
| S | demotion of K within the adopted rules' largest ceiling | r5s, r2, r8: K kept 284 of 363 (frozen 288; difference -2.34 to -0.32 points); demotion bound 2.34 points | passes |
