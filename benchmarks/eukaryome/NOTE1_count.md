# EUKARYOME long-read test, note 1: a count in section 2

Written 15 October 2026, during the run and before any score of it was seen. Recorded by `LOCK_NOTE1.sha256`.
It changes no rule, role, prediction, condition or outcome of `PROTOCOL_EUK.md`.

## The correction

Section 2 ("One rank up") says that 50 *Nanoglomus* reads lie inside a v4 culture named only to genus,
*Glomus* sp. 5089_2351A (42) or *Glomus* sp. 4406 (8). The two per-culture counts are right, but the 8 reads
inside *Glomus* sp. 4406 also lie inside 5089_2351A, so they were counted twice. The number of distinct reads is
**42**: 34 inside 5089_2351A only and 8 inside both (`genus_inside.reads.tsv`).

The decision of section 2, that *Nanoglomus* has no genus outcome, rests on there being any such read, not on
their number, and is unchanged.

## Code

`stage3/eukaryome_overlap.py genus_inside` printed per-culture counts only; it now also prints and writes the
number of distinct reads per genus. Its read table, from which the count above is taken, is unchanged.
