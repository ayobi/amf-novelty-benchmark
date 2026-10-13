# Thresholds on whole long reads: EUKARYOME Glomeromycota (exploratory)

For L. Tedersoo's follow-up to II.2-3 (13 October 2026): "how would the entire 4.5-kb fragment work here?". It fixes
no rule and is not part of any locked protocol.

Built by `stage3/longread_thresholds.py` (run in the Phase 0 workspace, on the windows of its Amendment 1) from the
EUKARYOME v2.1 long-read set (`General_EUK_longread_v2.1.fasta`, CC BY 4.0; SHA-256 536e3e10...). Reads: the
Glomeromycota reads whose species field gives a named species (named genus plus a lower-case epithet; no placeholder
marks), with at least two reads per species and all three windows cut. Each read is rebuilt in rRNA orientation as
SSU + ITS + LSU window. EUKARYOME v2.1 spells one species two ways; `Rhizoglomus_irregularis` was merged into
`Rhizoglomus_irregulare`. *Rhizoglomus* is read as *Rhizophagus* for the genus classes, as in the paper.

Run of 13 October 2026 (from its log): 347 reads of 43 species; left out 32 reads of species with a single read and
5 reads without the Kruger primer sites (Kruger region only). Lengths: full ITS median 483 nt (400-1,640); Kruger
region 342 reads, median 1,590 nt (218-2,730); whole read median 4,144 nt (1,059-6,976). Later runs write these to
`longread.limits.tsv`.

Regions: the full ITS window; the Kruger region, from the SSUmAf1 site to the LSUmBr site; the whole read. For a
cut-off t, a read is split from its species if no other read of its species lies within t, and joined to another
species if a read of another named species of its genus lies within t (edlib infix edit distance as % of the
shorter sequence, as in the paper). These are environmental reads labelled by EUKARYOME's classification, not
cultures: there is no genome-level view, and a label error counts as an overlap.

Files: `longread.summary.tsv` (database thresholds, the within-species spread and the best single cut-off per
region), `longread.sweep.tsv` (every cut-off from 0 to 12%), `longread.pairs.tsv` (the species pairs behind the
joins at 97%, 98.5% and the best cut-off), `longread.nearest.tsv` (nearest read per class and region),
`longread.species.tsv` (reads per species).

  python3 stage3/longread_thresholds.py --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta \
      --windows windows --out amf_longread --procs 4
