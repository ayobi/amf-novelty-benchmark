# The short-5.8S motif in UNITE (exploratory)

An independent test of the motif of `../README.md`, which was written from the culture copies. It fixes no rule.

Built by `stage3/short58_motif_db.py` from the UNITE public release of 19 February 2025
(`UNITE_public_19.02.2025.fasta.gz`, 2,069,189 sequences; SHA-256 abe98f43..., the file in
`external_refs/CHECKSUMS.sha256`). None of the culture copies (GenBank PX207096-PX215164) is in it. A motif hit
counts as in the 5.8S when fITS7 starts 30-45 nt before it, as in the culture copies (37 nt). SHs are the SH codes
in the UNITE records; sh.tsv counts only the sequences labelled *Rhizophagus* or *Rhizoglomus*.

Results (sequences labelled *Rhizophagus* or *Rhizoglomus* that cover the site, leaving out *Rhizophagus
diaphanus*, the old name of *Oehlia diaphana*, whose normal 5.8S carries the motif):

- 259 of 1,355 sequences (19%) carry the short-type hairpin.
- Of the 521 SHs holding these sequences, 111 hold only short-type sequences (249 sequences; 32 SHs with two or
  more), 407 only main-type sequences, and 3 hold both.
- Outside these genera the motif is in the 5.8S of 984 sequences, 940 of them in SHs with no *Rhizophagus*- or
  *Rhizoglomus*-labelled sequence: unplaced Glomeraceae 502, *Glomus* 265 (246 *Glomus* sp., 19 *G. invermaium*,
  now *Rhizophagus invermaius*), unplaced Glomeromycota 172, *Oehlia* (as *R. diaphanus*) and others. Whether these
  are unlabelled short-type *Rhizophagus* or other lineages with the same hairpin cannot be told from the motif.

Files: `genus.tsv`, `species.tsv`, `sh.tsv`, `hits.tsv` (accessions and positions only, no sequences),
`manifest.txt`.

  python3 stage3/short58_motif_db.py --fasta UNITE_public_19.02.2025.fasta.gz --out benchmarks/copy_checks/short58_motif/unite
