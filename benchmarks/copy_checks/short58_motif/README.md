# A sequence motif for the short-5.8S paralog (exploratory)

For L. Tedersoo's suggestion of a "microbarcode" that flags short-5.8S copies in databases (14 October 2026). It
fixes no rule and is not part of any locked protocol.

In *Rhizophagus* and [*Rhizoglomus*] the short type differs from the main type in one hairpin near the 3' end of
the 5.8S (main ...CACTCTCTGG(T/C)AACCCGGGGAGTATG..., short ...CACTCTCTT(T/C)GGAGAGTATG...). The motif
YTYGGRGARTATG (IUPAC) was written from these copies and is tested on the same copies, so the counts show that such a
flag can work, not its error rate on new data.

Built by `stage3/short58_motif.py` from `reference/culture_ref_v4` (903 retained copies). Short copies are those of
`stage3/copy_paralog_effects.py`: the short-5.8S class and every copy whose 5.8S is identical to one of it.

- `motif.copies.tsv`: every copy, its class, 5.8S length, whether its 5.8S carries the motif and where, and how many
  times the motif occurs elsewhere in the whole copy.
- `motif.summary.tsv`: copies and motif carriers per genus and class (*Rhizoglomus* counted with *Rhizophagus*).
- `primers.tsv`: start of the motif and of the ITS2 primers fITS7 (GTGARTCATCGAATCTTTG) and ITS3
  (GCATCGATGAAGAACGCAGC) in the 5.8S of the *Rhizophagus* copies, with up to two mismatches.

Result: all 41 short copies of *Rhizophagus* and [*Rhizoglomus*] carry it and none of their other 113 copies; it
occurs nowhere else in any copy. It misses the 3 short copies of *Septoglomus viscosum* (a different deletion), and
it is the normal 5.8S of *Oehlia diaphana* (29 of 29 copies) and is also in PX214775.1 (*Glomus* sp. 5089_2351A).
It starts at 5.8S position 121-122, downstream of fITS7 (84-85) and ITS3 (31).

  python3 stage3/short58_motif.py --ref reference/culture_ref_v4 --out benchmarks/copy_checks/short58_motif
