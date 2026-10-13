# Supplementary Data S1: alignments and trees of the culture copies

One folder per order. For each region (SSU flank, full ITS, LSU flank) and the linked copy:

- `<region>.alignment.fasta`: MAFFT L-INS-i alignment of the unique sequences of the order.
- `<region>.tree.nwk`: IQ-TREE maximum-likelihood tree (ModelFinder; branch labels SH-aLRT/UFBoot from 1,000 replicates each). Tip label: `accession|species|culture|copy class|n copies`, where the accession is the first GenBank accession carrying that exact sequence.
- `<region>.members.tsv`: every copy each tip stands for (identical sequences share one tip).
- `linked.partitions.nex`: the region partitions of the linked alignment.

Copy classes: main type, short5.8S, samelength (the paper's divergent classes). Built by `stage3/copy_trees.py` under `benchmarks/copy_checks/PROTOCOL_COPY_CHECKS.md`; relabelled by `stage3/copy_tree_outputs.py`.
