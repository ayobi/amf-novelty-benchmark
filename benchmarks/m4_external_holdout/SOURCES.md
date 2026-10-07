# Sources used for exclusion names and software installation

Consulted 29 September 2026. These sources justify nomenclatural mappings and software versions, not classifier accuracy. Benchmark numbers are computed from the supplied frozen culture collection.

- Rhizophagus prolifer / proliferus: Index Fungorum name record (orthographic variant): https://www.indexfungorum.org/Names/NamesRecord.asp?RecordID=803201; NCBI record https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2650738 .
- Rhizoglomus vesiculiferum and the earlier Rhizophagus combination: Błaszkowski et al. 2018, Nova Hedwigia 107:501–518, https://doi.org/10.1127/nova_hedwigia/2018/0488; primary author-uploaded paper https://www.researchgate.net/publication/325174120_A_new_genus_Oehlia_with_Oehlia_diaphana_comb_nov_and_an_emended_description_of_Rhizoglomus_vesiculiferum_comb_nov_in_the_Glomeromycotina; name record https://www.indexfungorum.org/Names/NamesRecord.asp?RecordID=824697 . The frozen V18 header uses Rhizophagus vesiculiferus.
- Septoglomus viscosum / Viscospora viscosa: Index Fungorum https://www.indexfungorum.org/names/NamesRecord.asp?RecordID=518471 . V18 includes both genus strings in a composite header. This mapping removes the same named species without adopting a new project-wide genus classification.
- MAFFT 7.526 packages: https://anaconda.org/conda-forge/mafft/files?version=7.526 .
- EPA-ng 0.3.8 packages: https://bioconda.github.io/recipes/epa-ng/README.html .
- RAxML-NG 1.2.2 packages: https://bioconda.github.io/recipes/raxml-ng/README.html .
- Fixed-backbone source: https://github.com/c383d893/AMF-LSU-Database-and-Pipeline2, commit f4a0014336f49be6aea48e8374dae052360a88bd.

The accompanying environment specifies versions available from the published package channels. The pilot used the same versions as M2b, with explicit local executable paths. A conda environment solve on the user's PC has not been performed here; the runner checks actual executable versions before starting work.
