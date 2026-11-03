# Measurement E: locked prediction

Written and committed on 2 November 2026, before any round of measurement E beyond the
*Funneliformis mosseae* pilot has run. Protocol v3, section 5.

**Measurement.** E reruns the 37 species rounds of the M4 harness (release 0.3.2) with each test
species removed from all three external databases. V18 loses the species' records as in M4. MaarjAM
and UNITE lose every record whose nearest named species of culture reference v4 is the withheld one
and that lies within that species' limit (NS31-AM1 1.71795%, full ITS 5.98803%). A VT whose v4
holders all belong to the species goes whole. Records nearer a congener stay, because they would
exist even if the species had never been sequenced.

**Prediction.** With the ITS distance rule (N6), false-known calls among the 673 main-type query rows
of the 34 first-outcome test species (the paper's denominator) rise from 53 to at most 71.

**Basis.** 18 of those rows had a novelty call resting on an unmatched VT, so removing the species'
VT records can turn them known. Another 34 rest on names that share no species (N3) and could also
flip, which is why the bound can fail, and why it is a test.

**Pilot.** The *F. mosseae* round is E's first and counts as development. It has run in part: no
false known under the paper's arm, and the genus returned for all 130 copies.
