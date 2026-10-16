#!/usr/bin/env python3
"""Presentation layer for the generated numbers and tables (version 19, after L. Tedersoo's comments).

The analysis scripts write names as deposited and the terms of earlier drafts. This script rewrites the generated
files in place so that the paper reads consistently:
  * genus names: species deposited in Rhizophagus are given in Rhizoglomus (Sieverding et al. 2015);
  * terms: 'sequence' for a deposited read ('copy' only for rRNA gene copies of a genome), 'assignment' for 'call',
    'matched' for 'carried', 'within-species limit' for 'within-species spread', 'SSU/LSU region' for 'flank',
    'near-best names' for 'candidate names', 'rRNA gene' for 'rDNA', and the primer LR0R (zero) for 'LROR';
  * a few captions and notes rewritten as plain titles (REWRITES, keyed on the text the scripts write).
It changes no number. It is idempotent, so it can be run after any numbers script:

  python3 analysis/display_names.py --paper .

Comment lines (starting with %) are left alone. A REWRITES entry whose old text is not found is reported, so a
change in a generator's wording shows up instead of passing silently.
"""
import argparse
import glob
import os
import re
import sys

NAMES = [  # (pattern, replacement), applied in order
    (r'\{\[\}Rhizoglomus\{\]\} vesiculiferum', 'Rhizoglomus vesiculiferum'),
    (r'\[Rhizoglomus\] vesiculiferum', 'Rhizoglomus vesiculiferum'),
    (r'\\spn\{\[Rhizoglomus\]\}', r'\\spn{Rhizoglomus}'),
    (r'Rhizophagus irregularis', 'Rhizoglomus irregulare'),
    (r'Rhizophagus intraradices', 'Rhizoglomus intraradices'),
    (r'Rhizophagus proliferus\b', 'Rhizoglomus proliferum'),
    (r'Rhizophagus prolifer\b', 'Rhizoglomus proliferum'),
    (r'Rhizophagus silesianus\b', 'Rhizoglomus silesianum'),
    (r'Rhizophagus vesiculifer\b', 'Rhizoglomus vesiculiferum'),
    (r'\birregularis\b', 'irregulare'),
    (r'\bprolifer\b', 'proliferum'),
    (r'\bsilesianus\b', 'silesianum'),
    (r'\bvesiculifer\b', 'vesiculiferum'),
    (r'\bRhizophagus(?! diaphanus)\b', 'Rhizoglomus'),
]

KEEP = ['recurrent deleted copies', 'recurrent deleted copy', 'a copy is transcribed', 'copy number']

TERMS = [
    (r'LROR', 'LR0R'),
    (r'\brDNA\b', 'rRNA gene'),
    (r'\bSSU flank\b', 'SSU region'), (r'\bLSU flank\b', 'LSU region'), (r'\bflank\b', 'region'),
    (r'\bLinked copy\b', 'Concatenated sequence'), (r'\blinked copy\b', 'concatenated sequence'),
    (r'\bWhole copy\b', 'Whole sequence'), (r'\bwhole copy\b', 'whole sequence'),
    (r'\bfiled as\b', 'classified as'),
    (r'\bCopies\b', 'Sequences'), (r'\bcopies\b', 'sequences'), (r'\bCopy\b', 'Sequence'), (r'\bcopy\b', 'sequence'),
    (r'\bCandidate names\b', 'Near-best names'), (r'\bcandidate names\b', 'near-best names'),
    (r'\bcandidate name\b', 'near-best name'),
    (r'\bno call\b', 'no name'),
    (r'\bCalls\b', 'Assignments'), (r'\bcalls\b', 'assignments'), (r'\bcall\b', 'assignment'),
    (r'\bcalled\b', 'assigned to'),
    (r'\bnot carried\b', 'unmatched'), (r'\bcarried\b', 'matched'), (r'\bcarries\b', 'has'),
    (r'(?<![A-Za-z\\])spreads\b', 'limits'), (r'(?<![A-Za-z\\])spread\b', 'limit'),
]

REWRITES = {
    'tables/table2_divergent.tex': [
        ('\\caption{The two divergent copy classes mislead each marker differently.}',
         '\\caption{Names given to the two classes of divergent sequences in the three systems, compared with the main '
         'type of their culture.}'),
    ],
    'tables/table4_benchmark.tex': [
        ('\\caption{Calls for withheld cultures and species, with the external databases intact.}',
         '\\caption{Assignment of the sequences of withheld cultures and species, with the external databases intact.}'),
        ('\\item \\emph{Known kept}: main-type copies of species held by another remaining culture that receive a known species or\n'
         'species-group call containing their species (148 rounds, one culture withheld in each). \\emph{False known}: known\n'
         'calls for copies whose species has no remaining culture: in 50 rounds with one species withheld in each, and, in the\n'
         'rounds with one culture withheld, for the main-type copies of species with a single culture. Averaged over species and\n'
         'over cultures, the final false-known rate with a species withheld is ',
         '\\item \\emph{Known kept}: main-type sequences of species represented by another remaining culture that were\n'
         'assigned to a known species or species group that includes their own species (148 rounds, one culture withheld in\n'
         'each). \\emph{False known}: sequences assigned to a known species although their species has no remaining culture,\n'
         'in 50 rounds with one species withheld in each and, in the rounds with one culture withheld, among the main-type\n'
         'sequences of species with a single culture. Averaged over species and over cultures, respectively, the final\n'
         'false-known rate with a species withheld is '),
    ],
    'tables/table5_external.tex': [
        ('\\caption{Withheld species also removed from the LSU reference.}',
         '\\caption{False-known assignments when the withheld species were also removed from the LSU reference.}'),
        ('\\item False-known calls among the main-type copies of the withheld species.\n'
         '$^{\\dagger}$\\,Development species: these two informed the ITS distance rule. \\emph{Intact}: copies placed on\n'
         'the full LSU reference (V18); \\emph{species removed}: V18 realigned without the species\' records, the model\n'
         'refitted and all copies placed again. The test species pool the species whose V18 records could be removed;',
         '\\item Numbers of main-type sequences of each withheld species that were assigned to a known species.\n'
         '\\emph{Intact}: sequences placed on the full LSU reference (V18); \\emph{species removed}: V18 realigned without\n'
         'the species\' records, the model refitted and all sequences placed again. $^{\\dagger}$\\,Development species, used\n'
         'to develop the ITS distance rule. The test species, i.e., all other species whose V18 records could be removed,\n'
         'are pooled;'),
    ],
    'tables/tableS1_m4_folds.tex': [
        ('\\caption{LSU-removal test by species: false-known calls among the main-type copies of each withheld species.}',
         '\\caption{False-known assignments for each withheld species in the LSU-removal test.}'),
        ('False-known calls with the LSU reference (V18) intact or with the species removed from it, each without and with the ITS distance rule.',
         'Main-type sequences of each withheld species that were assigned to a known species, with the LSU reference (V18) intact or with the species removed from it, each without and with the ITS distance rule.'),
        ('species not in V18 by label (removal changes nothing)', 'species not held by V18 under its name (removal changes nothing)'),
    ],
    'tables/table_paralog.tex': [
        ('\\caption{Paralog names: SHs and a VT reached by divergent copies and by at most one main-type copy.}',
         '\\caption{SHs and a VT that are reached by divergent sequences but by at most one main-type sequence (paralog '
         'names).}'),
        ('Name & Taxon of the best hit & Divergent copies & Main & Identity\\\\',
         'Taxon code & Taxon of the best hit & Divergent sequences & Main$^{a}$ & Identity$^{b}$\\\\'),
        ('\\item \\emph{Main}: main-type copies of any culture with that name. ',
         '\\item $^{a}$\\,Main-type sequences of any culture with that name. $^{b}$\\,BLAST identity of the divergent '
         'sequences to the best hit. Taxa are those of the best UNITE record (MaarjAM\'s for the VT). '),
        ('Taxa are those of the best UNITE record; the VT row is MaarjAM\'s.', ''),
    ],
    'tables/table_complexes.tex': [
        ('\\caption{Species complexes: pairs of named species whose nearest copies lie within the within-species spread on all three markers.}',
         '\\caption{Pairs of named species that form species complexes, i.e., whose closest sequences lie within the '
         'within-species limit on all three markers.}'),
    ],
    'tables/table3_sets.tex': [
        ('\\caption{Candidate names on the culture reference (version~4).}',
         '\\caption{Near-best names of the culture sequences in the three systems (culture reference version~4).}'),
    ],
    'tables/table2_edits.tex': [
        ('\\caption{Curation edits from the deposited copies to version~4 of the culture reference.}',
         '\\caption{Curation edits made to the deposited sequences to build version~4 of the culture reference.}'),
    ],
    'tables/table_stage1.tex': [
        ('\\caption{Barcode fragments cut from the 904 deposited rDNA copies (amfsplit v0.3.1; primer-trimmed fragments,\n'
         'primers found with up to 15\\% of their length in mismatches or gaps).}',
         '\\caption{Marker fragments cut from the 904 deposited rRNA gene sequences (primer-trimmed fragments; primers\n'
         'found with up to 15\\% of their length in mismatches or gaps; software version 0.3.1).}'),
    ],
    'tables/table_euk.tex': [
        ('\\caption{Field long reads: the rules as fixed, on \\EkTestN{} EUKARYOME reads.}',
         '\\caption{Assignment of the \\EkTestN{} EUKARYOME sequences with the rules as fixed in advance.}'),
        ('\\item With the ITS distance rule; in brackets, the same without it. \\emph{LSU reference}: \\emph{removed}, each absent species or genus that V18 holds removed in its own round, and intact for the reads of taxa it does not name; \\emph{intact}, every read placed on the full reference. \\emph{False known}: reads of taxa the reference does not hold called a known species or species group. \\emph{Absent species}: reads of species the reference does not hold, in genera it holds (\\spn{Dentiscutata erythropus}, whose genus is open, not counted). \\emph{Absent genera}: reads of genera the reference does not hold and of EUKARYOME placeholder lineages, so any genus given is wrong (\\spn{Nanoglomus} not counted, Section~\\ref{sec:euk-methods}). \\emph{Known kept}: reads of species the reference holds whose call contains their species. \\emph{Second}: VT coverage counted over the part of each MaarjAM reference that the read spans.',
         '\\item Counts with the ITS distance rule; in brackets, without it. \\emph{LSU reference}: \\emph{removed}, each absent species or genus that V18 holds was removed in its own round (sequences of taxa that V18 does not name were placed on the full reference); \\emph{intact}, every sequence was placed on the full reference. \\emph{False known}: sequences of taxa absent from the reference that were assigned to a known species or species group. \\emph{Absent species}: sequences of species absent from the reference, in genera present in it (\\spn{Dentiscutata erythropus}, whose genus is uncertain, not counted). \\emph{Absent genera}: sequences of genera absent from the reference and of EUKARYOME placeholder lineages, for which any genus given is wrong (\\spn{Nanoglomus} not counted; Section~\\ref{sec:euk-methods}). \\emph{Known kept}: sequences of species present in the reference whose assignment includes their own species. \\emph{Second}: VT coverage counted over the part of each MaarjAM reference that the sequence spans.'),
    ],
    'tables/tableS2_pcgenes.tex': [
        ('\\caption{Protein-coding genes on the rDNA species complexes:', '\\caption{Distances in protein-coding genes between the members of the species complexes:'),
    ],
}


def process(path, rel, report):
    text = open(path).read()
    for old, new in REWRITES.get(rel, []):
        if old in text:
            text = text.replace(old, new)
        elif new and new not in text:
            report.append(f'{rel}: rewrite not applied (old text not found): {old[:70]!r}')
    out = []
    for line in text.split('\n'):
        if line.lstrip().startswith('%'):
            out.append(line)
            continue
        for pat, rep in NAMES:
            line = re.sub(pat, rep, line)
        keep = {}
        for i, k in enumerate(KEEP):
            tok = f'\x00K{i}\x00'
            if k in line:
                line = line.replace(k, tok)
                keep[tok] = k
        for pat, rep in TERMS:
            line = re.sub(pat, rep, line)
        for tok, k in keep.items():
            line = line.replace(tok, k)
        out.append(line)
    new = '\n'.join(out)
    if new != open(path).read():
        open(path, 'w').write(new)
        return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--paper', default='.')
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.paper, 'numbers_*.tex')) + glob.glob(os.path.join(a.paper, 'tables', '*.tex')))
    report, changed = [], []
    for f in files:
        rel = os.path.relpath(f, a.paper)
        if process(f, rel, report):
            changed.append(rel)
    print(f'{len(changed)} of {len(files)} files changed' + (': ' + ', '.join(changed) if changed else ''))
    for r in report:
        print('  NOTE', r)
    return 0


if __name__ == '__main__':
    sys.exit(main())
