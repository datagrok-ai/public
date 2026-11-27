# Sequences

Playlist: **Sequences**. Features: `public/packages/Tutorials/bdd/features/howto/09-sequences/` (feature name
`How-to · Sequences`). For a biologist or peptide chemist new to
Datagrok: one idea per clip, everything computed in the browser (no Python, containers or outside services).

## Plan

| NN | Clip | One idea | Data |
|---|---|---|---|
| 01 | Open, inspect and convert sequences | A FASTA file from Browse becomes a Macromolecule column, drawn residue by residue; a residue's tooltip; the context panel for the column and for one sequence; Bio > Transform > Convert Sequence Notation to separator text (merged on 2026-11-22 from "Open sequences in a table" and "Convert between notations", 0:29) | Bio samples: FASTA.fasta |
| 03 | Align sequences and split them into positions | Bio > Analyze > MSA (Kalign, in the browser): motifs fall into columns; the header's conservation and logo; a click compares every row with one peptide; then Bio > Transform > Split to Monomers: one column per position, sort by a position (merged on 2026-11-22 from "Align sequences", 0:43, and "Split aligned sequences into positions", 0:32) | FASTA_PT |
| 04 | Composition with a WebLogo | Bio > Analyze > Composition: a logo per position; a letter's count; a click selects the peptides with that residue there | FASTA_PT_activity |
| 05 | Search sequences by motif, similarity and diversity | Bio > Search > Subsequence Search (a filter on the sequence column), Similarity Search and Diversity Search (merged on 2026-11-22 from "Find sequences by a motif", 0:24, and "Similar and diverse sequences") | FASTA_PT_activity |
| 07 | Map the sequence space and find the activity cliffs | Bio > Analyze > Sequence Space: similar peptides land together; then Bio > Analyze > Activity Cliffs: similar peptides far apart in activity, the pairs listed (merged on 2026-11-22 from "Map the sequence space", 0:34, and "Find activity cliffs", 0:38) | FASTA_PT_activity |
| 09 | Sequences to molecules | Bio > Transform > To Atomic Level; each monomer's atoms lit; then Chem measures the molecules | HELM sample |
| 10 | Monomer libraries | Bio > Manage > Monomers and Monomer Libraries: what a monomer is, where the libraries live | HELM sample |
| 11 | Peptide SAR | The Peptides package's SAR: residues by position, most potent residues, clusters; a residue selects its peptides | Demo peptides.csv (200 rows) |

Added after the plan:

| NN | Clip | One idea | Data |
|---|---|---|---|
| 12 | Number antibody chains | Bio > Annotate > Apply Numbering Scheme (IMGT, in the browser): regions on every chain, an aligned column, CDR3 extracted | Bio antibodies sample |
| 13 | Scan sequences for liabilities | Bio > Annotate > Scan Liabilities on the CDR3s: every hit marked, named under the pointer | Bio antibodies sample |

## 01 - Open, inspect and convert sequences (1:06)

Title: Open sequences in a table and convert their notation in Datagrok

Description: Datagrok recognizes sequences when it opens a file. This video opens Bio's sample FASTA file from
Browse (Files > App Data > Bio > samples): the sequence column is drawn residue by residue, each residue in its
color, and the pointer on a residue shows its name and structure. The context panel shows what Datagrok knows about
the column (FASTA notation, a peptide alphabet) and, for one sequence, how often each residue occurs. Bio >
Transform > Convert Sequence Notation then writes the column in another notation as a new column, here separator
text with dashes between the residues; the same dialog also writes HELM and BILN.

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/01-sequences-in-a-table.feature`

## 02 - Align sequences and split them into positions (1:14)

Title: Align sequences and split them into positions in Datagrok

Description: Bio > Analyze > MSA aligns a column of sequences right in your browser, using Kalign. In this video,
100 related peptides of different lengths, whose shared motifs sit at different places, are aligned into a new
column where the motifs fall into columns. The column header shows how conserved each position is and which
residues occur there, and a click on a peptide makes every other row show only what differs from it. Bio >
Transform > Split to Monomers, run on the aligned column, then gives every position a column of its own: position 8
is M in every peptide, and a double-click on a position's header sorts the peptides by the residue there. The data
is Bio's FASTA_PT sample (Files > App Data > Bio > samples).

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/02-align-sequences.feature`

## 03 - See the composition with a WebLogo (0:40)

Title: See sequence composition with a WebLogo in Datagrok

Description: Bio > Analyze > Composition adds a WebLogo: for every position of the sequences, the residues found
there, the most common on top. The pointer on a letter shows how many sequences have it, and a click selects them,
so every other viewer highlights the same rows. In this video, the peptides with N at position 3 all sit near the
top of the activity histogram, and nearly all of those with V at position 14 at its low end. The data is Bio's
FASTA_PT_activity sample (99 peptides with an activity).

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/03-composition-with-a-weblogo.feature`

## 04 - Search sequences by motif, similarity and diversity (1:04)

Title: Search sequences in Datagrok: by motif, by similarity, for diversity

Description: The three searches of the Bio > Search menu. Subsequence Search adds a filter for the sequence column
to the filter panel: type a motif and only the sequences that contain it stay, here 10 of 99 peptides, all from one
family; clear the box and all of them are back. Similarity Search lists the nearest neighbours of the current
sequence under the table, with the residues that differ in color, and follows as you click other rows. Diversity
Search lists the ten sequences most unlike each other, from all three families; clicking one finds it in the
table. The data is Bio's FASTA_PT_activity sample (99 peptides in three families).

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/04-search-sequences.feature`

## 05 - Map the sequence space and find the activity cliffs (1:13)

Title: Sequence space and activity cliffs for peptides in Datagrok

Description: Bio > Analyze > Sequence Space places every sequence on a scatter plot so that similar ones land close
together, and colors the groups it finds; this video compares sequences by Monomer chemical distance, which counts
how alike the monomers are as molecules. The pointer on a dot shows its peptide, and a click finds it in the table.
Bio > Analyze > Activity Cliffs then looks for pairs of nearly identical sequences with very different activity:
it colors the peptides by activity, links the pairs that are alike and makes the dots on steep cliffs bigger. Its
cliffs button lists the 21 pairs with the residue that differs, and a click on a pair zooms the plot to it and shows
both peptides in the table, here one residue apart and 3.84 against -2.81 in activity. The data is Bio's
FASTA_PT_activity sample (99 peptides in three families); everything is computed in the browser.

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/05-map-the-sequence-space.feature`

## 06 - Turn sequences into molecules (1:02)

Title: Turn peptide sequences into molecules in Datagrok

Description: Bio > Transform > To Atomic Level draws every sequence in a column as a full molecule, right in the
browser; dragging the new column's edge draws the molecules larger. With Highlight monomers on, each atom is
colored by the monomer it comes from, and the pointer on a
monomer of the sequence lights up its atoms. The result is an ordinary molecule column, so the chemistry tools work
on it: Chem > Calculate > Chemical Properties adds the molecular weight and LogP of every peptide. The data is
Bio's HELM sample of peptides with non-natural monomers.

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/06-sequences-to-molecules.feature`

## 07 - Monomers and monomer libraries (0:42)

Title: Monomers and monomer libraries in Datagrok

Description: Every sequence in Datagrok is built from monomers, and monomer libraries hold what each one is. Bio >
Manage > Monomers shows a library as a table of structures, symbols, names and attachment points, next to an editor
that follows the row you click; a click on an attachment point such as R3 shows the cap it gets when nothing bonds
there. Bio > Manage > Monomer Libraries lists the library files in use and, where two files define the same symbol
differently, lets you choose which one counts.

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/07-monomer-libraries.feature`

## 08 - Analyze peptide SAR (0:44)

Title: Peptide SAR analysis in Datagrok

Description: Bio > Analyze > SAR, from the Peptides package, shows which residue at which position drives the
activity of a peptide series. The dialog picks the sequence and activity columns; this video turns IC50 into pIC50
with the -lg scaling and sets how alike peptides must be to cluster. The result lays out a map of every residue at
every position (circle size: peptides in activity cliffs, color: how the activity changes), the most potent
residues, and the clusters with their logos; a click on a residue selects its peptides and lists them. The data is
the demo file Files > Demo > bio > peptides.csv.

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/08-peptide-sar.feature`

## 09 - Number antibody chains (1:02)

Title: Number antibody chains in IMGT and extract CDR3 in Datagrok

Description: Bio > Annotate > Apply Numbering Scheme numbers antibody chains in your browser, in IMGT, Kabat,
Chothia or another scheme. The frameworks and CDRs are then colored on every chain, the pointer on a residue names
its region, and a new column holds the chains aligned in the scheme's positions, with the regions drawn in its
header; the mouse wheel with Shift held over that header scrolls it to CDR3. A right-click on a chain extracts any
region, here CDR3, as a column of its own that keeps its IMGT positions. The data is Bio's antibodies sample (Files > App Data > Bio > samples > antibodies.csv).

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/09-number-antibody-chains.feature`

## 10 - Scan sequences for liabilities (0:43)

Title: Scan antibody sequences for liabilities in Datagrok

Description: Bio > Annotate > Scan Liabilities looks for sequence motifs that can cause trouble in development,
such as deamidation, isomerization, oxidation and N-glycosylation sites, and marks every hit on the sequences. In
this video the CDR3 of each heavy chain is extracted first from the right-click menu of a numbered chain, then
scanned; the pointer on a marked residue names the liability. Everything runs in the browser. The data is Bio's
antibodies sample.

Feature: `public/packages/Tutorials/bdd/features/howto/09-sequences/10-scan-liabilities.feature`

