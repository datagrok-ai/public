# Review: 09-sequences (Sequences)

Fresh-eyes review (2026-11-22, ~23:20-23:40) of the 14 delivered clips: feature, descriptions.md block, and frames
of each delivered take at its key lines and results (every delivered video matched its latest take). Judged on:
is the one idea clear to someone new to Datagrok, is each line true when it shows (data claims checked against the
frame or the data), plain wording, jumps. Nothing was changed: no refilm, no text edit.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Open sequences in a table (0:45) | OK. Residues colored; Details reads FASTA and a peptide alphabet; Composition and the molecule pane show for the clicked sequence. |
| 02 | Convert between sequence notations (0:29) | OK. The new separator(HELM) column is narrow (the command's default width) but shows enough to see the dashes; the pointer on hHis shows the same monomer. |
| 03 | Align sequences (0:43) | OK. The msa(sequence) column is narrow and shows its first ~15 positions only, but there the aligned block (MTMFC) visibly lines up across rows, so "the shared motifs fall into columns" holds on screen. |
| 04 | See the composition with a WebLogo (0:40) | OK. N at 3 selects 31 peptides, all in the top bin of activity (about 4.4); V at 14 selects 34, in the lower and middle bins. The second line ("nearly all at the low end") shows about 2 s before the click that proves it: it sits on the hover, not the click. Minor. |
| 05 | Find sequences by a motif (0:24) | OK. NEMTK: Filtered 10, every cluster value 1. |
| 06 | Find similar and diverse sequences (0:41) | OK. Three stacked bands (the table, the similarity list, the diversity list) are the two commands' default docking; each is labelled and the table keeps 6 rows in view. "From all three families" is asserted by the feature (3 cluster values), not visible in the list itself. |
| 07 | Map the sequence space (0:34) | OK. The Embed_X_1 / Embed_Y_1 / Cluster (DBSCAN) columns are the command's output; three clusters on the plot. |
| 08 | Find activity cliffs (0:38) | OK. 21 CLIFFS on the button; the clicked pair is rows 70/71, 3.84 against -2.81 (delta 6.65, similarity 0.95). The `zoomAndFilter` Given after the first visible step changes nothing on screen when it runs; it changes what follows: with the default, the click's zoom would also filter the table to the zoomed points. The lines hold either way. |
| 09 | Turn sequences into molecules (1:03) | OK. Over a minute for two linked steps (to atomic level, then Chemical Properties on the result); every atom colored by its monomer; MW and LogP for all 5 peptides (2203.33, -4.38 ...). |
| 10 | Monomers and monomer libraries (0:42) | OK. |
| 11 | Analyze peptide SAR (0:44) | OK. The `user listens for "peptides-sar-ready"` Given after the first visible step is a test hook: nothing on screen changes. The result view is busy (five panels), and the lines name where to look (bottom left, top right); a residue click selects 16 peptides and lists them. |
| 12 | Number antibody chains (1:02) | OK. Over a minute: number, scroll to CDR3, extract CDR3; the last overlaps 13's first step but belongs here. "CDR3 alone, with its IMGT positions: 105 to 117" is IMGT's CDR3 range; the header shows 105, 110, 111A. |
| 13 | Scan sequences for liabilities (0:43) | OK. The liabilities are underlined on the CDR3 column. |
| 14 | Split aligned sequences into positions (0:32) | OK. Checked against Bio's MSA.csv: of the first 50 peptides, 23 have hHis at position 2 (26 have a gap there), and the sort puts them first. |

Verdicts: 14 OK.

## Data claims

All true. Counted or read: 04 (31 and 34 selected, where they sit on the histogram), 05 (10 of 99, one cluster), 06
(3 clusters, by the feature's check), 08 (21 pairs; 3.84 against -2.81), 09 (MW and LogP for every row), 11 (16
selected), 12 (IMGT 105-117), 14 (23 with hHis, from the file).

## Titles and descriptions

Plain and specific, each naming its menu path and sample file; no banned words, no questions.

## Systematic

- New result columns open at the commands' default widths (02, 03): narrow, but each clip's line is about the part
  that shows.
- Two off-camera Givens after the first visible step (08, 11) change nothing visible when they run; 08's changes
  later behaviour (no table filter on zoom), which a person repeating the clip would see differently.

## Merge round (2026-11-22, ~23:15-23:35)

The lead: clips under ~35 s are too thin. Two merges, each filmed, reviewed by frames and delivered:

| NN | Clip | What changed |
|---|---|---|
| 01 | Open, inspect and convert sequences (1:06) | 01 "Open sequences in a table" (0:45) and 02 "Convert between sequence notations" (0:29) on the same FASTA table: open, residue tooltip (Leucine), column Details, one sequence's Composition, then the panel closed on camera and Convert Sequence Notation to separator text; the separator column reads M-I-E-V-F-L-F-G-I-V-L-G-L-I-F-I for row 2 and the pointer on its 6th monomer names Leucine again. The Molecular Structure pane was dropped (unreadable at this size). Converting to HELM was filmed first and dropped: in a column this narrow HELM is drawn as a chain of dots and the rows grow taller. The drag on the column edge is 200 px (was 300) to leave room for the new column. |
| 05 | Search sequences by motif, similarity and diversity (1:04) | 05 "Find sequences by a motif" (0:24) and 06 "Find similar and diverse sequences" (0:41): NEMTK keeps 10 of 99, all one family; clearing the box (filmed) brings 99 back; then Similarity Search and Diversity Search as before, two diversity clicks instead of three. The filter panel stays on the left, empty, during the last two searches. |

Absorbed: 02-convert-between-notations.feature and 06-similar-and-diverse-sequences.feature (and their generated
specs), "02 - Convert between sequence notations.mp4", "06 - Find similar and diverse sequences.mp4"; superseded by
the new names: "01 - Open sequences in a table.mp4", "05 - Find sequences by a motif.mp4".

## Merge round 2 (2026-11-22, ~23:35-23:50)

| NN | Clip | What changed |
|---|---|---|
| 03 | Align sequences and split them into positions (1:14) | 03 "Align sequences" (0:43) unchanged, then Split to Monomers on msa(sequence) (picked in the dialog on camera): 44 position columns, the grid scrolls to them (platform), "Position 8 holds M in every peptide" (asserted: every value of column 8 is M; the header tooltip at 100% agrees), the Methionine tooltip, a double-click on position 14 sorts (W, N, M... grouped). Replaces 14 "Split aligned sequences into positions" (0:32, on Bio's MSA sample of non-natural monomers). The pointer-away step after the hover was dropped (it outlined the whole grid for 2.7 s); the tooltip now fades under the next step's outline. |
| 07 | Map the sequence space and find the activity cliffs (1:13) | 07 (0:34) and 08 "Find activity cliffs" (0:38) on one table: Sequence Space (Monomer chemical distance), the dot tooltip, a click finds row 15; then Activity Cliffs docks under the first map (both stay: same peptides, by cluster and by activity); 21 cliffs, the first pair rows 70/71, QKYLDCLWMDYRPVNH vs ...WGD..., 3.84 against -2.81 (read on the frame). Off camera after the first visible step, as before: the cliffs plot set not to filter by zoom. While the cliffs plot is added the grid briefly scrolls right (new columns), and is back by the pair click. |

Absorbed: 08-find-activity-cliffs.feature and 14-split-to-monomers.feature (and their generated specs);
"08 - Find activity cliffs.mp4", "14 - Split aligned sequences into positions.mp4"; superseded by the new names:
"03 - Align sequences.mp4", "07 - Map the sequence space.mp4".
