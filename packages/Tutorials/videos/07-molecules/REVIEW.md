# Review: 07-molecules (Molecules)

Fresh-eyes review (2026-11-22, ~21:40-21:50): feature lines, descriptions.md, and a contact sheet at the key lines and
results of each delivered clip. Every delivered video matched its latest guide take, so the frames come from the
takes. The m:ss in the descriptions.md headings were the older renders: all refreshed to the delivered lengths.
08-r-group-analysis.feature was never delivered (the maker's trellis bug) and is not reviewed here.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Molecules in a table | refilmed (0:36) after the coordinator's no-jumps scan: Browse closed and the toolbox hidden off camera right after the file opened (the take dissolved from the Browse tree straight into a full-width grid, where a person repeating it sees the toolbox). Those two `Given` steps are gone: the toolbox now stays on the left as opening the file leaves it, and the grid still shows the structures. Text fixed: the description says how the context panel is opened (the status bar's Context Panel toggle). The ribbon shows Bio next to Chem on smiles (known product bug). |
| 02 | Draw or edit a molecule | OK. A chemist may notice that row 1's other columns (RingCount 3, NOCount 4...) still describe the old molecule after it becomes aspirin: they are plain numbers from the file, not computed. |
| 03 | Find molecules that contain a fragment | OK. |
| 04 | Find similar and diverse molecules | OK. Clicking a card makes its row current but the similarity viewer keeps the old reference (by design); the line only claims the row becomes current. |
| 05 | Calculate properties for every molecule | OK. "Computed in your browser" checked: with the Chemical Engine on Crux the properties run on Chem's Crux web workers. The Chemical Properties menu item's tooltip reads "Dynamically discovers and executes tagged biochemical calculators" (developer wording, known). |
| 06 | Map the chemical space | OK, minor: the dot's tooltip arrives at the very end of the "Point at a dot" step (21.5 of 15.8-21.8 s) and stays up while the next line "Color the map by any column" shows, so the tooltip is read under the wrong line; its structure drawing is faint. Not refilmed (a later tooltip render, not a step problem). |
| 07 | Find activity cliffs | refilmed (0:42). Text: "Pick the molecules and the activity, then OK" showed no pick (the dialog fills them in): now "Check the molecule and activity columns, then OK"; "Each row of the new table is one pair, with both activities" was false (the table shows the activity difference, SALI and similarity): now "Each row is one pair: both molecules and how far their activities differ". Jump (coordinator's scan): the pairs grid was enlarged by an off-camera drag; now the drag is filmed with the line "Drag its top border up to see more pairs". Needs work (lead, data/product): the line the clip hovers is "the first free line" and cannot be chosen; this take's tooltip shows two molecules with the same activity (7.68 and 7.68), the earlier take's 6.66 and 6.47, right after "Lines join similar molecules whose activity differs a lot"; and the cliffs table's top three rows are one structure twice (similarity 1.00, SALI Infinity, activity difference 0.19-0.45). The first cliffs a chemist sees are not cliffs. Needs either a library step to hover a named pair's line, or a product/data fix (identical fingerprints ranked first by SALI; equal-activity pairs drawn). |
| 08 | Highlight a fragment in every molecule | OK. The sketcher opens against the right edge of the window (fits, but tight). |
| 09 | What the Chem menu does | OK. |

Verdicts: 7 OK (06 with a minor late tooltip), 2 refilmed (01 with a description fix; 07 with two lines fixed and
the drag filmed), of which 07 still needs work: its tooltip and table examples contradict the clip's claim
(lead/product).

## Titles and descriptions

Plain and specific, written for a chemist; each names the menu path and the demo file. No hype, no questions.

## Jumps

None left. The coordinator's scan flagged 01 (Browse and toolbox closed off camera after opening the file) and 07
(the pairs grid enlarged off camera): both refilmed, 01 keeps the panels as they come, 07 films the drag. 07 still
hides the tooltip off camera before the next click, which moves nothing. All starts are off camera (smiles or
activity_cliffs opened, the grid narrowed); every new viewer, dialog and table appears on a filmed menu pick or click.

## Systematic

- Several first lines are bare menu paths ("Chem > Search > Substructure Search", "Chem > Calculate > Chemical
  Properties", "Chemistry > Highlight > Sketch"). They say where, not what for; with nested Chem menus the path helps
  a viewer follow the flyouts, so they were left. If the lead wants them changed, re-render from the take (edit the
  caption in the take's steps.json and the feature) without refilming.
- Hover tooltips of the scatter plot on molecules (06, 07) render late; a hover step on such a plot would need a
  longer rest after the tooltip appears.

## Merge round (2026-11-22, ~23:15-23:35)

The lead: clips under ~35 s are too thin. Two changes, each filmed, reviewed by frames and delivered:

| NN | Clip | What changed |
|---|---|---|
| 01 | View and edit molecules in a table (1:00) | 01 "Molecules in a table" (0:36) and 02 "Draw or edit a molecule" (0:28) in one clip, at the 1600x900 window of 01: open from Browse, Properties in the context panel, another molecule, then the sketcher on row 1, SMILES, OK. The panel follows the edit: it shows aspirin's C9H8O4 / 180.16 after OK (claimed). A hover on the column header and on a structure was tried and dropped: neither shows a tooltip. |
| 05 | Calculate properties for every molecule (0:57) | Extended from 0:29: two lines on what MW, LogP and PSA are and the Rule of Five limits, then the Search box (Ctrl+F, opens in the toolbox) "MW > 500" + Enter keeps 211 of 1000 (claimed exactly, and the status bar's "Filtered: 211" is pointed at), then a header double-click sorts by MW, heaviest first (1328.04 on top). Filtering before sorting: sorted first, the filter left the top rows unchanged on screen. |

Absorbed: 02-draw-or-edit-a-molecule.feature (and its generated spec), "02 - Draw or edit a molecule.mp4"; the old
"01 - Molecules in a table.mp4" is superseded by "01 - View and edit molecules in a table.mp4".

## Merge round 2 (2026-11-22, ~23:35-23:50)

The 31-37 s clips still left, merged in pairs; each filmed, reviewed by frames, delivered:

| NN | Clip | What changed |
|---|---|---|
| 03 | Find and highlight a substructure (1:12) | 03 "Find molecules that contain a fragment" (0:31) and 08 "Highlight a fragment in every molecule" (0:37), at 1280x720, grid narrowed to molregno + structures: Substructure Search C(=O)O keeps 314 of 1000 (claimed; the match drawn in red), the card's "Click to edit", the filter panel's reset icon brings 1000 back (filmed), then header click, Context Panel toggle, Chemistry > Highlight > Sketch with pyridine c1ccncc1: rows 1 and 3 show the ring in green, Rows: 1,000. A different fragment than the filter's, so the two results are not confused. Minor: once the context panel opens, the grid overlaps the filter panel's "Sketch" link ("Ske"), platform layout. |
| 04 | Find similar and diverse molecules, and map the chemical space (1:05) | 04 (0:33) and 06 "Map the chemical space" (0:31): Similarity Search, follow a row, click a card, Diversity Search, then Chem Space docks under the grid (four panels at 1080p, all readable); the pointer on a dot shows the structure tooltip, a click on row 15's dot makes it current and the similarity viewer follows (claimed by the scores reading), color by NumAromaticCarbocycles. |

Absorbed: 09-highlight-a-substructure.feature and 06-chemical-space.feature (and their generated specs);
"08 - Highlight a fragment in every molecule.mp4", "06 - Map the chemical space.mp4"; superseded by the new names:
"03 - Find molecules that contain a fragment.mp4", "04 - Find similar and diverse molecules.mp4".
