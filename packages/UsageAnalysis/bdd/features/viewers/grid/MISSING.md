# Grid — what the features cannot say yet

Written during the Grid / Filter panel / Legend gap round of 23 Oct 2026, against the grid's
"Automation surface" (`core/client/d4/lib/src/viewers/grid/CLAUDE.md`), the TestTrack cases
`files/TestTrack/Viewers/Grid/grid-ui.md`, `grid.md`, `grid-rows-select-filter-navigate.md` and the
library's rules (`public/libraries/bdd/CLAUDE.md`). Added in that round: `grid-pin-if-editable.feature`,
`grid-extracted-rows-project.feature` (GROK-19717), and the exact row count in `grid-columns.feature`
(326 after the popup's min handle reaches bin 10). Control+F and the block drag across the cells stay
out of scope, as `grid-rows.feature` records.

## Still open

- **The header histogram and the stats rows** (`grid-ui.md` "Column Stats — Visual Verification").
  The grid reports `stats rows` (which rows are drawn), not what they show: needs a
  `text of <stat> of <col>` reading (the string the stats row paints, `grid_special_rows.dart`) to
  claim "the Min row shows the minimum AGE". Same for the header histogram (`showColumnHistograms`):
  nothing says it was drawn for a column; a `header histogram <col>` area would.
- **Form designer field positions** (`grid-ui.md` "Summary Columns — Form Designer"): a field box
  dragged in the designer "settles where it was dropped" and comes back on Edit. The designer's field
  boxes carry no names and report no positions (`tile-viewer.ts` reads only which fields exist).
  Needs a name per field box (`form-field-<column>`) and its position in the designer's status.
- **The hamburger popup's categorical filter** (`grid-ui.md` "Context Panel — Column Hamburger Menu
  Inline Filter", steps 4-5): `grid-columns.feature` drives the histogram filter of the popup; for a
  string column the popup holds a categorical filter widget whose rows have no hit areas of their own
  (the filter panel's `category <label> of <caption>` areas are the panel's, not the popup's). Needs
  the popup's filter to publish `category <label>` areas like the panel's card does.

## Closed issues the regression round of 26 Oct 2026 could not claim

- **A cut-off cell's whole value in its tooltip (public#989, public#996)** regressed: the tooltip is empty;
  the steps and a video are on core#441. Claim it in `grid-viewer.feature` once fixed.

## Closed bugs round of 5 Nov 2026

Claimed in this round: public#3907 (`grid-column-groups.feature`, the last scenario: a member hidden before
its group is hidden and shown again stays hidden), core#417 (`grid-columns.feature`, the Sort dialog: its
column inputs offer a column renamed while it is open), public#3452 (`grid-pinned-molecules.feature`: every
row of a pinned SMILES column draws its molecule), public#3856 (`grid-appearance.feature`: a header colour
set over conditional coding shows next to it) and public#1812 (`filter-panel/filter-panel-linked-tables.feature`:
a column pinned in the target of a selection-to-filter link keeps its header and follows the filter).

- **Column menu icon and settings icon overlap (public#3787).** The grid places the column menu icon clear of
  the settings icon (`GridCore._columnMenuIconBounds`, `grid_core.dart:310`), but neither icon is a hit area
  of the grid's status, so a feature cannot say the two boxes do not intersect. Needs `column menu icon <col>`
  and `settings icon` hit areas, and a step `the {string} area of {widget} should not overlap the {string}
  area`.
- **A pinned column with an image renderer paints its cells (public#2691).** The pinned-molecule claim above
  covers a JS renderer in a pinned column; the curve (image) renderer of the issue needs a curves dataset in
  UsageAnalysis and the Curves package on the stand. Claim it in Curves' own features with
  `the "cell <r> of <col>" area of grid should be painted` after Pin > Pin Column.
- **Unselecting every source row empties a pinned column of the linked target (public#1812, part 2).** With
  no source row selected the selection-to-filter link stops filtering (`filter-panel-linked-tables.feature`:
  175 rows), so the empty state the issue describes is not reachable on this pair; the part that is
  (pinned header present, column follows the filter) is claimed.
- **A linked table shown inside the grid draws its structures in its own sub-column (core#369).** Needs
  `SPGI_v2_infinity` (12 MB) and a hit area per nested sub-column cell (`cell <r> of <linked table>/<col>`);
  the grid reports nested linked-table blocks as one column. Claim it once the nested cells are hit areas,
  on spgi-100 and SPGI-linked1 if the nested block shows for them.
