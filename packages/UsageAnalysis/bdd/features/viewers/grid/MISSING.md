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
- **No tooltip over a cell** (Tooltip > None, `grid.md` "Tooltip Settings"). The tooltip steps read
  what is visible; "none shown" after a hover passes before the 200 ms tooltip debounce
  (`_initCellTooltip`) has run. Needs a signal that the grid's tooltip request was answered
  (`GRID_SHOW_TOOLTIP` handled / declined), so the negative is read once the grid is quiet.
- **Form designer field positions** (`grid-ui.md` "Summary Columns — Form Designer"): a field box
  dragged in the designer "settles where it was dropped" and comes back on Edit. The designer's field
  boxes carry no names and report no positions (`tile-viewer.ts` reads only which fields exist).
  Needs a name per field box (`form-field-<column>`) and its position in the designer's status.
- **The hamburger popup's categorical filter** (`grid-ui.md` "Context Panel — Column Hamburger Menu
  Inline Filter", steps 4-5): `grid-columns.feature` drives the histogram filter of the popup; for a
  string column the popup holds a categorical filter widget whose rows have no hit areas of their own
  (the filter panel's `category <label> of <caption>` areas are the panel's, not the popup's). Needs
  the popup's filter to publish `category <label>` areas like the panel's card does.
- **The Order or Hide Columns dialog's list** reports its rows as `cell <r> of __name` and
  `cell <r> of x` (`grid-persistence.feature`): the internal column names of the dialog's grid, not
  "Name" and "Visible". A feature reading them depends on a private name; the dialog's grid should
  name its columns as the person sees them.
