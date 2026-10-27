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

## Edge cases (`grid-edge-cases.feature`)

Added with `grid-edge-cases.feature` (edge-case round of 2026-10-26, localhost). Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Steps requested: typing into the property panel.** `user sets "<caption>" property of grid to
   "<v>"` goes through the API, so it cannot say what the panel's editor accepts. Row Height typed as
   0 or 501 (the property declares 1..500) and Frozen Columns typed as -1 need a step that types into
   a property row of the context panel (`user types "<v>" into the "<caption>" property of grid`).
2. **Questions for owners**:
   - Frozen Columns above the column count or below 0: clamp, or refuse (core#1137 is the error only).
   - Add New Row On Last Row Edit with a rejected value in the last row: is a row added?
   - Allow Edit off: does Control+V show the read-only notice, and may Shift+Delete remove rows?
   - Where empty values go in a sort, ascending and descending.
   - Sort, pinned rows and column groups after a column is renamed (all are kept by name).
   - A group name of spaces only; Editable by with spaces only or another letter case.
   - Shift+Enter ("select same") with Allow Row Selection off, and over rows the filter hides.
   - More than 10 rows pinned (10 are drawn, the rest kept) and what Unpin All Rows leaves.
   - A pinned row deleted and brought back by Control+Z; a row pinned by an empty cell.
   - Design a Form... closed with its tab cross instead of Close and Apply; Default HTML Form with no
     column ticked; a malformed column Format (`0.0.0`).
3. **API or Dart tests to write**: a layout applied to a table missing its sort, pinned and frozen
   columns (ApiTests); a multi-line external text pasted into the last row (Dart regTest dispatching a
   `ClipboardEvent`).

### Won't do

- **Header height dragged to 0 and a row resizer dragged past 500 px**: no source says what the drag
  may reach.
