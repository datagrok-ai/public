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

1. **Decided and claimed** in `grid-edge-cases.feature` (localhost 2026-10-29, core master 575f6ae603).
   Known failures, each with its ticket and what the grid does today:
   - Frozen Columns typed in the context panel as 20 or -1 is not taken (the value is written), and is
     brought down to the column count when columns are removed (it stays at 6 with three data columns
     left): core#1137.
   - Allow Edit off: Shift+Delete keeps the rows (the selected row is deleted): core#1193.
   - Shift+Enter ("select same") selects nothing while Allow Row Selection is off (553 rows are
     selected), and only rows that pass the filter otherwise (471 selected rows are filtered out):
     core#1194 (`grid_row_selection.dart` selects over the whole table and does not check the setting).
   - A column group named with spaces only is not made (the group tag is set to the spaces): core#1195.
   - Twelve rows pinned show ten and the grid says that pinned rows are hidden (claimed as a warning
     balloon naming 10; the message's form is not decided), and Unpin All Rows also unpins the rows
     pinned beyond the ten shown (after Unpin All a layout pins 2 rows back): core#1191.
   - A pinned row deleted and brought back by Control+Z is pinned again (it comes back unpinned):
     core#1196.
   - Default HTML Form with no column ticked cannot be confirmed (OK is enabled): core#1197.
   - A malformed format (`0.0.0`) marks the Custom field invalid and is not taken (the field is not
     marked and the column's format tag becomes `0.0.0`): core#1198.
   - Row Height typed in the context panel as 0 or 501 is not taken (the property declares 1..500;
     the value is written): core#1200.

   Green: a value refused in the last row (`abc` in AGE under Add New Row On Last Row Edit) keeps the
   cell as it was and keeps the added row; a paste into a read-only grid changes nothing and shows no
   notice (a paste needs the cell in edit mode, and the double-click that opens it already says the
   grid is read-only); empty values go last in a sort both ways; the sort, a pinned row (through a
   layout) and a column group follow a renamed column.
2. **Steps requested**:
   - Design a Form... closed by the cross of its view tab works as Cancel (no form column). The tab's
     close icon is not found by `user clicks on close icon of {element}` (the tab reports one icon,
     "i"), and the view tab's context menu offers only Close others. Wanted: the close icon of a view
     tab, then `When user clicks on close icon of "Summary" view tab` and a claim on `column order`.
3. **Product suggestion**: Editable by could take its users from a user picker instead of free text
   (spaces only, another letter case). No scenario.
4. **API or Dart tests to write**: a layout applied to a table missing its sort, pinned and frozen
   columns (ApiTests); a multi-line external text pasted into the last row (Dart regTest dispatching a
   `ClipboardEvent`).

### Won't do

- **Header height dragged to 0 and a row resizer dragged past 500 px**: no source says what the drag
  may reach.
