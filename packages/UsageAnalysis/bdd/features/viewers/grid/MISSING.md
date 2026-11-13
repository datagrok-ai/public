# Grid — what the features cannot say yet

## Resolved

- **core#1195, core#1194 (Allow Row Selection)**: fixed in this round; a blank group name ungroups
  (`grid-column-groups.feature`), Shift+Enter selects nothing while Allow Row Selection is off (`grid-rows.feature`).
- **core#1137 (Pin / Unpin Column), core#1191**: fixed in this round; Dart tests "d4 | Viewers | Grid | Pin and Unpin
  Column work with Frozen Columns above the column count" and "… | Rows pinned past the cap leave no values for a layout".

## Still open

- **Header histogram and stats rows** (`grid-ui.md` "Column Stats"): need `header histogram <col>` areas and
  `text of <stat> of <col>` readings from the grid (`grid_special_rows.dart`).
- **Form designer field positions** (`grid-ui.md` "Summary Columns — Form Designer"): the field boxes need names and
  positions in the designer's status.
- **The column popup's categorical filter** (`grid-ui.md` "Column Hamburger Menu Inline Filter"): needs
  `category <label>` areas like the panel's card.
- **public#989, public#996**: a cut-off cell's tooltip is empty again (core#441); claim in `grid-viewer.feature`.
- **core#369**: nested linked-table sub-columns need `cell <r> of <table>/<col>` hit areas.
- **core#1198**: Format > Custom... takes a malformed pattern (0.0.0) and the cells fall back silently; claim when fixed:
  `Custom input in Format AGE dialog should be invalid` in `grid-appearance.feature`'s custom-format scenario.
- **core#1196**: a pinned row deleted and brought back by Control+Z is no longer pinned; claim when fixed in
  `grid-pinning.feature`: Pin Row, Shift+Delete, Control+Z, then `pinned rows` 1.
- **core#1197**: Default HTML Form's column dialog takes OK with no column ticked; claim when fixed:
  `OK button in "Select columns..." dialog should be disabled`.
- **core#1137 (typed values)**: Frozen Columns declares no minimum, so -1 typed in the panel is kept; claim in the Dart
  test "d4 | Viewers | Fixed issues | A property keeps a value inside its declared bounds" once `min: 0` is declared.
- **core#1193**: Shift+Delete removes the rows of a grid with Allow Edit off — awaiting a product decision.
- **core#1194 (the filter)**: Shift+Enter also selects the matching rows the filter hides — awaiting a product decision.

## Won't do

- **Header height dragged to 0, a row resizer dragged past 500 px**: no source says what the drag may reach.
- **Edge cases of the 2026-10-26 round dropped in review**: Row Height clamped to 1..500 (Dart "A property keeps a
  value inside its declared bounds"); empty values last in a sort (ddt `sorted_order_keys_test.dart`); the sort, a
  pinned row and a column group following a renamed column (Dart "d4 | Viewers | Grid | Sort and pinned rows follow a
  column rename"; a group is the column's own tag); an html summary column after its source column is removed
  (PowerGrid `summary-columns.feature`, GROK-19942); a value refused in the last row (no source says what it should do).
- **core#417, the Sort dialog after a rename**: Dart "d4 | Inputs | Column input | Renamed column"; no user can rename a
  column while the dialog is open.
- **public#3787, the column-menu and settings icons overlapping**: cosmetic.
- **public#2691, an image renderer in a pinned column**: a renderer's drawing, for Curves' package tests.
- **public#1812 part 2, no source row selected**: not reachable, the selection-to-filter link stops filtering.
- **Design a Form... closed by its tab's cross, Editable by as a user picker, a layout over a table missing its sort,
  pinned or frozen columns, a multi-line external paste**: not planned.
