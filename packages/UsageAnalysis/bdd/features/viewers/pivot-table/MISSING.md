# Pivot table — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `pivot_viewer/CLAUDE.md` "Automation
surface" and `pivot-table.md`. The Table switch is in `pivot-table-table-switch.feature`.

## Still open

- The + picker of the Aggregate row offers the type picked last, and Cancel adds nothing (scenario 6): needs a step that
  picks the type in the picker's aggregation selector (`pivot_grid.dart` `showColumnComboBox`).
- The tag-editor header duplicated on a Table switch (GROK-14995): needs a `tag rows` reading
  (`pivot_viewer_status.dart`).
- core#467 — a key column header dragged inside the inner grid logs a NullError on every repaint: the known failure of
  `pivot-table-inner-grid-drag.feature`.
- The Group by picker of a pivot docked at the right edge opens inside the window (public#3690): needs a `column picker
  popup` kind for `.d4-column-grid:has(.d4-column-selector-backdrop)`.

## Won't do

- ADD with every row filtered out publishes an empty table: allowed, the user's choice; no rows is the Dart Data
  Manipulation loop.
- A layout grouping by a column the table lacks: xamgle `d4 | Viewers | Invalid settings | Pivot table reports a group-by
  column the table lacks and draws once it is fixed`.
- The Group by column removed (another is picked), thousands of categories, column names with `, [ ]` or quotes in saved
  parameters: the rules every viewer shares (`../MISSING.md`) / low value.
