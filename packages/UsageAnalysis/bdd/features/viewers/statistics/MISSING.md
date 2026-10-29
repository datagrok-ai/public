# Statistics — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `stats_viewer/CLAUDE.md` "Automation
surface" and `statistics-ui.md` / `statistics.md`. Back Color is in `statistics-back-color.feature`.

## Needs a reading or an area

- **Reorder stat columns** (`statistics-ui.md`): dragging the `max` header left past `min`. The
  status re-exports the inner grid's cells and headers but not its `column order`, so a reorder has
  nothing to be claimed on. Needs `column order` re-exported from the inner grid in statistic names
  (`features/` status part of `stats_viewer`):

  ```gherkin
  Scenario: Dragging a statistic's header moves its column
    Then the "column order" reading of statistics viewer should be "values, nulls, unique, min, max, avg, med, stdev"
    When user drags the "header max" area of statistics viewer to the "header min" area
    Then the "column order" reading of statistics viewer should be "values, nulls, unique, max, min, avg, med, stdev"
    And the "max of AGE" reading of statistics viewer should be "89"
    And no errors should have been logged
  ```

- **Resize stat columns** (`statistics-ui.md`): the inner grid's `column resizer <col>` areas and
  `column width of <col>` readings are not re-exported. Needs `column resizer <stat>` /
  `column width of <stat>` (and the name column's), so the drag is by name:

  ```gherkin
  Scenario: Dragging a statistic's edge widens it and keeps its numbers whole
    When user drags the "column resizer avg" area of statistics viewer by 40 pixels to the right
    Then the "column width of avg" reading of statistics viewer should be higher than before
    And the "avg of HEIGHT" reading of statistics viewer should not be ""
    And no errors should have been logged
  ```

- **Full-screen shortcut** (`statistics.md`, "Full-screen shortcut"): the panel's expand icon and
  Alt+F. Nothing reports that a viewer is in full screen (`General > Full Screen` is only listed by
  `viewer-chrome.feature`). Needs a `full screen` reading on the base viewer status (or a library
  state `{widget} should be in full screen`, from the dock manager's maximized panel):

  ```gherkin
  Scenario: Full screen and back with Alt+F
    When user picks "General > Full Screen" from the viewer menu of statistics viewer
    Then statistics viewer should be in full screen
    When user presses Alt+F
    Then statistics viewer should not be in full screen
    And the "columns shown" reading of statistics viewer should be 11
  ```

## Edge cases (`statistics-edge-cases.feature`)

Written for `statistics-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data
shapes, columns removed or renamed under the open viewer, layouts naming a missing column, property
bounds, closing the viewer. Nothing here asks for a library, binding or code change made by that
feature; it lists what the feature could not claim and why.

### Edge cases: Still open

#### 1. Steps requested

None.

#### 2. Decided (the rules of `viewer-rules-edge-cases.feature`)

1. **Numeric statistics over no values are empty, not 0.** `statistics-edge-cases.feature` claims
   `avg of AGE` empty with every row filtered out, and an all-empty numeric column with avg and
   stdev empty: a known failure with no ticket yet (stdev reads 0.00). A constant column (min, max,
   avg 5.00, stdev 0.00) is claimed green.
2. **A Histograms > <column> column** is dropped when that column is removed: claimed green in
   `statistics-edge-cases.feature`. That it follows a rename (`[SEX]` is built once by name,
   `stats_viewer_core.dart:111-122`) is not written: `histogram columns` is a count, which a column
   kept under the old name satisfies too; it needs the histogram columns' names as a reading.

#### 3. Questions for owners

1. **A column added while Columns is an explicit list**: does it get a row? (With Columns empty it
   does; claimed.)

#### 3. API or Dart tests to write

None.

### Edge cases: Won't do

Nothing.
