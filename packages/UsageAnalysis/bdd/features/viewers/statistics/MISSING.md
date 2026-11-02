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

## Edge cases (now in `statistics-menu.feature`)

Written for the guard suite of edge cases of 26 Oct 2026 (its check is now an `@edge-cases`
scenario of its own, after the master scenarios of `statistics-menu.feature`): data
shapes, columns removed or renamed under the open viewer, layouts naming a missing column, property
bounds, closing the viewer. Nothing here asks for a library, binding or code change made by that
feature; it lists what the feature could not claim and why.

### Edge cases: Still open

#### 1. Steps requested

None.

#### 2. Decided

1. **Numeric statistics over no values are empty, not 0.** `avg of AGE` empty with every row
   filtered out, an all-empty numeric column with avg and stdev empty (core#1125, fixed), and a
   constant column (min, max, avg 5.00, stdev 0.00) are requested from ApiTests below.
2. **A Histograms > <column> column** is dropped when that column is removed: claimed green in
   `statistics-menu.feature`. That it follows a rename (`[SEX]` is built once by name,
   `stats_viewer_core.dart:111-122`) is not written: `histogram columns` is a count, which a column
   kept under the old name satisfies too; it needs the histogram columns' names as a reading.

#### 3. Questions for owners

None. A column added while Columns is an explicit list gets no row (2026-10-30): requested from
ApiTests below.

#### 4. API or Dart tests to write

1. ApiTests: a filter that passes no row shows 0 rows, 0 values for AGE and RACE, an empty avg of
   AGE and no error; resetting it brings back 1000 rows and avg 45.68.
2. ApiTests: Row Source Selected with nothing selected shows 0 rows and 0 values for AGE, no error.
3. ApiTests: a one-row table shows 1 row, 1 value and 0 nulls for its column.
4. ApiTests: a calculated column added while the viewer is open gets a row (12 columns, 1000 values).
5. ApiTests: removing RACE takes its row away (10 columns) and keeps avg of AGE at 45.68.
6. ApiTests: AGE renamed to AGE_YEARS keeps its row under the new name (avg 45.68), with no row AGE.
7. ApiTests: an int column with no values has 0 values, 1000 nulls and empty avg, stdev and min
   (core#1125, fixed).
8. ApiTests: a constant column of 5 has min, max and avg 5.00 and stdev 0.00.
9. ApiTests: with Columns set to AGE, HEIGHT, a calculated column added afterwards gets no row.

Moved out of BDD (no user gesture); requested from ApiTests.

### Edge cases: Won't do

Nothing.
