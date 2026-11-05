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
