@journey @viewers @realizes:viewers.heat-map
Feature: Heat map columns, sorting and column widths
  Which columns the map shows, in which row order, and how wide. The Columns property (core#376)
  says how many columns are visible ("11 visible") and its "..." opens the Order or Hide Columns
  dialog, whose unchecked column leaves the map: `columns shown`, `column order` and the column's
  header area all move, not only the text in the property row. A double-click on a header sorts
  the map's own rows (public#182): descending first, then ascending, then back to the table's
  order, read as `sort column` / `sort direction` and as `row order`, the first ten table rows the
  map draws top to bottom, so a sort that is reported but not drawn fails. The heat map gives every
  column the same width and does not let a border drag change it (GROK-18672: resizing was removed
  from heat map mode); the same drag on the same viewer with Is Heatmap off widens the column, so
  the drag is shown to land, and turning the heat map back on gives that column the uniform width
  again. Max Heatmap Columns 0 (GROK-17776) used to break the layout with a
  division by a zero width: it now leaves no column on the map and logs no error, and 100 brings
  all eleven back. Since core#928 the empty map says why: "No columns to show: Max Heatmap Columns is 0",
  and 100 takes the hint away (core#928). The hint is read as text shown in the viewer: the grid's
  status does not publish it as its `error` (see MISSING.md).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "columns shown" reading of heat map viewer should be 11
    And heat map viewer should be painted

  Scenario: The Columns property counts the visible columns and its dialog hides one from the map (core#376)
    When user clicks on settings icon of heat map viewer
    And user expands "Columns" category in context panel
    Then "Columns" property in context panel should have value "11 visible"
    When user clicks on first button in "Columns" property in context panel
    Then Order or Hide Columns dialog should be visible
    And the "AGE" column should be checked in the column list of Order or Hide Columns dialog
    When user takes a snapshot of heat map viewer
    And user toggles the "AGE" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then the "columns shown" reading of heat map viewer should be 10
    And the "column order" reading of heat map viewer should be "USUBJID, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And heat map viewer should not have a "header AGE" area
    And heat map viewer should not have a "column AGE" area
    And heat map viewer should have repainted
    And "Columns" property in context panel should have value "10 visible"
    When user clicks on first button in "Columns" property in context panel
    And user toggles the "AGE" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then the "columns shown" reading of heat map viewer should be 11
    And the "column order" reading of heat map viewer should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And heat map viewer should have a "column AGE" area
    And "Columns" property in context panel should have value "11 visible"
    And no errors should have been logged

  Scenario: A header double-click sorts descending, then ascending, then puts the table's order back (public#182)
    Then the "sort column" reading of heat map viewer should be ""
    And the "row order" reading of heat map viewer should be "1, 2, 3, 4, 5, 6, 7, 8, 9, 10"
    When user double-clicks on the "header AGE" area of heat map viewer
    Then the "sort column" reading of heat map viewer should be "AGE"
    And the "sort direction" reading of heat map viewer should be "descending"
    And the "row order" reading of heat map viewer should be "600, 720, 597, 469, 577, 639, 645, 848, 620, 440"
    And "AGE" column should have its maximum in row 600
    And heat map viewer should have repainted
    When user double-clicks on the "header AGE" area of heat map viewer
    Then the "sort direction" reading of heat map viewer should be "ascending"
    And the "row order" reading of heat map viewer should be "692, 816, 89, 460, 479, 988, 62, 672, 795, 893"
    And heat map viewer should have repainted
    When user double-clicks on the "header AGE" area of heat map viewer
    Then the "sort column" reading of heat map viewer should be ""
    And the "sort direction" reading of heat map viewer should be ""
    And the "row order" reading of heat map viewer should be "1, 2, 3, 4, 5, 6, 7, 8, 9, 10"
    And heat map viewer should have repainted
    And no errors should have been logged

  Scenario: A column border drag does not resize a heat map column, and does once Is Heatmap is off (GROK-18672)
    # The grid still reports `column resizer <col>` in heat map mode (grid_status.dart), a drag
    # target the heat map ignores (`ColumnResizing.check` refuses `isHeatmap`): see MISSING.md.
    Then the "column width of AGE" reading of heat map viewer should be 30
    When user drags the "column resizer AGE" area of heat map viewer by 40 pixels to the right
    Then the "column width of AGE" reading of heat map viewer should be 30
    And the "column width of SEX" reading of heat map viewer should be 30
    And the "column AGE" and "column SEX" areas of heat map viewer should be the same width
    When user sets "isHeatmap" property of heat map viewer to "false"
    Then the "is heatmap" reading of heat map viewer should be "false"
    When user remembers the "column width of AGE" reading of heat map viewer
    And user drags the "column resizer AGE" area of heat map viewer by 40 pixels to the right
    Then the "column width of AGE" reading of heat map viewer should be higher than remembered
    When user sets "isHeatmap" property of heat map viewer to "true"
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "row height" reading of heat map viewer should be between 0 and 8
    And the "column width of AGE" reading of heat map viewer should be 30
    And the "column AGE" and "column SEX" areas of heat map viewer should be the same width
    And no errors should have been logged

  @realizes:core#928
  Scenario: Max Heatmap Columns 0 leaves no column on the map and says why (GROK-17776, core#928)
    When user sets "maxHeatmapColumns" property of heat map viewer to "0"
    Then the "max heatmap columns" reading of heat map viewer should be 0
    And the "columns shown" reading of heat map viewer should be 0
    And the "column order" reading of heat map viewer should be ""
    And heat map viewer should not have a "column AGE" area
    And heat map viewer should show the text "No columns to show: Max Heatmap Columns is 0"
    And heat map viewer should have repainted
    When user sets "maxHeatmapColumns" property of heat map viewer to "100"
    Then the "columns shown" reading of heat map viewer should be 11
    And heat map viewer should have a "column AGE" area
    And heat map viewer should not show the text "No columns to show"
    And heat map viewer should be painted
    And no errors should have been logged
