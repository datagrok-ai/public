@viewers @edge-cases @realizes:viewers.scatter-plot
Feature: Scatter plot edge cases: empty sizes, no rows, one row, columns removed or renamed, a layout with a missing column
  Guards for the data shapes and the changes under an open scatter plot. Rows with an empty value in
  the size column are not drawn until Show Markers With Empty Size is on (help scatter-plot.md, "X,
  Y, Colors, Sizes, Markers"): demog-1000 has 128 empty HEIGHTs, so switching the option on adds
  their hollow circles to a plot sized by HEIGHT and switching it off takes them away again (the ink
  of the plot; the `rows shown` reading counts these rows either way, see MISSING.md). A plot over
  rows the filter removed entirely, or over a one-row table, draws what passes and logs no error. A
  column the plot uses that is removed is dropped from its settings (the viewer re-runs its
  auto-setup, which forgets names the table lacks); a column that is renamed is followed by every
  column setting that named it. A layout that names a column the table no longer has is applied
  without an unhandled error. A plot that filters the table by its zoom gives the rows back when it
  is closed (a filter that outlives its viewer is a defect, core#233). Data: demog-1000, X = AGE, Y
  = WEIGHT (no empty values in either, so every row is drawable).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | WEIGHT |
    Then scatter plot viewer should show 1000 rows

  Scenario: Rows with an empty size are hidden until Show Markers With Empty Size is on
    When user sets "Size" property of scatter plot viewer to "HEIGHT"
    And user sets "Show Markers With Empty Size" property of scatter plot viewer to "true"
    Then scatter plot viewer should have more ink than before
    When user sets "Show Markers With Empty Size" property of scatter plot viewer to "false"
    Then scatter plot viewer should have less ink than before
    And no errors should have been logged

  Scenario: A filter that removes every row leaves the plot empty without an error
    When user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And scatter plot viewer should show 0 rows
    And scatter plot viewer should report no error
    When user resets the filter
    Then scatter plot viewer should show 1000 rows
    And no errors should have been logged

  Scenario: A one-row table draws its one marker
    Given user opens a table "one-row" with:
      | name | x | y |
      | a    | 1 | 2 |
    And user adds a scatter plot viewer with:
      | X | x |
      | Y | y |
    Then scatter plot viewer should show 1 rows
    And scatter plot viewer should have a "marker of row 1" area
    And no errors should have been logged

  Scenario: Removing the X column drops it from the plot's settings without an error
    When user removes "AGE" column
    Then the table should not have a column "AGE"
    And "X" property of scatter plot viewer should not be "AGE"
    And no errors should have been logged

  Scenario: Removing the Color column drops it from the plot's settings without an error
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 4 items
    When user removes "RACE" column
    Then "Color" property of scatter plot viewer should not be "RACE"
    And scatter plot viewer should show 1000 rows
    And no errors should have been logged

  Scenario: Renaming the X and Color columns is followed by the plot's settings
    When user sets "Color" property of scatter plot viewer to "RACE"
    And user renames "AGE" column to "AGE_YEARS"
    And user renames "RACE" column to "ORIGIN"
    Then "X" property of scatter plot viewer should be "AGE_YEARS"
    And "Color" property of scatter plot viewer should be "ORIGIN"
    And scatter plot viewer should show 1000 rows
    And the legend of scatter plot viewer should list 4 items
    And no errors should have been logged

  Scenario: A layout naming a column the table no longer has is applied without an error
    When user sets "Color" property of scatter plot viewer to "RACE"
    And user saves the layout of the current table view
    And user removes "RACE" column
    And user loads the saved layout
    Then scatter plot viewer should be visible
    And scatter plot viewer should show 1000 rows
    And no errors should have been logged

  Scenario: Closing the plot gives back the rows its zoom filtered out
    When user sets "Zoom and Filter" property of scatter plot viewer to "filter by zoom"
    And user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then fewer than 1000 rows should pass the filter
    When user clicks on close icon of scatter plot viewer
    Then the open tableview should have 0 scatter plot viewers
    And all rows should pass the filter
    And no errors should have been logged
