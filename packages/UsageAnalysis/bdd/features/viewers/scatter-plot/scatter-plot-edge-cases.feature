@viewers @edge-cases @realizes:viewers.scatter-plot
Feature: Scatter plot edge cases: empty sizes, one row, the X column removed, a filter by zoom left on close
  Guards for the data shapes and the changes under an open scatter plot. Rows with an empty value in
  the size column are not drawn until Show Markers With Empty Size is on (help scatter-plot.md, "X,
  Y, Colors, Sizes, Markers"): demog-1000 has 128 empty HEIGHTs, so switching the option on adds
  their hollow circles to a plot sized by HEIGHT and switching it off takes them away again (the ink
  of the plot; the `rows shown` reading counts these rows either way, see MISSING.md). A plot over a
  one-row table draws its marker and logs no error. A removed X column is dropped from the settings
  (the viewer re-runs its auto-setup, which forgets names the table lacks). A plot that filters the
  table by its zoom gives the rows back when it is closed (a filter that outlives its viewer is a
  defect, core#233). No rows, a removed Color or Size column, renamed columns and a layout naming a
  missing column are rows of `viewer-rules-edge-cases.feature`. Data: demog-1000, X = AGE, Y =
  WEIGHT (no empty values in either, so every row is drawable).

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

  Scenario: Closing the plot gives back the rows its zoom filtered out
    When user sets "Zoom and Filter" property of scatter plot viewer to "filter by zoom"
    And user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then fewer than 1000 rows should pass the filter
    When user clicks on close icon of scatter plot viewer
    Then the open tableview should have 0 scatter plot viewers
    And all rows should pass the filter
    And no errors should have been logged

  Scenario Outline: Moving Average Window typed in the context panel as <typed> is not taken
    When user sets "Show Moving Average Line" property of scatter plot viewer to "true"
    And user clicks on settings icon of scatter plot viewer
    Then context panel should be visible
    Given "Statistics" category in context panel is expanded
    When user enters "<typed>" into "Moving Average Window" property in context panel
    Then "Moving Average Window" property of scatter plot viewer should be "10"
    And no errors should have been logged

    @known-failure:63 @realizes:core#1200
    Examples: The typed value is taken although the property declares 2..100
      | typed |
      | 1     |
      | 101   |
