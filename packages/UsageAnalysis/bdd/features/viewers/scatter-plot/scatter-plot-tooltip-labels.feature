@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot tooltip and marker labels
  Hovering a marker makes its row the hovered one and shows the table's tooltip over it; a custom
  column list replaces what the tooltip shows and "do not show" takes it away without taking the
  hover away. The tooltip's own Show Column Names (the `showLabels` property — the Labels category
  carries a caption of the same name) is set to Always so the names are read, since Auto hides
  them for one number and one short category. Marker labels follow Label Columns and Show Labels
  For: every row, then only the selected ones, and none once the selection is cleared on
  marker-free space. Labels on a datetime X axis (STARTED, with RACE on Y) survive a RACE card that
  keeps one category (Caucasian, 896 rows: the md keeps the first one listed, but a box over the inner
  80% of the plot selects none of the 15 Asian rows, which sit at the axis's edge) and a box selection over what is left (GROK-17251; the md's
  spgi-100 date column and Stereo Category become demog-1000's STARTED and RACE). One journey on
  demog-1000, X = WEIGHT, Y = HEIGHT. Closed issues: labels for selected rows the filter took out
  are not drawn (public#3189), and every label of the selected rows comes back after a wheel zoom in
  and out (public#3176). Row 11 is AGE 46, SEX F, WEIGHT 118.9, HEIGHT 168.974; every
  scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows
    And properties of scatter plot viewer should be:
      | Show Tooltip  | inherit from table |
      | Label Columns |                    |

  Scenario: Hovering a marker makes its row hovered and shows the table's tooltip
    Then the "hovered row" reading of scatter plot viewer should be 0
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the "hovered row" reading of scatter plot viewer should be 11
    And tooltip should be visible
    And the tooltip should show some columns
    And the tooltip should show "HEIGHT" as "168.974"
    When user moves the pointer away from scatter plot viewer
    Then tooltip should be hidden
    And no errors should have been logged

  Scenario: A custom column list replaces what the tooltip shows
    When user sets properties of scatter plot viewer:
      | Show Tooltip | show custom tooltip |
      | Row Tooltip  | AGE\nSEX            |
      | Data Values  | Do not add          |
      | showLabels   | Always              |
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the "hovered row" reading of scatter plot viewer should be 11
    And the tooltip should show columns "AGE, SEX"
    And the tooltip should show "AGE" as "46"
    And the tooltip should show "SEX" as "F"
    When user moves the pointer away from scatter plot viewer
    And user sets "Show Tooltip" property of scatter plot viewer to "do not show"
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the "hovered row" reading of scatter plot viewer should be 11
    And tooltip should be hidden
    When user moves the pointer away from scatter plot viewer
    And user sets properties of scatter plot viewer:
      | Show Tooltip | inherit from table |
      | Row Tooltip  |                    |
      | Data Values  | Merge              |
      | showLabels   | Auto               |
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then tooltip should be visible
    And the tooltip should show some columns
    When user moves the pointer away from scatter plot viewer
    Then no errors should have been logged

  Scenario: Labels are drawn for every row, then only for the selected ones
    Then no rows should be selected
    And the "labels shown" reading of scatter plot viewer should be 0
    When user sets "Label Columns" property of scatter plot viewer to "AGE"
    Then the "labels shown" reading of scatter plot viewer should be at least 1
    When user sets "Show Labels For" property of scatter plot viewer to "Selected"
    Then the "labels shown" reading of scatter plot viewer should be 0
    When user drags a selection box over the "view" area of scatter plot viewer
    Then some rows should be selected
    And the "labels shown" reading of scatter plot viewer should be higher than before
    When user clicks on the "empty space" area of scatter plot viewer
    Then no rows should be selected
    And the "labels shown" reading of scatter plot viewer should be 0
    When user sets properties of scatter plot viewer:
      | Label Columns   |     |
      | Show Labels For | All |
    Then the "labels shown" reading of scatter plot viewer should be 0
    And no errors should have been logged

  Scenario: Labels on a datetime axis survive a filter that keeps one category and a selection
    When user sets properties of scatter plot viewer:
      | X               | STARTED  |
      | Y               | RACE     |
      | Label Columns   | USUBJID  |
      | Show Labels For | Selected |
    Then scatter plot viewer should show 1000 rows
    When user adds a categorical filter on "RACE" keeping "Caucasian"
    Then 896 rows should pass the filter
    And scatter plot viewer should show 896 rows
    When user drags a selection box over the "view" area of scatter plot viewer
    Then some rows should be selected
    And every selected row should pass the filter
    And the "labels shown" reading of scatter plot viewer should be at least 1
    And scatter plot viewer should be painted
    And no errors should have been logged
    When user adds a categorical filter on "RACE" keeping "Asian, Black, Caucasian, Other"
    And user clicks on the "empty space" area of scatter plot viewer
    Then all rows should pass the filter
    And no rows should be selected
    When user sets properties of scatter plot viewer:
      | X               | WEIGHT |
      | Y               | HEIGHT |
      | Label Columns   |        |
      | Show Labels For | All    |
    Then scatter plot viewer should show 872 rows
    And no errors should have been logged

  @realizes:public#3189
  Scenario: Labels for selected rows are not drawn for the rows the filter took out (public#3189)
    When user sets properties of scatter plot viewer:
      | Label Columns   | AGE      |
      | Show Labels For | Selected |
    And user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    When user selects rows where "SEX" is "M"
    Then the "labels shown" reading of scatter plot viewer should be 373
    When user selects rows where "SEX" is "F"
    Then 553 rows should be selected
    And 0 selected rows should pass the filter
    And the "labels shown" reading of scatter plot viewer should be 0
    When user clears the row selection
    And user adds a categorical filter on "SEX" keeping "F, M"
    And user sets properties of scatter plot viewer:
      | Label Columns   |     |
      | Show Labels For | All |
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  @realizes:public#3176
  Scenario: Every label of the selected rows comes back after a zoom in and out (public#3176)
    When user sets properties of scatter plot viewer:
      | Zoom and Filter | no action |
      | Label Columns   | AGE       |
      | Show Labels For | Selected  |
    And user selects the first 5 rows
    Then the "labels shown" reading of scatter plot viewer should be 5
    When user remembers the value range of scatter plot viewer
    And user scrolls the mouse wheel up 4 times over the "view" area of scatter plot viewer
    Then scatter plot viewer should show a narrower value range than before
    When user scrolls the mouse wheel down 4 times over the "view" area of scatter plot viewer
    Then the "labels shown" reading of scatter plot viewer should be 5
    When user clears the row selection
    And user picks "Reset View" from the context menu of scatter plot viewer
    And user sets properties of scatter plot viewer:
      | Label Columns   |                |
      | Show Labels For | All            |
      | Zoom and Filter | filter by zoom |
    Then the "labels shown" reading of scatter plot viewer should be 0
    And all rows should pass the filter
    And no errors should have been logged
