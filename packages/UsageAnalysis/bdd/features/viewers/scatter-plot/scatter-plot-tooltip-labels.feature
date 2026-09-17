@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot tooltip and marker labels
  Hovering a marker makes its row the hovered one and shows the table's tooltip over it; a custom
  column list replaces what the tooltip shows and "do not show" takes it away without taking the
  hover away. The tooltip's own Show Column Names (the `showLabels` property — the Labels category
  carries a caption of the same name) is set to Always so the names are read, since Auto hides
  them for one number and one short category. Marker labels follow Label Columns and Show Labels
  For: every row, then only the selected ones, and none once the selection is cleared on
  marker-free space. One journey on demog-1000, X = WEIGHT, Y = HEIGHT; row 11 is AGE 46, SEX F,
  WEIGHT 118.9, HEIGHT 168.974; every scenario puts back what it changed.

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
