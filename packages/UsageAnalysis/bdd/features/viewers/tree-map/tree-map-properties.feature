@journey @viewers @realizes:viewers.tree-map
Feature: Tree map selection band, Row Source Selected, a filtering bar chart, Auto Layout and the categorical gate
  The tree map properties and ticket regressions the two other tree map features leave out.
  Show Selected Rows decides whether the band over a rectangle's selected share is drawn, not
  whether the rows are selected: `selection of <path>` is the band's own area, reported under the
  renderer's very condition, so switching it off takes the area away while `selected rows of
  <path>` and the table's selection stay. Row Source Selected builds the map from the selection
  itself, so it draws no band at all. A bar chart whose On Click is Filter narrows the tree map
  next to it to the clicked bar's rows (public#935: SEX = M is 447 rows of demog-1000, 416 of
  them Caucasian). Auto Layout keeps the on-viewer selectors only while the map is over 200 by
  200 (public#2131), so a 180 by 180 map loses them and keeps them once Auto Layout is off. A table
  without a categorical column gets the descriptor's own gate, "Requires categorical columns".
  The split selectors themselves are not driven here: see MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tree map viewer with:
      | splitByColumnNames | RACE |
    Then the "split columns" reading of tree map viewer should be "RACE"
    And the "leaves" reading of tree map viewer should be 4
    And the "rows shown" reading of tree map viewer should be 1000
    And the "rows of Asian" reading of tree map viewer should be 15
    And tree map viewer should be painted

  Scenario: Show Selected Rows off takes the selection band away and keeps the rows selected
    Given user clears the row selection
    When user selects rows where "RACE" is "Asian"
    Then 15 rows should be selected
    And the "selected rows of Asian" reading of tree map viewer should be 15
    And tree map viewer should have a "selection of Asian" area
    And tree map viewer should not have a "selection of Caucasian" area
    When user sets "showSelectedRows" property of tree map viewer to "false"
    Then tree map viewer should not have a "selection of Asian" area
    And the "leaf Asian" area of tree map viewer should have repainted
    And 15 rows should be selected
    And only rows where "RACE" is "Asian" should be selected
    And the "selected rows of Asian" reading of tree map viewer should be 15
    And tree map viewer should have a "leaf Asian" area
    When user sets "showSelectedRows" property of tree map viewer to "true"
    Then tree map viewer should have a "selection of Asian" area
    And the "leaf Asian" area of tree map viewer should have repainted
    When user clears the row selection
    Then tree map viewer should not have a "selection of Asian" area
    And no errors should have been logged

  Scenario: Row Source Selected builds the map from the selection and draws no band
    Given user clears the row selection
    When user selects rows where "RACE" is one of "Asian, Black"
    Then 42 rows should be selected
    And tree map viewer should have a "selection of Asian" area
    And tree map viewer should have a "selection of Black" area
    When user sets "rowSource" property of tree map viewer to "Selected"
    Then the "row source" reading of tree map viewer should be "Selected"
    And the "rows shown" reading of tree map viewer should be 42
    And the "leaves" reading of tree map viewer should be 2
    And the "selected rows of Asian" reading of tree map viewer should be 15
    And the "selected rows of Black" reading of tree map viewer should be 27
    And tree map viewer should have a "leaf Asian" area
    And tree map viewer should not have a "selection of Asian" area
    And tree map viewer should not have a "selection of Black" area
    When user sets "rowSource" property of tree map viewer to "Filtered"
    Then the "rows shown" reading of tree map viewer should be 1000
    And tree map viewer should have a "selection of Asian" area
    When user clears the row selection
    Then no errors should have been logged

  Scenario: A bar chart with On Click Filter narrows the tree map to the clicked bar (public#935)
    When user adds a bar chart viewer with:
      | Split           | SEX      |
      | Value Aggr Type | count    |
      | On Click        | Filter   |
    Then the "bars" reading of bar chart viewer should be 2
    And all rows should pass the filter
    When user clicks on the "bar M" area of bar chart viewer
    Then the filter should pass exactly the rows where "SEX" is "M"
    And 447 rows should pass the filter
    And the "rows shown" reading of tree map viewer should be 447
    And the "rows of Caucasian" reading of tree map viewer should be 416
    And the "rows of Asian" reading of tree map viewer should be 8
    And the "leaves" reading of tree map viewer should be 4
    When user clicks on empty plot space of bar chart viewer
    Then all rows should pass the filter
    And the "rows shown" reading of tree map viewer should be 1000
    And the "rows of Caucasian" reading of tree map viewer should be 896
    When user clicks on close icon of bar chart viewer
    Then the open tableview should have 0 bar chart viewers
    And no errors should have been logged

  Scenario: Auto Layout hides the selectors on a small map and keeps them once it is off (public#2131)
    Then the "selection panel shown" reading of tree map viewer should be "true"
    And tree map viewer should have a "split selector 1" area
    And tree map viewer should have a "color selector" area
    When user resizes tree map viewer to 180 by 180
    Then the "selection panel shown" reading of tree map viewer should be "false"
    And tree map viewer should not have a "split selector 1" area
    And tree map viewer should not have a "color selector" area
    And the "leaves" reading of tree map viewer should be 4
    When user sets "autoLayout" property of tree map viewer to "false"
    Then the "selection panel shown" reading of tree map viewer should be "true"
    And tree map viewer should have a "split selector 1" area
    And tree map viewer should have a "color selector" area
    When user sets "autoLayout" property of tree map viewer to "true"
    Then the "selection panel shown" reading of tree map viewer should be "false"
    And tree map viewer should not have a "split selector 1" area
    When user restores the size of tree map viewer
    Then the "selection panel shown" reading of tree map viewer should be "true"
    And tree map viewer should have a "split selector 1" area
    And no errors should have been logged

  Scenario: A table with no categorical column gets the categorical gate
    Given user opens a table "tree-map-numbers" with:
      | x | y   |
      | 1 | 2.5 |
      | 2 | 3.5 |
      | 3 | 4.5 |
    And user adds a tree map viewer
    Then the "error" reading of tree map viewer should be "Requires categorical columns"
    And no errors should have been logged
