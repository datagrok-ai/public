@journey @viewers @realizes:viewers.tree-map
Feature: Tree map colour, size, selection and filtering
  What the rectangles are coloured by, what they are sized by, what a click on one selects, and
  what a filter on the table leaves in them.
  The spec this replaces proved the colour and the size scenarios with `waitForCanvasChange
  (minDelta: 300)` — a repaint, which any of the four property writes would have produced. The
  map now reports the two numbers the layout is actually built from: `area of <path>` is the
  score the squarified layout sizes the rectangle by (the row count by default, the aggregation
  under a Size By) and `color score of <path>` is the score the leaf is coloured by. So "Color
  Aggr Type = max" is 45.67 becoming 89, and "Size by WEIGHT" is 896 becoming 71221.11, while
  `rows of Caucasian` stays 896 throughout — which is the point of the two being separate.
  `color score of` is reported for leaves only: `_refreshColors` scores no group.
  The later scenarios translate the md's manual part on the same map: Shift-click adds a rectangle's
  rows to the selection and Control-click takes them off; Row Source All, Selected, Filtered
  Selected and Filtered build the map from those rows (`rows shown`); the viewer's Filter narrows it
  and leaves the table alone; the outer margins inset the layout box (`view`); and a saved layout
  brings the two levels and the colour column back after the viewer was closed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tree map viewer with:
      | splitByColumnNames | RACE |
    Then the "split columns" reading of tree map viewer should be "RACE"
    And the "leaves" reading of tree map viewer should be 4
    And the "rows shown" reading of tree map viewer should be 1000
    And the "rows of Caucasian" reading of tree map viewer should be 896
    And tree map viewer should be painted

  Scenario: Colouring by AGE scores every leaf and the aggregation moves the score
    Then the "color column" reading of tree map viewer should be ""
    And the "color of Caucasian" reading of tree map viewer should be "#2ca02c"
    When user picks "AGE" in the "color" column selector of tree map viewer
    Then the "color column" reading of tree map viewer should be "AGE"
    And the "color aggregation" reading of tree map viewer should be "avg"
    And the "color score of Caucasian" reading of tree map viewer should be between 45.6 and 45.7
    And the "color score of Asian" reading of tree map viewer should be between 38.0 and 38.1
    And tree map viewer should have repainted
    When user sets "colorAggrType" property of tree map viewer to "max"
    Then the "color score of Caucasian" reading of tree map viewer should be 89
    And the "color score of Asian" reading of tree map viewer should be 64
    And the "color of Caucasian" reading of tree map viewer should be "#ff0000"
    And the "rows of Caucasian" reading of tree map viewer should be 896
    And tree map viewer should have repainted
    When user sets properties of tree map viewer:
      | colorColumnName |     |
      | colorAggrType   | avg |
    Then the "color column" reading of tree map viewer should be ""
    And no errors should have been logged

  Scenario: Sizing by WEIGHT replaces the row count with the aggregation
    Then the "size column" reading of tree map viewer should be ""
    And the "area of Caucasian" reading of tree map viewer should be 896
    When user sets "sizeColumnName" property of tree map viewer to "WEIGHT"
    Then the "size aggregation" reading of tree map viewer should be "sum"
    And the "area of Caucasian" reading of tree map viewer should be between 71221 and 71222
    And the "area of Asian" reading of tree map viewer should be between 1057 and 1058
    And the "rows of Caucasian" reading of tree map viewer should be 896
    And tree map viewer should have repainted
    When user sets "sizeAggrType" property of tree map viewer to "max"
    Then the "area of Caucasian" reading of tree map viewer should be 165
    And the "area of Asian" reading of tree map viewer should be between 91.7 and 91.9
    And the "rows of Caucasian" reading of tree map viewer should be 896
    And tree map viewer should have repainted
    When user sets properties of tree map viewer:
      | sizeColumnName |     |
      | sizeAggrType   | sum |
    Then the "area of Caucasian" reading of tree map viewer should be 896
    And no errors should have been logged

  Scenario: Clicking a rectangle selects exactly the rows it holds, and the band is drawn over it
    Given user clears the row selection
    Then the "selected rows of Caucasian" reading of tree map viewer should be 0
    And tree map viewer should not have a "selection of Caucasian" area
    When user clicks on the "leaf Caucasian" area of tree map viewer
    Then 896 rows should be selected
    And only rows where "RACE" is "Caucasian" should be selected
    And the "selected rows of Caucasian" reading of tree map viewer should be 896
    And the "selected rows of Asian" reading of tree map viewer should be 0
    And tree map viewer should have a "selection of Caucasian" area
    And tree map viewer should not have a "selection of Asian" area
    When user clears the row selection
    Then no rows should be selected
    And tree map viewer should not have a "selection of Caucasian" area
    And no errors should have been logged

  Scenario: A filter on the table reshapes the map without re-splitting it
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of tree map viewer should be 447
    And the "rows of Caucasian" reading of tree map viewer should be 416
    And the "rows of Other" reading of tree map viewer should be 14
    And the "rows of Black" reading of tree map viewer should be 9
    And the "rows of Asian" reading of tree map viewer should be 8
    And the "leaves" reading of tree map viewer should be 4
    And the "split columns" reading of tree map viewer should be "RACE"
    And tree map viewer should have repainted
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And the "rows shown" reading of tree map viewer should be 1000
    And the "rows of Caucasian" reading of tree map viewer should be 896
    And no errors should have been logged

  Scenario: Shift-click adds a rectangle's rows to the selection and Control-click takes them off
    Given user clears the row selection
    When user clicks on the "leaf Asian" area of tree map viewer
    Then 15 rows should be selected
    When user clicks on the "leaf Black" area of tree map viewer holding Shift
    Then 42 rows should be selected
    And only rows where "RACE" is one of "Asian, Black" should be selected
    And the "selected rows of Black" reading of tree map viewer should be 27
    When user clicks on the "leaf Asian" area of tree map viewer holding Control
    Then 27 rows should be selected
    And only rows where "RACE" is "Black" should be selected
    And the "selected rows of Asian" reading of tree map viewer should be 0
    When user clears the row selection
    Then no errors should have been logged

  Scenario: Row Source decides which rows the map is built from
    When user adds a categorical filter on "SEX" keeping "M"
    Then the "rows shown" reading of tree map viewer should be 447
    When user sets "rowSource" property of tree map viewer to "All"
    Then the "rows shown" reading of tree map viewer should be 1000
    And the "rows of Caucasian" reading of tree map viewer should be 896
    When user selects rows where "RACE" is "Asian"
    And user sets "rowSource" property of tree map viewer to "Selected"
    Then the "rows shown" reading of tree map viewer should be 15
    And the "leaves" reading of tree map viewer should be 1
    When user sets "rowSource" property of tree map viewer to "FilteredSelected"
    Then the "rows shown" reading of tree map viewer should be 8
    When user sets "rowSource" property of tree map viewer to "Filtered"
    Then the "rows shown" reading of tree map viewer should be 447
    When user clears the row selection
    And user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then the "rows shown" reading of tree map viewer should be 1000
    And no errors should have been logged

  Scenario: The viewer's Filter formula narrows the map and leaves the table alone
    When user sets "filter" property of tree map viewer to "${AGE} > 40"
    Then the "rows shown" reading of tree map viewer should be 635
    And the "rows of Caucasian" reading of tree map viewer should be 562
    And the "rows of Asian" reading of tree map viewer should be 6
    And all rows should pass the filter
    When user sets "filter" property of tree map viewer to ""
    Then the "rows shown" reading of tree map viewer should be 1000
    And no errors should have been logged

  Scenario: The outer margins inset the box the rectangles are laid out in
    When user sets properties of tree map viewer:
      | outerMarginLeft   | 30 |
      | outerMarginTop    | 30 |
      | outerMarginRight  | 30 |
      | outerMarginBottom | 30 |
    Then the "view" area of tree map viewer should be narrower than before
    And the "view" area of tree map viewer should be shorter than before
    And the "rows of Caucasian" reading of tree map viewer should be 896
    When user sets properties of tree map viewer:
      | outerMarginLeft   | 0 |
      | outerMarginTop    | 0 |
      | outerMarginRight  | 0 |
      | outerMarginBottom | 0 |
    Then the "view" area of tree map viewer should be wider than before
    And no errors should have been logged

  Scenario: A saved layout brings the two levels and the colour column back
    When user sets properties of tree map viewer:
      | splitByColumnNames | RACE, SEX |
      | colorColumnName    | AGE       |
    Then the "levels" reading of tree map viewer should be 2
    And the "color column" reading of tree map viewer should be "AGE"
    When user saves the layout of the current table view
    And user clicks on close icon of tree map viewer
    Then the open tableview should have 0 tree map viewers
    When user loads the saved layout
    Then the open tableview should have 1 tree map viewer
    And the "split columns" reading of tree map viewer should be "RACE, SEX"
    And the "color column" reading of tree map viewer should be "AGE"
    And the "rows of Caucasian | F" reading of tree map viewer should be 480
    And the "color score of Caucasian | F" reading of tree map viewer should be a finite number
    And no errors should have been logged
