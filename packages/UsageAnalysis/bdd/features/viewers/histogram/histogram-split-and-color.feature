@journey @viewers @realizes:viewers.histogram
Feature: Histogram split, stacking and color coding
  Color coding by a numerical column and a split are alternatives: setting Split disables the three
  Color properties in the context panel and replaces the bars with one spline per category, which
  answer to hover and to a click as a whole category; Split Stack brings the bars back as segments,
  one per category, and Normalize Values decides what the vertical axis is scaled to. A range
  narrowed under a split still leaves the table filter valid. One journey on demog-1000 with a
  histogram of AGE; every scenario puts back what it changed.
  Closed issues: a split with a single bin can be hovered without an error (core#998); the split
  lines take the split column's category colours (public#130).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | Value | AGE |
    And user resizes histogram viewer to 500 by 400
    And user clicks on settings icon of histogram viewer
    Then histogram viewer should show 1000 rows
    And histogram viewer should have a "bin 8" area

  Scenario: Color coding by a numerical column
    Then "Color" property in context panel should be enabled
    And "Color Aggr Type" property in context panel should be disabled
    When user sets "Color" property of histogram viewer to "HEIGHT"
    Then histogram viewer should have repainted by at least 2000 pixels
    And the "bin 1" area of histogram viewer should contain the color "#FF0000"
    And the "bin 17" area of histogram viewer should contain the color "#0000FF"
    And the "bin 1" and "bin 17" areas of histogram viewer should be painted in different colors
    And "Color Aggr Type" property in context panel should be enabled
    And "Invert Color Scheme" property in context panel should be enabled
    When user sets "Color Aggr Type" property of histogram viewer to "stdev"
    Then histogram viewer should have repainted by at least 500 pixels
    When user sets properties of histogram viewer:
      | Color Aggr Type     | avg  |
      | Invert Color Scheme | true |
    Then histogram viewer should have repainted by at least 2000 pixels
    And the "bin 1" area of histogram viewer should contain the color "#0000FF"
    And the "bin 17" area of histogram viewer should contain the color "#FF0000"
    When user sets properties of histogram viewer:
      | Invert Color Scheme | false |
      | Color               |       |
    Then histogram viewer should have repainted by at least 2000 pixels
    And "Color Aggr Type" property in context panel should be disabled
    And no errors should have been logged

  Scenario: A split disables color coding and draws lines instead of bars
    When user sets "Color" property of histogram viewer to "HEIGHT"
    And user sets "Split" property of histogram viewer to "RACE"
    Then "Color" property in context panel should be disabled
    And "Color Aggr Type" property in context panel should be disabled
    And "Invert Color Scheme" property in context panel should be disabled
    And histogram viewer should not have a "bin 8" area
    And histogram viewer should have a "line Caucasian" area
    And histogram viewer should have a "line Asian" area
    When user sets "Split" property of histogram viewer to ""
    Then "Color" property in context panel should be enabled
    And histogram viewer should have a "bin 8" area
    When user sets "Color" property of histogram viewer to ""
    Then no errors should have been logged

  Scenario: A split line answers to the pointer as a whole category
    Given user listens for "d4-histogram-mouse-over-line" event on histogram viewer
    And user listens for "d4-histogram-select-line" event on histogram viewer
    When user sets "Split" property of histogram viewer to "SEX"
    Then histogram viewer should have a "line F point" area
    And histogram viewer should have a "line M point" area
    When user hovers over the "line F point" area of histogram viewer
    Then "d4-histogram-mouse-over-line" event should have fired on histogram viewer
    When user clicks on the "line F point" area of histogram viewer
    Then "d4-histogram-select-line" event should have fired on histogram viewer
    And only rows where "SEX" is "F" should be selected
    And 553 rows should be selected
    When user clears the row selection
    And user moves the pointer away from histogram viewer
    And user sets "Split" property of histogram viewer to ""
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Split Stack brings the bars back as one segment per category
    When user sets properties of histogram viewer:
      | Split       | SEX  |
      | Split Stack | true |
    Then histogram viewer should have a "bin 8" area
    And histogram viewer should have a "bin 8 | F" area
    And histogram viewer should have a "bin 8 | M" area
    And the "bin 8 | F" and "bin 8 | M" areas of histogram viewer should be painted in different colors
    And the legend of histogram viewer should list 2 items
    When user sets "Show Values" property of histogram viewer to "true"
    Then histogram viewer should have a "bin labels" area
    And the "bin labels" area of histogram viewer should be painted
    When user sets properties of histogram viewer:
      | Show Values |       |
      | Split Stack | false |
      | Split       |       |
    Then histogram viewer should have a "bin 8" area
    And no errors should have been logged

  Scenario: Normalize Values decides what the vertical axis is scaled to
    Given histogram viewer should have a "y axis" area
    When user sets properties of histogram viewer:
      | Split           | SEX  |
      | Normalize Values | true |
    Then histogram viewer should have repainted by at least 500 pixels
    And histogram viewer should not have a "y axis" area
    When user sets "Normalize Values" property of histogram viewer to "false"
    Then histogram viewer should have a "y axis" area
    And the "y axis max" reading of histogram viewer should be 53
    And histogram viewer should have repainted by at least 500 pixels
    When user sets properties of histogram viewer:
      | Normalize Values | true |
      | Split            |      |
    Then histogram viewer should have a "bin 8" area
    And no errors should have been logged

  Scenario: Distribution lines and markers under a split
    When user sets "Split" property of histogram viewer to "SEX"
    And user sets "Show Distribution Lines" property of histogram viewer to "true"
    Then histogram viewer should have repainted by at least 500 pixels
    And histogram viewer should have more ink than before
    When user sets "Show Markers" property of histogram viewer to "false"
    Then histogram viewer should have less ink than before
    When user sets "Spline Tension" property of histogram viewer to "5"
    Then histogram viewer should have repainted by at least 500 pixels
    When user sets properties of histogram viewer:
      | Spline Tension          | 0     |
      | Show Markers            | true  |
      | Show Distribution Lines | false |
      | Split                   |       |
    Then histogram viewer should have a "bin 8" area
    And no errors should have been logged

  Scenario: A range narrowed under a split keeps the table filter valid
    When user sets properties of histogram viewer:
      | Split             | RACE |
      | Show Range Inputs | true |
    Then all rows should pass the filter
    When user enters "30" into the "range min input" area of histogram viewer
    Then 861 rows should pass the filter
    And the filter should pass exactly the rows where "AGE" is between 30 and 89
    And histogram viewer should show 861 rows
    And the legend of histogram viewer should list 4 items
    And histogram viewer should have a "line Caucasian" area
    When user enters "18" into the "range min input" area of histogram viewer
    Then all rows should pass the filter
    When user sets properties of histogram viewer:
      | Show Range Inputs | false |
      | Split             |       |
    Then histogram viewer should have a "bin 8" area
    And no errors should have been logged

  @realizes:core#998
  Scenario: A split with a single bin can be hovered without errors (core#998)
    When user sets "Split" property of histogram viewer to "SEX"
    Then histogram viewer should have a "line F" area
    When user sets "Bins" property of histogram viewer to "1"
    Then the "bins shown" reading of histogram viewer should be 1
    And histogram viewer should not have a "line F" area
    When user hovers over the "view" area of histogram viewer
    And user moves the pointer away from histogram viewer
    Then the "error" reading of histogram viewer should be ""
    When user sets properties of histogram viewer:
      | Bins  | 20 |
      | Split |    |
    Then the "bins shown" reading of histogram viewer should be 20
    And no errors should have been logged

  @realizes:public#130
  Scenario: Split lines take the split column's category colours (public#130)
    When user colors "SEX" column categorically:
      | F | #00FF00 |
      | M | #FF00FF |
    And user sets "Split" property of histogram viewer to "SEX"
    Then histogram viewer should have a "line F" area
    And the "line F" area of histogram viewer should contain the color "#00FF00"
    And the "line M" area of histogram viewer should contain the color "#FF00FF"
    When user sets "Split" property of histogram viewer to ""
    And user removes the coloring of "SEX" column
    Then no errors should have been logged
