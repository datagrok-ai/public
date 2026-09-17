@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot split columns and the inner viewer
  The grid the split columns build and the viewer that fills every cell of it: a clean add, the
  eight cells of SEX by RACE, each built-in inner type picked in the control panel and announcing
  itself, an excursion through a no-data state that must leave no stale frame, the control panel
  that hides without touching the type, gridlines under the three modes and two inner types, and an
  inner-viewer setting that redraws every cell. One journey on demog-1000; every scenario puts back
  what it changed.

  Not translated here: the (+) column picker's hover preview and its blank row (the popup is a
  canvas grid with no DOM rows, so both can only be driven by a measured layout model — the split
  columns are set through the property instead), and the Multi Curve and To Script menu entries,
  which a plain stand does not carry.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8

  Scenario: The viewer is added clean
    Then the open tableview should have 1 trellis plot viewer
    And the "error" reading of trellis plot viewer should be ""
    And the "cells drawn" reading of trellis plot viewer should be 8
    And the "blank cells" reading of trellis plot viewer should be 0
    And trellis plot viewer should show 1000 rows
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Two split columns give eight cells
    Then the "x categories" reading of trellis plot viewer should be 2
    And the "y categories" reading of trellis plot viewer should be 4
    And the cells of trellis plot viewer should be 2 wide and 4 tall
    And the "x categories shown" reading of trellis plot viewer should be 2
    And the "y categories shown" reading of trellis plot viewer should be 4
    And trellis plot viewer should have a "cell F | Caucasian" area
    And trellis plot viewer should have a "cell M | Asian" area
    And trellis plot viewer should have a "cell 1,1" area
    And trellis plot viewer should not have a "cell F | Klingon" area
    And the "distinct cell signatures" reading of trellis plot viewer should be 8
    And no errors should have been logged

  Scenario Outline: <type> picked in the control panel announces itself and draws a picture per cell
    Given user listens for "d4-trellis-plot-viewer-type-changed" event on trellis plot viewer
    When user picks "<type>" in the viewer selector of trellis plot viewer
    Then "d4-trellis-plot-viewer-type-changed" event should have fired on trellis plot viewer
    And the "inner viewer type" reading of trellis plot viewer should be "<type>"
    And the "cells drawn" reading of trellis plot viewer should be 8
    And the "distinct cell signatures" reading of trellis plot viewer should be at least <pictures>
    And no errors should have been logged

    Examples:
      | type         | pictures |
      | Bar chart    | 8        |
      | Histogram    | 8        |
      | Line chart   | 8        |
      | Pie chart    | 2        |
      | Box plot     | 8        |
      | Density plot | 2        |
      | Scatter plot | 8        |

  Scenario: A no-data excursion leaves no stale frame
    When user sets "Viewer Type" property of trellis plot viewer to "Pie chart"
    Then the "cells drawn" reading of trellis plot viewer should be 8
    And the "distinct cell signatures" reading of trellis plot viewer should be at least 2
    When user sets properties of trellis plot viewer:
      | Y Column Names |           |
      | Viewer Type    | Bar chart |
    Then the "y categories" reading of trellis plot viewer should be 1
    When user sets properties of trellis plot viewer:
      | Y Column Names | RACE      |
      | Viewer Type    | Pie chart |
    Then the "cells drawn" reading of trellis plot viewer should be 8
    And the "distinct cell signatures" reading of trellis plot viewer should be at least 2
    And the "blank cells" reading of trellis plot viewer should be 0
    And the "cell signature F | Caucasian" and "cell signature M | Asian" readings of trellis plot viewer should differ
    And no errors should have been logged
    When user sets "Viewer Type" property of trellis plot viewer to "Scatter plot"
    Then the "inner viewer type" reading of trellis plot viewer should be "Scatter plot"

  Scenario: Hiding the control panel keeps the type
    Given user sets "Viewer Type" property of trellis plot viewer to "Scatter plot"
    And trellis plot viewer should have a "control panel" area
    When user sets "Show Control Panel" property of trellis plot viewer to "false"
    Then trellis plot viewer should not have a "control panel" area
    And the "inner viewer type" reading of trellis plot viewer should be "Scatter plot"
    And the "cells drawn" reading of trellis plot viewer should be 8
    When user sets "Show Control Panel" property of trellis plot viewer to "true"
    Then trellis plot viewer should have a "control panel" area
    And no errors should have been logged

  Scenario: Gridlines follow the mode and the inner type
    When user sets "Show Gridlines" property of trellis plot viewer to "always"
    Then the "gridlines" reading of trellis plot viewer should be "true"
    When user sets "Show Gridlines" property of trellis plot viewer to "never"
    Then the "gridlines" reading of trellis plot viewer should be "false"
    When user sets properties of trellis plot viewer:
      | Show Gridlines | auto         |
      | Viewer Type    | Scatter plot |
    Then the "gridlines" reading of trellis plot viewer should be "true"
    When user sets "Viewer Type" property of trellis plot viewer to "Bar chart"
    Then the "gridlines" reading of trellis plot viewer should be "false"
    When user sets properties of trellis plot viewer:
      | Viewer Type    | Scatter plot |
      | Show Gridlines | always       |
    Then the "gridlines" reading of trellis plot viewer should be "true"
    And no errors should have been logged

  Scenario: Changing an inner viewer setting redraws every cell
    When user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user sets "xColumnName" inner property of trellis plot viewer to "AGE"
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should not be as remembered
    And the "cell signature M | Caucasian" reading of trellis plot viewer should differ from before
    And the "cells drawn" reading of trellis plot viewer should be 8
    And the "blank cells" reading of trellis plot viewer should be 0
    When user sets "xColumnName" inner property of trellis plot viewer to "HEIGHT"
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should be as remembered
    And no errors should have been logged
