@journey @viewers @realizes:viewers.tooltips
Feature: The grid's tooltip and the columns it cannot show
  With Show Visible Columns In Tooltip off (the default) the grid's row tooltip lists only the
  columns that fall off the grid, so hovering a cell shows nothing while every column fits and
  shows the pushed-out columns once a widened column pushes them past the right edge. Switched on
  in the grid's properties, the tooltip lists every column, the visible ones included, whether or
  not some fall off. Translated from the TestTrack case Tooltips/default-tooltip, on the first rows
  of energy_uk (three columns: value, source, target), written into the feature.

  By design since GROK-19901 (March 2026) the grid's own Show Tooltip defaults to "show custom
  tooltip" with an empty Row Tooltip, and then the grid shows no row tooltip at all, the pushed-out
  columns included; the first scenario claims that default, then the table is opened again and its
  grid switched to "inherit from table" in its properties, as the case now says. The tooltip's Show
  Column Names is set to Always, so the tooltip names the columns it lists (Auto leaves the names
  out of a grid's row tooltip).

  The case's other way of pushing the last column out of sight, widening the context panel, moves
  nothing the grid does not already show by the widened column, and is not translated.

  Background:
    Given user is logged in
    And user opens a table "energy_uk" with:
      | value              | source               | target         |
      | 124.72899627685547 | Agricultural 'waste' | Bio-conversion |
      | 0.597000002861023  | Bio-conversion       | Liquid         |
      | 26.86199951171875  | Bio-conversion       | Losses         |
      | 280.3219909667969  | Bio-conversion       | Solid          |
      | 81.14399719238281  | Bio-conversion       | Gas            |

  Scenario: By default the grid shows no row tooltip, even for the columns pushed off it
    Then "Show Tooltip" property of grid should be "show custom tooltip"
    And "Row Tooltip" property of grid should be ""
    And "Show Visible Columns In Tooltip" property of grid should be "false"
    When user drags the "column resizer value" area of grid by 1000 pixels to the right
    Then grid should not have a "cell 1 of target" area
    When user hovers over the "cell 1 of value" area of grid
    Then the mouse-over row of the table should be 1
    And tooltip should be hidden
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: The grid takes the table's tooltip, and Show Visible Columns In Tooltip is off
    When user closes all views
    And user opens a table "energy_uk" with:
      | value              | source               | target         |
      | 124.72899627685547 | Agricultural 'waste' | Bio-conversion |
      | 0.597000002861023  | Bio-conversion       | Liquid         |
      | 26.86199951171875  | Bio-conversion       | Losses         |
      | 280.3219909667969  | Bio-conversion       | Solid          |
      | 81.14399719238281  | Bio-conversion       | Gas            |
    Then "Show Tooltip" property of grid should be "show custom tooltip"
    When user hovers over the "cell 1 of value" area of grid
    And user clicks on "Edit properties (F4)" icon in grid
    Then context panel should be visible
    Given "Tooltip" category in context panel is expanded
    When user selects "inherit from table" in "Show Tooltip" property in context panel
    Then "Show Tooltip" property of grid should be "inherit from table"
    When user selects "Always" in "Show Column Names" property in context panel
    Then "Show Column Names" property of grid should be "Always"
    And "Show Visible Columns In Tooltip" property in context panel should be unchecked
    And no errors should have been logged

  Scenario: Off, no tooltip while every column fits
    Then the "column order" reading of grid should be "value, source, target"
    And grid should have a "cell 1 of target" area
    When user hovers over the "cell 1 of value" area of grid
    Then the mouse-over row of the table should be 1
    And tooltip should be hidden
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Switched on, the tooltip lists every column while every column fits
    When user checks "Show Visible Columns In Tooltip" property in context panel
    Then "Show Visible Columns In Tooltip" property of grid should be "true"
    When user hovers over the "cell 1 of value" area of grid
    Then tooltip should be visible
    And the tooltip should show columns "value, source, target"
    And the tooltip should show "source" as "Agricultural 'waste'"
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Switched on, the tooltip stays the same with columns pushed off the grid
    When user remembers the "column width of value" reading of grid
    And user drags the "column resizer value" area of grid by 1000 pixels to the right
    Then the "column width of value" reading of grid should be higher than remembered
    And grid should not have a "cell 1 of source" area
    And grid should not have a "cell 1 of target" area
    When user hovers over the "cell 1 of value" area of grid
    Then tooltip should be visible
    And the tooltip should show columns "value, source, target"
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Off again, the tooltip lists only the columns pushed off the grid
    When user unchecks "Show Visible Columns In Tooltip" property in context panel
    Then "Show Visible Columns In Tooltip" property of grid should be "false"
    When user hovers over the "cell 1 of value" area of grid
    Then tooltip should be visible
    And the tooltip should show columns "source, target"
    And the tooltip should show "target" as "Bio-conversion"
    When user moves the pointer away from grid
    Then no errors should have been logged
