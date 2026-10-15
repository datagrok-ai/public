@journey @viewers @realizes:viewers.tooltips
Feature: A viewer's own tooltip properties
  A viewer's Tooltip category holds Show Tooltip — "inherit from table" by default, "show custom
  tooltip" or "do not show", the latter two picked here — and Row Tooltip, greyed out and empty until
  the viewer shows a custom tooltip. Tooltip > Hide on a viewer with a custom tooltip sets it to
  "do not show" and hides its tooltip alone: the reference scatter plot keeps the table's tooltip
  and the box plot its own custom one. Tooltip > Show Custom brings the same tooltip back, and so
  does "show custom tooltip" after "do not show" in the properties. Hide on the reference viewer
  switches off every tooltip of the table — the reference's and both custom ones (by design,
  confirmed in the team thread on 2026-10-07; the TestTrack case was updated) — and Show Custom on
  it brings them all back. Translated from the TestTrack
  case Tooltips/tooltip-properties on demog-1000, with two scatter plots (the second keeps the table's
  tooltip, for reference) and a box plot. Both scatter plots are added before any property
  changes: a viewer added later copies the settings of the one of its type already in the view.

  Show Column Names is Always on the plots, so the tooltip's columns are read by name. Between two
  hovers the pointer rests above the grid.

  A custom tooltip with an empty Row Tooltip lists only the viewer's own data columns, not the
  table's tooltip columns (by design, confirmed 2026-10-07; the TestTrack case was updated): the
  scatter plot its axis columns (HEIGHT and WEIGHT, from Data Values = Merge), the box plot nothing,
  so the box plot is then given its own Row Tooltip (AGE, SEX). The grid, whose custom tooltip with
  an empty Row Tooltip shows nothing by design (claimed in grid-visible-columns-tooltip), is given
  its own Row Tooltip too (RACE, SEX) and shows exactly those columns.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer
    And user adds a scatter plot viewer
    And user adds a box plot viewer
    And user sets properties of scatter plot viewer:
      | showLabels | Always |
    And user sets properties of second scatter plot viewer:
      | showLabels | Always |
    And user sets properties of box plot viewer:
      | showLabels | Always |

  Scenario: Show Tooltip inherits from the table, and Row Tooltip is greyed out and empty
    When user clicks on settings icon of first scatter plot viewer
    Then context panel should be visible
    Given "Tooltip" category in context panel is expanded
    Then "Show Tooltip" property in context panel should have value "inherit from table"
    And "Row Tooltip" property in context panel should be disabled
    And "Row Tooltip" property of scatter plot viewer should be ""
    And no errors should have been logged

  Scenario: A custom tooltip with no columns of its own lists the viewer's own data columns
    When user selects "show custom tooltip" in "Show Tooltip" property in context panel
    Then "Show Tooltip" property of scatter plot viewer should be "show custom tooltip"
    And "Row Tooltip" property in context panel should be enabled
    And "Show Tooltip" property of second scatter plot viewer should be "inherit from table"
    And properties of scatter plot viewer should be:
      | X           | HEIGHT |
      | Y           | WEIGHT |
      | Data Values | Merge  |
    When user sets "Show Tooltip" property of box plot viewer to "show custom tooltip"
    And user hovers over the "marker" area of box plot viewer
    Then the mouse-over row of the table should be 952
    And tooltip should be hidden
    When user moves the pointer away from grid
    And user sets "Row Tooltip" property of box plot viewer to "AGE\nSEX"
    And user sets properties of grid:
      | Show Column Names | Always    |
      | Row Tooltip       | RACE\nSEX |
    Then "Show Tooltip" property of grid should be "show custom tooltip"
    When user hovers over the "cell 11 of AGE" area of grid
    Then the tooltip should show columns "RACE, SEX"
    When user moves the pointer away from grid
    Then "Show Tooltip" property of second scatter plot viewer should be "inherit from table"
    When user hovers over the "marker of row 11" area of second scatter plot viewer
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns "AGE, SEX"
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Tooltip > Hide on a viewer with a custom tooltip hides that tooltip alone
    When user picks "Tooltip > Hide" from the context menu of scatter plot viewer
    Then "Show Tooltip" property of scatter plot viewer should be "do not show"
    When user opens the context menu of scatter plot viewer
    Then the open menu should list "Tooltip > Show Custom"
    And the open menu should not list "Tooltip > Hide"
    When user closes the context menu
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the "hovered row" reading of scatter plot viewer should be 11
    And tooltip should be hidden
    When user moves the pointer away from grid
    And user hovers over the "marker of row 11" area of second scatter plot viewer
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns "AGE, SEX"
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Tooltip > Show Custom brings the same custom tooltip back
    When user picks "Tooltip > Show Custom" from the context menu of scatter plot viewer
    Then "Show Tooltip" property of scatter plot viewer should be "show custom tooltip"
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: "do not show" in the properties hides the custom tooltip too
    When user clicks on settings icon of first scatter plot viewer
    Given "Tooltip" category in context panel is expanded
    When user selects "do not show" in "Show Tooltip" property in context panel
    Then "Show Tooltip" property of scatter plot viewer should be "do not show"
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the "hovered row" reading of scatter plot viewer should be 11
    And tooltip should be hidden
    When user moves the pointer away from grid
    And user selects "show custom tooltip" in "Show Tooltip" property in context panel
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Hide on the reference viewer switches off every tooltip, and Show Custom brings them all back
    When user picks "Tooltip > Hide" from the context menu of second scatter plot viewer
    And user hovers over the "marker of row 11" area of second scatter plot viewer
    Then the "hovered row" reading of second scatter plot viewer should be 11
    And tooltip should be hidden
    When user moves the pointer away from grid
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the "hovered row" reading of scatter plot viewer should be 11
    And tooltip should be hidden
    When user moves the pointer away from grid
    And user hovers over the "marker" area of box plot viewer
    Then the mouse-over row of the table should be 952
    And tooltip should be hidden
    When user moves the pointer away from grid
    And user opens the context menu of second scatter plot viewer
    Then the open menu should list "Tooltip > Show Custom"
    When user picks "Tooltip > Show Custom" from the open menu
    And user hovers over the "marker of row 11" area of second scatter plot viewer
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    And user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns "AGE, SEX"
    When user moves the pointer away from grid
    Then no errors should have been logged
