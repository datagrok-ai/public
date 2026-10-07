@journey @viewers @realizes:viewers.tooltips
Feature: Hiding the table's default tooltip from a viewer's context menu
  Tooltip > Hide on any viewer switches off every tooltip of the table: the grid, the scatter plot
  and the box plot, which show its row tooltip, stop showing anything on hover, and so do the
  histogram's bin and the bar chart's bar, whose own tooltips show before Hide (by design, confirmed
  in the team thread on 2026-10-07; the TestTrack case was updated). The Tooltip group then
  offers Show Custom in place of Hide, and Show Custom brings the row tooltip back on all three and
  the histogram's and the bar chart's tooltips with it.
  Translated from the TestTrack case Tooltips/default-tooltip-visibility on demog-1000, with a scatter
  plot, a box plot, a histogram, a line chart, a bar chart and a trellis plot.

  The grid is switched to "inherit from table" with Show Visible Columns In Tooltip on, as the case
  says: the grid's own default (by design since GROK-19901) is a custom tooltip with no columns,
  which shows nothing whatever Hide does (claimed in grid-visible-columns-tooltip).
  Show Column Names is Always on the grid and on the two plots, so the tooltip's columns are read
  by name.

  The histogram and the bar chart report no hovered bin or bar, so their "no tooltip after Hide" is
  read on the same areas whose tooltips are claimed shown before Hide and after Show Custom.

  Between two hovers the pointer rests above the grid, on the ribbon: above any other viewer of
  this crowded view lies another viewer, whose own tooltip would answer the next claim.

  Not translated: the line chart's and the trellis plot's own tooltips going and coming back — neither reports a
  place a pointer can rest on to raise one (MISSING.md); they are in the view as the case has them.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user sets properties of grid:
      | Show Tooltip                    | inherit from table |
      | Show Column Names               | Always             |
      | Show Visible Columns In Tooltip | true               |
    And user adds a scatter plot viewer
    And user adds a box plot viewer
    And user adds a histogram viewer
    And user adds a line chart viewer
    And user adds a bar chart viewer
    And user adds a trellis plot viewer
    And user sets properties of scatter plot viewer:
      | showLabels | Always |
    And user sets properties of box plot viewer:
      | showLabels | Always |

  Scenario: Before anything is hidden, the grid and both plots show the table's tooltip
    When user hovers over the "cell 11 of AGE" area of grid
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "bin 1" area of histogram viewer
    Then tooltip should be visible
    And tooltip should contain text "AGE"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the first "bar" area of bar chart viewer
    Then tooltip should be visible
    And tooltip should contain text "count"
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Tooltip > Hide on the scatter plot hides the row tooltip on the grid and both plots
    When user opens the context menu of scatter plot viewer
    Then the open menu should list "Tooltip > Hide"
    When user picks "Tooltip > Hide" from the open menu
    And user hovers over the "cell 11 of AGE" area of grid
    Then the mouse-over row of the table should be 11
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
    Then no errors should have been logged

  Scenario: Hide switches off the histogram's and the bar chart's own tooltips too
    When user hovers over the "bin 1" area of histogram viewer
    Then tooltip should be hidden
    When user moves the pointer away from grid
    And user hovers over the first "bar" area of bar chart viewer
    Then tooltip should be hidden
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: The Tooltip group offers Show Custom in place of Hide, and Show Custom brings the tooltip back
    When user opens the context menu of box plot viewer
    Then the open menu should list "Tooltip > Show Custom"
    And the open menu should not list "Tooltip > Hide"
    When user picks "Tooltip > Show Custom" from the open menu
    And user hovers over the "cell 11 of AGE" area of grid
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "bin 1" area of histogram viewer
    Then tooltip should be visible
    And tooltip should contain text "AGE"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the first "bar" area of bar chart viewer
    Then tooltip should be visible
    And tooltip should contain text "count"
    When user moves the pointer away from grid
    Then no errors should have been logged
