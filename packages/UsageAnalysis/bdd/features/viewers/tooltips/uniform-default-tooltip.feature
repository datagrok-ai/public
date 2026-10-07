@journey @viewers @realizes:viewers.tooltips
Feature: Viewers show the same default tooltip
  A viewer that inherits the table's tooltip shows the same columns as every other one: a marker of
  the scatter plot, a marker of the box plot and a cell of the grid all list the table's default
  tooltip columns for their row. Translated from the TestTrack case Tooltips/uniform-default-tooltip
  on demog-1000.

  Each viewer's Show Column Names is set to Always, so the tooltip names the columns it lists. The
  grid is switched to "inherit from table" with Show Visible Columns In Tooltip on, as the case
  says: its own default (by design since GROK-19901) is a custom tooltip with no columns, which
  shows nothing (claimed in grid-visible-columns-tooltip), and without the visible columns it lists
  only the columns off the screen.

  The grid is reconfigured, so what is shared is the table's tooltip, not the three viewers' own
  defaults; edit-tooltip claims the same sharing for a tooltip edited in the dialog.

  Not translated: that the columns come in the same order on every viewer — the tooltip steps claim
  the set of columns, not their order, and the scatter plot puts its axis columns first
  (MISSING.md).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer
    And user adds a box plot viewer
    And user sets properties of scatter plot viewer:
      | showLabels | Always |
    And user sets properties of box plot viewer:
      | showLabels | Always |
    And user sets properties of grid:
      | Show Tooltip                    | inherit from table |
      | Show Column Names               | Always             |
      | Show Visible Columns In Tooltip | true               |

  Scenario: The scatter plot, the box plot and the grid list the same columns
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then tooltip should be visible
    And the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "marker" area of box plot viewer
    Then tooltip should be visible
    And the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "cell 11 of AGE" area of grid
    Then tooltip should be visible
    And the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then no errors should have been logged
