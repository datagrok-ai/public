@journey @viewers @realizes:viewers.pivot-table
Feature: Pivot table — the inner grid that shows the aggregation
  The aggregated grid is the pivot's own: its cells, headers and resizers are regions the pivot
  reports under a "grid " prefix, so the header menu is picked by path and a cell's colour is read
  back from the cell itself rather than from the look. Colour coding paints the value column, Grid >
  Hide takes a column away and the Order or Hide Columns dialog puts it back (GROK-16299), the viewer
  picker of the Aggregate row adds one in-cell viewer column per pivot category (GROK-15004) and none
  at all when two columns are pivoted, and a resizer widens a column. One journey on demog-1000 with
  the automatic cross tab, where None avg(AGE) runs from 36.92 for Indigestion to 52.30 for RA, so a
  linear scheme paints those two cells the ends of its scale.
  Not translated: dragging a key column header inside the inner grid to re-sort the rows (the manual
  checklist's second item). The gesture does not re-sort: it turns the dragged key column into a
  pivot column and the grid then logs a NullError from `col2screen` on every repaint — reported, not
  asserted.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And the "aggregated columns" reading of pivot table viewer should be 6
    And the "text of grid cell 5 of None avg(AGE)" reading of pivot table viewer should be "52.30"

  Scenario: Linear colour coding paints the value column from its own numbers
    Then the "color coding of None avg(AGE)" reading of pivot table viewer should be "Off"
    And the "color of grid cell 5 of None avg(AGE)" reading of pivot table viewer should be "#ffffff"
    When user picks "Grid > Color Coding > Linear" from the context menu of the "grid header None avg(AGE)" area of pivot table viewer
    Then the "color coding of None avg(AGE)" reading of pivot table viewer should be "Linear"
    And the "color of grid cell 5 of None avg(AGE)" reading of pivot table viewer should not be "#ffffff"
    And the "color of grid cell 2 of None avg(AGE)" reading of pivot table viewer should be "#0000ff"
    And the "color of grid cell 5 of None avg(AGE)" and "color of grid cell 2 of None avg(AGE)" readings of pivot table viewer should differ
    And the "grid cell 2 of None avg(AGE)" area of pivot table viewer should contain the color "#0000FF"
    And the "color coding of High avg(AGE)" reading of pivot table viewer should be "Off"
    When user picks "Grid > Color Coding > Off" from the context menu of the "grid header None avg(AGE)" area of pivot table viewer
    Then the "color coding of None avg(AGE)" reading of pivot table viewer should be "Off"
    And the "color of grid cell 5 of None avg(AGE)" reading of pivot table viewer should be "#ffffff"
    And no errors should have been logged

  Scenario: Grid > Hide takes a column away and Order or Hide Columns brings it back
    Then the "column visible of None avg(AGE)" reading of pivot table viewer should be "true"
    And the "columns shown" reading of pivot table viewer should be 7
    When user picks "Grid > Hide" from the context menu of the "grid header None avg(AGE)" area of pivot table viewer
    Then the "column visible of None avg(AGE)" reading of pivot table viewer should be "false"
    And the "columns shown" reading of pivot table viewer should be 6
    And pivot table viewer should not have a "grid header None avg(AGE)" area
    And the "aggregated columns" reading of pivot table viewer should be 6
    When user picks "Grid > Order or Hide Columns..." from the context menu of the "grid header Low avg(AGE)" area of pivot table viewer
    Then "Order or Hide Columns" dialog should be visible
    When user clicks the plain checkbox in the "Order or Hide Columns" dialog
    Then the "column visible of None avg(AGE)" reading of pivot table viewer should be "true"
    And the "columns shown" reading of pivot table viewer should be 7
    When user clicks on CLOSE button in "Order or Hide Columns" dialog
    Then "Order or Hide Columns" dialog should be absent
    And pivot table viewer should have a "grid header None avg(AGE)" area
    And no errors should have been logged

  Scenario: The viewer picker adds one in-cell viewer column per pivot category
    Then the "viewer columns" reading of pivot table viewer should be 0
    And the "pivot" reading of pivot table viewer should be "SEVERITY"
    When user adds a "Scatter plot" viewer column to pivot table viewer
    Then the "viewer columns" reading of pivot table viewer should be 5
    And the "viewer column names" reading of pivot table viewer should be "Critical, High, Low, Medium, None"
    And the "aggregate" reading of pivot table viewer should contain "Scatter plot"
    When user clicks on the "remove aggregate chip Scatter plot" area of pivot table viewer
    Then the "viewer columns" reading of pivot table viewer should be 0
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And no errors should have been logged

  Scenario: With two pivot columns the picker adds no viewer column at all
    When user sets "Pivot Column Names" property of pivot table viewer to "SEVERITY, SEX"
    Then the "aggregated columns" reading of pivot table viewer should be 10
    When user adds a "Scatter plot" viewer column to pivot table viewer
    Then the "aggregate" reading of pivot table viewer should contain "Scatter plot"
    And the "viewer columns" reading of pivot table viewer should be 0
    And the "viewer column names" reading of pivot table viewer should be ""
    When user clicks on the "remove aggregate chip Scatter plot" area of pivot table viewer
    And user sets "Pivot Column Names" property of pivot table viewer to "SEVERITY"
    Then the "aggregated columns" reading of pivot table viewer should be 6
    And the "viewer columns" reading of pivot table viewer should be 0
    And no errors should have been logged

  Scenario: A column resizer of the inner grid widens the column it belongs to
    When user drags the "grid column resizer None avg(AGE)" area of pivot table viewer by 60 pixels to the right
    Then the "column width of None avg(AGE)" reading of pivot table viewer should be higher than before
    And the "text of grid cell 5 of None avg(AGE)" reading of pivot table viewer should be "52.30"
    When user drags the "grid column resizer None avg(AGE)" area of pivot table viewer by 60 pixels to the left
    Then the "column width of None avg(AGE)" reading of pivot table viewer should be lower than before
    And no errors should have been logged
