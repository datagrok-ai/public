@viewers @realizes:viewers.pivot-table
Feature: Pivot table — a key column header dragged inside the inner grid
  The manual checklist's second inner-grid item (`pivot-table-ui.md`): drag a key column header of
  the aggregated grid onto a value column. On demog-1000 with the automatic cross tab the gesture
  does not re-sort the rows: it turns DIS_POP into a pivot column next to SEVERITY, and the grid
  logs a NullError from `GridCore.col2screen` on every repaint (core#467). The scenario claims only
  that nothing is logged, the known failure; what the drop should do instead is the owner's call.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And the "pivot" reading of pivot table viewer should be "SEVERITY"

  # core#467: dragging a column header inside the pivot's inner grid logs errors
  @known-failure:20
  Scenario: Dragging the DIS_POP header onto None avg(AGE) logs no error
    When user drags the "grid header DIS_POP" area of pivot table viewer onto the "grid header None avg(AGE)" area of pivot table viewer
    Then no errors should have been logged
