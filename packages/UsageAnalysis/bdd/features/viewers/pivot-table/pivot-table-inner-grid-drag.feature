@viewers @realizes:viewers.pivot-table
Feature: Pivot table — a key column header dragged inside the inner grid
  The manual checklist's second inner-grid item (`pivot-table-ui.md`): drag a key column header of
  the aggregated grid onto a value column. On demog-1000 with the automatic cross tab the gesture
  turned DIS_POP into a pivot column next to SEVERITY, and the grid logged a NullError from
  `GridCore.col2screen` on every repaint (core#467, fixed in the core on 2026-11-10: a header dragged
  inside the grid only reorders its columns). The scenario claims that nothing is logged.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And the "pivot" reading of pivot table viewer should be "SEVERITY"

  @realizes:core#467
  Scenario: Dragging the DIS_POP header onto None avg(AGE) logs no error
    When user drags the "grid header DIS_POP" area of pivot table viewer onto the "grid header None avg(AGE)" area of pivot table viewer
    Then no errors should have been logged
