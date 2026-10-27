@viewers @realizes:charts.viewer.sunburst
Feature: Sunburst regressions of closed issues
  Closed Sunburst issues, each claimed by the viewer's own report and its picture. public#2292: the
  sunburst did not respond to its Table property; picked in the context panel, spgi-100 rebinds it
  and it draws that table's segments, and demog-1000 picked back brings the first segments back.
  public#2262 (rows selected in the grid drawn on their segments) reproduces on the stand and is
  written up in MISSING.md, not here.

  Background:
    Given user is logged in
    And the package autostarts have completed

  @realizes:public#2292
  Scenario: The Table property moves the sunburst to another table and back (public#2292)
    Given user opens spgi dataset
    And user opens demog-1000 dataset
    And user adds a sunburst viewer
    Then sunburst viewer should be bound to table "demog-1000"
    And sunburst viewer should be painted
    When user remembers the "segment names" reading of sunburst viewer
    And user clicks on settings icon of sunburst viewer
    And user selects "spgi-100" in "Table" property in context panel
    Then sunburst viewer should be bound to table "spgi-100"
    And the "segment names" reading of sunburst viewer should not be as remembered
    And sunburst viewer should report no error
    When user selects "demog-1000" in "Table" property in context panel
    Then sunburst viewer should be bound to table "demog-1000"
    And the "segment names" reading of sunburst viewer should be as remembered
    And no errors should have been logged

