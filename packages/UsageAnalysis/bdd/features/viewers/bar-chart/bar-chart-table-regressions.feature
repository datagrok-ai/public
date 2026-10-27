@viewers @realizes:viewers.bar-chart
Feature: Bar chart regressions — the Table property switched in the context panel
  After the Table property of a viewer is switched in the context panel, the column properties of the
  same panel offer the new table's columns, not the old ones, without the panel being closed and opened
  again (public#405). A bar chart added to demog-1000 is moved to cars from its gear: the picker of
  its Split property offers cars' one categorical column, model, and not demog-1000's RACE. (On
  spgi-100 the picker of a context-panel column property does not hand the focus back to its selector,
  which the library's reading of the picker waits for, so the small cars table is the second table.)

  Background:
    Given user is logged in
    And user opens cars dataset
    And user opens demog-1000 dataset
    And the context panel is open

  # public#405
  @realizes:public#405
  Scenario: After a Table switch the column properties offer the new table's columns (public#405)
    Given user adds a bar chart viewer
    When user clicks on settings icon of bar chart viewer
    Then the context panel should show "Bar chart"
    When user selects "cars" in "Table" property in context panel
    Then bar chart viewer should be bound to table "cars"
    Given "X Axis" category in context panel is expanded
    Then "Split" property in context panel should offer the columns "model"
    And "Split" property in context panel should not offer the column "RACE"
    And no errors should have been logged
