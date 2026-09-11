@journey @eda @realizes:eda.viewer.pareto-front
Feature: Pareto front objectives
  The objectives of the Pareto front viewer, chosen in its properties in the context panel. Translated
  from files/TestTrack/EDA/pareto-front-viewer.md, steps 1 to 4 and 7, and the package's
  playwright/pareto-front-viewer.test.ts.

  Minimize and Maximize each open the platform's column picker, which lists the columns the property
  offers in the table's order. cars has sixteen numeric columns after its string model column, so a
  picker that offered strings would start with model. The old spec read the offer from the property's
  `choices`, which the viewer never sets: its claim that model and turbo were absent passed whatever
  the offer was. The conflict warning is text the viewer puts in its own element, not a picture on
  its canvas. The viewer reports no readings of its own (no getWidgetStatus): an axis or a label chosen
  by hand is claimed by the automatic choice the viewer turns off.

  One journey over one viewer, on purpose: for two seconds after a property is edited the platform
  ignores a change of the current object (AppEvents.propertyEdited), so a feature that edited a
  property and a next one that opens the settings of a new viewer at once would edit the old one.

  The case's cars-with-missing.csv is on no stand and in no repository. An int column whose every
  value is null stands in for its empty turbo; the property offers every numeric column, empty or
  not, while the viewer computes only over the non-empty ones, so that scenario is @known-failure
  until the offer leaves the empty column out.

  Background:
    Given user is logged in
    And user opens cars dataset
    And the context panel is open
    When user picks "ML > Pareto Front..." from the top menu
    Then pareto front viewer should be visible

  Scenario: The viewer's properties come in their categories, and only numeric columns are offered
    When user clicks on settings icon of pareto front viewer
    Then "Objectives" category in context panel should be visible
    And "Axes" category in context panel should be visible
    And "Labels" category in context panel should be visible
    And "Legend" category in context panel should be visible
    And "Description" category in context panel should be visible
    Given "Objectives" category in context panel is expanded
    Then "Minimize" property in context panel should contain text "2 / 16"
    And "Maximize" property in context panel should contain text "0 / 16"
    When user clicks on "..." button in "Maximize" property in context panel
    Then "Select columns..." dialog should be visible
    And the "rows" reading of grid viewer in "Select columns..." dialog should be 16
    And the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "diesel"
    When user clicks on CANCEL button in "Select columns..." dialog
    Then "Maximize" property of pareto front viewer should be ""
    And no errors should have been logged

  Scenario: Maximizing what is minimized is refused with a warning, and choosing again clears it
    When user clicks on "..." button in "Maximize" property in context panel
    And user clicks on All label in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Maximize" property in context panel should contain text "16 / 16"
    And "Maximize" property of pareto front viewer should contain "price"
    And pareto front viewer should contain text "Cannot minimize and maximize features at the same time"
    And pareto front viewer should contain text "highway.mpg"
    When user clicks on "..." button in "Maximize" property in context panel
    And user clicks on None label in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Maximize" property of pareto front viewer should be ""
    And pareto front viewer should not contain text "Cannot minimize and maximize"
    And no errors should have been logged

  Scenario: An axis and the labels chosen by hand turn their automatic choice off
    Given "Axes" category in context panel is expanded
    When user selects "horsepower" in "X Axis" property in context panel
    Then "X Axis" property of pareto front viewer should be "horsepower"
    And "Auto Axes Selection" property of pareto front viewer should not be "true"
    Given "Labels" category in context panel is expanded
    When user clicks on "..." button in "Label Columns" property in context panel
    And user clicks on None label in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Label Columns" property of pareto front viewer should be ""
    And "Auto Labels Selection" property of pareto front viewer should not be "true"
    And no errors should have been logged

  @known-failure
  Scenario: An empty column is not offered as an objective
    When user adds a calculated column "empty" with formula "If(true, null, 0)"
    Then "empty" column should have type "int"
    And "empty" column should have missing values
    When user clicks on "..." button in "Maximize" property in context panel
    Then the "rows" reading of grid viewer in "Select columns..." dialog should be 16
