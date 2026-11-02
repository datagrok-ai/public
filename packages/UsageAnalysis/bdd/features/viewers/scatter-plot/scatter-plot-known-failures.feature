@viewers @realizes:viewers.scatter-plot
Feature: Scatter plot — the known failures
  Temporary: the known failures of the scatter plot while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into scatter-plot-lines.feature as an ordinary @edge-cases scenario.

  A Moving Average Window typed in the context panel outside the 2..100 the property declares (1,
  101) should not be taken; the property keeps its value. The panel's number input marks such a
  value invalid but writes it into the property (core#1200). One Examples row per value, on
  demog-1000 with X = WEIGHT and Y = HEIGHT, as in scatter-plot-lines.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows

  @edge-cases @known-failure:28 @realizes:core#1200
  Scenario Outline: A Moving Average Window of <typed> typed in the context panel is not taken (core#1200)
    Given user sets "Moving Average Window" property of scatter plot viewer to "10"
    When user sets "Show Moving Average Line" property of scatter plot viewer to "true"
    And user clicks on settings icon of scatter plot viewer
    Then context panel should be visible
    Given "Statistics" category in context panel is expanded
    When user enters "<typed>" into "Moving Average Window" property in context panel
    Then "Moving Average Window" property of scatter plot viewer should be "10"

    Examples: The typed value is taken although the property declares 2..100
      | typed |
      | 1     |
      | 101   |
