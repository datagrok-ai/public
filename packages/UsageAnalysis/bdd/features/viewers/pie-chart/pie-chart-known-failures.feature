@viewers @realizes:viewers.pie-chart
Feature: Pie chart — the known failures
  Temporary: the known failures of the pie chart while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into pie-chart.feature as an ordinary @edge-cases scenario.

  A value typed in the context panel outside the property's declared bounds should not be taken;
  the property keeps its value: Max Radius 600 (150 stays), Start Angle -30 and 400 (0 stays). The
  panel's number input marks such a value invalid but writes it into the property (core#1200). One
  Examples row per value, on demog-1000 with a pie chart of RACE, as in pie-chart.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pie chart viewer with:
      | Category | RACE |
    Then the "slices" reading of pie chart viewer should be 4

  @edge-cases @known-failure:26 @realizes:core#1200
  Scenario Outline: <caption> <typed> typed in the context panel outside its declared bounds is not taken (core#1200)
    Given user sets "<caption>" property of pie chart viewer to "<kept>"
    When user clicks on settings icon of pie chart viewer
    Then context panel should be visible
    Given "Misc" category in context panel is expanded
    When user enters "<typed>" into "<caption>" property in context panel
    Then "<caption>" property of pie chart viewer should be "<kept>"

    Examples: The value is written into the property although the input marks it invalid
      | caption     | typed | kept |
      | Max Radius  | 600   | 150  |
      | Start Angle | -30   | 0    |
      | Start Angle | 400   | 0    |
