@viewers @realizes:viewers.trellis-plot
Feature: Trellis plot — the known failures
  Temporary: the known failures of the trellis plot while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into trellis-plot-tiles-and-layout.feature as an ordinary @edge-cases scenario.

  A Tiles Per Row typed in the context panel outside the property's declared bounds (1 to 10; 11
  here) should not be taken; the property keeps its value. The panel's number input marks such a
  value invalid but writes it into the property (core#1200). On demog-1000 with a trellis split by
  DIS_POP alone, as in trellis-plot-tiles-and-layout.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names  | DIS_POP      |
      | Y Column Names  |              |
      | Viewer Type     | Scatter plot |
      | Pack Categories | false        |
    Then the "x categories" reading of trellis plot viewer should be 6

  @edge-cases @known-failure:29 @realizes:core#1200
  Scenario: A Tiles Per Row of 11 typed in the context panel outside its declared bounds is not taken (core#1200)
    Given user sets "Tiles Per Row" property of trellis plot viewer to "4"
    When user clicks on settings icon of trellis plot viewer
    Then context panel should be visible
    Given "Misc" category in context panel is expanded
    When user enters "11" into "Tiles Per Row" property in context panel
    Then "Tiles Per Row" property of trellis plot viewer should be "4"
