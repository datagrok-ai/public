@viewers @realizes:viewers.histogram
Feature: Histogram known failures
  Temporary: the known failures of the histogram while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into histogram.feature as an ordinary @edge-cases scenario.

  Each scenario is a test of its own on demog-1000 with a histogram of AGE.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | Value | AGE |
    Then histogram viewer should show 1000 rows
    And the "bins shown" reading of histogram viewer should be 20

  @edge-cases @known-failure:23 @realizes:core#1200
  Scenario: A Bins value typed in the context panel outside its declared bounds is not taken (core#1200)
    When user clicks on settings icon of histogram viewer
    Then context panel should be visible
    Given "X Axis" category in context panel is expanded
    When user enters "201" into "Bins" property in context panel
    Then "Bins" property of histogram viewer should be "20"
