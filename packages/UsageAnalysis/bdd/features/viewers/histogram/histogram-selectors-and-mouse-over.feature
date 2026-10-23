@journey @viewers @realizes:viewers.histogram
Feature: Histogram split selector and mouse-over row
  The on-plot split selector comes and goes with Show Split Selector; the mouse-over row marker
  follows a hover in the grid and Show Mouse Over Row switches it off. One journey on demog-1000
  with a 500 by 400 histogram of AGE (the in-plot selectors need a root over 200 by 200); every
  scenario puts back what it changed. The bins slider against the Bins property (github-2296) is in
  `histogram.feature`; X Axis Height is in MISSING.md (core#1086).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | Value | AGE |
    And user resizes histogram viewer to 500 by 400
    Then the "bins shown" reading of histogram viewer should be 20

  Scenario: Show Split Selector takes the split selector off the plot and brings it back
    When user hovers over histogram viewer
    Then Split column input in histogram viewer should be visible
    When user sets "Show Split Selector" property of histogram viewer to "false"
    And user moves the pointer away from grid
    And user hovers over histogram viewer
    Then Split column input in histogram viewer should be hidden
    And Value column input in histogram viewer should be visible
    When user sets "Show Split Selector" property of histogram viewer to "true"
    And user moves the pointer away from grid
    And user hovers over histogram viewer
    Then Split column input in histogram viewer should be visible
    And no errors should have been logged

  Scenario: Show Mouse Over Row decides whether a grid hover puts the marker on the baseline
    Then histogram viewer should not have a "mouse over row marker" area
    When user takes a snapshot of histogram viewer
    And user hovers over the "cell 8 of AGE" area of grid
    Then histogram viewer should have a "mouse over row marker" area
    And the "mouse over row marker" area of histogram viewer should have repainted
    When user moves the pointer away from grid
    Then histogram viewer should not have a "mouse over row marker" area
    When user sets "Show Mouse Over Row" property of histogram viewer to "false"
    And user takes a snapshot of histogram viewer
    And user hovers over the "cell 8 of AGE" area of grid
    Then the mouse-over row of the table should be 8
    And histogram viewer should not have a "mouse over row marker" area
    And histogram viewer should not have repainted
    When user moves the pointer away from grid
    And user sets "Show Mouse Over Row" property of histogram viewer to "true"
    Then no errors should have been logged
