@journey @viewers @realizes:viewers.histogram
Feature: Histogram split selector, bins slider, x axis height and mouse-over row
  The on-plot split selector comes and goes with Show Split Selector; the bins slider and the Bins
  property agree on the bin count it sets (github-2296) — the 169 bins a 40-pixel drag gives depend on
  the slider's pixel scale at this 500-pixel width, so the number holds for this size only; X Axis Height makes the axis box taller (it
  applies with Auto Layout off — Auto Layout sizes the axis itself, histogram_core.dart);
  the mouse-over row marker follows a hover in the grid and Show Mouse Over Row switches it off. One journey on demog-1000 with a 500 by 400
  histogram of AGE (the in-plot selectors need a root over 200 by 200); every scenario puts back
  what it changed.

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

  Scenario: The bins slider and the Bins property agree on the bin count (github-2296)
    When user moves the pointer away from grid
    And user hovers over histogram viewer
    And user drags the "bins slider" area of histogram viewer by 40 pixels to the right
    Then "Bins" property of histogram viewer should be "169"
    And the "bins shown" reading of histogram viewer should be 169
    When user sets "Bins" property of histogram viewer to "20"
    Then the "bins shown" reading of histogram viewer should be 20
    And no errors should have been logged

  Scenario: With Auto Layout off, X Axis Height makes the x axis box taller and gives the height back
    When user sets "Auto Layout" property of histogram viewer to "false"
    And user sets "X Axis Height" property of histogram viewer to "60"
    Then the "x axis" area of histogram viewer should be taller than before
    And the "x axis" area of histogram viewer should be at least 60 pixels tall
    When user sets "X Axis Height" property of histogram viewer to "30"
    Then the "x axis" area of histogram viewer should be shorter than before
    When user sets "Auto Layout" property of histogram viewer to "true"
    Then histogram viewer should have an "x axis" area
    And no errors should have been logged

  Scenario: Show Mouse Over Row decides whether a grid hover puts the marker on the baseline
    Given histogram viewer should not have a "mouse over row marker" area
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
