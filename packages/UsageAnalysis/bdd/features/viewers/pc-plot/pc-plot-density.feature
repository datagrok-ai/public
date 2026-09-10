@journey @viewers @realizes:viewers.pc-plot
Feature: PC plot density overlay
  Show Density draws a distribution per axis inside the box the plot reports as `density <col>`: a
  circle cloud, a box plot or a violin, each with its own parts — the interquartile range, the
  median, the mean cross, the two whisker dashes and, on the box plot, the circles — and a bin
  count the violin's outline follows. Every claim here is the ink of that one box, not of the
  canvas, so a change to one axis's shape is read where it was drawn.
  The journey runs with Show All Lines off, so the 1000 polylines do not paint over the boxes: the
  old spec had to build that sparse canvas to see the overlay at all, and it stays because a box
  drawn over a thousand overlapping lines adds no measurable ink.
  demog-1000 has exactly 1000 rows and `autoStyle` sets `showDensity = rowCount > 1000`, so a fresh
  plot has NO density and `density style` reads "" — the old spec's claim that the plot opens with
  density on holds only on the larger demog.csv. The style default (circles) survives regardless.
  Every scenario puts the density settings back.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pc plot viewer with:
      | Column Names   | AGE, HEIGHT, WEIGHT |
      | Show All Lines | false               |
    Then pc plot viewer should show 1000 rows
    And the "lines drawn" reading of pc plot viewer should be 1
    And "Show Density" property of pc plot viewer should be "false"
    And the "density style" reading of pc plot viewer should be ""
    And pc plot viewer should not have a "density \"AGE\"" area

  Scenario: The density style defaults to circles and Show Density draws them
    Then "Density Style" property of pc plot viewer should be "circles"
    When user sets "Show Density" property of pc plot viewer to "true"
    Then the "density style" reading of pc plot viewer should be "circles"
    And pc plot viewer should have a "density \"AGE\"" area
    And pc plot viewer should have a "density \"HEIGHT\"" area
    And pc plot viewer should have a "density \"WEIGHT\"" area
    And the "density \"AGE\"" area of pc plot viewer should be painted
    And pc plot viewer should have more ink than before
    When user sets "Show Density" property of pc plot viewer to "false"
    Then pc plot viewer should not have a "density \"AGE\"" area
    And pc plot viewer should have less ink than before
    And no errors should have been logged

  Scenario: Each density style draws a different shape in the same box
    When user sets "Show Density" property of pc plot viewer to "true"
    Then the "density \"AGE\"" area of pc plot viewer should be painted
    When user sets "Density Style" property of pc plot viewer to "box plot"
    Then the "density style" reading of pc plot viewer should be "box plot"
    And the "density \"AGE\"" area of pc plot viewer should have repainted
    And the "density \"AGE\"" area of pc plot viewer should be painted
    When user sets "Density Style" property of pc plot viewer to "violin plot"
    Then the "density style" reading of pc plot viewer should be "violin plot"
    And the "density \"AGE\"" area of pc plot viewer should have repainted
    And the "density \"AGE\"" area of pc plot viewer should be painted
    When user sets "Density Style" property of pc plot viewer to "circles"
    Then the "density \"AGE\"" area of pc plot viewer should have repainted
    And no errors should have been logged

  Scenario: The box plot's parts each add their own ink
    When user sets "Density Style" property of pc plot viewer to "box plot"
    Then the "density \"AGE\"" area of pc plot viewer should be painted
    When user sets "Show Circles" property of pc plot viewer to "false"
    Then the "density \"AGE\"" area of pc plot viewer should have less ink than before
    When user sets "Show Interquartile Range" property of pc plot viewer to "false"
    Then the "density \"AGE\"" area of pc plot viewer should have less ink than before
    When user sets "Show Median" property of pc plot viewer to "false"
    Then the "density \"AGE\"" area of pc plot viewer should have less ink than before
    When user sets "Show Mean Cross" property of pc plot viewer to "false"
    Then the "density \"AGE\"" area of pc plot viewer should have less ink than before
    When user sets "Show Upper Dash" property of pc plot viewer to "false"
    Then the "density \"AGE\"" area of pc plot viewer should have less ink than before
    When user sets properties of pc plot viewer:
      | Show Upper Dash          | true |
      | Show Mean Cross          | true |
      | Show Median              | true |
      | Show Interquartile Range | true |
      | Show Circles             | true |
    Then the "density \"AGE\"" area of pc plot viewer should have more ink than before
    And no errors should have been logged

  Scenario: The bin count changes the violin's outline
    When user sets "Density Style" property of pc plot viewer to "violin plot"
    Then "Bins" property of pc plot viewer should be "100"
    When user sets "Bins" property of pc plot viewer to "200"
    Then "Bins" property of pc plot viewer should be "200"
    And the "density \"AGE\"" area of pc plot viewer should have repainted
    When user sets "Bins" property of pc plot viewer to "20"
    Then the "density \"AGE\"" area of pc plot viewer should have repainted
    When user sets "Bins" property of pc plot viewer to "100"
    Then the "density \"AGE\"" area of pc plot viewer should have repainted
    And no errors should have been logged

  Scenario: The density survives a normalization double-toggle (github-1546)
    When user sets "Density Style" property of pc plot viewer to "box plot"
    And user sets "Normalize Each Column" property of pc plot viewer to "false"
    Then the "normalization" reading of pc plot viewer should be "global"
    And the "density \"AGE\"" area of pc plot viewer should be painted
    When user sets "Normalize Each Column" property of pc plot viewer to "true"
    And user sets "Normalize Each Column" property of pc plot viewer to "false"
    And user sets "Normalize Each Column" property of pc plot viewer to "true"
    Then the "normalization" reading of pc plot viewer should be "per column"
    And the "density style" reading of pc plot viewer should be "box plot"
    And the "density \"AGE\"" area of pc plot viewer should be painted
    And the "error" reading of pc plot viewer should be ""
    And no errors should have been logged

  Scenario: A logarithmic AGE axis keeps the density painted
    When user sets "Log Columns" property of pc plot viewer to "AGE"
    Then "Log Columns" property of pc plot viewer should be "AGE"
    And pc plot viewer should have repainted
    And the "density \"AGE\"" area of pc plot viewer should be painted
    And the "density \"WEIGHT\"" area of pc plot viewer should be painted
    When user sets "Log Columns" property of pc plot viewer to ""
    Then "Log Columns" property of pc plot viewer should be ""
    And the "density \"AGE\"" area of pc plot viewer should be painted
    And pc plot viewer should have repainted
    When user sets properties of pc plot viewer:
      | Show Density   | false    |
      | Density Style  | circles  |
      | Show All Lines | true     |
    Then the "density style" reading of pc plot viewer should be ""
    And pc plot viewer should not have a "density \"AGE\"" area
    And the "lines drawn" reading of pc plot viewer should be 1000
    And no errors should have been logged
