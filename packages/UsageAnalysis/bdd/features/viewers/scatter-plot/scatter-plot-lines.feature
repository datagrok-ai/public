@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot trend lines
  The lines the plot draws over the markers: a regression line, one fit per color category under
  Regression Per Category (on by default) with the equation table, the R key as the shortcut for
  the same property, a regression over a datetime axis with a time unit, and the moving average
  ladder — the line, its window, the split by category and the deviation band. Read from what the
  plot reports: the `regression lines` and `formula lines` readings and the `regression line`,
  `regression stats` and `moving average` hit areas. One journey on demog-1000, X = WEIGHT,
  Y = HEIGHT (872 of the 1000 rows have a HEIGHT); every scenario puts back what it changed.
  The Formula Lines dialog is PowerPack's and has no scenarios here.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows
    And the "formula lines" reading of scatter plot viewer should be 0

  Scenario: A regression line is drawn, on a logarithmic axis too
    Then "Show Regression Line" property of scatter plot viewer should be "false"
    And the "regression lines" reading of scatter plot viewer should be 0
    And scatter plot viewer should not have a "regression line" area
    When user sets "Show Regression Line" property of scatter plot viewer to "true"
    Then the "regression lines" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "regression line" area
    And the "regression line" area of scatter plot viewer should be painted
    When user sets "Y Axis Type" property of scatter plot viewer to "logarithmic"
    Then the "regression lines" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "regression line" area
    And no error or warning balloon should have been shown
    When user sets properties of scatter plot viewer:
      | Y Axis Type          | linear |
      | Show Regression Line | false  |
    Then the "regression lines" reading of scatter plot viewer should be 0
    And scatter plot viewer should not have a "regression line" area
    And no errors should have been logged

  Scenario: Regression Per Category fits one line per color, and the equation table is drawn
    Then "Regression Per Category" property of scatter plot viewer should be "true"
    When user sets properties of scatter plot viewer:
      | Color                | RACE |
      | Show Regression Line | true |
    Then the "regression lines" reading of scatter plot viewer should be 4
    And "Show Regression Line Equation" property of scatter plot viewer should be "true"
    When user sets "Regression Per Category" property of scatter plot viewer to "false"
    Then the "regression lines" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "regression stats" area
    And the "regression stats" area of scatter plot viewer should be painted
    When user sets "Show Regression Line Equation" property of scatter plot viewer to "false"
    Then scatter plot viewer should not have a "regression stats" area
    When user sets properties of scatter plot viewer:
      | Show Regression Line Equation | true |
      | Regression Per Category       | true |
    Then the "regression lines" reading of scatter plot viewer should be 4
    When user sets properties of scatter plot viewer:
      | Show Regression Line | false |
      | Color                |       |
    Then the "regression lines" reading of scatter plot viewer should be 0
    And no errors should have been logged

  Scenario: A regression over a datetime axis takes a time unit
    When user sets properties of scatter plot viewer:
      | X                    | STARTED |
      | Show Regression Line | true    |
    Then the "regression lines" reading of scatter plot viewer should be 1
    And scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    When user sets "X Map" property of scatter plot viewer to "year"
    Then "X Map" property of scatter plot viewer should be "year"
    And scatter plot viewer should have repainted
    And no error or warning balloon should have been shown
    When user sets properties of scatter plot viewer:
      | X Map                |        |
      | X                    | WEIGHT |
      | Show Regression Line | false  |
    Then the "regression lines" reading of scatter plot viewer should be 0
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: The moving average line, its window, its split and its deviation band
    Then scatter plot viewer should not have a "moving average" area
    When user sets properties of scatter plot viewer:
      | Color                       | RACE  |
      | Moving Average Per Category | false |
    And user sets "Show Moving Average Line" property of scatter plot viewer to "true"
    Then scatter plot viewer should have a "moving average" area
    And scatter plot viewer should have more ink than before
    When user sets "Moving Average Window" property of scatter plot viewer to "200"
    Then scatter plot viewer should have repainted
    When user sets "Moving Average Per Category" property of scatter plot viewer to "true"
    Then scatter plot viewer should have repainted
    When user sets "Show Moving Average Deviation" property of scatter plot viewer to "true"
    Then scatter plot viewer should have more ink than before
    When user sets properties of scatter plot viewer:
      | Show Moving Average Deviation | false |
      | Moving Average Window         | 10    |
      | Show Moving Average Line      | false |
      | Color                         |       |
    Then scatter plot viewer should not have a "moving average" area
    And scatter plot viewer should have less ink than before
    And the "formula lines" reading of scatter plot viewer should be 0
    And no errors should have been logged
