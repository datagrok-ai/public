@journey @viewers @realizes:viewers.line-chart
Feature: Line chart regression, moving average and formula lines
  The four analytical overlays the chart can draw over its series, and whether they survive a
  layout round-trip through the server.
  Every one of the five overlay scenarios in the old spec asserted exactly one thing: that the
  page's error count had not gone up while the property was being set. Nothing looked at the chart.
  The formula-lines scenario went one better — it wrote two lines into the `formulaLines` property
  as JSON, then parsed *that same JSON string* back out of the property and counted 2, which is a
  test of `JSON.parse`. The chart now counts the fits it drew (`regression lines`,
  `moving average lines`), says whether the deviation band is on, counts the formula items that are
  active, and gives each drawn formula line and band a hit area named after its title — so an item
  that parsed but never rendered has no area.
  Fixture: spgi-100 on `CAST Idea ID` × `Chemical Space X`; `Stereo Category` has 5 categories, so
  a split turns one fit into five.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | CAST Idea ID     |
      | yColumnNames | Chemical Space X |
    Then 100 rows should pass the filter
    And the "regression lines" reading of line chart viewer should be 0
    And the "moving average lines" reading of line chart viewer should be 0
    And the "formula lines" reading of line chart viewer should be 0
    And line chart viewer should not have a "regression line" area
    And line chart viewer should not have a "moving average" area
    And line chart viewer should report no error

  Scenario: The regression line is one fit over the whole series
    When user sets "showRegressionLine" property of line chart viewer to "true"
    Then the "regression lines" reading of line chart viewer should be 1
    And line chart viewer should have a "regression line" area
    And line chart viewer should have more ink than before
    When user sets "showRegressionLine" property of line chart viewer to "false"
    Then the "regression lines" reading of line chart viewer should be 0
    And line chart viewer should not have a "regression line" area
    And line chart viewer should have less ink than before
    And no errors should have been logged

  Scenario: The R key toggles the regression line the chart says it binds
    When user clicks on the "plot" area of line chart viewer
    And user presses R in line chart viewer
    Then the "regression lines" reading of line chart viewer should be 1
    And line chart viewer should have a "regression line" area
    When user presses R in line chart viewer
    Then the "regression lines" reading of line chart viewer should be 0
    And no errors should have been logged

  Scenario: The moving average is a series of its own, and the deviation band widens it
    When user sets "showMovingAverageLine" property of line chart viewer to "true"
    Then the "moving average lines" reading of line chart viewer should be 1
    And the "moving average deviation" reading of line chart viewer should be "false"
    And line chart viewer should have a "moving average" area
    And line chart viewer should have more ink than before
    When user sets "showMovingAverageDeviation" property of line chart viewer to "true"
    Then the "moving average deviation" reading of line chart viewer should be "true"
    And line chart viewer should have more ink than before
    When user sets properties of line chart viewer:
      | showMovingAverageDeviation | false |
      | showMovingAverageLine      | false |
    Then the "moving average lines" reading of line chart viewer should be 0
    And line chart viewer should not have a "moving average" area
    And no errors should have been logged

  Scenario: A split gives every category its own fit and its own moving average
    When user sets properties of line chart viewer:
      | showRegressionLine    | true |
      | showMovingAverageLine | true |
    Then the "regression lines" reading of line chart viewer should be 1
    And the "moving average lines" reading of line chart viewer should be 1
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category"
    Then the "lines" reading of line chart viewer should be 5
    And the "regression lines" reading of line chart viewer should be 5
    And the "moving average lines" reading of line chart viewer should be 5
    And line chart viewer should have repainted
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then the "regression lines" reading of line chart viewer should be 1
    When user sets properties of line chart viewer:
      | showRegressionLine    | false |
      | showMovingAverageLine | false |
    Then no errors should have been logged

  Scenario: A formula line and a formula band each get a hit area named after their title
    When user sets "formulaLines" property of line chart viewer to '[{"type": "line", "formula": "${Chemical Space X} = 5", "title": "const-line", "color": "#FF0000"}, {"type": "band", "formula": "${Chemical Space X} in(2, 8)", "title": "const-band", "color": "#00FF00"}]'
    Then the "formula lines" reading of line chart viewer should be 2
    And line chart viewer should have a 'formula line "const-line"' area
    And line chart viewer should have a 'formula band "const-band"' area
    And line chart viewer should have more ink than before
    And the 'formula band "const-band"' area of line chart viewer should be at least 10 pixels tall
    When user sets "formulaLines" property of line chart viewer to ""
    Then the "formula lines" reading of line chart viewer should be 0
    And line chart viewer should not have a 'formula line "const-line"' area
    And line chart viewer should not have a 'formula band "const-band"' area
    And no errors should have been logged

  Scenario: The formula lines come back from a saved layout
    When user sets "formulaLines" property of line chart viewer to '[{"type": "line", "formula": "${Chemical Space X} = 5", "title": "const-line", "color": "#FF0000"}, {"type": "band", "formula": "${Chemical Space X} in(2, 8)", "title": "const-band", "color": "#00FF00"}]'
    Then the "formula lines" reading of line chart viewer should be 2
    When user saves the layout of the current table view to the server
    And user sets "formulaLines" property of line chart viewer to ""
    Then the "formula lines" reading of line chart viewer should be 0
    When user loads the saved layout
    Then the "formula lines" reading of line chart viewer should be 2
    And line chart viewer should have a 'formula line "const-line"' area
    And line chart viewer should have a 'formula band "const-band"' area
    When user sets "formulaLines" property of line chart viewer to ""
    Then the "formula lines" reading of line chart viewer should be 0
    And no errors should have been logged
