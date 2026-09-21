@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart stack aggregation and datetime split
  Stacking needs an additive aggregation: with avg a Stack column draws no segments and no legend,
  sum and count build the stack, removing the Stack collapses it. A datetime Split column enables
  the Split Map, whose year and month re-categorize the bars. One journey on demog-1000 with a bar
  chart of AGE by RACE.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split           | RACE |
      | Value           | AGE  |
      | Value Aggr Type | avg  |
    Then the "bars" reading of bar chart viewer should be 4

  Scenario: A Stack under avg draws no segments and no legend
    When user sets "Stack" property of bar chart viewer to "SEX"
    Then the "stack segments" reading of bar chart viewer should be 0
    And legend of bar chart viewer should be hidden
    And no errors should have been logged

  Scenario: Sum builds the stack and its legend
    When user sets "Value Aggr Type" property of bar chart viewer to "sum"
    Then the "stack segments" reading of bar chart viewer should be 8
    And legend of bar chart viewer should be visible
    And legend of bar chart viewer should have 2 items
    And no errors should have been logged

  Scenario: Count keeps the stack; removing the Stack collapses it
    When user sets "Value Aggr Type" property of bar chart viewer to "count"
    Then the "stack segments" reading of bar chart viewer should be 8
    And legend of bar chart viewer should have 2 items
    When user sets "Stack" property of bar chart viewer to ""
    Then the "stack segments" reading of bar chart viewer should be 0
    And legend of bar chart viewer should be hidden
    And no errors should have been logged

  Scenario: A string Split has no Split Map
    When user clicks on settings icon of bar chart viewer
    Then "Split Map" property in context panel should be disabled
    And the "bars" reading of bar chart viewer should be 4
    And no errors should have been logged

  Scenario: A datetime Split enables the Split Map, which re-categorizes the bars
    When user sets "Split" property of bar chart viewer to "STARTED"
    Then "Split Map" property in context panel should be enabled
    And Split column input in bar chart viewer should contain text "STARTED"
    And Split column input in bar chart viewer should not contain text "RACE"
    When user sets "Split Map" property of bar chart viewer to "month"
    Then the "bars" reading of bar chart viewer should be higher than before
    And bar chart viewer should have repainted
    When user sets "Split Map" property of bar chart viewer to "year"
    Then the "bars" reading of bar chart viewer should be lower than before
    And bar chart viewer should have repainted
    And "Split" property of bar chart viewer should be "STARTED"
    And no errors should have been logged
    When user sets properties of bar chart viewer:
      | Split           | RACE |
      | Value Aggr Type | avg  |
    Then the "bars" reading of bar chart viewer should be 4
