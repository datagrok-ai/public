@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart sorting and orientation
  Vertical bars sorted by value descend from left to right and swap sides when the order flips;
  stacking holds on a value whose sums include negatives and the negative bar hangs below the
  baseline error-free (GROK-19480); vertical bars with small counts still fill the view
  (github-3417). One journey on spgi-100 with a bar chart of the Chemical Space X sum by Primary
  Series Name.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a bar chart viewer with:
      | Split           | Primary Series Name |
      | Value           | Chemical Space X    |
      | Value Aggr Type | sum                 |
    Then the "bars" reading of bar chart viewer should be 5
    And the bars of bar chart viewer should lie one under another

  Scenario: Vertical and descending by value puts the tallest bar on the left
    When user sets properties of bar chart viewer:
      | Orientation    | vertical |
      | Bar Sort Type  | by value |
      | Bar Sort Order | desc     |
    Then bar chart viewer should have repainted
    And the bars of bar chart viewer should descend from left to right
    And the "bar Triazoles" area of bar chart viewer should hang below the "bar Aminopiperidines" area
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Stacking holds on a negative-sum value
    When user sets "Legend Visibility" property of bar chart viewer to "Always"
    Then legend of bar chart viewer should be hidden
    When user sets "Stack" property of bar chart viewer to "Stereo Category"
    Then legend of bar chart viewer should be visible
    And legend of bar chart viewer should have 5 items
    And the "stack segments" reading of bar chart viewer should be higher than before
    And bar chart viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user sets properties of bar chart viewer:
      | Stack             |      |
      | Legend Visibility | Auto |
    Then legend of bar chart viewer should be hidden
    And the "stack segments" reading of bar chart viewer should be 0
    And no errors should have been logged

  Scenario: Ascending swaps the tall side, horizontal puts the bars back under each other
    When user sets "Bar Sort Order" property of bar chart viewer to "asc"
    Then bar chart viewer should have repainted
    And the bars of bar chart viewer should ascend from left to right
    And the "bar Triazoles" area of bar chart viewer should hang below the "bar Aminopiperidines" area
    When user sets "Orientation" property of bar chart viewer to "horizontal"
    Then bar chart viewer should have repainted
    And the bars of bar chart viewer should lie one under another
    And "Bar Sort Order" property of bar chart viewer should be "asc"
    And no errors should have been logged

  Scenario: Vertical bars with small counts fill the view
    When user sets properties of bar chart viewer:
      | Value           | CAST Idea ID |
      | Value Aggr Type | count        |
    And user sets "Orientation" property of bar chart viewer to "vertical"
    Then the tallest bar of bar chart viewer should start within the top 40% of the view
    And the "bar Triazoles" area of bar chart viewer should be painted
    And no errors should have been logged
    When user sets "Orientation" property of bar chart viewer to "horizontal"
    Then bar chart viewer should have repainted
    And the bars of bar chart viewer should lie one under another
    And no errors should have been logged
    When user sets "Orientation" property of bar chart viewer to "auto"
    Then properties of bar chart viewer should be:
      | Orientation     | auto         |
      | Value           | CAST Idea ID |
      | Value Aggr Type | count        |
    And the bars of bar chart viewer should lie one under another
    And bar chart viewer should be painted
    And no errors should have been logged
    When user sets properties of bar chart viewer:
      | Value           | Chemical Space X |
      | Value Aggr Type | sum              |
      | Bar Sort Type   | by category      |
      | Bar Sort Order  | asc              |
