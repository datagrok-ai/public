@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart stacking and relative values
  Relative Values normalizes stacked bars to one length and is inert without a Stack column: the
  bars keep their absolute lengths and no legend shows; a Stack column set later activates it at
  once, repeatably, and the stacked bars render on a value whose sums include negatives
  (github-2659, GROK-19480). One journey on spgi-100 with a bar chart of the Chemical Space X sum
  by Primary Series Name, stacked by Scaffold Names (six categories, one of them empty).

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a bar chart viewer with:
      | Split           | Primary Series Name |
      | Value           | Chemical Space X    |
      | Value Aggr Type | sum                 |
    Then the "bars" reading of bar chart viewer should be 5
    And the "stack segments" reading of bar chart viewer should be 0
    And "Relative Values" property of bar chart viewer should be "false"

  Scenario: The stacked bars render on negative sums without errors
    When user sets "Stack" property of bar chart viewer to "Scaffold Names"
    Then the "stack segments" reading of bar chart viewer should be higher than before
    And legend of bar chart viewer should be visible
    And legend of bar chart viewer should have 6 items
    And bar chart viewer should have a "bar Triazoles | TRISUBSTITUTED" area
    And the "bar Triazoles | TRISUBSTITUTED" area of bar chart viewer should be painted
    And bar chart viewer should have a "bar Pyrrolidines | AMINOPYRROLIDINE" area
    And the bars of bar chart viewer should differ in length
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Relative Values with a Stack column normalizes the bars
    When user sets "Relative Values" property of bar chart viewer to "true"
    Then bar chart viewer should have repainted by at least 2000 pixels
    And the bars of bar chart viewer should be of equal length
    And the "bar Pyrrolidines | AMINOPYRROLIDINE" area of bar chart viewer should have more ink than before
    And bar chart viewer should be painted
    And no errors should have been logged

  Scenario: Relative Values off restores the absolute lengths
    When user sets "Relative Values" property of bar chart viewer to "false"
    Then bar chart viewer should have repainted by at least 2000 pixels
    And the bars of bar chart viewer should differ in length
    And "Stack" property of bar chart viewer should be "Scaffold Names"
    And no errors should have been logged

  Scenario: Removing the Stack collapses the bars to single segments
    When user sets "Stack" property of bar chart viewer to ""
    Then the "stack segments" reading of bar chart viewer should be 0
    And legend of bar chart viewer should be hidden
    And bar chart viewer should be painted
    And no errors should have been logged

  Scenario: Relative Values alone is inert
    When user sets "Relative Values" property of bar chart viewer to "true"
    Then "Stack" property of bar chart viewer should be ""
    And bar chart viewer should not have repainted
    And the bars of bar chart viewer should differ in length
    And legend of bar chart viewer should be hidden
    And no errors should have been logged

  Scenario: A Stack column activates it at once, and again after it is removed
    When user sets "Stack" property of bar chart viewer to "Scaffold Names"
    Then the bars of bar chart viewer should be of equal length
    And legend of bar chart viewer should have 6 items
    And the "stack segments" reading of bar chart viewer should be higher than before
    When user sets "Stack" property of bar chart viewer to ""
    Then the bars of bar chart viewer should differ in length
    And legend of bar chart viewer should be hidden
    When user sets "Stack" property of bar chart viewer to "Scaffold Names"
    Then the bars of bar chart viewer should be of equal length
    And legend of bar chart viewer should have 6 items
    When user sets properties of bar chart viewer:
      | Stack           |       |
      | Relative Values | false |
    Then "Relative Values" property of bar chart viewer should be "false"
    And the bars of bar chart viewer should differ in length
    And legend of bar chart viewer should be hidden
    And bar chart viewer should be painted
    And no errors should have been logged
