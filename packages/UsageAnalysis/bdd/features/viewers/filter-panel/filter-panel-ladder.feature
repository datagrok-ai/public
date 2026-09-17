@journey @viewers @realizes:viewers.filters
Feature: Filter panel core ladder
  The ladder every other filter panel feature stands on: a card added from the header picker
  filters nothing, a click on a category name keeps that category alone while a click on its
  checkbox adds one, a second criterion composes with the first, the header counter counts the
  cards that restrict rows and its tooltip names them, the master toggle stashes and gives back
  every card's state, Escape on the panel does the same, the header search hides cards and nothing
  else, the reset icon clears the criteria and keeps the cards, and a card's own checkbox suspends
  its criterion without losing it. One journey on demog-1000: RACE is Black 27, Other 62,
  Caucasian 896, Asian 15, and 708 of the 1000 rows are aged 30 to 60 — 19 of them Black.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    And counter of filter panel should be hidden

  Scenario: A card added from the header picker filters nothing
    When user adds a card for "RACE" to the filter panel
    Then "RACE" filter card should be visible
    And the "cards" reading of filter panel should be "RACE"
    And the "type of RACE" reading of filter panel should be "categorical"
    And the "categories of RACE" reading of filter panel should be "Asian, Black, Caucasian, Other"
    And all rows should pass the filter
    And the "filters" reading of filter panel should be 0
    And the "filtering of RACE" reading of filter panel should be "false"
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: A click on a category name keeps that category alone
    When user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    And the filter should pass exactly the rows where "RACE" is "Black"
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And the "rows shown" reading of filter panel should be 27
    And the "filters" reading of filter panel should be 1
    And counter of filter panel should have text "1"
    And no errors should have been logged

  Scenario: A click on a checkbox adds a category to the ones kept
    When user clicks on the "checkbox Other of RACE" area of filter panel
    Then 89 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Black, Other"
    And counter of filter panel should have text "1"
    When user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And no errors should have been logged

  Scenario: A second criterion composes with the first
    When user adds a range filter on "AGE" from 30 to 60
    Then "AGE" filter card should be visible
    And 19 rows should pass the filter
    And the "min of AGE" reading of filter panel should be 30
    And the "max of AGE" reading of filter panel should be 60
    And the "filters" reading of filter panel should be 2
    And counter of filter panel should have text "2"
    And no errors should have been logged

  Scenario: The counter's tooltip names the cards that restrict rows
    When user adds a card for "SEX" to the filter panel
    Then the "filters" reading of filter panel should be 2
    When user hovers over counter of filter panel
    Then tooltip should contain the text "RACE"
    And tooltip should contain the text "Black"
    And tooltip should contain the text "AGE"
    And tooltip should not contain the text "SEX"
    And no errors should have been logged

  Scenario: The master toggle stashes every card's state and gives it back
    When user hovers over filter panel
    And user unchecks master of filter panel
    Then all rows should pass the filter
    And the "active" reading of filter panel should be "false"
    And the "cards" reading of filter panel should be "SEX, AGE, RACE"
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And counter of filter panel should be hidden
    When user checks master of filter panel
    Then 19 rows should pass the filter
    And the "active" reading of filter panel should be "true"
    And counter of filter panel should have text "2"
    And no errors should have been logged

  Scenario: The header search hides cards and leaves the rows alone
    When user hovers over filter panel
    And user clicks on search icon of filter panel
    And user types "RACE" into search of filter panel
    Then "AGE" filter card should be hidden
    And "SEX" filter card should be hidden
    And "RACE" filter card should be visible
    And 19 rows should pass the filter
    And counter of filter panel should have text "2"
    When user clears search of filter panel
    Then "AGE" filter card should be visible
    And "SEX" filter card should be visible
    And 19 rows should pass the filter
    And no errors should have been logged

  Scenario: The reset icon clears the criteria and keeps the cards
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the "filters" reading of filter panel should be 0
    And counter of filter panel should be hidden
    And the "cards" reading of filter panel should be "SEX, AGE, RACE"
    And "RACE" filter card should be enabled
    And the "filtering of RACE" reading of filter panel should be "false"
    And no errors should have been logged

  Scenario: A card's checkbox suspends its criterion and keeps it
    When user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user hovers over "RACE" filter card
    And user unchecks checkbox of "RACE" filter card
    Then all rows should pass the filter
    And "RACE" filter card should be disabled
    And the "enabled of RACE" reading of filter panel should be "false"
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And counter of filter panel should be hidden
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then "RACE" filter card should be enabled
    And the "enabled of RACE" reading of filter panel should be "true"
    And all rows should pass the filter
    And no errors should have been logged
