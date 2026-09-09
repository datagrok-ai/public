@journey @viewers @realizes:viewers.filters
Feature: Filter panel persistence
  What a saved layout and a saved project carry: the cards, the criterion each one holds and the
  rows they leave — and what they must not carry, the header search that was open when the layout
  was saved (GROK-16677). A project reopened from the server brings its filter panel back with the
  same criterion (GROK-19152). One journey on demog-1000, where RACE is Caucasian in 896 rows and
  633 of those are aged 30 to 60.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    Then the filter panel should have 0 filters
    And all rows should pass the filter

  Scenario: A saved layout brings the cards and their criteria back
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Caucasian of RACE" area of filter panel
    And user adds a range filter on "AGE" from 30 to 60
    Then 633 rows should pass the filter
    And the "filters" reading of filter panel should be 2
    When user saves the layout of the current table view to the server
    And user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the "filters" reading of filter panel should be 0
    When user loads the saved layout
    Then filter panel should be visible
    And "RACE" filter card should be visible
    And "AGE" filter card should be visible
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And the "min of AGE" reading of filter panel should be 30
    And the "max of AGE" reading of filter panel should be 60
    And 633 rows should pass the filter
    And no errors should have been logged

  Scenario: A layout saved while the header search was open carries no residue
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    When user clicks on search icon of filter panel
    And user types "RACE" into search of filter panel
    Then "AGE" filter card should be hidden
    And "RACE" filter card should be visible
    When user saves the layout of the current table view to the server
    And user clears search of filter panel
    And user clicks on the "category Asian of RACE" area of filter panel
    Then 15 rows should pass the filter
    When user loads the saved layout
    Then "RACE" filter card should be visible
    And "AGE" filter card should be visible
    And all rows should pass the filter
    And the "filtering of RACE" reading of filter panel should be "false"
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: A project reopens with its filter panel and its criterion
    When user clicks on the "category Caucasian of RACE" area of filter panel
    Then 896 rows should pass the filter
    And counter of filter panel should have text "1"
    When user saves the current view as project "bdd filter panel round trip"
    And user closes all views
    And user opens the "bdd filter panel round trip" project
    Then filter panel should be visible
    And "RACE" filter card should be visible
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And 896 rows should pass the filter
    And no errors should have been logged
