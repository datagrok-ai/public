@journey @viewers @realizes:viewers.filters
Feature: Hierarchical filter card
  The hierarchical card filters through a tree of levels: a click on a row keeps that branch alone,
  a branch opens to the next level whose rows carry the count of the rows they still pass, the
  card's search hides nodes without touching the table (GROK-19968), reordering the levels rebuilds
  the roots and clears the criterion (GROK-16528), and a saved layout brings the levels and the
  criterion back. One journey on demog-1000: 553 rows are F, 480 of them Caucasian, and 285 of
  those have no SEVERITY.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user picks "Add Filter | Hierarchical" from the filter panel menu
    And user configures the hierarchical filter with columns "SEX, RACE"
    Then "SEX / RACE" filter card should be visible
    And the filter panel should have 1 filter
    And all rows should pass the filter
    And counter of filter panel should be hidden

  Scenario: A click on a branch keeps only its rows
    Then the hierarchical filter card should list the "F" row
    And the "F" row of the hierarchical filter card should count 553 rows
    When user clicks on the "F" row of the hierarchical filter card
    Then 553 rows should pass the filter
    And the filter should pass exactly the rows where "SEX" is "F"
    And the "filters" reading of filter panel should be 1
    And counter of filter panel should have text "1"
    And the "M" row of the hierarchical filter card should count 0 rows
    And no errors should have been logged

  Scenario: A branch opens to the next level and a click there narrows to the leaf
    When user expands the "F" row of the hierarchical filter card
    Then the hierarchical filter card should list the "F / Caucasian" row
    And the hierarchical filter card should list the "F / Other" row
    When user clicks on the "F / Caucasian" row of the hierarchical filter card
    Then 480 rows should pass the filter
    And the "F / Caucasian" row of the hierarchical filter card should count 480 rows
    And the "F / Other" row of the hierarchical filter card should count 0 rows
    And counter of filter panel should have text "1"
    And no errors should have been logged

  Scenario: A third level narrows inside the leaf of the second
    When user configures the hierarchical filter with columns "SEX, RACE, SEVERITY"
    Then "SEX / RACE / SEVERITY" filter card should be visible
    And all rows should pass the filter
    When user expands the "F" row of the hierarchical filter card
    And user expands the "F / Caucasian" row of the hierarchical filter card
    Then the hierarchical filter card should list the "F / Caucasian / None" row
    When user clicks on the "F / Caucasian / None" row of the hierarchical filter card
    Then 285 rows should pass the filter
    And the "F / Caucasian / None" row of the hierarchical filter card should count 285 rows
    And the "F / Caucasian / Low" row of the hierarchical filter card should count 0 rows
    And counter of filter panel should have text "1"
    And no errors should have been logged

  Scenario: Reordering the levels rebuilds the roots and clears the criterion
    When user configures the hierarchical filter with columns "RACE, SEX"
    Then "RACE / SEX" filter card should be visible
    And the hierarchical filter card should list the "Caucasian" row
    And the hierarchical filter card should not list the "F" row
    And all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: The card's search hides rows and leaves the table alone
    When user clicks on the "Asian" row of the hierarchical filter card
    Then 15 rows should pass the filter
    When user hovers over "RACE / SEX" filter card
    And user clicks on search icon of "RACE / SEX" filter card
    And user types "Cau" into the search of the "RACE / SEX" filter card
    Then the hierarchical filter card should not list the "Asian" row
    And the hierarchical filter card should list the "Caucasian" row
    And 15 rows should pass the filter
    And no rows where "RACE" is "Caucasian" should pass the filter
    When user clears the search of the "RACE / SEX" filter card
    Then the hierarchical filter card should list the "Asian" row
    And 15 rows should pass the filter
    And no errors should have been logged

  Scenario: A saved layout brings the levels and the criterion back
    When user configures the hierarchical filter with columns "SEX, RACE"
    And user clicks on the "F" row of the hierarchical filter card
    Then 553 rows should pass the filter
    When user saves the layout of the current table view to the server
    And user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    When user loads the saved layout
    Then "SEX / RACE" filter card should be visible
    And the filter panel should have 1 filter
    And 553 rows should pass the filter
    And counter of filter panel should have text "1"
    And no errors should have been logged
